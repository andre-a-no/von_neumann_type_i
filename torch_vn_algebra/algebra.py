"""
Type I von Neumann algebra M = ⊕_{c=1}^{C} M_{n_c}(ℂ) with subspace dimensions k_c.

Storage shape: (batch, C, k_max, k_max)
- batch: parallel Monte Carlo samples
- C: number of factors (channels)
- k_max = max(k_c): maximum subspace dimension across channels

Each channel c has an active block of size k_c × k_c in the top-left corner.
The remaining (k_max - k_c) dimensions are padding.
"""

import torch
import numpy as np
from typing import Optional, List, Tuple, Callable, Union
from .hilbert_space import HilbertSpace
from . import cost
from . import functional as _fc
from .functional import cache_enabled as _fc_cache_enabled
from .functional import divided_differences as _divided_differences   # noqa: F401 (tests)

def _as_tensor(x) -> torch.Tensor:
    """torch.as_tensor that keeps Python and NumPy numbers in double precision (torch.as_tensor turns a list of
    Python floats into float32, which would round data meant for a double-precision algebra)."""
    if isinstance(x, torch.Tensor):
        return x
    try:
        return torch.as_tensor(np.asarray(x))
    except (TypeError, ValueError, RuntimeError):       # e.g. a list of tensors that require gradients
        return torch.as_tensor(x)


def _same_device(a, b) -> bool:
    """Device equality that treats 'cpu' and torch.device('cpu'), 'cuda' and 'cuda:<current>' as equal."""
    a, b = torch.device(a), torch.device(b)
    if a.type != b.type:
        return False
    if a.type == 'cuda':
        cur = torch.cuda.current_device() if torch.cuda.is_available() else 0
        return (a.index if a.index is not None else cur) == (b.index if b.index is not None else cur)
    return True


def _gram_schmidt_preferred(Z: torch.Tensor) -> bool:
    """On a GPU, torch.linalg.qr factorises a batch of small matrices one matrix at a time (tens of
    microseconds each); Gram-Schmidt costs O(n) batched kernels regardless of the batch, so it is used
    when the batch is large compared with n. On the CPU, LAPACK QR is fast and is kept."""
    return Z.is_cuda and Z.dim() == 3 and Z.shape[-1] <= 256 and Z.shape[0] >= Z.shape[-1] \
        and Z.shape[-2] <= 4096


def _orthonormal_columns(Z: torch.Tensor) -> torch.Tensor:
    """
    Q of the thin QR decomposition Z = QR with diag(R) > 0, for a batch of matrices of full column rank, by
    classical Gram-Schmidt with one reorthogonalisation ("twice is enough"; orthogonality to working
    precision). With diag(R) > 0 the factor Q is unique, so for a Ginibre Z it is Haar distributed
    (the same matrix as Mezzadri's QR recipe with the phase correction).
    """
    n = Z.shape[-1]
    if Z.requires_grad and torch.is_grad_enabled():      # differentiable variant: no in-place writes
        cols = []
        for j in range(n):
            v = Z[..., j:j + 1]
            if j:
                P = torch.cat(cols, dim=-1)
                v = v - P @ (P.mH @ v)
                v = v - P @ (P.mH @ v)
            cols.append(v / torch.linalg.vector_norm(v, dim=-2, keepdim=True))
        return torch.cat(cols, dim=-1)
    Q = torch.empty_like(Z)
    for j in range(n):
        v = Z[..., j:j + 1]
        if j:
            P = Q[..., :j]
            v = v - P @ (P.mH @ v)
            v = v - P @ (P.mH @ v)
        Q[..., j:j + 1] = v / torch.linalg.vector_norm(v, dim=-2, keepdim=True)
    return Q


class TypeIAlgebra:
    """
    Type I von Neumann algebra M = ⊕_{c=1}^{C} M_{n_c}(ℂ).
    
    Operator storage shape: (batch, C, k_max, k_max)
    """

    def __init__(
        self,
        n_factors: List[int],
        k_factors: List[int],
        hilbert: Optional[HilbertSpace] = None,
        batch_size: int = 1,
        complex_valued: bool = True,
        device: Optional[torch.device] = None,
        precision: str = 'single',
        charges: Optional[List] = None
    ):
        """
        n_factors, k_factors: sizes n_c of the factors and k_c <= n_c of the active subspaces.
        precision: 'single' (float32 / complex64, default) or 'double' (float64 / complex128).
        charges: labels of the sectors (e.g. particle numbers); default 0, 1, ..., C-1. They are
            used to merge sectors in symmetric tensor products (composite.FusedProduct).
        """
        C = len(n_factors)
        assert len(k_factors) == C, "n_factors and k_factors must have same length"

        for n, k in zip(n_factors, k_factors):
            if not (n > 0 and 0 <= k <= n):
                raise ValueError(f"need 0 <= k_c <= n_c and n_c > 0, got n={n}, k={k}")
        if C == 0 or max(k_factors) == 0:
            raise ValueError("the algebra needs at least one sector with k_c > 0")

        self.n_factors = list(n_factors)
        self.k_factors = list(k_factors)
        k_factors = self.k_factors
        self.C = C
        self.k_max = max(k_factors) if k_factors else 0
        self.total_subspace_dim = sum(k_factors)
        # kept for compatibility; lambda_max / lambda_min are now exact for every block size
        self.exact_eig_max_dim = 256
        self.charges = list(range(C)) if charges is None else list(charges)
        assert len(self.charges) == C, "one charge per sector"

        # padded storage: warn when most of the (C, k_max, k_max) layout is padding
        padded = C * self.k_max ** 2
        used = sum(k * k for k in k_factors)
        if padded > 2 ** 24 and used < 0.5 * padded:
            import warnings
            warnings.warn(f"sector sizes are very uneven: {100 * (1 - used / padded):.0f}% of every operator "
                          f"(C x k_max^2 = {padded:.2e} entries) is padding; consider working in single sectors "
                          f"(e.g. SpinChain.sector_algebra).", cost.CostWarning, stacklevel=2)

        if hilbert is None:
            self.hilbert = HilbertSpace(
                n=self.total_subspace_dim,
                k=self.total_subspace_dim,
                batch_size=batch_size,
                n_channels=C,
                complex_valued=complex_valued,
                device=device,
                precision=precision
            )
        else:
            self.hilbert = hilbert

        # sectors grouped by block size, with an index usable on the stored tensors: a slice when the sectors of
        # a group are contiguous (a view, no copy), otherwise a cached device index tensor
        groups = {}
        for c, k_c in enumerate(self.k_factors):
            if k_c:
                groups.setdefault(k_c, []).append(c)
        self._size_groups_cache = list(groups.items())
        self._const_cache = {}

    def _sel(self, idx: List[int]):
        """Index for the sector axis selecting the sectors idx (sorted): a slice if they are contiguous, else a
        cached index tensor on the algebra's device (built once, so indexing does not copy a list to the device)."""
        if idx[-1] - idx[0] + 1 == len(idx):
            return slice(idx[0], idx[-1] + 1)
        key = ('sel', tuple(idx))
        t = self._const_cache.get(key)
        if t is None:
            t = self._const_cache[key] = torch.tensor(idx, dtype=torch.long, device=self.hilbert.device)
        return t

    def _const(self, name: str, dtype: Optional[torch.dtype] = None) -> torch.Tensor:
        """Constant tensors of the sector structure, built once per (name, dtype) on the algebra's device:
        'k' (C,) block sizes; 'rows' (C, k_max) active rows; 'block_mask' (C, k_max, k_max); 'eye' (C, k_max, k_max)
        identity of every block; 'w_norm' (C,) 1/k_c; 'w_vn' (C,) 1/(C k_c) (zero for empty sectors)."""
        dtype = dtype or self.dtype
        key = (name, dtype)
        t = self._const_cache.get(key)
        if t is not None:
            return t
        dev = self.hilbert.device
        k = torch.tensor(self.k_factors, device=dev)
        rows = torch.arange(self.k_max, device=dev) < k[:, None]
        if name == 'k':
            t = k.to(dtype)
        elif name == 'rows':
            t = rows
        elif name == 'block_mask':
            t = rows[:, :, None] & rows[:, None, :]
        elif name == 'eye':
            t = torch.diag_embed(rows.to(dtype))
        elif name in ('w_norm', 'w_vn'):
            kk = k.to(dtype)
            w = torch.where(k > 0, 1 / kk.clamp(min=1), torch.zeros_like(kk))
            t = w if name == 'w_norm' else w / sum(1 for kc in self.k_factors if kc)
        else:
            raise KeyError(name)
        self._const_cache[key] = t
        return t

    @property
    def complex_valued(self) -> bool:
        return self.hilbert.complex_valued

    @property
    def precision(self) -> str:
        return self.hilbert.precision

    @property
    def dtype(self) -> torch.dtype:
        return self.hilbert.dtype

    def like(self, complex_valued: Optional[bool] = None, precision: Optional[str] = None,
             device=None) -> 'TypeIAlgebra':
        """The same sector structure with another number field / precision / device."""
        alg = TypeIAlgebra(self.n_factors, self.k_factors,
                           complex_valued=self.complex_valued if complex_valued is None else complex_valued,
                           precision=self.precision if precision is None else precision,
                           device=self.hilbert.device if device is None else device, charges=self.charges)
        alg.exact_eig_max_dim = self.exact_eig_max_dim
        return alg

    def bytes_per_operator(self, batch_size: int = 1) -> int:
        """Memory of one materialised operator of this algebra (padded layout)."""
        return batch_size * self.C * self.k_max ** 2 * self.dtype.itemsize

    # ========================================================================
    # Random Unitary / Orthogonal / SU(n)
    # ========================================================================

    UNITARY_MEASURES = ('haar', 'haar_su', 'coe', 'cse', 'diag')

    def random_isometry(self, n: int, m: int, batch_size: Optional[int] = None) -> torch.Tensor:
        """
        Haar-random isometry C^m -> C^n (R^m -> R^n for a real algebra), m <= n: the first m columns of a Haar
        unitary, drawn directly as the Q factor of an n x m Ginibre matrix with diag(R) > 0.
        Returns a tensor of shape (n, m), or (batch_size, n, m) if batch_size is given.
        """
        if not 0 < m <= n:
            raise ValueError(f"need 0 < m <= n, got n={n}, m={m}")
        shape = () if batch_size is None else (batch_size,)
        Z = torch.randn(*shape, n, m, dtype=self.hilbert.dtype, device=self.hilbert.device)
        if _gram_schmidt_preferred(Z):
            return _orthonormal_columns(Z)
        Q, R = cost.batched_call(torch.linalg.qr, Z, f"random_isometry(n={n}, m={m})", kind='qr') if Z.dim() == 3 \
            else torch.linalg.qr(Z)
        d = torch.diagonal(R, dim1=-2, dim2=-1)
        d = d.sgn() if self.hilbert.complex_valued else d.sign()
        return Q * d.unsqueeze(-2)

    def random_unitary(self, n: int, measure: str = 'haar',
                       batch_size: Optional[int] = None) -> torch.Tensor:
        """
        Random n x n unitary (orthogonal if the algebra is real).

        measure:
            'haar'    - Haar measure on U(n) (complex) or O(n) (real)
            'haar_su' - Haar measure on SU(n) (complex) or SO(n) (real)
            'coe'     - Circular Orthogonal Ensemble, S = W^T W (complex only)
            'cse'     - Circular Symplectic Ensemble, S = W^R W with
                        W^R = J W^T J^T (complex only, n even)
            'diag'    - diagonal matrix of i.i.d. uniform phases (complex only)

        Returns a tensor of shape (n, n), or (batch_size, n, n) if batch_size is given.
        """
        shape = () if batch_size is None else (batch_size,)
        device = self.hilbert.device
        complex_ = self.hilbert.complex_valued
        if measure in ('coe', 'cse', 'diag') and not complex_:
            raise ValueError(f"measure '{measure}' requires a complex-valued algebra")

        if measure == 'haar':
            # Mezzadri's recipe: QR of a Ginibre matrix with the phases of diag(R) fixed.
            # torch.randn with a complex dtype already draws circular N(0, 1) entries.
            dtype = self.hilbert.dtype
            cost.check_memory(3 * cost.tensor_bytes((*shape, n, n), dtype), device, f"random_unitary(n={n})")
            Z = torch.randn(*shape, n, n, dtype=dtype, device=device)
            if _gram_schmidt_preferred(Z):
                return _orthonormal_columns(Z)
            Q, R = cost.batched_call(torch.linalg.qr, Z, f"random_unitary(n={n}, batch={shape[0] if shape else 1})",
                                     kind='qr')
            d = torch.diagonal(R, dim1=-2, dim2=-1)
            d = d.sgn() if complex_ else d.sign()
            return Q * d.unsqueeze(-2)
        elif measure == 'haar_su':
            Q = self.random_unitary(n, 'haar', batch_size)
            det = torch.linalg.det(Q)
            if complex_:
                # Multiply by det^{-1/n}: a Haar U(n) matrix times a phase is Haar on SU(n).
                phase = torch.exp(-1j * torch.angle(det) / n).to(Q.dtype)
                return Q * phase[..., None, None]
            # Real case: flip the first column when det = -1 to land in SO(n).
            s_ = torch.ones(*det.shape, n, dtype=Q.dtype, device=Q.device)
            s_[..., 0] = det.sign()
            return Q * s_.unsqueeze(-2)
        elif measure == 'coe':
            W = self.random_unitary(n, 'haar', batch_size)
            return W.transpose(-2, -1) @ W
        elif measure == 'cse':
            if n % 2 != 0:
                raise ValueError("CSE requires even dimension")
            W = self.random_unitary(n, 'haar', batch_size)
            h = n // 2
            Wt = W.transpose(-2, -1)                     # W^R = J W^T J^T with J = [[0, 1], [-1, 0]]: a signed
            W_dual = torch.cat([torch.cat([Wt[..., h:, h:], -Wt[..., h:, :h]], dim=-1),      # block permutation
                                torch.cat([-Wt[..., :h, h:], Wt[..., :h, :h]], dim=-1)], dim=-2)
            return W_dual @ W
        elif measure == 'diag':
            phases = torch.exp(2j * torch.pi * torch.rand(*shape, n, device=device, dtype=self.hilbert.real_dtype))
            return torch.diag_embed(phases.to(self.hilbert.dtype))
        else:
            raise ValueError(f"Unknown measure: {measure}. Choose from {self.UNITARY_MEASURES}")

    # ========================================================================
    # Operator Class
    # ========================================================================

    class Operator:
        """
        Operator in the von Neumann algebra.
        
        Storage: (batch, C, k_max, k_max)
        """

        def __init__(
            self,
            algebra: 'TypeIAlgebra',
            generator: Optional[Callable[[], torch.Tensor]] = None,
            matrix: Optional[torch.Tensor] = None,
            is_self_adjoint: Optional[bool] = None,
            is_normal: Optional[bool] = None,
            is_positive: Optional[bool] = None,
            is_invertible: Optional[bool] = None,
            is_projection: Optional[bool] = None,
            eigenvalues: Optional[List[torch.Tensor]] = None
        ):
            self.algebra = algebra
            self._generator = generator
            self._matrix = matrix
            self._is_materialized = matrix is not None

            self._is_self_adjoint = is_self_adjoint
            self._is_self_adjoint_set = is_self_adjoint is not None
            self._is_normal = is_normal
            self._is_normal_set = is_normal is not None
            self._is_positive = is_positive
            self._is_positive_set = is_positive is not None
            self._is_invertible = is_invertible
            self._is_invertible_set = is_invertible is not None
            self._is_projection = is_projection
            self._is_projection_set = is_projection is not None
            self._eigenvalues = eigenvalues

            self._lambda_max = None
            self._lambda_min = None
            self._trace = None
            self._inverse = None
            self._abs = None
            self._spectral_cache = None          # {'kind', 'key', 'groups': [(w, V), ...]} when caching is on
            self._cache_on = False
            self._batch = None if matrix is None else matrix.shape[0]   # known batch size of a lazy operator

            self._update_memory()

        def _update_memory(self):
            if self._matrix is not None:
                self._bytes_held = self._matrix.element_size() * self._matrix.numel()
            else:
                self._bytes_held = 0

        @property
        def matrix(self) -> torch.Tensor:
            if not self._is_materialized:
                self._matrix = self._generator()
                self._is_materialized = True
                # drop the recipe: it holds references to the parent operators (and their
                # cached matrices), which would otherwise stay alive as long as this one
                self._generator = None
                self._validate_properties()
                self._update_memory()
            return self._matrix

        def _eps(self) -> float:
            return 1e-5 if self._matrix.dtype in (torch.float32, torch.complex64) else 1e-10

        def _tol_t(self) -> torch.Tensor:
            """Tolerance for property checks as a 0-dim tensor on the device (no synchronisation): eps times the
            largest |entry|, i.e. relative to the size of the operator (an absolute floor would misclassify
            operators of small norm, e.g. 1e-6 times a nilpotent as self-adjoint)."""
            mat = self._matrix
            scale = mat.abs().amax() if mat.numel() else mat.new_zeros((), dtype=mat.real.dtype)
            return self._eps() * torch.clamp(scale, min=torch.finfo(scale.dtype).tiny)

        def _tol(self) -> float:
            """The same tolerance as a Python number (one host synchronisation)."""
            return float(self._tol_t())

        def _max_abs(self, t: torch.Tensor) -> float:
            return t.abs().max().item() if t.numel() else 0.0

        def _validate_properties(self):
            """Check the property flags passed by the user (one host synchronisation for all of them): run when a
            lazy operator is materialised, and by alg.operator(..., validate=True); flags of an eager operator are
            otherwise trusted. Flags set by the library itself are mathematically guaranteed and not checked."""
            checks = []
            want = lambda flag: getattr(self, f'_is_{flag}_set') and getattr(self, f'_is_{flag}')
            if not any(want(f) for f in ('self_adjoint', 'normal', 'positive', 'invertible', 'projection')):
                return
            mat = self._matrix
            if mat.numel() == 0:
                return
            tol = self._tol_t()
            mat_d = mat.conj().transpose(-2, -1)
            if want('self_adjoint'):
                checks.append(("Operator marked self-adjoint but not Hermitian", (mat - mat_d).abs().amax() > tol))
            if want('normal'):
                checks.append(("Operator marked normal but not normal",
                               (mat @ mat_d - mat_d @ mat).abs().amax() > tol * mat.abs().amax()))
            if want('positive'):
                checks.append(("Operator marked positive but has negative eigenvalues", (self.lambda_min < -tol).any()))
            if want('invertible'):
                checks.append(("Operator marked invertible but has zero eigenvalue",
                               (self._block_svdvals_min() <= tol).any()))
            if want('projection'):
                checks.append(("Operator marked projection but P^2 != P", (mat @ mat - mat).abs().amax() > tol))
            failed = torch.stack([c for _, c in checks]).tolist()
            for (msg, _), bad in zip(checks, failed):
                if bad:
                    raise ValueError(msg)

        def _block_svdvals_min(self) -> torch.Tensor:
            """Smallest singular value over the active blocks, shape (batch,)."""
            B = self.matrix.shape[0]
            mins = [_fc.singular_values(x.detach())[..., -1].reshape(B, -1) for _, _, x in self._groups()]
            return torch.cat(mins, dim=-1).min(dim=-1)[0]

        @property
        def shape(self) -> Tuple[int, ...]:
            if not self._is_materialized and self._batch is not None:
                return (self._batch, self.algebra.C, self.algebra.k_max, self.algebra.k_max)
            return self.matrix.shape

        # ----- type checks -----
        @property
        def is_self_adjoint(self) -> bool:
            if self._is_self_adjoint is not None:
                return self._is_self_adjoint
            mat = self.matrix
            if self._is_self_adjoint is not None:          # set while materialising
                return self._is_self_adjoint
            self._is_self_adjoint = bool((mat - mat.conj().transpose(-2, -1)).abs().amax() <= self._tol_t()) \
                if mat.numel() else True
            return self._is_self_adjoint

        @is_self_adjoint.setter
        def is_self_adjoint(self, value: bool):
            self._is_self_adjoint = value
            self._is_self_adjoint_set = True

        @property
        def is_normal(self) -> bool:
            if self._is_normal is not None:
                return self._is_normal
            if self._is_self_adjoint:
                self._is_normal = True
                return True
            mat = self.matrix
            if self._is_normal is not None:                # set while materialising
                return self._is_normal
            if mat.numel() == 0:
                self._is_normal = True
                return True
            mat_d = mat.conj().transpose(-2, -1)
            diff = mat @ mat_d - mat_d @ mat
            self._is_normal = bool(diff.abs().amax() <= self._tol_t() * mat.abs().amax())
            return self._is_normal

        @is_normal.setter
        def is_normal(self, value: bool):
            self._is_normal = value
            self._is_normal_set = True

        @property
        def is_positive(self) -> bool:
            if self._is_positive is not None:
                return self._is_positive
            if not self.is_self_adjoint:
                self._is_positive = False
                return False
            self.matrix
            self._is_positive = bool((self.lambda_min >= -self._tol_t()).all())
            return self._is_positive

        @is_positive.setter
        def is_positive(self, value: bool):
            self._is_positive = value
            self._is_positive_set = True

        @property
        def is_invertible(self) -> bool:
            if self._is_invertible is not None:
                return self._is_invertible
            self.matrix
            self._is_invertible = bool((self._block_svdvals_min() > self._tol_t()).all())
            return self._is_invertible

        @is_invertible.setter
        def is_invertible(self, value: bool):
            self._is_invertible = value
            self._is_invertible_set = True

        @property
        def is_projection(self) -> bool:
            if self._is_projection is not None:
                return self._is_projection
            if not self.is_self_adjoint:
                self._is_projection = False
                return False
            mat = self.matrix
            self._is_projection = bool((mat @ mat - mat).abs().amax() <= self._tol_t()) if mat.numel() else True
            return self._is_projection

        @is_projection.setter
        def is_projection(self, value: bool):
            self._is_projection = value
            self._is_projection_set = True

        # ----- lambda_max and lambda_min (global extremes) -----
        @property
        def lambda_max(self) -> torch.Tensor:
            if self._lambda_max is not None:
                return self._lambda_max
            if not self.is_self_adjoint:
                raise RuntimeError("lambda_max requires self-adjoint operator")
            # exact for every block size: shifted power iteration was unreliable for symmetric spectra (+-l)
            # and for nearly degenerate edges; eigenvalues() skips the padding and forms no eigenvectors, and
            # uses the spectrum of operator_from_eigenvalues without any decomposition
            all_eigvals = torch.cat([w for w in self.eigenvalues() if w.numel()], dim=-1)   # (batch, sum k_c)
            self._lambda_max = all_eigvals.max(dim=1)[0]
            self._lambda_min = all_eigvals.min(dim=1)[0]
            return self._lambda_max

        @property
        def lambda_min(self) -> torch.Tensor:
            if self._lambda_min is not None:
                return self._lambda_min
            _ = self.lambda_max
            return self._lambda_min

        # ----- other methods (trace, norms, etc.) -----
        @property
        def trace(self) -> torch.Tensor:
            if self._trace is not None:
                return self._trace
            diag = torch.diagonal(self.matrix, dim1=-2, dim2=-1)
            self._trace = torch.sum(diag, dim=(-2, -1))
            return self._trace

        @trace.setter
        def trace(self, value: torch.Tensor):
            self._trace = value

        def trace_norm(self) -> torch.Tensor:
            """||A||_1 = sum of the singular values of the active blocks (for self-adjoint A: sum of |eigenvalues|,
            computed with the batched Hermitian eigensolver)."""
            B = self.matrix.shape[0]
            if self.is_self_adjoint:            # one check (cached) instead of an SVD
                vals = [cost.batched_call(torch.linalg.eigvalsh, x, f"eigvalsh (blocks of size {k})", kind='eigh').abs()
                        for k, _, x in self._groups()]
            else:
                vals = [_fc.singular_values(x) for k, _, x in self._groups()]
            total = torch.zeros(B, dtype=self.matrix.real.dtype, device=self.matrix.device)
            for v in vals:
                total = total + v.reshape(B, -1).sum(-1)
            return total

        def frobenius_norm(self) -> torch.Tensor:
            """(sum_c ||A_c||_F^2)^{1/2}; the padding is zero."""
            mat = self.matrix
            s = mat.abs().amax(dim=(-3, -2, -1)).clamp(min=torch.finfo(mat.real.dtype).tiny)      # no overflow
            return s * torch.linalg.vector_norm(mat / s[:, None, None, None].to(mat.dtype), dim=(-3, -2, -1))

        def _groups(self):
            """
            Sectors grouped by block size: [(k, idx, blocks)] with idx a list of sector indices of size k and
            blocks the active blocks of these sectors, shape (batch * len(idx), k, k). Decompositions run as
            one batched call per distinct block size (one call in total when all sectors have the same size).
            """
            mat = self.matrix
            alg = self.algebra
            out = []
            for k, idx in alg._size_groups():
                x = mat[:, alg._sel(idx), :k, :k]          # a view when the sectors are contiguous
                out.append((k, idx, x.reshape(-1, k, k)))
            return out

        def operator_norm(self) -> torch.Tensor:
            # largest singular value over the active blocks (exact; the padding is skipped)
            if self._is_self_adjoint:                      # from the extreme eigenvalues (cached, no SVD)
                return torch.maximum(self.lambda_max.abs(), self.lambda_min.abs())
            B = self.matrix.shape[0]
            norms = [_fc.singular_values(x)[..., 0].reshape(B, -1) for _, _, x in self._groups()]
            if not norms:
                return torch.zeros(B, dtype=self.matrix.real.dtype, device=self.matrix.device)
            return torch.cat(norms, dim=-1).max(dim=-1)[0]

        def _blockwise(self, f: Callable[[torch.Tensor], torch.Tensor], what: str) -> torch.Tensor:
            """f applied to the active block of every sector, one call per distinct block size (f acts on a
            stack (n, k, k)); the padding stays zero, so that decompositions and their gradients never see the
            zero singular values of the padding."""
            mat = self.matrix
            B = mat.shape[0]
            sel = self.algebra._sel
            groups = self._groups()
            if len(groups) == 1 and groups[0][0] == self.algebra.k_max and len(groups[0][1]) == self.algebra.C:
                return f(groups[0][2]).reshape(mat.shape)    # no padding: no scatter
            out = torch.zeros_like(mat)
            for k, idx, x in groups:
                out[:, sel(idx), :k, :k] = f(x).reshape(B, len(idx), k, k)
            return out

        def inverse(self, tol: Optional[float] = None) -> 'TypeIAlgebra.Operator':
            """(Pseudo-)inverse by SVD of every block; singular values below tol are treated as zero (default:
            relative tolerance k_c * eps * largest singular value of the block, as numpy.linalg.pinv per block)."""
            if self.is_self_adjoint:                        # Hermitian eigensolver instead of an SVD
                def f(w):                                   # w: eigenvalues of one block, (..., k_c)
                    thr = (w.shape[-1] * torch.finfo(w.dtype).eps * w.abs().amax(dim=-1, keepdim=True)
                           if tol is None else tol)
                    big = w.abs() > thr
                    return torch.where(big, 1.0 / torch.where(big, w, torch.ones_like(w)), torch.zeros_like(w))
                op = self.apply_function(f, per_block=True)
                op._is_self_adjoint = True
                op._is_normal = True
                op._is_positive = self._is_positive
                op._is_invertible = self._is_invertible
                return op

            def inv_block(A):
                return _fc.pinv(A, tol)
            op = TypeIAlgebra.Operator(self.algebra, generator=lambda: self._blockwise(inv_block, "inverse"))
            # carry over what is known (A^+ is normal / invertible with A); nothing is computed here
            op._is_normal = self._is_normal
            op._is_invertible = self._is_invertible
            op._batch = self._batch
            return op

        @property
        def inv(self) -> 'TypeIAlgebra.Operator':
            if self._inverse is None:
                self._inverse = self.inverse()
            return self._inverse

        def abs(self) -> 'TypeIAlgebra.Operator':
            """|A| = (A^* A)^{1/2}: for self-adjoint A by the functional calculus, otherwise from the singular
            value decomposition (on GPUs for blocks > 32 via the Hermitian eigensolver of [[0, A], [A^*, 0]]);
            gradients finite at degenerate singular values."""
            if self.is_self_adjoint:
                op = self.apply_function(torch.abs)
                op._is_self_adjoint = True
                op._is_positive = True
                op._is_normal = True
                return op

            def abs_block(A):
                return _fc.abs_(A)
            op = TypeIAlgebra.Operator(self.algebra, generator=lambda: self._blockwise(abs_block, "abs"))
            op._is_self_adjoint = True
            op._is_positive = True
            op._is_normal = True
            op._batch = self._batch
            return op

        def polar(self) -> Tuple['TypeIAlgebra.Operator', 'TypeIAlgebra.Operator']:
            """
            Polar decomposition A = W |A|, block by block: W unitary for invertible A (a partial isometry on the
            kernel of a singular A), |A| = (A^* A)^{1/2}. |A| is differentiable; W carries no gradient.
            """
            memo = {}                                       # one SVD per block size for both factors

            def factors():
                if not memo:
                    B, alg = self.matrix.shape[0], self.algebra
                    w_out, p_out = torch.zeros_like(self.matrix), torch.zeros_like(self.matrix)
                    for k, idx, x in self._groups():
                        W_, P_ = _fc.polar(x)
                        w_out[:, alg._sel(idx), :k, :k] = W_.reshape(B, len(idx), k, k)
                        p_out[:, alg._sel(idx), :k, :k] = P_.reshape(B, len(idx), k, k)
                    memo['W'], memo['P'] = w_out, p_out
                return memo
            W = TypeIAlgebra.Operator(self.algebra, generator=lambda: factors()['W'].detach())
            P = TypeIAlgebra.Operator(self.algebra, generator=lambda: factors()['P'])
            W._is_normal = None
            P._is_self_adjoint = P._is_positive = P._is_normal = True
            W._batch = P._batch = self._batch
            return W, P

        def _require_positive(self, what: str):
            """Positivity check that shares its decomposition with the functional calculus that follows: an
            unknown sign is decided from the eigendecomposition (computed once, kept for apply_function)."""
            if self._is_positive is True:
                return
            if self._is_positive is None and self.is_self_adjoint:
                dec = self._decompositions('hermitian', keep_once=True)
                wmin = torch.stack([wv[0].amin() for _, _, _, wv in dec]).amin()
                self._is_positive = bool(wmin >= -self._tol_t())
            if not self._is_positive:
                raise RuntimeError(what)

        def sqrt(self) -> 'TypeIAlgebra.Operator':
            """A^{1/2} of a positive operator (spectral calculus on the active blocks)."""
            self._require_positive("sqrt requires positive operator")
            op = self.apply_function(lambda w: torch.sqrt(torch.clamp(w, min=0)))
            op._is_self_adjoint = True
            op._is_positive = True
            op._is_normal = True
            return op

        def trace_a_log_a(self, tol: float = 1e-12) -> torch.Tensor:
            if self._is_positive is False or not self.is_self_adjoint:
                raise RuntimeError("trace_a_log_a requires positive operator")
            # eigenvalues of the active blocks only (the zero padding must not contribute), computed once and
            # also used for the positivity check if that is not known; values below tol count as exact zeros
            B = self.matrix.shape[0]
            out = torch.zeros(B, dtype=self.matrix.real.dtype, device=self.matrix.device)
            mins = []
            for k_c, _, x in self._groups():
                S = cost.batched_call(torch.linalg.eigvalsh, 0.5 * (x + x.conj().transpose(-2, -1)),
                                      f"eigvalsh (blocks of size {k_c})")
                mins.append(S.detach().amin())
                S = torch.where(S > tol, S, torch.zeros_like(S))
                out = out + torch.special.xlogy(S, S).reshape(B, -1).sum(-1)
            if self._is_positive is None:
                self._is_positive = bool(torch.stack(mins).amin() >= -self._tol_t())
            if not self._is_positive:
                raise RuntimeError("trace_a_log_a requires positive operator")
            return out

        def entropy(self) -> torch.Tensor:
            return -self.trace_a_log_a()

        @property
        def michelson_contrast(self) -> torch.Tensor:
            """
            Michelson contrast Delta(A) = (lambda_max - lambda_min) / (lambda_max + lambda_min)
            of a positive operator, with the extreme eigenvalues taken over the whole
            algebra (all channels). Returns a tensor of shape (batch,); Delta = 0 for A = 0.
            """
            if self._is_positive is not True and not self.is_positive:
                raise RuntimeError("Michelson contrast requires a positive operator.")
            lmax, lmin = self.lambda_max, self.lambda_min
            denom = lmax + lmin
            # A = 0 (up to rounding): Delta = 0; the cutoff is relative, Delta does not depend on the scale of A
            ok = denom > 4 * torch.finfo(denom.dtype).eps * (lmax.abs() + lmin.abs())
            safe = torch.where(ok, denom, torch.ones_like(denom))
            return torch.where(ok, (lmax - lmin) / safe, torch.zeros_like(denom))

        def _compatible(self, other, op: str) -> None:
            """Binary operations need operators of the same algebra (or of algebras with the same structure)."""
            if not isinstance(other, TypeIAlgebra.Operator):
                raise TypeError(f"'{op}' needs two operators; for a multiple of the identity use alg.identity() * x")
            a, b = self.algebra, other.algebra
            if a is not b and (a.k_factors != b.k_factors or a.n_factors != b.n_factors or a.charges != b.charges
                               or a.dtype != b.dtype or _same_device(a.hilbert.device, b.hilbert.device) is False):
                raise ValueError(f"'{op}' of operators from different algebras ({a} with {a.dtype} and {b} with "
                                 f"{b.dtype}); convert one with op.cast(alg)")

        def _derived(self, generator, other=None, **flags) -> 'TypeIAlgebra.Operator':
            """Lazy result with the property flags that are known to hold (not checked; None = unknown)."""
            op = TypeIAlgebra.Operator(self.algebra, generator=generator)
            for name, value in flags.items():
                setattr(op, f'_is_{name}', value)
            # known batch of the result: unknown (None) as soon as one operand's batch is unknown
            op._batch = self._batch if other is None else (
                None if self._batch is None or other._batch is None else max(self._batch, other._batch))
            return op

        def __add__(self, other: 'TypeIAlgebra.Operator') -> 'TypeIAlgebra.Operator':
            self._compatible(other, '+')
            sa = True if (self._is_self_adjoint and other._is_self_adjoint) else None
            pos = True if (self._is_positive and other._is_positive) else None
            return self._derived(lambda: self.matrix + other.matrix, other, self_adjoint=sa, positive=pos,
                                 normal=True if sa else None)

        def __matmul__(self, other: 'TypeIAlgebra.Operator') -> 'TypeIAlgebra.Operator':
            self._compatible(other, '@')
            inv = True if (self._is_invertible and other._is_invertible) else None
            return self._derived(lambda: torch.matmul(self.matrix, other.matrix), other, invertible=inv)

        def __mul__(self, scalar: Union[float, complex, torch.Tensor]) -> 'TypeIAlgebra.Operator':
            """Multiplication by a number, or by a tensor of shape (batch,) holding one number per sample."""
            if isinstance(scalar, torch.Tensor) and scalar.dim() > 0:
                if scalar.dim() != 1:
                    raise ValueError("multiply by a number or by a tensor of shape (batch,)")
                scalar = scalar.reshape(-1, 1, 1, 1)
            if not self.algebra.complex_valued and (
                    isinstance(scalar, complex) or (isinstance(scalar, torch.Tensor) and torch.is_complex(scalar))):
                raise ValueError("complex scalar in a real algebra")
            def scaled():
                m = self.matrix
                if isinstance(scalar, torch.Tensor):     # keep the algebra's dtype (no promotion to float64)
                    sc = scalar.to(m.device)
                    sc = sc.to(m.dtype if torch.is_complex(sc) else m.real.dtype)
                    return m * sc
                return m * scalar
            real = isinstance(scalar, (int, float)) or (isinstance(scalar, torch.Tensor) and not torch.is_complex(scalar))
            nonneg = isinstance(scalar, (int, float)) and scalar >= 0
            nonzero = isinstance(scalar, (int, float, complex)) and scalar != 0
            if isinstance(scalar, (int, float, complex)) and scalar == 0:     # 0 * A = 0, whatever A is
                return self._derived(scaled, normal=True, self_adjoint=True, positive=True, invertible=None)
            out = self._derived(scaled, normal=self._is_normal,
                                self_adjoint=self._is_self_adjoint if real else None,
                                positive=True if (self._is_positive and nonneg) else None,
                                invertible=self._is_invertible if nonzero else None)
            if isinstance(scalar, torch.Tensor) and scalar.dim() > 0 and self._batch is not None:
                out._batch = max(self._batch, scalar.shape[0])
            return out

        __rmul__ = __mul__

        def __sub__(self, other: 'TypeIAlgebra.Operator') -> 'TypeIAlgebra.Operator':
            self._compatible(other, '-')
            sa = True if (self._is_self_adjoint and other._is_self_adjoint) else None
            return self._derived(lambda: self.matrix - other.matrix, other, self_adjoint=sa, normal=True if sa else None)

        def __neg__(self) -> 'TypeIAlgebra.Operator':
            return self._derived(lambda: -self.matrix, self_adjoint=self._is_self_adjoint, normal=self._is_normal,
                                 invertible=self._is_invertible)

        def adjoint(self) -> 'TypeIAlgebra.Operator':
            return self._derived(lambda: self.matrix.conj().transpose(-2, -1), self_adjoint=self._is_self_adjoint,
                                 normal=self._is_normal, positive=self._is_positive,
                                 invertible=self._is_invertible, projection=self._is_projection)

        def cast(self, algebra: 'TypeIAlgebra') -> 'TypeIAlgebra.Operator':
            """
            This operator as an element of `algebra` (same sectors, other field / precision / device),
            e.g. alg.like(precision='double'). Casting a complex operator to a real algebra drops the
            imaginary part and raises if it is not negligible.
            """
            if list(algebra.k_factors) != list(self.algebra.k_factors):
                raise ValueError("cast: the sector structure must agree")
            mat = self.matrix
            if torch.is_complex(mat) and not algebra.complex_valued:
                if mat.imag.abs().max() > 1e-6 * max(1.0, mat.abs().max().item()):
                    raise ValueError("operator has a non-negligible imaginary part; the target algebra is real")
                mat = mat.real
            out = type(self).__new__(type(self))
            out.__dict__.update(self.__dict__)
            out.algebra = algebra
            out._matrix = mat.to(dtype=algebra.dtype, device=algebra.hilbert.device)
            out._is_materialized = True
            out._generator = None
            out._lambda_max = out._lambda_min = out._trace = out._inverse = out._abs = None
            out._eigenvalues = None
            out._spectral_cache = None
            out._update_memory()
            return out

        # ----- general functional calculus (spectral theorem, block by block) -----
        def eigh(self) -> List[Tuple[torch.Tensor, torch.Tensor]]:
            """
            Spectral decomposition of a self-adjoint operator, one entry per channel:
            (eigenvalues (batch, k_c), eigenvectors (batch, k_c, k_c)), A_c = V diag(w) V*.
            Channels with k_c = 0 give empty tensors.
            """
            if not self.is_self_adjoint:
                raise RuntimeError("eigh requires a self-adjoint operator")
            if self._caching() and not (self.matrix.requires_grad and torch.is_grad_enabled()):
                return self._from_cache(lambda wv: wv)
            return self._per_sector(lambda x: cost.batched_call(torch.linalg.eigh, x, f"eigh (blocks of size {x.shape[-1]})"))

        def _from_cache(self, pick):
            """Per-sector results pick((w, V)) from the cached Hermitian decomposition."""
            B = self.matrix.shape[0]
            res = [None] * self.algebra.C
            for k, idx, x, (w, V) in self._decompositions('hermitian'):
                y = pick((w, V))
                parts = y if isinstance(y, tuple) else (y,)
                parts = [p.reshape(B, len(idx), *p.shape[1:]) for p in parts]
                for j, c in enumerate(idx):
                    r = tuple(p[:, j] for p in parts)
                    res[c] = r if isinstance(y, tuple) else r[0]
            mat = self.matrix
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c == 0:
                    e = torch.linalg.eigh(mat[:, c, :0, :0])
                    res[c] = pick((e[0], e[1]))
            return res

        def _per_sector(self, fn):
            """fn on the stacked blocks of each size group, split back into one result per sector (tuples
            of tensors with leading dimensions (batch, ...)); empty sectors give empty tensors."""
            mat = self.matrix
            B = mat.shape[0]
            res = [None] * self.algebra.C
            for k, idx, x in self._groups():
                y = fn(x)
                parts = y if isinstance(y, tuple) else (y,)
                parts = [p.reshape(B, len(idx), *p.shape[1:]) for p in parts]
                for j, c in enumerate(idx):
                    r = tuple(p[:, j] for p in parts)
                    res[c] = r if isinstance(y, tuple) else r[0]
            for c, k_c in enumerate(self.algebra.k_factors):
                if k_c == 0:                                  # keep the old behaviour for empty sectors
                    res[c] = fn(mat[:, c, :0, :0])
            return res

        def eigenvalues(self) -> List[torch.Tensor]:
            """
            Eigenvalues of a self-adjoint operator, one real tensor (batch, k_c) per channel, ascending.
            Cheaper than eigh() (no eigenvectors are formed or stored); differentiable.
            Channels with k_c = 0 give empty tensors.
            """
            if not self.is_self_adjoint:
                raise RuntimeError("eigenvalues requires a self-adjoint operator")
            if self._eigenvalues is not None and self._batch is not None:   # spectrum known by construction
                return [torch.broadcast_to(e.real if torch.is_complex(e) else e, (self._batch, k)).sort(dim=-1).values
                        for e, k in zip(self._eigenvalues, self.algebra.k_factors)]
            if self._caching() and not (self.matrix.requires_grad and torch.is_grad_enabled()):
                return self._from_cache(lambda wv: wv[0])
            return self._per_sector(lambda x: cost.batched_call(torch.linalg.eigvalsh, x,
                                                                f"eigvalsh (blocks of size {x.shape[-1]})"))

        # ----- spectral cache (opt-in) -----
        def cache_spectrum(self, enable: bool = True) -> 'TypeIAlgebra.Operator':
            """
            Keep the eigendecomposition after its first use, so that further functions of this operator
            (apply_function, sqrt, log, power, eigenvalues, propagators for many times, ...) cost matrix products
            only. Costs one extra copy of the operator. Also switched on for all operators by
            `torch_vn_algebra.caching()`. Returns self. The cache is invalidated by in-place changes of the matrix
            (PyTorch version counter), but not by changes through `.data`.
            """
            self._cache_on = bool(enable)
            if not enable:
                self._spectral_cache = None
            return self

        def clear_cache(self):
            self._spectral_cache = None

        def cache_info(self) -> dict:
            c = self._spectral_cache
            nbytes = 0 if c is None else sum(t.element_size() * t.numel() for g in c['groups'] for t in g)
            return {'enabled': self._caching(), 'cached': c is not None and self._cache_valid(c),
                    'kind': None if c is None else c['kind'], 'bytes': nbytes}

        def _caching(self) -> bool:
            return self._cache_on or _fc_cache_enabled()

        def _cache_key(self):
            m = self.matrix
            return (id(m), m._version)

        def _cache_valid(self, c) -> bool:
            return c['key'] == self._cache_key()

        def _decompositions(self, kind: str, keep_once: bool = False):
            """[(k, idx, x, (w, V))] per block-size group, from the cache if it is valid (no gradients). keep_once:
            keep the result for the next call even without caching (a positivity check followed by the functional
            calculus decomposes once)."""
            groups = self._groups()
            c = self._spectral_cache
            if c is not None and c['kind'] == kind and self._cache_valid(c):
                if c.get('once') and not self._caching():
                    self._spectral_cache = None              # used: drop the one-off copy
                return [(k, idx, x, wv) for (k, idx, x), wv in zip(groups, c['groups'])]
            out = []
            with torch.no_grad():
                for k, idx, x in groups:
                    xd = x.detach()
                    wv = (_fc.hermitian_eig(xd, f"eigh (blocks of size {k})") if kind == 'hermitian'
                          else _fc.normal_eig(xd, f"eigh (blocks of size {k})"))
                    out.append((k, idx, x, wv))
            if self._caching() or keep_once:
                self._spectral_cache = {'kind': kind, 'key': self._cache_key(), 'groups': [o[3] for o in out],
                                        'once': not self._caching()}
            return out

        def apply_function(self, f: Callable[[torch.Tensor], torch.Tensor],
                           per_block: bool = False) -> 'TypeIAlgebra.Operator':
            """
            f(A) by the spectral theorem, f(A)_c = V_c f(w_c) V_c*, for self-adjoint A (f on the real spectrum)
            and for normal A such as unitaries (f on the complex spectrum; f must accept complex tensors).
            `f` must act elementwise: it is called once per block size with the eigenvalues of all n sectors of
            that size, a tensor (batch, n * k_c), and returns a tensor of the same shape (per-sample parameters
            of shape (batch, 1) broadcast; the operator must then have the full batch, e.g. (A * t) instead of a
            batch-1 A). With per_block=True f instead receives one block per row, (batch * n, k_c), and may
            couple the eigenvalues of a block (e.g. thresholds relative to the block); per-sample parameters do
            not broadcast then. A complex-valued result requires a complex algebra. Only the Hermitian
            eigensolver is used (for normal A, on H + alpha K with A = H + iK), one batched call per block size.
            Gradients with respect to A use the Daleckii-Krein formula (divided differences of f): finite also
            at degenerate spectra; for normal non-self-adjoint A they assume f holomorphic near the spectrum.
            Gradients with respect to parameters of f (e.g. t in exp(-itw)) are propagated as usual. For
            non-normal A use holomorphic_function.
            """
            alg = self.algebra

            def generator():
                if self.is_self_adjoint:
                    kind = 'hermitian'
                elif self.is_normal:
                    kind = 'normal'
                else:
                    raise RuntimeError("apply_function requires a self-adjoint or normal operator "
                                       "(use holomorphic_function for other operators)")
                B = self.matrix.shape[0]
                parts = []
                for k_c, idx, x, wv in self._decompositions(kind):
                    n = len(idx)

                    def f_values(w, n=n, k_c=k_c):
                        # f sees the eigenvalues as (batch, n * k_c) (one call for all sectors of this size; an
                        # elementwise f with per-sample parameters of shape (batch, 1) broadcasts as for (batch, k_c))
                        if per_block:
                            return f(w)
                        return f(w.reshape(B, n * k_c)).reshape(-1, k_c)
                    y = _fc.spectral_apply(x, f_values, kind, decomposition=wv)
                    parts.append((k_c, idx, y))
                cplx = any(torch.is_complex(y) for _, _, y in parts)
                out_dtype = self.matrix.dtype
                if cplx and not torch.is_complex(self.matrix):
                    big = torch.stack([y.detach().imag.abs().amax() for _, _, y in parts if torch.is_complex(y)]).amax()
                    if bool(big > self._tol_t()):                # one synchronisation (validation)
                        raise ValueError("complex-valued f(A) requires a complex-valued algebra")
                out = torch.zeros(self.matrix.shape, dtype=out_dtype, device=self.matrix.device)
                for k_c, idx, y in parts:
                    y = y.real if torch.is_complex(y) and not torch.is_complex(out) else y
                    out[:, alg._sel(idx), :k_c, :k_c] = y.to(out_dtype).reshape(B, len(idx), k_c, k_c)
                result._is_normal = True                       # f(A) of a normal A is normal
                if kind == 'hermitian' and not cplx:
                    result._is_self_adjoint = True
                return out
            result = self._derived(generator)
            return result

        def holomorphic_function(self, f: Callable[[torch.Tensor], torch.Tensor], center=None, radius=None,
                                 nodes: int = 32, max_nodes: int = 2048,
                                 rtol: Optional[float] = None, atol: float = 0.0) -> 'TypeIAlgebra.Operator':
            """
            f(A) for any operator A and f holomorphic near the spectrum (holomorphic functional calculus),
            f(A)_c = (2 pi i)^{-1} oint f(z) (z - A_c)^{-1} dz on a circle, by the trapezoidal rule with batched
            LU solves; the number of nodes is doubled until the result changes by less than rtol (default 1e-10 in
            double, 1e-5 in single precision) relative to its size; ValueError if the estimate stays above
            sqrt(rtol). Default circle: just around the spectrum of every block (functional.default_contour);
            f must be holomorphic on the closed disk and accept complex tensors. For functions with singularities
            (log, powers, resolvents) pass center and radius: numbers, or callables mapping a stack of blocks
            (n, k, k) to tensors (n,). The estimated error is stored as result.contour_error after evaluation.
            The error is relative to max(|f(X_c)|, atol): when f(X) may vanish (sin at X = 0, z^3 at a
            nilpotent block) pass atol > 0, a scale for f near the spectrum; the absolute accuracy accepted where
            |f(X_c)| < atol is rtol * atol. Differentiable through the solves. For self-adjoint or normal A, apply_function is exact and cheaper;
            for exp, expm is more accurate when the spectrum is spread.
            """
            alg = self.algebra
            result = None

            def generator():
                out = torch.zeros_like(self.matrix)
                B = out.shape[0]
                error = 0.0
                for k_c, idx, x in self._groups():
                    y, err = _fc.holomorphic(x, f, center, radius, nodes, max_nodes, rtol, atol=atol)
                    error = max(error, err)
                    out[:, alg._sel(idx), :k_c, :k_c] = y.to(out.dtype).reshape(B, len(idx), k_c, k_c)
                result.contour_error = error
                return out
            result = self._derived(generator)
            return result

        def expm(self, scale: Union[float, complex] = 1.0) -> 'TypeIAlgebra.Operator':
            """exp(scale * A) for any operator (torch.linalg.matrix_exp on each active block)."""
            alg = self.algebra
            if isinstance(scale, complex) and not alg.hilbert.complex_valued:
                raise ValueError("complex scale requires a complex-valued algebra")
            if isinstance(scale, torch.Tensor) and scale.dim() > 0:
                raise ValueError("expm: scale must be a number (for one value per sample use (A * t).expm())")

            def generator():
                return self._blockwise(lambda x: cost.batched_call(
                    torch.linalg.matrix_exp, scale * x, f"expm (blocks of size {x.shape[-1]})", kind='expm'), 'expm')
            return self._derived(generator, invertible=True)

        def log(self, eps: float = 1e-12) -> 'TypeIAlgebra.Operator':
            """log(A) of a positive operator; eigenvalues are clamped to eps."""
            self._require_positive("log requires a positive operator")
            op = self.apply_function(lambda w: torch.log(torch.clamp(w, min=eps)))
            op._is_self_adjoint = True
            op._is_normal = True
            return op

        def power(self, p: float) -> 'TypeIAlgebra.Operator':
            """A^p of a positive operator (negative p: on the support only)."""
            self._require_positive("power requires a positive operator")
            def f(w):                                       # w: eigenvalues of one block, (..., k_c)
                w = torch.clamp(w, min=0.0)
                if p >= 0:
                    return w ** p
                # support: eigenvalues above the rounding level of the block (as for inverse / pinv)
                thr = w.shape[-1] * torch.finfo(w.dtype).eps * w.amax(dim=-1, keepdim=True)
                big = w > thr
                return torch.where(big, torch.where(big, w, torch.ones_like(w)) ** p, torch.zeros_like(w))
            op = self.apply_function(f, per_block=p < 0)
            op._is_self_adjoint = True
            op._is_positive = True
            return op

        def _weighted_trace(self, name: str) -> torch.Tensor:
            """sum_c w_c Tr(A_c), in the number field of the algebra (complex for complex algebras); the weights are
            cached on the algebra."""
            tr_c = torch.diagonal(self.matrix, dim1=-2, dim2=-1).sum(-1)          # (batch, C), padding = 0
            return (tr_c * self.algebra._const(name, self.matrix.real.dtype)).sum(-1)

        def Tr_blunt(self) -> torch.Tensor:
            """sum_c Tr(A_c)."""
            return torch.diagonal(self.matrix, dim1=-2, dim2=-1).sum((-2, -1))

        def Tr_norm(self) -> torch.Tensor:
            """sum_c Tr(A_c) / k_c (normalised trace of every factor, summed)."""
            return self._weighted_trace('w_norm')

        def tau_vN(self) -> torch.Tensor:
            """(1/C) sum_c Tr(A_c) / k_c over the non-empty sectors: the faithful normal tracial state with equal
            sector weights."""
            return self._weighted_trace('w_vn')

        def __repr__(self) -> str:
            if not self._is_materialized:   # do not trigger materialisation
                return f"Operator(C={self.algebra.C}, k_max={self.algebra.k_max}, lazy)"
            mem = self._bytes_held / 1024
            return f"Operator({tuple(self._matrix.shape)}, materialized, {mem:.1f}KB)"

    # ========================================================================
    # Factory Methods
    # ========================================================================

    def operator_from_eigenvalues(
        self,
        eigenvalue_sampler: Callable[[int], torch.Tensor],
        batch_size: int = 1,
        unitary_measure: str = 'haar',
        force_self_adjoint: Optional[bool] = None,
        force_positive: Optional[bool] = None,
        force_normal: Optional[bool] = None,
        force_invertible: Optional[bool] = None,
        force_projection: Optional[bool] = None,
    ) -> Operator:
        cost.check_memory(4 * cost.tensor_bytes((batch_size, self.C, self.k_max, self.k_max), self.hilbert.dtype),
                          self.hilbert.device, f"operator_from_eigenvalues(batch={batch_size})")
        if unitary_measure not in self.UNITARY_MEASURES:
            raise ValueError(f"Unknown measure: {unitary_measure}. Choose from {self.UNITARY_MEASURES}")
        all_eig = []
        flags = []                                   # per sector: non-real, negative, zero, not 0/1
        for c, k_c in enumerate(self.k_factors):
            if k_c == 0:
                all_eig.append(torch.zeros(0, dtype=self.hilbert.dtype, device=self.hilbert.device))
                continue
            eig = eigenvalue_sampler(k_c)
            eig = _as_tensor(eig)
            if eig.shape[-1:] != (k_c,) or eig.dim() > 2 or (eig.dim() == 2 and eig.shape[0] not in (1, batch_size)):
                raise ValueError(f"eigenvalue_sampler({k_c}) must return shape ({k_c},) or ({batch_size}, {k_c}), "
                                 f"got {tuple(eig.shape)}")
            # the flags are computed where the samples are (no device round trip for samplers on the CPU)
            re, im = (eig.real, eig.imag.abs()) if torch.is_complex(eig) else (eig, torch.zeros_like(eig))
            flags.append(torch.stack([torch.any(im > 1e-12), torch.any(re < -1e-12),
                                      torch.any(torch.abs(eig) < 1e-12),
                                      ~torch.all((torch.abs(re) < 1e-12) | (torch.abs(re - 1) < 1e-12))]))
            all_eig.append(eig)
        flags = torch.stack([f.to(flags[0].device) for f in flags]).tolist()   # one host synchronisation at most
        it = iter(flags)
        flags = [next(it) if k_c else [False, False, False, False] for k_c in self.k_factors]
        for c, k_c in enumerate(self.k_factors):
            if k_c and torch.is_complex(all_eig[c]) and not self.hilbert.complex_valued and flags[c][0]:
                raise ValueError("eigenvalues: complex values in a real algebra (create it with complex_valued=True)")
        all_eig = [e if not k_c else self._to_field(e, "eigenvalues") for e, k_c in zip(all_eig, self.k_factors)]
        nonreal = [f[0] for f in flags]
        # U diag(eig) U^* is always normal; a non-real eigenvalue rules out self-adjoint, positive and projection
        inferred_normal = True
        inferred_self_adjoint = not any(nonreal)
        inferred_positive = not any(f[0] or f[1] for f in flags)
        inferred_invertible = not any(f[2] for f in flags)
        inferred_projection = not any(f[0] or f[3] for f in flags)

        if force_self_adjoint is not None:
            if force_self_adjoint and not inferred_self_adjoint:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig) and torch.any(eig.imag.abs() > 1e-12):
                        raise ValueError(f"force_self_adjoint=True but eigenvalues for channel {c} have non-zero imaginary part")
            inferred_self_adjoint = force_self_adjoint

        if force_positive is not None:
            if force_positive and not inferred_positive:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig):
                        if torch.any(eig.imag.abs() > 1e-12):
                            raise ValueError(f"force_positive=True but eigenvalues for channel {c} are not real")
                        if torch.any(eig.real < -1e-12):
                            bad_indices = torch.where(eig.real < -1e-12)[0].tolist()
                            raise ValueError(
                                f"force_positive=True but eigenvalues for channel {c} are negative (real part): "
                                f"indices {bad_indices}"
                            )
                    else:
                        if torch.any(eig < -1e-12):
                            bad_indices = torch.where(eig < -1e-12)[0].tolist()
                            raise ValueError(
                                f"force_positive=True but eigenvalues for channel {c} are negative: "
                                f"indices {bad_indices}"
                            )
            inferred_positive = force_positive

        if force_normal is not None:
            inferred_normal = force_normal

        if force_invertible is not None:
            if force_invertible and not inferred_invertible:
                for c, eig in enumerate(all_eig):
                    if torch.any(torch.abs(eig) < 1e-12):
                        raise ValueError(f"force_invertible=True but eigenvalues for channel {c} contain zeros")
            inferred_invertible = force_invertible

        if force_projection is not None:
            if force_projection and not inferred_projection:
                for c, eig in enumerate(all_eig):
                    if torch.is_complex(eig):
                        bad_mask = ((torch.abs(eig.real) > 1e-12) & (torch.abs(eig.real - 1) > 1e-12)) \
                            | (eig.imag.abs() > 1e-12)
                    else:
                        bad_mask = (torch.abs(eig) > 1e-12) & (torch.abs(eig - 1) > 1e-12)
                    if torch.any(bad_mask):
                        raise ValueError(f"force_projection=True but eigenvalues for channel {c} are not 0 or 1")
            inferred_projection = force_projection

        def generator():
            # one random_unitary call per distinct block size, for all sectors of that size at once
            mat = torch.zeros(batch_size, self.C, self.k_max, self.k_max, dtype=self.hilbert.dtype,
                              device=self.hilbert.device)
            for k, idx in self._size_groups():
                eig = torch.stack([torch.broadcast_to(all_eig[c], (batch_size, k)) for c in idx], dim=1)
                U = self.random_unitary(k, measure=unitary_measure, batch_size=batch_size * len(idx))
                U = U.reshape(batch_size, len(idx), k, k)
                A = (U * eig.unsqueeze(-2)) @ U.mH
                # real spectrum: remove round-off so that the block is exactly Hermitian (decided on the host)
                real = [not nonreal[c] for c in idx]
                if all(real):
                    A = 0.5 * (A + A.mH)
                elif any(real):
                    m = torch.tensor(real, device=A.device)
                    A = torch.where(m[:, None, None], 0.5 * (A + A.mH), A)
                if k == self.k_max and len(idx) == self.C:
                    return A                                 # one group without padding: no copy
                mat[:, self._sel(idx), :k, :k] = A
            return mat

        op = self.Operator(self, generator=generator)
        op._batch = batch_size
        op._eigenvalues = all_eig
        op._is_self_adjoint = inferred_self_adjoint
        op._is_positive = inferred_positive
        op._is_normal = inferred_normal
        op._is_invertible = inferred_invertible
        op._is_projection = inferred_projection
        return op

    def _to_field(self, t: torch.Tensor, what: str) -> torch.Tensor:
        """Convert input data to this algebra's dtype and device; complex data with a non-zero imaginary
        part in a real algebra is an error rather than a silent truncation."""
        t = _as_tensor(t)
        if torch.is_complex(t) and not self.hilbert.complex_valued:
            if t.numel() and bool(t.imag.abs().amax() > 1e-6 * torch.clamp(t.abs().amax(), min=1)):
                raise ValueError(f"{what}: complex values in a real algebra (create it with complex_valued=True)")
            t = t.real
        return t.to(device=self.hilbert.device, dtype=self.hilbert.dtype)

    def operator(self, matrix: torch.Tensor, validate: bool = False, **properties) -> Operator:
        """
        Wrap a tensor of shape (batch, C, k_max, k_max) as an operator of this algebra.
        Keyword arguments (is_self_adjoint, is_positive, ...) tag known properties. They are trusted (no
        synchronisation, no decomposition) unless validate=True, which checks them once; a wrong flag selects a
        wrong algorithm (e.g. a non-Hermitian matrix tagged is_self_adjoint is treated as its Hermitian part).
        """
        expected = (self.C, self.k_max, self.k_max)
        if tuple(matrix.shape[1:]) != expected:
            raise ValueError(f"expected shape (batch, {self.C}, {self.k_max}, {self.k_max}), got {tuple(matrix.shape)}")
        mat = self._to_field(matrix, "matrix")
        if any(k_c < self.k_max for k_c in self.k_factors):     # the padding must be zero (one mask product)
            mat = mat * self._const('block_mask', torch.bool)
        op = self.Operator(self, matrix=mat, **properties)
        if validate:
            op._validate_properties()
        return op

    def from_blocks(self, blocks: List[torch.Tensor]) -> Operator:
        """Operator from a list of C blocks of shape (batch, k_c, k_c) (or (k_c, k_c))."""
        if len(blocks) != self.C:
            raise ValueError(f"from_blocks: need {self.C} blocks, got {len(blocks)}")
        blocks = [_as_tensor(b) for b in blocks]
        blocks = [b if b.dim() == 3 else b.unsqueeze(0) for b in blocks]
        batch = max(b.shape[0] for b in blocks)
        mat = torch.zeros(batch, self.C, self.k_max, self.k_max, dtype=self.hilbert.dtype, device=self.hilbert.device)
        for c, (b, k_c) in enumerate(zip(blocks, self.k_factors)):
            if tuple(b.shape[-2:]) != (k_c, k_c):
                raise ValueError(f"block {c} has shape {tuple(b.shape[-2:])}, expected ({k_c}, {k_c})")
            mat[:, c, :k_c, :k_c] = self._to_field(b, f"block {c}")
        return self.Operator(self, matrix=mat)

    def random_unitary_operator(self, batch_size: int = 1, measure: str = 'haar') -> Operator:
        """
        Block-diagonal unitary U = (+)_c U_c in the algebra, with independent
        U_c drawn from `measure` (see random_unitary) for every channel and sample.
        """
        if measure not in self.UNITARY_MEASURES:
            raise ValueError(f"Unknown measure: {measure}. Choose from {self.UNITARY_MEASURES}")

        def generator():
            mat = torch.zeros(batch_size, self.C, self.k_max, self.k_max,
                              dtype=self.hilbert.dtype, device=self.hilbert.device)
            for k, idx in self._size_groups():          # one call per distinct block size
                U = self.random_unitary(k, measure=measure, batch_size=batch_size * len(idx))
                if k == self.k_max and len(idx) == self.C:
                    return U.reshape(batch_size, self.C, k, k)
                mat[:, self._sel(idx), :k, :k] = U.reshape(batch_size, len(idx), k, k)
            return mat
        op = self.Operator(self, generator=generator)
        op._batch = batch_size
        op._is_normal = True
        op._is_invertible = True
        return op

    def identity(self, batch_size: int = 1) -> Operator:
        """The unit of the algebra. Its matrix is a read-only expanded view of a cached (C, k_max, k_max) tensor
        (no memory per sample); clone it before modifying it in place."""
        def generator():
            return self._const('eye').unsqueeze(0).expand(batch_size, -1, -1, -1)
        op = self.Operator(self, generator=generator)
        op._batch = batch_size
        op._is_self_adjoint = True
        op._is_normal = True
        op._is_positive = True
        op._is_invertible = True
        op._is_projection = True
        return op

    def zero(self, batch_size: int = 1) -> Operator:
        shape = (batch_size, self.C, self.k_max, self.k_max)
        def generator():
            return torch.zeros(shape, dtype=self.hilbert.dtype, device=self.hilbert.device)
        op = self.Operator(self, generator=generator)
        op._batch = batch_size
        op._is_self_adjoint = True
        op._is_normal = True
        op._is_positive = True
        op._is_invertible = False
        op._is_projection = True
        return op

    def _size_groups(self) -> List[Tuple[int, List[int]]]:
        """Non-empty sectors grouped by block size, in order of first appearance: [(k, [c, ...]), ...]."""
        return self._size_groups_cache

    def central(self, scalars: List[Union[float, complex, torch.Tensor]], batch_size: int = 1) -> Operator:
        """
        Central element sum_c s_c 1_c: one number per sector, or a tensor (batch,) per sector for one value per
        sample. Tensors keep their gradients; complex values need a complex algebra.
        """
        if len(scalars) != self.C:
            raise ValueError(f"central: need one number per sector ({self.C}), got {len(scalars)}")
        # Python numbers in the algebra's precision (torch.as_tensor(1/3) would be float32)
        ts = [x if isinstance(x, torch.Tensor) else
              torch.tensor(x, dtype=(torch.complex128 if self.hilbert.real_dtype == torch.float64 else torch.complex64)
                           if isinstance(x, complex) else self.hilbert.real_dtype)
              for x in scalars]
        if any(t.dim() > 1 for t in ts):
            raise ValueError("central: each value must be a number or a tensor of shape (batch,)")
        B = max([batch_size] + [t.shape[0] for t in ts if t.dim() == 1])
        dt = self.hilbert.dtype
        for t in ts:
            dt = torch.promote_types(dt, t.dtype)
        V = torch.stack([t.to(device=self.hilbert.device, dtype=dt).reshape(-1).expand(B) for t in ts], dim=1)
        V = self._to_field(V, "central values")                              # (B, C)

        def generator():
            return V[:, :, None, None] * self._const('eye')
        op = self.Operator(self, generator=generator)
        op._batch = B
        act = V.detach()[:, [c for c, k_c in enumerate(self.k_factors) if k_c]]
        re = act.real if torch.is_complex(act) else act
        imag0 = (act.imag == 0).all() if torch.is_complex(act) else torch.ones((), dtype=torch.bool, device=act.device)
        real, pos, inv, proj = torch.stack([imag0, (re >= 0).all(), (act != 0).all(),
                                            ((re == 0) | (re == 1)).all()]).tolist()   # one synchronisation
        op._is_self_adjoint = real
        op._is_normal = True                               # every central element is normal
        op._is_positive = real and pos
        op._is_invertible = inv
        op._is_projection = real and proj
        return op

    def __repr__(self) -> str:
        return f"TypeIAlgebra(C={self.C}, n={self.n_factors}, k={self.k_factors})"