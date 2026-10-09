import torch
from typing import Optional, Tuple, List


class HilbertSpace:
    """
    Hilbert space H of dimension n with a distinguished k-dimensional subspace H0.
    
    The space has a direct sum decomposition over n_channels factors (central decomposition
    of the von Neumann algebra). Vectors have shape (batch, n_channels, n, 1).
    Norm is defined as sqrt( sum_{c=1}^{n_channels} sum_{i=1}^{n} |ψ_{c,i}|^2 ).
    
    Attributes:
        n: Ambient dimension (size of each factor space)
        k: Subspace dimension (support of the algebra, same for all channels)
        batch_size: Batch dimension
        n_channels: Number of factors in the direct sum
        complex_valued: bool
        device: torch.device
        dtype: torch.dtype
    """
    
    def __init__(
        self,
        n: int,
        k: int,
        batch_size: int = 1,
        n_channels: int = 1,
        complex_valued: bool = True,
        device: Optional[torch.device] = None
    ):
        assert k <= n, f"k={k} > n={n}"
        assert k >= 0
        assert n > 0
        assert batch_size > 0
        assert n_channels > 0
        
        self.n = n
        self.k = k
        self.batch_size = batch_size
        self.n_channels = n_channels
        self.complex_valued = complex_valued
        
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device
        
        self.dtype = torch.complex64 if complex_valued else torch.float32
    
    @property
    def dim(self) -> int:
        return self.n
    
    @property
    def subspace_dim(self) -> int:
        return self.k
    
    @property
    def shape_ambient(self) -> Tuple[int, int, int, int]:
        return (self.batch_size, self.n_channels, self.n, self.n)
    
    @property
    def shape_ket(self) -> Tuple[int, int, int, int]:
        return (self.batch_size, self.n_channels, self.n, 1)
    
    @property
    def shape_bra(self) -> Tuple[int, int, int, int]:
        return (self.batch_size, self.n_channels, 1, self.n)
    
    def inner_product(self, bra: torch.Tensor, ket: torch.Tensor) -> torch.Tensor:
        """
        Global inner product ⟨ψ|φ⟩, summing over channels and spatial dimensions.
        Returns tensor of shape (batch, 1, 1, 1).
        """
        # Input shapes: (batch, n_channels, 1, n) and (batch, n_channels, n, 1)
        # or flattened variants.
        if bra.dim() == 4 and bra.shape[-2] == 1:
            bra_flat = bra.squeeze(-2)  # (batch, n_channels, n)
        else:
            bra_flat = bra
        
        if ket.dim() == 4 and ket.shape[-1] == 1:
            ket_flat = ket.squeeze(-1)  # (batch, n_channels, n)
        else:
            ket_flat = ket
        
        # Sum over n (last dimension) and over channels (dim=1)
        if self.complex_valued:
            prod = bra_flat.conj() * ket_flat                 # (batch, n_channels, n)
            inner = prod.sum(dim=(-1, -2), keepdim=True)      # (batch, 1, 1)
        else:
            prod = bra_flat * ket_flat
            inner = prod.sum(dim=(-1, -2), keepdim=True)
        
        # Expand to (batch, 1, 1, 1) for consistency
        return inner.unsqueeze(-1)   # (batch, 1, 1, 1)
    
    def norm(self, ket: torch.Tensor) -> torch.Tensor:
        """
        Global norm: sqrt( sum_c sum_i |ψ_{c,i}|^2 ).
        Returns (batch, 1, 1, 1).
        """
        # ket shape: (batch, n_channels, n, 1)
        norm_sq = (ket.conj() * ket).real.sum(dim=(-2, -3), keepdim=True)  # (batch, 1, 1, 1)
        norm_sq = torch.clamp(norm_sq, min=0.0)
        return torch.sqrt(norm_sq)
    
    def normalize(self, ket: torch.Tensor, eps: float = 1e-12) -> torch.Tensor:
        """Normalize ket to unit global norm."""
        norm = self.norm(ket)
        # Replace zero norms with 1 to avoid division by zero
        norm = torch.where(norm < eps, torch.ones_like(norm), norm)
        return ket / norm
    
    class Basis:
        """
        Orthonormal basis for the k-dimensional subspace H0.
        
        The basis matrix V has shape (batch, n_channels, n, k). Each channel may have
        its own orthonormal basis (columns orthonormal per channel).
        """
        
        def __init__(
            self,
            parent: 'HilbertSpace',
            V: Optional[torch.Tensor] = None,
            initialization: str = 'random'
        ):
            self.parent = parent
            self.n = parent.n
            self.k = parent.k
            self.batch_size = parent.batch_size
            self.n_channels = parent.n_channels
            self.device = parent.device
            self.dtype = parent.dtype
            self.complex_valued = parent.complex_valued
            
            if V is not None:
                self.V = V.to(device=self.device, dtype=self.dtype)
                expected = (self.batch_size, self.n_channels, self.n, self.k)
                assert self.V.shape == expected, f"Expected {expected}, got {self.V.shape}"
            else:
                self.V = self._initialize(initialization)
        
        def _initialize(self, method: str) -> torch.Tensor:
            if method == 'standard':
                V = torch.zeros(self.batch_size, self.n_channels, self.n, self.k,
                                dtype=self.dtype, device=self.device)
                for i in range(min(self.k, self.n)):
                    V[..., i, i] = 1.0
                return V
            
            elif method == 'random':
                random_mat = torch.randn(self.batch_size, self.n_channels, self.n, self.k,
                                         device=self.device)
                if self.complex_valued:
                    random_imag = torch.randn(self.batch_size, self.n_channels, self.n, self.k,
                                              device=self.device)
                    random_mat = random_mat + 1j * random_imag
                random_mat = random_mat.to(self.dtype)
                
                # QR per batch and channel
                V_list = []
                for b in range(self.batch_size):
                    V_channel_list = []
                    for c in range(self.n_channels):
                        q, r = torch.linalg.qr(random_mat[b, c], mode='reduced')
                        signs = torch.diag(torch.sign(torch.diag(r).real))
                        q = torch.matmul(q, signs.to(self.dtype))
                        V_channel_list.append(q)
                    V_list.append(torch.stack(V_channel_list, dim=0))
                return torch.stack(V_list, dim=0)
            
            elif method == 'haar':
                full_mat = torch.randn(self.batch_size, self.n_channels, self.n, self.n,
                                       device=self.device)
                if self.complex_valued:
                    full_imag = torch.randn(self.batch_size, self.n_channels, self.n, self.n,
                                            device=self.device)
                    full_mat = full_mat + 1j * full_imag
                full_mat = full_mat.to(self.dtype)
                
                V_list = []
                for b in range(self.batch_size):
                    V_channel_list = []
                    for c in range(self.n_channels):
                        q, _ = torch.linalg.qr(full_mat[b, c])
                        V_channel_list.append(q[:, :self.k])
                    V_list.append(torch.stack(V_channel_list, dim=0))
                return torch.stack(V_list, dim=0)
            else:
                raise ValueError(f"Unknown method: {method}")
        
        @property
        def projection(self) -> torch.Tensor:
            V_conj_T = torch.conj(torch.transpose(self.V, -2, -1))
            return torch.matmul(self.V, V_conj_T)
        
        @property
        def restriction(self) -> torch.Tensor:
            return torch.conj(torch.transpose(self.V, -2, -1))
        
        def embed_vector(self, v_sub: torch.Tensor) -> torch.Tensor:
            """Embed vector from subspace coordinates to ambient space."""
            return torch.matmul(self.V, v_sub)
        
        def restrict_vector(self, v_amb: torch.Tensor) -> torch.Tensor:
            V_conj_T = torch.conj(torch.transpose(self.V, -2, -1))
            return torch.matmul(V_conj_T, v_amb)
        
        def embed_operator(self, A_sub: torch.Tensor) -> torch.Tensor:
            tmp = torch.matmul(self.V, A_sub)
            V_conj_T = torch.conj(torch.transpose(self.V, -2, -1))
            return torch.matmul(tmp, V_conj_T)
        
        def restrict_operator(self, A_amb: torch.Tensor) -> torch.Tensor:
            V_conj_T = torch.conj(torch.transpose(self.V, -2, -1))
            tmp = torch.matmul(V_conj_T, A_amb)
            return torch.matmul(tmp, self.V)
        
        def ket(self, idx: int) -> torch.Tensor:
            assert 0 <= idx < self.k
            return self.V[..., idx].unsqueeze(-1)
        
        def bra(self, idx: int) -> torch.Tensor:
            ket = self.ket(idx)
            return torch.conj(torch.transpose(ket, -2, -1))
        
        def outer_product(self, i: int, j: int) -> torch.Tensor:
            return torch.matmul(self.ket(i), self.bra(j))
        
        def gram_matrix(self) -> torch.Tensor:
            V_conj_T = torch.conj(torch.transpose(self.V, -2, -1))
            return torch.matmul(V_conj_T, self.V)
        
        def is_orthonormal(self, tol: float = 1e-6) -> bool:
            gram = self.gram_matrix()
            identity = torch.eye(self.k, dtype=self.dtype, device=self.device)
            identity_exp = identity.unsqueeze(0).unsqueeze(0).expand(
                self.batch_size, self.n_channels, self.k, self.k
            )
            return torch.allclose(gram, identity_exp, atol=tol)
        
        def orthonormalize(self):
            V_new_list = []
            for b in range(self.batch_size):
                V_channel_list = []
                for c in range(self.n_channels):
                    q, r = torch.linalg.qr(self.V[b, c], mode='reduced')
                    signs = torch.diag(torch.sign(torch.diag(r).real))
                    q = torch.matmul(q, signs.to(self.dtype))
                    V_channel_list.append(q)
                V_new_list.append(torch.stack(V_channel_list, dim=0))
            with torch.no_grad():
                self.V.copy_(torch.stack(V_new_list, dim=0))
        
        def random_subspace_vector(self, normalize: bool = True) -> torch.Tensor:
            """
            Generate a random vector in the subspace (in ambient coordinates).
            If normalize=True, the global norm (sum over channels) becomes 1.
            """
            if self.complex_valued:
                real = torch.randn(self.batch_size, self.n_channels, self.k, 1, device=self.device)
                imag = torch.randn(self.batch_size, self.n_channels, self.k, 1, device=self.device)
                v_sub = torch.complex(real, imag)
            else:
                v_sub = torch.randn(self.batch_size, self.n_channels, self.k, 1, device=self.device)
            
            if normalize:
                # Compute global norm: sum over channels (dim=-3) and over k (dim=-2)
                global_norm_sq = (v_sub.conj() * v_sub).real.sum(dim=(-2, -3), keepdim=True)
                global_norm_sq = torch.clamp(global_norm_sq, min=0.0)
                global_norm = torch.sqrt(global_norm_sq)  # (batch, 1, 1, 1)
                global_norm = torch.where(global_norm < 1e-12, torch.ones_like(global_norm), global_norm)
                v_sub = v_sub / global_norm
            
            return self.embed_vector(v_sub)
        
        def __repr__(self) -> str:
            return f"Basis(n={self.n}, k={self.k}, batch={self.batch_size}, channels={self.n_channels})"
    
    def standard_basis(self) -> Basis:
        return self.Basis(self, initialization='standard')
    
    def random_basis(self) -> Basis:
        return self.Basis(self, initialization='random')
    
    def haar_basis(self) -> Basis:
        return self.Basis(self, initialization='haar')
    
    def __repr__(self) -> str:
        return f"HilbertSpace(n={self.n}, k={self.k}, batch={self.batch_size}, channels={self.n_channels}, complex={self.complex_valued}, device={self.device})"
    
    def __str__(self) -> str:
        return f"Hilbert space H (dim={self.n}) with subspace H0 (dim={self.k}), batch={self.batch_size}, channels={self.n_channels}"


def tensor_product_hilbert(spaces: List[HilbertSpace]) -> HilbertSpace:
    n_total = 1
    k_total = 1
    batch_size = max(s.batch_size for s in spaces)
    n_channels = max(s.n_channels for s in spaces)
    for s in spaces:
        n_total *= s.n
        k_total *= s.k
    return HilbertSpace(
        n=n_total, k=k_total,
        batch_size=batch_size, n_channels=n_channels,
        complex_valued=spaces[0].complex_valued, device=spaces[0].device
    )


def direct_sum_hilbert(spaces: List[HilbertSpace]) -> HilbertSpace:
    n_total = sum(s.n for s in spaces)
    k_total = sum(s.k for s in spaces)
    batch_size = max(s.batch_size for s in spaces)
    n_channels = max(s.n_channels for s in spaces)
    return HilbertSpace(
        n=n_total, k=k_total,
        batch_size=batch_size, n_channels=n_channels,
        complex_valued=spaces[0].complex_valued, device=spaces[0].device
    )