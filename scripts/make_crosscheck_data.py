#!/usr/bin/env python3
"""
Reference values from independent, established software for tests/test_crosscheck_external.py.

The library torch_vn_algebra is NOT imported here: every reference value is produced by external code on inputs
generated with NumPy (fixed seeds) and stored together with those inputs, so that the tests compare the library
with the external software on identical inputs without needing the external package installed.

    QuTiP       Lindblad master equation (mesolve, liouvillian), partial traces, entropies, fidelity,
                relative entropy, trace distance                                   -> lindblad.json, composite.json
    QuSpin      XXZ chains in magnetisation and momentum sectors: spectra, ground states, correlations,
                entanglement entropies, time evolution; J1-J2 Heisenberg ring       -> xxz.json
    GAP         character tables (Irr), Frobenius-Schur indicators, classes, induced / permutation
                characters and their decomposition into irreducibles (via passagemath-gap / libgap)
                                                                                   -> groups.json
    SymPy       Clebsch-Gordan coefficients (sympy.physics.wigner), exact, j1, j2 <= 4  -> cg.json

Installation (any one of them may be missing; its section is then skipped):

    pip install qutip quspin sympy passagemath-gap

    python scripts/make_crosscheck_data.py                 # all sections
    python scripts/make_crosscheck_data.py lindblad xxz    # some sections

Conventions of the stored data (independent of the library):
* complex arrays are stored as {"re": [...], "im": [...]} (nested lists of the same shape);
* "label basis" of an L-site spin-1/2 chain: basis state x = sum_i b_i 2^i with b_i = 1 if site i is up;
* block-diagonal matrices of a direct sum (+)_c C^{k_c} use the order: sector 0 first, then sector 1, ...
"""
import json
import math
import pathlib
import sys
from fractions import Fraction

import numpy as np

OUT = pathlib.Path(__file__).resolve().parents[1] / 'tests' / 'data' / 'crosscheck'


def cplx(a):
    a = np.asarray(a, dtype=complex)
    return {"re": a.real.tolist(), "im": a.imag.tolist()}


def save(name, data):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    with open(path, 'w') as f:
        json.dump(data, f, separators=(',', ':'))
        f.write('\n')
    print(f"wrote {path} ({path.stat().st_size / 1024:.0f} kB)")


def rand_herm(rng, n, scale=1.0):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return scale * (X + X.conj().T) / 2


def rand_cplx(rng, n, scale=1.0):
    return scale * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))


def rand_state(rng, n, rank=None):
    G = rng.normal(size=(n, rank or n)) + 1j * rng.normal(size=(n, rank or n))
    R = G @ G.conj().T
    return R / np.trace(R).real


def block_diag(blocks):
    n = sum(b.shape[0] for b in blocks)
    out = np.zeros((n, n), dtype=complex)
    o = 0
    for b in blocks:
        k = b.shape[0]
        out[o:o + k, o:o + k] = b
        o += k
    return out


def mask_blocks(X, groups):
    """Keep only the entries (i, j) with i, j in the same group (groups: list of index lists)."""
    out = np.zeros_like(X)
    for g in groups:
        out[np.ix_(g, g)] = X[np.ix_(g, g)]
    return out


# ======================================================================================================
# QuTiP: Lindblad dynamics
# ======================================================================================================
def chain_ops_label_basis(L):
    """S^+_i, S^-_i, S^z_i of an L-site chain built with qutip.tensor (qutip order: site 0 = most significant
    factor, basis(2, 0) = up), returned in the label basis (index x, bit i of x = 1 iff site i is up)."""
    import qutip
    up, dn = qutip.basis(2, 0), qutip.basis(2, 1)
    sp1 = up * dn.dag()                                   # |up><down|
    sz1 = 0.5 * (up * up.dag() - dn * dn.dag())
    eye = qutip.qeye(2)
    # permutation qutip index -> label
    perm = np.zeros(2 ** L, dtype=int)                    # perm[label] = qutip index
    for x in range(2 ** L):
        perm[x] = sum((1 - ((x >> i) & 1)) << (L - 1 - i) for i in range(L))

    def site(op, i):
        M = qutip.tensor([op if j == i else eye for j in range(L)]).full()
        return M[np.ix_(perm, perm)]
    Sp = [site(sp1, i) for i in range(L)]
    Sm = [s.conj().T for s in Sp]
    Sz = [site(sz1, i) for i in range(L)]
    return Sp, Sm, Sz


def lindblad():
    import qutip
    rng = np.random.default_rng(20261015)
    opts = dict(atol=1e-13, rtol=1e-12, nsteps=10 ** 7, method='adams')
    data = {"software": f"qutip {qutip.__version__}", "cases": {}}

    # --- case 1: block-diagonal algebra M_2 (+) M_3 (+) M_1, constant and time-dependent H -------------
    dims = [2, 3, 1]
    H0 = block_diag([rand_herm(rng, k) for k in dims])
    H1 = block_diag([rand_herm(rng, k, 0.5) for k in dims])
    jumps = [block_diag([rand_cplx(rng, k, 0.4) for k in dims]) for _ in range(2)]
    rates = [0.7, 0.3]
    rho0 = block_diag([p * rand_state(rng, k) for p, k in zip([0.5, 0.3, 0.2], dims)])
    times = [0.0, 0.25, 0.5, 0.75, 1.0]
    omega = 1.3
    c_ops = [qutip.Qobj(math.sqrt(g) * L) for g, L in zip(rates, jumps)]
    res = qutip.mesolve(qutip.Qobj(H0), qutip.Qobj(rho0), times, c_ops, options=opts)
    res_t = qutip.mesolve([qutip.Qobj(H0), [qutip.Qobj(H1), lambda t: np.cos(omega * t)]], qutip.Qobj(rho0), times,
                          c_ops, options=opts)
    # channel exp(t L) at t = 0.8 on the matrix units E_ij of every block
    t_ch = 0.8
    S = (qutip.liouvillian(qutip.Qobj(H0), c_ops) * t_ch).expm()
    units, images = [], []
    off = 0
    for c, k in enumerate(dims):
        for i in range(k):
            for j in range(k):
                E = np.zeros((sum(dims), sum(dims)), dtype=complex)
                E[off + i, off + j] = 1.0
                img = qutip.vector_to_operator(S * qutip.operator_to_vector(qutip.Qobj(E))).full()
                units.append([c, i, j])
                images.append(img)
        off += k
    data["cases"]["blocks"] = dict(
        dims=dims, H0=cplx(H0), H1=cplx(H1), omega=omega, jumps=[cplx(L) for L in jumps], rates=rates,
        rho0=cplx(rho0), times=times, rho_t=[cplx(s.full()) for s in res.states],
        rho_t_timedep=[cplx(s.full()) for s in res_t.states],
        channel_t=t_ch, channel_units=units, channel_images=[cplx(m) for m in images],
        note="mesolve(H0) / mesolve(H0 + cos(omega t) H1) with c_ops sqrt(rate) L; exp(t liouvillian) on E_ij")

    # --- case 2: XXZ chain L = 4 with losses (jumps between sectors), in the full 2^L space ---------------
    L = 4
    Sp, Sm, Sz = chain_ops_label_basis(L)
    J, Delta, h = 1.0, 0.6, [0.3, -0.2, 0.5, 0.1]
    H = sum(J / 2 * (Sp[i] @ Sm[i + 1] + Sm[i] @ Sp[i + 1]) + J * Delta * Sz[i] @ Sz[i + 1] for i in range(L - 1))
    H = H + sum(h[i] * Sz[i] for i in range(L))
    jump_list = [("lowering", 0, Sm[0], 0.4), ("lowering", 2, Sm[2], 0.25),
                 ("collective_lowering", [1, 3], Sm[1] + Sm[3], 0.15),
                 ("dephasing", 1, Sz[1], 0.3), ("hop", [0, 1], Sp[0] @ Sm[1], 0.2)]
    pop = np.array([bin(x).count('1') for x in range(2 ** L)])
    sectors = [list(np.nonzero(pop == N)[0]) for N in range(L + 1)]
    rho0 = mask_blocks(rand_state(rng, 2 ** L), sectors)
    rho0 /= np.trace(rho0).real
    times = [0.0, 0.5, 1.0, 1.5, 2.0]
    c_ops = [qutip.Qobj(math.sqrt(g) * M) for _, _, M, g in jump_list]
    res = qutip.mesolve(qutip.Qobj(H), qutip.Qobj(rho0), times, c_ops, options=opts)
    states = [s.full() for s in res.states]
    off_sector = max(float(np.abs(s - mask_blocks(s, sectors)).max()) for s in states)
    data["cases"]["chain"] = dict(
        L=L, J=J, Delta=Delta, h=h, boundary='open', H=cplx(H),
        jumps=[dict(kind=kind, sites=sites, rate=g, matrix=cplx(M)) for kind, sites, M, g in jump_list],
        rho0=cplx(rho0), times=times, rho_t=[cplx(s) for s in states], max_off_sector=off_sector,
        note="label basis; qutip mesolve in the full 2^L space; max_off_sector = largest entry of rho(t) between "
             "different particle-number sectors")

    # --- case 3: a real algebra M_3(R) (+) M_2(R), no Hamiltonian, real jumps -------------------------------
    rng = np.random.default_rng(99)
    dims = [3, 2]
    jumps = [block_diag([0.5 * rng.normal(size=(k, k)) for k in dims]) for _ in range(3)]
    rates = [1.0, 0.5, 0.25]
    G = [rng.normal(size=(k, k + 1)) for k in dims]
    rho0 = block_diag([p * g @ g.T / np.trace(g @ g.T) for p, g in zip([0.7, 0.3], G)]).real
    times = [0.0, 0.5, 1.0]
    c_ops = [qutip.Qobj(math.sqrt(g) * L) for g, L in zip(rates, jumps)]
    res = qutip.mesolve(qutip.Qobj(np.zeros((5, 5))), qutip.Qobj(rho0), times, c_ops, options=opts)
    S = (qutip.liouvillian(None, c_ops) * 1.0).expm()
    data["cases"]["real_blocks"] = dict(
        dims=dims, jumps=[L.real.tolist() for L in jumps], rates=rates, rho0=rho0.tolist(), times=times,
        rho_t=[cplx(st.full()) for st in res.states],
        channel_t=1.0,
        channel_rho0=cplx(qutip.vector_to_operator(S * qutip.operator_to_vector(qutip.Qobj(rho0))).full()),
        note="real matrices, H = 0")

    # --- case 4: a channel as jump (amplitude damping of site 1 of an L = 3 ring): GKSL with jumps K_0, K_1 --
    L, gamma = 3, 0.3
    Sp, Sm, Sz = chain_ops_label_basis(L)
    Pup = Sp[1] @ Sm[1]
    K0 = (np.eye(2 ** L) - Pup) + math.sqrt(1 - gamma) * Pup
    K1 = math.sqrt(gamma) * Sm[1]
    J, Delta, hz = 1.0, 0.4, 0.3
    H = sum(J / 2 * (Sp[i] @ Sm[(i + 1) % L] + Sm[i] @ Sp[(i + 1) % L]) + J * Delta * Sz[i] @ Sz[(i + 1) % L]
            for i in range(L)) + hz * sum(Sz)
    pop = np.array([bin(x).count('1') for x in range(2 ** L)])
    sectors = [list(np.nonzero(pop == N)[0]) for N in range(L + 1)]
    rho0 = mask_blocks(rand_state(rng, 2 ** L), sectors)
    rho0 /= np.trace(rho0).real
    rate = 0.5
    times = [0.0, 1.0, 2.0]
    res = qutip.mesolve(qutip.Qobj(H), qutip.Qobj(rho0), times, [qutip.Qobj(math.sqrt(rate) * K) for K in (K0, K1)],
                        options=opts)
    data["cases"]["channel_jump"] = dict(
        L=L, boundary='periodic', J=J, Delta=Delta, h=hz, site=1, gamma=gamma, rate=rate, K0=K0.real.tolist(),
        K1=K1.real.tolist(), rho0=cplx(rho0), times=times, rho_t=[cplx(st.full()) for st in res.states],
        note="label basis; qutip mesolve with c_ops sqrt(rate) K_0, sqrt(rate) K_1 (Kraus operators of amplitude "
             "damping with probability gamma on one site)")
    save('lindblad.json', data)


# ======================================================================================================
# QuTiP: tensor products, partial traces, entropies
# ======================================================================================================
def composite():
    import qutip
    rng = np.random.default_rng(1729)
    data = {"software": f"qutip {qutip.__version__}", "cases": {}}
    configs = {
        # name: (k1, charges1, k2, charges2, product)
        "plain": ([2, 1], [0, 1], [2, 2], [0, 1], 'plain'),
        "fused": ([2, 1], [0, 1], [2, 2], [0, 1], 'fused'),
        "fused_repeated_charges": ([2, 1], [0, 0], [2, 2], [0, 1], 'fused'),
    }
    for name, (k1, q1, k2, q2, kind) in configs.items():
        n1, n2 = sum(k1), sum(k2)
        o1 = np.cumsum([0] + k1)
        o2 = np.cumsum([0] + k2)
        # sector of every basis vector of H1 (x) H2 (index i1 * n2 + i2), and its pair / fused label
        sec_of = {}
        for i1 in range(n1):
            c = int(np.searchsorted(o1, i1, side='right') - 1)
            for i2 in range(n2):
                d = int(np.searchsorted(o2, i2, side='right') - 1)
                key = (c, d) if kind == 'plain' else q1[c] + q2[d]
                sec_of.setdefault(key, []).append(i1 * n2 + i2)
        groups = list(sec_of.values())
        n = n1 * n2

        def state():
            R = mask_blocks(rand_state(rng, n), groups)
            return R / np.trace(R).real
        rho, sigma = state(), state()
        X = mask_blocks(rand_cplx(rng, n), groups)
        A = block_diag([rand_cplx(rng, k) for k in k1])
        B = block_diag([rand_cplx(rng, k) for k in k2])
        dims = [[n1, n2], [n1, n2]]
        Q = lambda M: qutip.Qobj(M, dims=dims)                                     # noqa: E731
        r, s = Q(rho), Q(sigma)
        data["cases"][name] = dict(
            k1=k1, charges1=q1, k2=k2, charges2=q2, product=kind,
            rho=cplx(rho), sigma=cplx(sigma), X=cplx(X), A=cplx(A), B=cplx(B),
            ptrace_rho_keep1=cplx(r.ptrace(0).full()), ptrace_rho_keep2=cplx(r.ptrace(1).full()),
            ptrace_X_keep1=cplx(Q(X).ptrace(0).full()), ptrace_X_keep2=cplx(Q(X).ptrace(1).full()),
            kron_AB=cplx(qutip.tensor(qutip.Qobj(A), qutip.Qobj(B)).full()),
            entropy_rho=float(qutip.entropy_vn(r)), entropy_sigma=float(qutip.entropy_vn(s)),
            entropy_rho1=float(qutip.entropy_vn(r.ptrace(0))), entropy_rho2=float(qutip.entropy_vn(r.ptrace(1))),
            mutual_information=float(qutip.entropy_mutual(r, 0, 1)),
            fidelity_squared=float(qutip.fidelity(r, s) ** 2),
            relative_entropy=float(qutip.entropy_relative(r, s)),
            relative_entropy_reduced1=float(qutip.entropy_relative(r.ptrace(0), s.ptrace(0))),
            trace_distance=float(qutip.tracedist(r, s)),
            purity=float((r * r).tr().real),
            note="H1 = (+)_c C^{k1_c}, H2 = (+)_d C^{k2_d}, H1 (x) H2 in Kronecker order (index i1 * n2 + i2); "
                 "qutip.fidelity is sqrt-fidelity, squared here; natural logarithms")
    save('composite.json', data)


# ======================================================================================================
# QuSpin: XXZ chains
# ======================================================================================================
def _quspin_static(L, bonds, J, Delta, h):
    static = [["+-", [[J / 2, i, j] for i, j in bonds]], ["-+", [[J / 2, i, j] for i, j in bonds]],
              ["zz", [[J * Delta, i, j] for i, j in bonds]]]
    if h is not None:
        static.append(["z", [[float(h[i]), i] for i in range(L)]])
    return static


def xxz():
    import quspin
    from quspin.basis import spin_basis_1d
    from quspin.operators import hamiltonian
    kw = dict(check_symm=False, check_herm=False, check_pcon=False)
    data = {"software": f"quspin {getattr(quspin, '__version__', '?')}", "cases": {}}
    rng = np.random.default_rng(7)

    # --- open chain L = 10 with random fields: sector spectra, ground state, correlations, entanglement ----
    L, J, Delta = 10, 1.0, 0.8
    h = rng.uniform(-0.5, 0.5, size=L).round(6)
    bonds = [(i, i + 1) for i in range(L - 1)]
    case = dict(L=L, boundary='open', J=J, Delta=Delta, h=h.tolist(), spectra={})
    for N in (3, 4, 5):
        basis = spin_basis_1d(L, Nup=N, pauli=0)
        H = hamiltonian(_quspin_static(L, bonds, J, Delta, h), [], basis=basis, dtype=np.float64, **kw)
        case["spectra"][str(N)] = np.linalg.eigvalsh(H.toarray()).tolist()
    N = 5
    basis = spin_basis_1d(L, Nup=N, pauli=0)
    H = hamiltonian(_quspin_static(L, bonds, J, Delta, h), [], basis=basis, dtype=np.float64, **kw)
    E, V = np.linalg.eigh(H.toarray())
    assert E[1] - E[0] > 1e-3, "ground state must be non-degenerate"
    psi = V[:, 0]
    op = lambda static: hamiltonian(static, [], basis=basis, dtype=np.float64, **kw)          # noqa: E731
    zz = np.zeros((L, L))
    pm = np.zeros((L, L))
    for i in range(L):
        for j in range(L):
            zz[i, j] = op([["zz", [[1.0, i, j]]]] if i != j else [["I", [[0.25, i]]]]).expt_value(psi).real
            if i != j:
                pm[i, j] = op([["+-", [[1.0, i, j]]]]).expt_value(psi).real
            else:
                pm[i, j] = op([["z", [[1.0, i]]], ["I", [[0.5, i]]]]).expt_value(psi).real
    ent = [float(basis.ent_entropy(psi, sub_sys_A=list(range(n)), density=False)["Sent_A"]) for n in range(1, L)]
    case["ground_state"] = dict(N=N, E0=float(E[0]), E1=float(E[1]), SzSz=zz.tolist(), SpSm=pm.tolist(),
                                entanglement=ent)
    # time evolution of the Neel state 1010... (site 0 up)
    bits = [1 - (i % 2) for i in range(L)]
    # QuSpin stores site i at bit L-1-i of the state integer
    s_int = sum(b << (L - 1 - i) for i, b in enumerate(bits))
    psi0 = np.zeros(basis.Ns, dtype=complex)
    psi0[basis.index(s_int)] = 1.0
    sz_ops = [op([["z", [[1.0, i]]]]) for i in range(L)]
    check = [sz.expt_value(psi0).real for sz in sz_ops]
    assert np.allclose(check, [b - 0.5 for b in bits]), "QuSpin site convention"
    times = [0.0, 0.5, 1.0, 2.0, 3.0]
    Hc = hamiltonian(_quspin_static(L, bonds, J, Delta, h), [], basis=basis, dtype=np.complex128, **kw)
    psit = Hc.evolve(psi0, 0.0, times)                                         # (Ns, T)
    case["evolution"] = dict(N=N, bits=bits, times=times,
                             Sz=[[float(sz.expt_value(psit[:, t]).real) for sz in sz_ops] for t in range(len(times))],
                             entanglement_half=[float(basis.ent_entropy(psit[:, t], sub_sys_A=list(range(L // 2)),
                                                                        density=False)["Sent_A"])
                                                for t in range(len(times))])
    data["cases"]["open10"] = case

    # --- periodic ring L = 12: momentum sectors -------------------------------------------------------------
    L, J, Delta, hz = 12, 1.0, 0.5, 0.2
    bonds = [(i, (i + 1) % L) for i in range(L)]
    case = dict(L=L, boundary='periodic', J=J, Delta=Delta, h=hz, sectors={})
    for N in (5, 6):
        blocks = {}
        for m in range(L):
            basis = spin_basis_1d(L, Nup=N, kblock=m, pauli=0)
            if basis.Ns == 0:
                blocks[str(m)] = []
                continue
            H = hamiltonian(_quspin_static(L, bonds, J, Delta, [hz] * L), [], basis=basis, dtype=np.complex128, **kw)
            blocks[str(m)] = np.linalg.eigvalsh(H.toarray()).tolist()
        full = spin_basis_1d(L, Nup=N, pauli=0)
        H = hamiltonian(_quspin_static(L, bonds, J, Delta, [hz] * L), [], basis=full, dtype=np.float64, **kw)
        case["sectors"][str(N)] = dict(momentum_spectra=blocks, spectrum=np.linalg.eigvalsh(H.toarray()).tolist())
    data["cases"]["ring12"] = case

    # --- J1-J2 Heisenberg ring L = 10 (SU(2)): full spectrum, ground-state correlations --------------------
    L, J1, J2 = 10, 1.0, 0.2
    b1 = [(i, (i + 1) % L) for i in range(L)]
    b2 = [(i, (i + 2) % L) for i in range(L)]
    static = _quspin_static(L, b1, J1, 1.0, None) + _quspin_static(L, b2, J2, 1.0, None)
    spec = []
    gs = None
    for N in range(L + 1):
        basis = spin_basis_1d(L, Nup=N, pauli=0)
        H = hamiltonian(static, [], basis=basis, dtype=np.float64, **kw)
        e, v = np.linalg.eigh(H.toarray())
        spec += e.tolist()
        if gs is None or e[0] < gs[0] - 1e-12:
            gs = (e[0], v[:, 0], basis, N)
    E0, psi, basis, N = gs
    assert N == L // 2
    corr = []
    for r in range(1, L // 2 + 1):
        ops = [["+-", [[0.5, 0, r]]], ["-+", [[0.5, 0, r]]], ["zz", [[1.0, 0, r]]]]
        corr.append(float(hamiltonian(ops, [], basis=basis, dtype=np.float64, **kw).expt_value(psi).real))
    data["cases"]["j1j2_ring10"] = dict(L=L, J1=J1, J2=J2, spectrum=sorted(spec), E0=float(E0),
                                        S0_dot_Sr=corr, note="<S_0 . S_r> in the (singlet) ground state, r = 1..L/2")

    # --- sweep: momentum blocks of all sectors of periodic rings L = 2..8 (J = 0.8, Delta = 1.3, h = 0.37) -----
    J, Delta, hz = 0.8, 1.3, 0.37
    sweep = {}
    for L in range(2, 9):
        bonds = [(i, (i + 1) % L) for i in range(L)] if L > 2 else [(0, 1)]
        for N in range(L + 1):
            blocks = {}
            for m in range(L):
                basis = spin_basis_1d(L, Nup=N, kblock=m, pauli=0)
                if basis.Ns == 0:
                    blocks[str(m)] = []
                    continue
                H = hamiltonian(_quspin_static(L, bonds, J, Delta, [hz] * L), [], basis=basis, dtype=np.complex128,
                                **kw)
                blocks[str(m)] = np.linalg.eigvalsh(H.toarray()).tolist()
            sweep[f"{L},{N}"] = blocks
    data["cases"]["momentum_sweep"] = dict(J=J, Delta=Delta, h=hz, blocks=sweep,
                                           note="key 'L,N'; periodic ring, L = 2 has the single bond (0, 1)")

    # --- sweep: Heisenberg chains L = 1..8, open / periodic, J2 = 0 / 0.35, uniform field h = 0.21 ----------
    sweep = {}
    for L in range(1, 9):
        for bc in ('open', 'periodic'):
            for J2 in (0.0, 0.35):
                if J2 and L < 4:
                    continue
                if bc == 'periodic':
                    b1 = [(i, (i + 1) % L) for i in range(L)] if L > 2 else [(i, i + 1) for i in range(L - 1)]
                    b2 = [(i, (i + 2) % L) for i in range(L)] if L > 4 else [(i, i + 2) for i in range(L - 2)]
                else:
                    b1 = [(i, i + 1) for i in range(L - 1)]
                    b2 = [(i, i + 2) for i in range(L - 2)]
                spec = []
                for N in range(L + 1):
                    basis = spin_basis_1d(L, Nup=N, pauli=0)
                    st = (_quspin_static(L, b1, 1.0, 1.0, [0.21] * L) if b1
                          else [["z", [[0.21, i] for i in range(L)]]])
                    if J2:
                        st = st + _quspin_static(L, b2, J2, 1.0, None)
                    Hm = hamiltonian(st, [], basis=basis, dtype=np.float64, **kw).toarray()
                    spec += np.linalg.eigvalsh(Hm).tolist()
                sweep[f"{L},{bc},{J2}"] = dict(spectrum=sorted(spec), nn_bonds=b1, nnn_bonds=b2 if J2 else [])
    data["cases"]["heisenberg_sweep"] = dict(J=1.0, h=0.21, chains=sweep,
                                             note="key 'L,boundary,J2'; full spectra (2^L levels) with the bonds given")
    save('xxz.json', data)


# ======================================================================================================
# GAP: character tables, indicators, induced characters
# ======================================================================================================
GAP_GROUPS = {
    # name: (GAP expression of a permutation group, [(kind, order of h, power of E(order)) ...]);
    # every group also gets the tensor square of its natural representation ("natural2")
    "S4": ("SymmetricGroup(4)", [("natural",), ("regular",), ("induced", 4, 1), ("induced", 3, 0)]),
    "S5": ("SymmetricGroup(5)", [("natural",), ("regular",), ("induced", 5, 1), ("induced", 6, 1),
                                 ("induced", 4, 0)]),
    "A5": ("AlternatingGroup(5)", [("natural",), ("regular",), ("induced", 5, 1), ("induced", 3, 1)]),
    "SL(2,3)": ("Image(SmallerDegreePermutationRepresentation(Image(IsomorphismPermGroup(SL(2,3)))))",
                [("natural",), ("regular",), ("induced", 6, 1), ("induced", 4, 1), ("induced", 3, 1)]),
    "D12": ("DihedralGroup(IsPermGroup, 12)", [("natural",), ("regular",), ("induced", 6, 1), ("induced", 2, 1)]),
    "Q8xC2": ("DirectProduct(QuaternionGroup(IsPermGroup, 8), CyclicGroup(IsPermGroup, 2))",
              [("natural",), ("regular",), ("induced", 4, 1), ("induced", 2, 1)]),
}

GAP_CODE = r'''function(G, specs)
  local tbl, cc, els, deg, irr, cls, reps, spec, h, n, H, lam, ind, out, hidx;
  tbl := CharacterTable(G);
  cc := ConjugacyClasses(tbl);
  els := AsSSortedList(G);
  deg := LargestMovedPoint(G);
  cls := List(els, g -> First([1..Length(cc)], i -> g in cc[i]));
  irr := List(Irr(tbl), chi -> List(ValuesOfClassFunction(chi), x -> [Conductor(x), CoeffsCyc(x, Conductor(x))]));
  reps := [];
  for spec in specs do
    if spec[1] = "natural" then
      ind := ClassFunction(tbl, List(cc, c -> deg - NrMovedPoints(Representative(c))));
      Add(reps, rec(kind := "natural", mults := List(Irr(tbl), chi -> ScalarProduct(ind, chi))));
    elif spec[1] = "regular" then
      Add(reps, rec(kind := "regular", mults := List(Irr(tbl), chi -> chi[1])));
    elif spec[1] = "natural2" then
      ind := ClassFunction(tbl, List(cc, c -> (deg - NrMovedPoints(Representative(c)))^2));
      Add(reps, rec(kind := "natural2", mults := List(Irr(tbl), chi -> ScalarProduct(ind, chi))));
    else
      n := spec[2];
      hidx := First([1..Length(els)], i -> Order(els[i]) = n);
      h := els[hidx];
      H := Subgroup(G, [h]);
      lam := First(Irr(H), chi -> h^chi = E(n)^spec[3]);
      ind := InducedClassFunction(lam, G);
      Add(reps, rec(kind := "induced", order := n, power := spec[3], h := hidx - 1,
                    mults := List(Irr(tbl), chi -> ScalarProduct(ind, chi))));
    fi;
  od;
  return rec(deg := deg, perms := List(els, g -> ListPerm(g, deg) - 1), classes := cls - 1, irr := irr,
             indicators := Indicator(tbl, 2), class_sizes := SizesConjugacyClasses(tbl), reps := reps,
             identification := IdGroup(G), structure := StructureDescription(G),
             generators := List(SmallGeneratingSet(G), g -> Position(els, g) - 1));
end'''


def _py(x):
    """Sage/GAP values -> plain Python (int, str, lists)."""
    if isinstance(x, (list, tuple)):
        return [_py(v) for v in x]
    if isinstance(x, str):
        return x
    try:
        return int(x)
    except (TypeError, ValueError):
        return str(x)


def groups():
    from sage.all__sagemath_gap import libgap
    f = libgap.eval(GAP_CODE)
    data = {"software": f"GAP {libgap.eval('GAPInfo.Version')} (passagemath-gap)",
            "note": "perms: 0-based images of the points 0..deg-1 (GAP ListPerm), elements in AsSSortedList order; "
                    "the group law is GAP's: (g*h)(i) = h(g(i)). irr: complex character values per class (rows = "
                    "Irr(CharacterTable(G))). indicators: Frobenius-Schur indicators. reps[].mults: multiplicities "
                    "<chi_rep, chi_irr> computed by GAP (ScalarProduct), for the natural permutation character on "
                    "deg points, its square (natural2: the tensor square of the permutation representation), the "
                    "regular character and InducedClassFunction(lambda, G) with lambda(h) = "
                    "E(order)^power on the cyclic group <h>, h = element number h (0-based) of the sorted list. "
                    "generators: SmallGeneratingSet(G) as 0-based element numbers.",
            "groups": {}}
    for name, (expr, specs) in GAP_GROUPS.items():
        G = libgap.eval(expr)
        r = f(G, [list(s) for s in specs] + [["natural2"]]).sage()
        irr = []
        for row in r['irr']:
            vals = []
            for cond, coeffs in row:
                z = sum(int(a) * np.exp(2j * np.pi * k / int(cond)) for k, a in enumerate(coeffs))
                vals.append([float(z.real), float(z.imag)])
            irr.append(vals)
        data["groups"][name] = dict(
            gap=expr, order=len(r['perms']), structure=str(r['structure']), id=[int(x) for x in r['identification']],
            deg=int(r['deg']), perms=[[int(x) for x in p] for p in r['perms']],
            classes=[int(c) for c in r['classes']], class_sizes=[int(s) for s in r['class_sizes']], irr=irr,
            indicators=[int(x) for x in r['indicators']], generators=[int(x) for x in r['generators']],
            reps=[{str(k): _py(v) for k, v in rep.items()} for rep in r['reps']])
    save('groups.json', data)


# ======================================================================================================
# SymPy: Clebsch-Gordan coefficients
# ======================================================================================================
def cg():
    import sympy
    from sympy import Rational, nsimplify
    from sympy.physics.wigner import clebsch_gordan
    TWO_J_MAX = 8
    table = {}
    for a in range(TWO_J_MAX + 1):
        for b in range(TWO_J_MAX + 1):
            for e in range(abs(a - b), a + b + 1, 2):
                vals = []
                for x in range(-a, a + 1, 2):
                    for y in range(-b, b + 1, 2):
                        if abs(x + y) > e:
                            continue
                        v = clebsch_gordan(Rational(a, 2), Rational(b, 2), Rational(e, 2), Rational(x, 2),
                                           Rational(y, 2), Rational(x + y, 2))
                        sq = Fraction(str(nsimplify(v ** 2)))
                        sign = 0 if v == 0 else (1 if v > 0 else -1)
                        vals.append(f"{'-' if sign < 0 else ''}{sq}")
                table[f"{a},{b},{e}"] = vals
    save('cg.json', {"software": f"sympy {sympy.__version__} (sympy.physics.wigner.clebsch_gordan)",
                     "note": "key '2j1,2j2,2J'; values for m1 = -j1..j1, m2 = -j2..j2 (outer loop m1) with "
                             "|m1+m2| <= J, "
                             "M = m1 + m2; entry 's' = sign * sqrt(|s|) as an exact fraction",
                     "two_j_max": TWO_J_MAX, "table": table})


SECTIONS = dict(lindblad=lindblad, composite=composite, xxz=xxz, groups=groups, cg=cg)

if __name__ == '__main__':
    names = sys.argv[1:] or list(SECTIONS)
    for name in names:
        try:
            SECTIONS[name]()
        except ImportError as exc:
            print(f"skipped {name}: {exc}")
