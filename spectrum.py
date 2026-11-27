"""Zero-point focusing of the Weyl-ordered clock in arXiv:2611.34924.

One radiative polarization z1 in one angular harmonic, corner data fixed
(the setting of the paper's Appendix D.3).  With x = v h on [eps, 1],
x(eps) = x(1) = 0, the symplectic form is (D.42a)-(D.43)
    Omega(x, y) = kappa * int [(dx) y - (dy) x] dv = 2 kappa int (dx) y dv,
with kappa ~ rho_S^2 / G, so hbar/kappa ~ 1/s.  We set hbar = 1, kappa = 1.
The sphere-averaged focusing (4.22) at v is, for this harmonic,
    Q_v[x] = (1/4) int_eps^v (1/t) (dx - x/t)^2 dt     ( = (1/4) int t (dh)^2 ).
Classical coordinates: coefficients c_n of x = sum c_n sin(n pi tau).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

EPS = 0.1
L = 1.0 - EPS

def omega_matrix(N):
    m = np.arange(1, N + 1)[:, None]; n = np.arange(1, N + 1)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        I = np.where(m == n, 0.0, n * (1 - (-1.0) ** (m + n)) / (np.pi * (n**2 - m**2)))
    # int e_m' e_n dt = m pi * I_mn ; Omega_mn = kappa*(int e_m' e_n - int e_n' e_m)
    A = m * np.pi * I
    return A - A.T

def q_matrix(N, v, npts=None):
    npts = npts or max(400, 8 * N)
    xs, ws = leggauss(npts)
    t = EPS + (v - EPS) * (xs + 1) / 2; w = ws * (v - EPS) / 2
    tau = (t - EPS) / L
    k = np.arange(1, N + 1)[:, None]
    e = np.sin(k * np.pi * tau); de = (k * np.pi / L) * np.cos(k * np.pi * tau)
    phi = de - e / t                    # (dx - x/t) for each basis fn
    return 0.25 * (phi * (w / t)) @ phi.T

def comm_matrix(N):
    # [c_m, c_n] = i J_mn with J = Omega^{-1} (up to a sign that does not affect nu)
    return np.linalg.inv(omega_matrix(N))

def sympl_eigs(Qm, J):
    ev = np.linalg.eigvals(J @ (2 * Qm))
    nu = np.sort(np.abs(ev.imag))[::2]   # +-i nu pairs
    return nu[nu > 1e-12]

def williamson_vacuum_cov(Qm, J):
    H = 2 * Qm
    w, U = np.linalg.eigh(H)
    Hh = U @ np.diag(np.sqrt(w)) @ U.T; Hmh = U @ np.diag(1 / np.sqrt(w)) @ U.T
    K = Hh @ J @ Hh
    K = (K - K.T) / 2
    import scipy.linalg as sl
    T, O = sl.schur(K, output="real")
    n = K.shape[0]; d = np.zeros(n)
    for i in range(0, n, 2):
        om = abs(T[i, i + 1]); d[i] = d[i + 1] = np.sqrt(om)
    S = Hmh @ O @ np.diag(d)
    return 0.5 * S @ S.T

if __name__ == "__main__":
    print("(P1) symplectic eigenvalues of the full-cone focusing Q_1  (hbar=kappa=1)")
    for N in [20, 40, 80, 160]:
        J = comm_matrix(N); nu = sympl_eigs(q_matrix(N, 1.0), J)
        print(f"  N={N:4d}: nu_1..4 = {np.round(nu[:4],4)}  nu_last/k_last = {nu[-1]/len(nu):.4f}  zero-point = {nu.sum()/2:10.3f}  zp/N^2 = {nu.sum()/2/N**2:.5f}")
    print("\n(P2) zero-point of partial focusing Q_v")
    for v in [0.3, 0.55, 0.8]:
        row = []
        for N in [20, 40, 80, 160]:
            J = comm_matrix(N); nu = sympl_eigs(q_matrix(N, v), J)
            row.append(f"N={N}: {nu.sum()/2:9.3f}")
        print(f"  v={v}: " + "  ".join(row))
    print("\n(P3) inf over states of <:Q_v:> with normal ordering w.r.t. the vacuum of Q_1")
    for N in [20, 40, 80, 160]:
        J = comm_matrix(N); Q1 = q_matrix(N, 1.0); Sig = williamson_vacuum_cov(Q1, J)
        out = []
        for v in [0.3, 0.55, 0.8, 1.0]:
            Qv = q_matrix(N, v); nu = sympl_eigs(Qv, J)
            out.append(f"v={v}: {nu.sum()/2 - np.trace(Qv @ Sig):9.3f}")
        print(f"  N={N:4d}: " + "  ".join(out))
