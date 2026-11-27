"""Cutoff-free lower bound on the Weyl-ordered partial focusing Q_v.

Split [eps, v] into N intervals I_n = [c-l/2, c+l/2].  On each, with
s = (t - c)/l + 1/2 in [0,1] and bump b(s) = sin^2(pi s):
    f(t) = t * b'(s),  g(t) = t * b''(s)/(2 pi^2)     (int f/t = int g/t = 0)
    A_n = int f phi,  B_n = int g phi,   phi = dx - x/t = t dh.
Exact bracket (kappa = hbar = 1):   {A,B} = (1/2)[ int f' g dt + int f G / t dt ],
    G(t) = int_{c-l/2}^t g(u)/u du   (compact, since int g/t = 0).
Cauchy-Schwarz in <u,w> = int u w / t on I_n gives the classical form inequality
    (1/4) int_{I_n} phi^2 / t  >=  (1/4) (A,B) Gram^{-1} (A,B)^T,
    Gram = [[int t f^2, int t f g], [int t f g, int t g^2]].
Weyl quantization preserves the order of quadratic forms, and the Weyl operator of the
2x2 form has spectrum >= nu_n/2 with nu_n = (1/4) |{A,B}| / sqrt(det Gram).
Pairs on disjoint intervals commute, hence  Q_v(Weyl) >= sum_n nu_n / 2  (times identity).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

EPS = 0.1
xs, ws = leggauss(200)

def interval_nu(c, l):
    s = (xs + 1) / 2; t = c - l / 2 + l * s; w = ws * l / 2
    b1 = np.pi * np.sin(2 * np.pi * s)                  # b'(s)
    b2 = 2 * np.pi**2 * np.cos(2 * np.pi * s) / (2 * np.pi**2)   # b''(s)/(2 pi^2)
    b2p = -2 * np.pi * np.sin(2 * np.pi * s)            # d/ds of b2
    b1p = 2 * np.pi**2 * np.cos(2 * np.pi * s)          # b''(s)
    f = t * b1; g = t * b2
    fp = b1 + t * b1p / l                              # df/dt
    # G(t) = int g/u du = int b2(s) ds * l  -> antiderivative of cos: sin(2 pi s)/(2 pi)
    G = l * np.sin(2 * np.pi * s) / (2 * np.pi)
    br = 0.5 * (np.sum(w * fp * g) + np.sum(w * f * G / t))
    gram = np.array([[np.sum(w * t * f * f), np.sum(w * t * f * g)],
                     [np.sum(w * t * f * g), np.sum(w * t * g * g)]])
    return 0.25 * abs(br) / np.sqrt(np.linalg.det(gram))

for v in [0.3, 0.55, 0.8, 1.0]:
    row = []
    for N in [10, 100, 1000, 10000]:
        l = (v - EPS) / N
        cs = EPS + l * (np.arange(N) + 0.5)
        lb = sum(interval_nu(c, l) for c in cs) / 2
        row.append(f"N={N}: {lb:.4g} (/N^2 = {lb/N**2:.4g})")
    print(f"v={v}: " + "  ".join(row))
