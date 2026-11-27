"""Cutoff-free lower bound on the model focusing operator Q_v (Lemma 1 of PROOF.md).

Split [eps, v] into N intervals I_n = [c - l/2, c + l/2], tau = (t - c)/l + 1/2 in [0, 1]:
    u1(tau) = sin^2(pi tau) sin(2 pi tau)              (odd about 1/2, mean 0, C^1, vanishes at ends)
    u2(tau) = sin^2(pi tau) (cos(2 pi tau) + 1/2)      (even, mean 0, C^1, vanishes at ends)
    f = t u1,  g = t u2,  A_n = int f phi,  B_n = int g phi,  phi = dx - x/t = t dh.
Exact bracket (kappa = hbar = 1):  {A,B} = (1/2)[ int f' g dt + int f G / t dt ],
    G(t) = int_{c-l/2}^t g(u)/u du = l * U2(tau),  U2 = antiderivative of u2 (vanishes at both ends).
Bessel inequality in <p,q> = int p q / t on I_n:
    (1/4) int_{I_n} phi^2 / t  >=  (1/4) (A,B) Gram^{-1} (A,B)^T,   Gram = Gram{t f, t g}.
Ground energy of z^T M z with [A,B] = i c is |c| sqrt(det M); here M = Gram^{-1}/4, so
    E0_n = |{A,B}| / (4 sqrt(det Gram)).
Lemma 1: Q_v >= sum_n E0_n  (as quadratic forms on every vector).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

EPS = 0.1
xs, ws = leggauss(400)
s = (xs + 1) / 2; wq = ws / 2
pi = np.pi
u1 = np.sin(pi*s)**2 * np.sin(2*pi*s)
u1p = 2*pi*np.sin(pi*s)*np.cos(pi*s)*np.sin(2*pi*s) + 2*pi*np.sin(pi*s)**2*np.cos(2*pi*s)
u2 = np.sin(pi*s)**2 * (np.cos(2*pi*s) + 0.5)
# U2(tau) = int_0^tau u2 ; closed form via sin^2 = (1 - cos 2 pi s)/2
U2 = (0.5*(np.sin(2*pi*s)/(2*pi)) - 0.5*(s/2 + np.sin(4*pi*s)/(8*pi))) \
     + 0.5*(0.5*s - np.sin(2*pi*s)/(4*pi))

def interval_E0(c, l):
    t = c - l/2 + l*s; w = wq * l
    f, g = t*u1, t*u2
    fp = u1 + t*u1p/l
    G = l*U2
    br = 0.5*(np.sum(w*fp*g) + np.sum(w*f*G/t))
    gram = np.array([[np.sum(w*t*f*f), np.sum(w*t*f*g)],
                     [np.sum(w*t*f*g), np.sum(w*t*g*g)]])
    return abs(br) / (4*np.sqrt(np.linalg.det(gram)))

# asymptotic constant: E0 ~ alpha / (8 c l) with alpha = I1 / sqrt(D)
I1 = np.sum(wq*u1p*u2)
D = np.sum(wq*u1*u1)*np.sum(wq*u2*u2) - np.sum(wq*u1*u2)**2
alpha = I1/np.sqrt(D)
print(f"checks: mean u1 = {np.sum(wq*u1):.2e}, mean u2 = {np.sum(wq*u2):.2e}, U2(1) = {U2[-1]:.2e}")
print(f"alpha = I1/sqrt(D) = {alpha:.6f};  predicted  sum_n E0_n ~ (alpha/8) N^2 ln(v/eps)/(v-eps)^2\n")
for v in [0.3, 0.55, 0.8, 1.0]:
    pred = alpha/8*np.log(v/EPS)/(v-EPS)**2
    row = []
    for N in [10, 100, 1000, 10000]:
        l = (v - EPS)/N
        cs = EPS + l*(np.arange(N) + 0.5)
        lb = sum(interval_E0(c, l) for c in cs)
        row.append(f"N={N}: /N^2 = {lb/N**2:.5f}")
    print(f"v={v}: " + "  ".join(row) + f"   predicted {pred:.5f}")
