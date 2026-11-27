# The Weyl-ordered gravitational clock of arXiv:2611.34924 is infinite on every state

**Paper:** V. Chandrasekaran, *An algebraic quantization of gravity from the null cone*, arXiv:2611.34924v1. Equation numbers below refer to that paper.

**Result.** The clock in (4.22)–(4.27) squares the shear operator, and those squares carry an infinite zero-point term. Once every radial mode on the null cone is kept, three things follow:

1. The clock's quadratic form is +∞ on every nonzero vector.
2. The causal-diamond space H_D = P_D H_ω is therefore {0}.
3. Both steps the paper uses to argue otherwise fail:
   - the Wigner-measure argument (5.9)/(C.15b);
   - the order-by-order expansion (5.14)/(C.15e).

Normal ordering removes the infinity, but then the focusing operator Ξ̂_y(v) is no longer positive. The bound (5.21) that ties the clock to the corner area depends on that positivity.

## Setting

We work in the sector the paper itself uses to count graviton states (App. D.3): corner data held at spherical values, and radiative data linearized about the shear-free incoming cone.

- **Corner data.** ϱ_S(y) = b on S, so the area is A = 4πb² and s = A/(4Għ).
- **Radiative data.** One polarization z₁ of the conformal metric, in one real normalized angular harmonic Y: h(v) = ∫_S ε₀ Y z₁. This obeys h(ε) = h(1) = 0 (D.42).
- **Symplectic form.** From (D.42a)–(D.43), with x = v h:
  - Ω(x, y) = κ ∫_ε^1 [(∂x) y − (∂y) x] dv = 2κ ∫_ε^1 (∂x) y dv
  - κ ∝ ϱ_S²/G, so ħ/κ ∝ 1/s.
- **Brackets.** For linear observables F_a = ∫ a x with ∫ a = 0:
  - {F_a, F_b} = (1/2κ) ∫ a(v) B(v) dv, where B(v) = ∫_ε^v b.
  - So observables supported in disjoint intervals commute.
- **Focusing operator.** By (4.22), dropping the nonnegative z₂ term and projecting onto Y (Cauchy–Schwarz on S, an inequality of classical quadratic forms):

  Ξ̄(v) := ∫_S ε₀ Ξ̂_y(v) ≥ Q_v,  where Q_v[x] = ¼ ∫_ε^v t (∂h)² dt = ¼ ∫_ε^v (1/t) φ² dt,  φ = ∂x − x/t.

Every operator below is the Weyl quantization of a classical quadratic form, as in (4.20)–(4.22). Weyl quantization preserves the order of quadratic forms. For a positive form on one canonical pair with {A, B} = c, the Weyl operator is at least ħν/2, where ν = |c|·√det M is the form's symplectic eigenvalue.

## Lemma 1 (local zero-point)

For every v ∈ (ε, 1] and every N, as operators:

  Q_v ≥ Σ_{n=1}^N (ħ ν_n / 2) · 1,  with  ν_n ≥ π/(4κ c_n ℓ) · (1 − O(ℓ/ε)).

Here ℓ = (v − ε)/N, and c_n are the midpoints of the N equal intervals I_n partitioning [ε, v]. Hence

  Q_v ≥ (π/8)(ħ/κ) · N² ln(v/ε)/(v − ε)² · (1 − O(1/N)) for all N,

so Q_v exceeds every multiple of the identity.

**Proof.**
1. **Local pairs.** On I_n write τ = (t − c_n)/ℓ + ½ ∈ [0, 1] and b(τ) = sin²(πτ). Define
   - f = t·b′(τ), g = t·b″(τ)/(2π²);
   - A_n = ∫ f φ, B_n = ∫ g φ.

   Both f/t and g/t integrate to zero, so A_n and B_n are linear observables F_a with ∫ a = 0, supported in I_n. Pairs on different intervals therefore commute.
2. **Their bracket.** Integrating by parts:
   - {A_n, B_n} = (1/2κ)[∫ f′g + ∫ f G/t], with G(t) = ∫^t g/u;
   - ∫ f′g = c_n²π² + O(ℓ²) and |∫ f G/t| ≤ ℓ²/2.
3. **Cauchy–Schwarz.** In the inner product ⟨u, w⟩ = ∫_{I_n} u w/t, A_n = ⟨t f, φ⟩ and B_n = ⟨t g, φ⟩. So, as classical forms,

   ¼ ∫_{I_n} φ²/t ≥ ¼ (A_n, B_n) Gram⁻¹ (A_n, B_n)ᵀ,

   where Gram is the Gram matrix of {t f, t g}. Its determinant satisfies √det Gram = c_n³ ℓ π/2 · (1 + O(ℓ/c_n)).
4. **Symplectic eigenvalue.** This gives ν_n = ¼ |{A_n, B_n}| / √det Gram = π/(4κ c_n ℓ) · (1 + O(ℓ/ε)).
5. **Sum.** Discarding the rest of [ε, v] only lowers the form. Weyl quantization preserves the order, and commuting pairs add their zero points. Summing, Σ_n 1/(c_n ℓ) = N² ln(v/ε)/(v − ε)² · (1 + O(1/N)). ∎

`partition_bound.py` evaluates the exact per-interval ν_n by quadrature. The total grows exactly as N² (N = 10 to 10⁴), with the predicted constant (π/8)·ln(v/ε)/(v − ε)²: 1.116 against 1.1163 at v = 1, ε = 0.1.

## Theorem 1 (the clock form is infinite; H_D = {0})

In this sector the quadratic form (4.27) of the clock satisfies Q[ψ] = +∞ for every ψ ≠ 0. Consequently:
- M, the closure of the form domain in (C.13), is {0};
- Ψ_L = L(L + û*)⁻¹ P_M ψ = 0;
- 1 − ⟨Ω_A|Ψ_L⟩ = 1 for every L;
- ‖P_D Ω_A‖ = 0 for every ℓ(D), contradicting (5.19);
- H_D = P_D H_ω = {0}.

**Proof.**
1. **The clock in this sector.** With ϱ̂_S = b and â = 4πb, (4.26) gives û* = (b/4π) ∫ ε₀ ∫_ε^1 dv e^{Ξ̂_y(v)}.
2. **Cut-off forms.** The form Q_{Ξ*} of (C.11) uses min(Ξ̂_y(v), Ξ*). Since e^x ≥ 1 + x,

   Q_{Ξ*}[ψ] ≥ (b/4π) ∫ dv ∫ ε₀ ⟨ψ| 1 + min(Ξ̂_y(v), Ξ*) |ψ⟩.
3. **Removing the cutoff.** As Ξ* → ∞, monotone convergence takes the right-hand side to (b/4π) ∫ dv [4π‖ψ‖² + ⟨ψ|Ξ̄(v)|ψ⟩].
4. **Conclusion.** By Lemma 1, ⟨ψ|Ξ̄(v)|ψ⟩ ≥ ⟨ψ|Q_v|ψ⟩ = +∞ for every v > ε and every ψ ≠ 0. The forms Q_{Ξ*} increase to Q (C.11)–(C.13), so Q[ψ] = +∞. ∎

The argument uses only the uncertainty principle, through the zero point of each local canonical pair. So it holds in every representation of the canonical commutation relations, for every choice of Gaussian centres and widths (4.7), and for every state.

## Theorem 2 (the paper's two arguments for (5.19) fail)

**(a) Wigner measure, (5.9)/(C.15b).** Let μ be the Wigner function of any pure Gaussian state, with covariance Σ. Then Q_v = +∞ μ-almost surely, so ū* = +∞ μ-almost surely. This contradicts (C.15b)'s claim that μ-almost every configuration has a conformal metric smooth up to S. It also gives ⟨f₀⟩_μ = 1, not → 0, in (5.9).

*Proof.*
1. **The trace is infinite.** Uncertainty, Σ + (iħ/2)J ≥ 0, applied pair by pair to the local pairs of Lemma 1, gives tr(Q_v Σ) ≥ Σ_n ħν_n/2 for every N. So tr(Q_v Σ) = ∞.
2. **Laplace transform.** Let Q^(N) = Σ_n P_n be the finite-rank form built from the N local pairs of Lemma 1, so Q_v ≥ Q^(N). For a Gaussian measure, E[e^{−τQ^(N)}] ≤ Π_k (1 + 2τλ_k)^{−1/2}, where λ_k are the eigenvalues of Σ^{1/2} Q^(N) Σ^{1/2}. A nonzero mean only lowers the left-hand side.
3. **Divergence.** Π_k (1 + 2τλ_k) ≥ 1 + 2τ Σ_k λ_k, and Σ_k λ_k = tr(Q^(N) Σ) ≥ Σ_n ħν_n/2 → ∞ as N → ∞. So E[e^{−τQ_v}] ≤ inf_N E[e^{−τQ^(N)}] = 0 for every τ > 0. Hence P(Q_v < ∞) = lim_{τ→0} E[e^{−τQ_v}] = 0.
4. **The clock.** Fubini over v and e^x ≥ 1 + x give ū* = ∞ almost surely. ∎

**(b) Order-by-order expansion, (5.14)/(C.15e).**
1. At shear-free data, ∂_I ū* = 0 along radiative directions, since focusing starts at second order.
2. So the O(ħ) Wick contraction (5.14) is ½ Σ^{IJ} ∂_I∂_J f₀ = L/(2(L + ū₀)²) · Σ^{IJ} ∂_I∂_J ū*.
3. The Hessian of ū* contains (b/4π) ∫ dv · 2Q_v, and its contraction with Σ is ∞ by step 1 of (a).
4. So the coefficient of ħ¹ in 1 − ⟨Ω_A|Ψ_L⟩ is +∞ for every L, not O(L⁻¹) as stated in (C.15e). ∎

## Proposition 3 (normal ordering breaks positivity)

Subtracting the zero point (normal ordering with respect to a Gaussian Ω) makes Ξ̄(1) finite. But the operator inequality Ξ̂_y(v) ≥ 0 asserted below (4.22) no longer holds. That inequality is the input to (5.21a) and so to the area bound (5.21d).

- **Numerics.** `spectrum.py` takes Ω to be the vacuum of the full-cone form Q₁, which is the most favourable choice for the paper. In a Galerkin basis of N sine modes it computes inf_ψ ⟨ψ|:Q_v:|ψ⟩ = (ħ/2)Σν(Q_v) − tr(Q_v Σ_Ω):

  | N | v = 0.3 | v = 0.55 | v = 0.8 | v = 1 |
  |---|---|---|---|---|
  | 20 | −3.27 | −2.64 | −2.15 | 0 |
  | 40 | −8.08 | −6.04 | −5.04 | 0 |
  | 80 | −16.85 | −13.46 | −11.23 | 0 |
  | 160 | −36.67 | −29.37 | −24.67 | 0 |

- **Reading.** :Q_v: has negative expectation values for every v < 1. Their magnitude grows with the mode cutoff, consistent with Q_v's sharp edge at v inside the interval.
- **Supporting result.** For a two-dimensional massless field, the optimal lower bound on ∫ w :T_uu: is −(ħ/48π) ∫ (w′)²/w (Flanagan 1997). That bound is −∞ for a weight with a sharp edge.

This proposition rests on truncated numerics and the cited bound. It is not proved here at the level of Theorems 1–2.

## Remark: the paper is internally inconsistent

- The clock (4.22) multiplies field operators, so it keeps the zero point. Theorems 1–2 apply.
- The state count (D.45a)–(D.46c) instead counts states with Σ_i n_i ν_i ≤ Q*, which subtracts the zero point. Proposition 3 applies.

Neither choice supports the full chain (4.22) → (5.21) → (5.30) → finite dim H_D. As a cross-check, `spectrum.py` finds the symplectic eigenvalues of Q₁ growing linearly in the radial mode number (ν_k ≈ 3.21k − 0.4 at κ = ħ = 1, ε = 0.1). That reproduces the k/s scaling of (D.41).

## Scope

- **Proven (Lemma 1, Theorems 1–2):** the frozen-corner sector with radiative data linearized about the shear-free incoming cone. These are the same approximations the paper uses for its mode count (App. D.3) and its expansion about ϕ₀ in (5.13)–(5.14).
- **Not covered:**
  - nonlinear radiative couplings (the e^{2z₁} factor and the z₁–z₂ cross terms in Ω);
  - fluctuating corner data;
  - other choices of operator ordering.
- **Assessment of the full theory:** these would have to cancel a positive divergence that is quadratic in the local mode cutoff, uniformly in the state. Nothing in the paper indicates such a cancellation.

## Files

- `partition_bound.py`: exact per-interval symplectic eigenvalues for Lemma 1, for N up to 10⁴.
- `spectrum.py`: Galerkin symplectic spectra of Q_v (zero points growing as N²) and the normal-ordered infima of Proposition 3.

Both run with numpy and scipy.
