# Zero-point focusing makes the null-cone clock of arXiv:2611.34924 infinite

**Paper:** V. Chandrasekaran, *An algebraic quantization of gravity from the null cone*, arXiv:2611.34924v1. Equation numbers below refer to that paper.

## Result

The gravitational clock (4.22)–(4.27) is built from products of shear operators, which keep their zero-point fluctuations. We work in the model the paper uses for its graviton count and for its expansion about the spherical cone:

- corner data held at spherical values;
- radiative data linearized about the shear-free incoming cone;
- every radial mode kept, as in the infinite product (4.7).

In this model:

1. **Lemma 1.** The focusing operator exceeds every constant. The zero-point focusing grows without bound as the cone is resolved more finely.
2. **Theorem 1.** The clock's quadratic form (4.27) is +∞ on every nonzero vector. The model's causal-diamond space P_D H is therefore {0}.
3. **Theorem 2.** Both steps the paper uses for (5.19) fail in the model:
   - the Wigner-measure argument (5.9)/(C.15b);
   - the order-by-order expansion (5.14)/(C.15e).
4. **Proposition 3.** If the zero point is removed by normal ordering, the focusing operator is no longer positive for v < 1. That positivity is asserted below (4.22).

The paper evaluates (5.8) and (C.15)–(C.16) "at fixed resolution". At a fixed radial resolution N, the zero point is finite and scales as N²/s. So the paper's estimates hold resolution by resolution but not uniformly in N. The limit N → ∞, which the infinite product (4.7) requires, is the case covered by Theorems 1–2.

## Setting

- **Corner data.** ϱ_S(y) = b on S, so the area is A = 4πb² and s = A/(4Għ).
- **Radiative data.** One polarization z₁ of the conformal metric, in one real normalized angular harmonic Y: h(v) = ∫_S ε₀ Y z₁. This obeys h(ε) = h(1) = 0 (D.42).
- **Symplectic form.** From (D.42a)–(D.43), with x = v h:
  - Ω(x, y) = κ ∫_ε^1 [(∂x) y − (∂y) x] dv = 2κ ∫_ε^1 (∂x) y dv
  - κ ∝ ϱ_S²/G, so ħ/κ ∝ 1/s.
- **Brackets.** Linear observables F_a = ∫ a x are admissible when ∫ a = 0, which is needed for a Hamiltonian vector field to exist. For such observables:
  - {F_a, F_b} = (1/2κ) ∫ a(v) B(v) dv, with B(v) = ∫_ε^v b. B vanishes outside supp b.
- **Focusing.** Write φ = ∂x − x/t = t∂h. Dropping the nonnegative z₂ term of (4.22) and projecting onto Y (Bessel's inequality on S), the sphere-integrated focusing satisfies

  Ξ̄(v) := ∫_S ε₀ Ξ̂_y(v) ≥ Q_v,  where Q_v = ¼ ∫_ε^v φ²/t dt.

  The ¼ comes from the ⅛ in (3.26c), because (∂q̄)·(∂q̄) = 2 g_{αβ} ∂z^α ∂z^β.

**Operator meaning.** The paper applies Weyl ordering "wherever the smeared products are defined" (4.20). It asserts Ξ̂ ≥ 0 below (4.22) because Ξ̂ is a sum of squares. Consistently with that, we read the operator form of Q_v as a sum of squares:

  Q̂_v[ψ] := sup over finite orthonormal families {e_k} ⊂ L²([ε, v], dt/t) of ¼ Σ_k ‖φ̂(e_k) ψ‖²,

where φ̂(e) = ∫ e φ̂ dt/t is the smeared field. This definition:
- takes values in [0, ∞];
- does not depend on a basis;
- needs no infinite-dimensional Weyl calculus.

## Lemma 1 (local zero point)

For every v ∈ (ε, 1], every N, and every vector ψ in a regular representation (one where the C¹, compactly supported smearings below define self-adjoint field operators):

  Q̂_v[ψ] ≥ ‖ψ‖² Σ_{n=1}^N E_n,  with E_n = (α/8)(ħ/κ)/(c_n ℓ) · (1 + O(ℓ/ε)),  α = 4√10 π/5.

Here ℓ = (v − ε)/N and c_n are the midpoints of the N equal intervals partitioning [ε, v]. Hence

  Q̂_v[ψ] ≥ ‖ψ‖² · (π/√10)(ħ/κ) · N² ln(v/ε)/(v − ε)² · (1 + O(1/N))

for every N, so Q̂_v[ψ] = +∞ for every ψ ≠ 0.

**Proof.**
1. **Smearings.** On the interval I_n write τ = (t − c_n)/ℓ + ½ ∈ [0, 1], and take
   - u₁ = sin²(πτ) sin 2πτ,
   - u₂ = sin²(πτ)(cos 2πτ + ½).

   Both are C¹, vanish with their derivatives at τ = 0 and 1, and have zero mean. Set f = t u₁, g = t u₂, A_n = ∫ f φ and B_n = ∫ g φ. Integrating by parts, A_n = F_a with a = −f′ − f/t, which is continuous with ∫ a = −∫ f/t = −ℓ∫u₁ = 0. The same holds for B_n. So A_n and B_n are admissible observables supported in I_n.
2. **Bracket.** {A_n, B_n} = (1/2κ)[∫ f′g + ∫ f G/t], where G(t) = ∫^t g/u = ℓ U₂(τ) and U₂ = ∫₀^τ u₂ vanishes at both ends. As ℓ → 0 this gives {A_n, B_n} = (c_n²/2κ) I₁ (1 + O(ℓ/c_n)), with I₁ = ∫ u₁′ u₂ = π/4.
3. **Bessel's inequality.** In ⟨p, q⟩ = ∫_{I_n} p q/t, A_n = ⟨t f, φ⟩ and B_n = ⟨t g, φ⟩. Choosing {e_k} as an orthonormal basis of span{t f, t g} in the definition of Q̂_v gives

   ¼ Σ_k ‖φ̂(e_k)ψ‖² = ⟨ψ| ¼ (Â_n, B̂_n) Gram⁻¹ (Â_n, B̂_n)ᵀ |ψ⟩.

   Here Gram is the Gram matrix of {t f, t g}, with √det Gram = c_n³ ℓ √D (1 + O(ℓ/c_n)). The constant D = ∫u₁² ∫u₂² − (∫u₁u₂)² = (5/32)(1/16).
4. **Ground energy.** The operator zᵀMz with M = ¼ Gram⁻¹ and [Â, B̂] = iħc is at least ħ|c|√det M = ħ|c|/(4√det Gram). This gives E_n = (α/8)(ħ/κ)/(c_n ℓ) · (1 + O(ℓ/ε)), with α = I₁/√D = 4√10 π/5.
5. **Adding intervals.** Families on different intervals together form a single orthonormal family, so the bounds add. The pairs do not need to commute, since each term is bounded separately. Finally, Σ_n 1/(c_n ℓ) = N² ln(v/ε)/(v − ε)² · (1 + O(1/N)). ∎

`partition_bound.py` evaluates the exact E_n by quadrature. Σ_n E_n/N² matches (α/8) ln(v/ε)/(v − ε)² to six digits at N = 10⁴ (2.82410 at v = 1, ε = 0.1, κ = ħ = 1). The closed form α = 4√10 π/5 follows from I₁ = π/4, ∫u₁² = 5/32, ∫u₂² = 1/16 and ∫u₁u₂ = 0.

## Theorem 1 (the model clock is infinite)

Assume the clock satisfies, as quadratic forms,

  û* ≥ (b/4π) ∫_ε^1 dv [4π + Ξ̄(v)].

This holds in the model for the clock (4.26) under any definition of e^{Ξ̂_y(v)} that obeys e^x ≥ 1 + x, through the cut-off forms min(Ξ̂_y(v), Ξ*) of (C.11) and monotone convergence. Note that (4.26) needs Ξ̂_y at a single point y, which is not itself an operator, so this assumption is what any definition has to supply.

Under that assumption, the clock's quadratic form (4.27) is +∞ on every ψ ≠ 0. Consequently:
- M, the closure of the form domain in (C.13), is {0};
- Ψ_L = L(L + û*)⁻¹ P_M ψ = 0, so 1 − ⟨Ω_A|Ψ_L⟩ = 1 for every L;
- ‖P_D Ω_A‖ = 0 for every ℓ(D);
- the model's causal-diamond space is P_D H = {0}.

**Proof.** In the model, ϱ̂_S = b and â = 4πb. Lemma 1 gives ⟨ψ|Ξ̄(v)|ψ⟩ ≥ Q̂_v[ψ] = +∞ for every v > ε. Integrating over v proves the claim. The rest follows from (C.13), whose footnote extends the projections by zero on M⊥. ∎

## Theorem 2 (the paper's arguments for (5.19) fail in the model)

**(a) Wigner measure, (5.9)/(C.15b).** Let μ be the Wigner function of a pure Gaussian state in which the smeared fields A_n and B_n have finite variance. This holds for Hadamard-type states and for the widths (Δq_k)² = 1/(2sk) used in (4.32). Then Q_v = +∞ μ-almost surely, and so ū* = +∞ μ-almost surely. This contradicts (C.15b)'s claim that μ-almost every configuration is smooth up to S. In (5.9) it gives ⟨f₀⟩_μ = 1 for every L, not → 0.

*Proof.*
1. **Lower bound.** Let Q^(N) = Σ_n ¼ (A_n, B_n) Gram_n⁻¹ (A_n, B_n)ᵀ. Then Q_v ≥ Q^(N) pointwise, by step 3 of Lemma 1.
2. **Uncertainty.** The covariance Σ of a quantum state satisfies Σ + (iħ/2)J ≥ 0. Applied pair by pair, this gives tr(Q^(N) Σ) ≥ Σ_n E_n, which tends to ∞ as N → ∞.
3. **Laplace transform.** For a Gaussian measure, E[e^{−τQ^(N)}] ≤ Π_k (1 + 2τλ_k)^{−1/2}, where λ_k are the eigenvalues of Σ^{1/2} Q^(N) Σ^{1/2}. A nonzero mean only lowers the left-hand side. Since Π_k (1 + 2τλ_k) ≥ 1 + 2τ tr(Q^(N) Σ), we get E[e^{−τQ_v}] ≤ inf_N E[e^{−τQ^(N)}] = 0 for every τ > 0.
4. **Conclusion.** P(Q_v < ∞) = lim_{τ→0} E[e^{−τQ_v}] = 0. Fubini over v and e^x ≥ 1 + x give ū* = +∞ almost surely. ∎

**(b) Order-by-order expansion, (5.14)/(C.15e).**
1. At shear-free data, ∂_I ū* = 0 along radiative directions, since focusing starts at second order.
2. Along corner directions the first derivative also vanishes, because â is maximal at the round sphere, and the corner Hessian is ≥ 0.
3. So the O(ħ) Wick term (5.14) is at least L/(2(L + ū₀)²) · (b/4π) ∫ dv · tr(Σ · 2Q_v), which is +∞ by step 2 of (a).
4. The coefficient of ħ¹ in 1 − ⟨Ω_A|Ψ_L⟩ is therefore +∞ for every L, not O(L⁻¹) as stated in (C.15e). ∎

## Fixed resolution

With N radial modes kept, the zero point is finite.
- Galerkin evaluation: the low symplectic eigenvalues of Q₁ converge to 2.801, 6.011, 9.214, 12.411, … (κ = ħ = 1, ε = 0.1). They grow linearly, reproducing the k/s scaling of (D.41).
- The truncated zero point grows as about 0.3 N² (ħ/κ), i.e. of order N²/s.

So the clock shift from vacuum focusing is negligible while N ≪ √s, and grows without bound beyond that. The paper's fixed-resolution estimates, (5.8) and (C.15)–(C.16), do not control this limit. Its Hilbert space (4.7) is an infinite product over all modes, which is the case covered by Theorems 1–2.

## Proposition 3 (normal ordering loses positivity)

Normal ordering removes the zero point. Take normal ordering with respect to Ω, the vacuum of the full-cone form Q₁. Then:
- ⟨Ω|:Q̂_v:|Ω⟩ = 0, and :Q̂_v: ≥ 0 only if :Q̂_v:Ω = 0.
- If the coefficient K₁₁ of a₁†a₁† in :Q̂_v: is nonzero, take ψ = Ω + η a₁†a₁†Ω with small η of the right sign. This gives ⟨ψ|:Q̂_v:|ψ⟩ < 0.

`spectrum.py` computes |K₁₁| in a Galerkin sine basis. It is converged in the cutoff and nonzero for every v < 1, and vanishes at v = 1:

| v | 0.3 | 0.55 | 0.8 | 1 |
|---|---|---|---|---|
| \|K₁₁\| (N = 40 to 240) | 0.10277 | 0.00352 | 0.32242 | 0 |

So, under normal ordering, the operator inequality Ξ̂_y(v) ≥ 0 asserted below (4.22) fails for v < 1.

This does not by itself break (5.21a). That bound needs only a lower bound on ∫ dv e^{Ξ̂}, and a bound :Ξ̂: ≥ −C would weaken it by a factor e^{−C} without destroying it. Whether such a C exists, uniformly in the state, is not settled here.

Appendix D.3 counts graviton states by Σ_i n_i ν_i ≤ Q*, i.e. with normal-ordered excitation energies. That counting applies to a different functional, Q[h] in (D.40), bounded through ∫Θ(n) in (D.38). The paper does not say which ordering the clock (4.22) is meant to carry.

## Scope

- **Proven:** Lemma 1 and Theorems 1–2, in the model above:
  - z₁ polarization, one angular harmonic;
  - radiative data linearized about the shear-free cone;
  - corner frozen at spherical values;
  - all radial modes kept;
  - the clock read as a sum of squares.
- **Not covered:**
  - fluctuating corner data. Through δK∧α in (3.42b) the local brackets then scale as 1/K̂(y);
  - the nonlinear couplings (the e^{2z₁} factor and the z₁–z₂ cross terms in Ω);
  - other definitions of e^{Ξ̂_y} at a point.
- **Not proven:** Proposition 3 is a converged numerical statement about one coefficient, plus an elementary argument from it.

## Files

- `partition_bound.py`: exact per-interval ground energies for Lemma 1 (N up to 10⁴), against the closed-form constant.
- `spectrum.py`: Galerkin symplectic spectra of Q₁, the truncated zero points, and the pair-creation coefficient K₁₁ for Proposition 3.

Both run with numpy and scipy.
