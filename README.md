# Zero-point focusing makes the null-cone clock of arXiv:2611.34924 infinite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24577548.svg)](https://doi.org/10.5281/zenodo.24577548)

V. Chandrasekaran, *An algebraic quantization of gravity from the null cone* (arXiv:2611.34924), derives a finite-dimensional Hilbert space for a causal diamond in 4D gravity. The derivation uses a gravitational clock: an affine parameter, averaged over the null generators, that is quantized as a product of shear operators ((4.22)–(4.27)). The causal diamond is the subspace where the clock reads at most ℓ(D).

This repository shows that the derivation fails. Each mode on the null cone carries a zero-point fluctuation, and the cone has infinitely many modes. Splitting the cone into N intervals and taking one canonical pair in each gives, for every N,

    Q_v ≥ (π/√10)(ħ/κ) · N² ln(v/ε)/(v−ε)²,

where Q_v is the shear part of the focusing operator up to affine position v, and ħ/κ ∝ 1/s with s = A/(4Għ).

The results below hold in the setting the paper uses for its own state counting (App. D.3) and for its expansion about the spherical cone: corner data fixed, radiation linearized, every radial mode kept.

- **Theorem 1.** The clock's quadratic form is +∞ on every nonzero vector. The causal-diamond space is therefore {0}.
- **Theorem 2.** Both of the paper's arguments that the diamond has states fail:
  - (C.15b) claims Wigner-typical configurations are smooth; they are almost surely not.
  - (C.15e) claims the O(ħ) correction (5.14) is O(1/L); it is +∞.
- **Proposition 3.** Normal ordering removes the infinity, but the partial focusing then takes negative values for v < 1. The paper asserts that focusing is positive (below (4.22)) and uses this to bound the area.

The full statements, proofs and scope are in [PROOF.md](PROOF.md).

## Running the checks

```
pip install numpy scipy
python partition_bound.py   # Lemma 1: exact per-interval ground energies, N up to 10^4
python spectrum.py          # Galerkin spectra of the focusing form; pair coefficient K_11
```

## Citation

V. S. Baraiya, *Zero-point focusing makes the null-cone clock of arXiv:2611.34924 infinite*, Zenodo, doi:[10.5281/zenodo.24577548](https://doi.org/10.5281/zenodo.24577548) (all versions; v1.0: doi:[10.5281/zenodo.24577549](https://doi.org/10.5281/zenodo.24577549)). See [CITATION.cff](CITATION.cff).

## License

Code (`*.py`): MIT License. Text (`README.md`, `PROOF.md`): CC BY 4.0. See [LICENSE](LICENSE).
