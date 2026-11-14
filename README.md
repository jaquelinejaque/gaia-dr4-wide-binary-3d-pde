# Zenodo Pre-Registration Deposit Package (Pre-Gaia DR4)

**Title:** Nonlinear Field-Equation Predictions of an Additive Cubic-Gradient Scalar Model (`v0`) for Gaia DR4 Wide Binaries and Tests on Dwarf Galaxies  
**Author:** Jaqueline Martins (Independent Researcher, United Kingdom)  
**Deposit Date (Frozen Pre-DR4):** 14 November 2026  
**Target Data Release:** ESA *Gaia* Data Release 4 (expected December 2026)  
**Licences:** CC-BY-4.0 (Manuscript, Tables, Figures) · MIT Licence (Numerical Solver & Audit Scripts)

---

## Package Structure

- `Martins_2026_Gaia_DR4_PreRegistration_Technical_Note.pdf` — Typeset A4 Pre-Registration & Technical Note for archival Zenodo deposit.
- `Martins_2026_Gaia_DR4_PreRegistration_Technical_Note.md` — Complete Markdown/LaTeX source of the manuscript.
- `SHA256SUMS.txt` — Cryptographic SHA-256 manifest locking all manuscript, solver, data, and figure files prior to *Gaia* DR4.
- `code/`
  - `pde_solver_v0_isolated_validation.py` — Standalone cylindrical Finite-Volume Newton–Raphson validation suite (Newtonian Poisson check + $200\times 400 \to 400\times 800$ Gauss spherical shell convergence).
  - `pde_v0_efe_dwarfs_binaries.py` — Standalone 3D axisymmetric nonlinear PDE solver (`ub_z > 0`) for the 8 classical Milky Way dwarf spheroidals, *Gaia* wide binaries ($1\text{--}30\text{ kAU}$), explicit two-body check ($0.85 + 0.65\,M_\odot$), and Monte Carlo sky-projected quadrupole ($\Delta_\gamma = -0.024$).
  - `pde_v13_poisson_aqual_solver.py` — Standalone 3D axisymmetric nonlinear AQUAL PDE solver for the `v1.3` exponential Poisson interpolating extension ($\mu(x) = 1 - e^{-x}$, $a_0 = 1.20 \times 10^{-10}\text{ m s}^{-2}$).
  - `audit_v13_sparc_dwarfs_efe.py` — Controlled 175-galaxy SPARC rotation-curve and Chae et al. (2021) External Field Effect audit (`RAR` vs `RAR + EFE` vs `v1.3 Poisson + EFE`).
- `data/`
  - `pde_v0_validation_results.json` — Machine-precision numerical verification outputs (Table 1).
  - `pde_v0_efe_dwarfs_binaries_results.json` — Exact `v0` nonlinear PDE outputs for wide binaries and Milky Way dwarf spheroidals (Tables 2 & 3).
  - `pde_v13_poisson_results.json` — Exact `v1.3` nonlinear PDE outputs for wide binaries and Milky Way dwarf spheroidals (Tables 2 & 3).
  - `audit_v13_results.json` & `gas_rich_dwarfs_v13_results.json` — Controlled SPARC subsample statistics (Table 4).
- `figures/`
  - `fig1_gaia_dr4_wide_binaries.pdf` / `.png` — Figure 1: Wide-binary relative acceleration profile $\gamma_{\mathrm{rel}}(s)$ and sky-projected quadrupole $\Delta_\gamma(s)$.
  - `fig2_dsphs_and_sparc_tension.pdf` / `.png` — Figure 2: Classical Milky Way dSph dispersion ratios $\gamma_\sigma$ and SPARC rotation-curve RMS scatter comparison.
  - `make_publication_figures.py` — Reproducible figure generation script (`matplotlib` STIX serif).

---

## Reproducing All Numerical Results From Scratch

Requirements: Python $\ge 3.10$ with `numpy`, `scipy`, and `matplotlib`.

```bash
python3 code/pde_solver_v0_isolated_validation.py
python3 code/pde_v0_efe_dwarfs_binaries.py
python3 code/pde_v13_poisson_aqual_solver.py
python3 figures/make_publication_figures.py
shasum -a 256 -c SHA256SUMS.txt
```

---

## Citation & Official Archive (CERN / Zenodo)

- **DOI:** [10.5281/zenodo.24210580](https://doi.org/10.5281/zenodo.24210580)
- **ORCID:** [0009-0009-7732-5915](https://orcid.org/0009-0009-7732-5915)

## Research Compute & Open-Science Patronage

Independent open-science reproduction and 3D nonlinear PDE compute for *Gaia* DR4 wide-binary testing is supported by community patronage:
- **Direct Research Support (Stripe UK — GBP / USD / EUR):** [https://donate.stripe.com/9B64gzeaI8Dr6By305dwc3G](https://donate.stripe.com/9B64gzeaI8Dr6By305dwc3G)
