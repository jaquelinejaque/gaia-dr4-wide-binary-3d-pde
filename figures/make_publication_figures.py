#!/usr/bin/env python3
"""
Generate publication-grade MNRAS-style figures (vector PDF + 300 DPI PNG)
for the Zenodo Pre-Registration & Technical Note (Martins 2026).
Zero decorative elements; clean serif British academic typography.
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(PKG_ROOT, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Palatino", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.labelsize": 10.5,
    "axes.titlesize": 10.5,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "axes.linewidth": 0.9,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.35,
})

with open(os.path.join(PKG_ROOT, "data", "pde_v0_efe_dwarfs_binaries_results.json"), encoding="utf-8") as fh:
    pde_data = json.load(fh)
with open(os.path.join(PKG_ROOT, "data", "pde_v13_poisson_results.json"), encoding="utf-8") as fh:
    v13_data = json.load(fh)

# ==============================================================================
# FIGURE 1: Gaia DR4 Wide-Binary Pre-Registered Predictions (Boost & Quadrupole)
# ==============================================================================
wb_std = pde_data["gaia_wide_binaries"]["Convencao_Preregistro_aext_1p9e-10"]["bins"]
wb_bar = pde_data["gaia_wide_binaries"]["Convencao_Barionica_MW_gN_1p2e-10"]["bins"]

s_kau = np.array([b["s_kAU"] for b in wb_std])
gam_ref_std = wb_std[0]["gamma_PDE_sph"]
gam_rel_std = np.array([b["gamma_PDE_sph"] / gam_ref_std for b in wb_std])

gam_ref_bar = wb_bar[0]["gamma_PDE_sph"]
gam_rel_bar = np.array([b["gamma_PDE_sph"] / gam_ref_bar for b in wb_bar])

quad_std = np.array([b["quadrupole_sky_perp_minus_par"] for b in wb_std])
quad_bar = np.array([b["quadrupole_sky_perp_minus_par"] for b in wb_bar])

# Reference curves from PREREGISTRO_DR4.md & Poisson AQUAL PDE
s_mond = np.array([1.5, 3.5, 7.5, 15.0, 20.0, 25.0])
gam_mond = np.array([1.00, 1.15, 1.35, 1.39, 1.40, 1.40])
s_v13 = np.array([1.5, 15.0, 25.0])
gam_v13_lo = np.array([1.000, 1.118, 1.119])
gam_v13_hi = np.array([1.000, 1.179, 1.180])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.1), dpi=300)

# Panel (a): Relative acceleration boost gamma_rel(s)
ax1.axhspan(1.05, 1.25, color="#d9e6f2", alpha=0.55, label=r"Frozen v0 confirmatory band $[1.05, 1.25]$")
ax1.axvspan(10.0, 30.0, color="#eeeeee", alpha=0.55)
ax1.axhline(1.00, color="#333333", linestyle="--", linewidth=1.1, label="Newtonian GR baseline (1.000)")
ax1.plot(s_mond, gam_mond, color="#8c2d19", linestyle=":", marker="s", markersize=4.5, linewidth=1.4,
         label=r"Canonical MOND + EFE ($\mu_{\mathrm{RAR}}$, $\sim 1.40$)")
ax1.plot(s_kau, gam_rel_std, color="#0b3c5d", linestyle="-", marker="o", markersize=5.5, linewidth=2.0,
         label=r"v0 3D PDE ($g_{\mathrm{ext}}=1.9\times 10^{-10}$, $\gamma_{\mathrm{rel}}=1.166$)")
ax1.plot(s_kau, gam_rel_bar, color="#328cc1", linestyle="-.", marker="^", markersize=4.8, linewidth=1.4,
         label=r"v0 3D PDE ($g_{N,\mathrm{MW}}=1.2\times 10^{-10}$, $\gamma_{\mathrm{rel}}=1.223$)")
ax1.fill_between(s_v13, gam_v13_lo, gam_v13_hi, color="#2e6f40", alpha=0.22,
                 label="v1.3 Poisson AQUAL PDE (1.118\u20131.179)")

ax1.text(22.2, 1.305, "Primary DR4 bin\n(10\u201330 kAU)", fontsize=8.5, ha="center", color="#444444",
         bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc", alpha=0.85, lw=0.6))
ax1.set_xlabel(r"Projected separation $s$ (kAU)")
ax1.set_ylabel(r"Relative acceleration ratio $\gamma_{\mathrm{rel}}(s) \equiv \gamma_a(s)/\gamma_a(1.5\,\mathrm{kAU})$")
ax1.set_title("(a) Pre-registered radial acceleration profile")
ax1.set_xlim(0.5, 27.5)
ax1.set_ylim(0.96, 1.50)
ax1.grid(True)
ax1.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#cccccc")

# Panel (b): Sky-projected quadrupole anisotropy
ax2.axvspan(10.0, 30.0, color="#eeeeee", alpha=0.55)
ax2.axhline(0.0, color="#333333", linestyle="--", linewidth=1.0, label="Isotropic baseline (0.000)")
ax2.axhspan(0.03, 0.05, color="#f4cccc", alpha=0.45, label="Superseded mental-algebra sign (+0.03 to +0.05)")
ax2.plot(s_kau, quad_std, color="#0b3c5d", linestyle="-", marker="o", markersize=5.5, linewidth=2.0,
         label=r"Exact 3D PDE ($g_{\mathrm{ext}}=1.9\times 10^{-10}$: $-0.024$ in 10–30 kAU)")
ax2.plot(s_kau, quad_bar, color="#328cc1", linestyle="-.", marker="^", markersize=4.8, linewidth=1.4,
         label=r"Exact 3D PDE ($g_{N,\mathrm{MW}}=1.2\times 10^{-10}$: $-0.022$ in 10–30 kAU)")

ax2.annotate("Corrected pre-DR4 prediction:\n" + r"$\gamma_{\perp}-\gamma_{\parallel} = -0.024$ ($\gamma_{\parallel} > \gamma_{\perp}$)",
             xy=(20.0, -0.0259), xytext=(8.5, -0.043),
             arrowprops=dict(arrowstyle="->", color="#0b3c5d", lw=1.0),
             fontsize=8.5, color="#0b3c5d", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#0b3c5d", lw=0.7))

ax2.set_xlabel(r"Projected separation $s$ (kAU)")
ax2.set_ylabel(r"Sky-projected quadrupole $\Delta_{\gamma} \equiv \gamma_{\perp} - \gamma_{\parallel}$")
ax2.set_title("(b) Quadrupole anisotropy relative to Galactic Centre")
ax2.set_xlim(0.5, 27.5)
ax2.set_ylim(-0.052, 0.058)
ax2.grid(True)
ax2.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#cccccc")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1_gaia_dr4_wide_binaries.pdf"))
fig.savefig(os.path.join(OUT_DIR, "fig1_gaia_dr4_wide_binaries.png"), dpi=300)
plt.close(fig)

# ==============================================================================
# FIGURE 2: Milky Way Classical dSphs & SPARC HI-Rich Dwarf vs Spiral Tension
# ==============================================================================
dsphs = pde_data["milky_way_8_dsphs"]
names = [d["galaxy"] + ("*" if d["tidal_system"] else "") for d in dsphs]
gam_lin = [d["pop_result"]["gamma_lin"] for d in dsphs]
gam_v0 = [d["pop_result"]["gamma_pde"] for d in dsphs]
gam_v0_hi = [d["gamma_pde_range_ML_1p2_to_2p0"][1] for d in dsphs]
# v1.3 PDE results loaded directly from verified data/pde_v13_poisson_results.json
v13_raw = {d["galaxy"]: d["gamma_sigma"] for d in v13_data["dsphs"]}
gam_v13 = [v13_raw[d["galaxy"]] for d in dsphs]

fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(10.4, 4.2), dpi=300)

x_idx = np.arange(len(names))
w = 0.26
bx1.axhspan(0.90, 1.10, color="#d9ead3", alpha=0.45, label=r"Observational consistency ($\pm 10\%$)")
bx1.axhline(1.00, color="#333333", linestyle="--", linewidth=1.0)
bx1.bar(x_idx - w, gam_lin, width=w, color="#b7b7b7", edgecolor="#555555", linewidth=0.6,
        label="v0 linearised approx. (invalid for $\\eta_{1/2}>1$)")
bx1.bar(x_idx, gam_v0, width=w, color="#0b3c5d", edgecolor="#062133", linewidth=0.6,
        label=r"v0 exact nonlinear PDE ($a_c=1.9\times 10^{-11}$)")
bx1.errorbar(x_idx, gam_v0, yerr=[np.zeros(len(names)), np.maximum(0, np.array(gam_v0_hi) - np.array(gam_v0))],
             fmt="none", ecolor="#0b3c5d", capsize=2.5, elinewidth=0.9)
bx1.bar(x_idx + w, gam_v13, width=w, color="#2e6f40", edgecolor="#193d23", linewidth=0.6,
        label=r"v1.3 Poisson AQUAL PDE ($a_0=1.2\times 10^{-10}$)")

bx1.set_xticks(x_idx)
bx1.set_xticklabels(names, rotation=28, ha="right")
bx1.set_ylabel(r"Predicted-to-observed dispersion $\gamma_{\sigma} \equiv \sigma_{\mathrm{PDE}}/\sigma_{\mathrm{obs}}$")
bx1.set_title("(a) Classical Milky Way dSphs (*tidally disrupted)")
bx1.set_ylim(0.0, 1.38)
bx1.grid(True, axis="y")
bx1.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#cccccc")

# Panel (b): SPARC RMS scatter comparison across galaxy regimes
categories = ["Canonical SPARC\n(144 galaxies)", "Canonical Dwarfs\n(40 gal., $V_{\\max}<80$)", "HI-Rich Dwarfs\n($f_{\\mathrm{gas}}\\geq 85\\%$, 8 gal.)", "Ultra-HI Dwarfs\n($f_{\\mathrm{gas}}\\geq 90\\%$, 4 gal.)"]
rms_rar_noefe = [0.1325, 0.1899, 0.1963, 0.1460]
rms_rar_efe   = [0.1391, 0.1702, 0.1845, 0.1545]
rms_v13_efe   = [0.1556, 0.1700, 0.1762, 0.1538]

cx = np.arange(len(categories))
bx2.bar(cx - w, rms_rar_noefe, width=w, color="#8c2d19", alpha=0.82, edgecolor="#4a170d", linewidth=0.6,
        label="Standard RAR (no EFE)")
bx2.bar(cx, rms_rar_efe, width=w, color="#d9822b", alpha=0.88, edgecolor="#6e3f11", linewidth=0.6,
        label="RAR + EFE control (Chae+2021)")
bx2.bar(cx + w, rms_v13_efe, width=w, color="#2e6f40", alpha=0.90, edgecolor="#193d23", linewidth=0.6,
        label=r"v1.3 Poisson $\mu(x)=1-e^{-x}$ + EFE")

bx2.set_xticks(cx)
bx2.set_xticklabels(categories, fontsize=8.4)
bx2.set_ylabel(r"Residual scatter $\mathrm{RMS}(\log_{10}(g_{\mathrm{obs}}/g_{\mathrm{pred}}))$ [dex]")
bx2.set_title("(b) SPARC rotation-curve trade-off (dwarfs vs spirals)")
bx2.set_ylim(0.10, 0.245)
bx2.grid(True, axis="y")
bx2.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, edgecolor="#cccccc")

fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "fig2_dsphs_and_sparc_tension.pdf"))
fig2.savefig(os.path.join(OUT_DIR, "fig2_dsphs_and_sparc_tension.png"), dpi=300)
plt.close(fig2)

print("Saved figures successfully in:", OUT_DIR)
