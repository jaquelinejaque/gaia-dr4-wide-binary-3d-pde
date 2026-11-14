#!/usr/bin/env python3
"""
Controlled SPARC Rotation-Curve & External Field Effect Audit (175 Galaxies)
Author: Jaqueline Martins (2026)

Evaluates at fixed canonical stellar mass-to-light ratios (Upsilon_disk = 0.5, Upsilon_bulge = 0.7 M_sun/L_sun)
and canonical acceleration scale a0 = 1.20e-10 m s^-2:
  1. Standard RAR (McGaugh et al. 2016): g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0)))
  2. Invalid explicit baryonic exponential: g_obs = g_bar / (1 - exp(-g_bar/a0)) [asymptotic sanity check]
  3. Implicit Milgrom (1999) Poisson interpolation (v1.3 isolated): (1 - exp(-g_obs/a0)) * g_obs = g_bar
  4. v1.3 Poisson + EFE and RAR + EFE control using environmental external fields from Chae et al. (2021)
  5. Redshift scaling sensitivity check a0(z) across SPARC distances D (Mpc) -> z = H0 * D / c
"""

import glob
import os
import re
import json
import numpy as np
from scipy.optimize import brentq, minimize_scalar

PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.environ.get("SPARC_DATA_DIR", os.path.join(PKG_ROOT, "sparc_data"))
KMS2_PER_KPC_TO_SI = 1e6 / 3.08567758e19
A0_CANON = 1.20e-10  # m/s^2
H0_KMS_MPC = 70.0
C_KMS = 299792.458

# 1. Carregar metadados oficiais SPARC_Lelli2016c.mrt
meta = {}
with open(f"{BASE}/data/sparc_raw/SPARC_Lelli2016c.mrt", encoding="utf-8") as fh:
    lines = fh.readlines()
i0 = max(i for i, l in enumerate(lines) if l.startswith("---")) + 1
for l in lines[i0:]:
    t = l.split()
    if len(t) >= 18:
        try:
            meta[t[0]] = {
                "T": int(t[1]),
                "D_Mpc": float(t[2]),
                "Inc": float(t[5]),
                "L36": float(t[7]),
                "Vflat": float(t[15]),
                "Q": int(t[17]),
            }
        except ValueError:
            pass

# 2. Carregar tabela ambiental real Chae et al. (2020/2021) (log10 e_N = log10(g_N,ext / a0))
chae_env = {}
env_path = f"{BASE}/out/sparc/chae2021_env.csv"
if os.path.exists(env_path):
    with open(env_path, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                gname = parts[0].strip()
                # extrair numero float com sinal (tratando minus unicode −)
                s_max = parts[1].replace("−", "-")
                s_min = parts[2].replace("−", "-")
                m1 = re.search(r"[-+]?\d+\.\d+", s_max)
                m2 = re.search(r"[-+]?\d+\.\d+", s_min)
                if m1 and m2:
                    log_emax = float(m1.group(0))
                    log_emin = float(m2.group(0))
                    # media geometrica entre clustering max e min (padrao Chae+2021)
                    log_emed = 0.5 * (log_emax + log_emin)
                    chae_env[gname] = {
                        "eN_max": 10.0 ** log_emax,
                        "eN_min": 10.0 ** log_emin,
                        "eN_med": 10.0 ** log_emed,
                    }

print(f"Metadados SPARC carregados: {len(meta)} galaxias | Chae+2021 EFE carregado: {len(chae_env)} galaxias")

# Mediana ambiental para galaxias sem entrada individual em Chae+2021
default_eN_med = float(np.median([v["eN_med"] for v in chae_env.values()]))
default_eN_max = float(np.median([v["eN_max"] for v in chae_env.values()]))
print(f"Mediana Chae+2021: eN_med = {default_eN_med:.5f} (log10={np.log10(default_eN_med):.3f}), eN_max = {default_eN_max:.5f}")


# --- Modelos de Aceleracao ---
def pred_rar_std(gbar, a0):
    """McGaugh 2016 RAR padrao sem EFE."""
    y = np.maximum(gbar / a0, 1e-30)
    return gbar / (1.0 - np.exp(-np.sqrt(y)))


def pred_explicit_baryonic_exp_invalid(gbar, a0):
    """Explicit baryonic exponential (invalid in the deep regime as g_obs -> a0): g_obs = g_bar / (1 - exp(-g_bar/a0))."""
    y = np.maximum(gbar / a0, 1e-30)
    return gbar / (1.0 - np.exp(-y))


def pred_poisson_mu_isolated(gbar, a0):
    """Poisson/Boltzmann correto em x = g_obs/a0 (sem EFE): (1 - exp(-g_obs/a0)) * g_obs = g_bar."""
    out = np.empty_like(gbar)
    for i, gb in enumerate(gbar):
        def f(log_g):
            g = 10.0 ** log_g
            return (1.0 - np.exp(-g / a0)) * g - gb
        out[i] = 10.0 ** brentq(f, -16.0, -5.0)
    return out


def pred_poisson_mu_efe_sub(gbar, a0, eN):
    """
    Poisson/Boltzmann v1.3 + EFE pela regra de subtracao de campo total menos ambiente (Jaqueline / AQUAL radial):
      1) Ambiente resolve: (1 - exp(-g_ext/a0)) * g_ext = g_N,ext = eN * a0
      2) Total resolve:    (1 - exp(-g_tot/a0)) * g_tot = g_bar + g_N,ext
      3) Subtracao fisica: g_obs = g_tot - g_ext
      (Para media esferica 3D de orientacoes com angulo medio, AQUAL adiciona fator ~1.1 no limiar EFE;
       testamos a subtracao pura sem parametros livres extras!)
    """
    gN_ext = eN * a0
    # Resolver g_ext uma unica vez por galaxia
    def f_ext(log_g):
        g = 10.0 ** log_g
        return (1.0 - np.exp(-g / a0)) * g - gN_ext
    g_ext = 10.0 ** brentq(f_ext, -16.0, -5.0)

    out = np.empty_like(gbar)
    for i, gb in enumerate(gbar):
        def f_tot(log_g):
            g = 10.0 ** log_g
            return (1.0 - np.exp(-g / a0)) * g - (gb + gN_ext)
        g_tot = 10.0 ** brentq(f_tot, -16.0, -5.0)
        out[i] = max(g_tot - g_ext, 1e-20)
    return out


def pred_rar_efe_chae(gbar, a0, eN):
    """
    RAR + EFE analitico de Chae et al. (2020, ApJ 904, 51, Eq. 6 / AQUAL EFE):
    Substitui nu(y) por nu_EFE(y, eN) via subtracao / Chae+2020:
    Usamos a mesma regra exata de subtracao na RAR para comparacao 100% justa:
      g_tot = rar(gbar + gN_ext), g_ext = rar(gN_ext) -> g_obs = g_tot - g_ext
    """
    gN_ext = eN * a0
    g_ext = pred_rar_std(np.array([gN_ext]), a0)[0]
    g_tot = pred_rar_std(gbar + gN_ext, a0)
    return np.maximum(g_tot - g_ext, 1e-20)


# --- Funcoes Cosmologicas para a0(z) ---
def a0_hubble(z, a0=A0_CANON):
    return a0 * np.sqrt(0.315 * (1.0 + z) ** 3 + 0.685)


def a0_desi_cpl(z, a0=A0_CANON, w0=-0.75, wa=-0.86):
    a = 1.0 / (1.0 + z)
    rho_ratio = (a ** (-3.0 * (1.0 + w0 + wa))) * np.exp(-3.0 * wa * (1.0 - a))
    return a0 * np.sqrt(rho_ratio)


# Carregar todos os 175 arquivos rotmod
files = sorted(glob.glob(f"{BASE}/data/sparc_raw/Rotmod_LTG/*_rotmod.dat"))
gal_records = []

for fp in files:
    name = os.path.basename(fp).replace("_rotmod.dat", "")
    mi = meta.get(name)
    if mi is None:
        continue
    d = np.loadtxt(fp)
    if d.ndim == 1:
        d = d[None, :]
    R, Vobs, eV, Vgas, Vdisk, Vbul = d[:, :6].T
    vmax = float(np.max(Vobs))
    z_gal = max((H0_KMS_MPC * mi["D_Mpc"]) / C_KMS, 0.0)
    env_info = chae_env.get(name, {"eN_med": default_eN_med, "eN_max": default_eN_max, "has_chae": False})
    if name in chae_env:
        env_info["has_chae"] = True

    gal_records.append({
        "name": name,
        "meta": mi,
        "R": R,
        "Vobs": Vobs,
        "eV": eV,
        "Vgas": Vgas,
        "Vdisk": Vdisk,
        "Vbul": Vbul,
        "Vmax": vmax,
        "z": z_gal,
        "env": env_info,
    })

print(f"Total arquivos processados: {len(gal_records)}")


def evaluate_subset(subset_name, filter_fn, point_err_cut=None):
    res = {
        "RAR_padrao": [],
        "explicit_baryonic_exp_invalid": [],
        "v13_Poisson_sem_EFE": [],
        "v13_Poisson_EFE_ChaeMed": [],
        "v13_Poisson_EFE_ChaeMax": [],
        "RAR_EFE_ChaeMed": [],
        "RAR_EFE_ChaeMax": [],
        "v13_Poisson_Hubble_z": [],
        "v13_Poisson_DESI_z": [],
    }
    ngal = 0
    z_vals = []

    for g in gal_records:
        if not filter_fn(g):
            continue
        R, Vobs, eV = g["R"], g["Vobs"], g["eV"]
        mask = (R > 0) & (Vobs > 0)
        if point_err_cut is not None:
            mask = mask & (eV / np.maximum(Vobs, 1e-9) < point_err_cut)
        if mask.sum() == 0:
            continue
        gobs = (Vobs[mask] ** 2 / R[mask]) * KMS2_PER_KPC_TO_SI
        gbar = (
            np.sign(g["Vgas"][mask]) * g["Vgas"][mask] ** 2
            + 0.5 * np.sign(g["Vdisk"][mask]) * g["Vdisk"][mask] ** 2
            + 0.7 * g["Vbul"][mask] ** 2
        ) / R[mask] * KMS2_PER_KPC_TO_SI
        ok = gbar > 0
        if ok.sum() == 0:
            continue

        ngal += 1
        z_vals.append(g["z"])
        go = gobs[ok]
        gb = gbar[ok]
        log_go = np.log10(go)

        eN_med = g["env"]["eN_med"]
        eN_max = g["env"]["eN_max"]

        res["RAR_padrao"].append(log_go - np.log10(pred_rar_std(gb, A0_CANON)))
        res["explicit_baryonic_exp_invalid"].append(log_go - np.log10(pred_explicit_baryonic_exp_invalid(gb, A0_CANON)))
        res["v13_Poisson_sem_EFE"].append(log_go - np.log10(pred_poisson_mu_isolated(gb, A0_CANON)))
        res["v13_Poisson_EFE_ChaeMed"].append(log_go - np.log10(pred_poisson_mu_efe_sub(gb, A0_CANON, eN_med)))
        res["v13_Poisson_EFE_ChaeMax"].append(log_go - np.log10(pred_poisson_mu_efe_sub(gb, A0_CANON, eN_max)))
        res["RAR_EFE_ChaeMed"].append(log_go - np.log10(pred_rar_efe_chae(gb, A0_CANON, eN_med)))
        res["RAR_EFE_ChaeMax"].append(log_go - np.log10(pred_rar_efe_chae(gb, A0_CANON, eN_max)))
        res["v13_Poisson_Hubble_z"].append(log_go - np.log10(pred_poisson_mu_isolated(gb, a0_hubble(g["z"]))))
        res["v13_Poisson_DESI_z"].append(log_go - np.log10(pred_poisson_mu_isolated(gb, a0_desi_cpl(g["z"]))))

    print(f"\n======================================================================")
    print(f"SUBCONJUNTO: {subset_name} | Galaxias: {ngal} | Pontos: {sum(len(x) for x in res['RAR_padrao'])}")
    print(f"Redshift z das galaxias: min={min(z_vals):.5f}, med={np.median(z_vals):.5f}, max={max(z_vals):.5f}")
    print(f"----------------------------------------------------------------------")
    print(f"{'Modelo':<30} | {'RMS (dex)':>10} | {'Std (dex)':>10} | {'Bias Mediano':>12} | {'Bias Medio':>11}")
    print(f"----------------------------------------------------------------------")
    summary = {}
    for k, arr_list in res.items():
        r = np.concatenate(arr_list)
        rms = float(np.sqrt(np.mean(r ** 2)))
        std = float(np.std(r))
        med = float(np.median(r))
        mean = float(np.mean(r))
        summary[k] = {"rms": rms, "std": std, "median_bias": med, "mean_bias": mean, "n_pts": int(len(r)), "n_gal": ngal}
        print(f"{k:<30} | {rms:10.4f} | {std:10.4f} | {med:+12.4f} | {mean:+11.4f}")
    return summary


s_all = evaluate_subset("1. SPARC Completo (175 galaxias, sem cortes)", lambda g: True, point_err_cut=None)
s_can = evaluate_subset("2. SPARC Canonico (Q<=2, i>30, eV/V<0.10)", lambda g: g["meta"]["Q"] <= 2 and g["meta"]["Inc"] > 30, point_err_cut=0.10)
s_dw_all = evaluate_subset("3. Anas SPARC Completo (Vmax < 80 km/s, sem cortes)", lambda g: g["Vmax"] < 80.0, point_err_cut=None)
s_dw_can = evaluate_subset("4. Anas SPARC Canonicas (Vmax < 80 km/s, Q<=2, i>30, eV/V<0.10)", lambda g: g["Vmax"] < 80.0 and g["meta"]["Q"] <= 2 and g["meta"]["Inc"] > 30, point_err_cut=0.10)
s_dw_chae = evaluate_subset("5. Anas com EFE medido em Chae+2021 (Vmax < 80 km/s, Q<=2, i>30)", lambda g: g["Vmax"] < 80.0 and g["meta"]["Q"] <= 2 and g["meta"]["Inc"] > 30 and g["env"]["has_chae"], point_err_cut=0.10)

# Checar tambem a inclinacao da BTFR (Lei de Tully-Fisher Barionica) para provar por que g_bar/(1-e^(-g_bar/a0)) quebra!
print("\n=== TESTE DE SANIDADE ASSINTOTICA (Borda Galactica g_bar = 1e-12 m/s^2 = 0.0083 * a0) ===")
gb_test = np.array([1e-12])
print(f"  g_bar                  = {gb_test[0]:.3e} m/s^2")
print(f"  RAR Padrao             = {pred_rar_std(gb_test, A0_CANON)[0]:.3e} m/s^2 (= sqrt(g_bar*a0) = {np.sqrt(gb_test[0]*A0_CANON):.3e})")
print(f"  v1.3 Poisson (x=g_obs) = {pred_poisson_mu_isolated(gb_test, A0_CANON)[0]:.3e} m/s^2 (= sqrt(g_bar*a0) = {np.sqrt(gb_test[0]*A0_CANON):.3e})")
print(f"  Explicit baryonic exp  = {pred_explicit_baryonic_exp_invalid(gb_test, A0_CANON)[0]:.3e} m/s^2 (freezes at a0 = 1.20e-10 -> V grows as sqrt(R))")

out_json = os.path.join(PKG_ROOT, "data", "audit_v13_results.json")
with open(out_json, "w", encoding="utf-8") as fh:
    json.dump({"all_175": s_all, "canonical_144": s_can, "dwarfs_all": s_dw_all, "dwarfs_can_40": s_dw_can, "dwarfs_chae": s_dw_chae}, fh, indent=2)
print(f"\nSaved JSON to: {out_json}")
