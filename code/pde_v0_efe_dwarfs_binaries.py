#!/usr/bin/env python3
"""
Nonlinear 3D Axisymmetric PDE Solver for the v0 Cubic-Gradient Scalar Model with EFE (ub_z > 0)
Method: Conservative Cylindrical Finite-Volume Discretisation + Sparse Newton-Raphson (Analytical Jacobian + Armijo Line Search)
Author: Jaqueline Martins (2026)

Targets solved without linearisation approximations:
  1. Weak-EFE Analytical Benchmark: numerical comparison against u = kappa / (ub_z * sqrt(2*R^2 + Z^2))
  2. Classical Milky Way Dwarf Spheroidals (Fornax, Leo I, Leo II, Sculptor, Carina, Sextans, Draco, Ursa Minor)
     + 2x mesh refinement (200x400 -> 400x800) on Fornax, Leo I, and Leo II
  3. Gaia Wide Binaries (1-30 kAU, M_tot = 1.5 M_sun, frozen a_c = 1.9e-11 m s^-2):
     - Radial acceleration ratio gamma_a(s) across frozen separation bins
     - Relative boost gamma_rel = gamma_a(10-30 kAU) / gamma_a(1-2 kAU)
     - Sky-projected quadrupole anisotropy Delta_gamma = gamma_perp - gamma_par
     - Explicit two-body check (M1 = 0.85 M_sun, M2 = 0.65 M_sun)
"""

import json
import math
import os
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# Constantes fisicas SI e constante congelada a_c (12/11/2026)
G_SI = 6.67430e-11
M_SUN = 1.98847e30
PC_M = 3.08567758149e16
KPC_M = 1.0e3 * PC_M
AU_M = 1.495978707e11
KAU_M = 1.0e3 * AU_M
A_C_FROZEN = 1.9e-11  # m/s^2 (congelado no SPARC, proibido reajustar)


class SparseNewtonAxisymPDE:
    """
    Resolve R(u) = -div( |grad u + ub_z e_z| * (grad u + ub_z e_z) ) - kappa * rho_nd = 0
    onde u = u_self e grad u aponta radialmente para dentro (u > 0 na origem, u -> 0 no infinito).
    """

    def __init__(self, nr, nz, xmax, kappa, ub_z=0.0, rho_nd=None):
        self.nr, self.nz, self.xmax = nr, nz, xmax
        self.dx = xmax / nr
        self.dz = 2.0 * xmax / nz
        self.R = (np.arange(nr) + 0.5) * self.dx
        self.Z = -xmax + (np.arange(nz) + 0.5) * self.dz
        self.RR, self.ZZ = np.meshgrid(self.R, self.Z, indexing="ij")
        self.x = np.sqrt(self.RR ** 2 + self.ZZ ** 2)
        self.rho_nd = 3.0 * (1.0 + self.x ** 2) ** (-2.5) if rho_nd is None else rho_nd
        self.kappa = kappa
        self.ub_z = ub_z
        self.Rf = np.arange(nr + 1) * self.dx
        # Volume de cada celula anular: 2*pi*R_i*dx*dz (para verificar massa total / dois corpos)
        self.dV = 2.0 * math.pi * self.R[:, None] * self.dx * self.dz
        self.bc_R = np.zeros(nz)
        self.bc_Z0 = np.zeros(nr)
        self.bc_Z1 = np.zeros(nr)

    def set_asymptotic_bc(self, mass_norm=1.0):
        eff_k = self.kappa * mass_norm
        if abs(self.ub_z) < 1e-12:
            self.bc_R = -math.sqrt(eff_k) * np.log(np.sqrt(self.xmax ** 2 + self.Z ** 2) / self.xmax)
            self.bc_Z0 = -math.sqrt(eff_k) * np.log(np.sqrt(self.R ** 2 + self.xmax ** 2) / self.xmax)
            self.bc_Z1 = self.bc_Z0.copy()
        else:
            self.bc_R = eff_k / (self.ub_z * np.sqrt(2.0 * self.xmax ** 2 + self.Z ** 2))
            self.bc_Z0 = eff_k / (self.ub_z * np.sqrt(2.0 * self.R ** 2 + self.xmax ** 2))
            self.bc_Z1 = self.bc_Z0.copy()

    def grad_self(self, u):
        dx, dz = self.dx, self.dz
        dR = np.empty_like(u)
        dR[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2.0 * dx)
        dR[0, :] = (u[1, :] - u[0, :]) / (2.0 * dx)
        dR[-1, :] = (self.bc_R - u[-2, :]) / (1.5 * dx)

        dZ = np.empty_like(u)
        dZ[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2.0 * dz)
        dZ[:, 0] = (u[:, 1] - self.bc_Z0) / (1.5 * dz)
        dZ[:, -1] = (self.bc_Z1 - u[:, -2]) / (1.5 * dz)
        return dR, dZ

    def eval_residual_and_jacobian(self, u, build_jac=True):
        nr, nz, dx, dz, R, Rf = self.nr, self.nz, self.dx, self.dz, self.R, self.Rf
        eps2 = 1e-24
        ub = self.ub_z
        dR_c, dZ_c = self.grad_self(u)

        # 1. Faces R internas
        wR_fR = (u[1:, :] - u[:-1, :]) / dx
        wZ_fR = 0.5 * (dZ_c[:-1, :] + dZ_c[1:, :]) + ub
        A_fR = np.sqrt(wR_fR ** 2 + wZ_fR ** 2 + eps2)
        FR_int = A_fR * wR_fR

        # 2. Face R externa (R = xmax)
        wR_bR = (self.bc_R - u[-1, :]) / (0.5 * dx)
        wZ_bR = dZ_c[-1, :] + ub
        A_bR = np.sqrt(wR_bR ** 2 + wZ_bR ** 2 + eps2)
        FR_bnd = A_bR * wR_bR

        # 3. Faces Z internas
        wZ_fZ = (u[:, 1:] - u[:, :-1]) / dz + ub
        wR_fZ = 0.5 * (dR_c[:, :-1] + dR_c[:, 1:])
        A_fZ = np.sqrt(wR_fZ ** 2 + wZ_fZ ** 2 + eps2)
        FZ_int = A_fZ * wZ_fZ

        # 4. Faces Z externas (Z = -xmax, Z = +xmax)
        wZ_bZ0 = (u[:, 0] - self.bc_Z0) / (0.5 * dz) + ub
        wR_bZ0 = dR_c[:, 0]
        A_bZ0 = np.sqrt(wR_bZ0 ** 2 + wZ_bZ0 ** 2 + eps2)
        FZ_b0 = A_bZ0 * wZ_bZ0

        wZ_bZ1 = (self.bc_Z1 - u[:, -1]) / (0.5 * dz) + ub
        wR_bZ1 = dR_c[:, -1]
        A_bZ1 = np.sqrt(wR_bZ1 ** 2 + wZ_bZ1 ** 2 + eps2)
        FZ_b1 = A_bZ1 * wZ_bZ1

        res = -self.kappa * self.rho_nd.copy()
        res[:-1, :] -= (Rf[1:nr, None] * FR_int) / (R[:-1, None] * dx)
        res[1:, :] += (Rf[1:nr, None] * FR_int) / (R[1:, None] * dx)
        res[-1, :] -= (Rf[nr] * FR_bnd) / (R[-1] * dx)

        res[:, :-1] -= FZ_int / dz
        res[:, 1:] += FZ_int / dz
        res[:, 0] += FZ_b0 / dz
        res[:, -1] -= FZ_b1 / dz

        if not build_jac:
            return res, None

        t_fR = A_fR + (wR_fR ** 2) / A_fR
        t_bR = A_bR + (wR_bR ** 2) / A_bR
        t_fZ = A_fZ + (wZ_fZ ** 2) / A_fZ
        t_bZ0 = A_bZ0 + (wZ_bZ0 ** 2) / A_bZ0
        t_bZ1 = A_bZ1 + (wZ_bZ1 ** 2) / A_bZ1

        I = np.arange(nr * nz).reshape(nr, nz)
        rows, cols, vals = [], [], []
        diag = np.zeros(nr * nz)

        c_out = t_fR * Rf[1:nr, None] / (R[:-1, None] * dx * dx)
        rows.append(I[:-1, :].ravel()); cols.append(I[1:, :].ravel()); vals.append((-c_out).ravel())
        np.add.at(diag, I[:-1, :].ravel(), c_out.ravel())

        c_in = t_fR * Rf[1:nr, None] / (R[1:, None] * dx * dx)
        rows.append(I[1:, :].ravel()); cols.append(I[:-1, :].ravel()); vals.append((-c_in).ravel())
        np.add.at(diag, I[1:, :].ravel(), c_in.ravel())

        cRb = t_bR * Rf[nr] / (R[-1] * dx * dx) * 2.0
        np.add.at(diag, I[-1, :].ravel(), cRb.ravel())

        cz = t_fZ / (dz * dz)
        rows.append(I[:, :-1].ravel()); cols.append(I[:, 1:].ravel()); vals.append((-cz).ravel())
        np.add.at(diag, I[:, :-1].ravel(), cz.ravel())
        rows.append(I[:, 1:].ravel()); cols.append(I[:, :-1].ravel()); vals.append((-cz).ravel())
        np.add.at(diag, I[:, 1:].ravel(), cz.ravel())

        cZ0 = t_bZ0 / (dz * dz) * 2.0
        np.add.at(diag, I[:, 0].ravel(), cZ0.ravel())
        cZ1 = t_bZ1 / (dz * dz) * 2.0
        np.add.at(diag, I[:, -1].ravel(), cZ1.ravel())

        # Derivadas cruzadas transversais
        cr_fR = (wR_fR[:, 1:-1] * wZ_fR[:, 1:-1] / A_fR[:, 1:-1]) / (4.0 * dz)
        s_out = Rf[1:nr, None] / (R[:-1, None] * dx) * cr_fR
        s_in = Rf[1:nr, None] / (R[1:, None] * dx) * cr_fR
        for di, sgn_cell, coeff in [(0, -1.0, s_out), (1, +1.0, s_in)]:
            row_idx = I[di : nr - 1 + di, 1:-1].ravel()
            for nb_i in [0, 1]:
                rows.append(row_idx); cols.append(I[nb_i : nr - 1 + nb_i, 2:].ravel()); vals.append((sgn_cell * coeff).ravel())
                rows.append(row_idx); cols.append(I[nb_i : nr - 1 + nb_i, :-2].ravel()); vals.append((-sgn_cell * coeff).ravel())

        cr_fZ = (wR_fZ[1:-1, :] * wZ_fZ[1:-1, :] / A_fZ[1:-1, :]) / (4.0 * dx * dz)
        for dj, sgn_cell in [(0, -1.0), (1, +1.0)]:
            row_idx = I[1:-1, dj : nz - 1 + dj].ravel()
            for nb_j in [0, 1]:
                rows.append(row_idx); cols.append(I[2:, nb_j : nz - 1 + nb_j].ravel()); vals.append((sgn_cell * cr_fZ).ravel())
                rows.append(row_idx); cols.append(I[:-2, nb_j : nz - 1 + nb_j].ravel()); vals.append((-sgn_cell * cr_fZ).ravel())

        rows.append(I.ravel()); cols.append(I.ravel()); vals.append(diag.ravel())
        J = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nr * nz, nr * nz))
        return res, J

    def solve_newton(self, max_iters=28, tol=1e-11, mass_norm=1.0):
        self.set_asymptotic_bc(mass_norm=mass_norm)
        eff_k = self.kappa * mass_norm
        if abs(self.ub_z) < 1e-12:
            u = -math.sqrt(eff_k) * np.log(np.sqrt(1.0 + self.x ** 2) / self.xmax)
        else:
            u = eff_k / (self.ub_z * np.sqrt(2.0 * self.RR ** 2 + self.ZZ ** 2 + 2.0) + math.sqrt(eff_k))
        res, _ = self.eval_residual_and_jacobian(u, build_jac=False)
        ires = float(np.max(np.abs(res)))
        for it in range(1, max_iters + 1):
            res, J = self.eval_residual_and_jacobian(u, build_jac=True)
            norm_r = float(np.max(np.abs(res)))
            if norm_r / max(ires, 1e-30) < tol and norm_r < 1e-10:
                return u, ires, norm_r, it
            delta = spla.spsolve(J.tocsc(), -res.ravel()).reshape(self.nr, self.nz)
            alpha = 1.0
            for _ in range(10):
                r_try, _ = self.eval_residual_and_jacobian(u + alpha * delta, build_jac=False)
                if float(np.max(np.abs(r_try))) < norm_r * (1.0 - 0.05 * alpha) or alpha < 0.01:
                    break
                alpha *= 0.5
            u = u + alpha * delta
        res, _ = self.eval_residual_and_jacobian(u, build_jac=False)
        return u, ires, float(np.max(np.abs(res))), max_iters


def run_weak_efe_benchmark():
    print("====================================================================================")
    print("SECAO 1 — BENCHMARK ANALITICO EFE FRACO (Funcao de Green Exata Anisotropica)")
    print("====================================================================================")
    kappa, ub, xmax = 0.01, 10.0, 35.0
    solver = SparseNewtonAxisymPDE(220, 440, xmax, kappa, ub)
    u, ires, fres, it = solver.solve_newton()
    dR, dZ = solver.grad_self(u)
    # Aceleracao radial interna para dentro: g_r = -(R*dR + Z*dZ)/x
    g_psi_r = -(solver.RR * dR + solver.ZZ * dZ) / solver.x
    g_N_r = kappa * solver.x / (1.0 + solver.x ** 2) ** 1.5

    # Medir na casca esferica x in [11.5, 12.5] (longe do nucleo Plummer r_p=1 e longe da borda xmax=35)
    shell = (solver.x > 11.5) & (solver.x < 12.5)
    # Media esferica ponderada por sin(theta) proporcional a R/x (elemento de area 2*pi*R*ds)
    w_sph = solver.RR[shell] / solver.x[shell]
    ratio_sph = float(np.sum(w_psi := (g_psi_r[shell] / g_N_r[shell]) * w_sph) / np.sum(w_sph))
    exact_sph = (math.pi / 4.0) / ub  # pi / (4 * ub) = 0.0785398

    iz_pos = np.argmin(np.abs(solver.Z - 12.0))
    ir_pos = np.argmin(np.abs(solver.R - 12.0))
    iz_mid = np.argsort(np.abs(solver.Z))[:2]
    g_par = -0.5 * (dZ[0, iz_pos] - dZ[0, solver.nz - 1 - iz_pos])
    g_perp = -0.5 * (dR[ir_pos, iz_mid[0]] + dR[ir_pos, iz_mid[1]])
    ratio_qp = float(g_perp / g_par)

    print(f"Newton iters={it} | Residuo: {ires:.3e} -> {fres:.3e} (queda {math.log10(ires/fres):.1f} ordens)")
    print(f"Media esferica <g_Psi / g_N> em x=12 : num={ratio_sph:.6f} vs analitico pi/(4*ub)={exact_sph:.6f} (erro={abs(ratio_sph-exact_sph)/exact_sph*100:.2f}%)")
    print(f"Anisotropia g_perp / g_par em x=12   : num={ratio_qp:.5f} vs limite pontual 1/sqrt(2)=0.70711 (g_perp < g_par!)")
    return {
        "iters": it,
        "drop_orders": math.log10(ires / fres),
        "spherical_boost_num": ratio_sph,
        "spherical_boost_exact": exact_sph,
        "quadrupole_perp_over_par": ratio_qp,
    }


def run_8_milky_way_dsphs():
    print("\n====================================================================================")
    print("SECAO 2 — AS 8 ANAS ESFEROIDAIS DA VIA LACTEA NA PDE NAO-LINEAR EXATA (SEM LINEARIZAR)")
    print("Parametros: Wolf+2010 / McGaugh & Wolf 2010 | M_bar,MW = 6.5e10 M_sun | a_c = 1.9e-11")
    print("====================================================================================")
    M_MW = 6.5e10 * M_SUN
    # Nome, L_V (10^6 L_sun), r_1/2_3D (pc), sigma_obs (km/s), e_sigma (km/s), D_GC (kpc), M/L_pop, tidal_flag
    dsphs = [
        ("Fornax",     20.40, 944.0, 11.7, 0.9, 149.0, 1.2, False),
        ("Leo I",       5.50, 388.0,  9.2, 1.4, 258.0, 1.2, False),
        ("Leo II",      0.74, 233.0,  6.6, 0.7, 236.0, 1.6, False),
        ("Sculptor",    2.30, 375.0,  9.2, 1.1,  86.0, 1.8, False),
        ("Carina",      0.38, 334.0,  6.6, 1.2, 107.0, 1.2, False),
        ("Sextans",     0.44, 925.0,  7.9, 1.3,  89.0, 1.8, True),
        ("Draco",       0.27, 291.0,  9.1, 1.2,  76.0, 1.8, True),
        ("Ursa Minor",  0.29, 240.0,  9.5, 1.2,  78.0, 1.8, True),
    ]
    plummer_r12_over_rp = 1.0 / math.sqrt(2.0 ** (2.0 / 3.0) - 1.0)  # 1.304766

    out_rows = []
    header = (
        f"{'Galaxia':11s} | {'M/L':3s} | {'kappa':6s} | {'ub_z':5s} | {'eta_1/2':7s} | "
        f"{'Geff_Lin':8s} | {'Geff_PDE':8s} | {'sig_obs':7s} | {'sig_PDE':7s} | {'gamma_Lin':9s} | {'gamma_PDE [ML 1.2-2.0]':22s}"
    )
    print(header)
    print("-" * len(header))

    for nm, LV6, r12_pc, sig_obs, esig, d_kpc, ml_pop, is_tidal in dsphs:
        rp_m = (r12_pc / plummer_r12_over_rp) * PC_M
        gN_MW = G_SI * M_MW / ((d_kpc * KPC_M) ** 2)
        ub_z = math.sqrt(gN_MW / A_C_FROZEN)

        # Rodar para M/L = 1.2, ml_pop e 2.0 na PDE nao-linear (pois kappa depende de M_star!)
        ml_cache = {}
        for ml in sorted(set([1.2, ml_pop, 2.0])):
            M_star = ml * LV6 * 1.0e6 * M_SUN
            kappa = G_SI * M_star / (A_C_FROZEN * rp_m ** 2)
            solver = SparseNewtonAxisymPDE(200, 400, 25.0, kappa, ub_z)
            u, ires, fres, it = solver.solve_newton()
            dR, dZ = solver.grad_self(u)
            g_psi_r = -(solver.RR * dR + solver.ZZ * dZ) / solver.x
            g_N_r = kappa * solver.x / (1.0 + solver.x ** 2) ** 1.5

            # 1) Estimador na casca r_1/2 (x = 1.3048, identico ao ponto de Wolf+2010)
            sh = (solver.x > 1.22) & (solver.x < 1.39)
            w_s = solver.RR[sh] / solver.x[sh]
            geff_half = 1.0 + float(np.sum((g_psi_r[sh] / g_N_r[sh]) * w_s) / np.sum(w_s))
            eta_half = float(np.sqrt(np.sum(((dR[sh] ** 2 + dZ[sh] ** 2)) * w_s) / np.sum(w_s))) / ub_z

            # 2) Teorema do Virial 3D exato: 3*M*sigma_los^2 = int rho(x) * (r . g_tot) dV
            r_dot_g_psi = -(solver.RR * dR + solver.ZZ * dZ)
            r_dot_g_N = kappa * (solver.x ** 2) / (1.0 + solver.x ** 2) ** 1.5
            W_N = float(np.sum(solver.rho_nd * r_dot_g_N * solver.dV))
            W_psi = float(np.sum(solver.rho_nd * r_dot_g_psi * solver.dV))
            geff_vir = 1.0 + W_psi / W_N

            # Dispersao Newtoniana Wolf em r_1/2: sigma_N^2 = G*(0.5*M_star)/(3*r_1/2)
            sig_N_kms = math.sqrt(G_SI * (0.5 * M_star) / (3.0 * r12_pc * PC_M)) / 1.0e3
            sig_pde_half = sig_N_kms * math.sqrt(geff_half)
            sig_pde_vir = sig_N_kms * math.sqrt(geff_vir)

            geff_lin_r12 = 1.0 + 0.75 * math.sqrt(A_C_FROZEN / gN_MW)
            sig_lin = sig_N_kms * math.sqrt(geff_lin_r12)

            ml_cache[ml] = {
                "kappa": kappa,
                "iters": it,
                "drop_orders": math.log10(ires / max(fres, 1e-30)),
                "eta_half": eta_half,
                "geff_lin": geff_lin_r12,
                "geff_pde_half": geff_half,
                "geff_pde_virial": geff_vir,
                "sig_N_kms": sig_N_kms,
                "sig_lin_kms": sig_lin,
                "sig_pde_kms": sig_pde_half,
                "sig_pde_vir_kms": sig_pde_vir,
                "gamma_lin": sig_lin / sig_obs,
                "gamma_pde": sig_pde_half / sig_obs,
                "gamma_pde_vir": sig_pde_vir / sig_obs,
            }

        c_pop = ml_cache[ml_pop]
        c_lo = ml_cache[1.2]
        c_hi = ml_cache[2.0]

        # Refinamento 2x (400x800) nas 3 anas distantes (Fornax, Leo I, Leo II)
        mesh_diff_pct = None
        if nm in ("Fornax", "Leo I", "Leo II"):
            s_fine = SparseNewtonAxisymPDE(400, 800, 25.0, c_pop["kappa"], ub_z)
            uf, _, _, _ = s_fine.solve_newton()
            dRf, dZf = s_fine.grad_self(uf)
            gp_f = -(s_fine.RR * dRf + s_fine.ZZ * dZf) / s_fine.x
            gn_f = c_pop["kappa"] * s_fine.x / (1.0 + s_fine.x ** 2) ** 1.5
            shf = (s_fine.x > 1.22) & (s_fine.x < 1.39)
            wf = s_fine.RR[shf] / s_fine.x[shf]
            geff_f = 1.0 + float(np.sum((gp_f[shf] / gn_f[shf]) * wf) / np.sum(wf))
            mesh_diff_pct = abs(math.sqrt(geff_f) - math.sqrt(c_pop["geff_pde_half"])) / math.sqrt(c_pop["geff_pde_half"]) * 100.0

        print(
            f"{nm:11s} | {ml_pop:3.1f} | {c_pop['kappa']:6.3f} | {ub_z:5.3f} | {c_pop['eta_half']:7.2f} | "
            f"{c_pop['geff_lin']:8.2f} | {c_pop['geff_pde_half']:8.2f} | {sig_obs:7.2f} | {c_pop['sig_pde_kms']:7.2f} | "
            f"{c_pop['gamma_lin']:9.2f} | {c_pop['gamma_pde']:5.2f} [{c_lo['gamma_pde']:.2f}-{c_hi['gamma_pde']:.2f}]"
            + (f" (malha 2x: {mesh_diff_pct:.3f}%)" if mesh_diff_pct is not None else "")
        )

        out_rows.append({
            "galaxy": nm,
            "tidal_system": is_tidal,
            "D_GC_kpc": d_kpc,
            "sigma_obs_kms": sig_obs,
            "e_sigma_kms": esig,
            "ub_z": ub_z,
            "ml_pop": ml_pop,
            "mesh_2x_change_pct": mesh_diff_pct,
            "pop_result": c_pop,
            "gamma_pde_range_ML_1p2_to_2p0": [c_lo["gamma_pde"], c_hi["gamma_pde"]],
            "gamma_lin_range_ML_1p2_to_2p0": [c_lo["gamma_lin"], c_hi["gamma_lin"]],
        })
    return out_rows


def run_gaia_wide_binaries():
    print("\n====================================================================================")
    print("SECAO 3 — BINARIOS LARGOS DO GAIA (1-30 kAU) NA PDE NAO-LINEAR EXATA + QUADRUPOLO")
    print("Parametros congelados (PREREGISTRO_DR4.md): M_tot = 1.5 M_sun | a_c = 1.9e-11 m/s^2")
    print("====================================================================================")
    M_tot = 1.5 * M_SUN
    results_wb = {}

    for env_tag, g_ext_SI in [
        ("Convencao_Preregistro_aext_1p9e-10", 1.9e-10),
        ("Convencao_Barionica_MW_gN_1p2e-10", 1.2e-10),
    ]:
        ub_z = math.sqrt(g_ext_SI / A_C_FROZEN)
        print(f"\n--- {env_tag} (g_ext = {g_ext_SI:.2e} m/s^2, ub_z = {ub_z:.4f}) ---")
        seps_kau = [1.5, 3.5, 7.5, 15.0, 20.0, 25.0]
        bin_rows = []

        # Para fonte compacta r_p = 0.08 * s (estrela pontual vista na escala orbital x_orb = s / r_p = 12.5, xmax = 40.0)
        x_orb = 12.5
        xmax = 40.0
        for s_kau in seps_kau:
            s_m = s_kau * KAU_M
            rp_m = s_m / x_orb
            kappa = G_SI * M_tot / (A_C_FROZEN * rp_m ** 2)

            solver = SparseNewtonAxisymPDE(220, 440, xmax, kappa, ub_z)
            u, ires, fres, it = solver.solve_newton()
            dR, dZ = solver.grad_self(u)
            g_psi_r = -(solver.RR * dR + solver.ZZ * dZ) / solver.x
            g_N_r = kappa * solver.x / (1.0 + solver.x ** 2) ** 1.5

            sh = (solver.x > x_orb - 0.5) & (solver.x < x_orb + 0.5)
            w_s = solver.RR[sh] / solver.x[sh]
            gamma_sph = 1.0 + float(np.sum((g_psi_r[sh] / g_N_r[sh]) * w_s) / np.sum(w_s))

            # Extrair gamma_par (theta ~ 0..25 deg, alinhado ao Centro Galactico) vs gamma_perp (theta ~ 65..90 deg, perpendicular)
            cos_th = np.abs(solver.ZZ[sh]) / solver.x[sh]
            sin_th = solver.RR[sh] / solver.x[sh]
            mask_par = cos_th > math.cos(math.radians(25.0))
            mask_perp = sin_th > math.sin(math.radians(65.0))

            gamma_par_3d = 1.0 + float(np.mean(g_psi_r[sh][mask_par] / g_N_r[sh][mask_par]))
            gamma_perp_3d = 1.0 + float(np.mean(g_psi_r[sh][mask_perp] / g_N_r[sh][mask_perp]))

            # Projecao real no ceu (GC esta a ~90 deg da linha de visada para binarios solares d < 200 pc!):
            # Eixo Z = direcao do Centro Galactico (no plano do ceu X_sky = Z).
            # Eixo Y_sky = perpendicular no plano do ceu; Eixo Z_los = linha de visada.
            # Para distribuicao isotropica 3D em (phi): R = r*sqrt(sin^2(th_sky)*cos^2(psi) + cos^2(...) ) ->
            # Fazendo media Monte Carlo exata sobre orientacoes 3D com angulo projetado theta_sky no ceu:
            rng = np.random.default_rng(42)
            # Amostrar pontos uniformes na esfera 3D: z_gc = cos(alpha), y_perp = sin(alpha)*cos(beta), x_los = sin(alpha)*sin(beta)
            n_mc = 200000
            z_gc = rng.uniform(-1.0, 1.0, n_mc)
            beta = rng.uniform(0.0, 2.0 * math.pi, n_mc)
            s_a = np.sqrt(1.0 - z_gc ** 2)
            y_perp = s_a * np.cos(beta)
            # Angulo projetado no ceu em relacao a direcao do GC (z_gc vs y_perp):
            th_sky = np.arctan2(np.abs(y_perp), np.abs(z_gc))  # [0, pi/2]
            # Interpolar boost radial 3D em funcao de |cos(alpha)| = |z_gc|
            order = np.argsort(cos_th)
            boost_3d_mc = np.interp(np.abs(z_gc), cos_th[order], (g_psi_r[sh] / g_N_r[sh])[order])
            gam_mc = 1.0 + boost_3d_mc
            gam_sky_par = float(np.mean(gam_mc[th_sky < math.pi / 4.0]))
            gam_sky_perp = float(np.mean(gam_mc[th_sky >= math.pi / 4.0]))
            quad_sky = gam_sky_perp - gam_sky_par

            # Formula escalar 1D de subtracao de ontem para comparacao direta:
            gN_int = G_SI * M_tot / (s_m ** 2)
            gam_1d_sub = 1.0 + (math.sqrt(A_C_FROZEN * (gN_int + g_ext_SI)) - math.sqrt(A_C_FROZEN * g_ext_SI)) / gN_int

            bin_rows.append({
                "s_kAU": s_kau,
                "gN_int_over_gext": gN_int / g_ext_SI,
                "gamma_1D_sub": gam_1d_sub,
                "gamma_PDE_sph": gamma_sph,
                "gamma_3D_par": gamma_par_3d,
                "gamma_3D_perp": gamma_perp_3d,
                "gamma_sky_par": gam_sky_par,
                "gamma_sky_perp": gam_sky_perp,
                "quadrupole_sky_perp_minus_par": quad_sky,
                "iters": it,
                "drop_orders": math.log10(ires / max(fres, 1e-30)),
            })

        # Teste Exato de Dois Corpos na PDE Nao-Linear em s = 1.5 kAU e s = 20.0 kAU (M1=0.85, M2=0.65 M_sun alinhados em Z)
        twobody_check = {}
        for s_tb in [1.5, 20.0]:
            s_m = s_tb * KAU_M
            # Escala da malha: separacao s = 8 unidades (Z1 = -3.467, Z2 = +4.533, baricentro em 0), nucleo compacto a_star = 0.35
            L_unit = s_m / 8.0
            a_star = 0.35
            q1, q2 = 0.85 / 1.5, 0.65 / 1.5
            Z1, Z2 = -q2 * 8.0, +q1 * 8.0
            kappa_tb = G_SI * M_tot / (A_C_FROZEN * L_unit ** 2)

            tb_solver = SparseNewtonAxisymPDE(200, 480, 28.0, kappa_tb, ub_z)
            r1_sq = tb_solver.RR ** 2 + (tb_solver.ZZ - Z1) ** 2
            r2_sq = tb_solver.RR ** 2 + (tb_solver.ZZ - Z2) ** 2
            rho1 = q1 * (3.0 / (a_star ** 3)) * (1.0 + r1_sq / (a_star ** 2)) ** (-2.5)
            rho2 = q2 * (3.0 / (a_star ** 2) / a_star) * (1.0 + r2_sq / (a_star ** 2)) ** (-2.5)
            # Normalizar cada massa exatamente na malha discreta para eliminar erro de discretizacao do nucleo
            rho1 *= (4.0 * math.pi * q1) / float(np.sum(rho1 * tb_solver.dV))
            rho2 *= (4.0 * math.pi * q2) / float(np.sum(rho2 * tb_solver.dV))
            tb_solver.rho_nd = rho1 + rho2

            u_tb, ir_tb, fr_tb, it_tb = tb_solver.solve_newton(mass_norm=1.0)
            _, dZ_tb = tb_solver.grad_self(u_tb)
            # Na estrela 2 (em Z2 > 0), a atracao mutua aponta para -Z, logo g_self_Z < 0 -> dZ_tb < 0 (pois grad u = g_inward).
            # Pela 3a lei de Newton (conservacao do momento da Lagrangiana!), F1_z + F2_z = 0!
            # Aceleracao relativa extra ao longo de Z: a_rel_psi = |<g_psi_Z>_2 - <g_psi_Z>_1|
            gZ_1 = float(np.sum(rho1 * dZ_tb * tb_solver.dV)) / (4.0 * math.pi * q1)
            gZ_2 = float(np.sum(rho2 * dZ_tb * tb_solver.dV)) / (4.0 * math.pi * q2)
            a_rel_psi_nd = abs(gZ_2 - gZ_1)  # Nota: auto-forca Newtoniana isolada integra zero por simetria
            # Para remover qualquer assimetria residual de malha em cada estrela isolada sob ub_z, subtrair estrela unica na mesma posicao!
            tb_s1 = SparseNewtonAxisymPDE(200, 480, 28.0, kappa_tb, ub_z, rho_nd=rho1)
            u_s1, _, _, _ = tb_s1.solve_newton(mass_norm=q1)
            _, dZ_s1 = tb_s1.grad_self(u_s1)
            gZ_1_self = float(np.sum(rho1 * dZ_s1 * tb_s1.dV)) / (4.0 * math.pi * q1)

            tb_s2 = SparseNewtonAxisymPDE(200, 480, 28.0, kappa_tb, ub_z, rho_nd=rho2)
            u_s2, _, _, _ = tb_s2.solve_newton(mass_norm=q2)
            _, dZ_s2 = tb_s2.grad_self(u_s2)
            gZ_2_self = float(np.sum(rho2 * dZ_s2 * tb_s2.dV)) / (4.0 * math.pi * q2)

            a_mut_psi_nd = abs((gZ_2 - gZ_2_self) - (gZ_1 - gZ_1_self))
            a_mut_N_nd = kappa_tb / (8.0 ** 2)
            twobody_check[f"{s_tb:.1f}_kAU"] = {
                "gamma_2body_aligned": 1.0 + a_mut_psi_nd / a_mut_N_nd,
                "iters": it_tb,
            }

        gam_ref = bin_rows[0]["gamma_PDE_sph"]  # 1.5 kAU (bin 1-2 kAU)
        gam_10_30 = 0.5 * (bin_rows[3]["gamma_PDE_sph"] + bin_rows[5]["gamma_PDE_sph"])
        gam_rel_pde = gam_10_30 / gam_ref
        gam_rel_1d = (0.5 * (bin_rows[3]["gamma_1D_sub"] + bin_rows[5]["gamma_1D_sub"])) / bin_rows[0]["gamma_1D_sub"]
        quad_10_30 = 0.5 * (bin_rows[3]["quadrupole_sky_perp_minus_par"] + bin_rows[5]["quadrupole_sky_perp_minus_par"])

        for r in bin_rows:
            print(
                f"  s={r['s_kAU']:4.1f} kAU | gN/gext={r['gN_int_over_gext']:6.3f} | "
                f"1D_sub={r['gamma_1D_sub']:.4f} | PDE_sph={r['gamma_PDE_sph']:.4f} (rel 1-2kAU = {r['gamma_PDE_sph']/gam_ref:.4f}) | "
                f"sky_perp={r['gamma_sky_perp']:.4f} vs sky_par={r['gamma_sky_par']:.4f} -> d_quad={r['quadrupole_sky_perp_minus_par']:+.4f}"
            )
        print(f"  >>> SINTESE EMENDA 3: gamma_rel(10-30 / 1-2) = {gam_rel_pde:.4f} (vs 1D subtracao = {gam_rel_1d:.4f})")
        print(f"  >>> SINTESE EMENDA 4: quadrupolo projetado (gamma_perp - gamma_par) em 10-30 kAU = {quad_10_30:+.4f}")
        print(f"  >>> CHECK 2 CORPOS ALINHADOS (M1=0.85, M2=0.65): 1.5kAU={twobody_check['1.5_kAU']['gamma_2body_aligned']:.4f} | 20kAU={twobody_check['20.0_kAU']['gamma_2body_aligned']:.4f} (ratio={twobody_check['20.0_kAU']['gamma_2body_aligned']/twobody_check['1.5_kAU']['gamma_2body_aligned']:.4f})")

        results_wb[env_tag] = {
            "g_ext_SI": g_ext_SI,
            "ub_z": ub_z,
            "bins": bin_rows,
            "gamma_rel_10_30_over_1_2_PDE": gam_rel_pde,
            "gamma_rel_10_30_over_1_2_1D": gam_rel_1d,
            "quadrupole_sky_10_30_perp_minus_par": quad_10_30,
            "twobody_aligned_check": twobody_check,
        }
    return results_wb


if __name__ == "__main__":
    t_all = time.time()
    res_bench = run_weak_efe_benchmark()
    res_dsphs = run_8_milky_way_dsphs()
    res_wb = run_gaia_wide_binaries()

    full_payload = {
        "elapsed_seconds": time.time() - t_all,
        "frozen_ac_SI": A_C_FROZEN,
        "weak_efe_benchmark": res_bench,
        "milky_way_8_dsphs": res_dsphs,
        "gaia_wide_binaries": res_wb,
    }
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_json = os.path.join(pkg_root, "data", "pde_v0_efe_dwarfs_binaries_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(full_payload, fh, indent=2)
    print(f"\nCOMPLETED in {time.time() - t_all:.1f}s | Saved to: {out_json}")
