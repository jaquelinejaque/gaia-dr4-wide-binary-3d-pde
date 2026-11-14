#!/usr/bin/env python3
"""
Standalone Sparse Newton-Raphson Axisymmetric PDE Solver for the v1.3 Extension
(AQUAL Lagrangian with Milgrom 1999 Exponential Poisson Interpolating Function):
  -div( (1 - exp(-|grad Phi|/a0)) * grad Phi ) = 4 * pi * G * rho
with canonical scale a0 = 1.20e-10 m s^-2 and External Field Effect (ub_z > 0).

Reproduces the v1.3 comparison columns in Table 2 (Gaia wide binaries) and
Table 3 (Classical Milky Way dwarf spheroidals) of Martins (2026).
"""

import json
import math
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

G_SI = 6.67430e-11
M_SUN = 1.98847e30
PC_M = 3.08567758149e16
KPC_M = 1.0e3 * PC_M
KAU_M = 1.0e3 * 1.495978707e11
A0_CANON = 1.20e-10


class PoissonAQUALPDESolver:
    def __init__(self, nr, nz, xmax, kappa, ub_z):
        self.nr, self.nz, self.xmax = nr, nz, xmax
        self.dx = xmax / nr
        self.dz = 2.0 * xmax / nz
        self.R = (np.arange(nr) + 0.5) * self.dx
        self.Z = -xmax + (np.arange(nz) + 0.5) * self.dz
        self.RR, self.ZZ = np.meshgrid(self.R, self.Z, indexing="ij")
        self.x = np.sqrt(self.RR ** 2 + self.ZZ ** 2)
        self.rho_nd = 3.0 * (1.0 + self.x ** 2) ** (-2.5)
        self.kappa = kappa
        self.ub_z = ub_z
        self.Rf = np.arange(nr + 1) * self.dx
        self.mu_e = 1.0 - math.exp(-ub_z)
        self.L_e = (ub_z * math.exp(-ub_z)) / self.mu_e
        lam = 1.0 + self.L_e
        self.bc_R = kappa / (self.mu_e * np.sqrt(lam * self.xmax ** 2 + self.Z ** 2))
        self.bc_Z0 = kappa / (self.mu_e * np.sqrt(lam * self.R ** 2 + self.xmax ** 2))
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

    @staticmethod
    def mu_and_dmu(A):
        ex = np.exp(-np.minimum(A, 500.0))
        return 1.0 - ex, ex

    def eval_res_jac(self, u, build_jac=True):
        nr, nz, dx, dz, R, Rf, ub = self.nr, self.nz, self.dx, self.dz, self.R, self.Rf, self.ub_z
        eps2 = 1e-24
        dR_c, dZ_c = self.grad_self(u)

        wR_fR = (u[1:, :] - u[:-1, :]) / dx
        wZ_fR = 0.5 * (dZ_c[:-1, :] + dZ_c[1:, :]) + ub
        A_fR = np.sqrt(wR_fR ** 2 + wZ_fR ** 2 + eps2)
        mu_fR, mp_fR = self.mu_and_dmu(A_fR)
        FR_int = mu_fR * wR_fR

        wR_bR = (self.bc_R - u[-1, :]) / (0.5 * dx)
        wZ_bR = dZ_c[-1, :] + ub
        A_bR = np.sqrt(wR_bR ** 2 + wZ_bR ** 2 + eps2)
        mu_bR, mp_bR = self.mu_and_dmu(A_bR)
        FR_bnd = mu_bR * wR_bR

        wZ_fZ = (u[:, 1:] - u[:, :-1]) / dz + ub
        wR_fZ = 0.5 * (dR_c[:, :-1] + dR_c[:, 1:])
        A_fZ = np.sqrt(wR_fZ ** 2 + wZ_fZ ** 2 + eps2)
        mu_fZ, mp_fZ = self.mu_and_dmu(A_fZ)
        FZ_int = mu_fZ * wZ_fZ

        wZ_bZ0 = (u[:, 0] - self.bc_Z0) / (0.5 * dz) + ub
        A_bZ0 = np.sqrt(dR_c[:, 0] ** 2 + wZ_bZ0 ** 2 + eps2)
        mu_bZ0, mp_bZ0 = self.mu_and_dmu(A_bZ0)
        FZ_b0 = mu_bZ0 * wZ_bZ0

        wZ_bZ1 = (self.bc_Z1 - u[:, -1]) / (0.5 * dz) + ub
        A_bZ1 = np.sqrt(dR_c[:, -1] ** 2 + wZ_bZ1 ** 2 + eps2)
        mu_bZ1, mp_bZ1 = self.mu_and_dmu(A_bZ1)
        FZ_b1 = mu_bZ1 * wZ_bZ1

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

        t_fR = mu_fR + mp_fR * (wR_fR ** 2) / A_fR
        t_bR = mu_bR + mp_bR * (wR_bR ** 2) / A_bR
        t_fZ = mu_fZ + mp_fZ * (wZ_fZ ** 2) / A_fZ
        t_bZ0 = mu_bZ0 + mp_bZ0 * (wZ_bZ0 ** 2) / A_bZ0
        t_bZ1 = mu_bZ1 + mp_bZ1 * (wZ_bZ1 ** 2) / A_bZ1

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
        np.add.at(diag, I[:, 0].ravel(), (t_bZ0 / (dz * dz) * 2.0).ravel())
        np.add.at(diag, I[:, -1].ravel(), (t_bZ1 / (dz * dz) * 2.0).ravel())

        rows.append(I.ravel()); cols.append(I.ravel()); vals.append(diag.ravel())
        J = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nr * nz, nr * nz))
        return res, J

    def solve(self):
        lam = 1.0 + self.L_e
        u = self.kappa / (self.mu_e * np.sqrt(lam * self.RR ** 2 + self.ZZ ** 2 + 1.5))
        res, _ = self.eval_res_jac(u, False)
        ires = float(np.max(np.abs(res)))
        for it in range(1, 30):
            res, J = self.eval_res_jac(u, True)
            nr_val = float(np.max(np.abs(res)))
            if nr_val / max(ires, 1e-30) < 1e-10 and nr_val < 1e-10:
                return u, ires, nr_val, it
            du = spla.spsolve(J.tocsc(), -res.ravel()).reshape(self.nr, self.nz)
            alpha = 1.0
            for _ in range(8):
                rt, _ = self.eval_res_jac(u + alpha * du, False)
                if float(np.max(np.abs(rt))) < nr_val * (1.0 - 0.05 * alpha) or alpha < 0.02:
                    break
                alpha *= 0.5
            u = u + alpha * du
        res, _ = self.eval_res_jac(u, False)
        return u, ires, float(np.max(np.abs(res))), 30


def solve_ub_from_gN(gN_SI):
    y = gN_SI / A0_CANON
    return brentq(lambda x: (1.0 - math.exp(-x)) * x - y, 1e-8, 100.0)


if __name__ == "__main__":
    M_MW = 6.5e10 * M_SUN
    r12_fac = 1.0 / math.sqrt(2.0 ** (2.0 / 3.0) - 1.0)
    dsph_out = []
    for nm, LV6, r12_pc, sig_obs, d_kpc, ml in [
        ("Fornax", 20.40, 944.0, 11.7, 149.0, 1.2),
        ("Leo I", 5.50, 388.0, 9.2, 258.0, 1.2),
        ("Leo II", 0.74, 233.0, 6.6, 236.0, 1.6),
        ("Sculptor", 2.30, 375.0, 9.2, 86.0, 1.8),
        ("Carina", 0.38, 334.0, 6.6, 107.0, 1.2),
        ("Sextans", 0.44, 925.0, 7.9, 89.0, 1.8),
        ("Draco", 0.27, 291.0, 9.1, 76.0, 1.8),
        ("Ursa Minor", 0.29, 240.0, 9.5, 78.0, 1.8),
    ]:
        rp_m = (r12_pc / r12_fac) * PC_M
        gN_MW = G_SI * M_MW / ((d_kpc * KPC_M) ** 2)
        ub = solve_ub_from_gN(gN_MW)
        M_star = ml * LV6 * 1.0e6 * M_SUN
        kappa = G_SI * M_star / (A0_CANON * rp_m ** 2)
        s = PoissonAQUALPDESolver(180, 360, 25.0, kappa, ub)
        u, ir, fr, it = s.solve()
        dR, dZ = s.grad_self(u)
        gtot_r = -(s.RR * dR + s.ZZ * dZ) / s.x
        gN_r = kappa * s.x / (1.0 + s.x ** 2) ** 1.5
        sh = (s.x > 1.22) & (s.x < 1.39)
        ws = s.RR[sh] / s.x[sh]
        geff = float(np.sum((gtot_r[sh] / gN_r[sh]) * ws) / np.sum(ws))
        sigN = math.sqrt(G_SI * (0.5 * M_star) / (3.0 * r12_pc * PC_M)) / 1.0e3
        sig_p = sigN * math.sqrt(geff)
        dsph_out.append({"galaxy": nm, "ML_pop": ml, "Geff_over_G": geff, "sigma_PDE_kms": sig_p, "sigma_obs_kms": sig_obs, "gamma_sigma": sig_p / sig_obs})
        print(f"{nm:11s}: Geff/G={geff:6.2f} | sig_PDE={sig_p:5.2f} vs obs={sig_obs:5.2f} -> gamma={sig_p/sig_obs:.3f}")

    wb_out = {}
    for env_lbl, ub_wb in [("Solar_gtot_1p9e-10", 1.9e-10 / A0_CANON), ("Solar_gN_bar_1p2e-10", solve_ub_from_gN(1.2e-10))]:
        gam_dict = {}
        for s_kau in [1.5, 15.0, 25.0]:
            s_m = s_kau * KAU_M
            rp_m = s_m / 12.5
            kappa = G_SI * (1.5 * M_SUN) / (A0_CANON * rp_m ** 2)
            sv = PoissonAQUALPDESolver(180, 360, 40.0, kappa, ub_wb)
            u, _, _, _ = sv.solve()
            dR, dZ = sv.grad_self(u)
            gtot_r = -(sv.RR * dR + sv.ZZ * dZ) / sv.x
            gN_r = kappa * sv.x / (1.0 + sv.x ** 2) ** 1.5
            sh = (sv.x > 12.0) & (sv.x < 13.0)
            ws = sv.RR[sh] / sv.x[sh]
            gam_dict[str(s_kau)] = float(np.sum((gtot_r[sh] / gN_r[sh]) * ws) / np.sum(ws))
        g_10_30 = 0.5 * (gam_dict["15.0"] + gam_dict["25.0"])
        wb_out[env_lbl] = {"bins": gam_dict, "gamma_rel": g_10_30 / gam_dict["1.5"]}
        print(f"{env_lbl}: gamma_rel = {wb_out[env_lbl]['gamma_rel']:.4f}")

    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_json = os.path.join(pkg_root, "data", "pde_v13_poisson_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump({"dsphs": dsph_out, "wide_binaries": wb_out}, fh, indent=2)
