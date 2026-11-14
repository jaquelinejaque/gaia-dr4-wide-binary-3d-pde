#!/usr/bin/env python3
"""
Isolated Validation & Convergence Suite for the Nonlinear Cylindrical PDE Solver (v0)
Equation: -div(|grad u| grad u) = kappa * rho_nd(x), with rho_nd(x) = 3 * (1 + x^2)^(-5/2) (normalised Plummer sphere, M_tot = 4*pi).
Author: Jaqueline Martins (2026)

Verification criteria (Table 1):
  1. Newtonian Machine Check (-div(grad u) = kappa * rho_nd): ratio vs exact analytical kappa / 2^(3/2)
  2. Diagnostic reproduction of unnormalised Plummer source (missing factor 3 -> ratio 1/sqrt(3) = 0.5774)
  3. Isolated Nonlinear Gauss Benchmark (200x400 baseline and 400x800 2x refined mesh):
     - Criterion 1: ratio = a_num / 2.659148 within 1.00 +- 0.02
     - Criterion 2: nonlinear residual reduction >= 6 orders of magnitude
     - Criterion 3: 2x mesh-refinement drift < 2.0%
"""

import math
import os
import time
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class ValidatedCylindricalPDESolver:
    def __init__(self, nr, nz, xmax, kappa, ub_z=0.0, bc_vals=None, newtonian=False, plume_factor=3.0):
        self.nr, self.nz, self.xmax = nr, nz, xmax
        self.dx = xmax / nr
        self.dz = 2.0 * xmax / nz
        self.R = (np.arange(nr) + 0.5) * self.dx
        self.Z = -xmax + (np.arange(nz) + 0.5) * self.dz
        self.RR, self.ZZ = np.meshgrid(self.R, self.Z, indexing="ij")
        self.x = np.sqrt(self.RR ** 2 + self.ZZ ** 2)
        self.rho_nd = plume_factor * (1.0 + self.x ** 2) ** (-2.5)
        self.kappa = kappa
        self.ub_z = ub_z
        self.newtonian = newtonian
        self.Rf = np.arange(nr + 1) * self.dx
        if bc_vals is None:
            self.bc_vals = np.zeros_like(self.x)
        else:
            self.bc_vals = bc_vals

    def _face_amplitudes(self, u):
        """
        Calcula |grad u + ub_z e_z| diretamente nas faces R e Z usando stencil consistente de 2a ordem,
        garantindo que o operador M(u)*u - rhs(u) seja 100% consistente com a lei de conservacao FV.
        """
        nr, nz, dx, dz = self.nr, self.nz, self.dx, self.dz
        if self.newtonian:
            return (
                np.ones((nr - 1, nz)),
                np.ones(nz),
                np.ones((nr, nz - 1)),
                np.ones(nr),
                np.ones(nr),
            )

        eps = 1e-14
        # Gradientes centrados nas celulas para componentes transversais nas faces
        dR_c = np.empty_like(u)
        dR_c[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2.0 * dx)
        dR_c[0, :] = (u[1, :] - u[0, :]) / (2.0 * dx)  # simetria em R=0 -> dR(0)=0, media na celula 0
        dR_c[-1, :] = (self.bc_vals[-1, :] - u[-2, :]) / (1.5 * dx)

        dZ_c = np.empty_like(u)
        dZ_c[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2.0 * dz)
        dZ_c[:, 0] = (u[:, 1] - self.bc_vals[:, 0]) / (1.5 * dz)
        dZ_c[:, -1] = (self.bc_vals[:, -1] - u[:, -2]) / (1.5 * dz)

        # Faces internas R (entre i e i+1, shape (nr-1, nz))
        dR_fR = (u[1:, :] - u[:-1, :]) / dx
        gZ_fR = 0.5 * (dZ_c[:-1, :] + dZ_c[1:, :]) + self.ub_z
        A_fR = np.sqrt(dR_fR ** 2 + gZ_fR ** 2 + eps ** 2)

        # Face contorno R=xmax (meia celula dx/2 ate a borda)
        dR_bR = (self.bc_vals[-1, :] - u[-1, :]) / (0.5 * dx)
        gZ_bR = dZ_c[-1, :] + self.ub_z
        A_bR = np.sqrt(dR_bR ** 2 + gZ_bR ** 2 + eps ** 2)

        # Faces internas Z (entre j e j+1, shape (nr, nz-1))
        dZ_fZ = (u[:, 1:] - u[:, :-1]) / dz + self.ub_z
        dR_fZ = 0.5 * (dR_c[:, :-1] + dR_c[:, 1:])
        A_fZ = np.sqrt(dR_fZ ** 2 + dZ_fZ ** 2 + eps ** 2)

        # Faces contorno Z=-xmax e Z=+xmax
        dZ_bZ0 = (u[:, 0] - self.bc_vals[:, 0]) / (0.5 * dz) + self.ub_z
        A_bZ0 = np.sqrt(dR_c[:, 0] ** 2 + dZ_bZ0 ** 2 + eps ** 2)

        dZ_bZ1 = (self.bc_vals[:, -1] - u[:, -1]) / (0.5 * dz) + self.ub_z
        A_bZ1 = np.sqrt(dR_c[:, -1] ** 2 + dZ_bZ1 ** 2 + eps ** 2)

        return A_fR, A_bR, A_fZ, A_bZ0, A_bZ1

    def assemble(self, u):
        nr, nz, dx, dz, R = self.nr, self.nz, self.dx, self.dz, self.R
        A_fR, A_bR, A_fZ, A_bZ0, A_bZ1 = self._face_amplitudes(u)

        I = np.arange(nr * nz).reshape(nr, nz)
        rows, cols, vals = [], [], []
        diag = np.zeros(nr * nz)
        rhs = (self.kappa * self.rho_nd).copy()

        # Faces R internas (i <-> i+1)
        c_out = A_fR * self.Rf[1:nr, None] / (R[:-1, None] * dx * dx)
        rows.append(I[:-1, :].ravel()); cols.append(I[1:, :].ravel()); vals.append((-c_out).ravel())
        np.add.at(diag, I[:-1, :].ravel(), c_out.ravel())

        c_in = A_fR * self.Rf[1:nr, None] / (R[1:, None] * dx * dx)
        rows.append(I[1:, :].ravel()); cols.append(I[:-1, :].ravel()); vals.append((-c_in).ravel())
        np.add.at(diag, I[1:, :].ravel(), c_in.ravel())

        # Faces Z internas (j <-> j+1)
        cz = A_fZ / (dz * dz)
        rows.append(I[:, :-1].ravel()); cols.append(I[:, 1:].ravel()); vals.append((-cz).ravel())
        np.add.at(diag, I[:, :-1].ravel(), cz.ravel())
        rows.append(I[:, 1:].ravel()); cols.append(I[:, :-1].ravel()); vals.append((-cz).ravel())
        np.add.at(diag, I[:, 1:].ravel(), cz.ravel())

        # Contorno Dirichlet R = xmax (distancia dx/2)
        cRb = A_bR * self.Rf[-1] / (R[-1] * dx * dx) * 2.0
        rhs[-1, :] += cRb * self.bc_vals[-1, :]
        np.add.at(diag, I[-1, :].ravel(), cRb.ravel())

        # Contorno Dirichlet Z = -xmax e Z = +xmax (distancia dz/2)
        cZ0 = A_bZ0 / (dz * dz) * 2.0
        rhs[:, 0] += cZ0 * self.bc_vals[:, 0]
        np.add.at(diag, I[:, 0].ravel(), cZ0.ravel())

        cZ1 = A_bZ1 / (dz * dz) * 2.0
        rhs[:, -1] += cZ1 * self.bc_vals[:, -1]
        np.add.at(diag, I[:, -1].ravel(), cZ1.ravel())

        rows.append(I.ravel()); cols.append(I.ravel()); vals.append(diag.ravel())
        M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nr * nz, nr * nz))
        return M, rhs.ravel()

    def nonlinear_residual(self, u):
        """Residuo discreto exato da PDE de Volumes Finitos: r(u) = M(u)*u - rhs(u)."""
        M, rhs = self.assemble(u)
        return (M @ u.ravel() - rhs).reshape(self.nr, self.nz)

    def solve(self, u0, max_iters=80, damp=0.65, tol=1e-11):
        u = u0.copy()
        init_res = float(np.max(np.abs(self.nonlinear_residual(u))))
        if self.newtonian:
            M, rhs = self.assemble(u)
            u = spla.spsolve(M.tocsc(), rhs).reshape(self.nr, self.nz)
            fres = float(np.max(np.abs(self.nonlinear_residual(u))))
            return u, init_res, fres, 1, 0.0

        du_last = 0.0
        for it in range(1, max_iters + 1):
            M, rhs = self.assemble(u)
            u_star = spla.spsolve(M.tocsc(), rhs).reshape(self.nr, self.nz)
            du = float(np.max(np.abs(u_star - u)))
            du_last = du
            u = (1.0 - damp) * u + damp * u_star
            res_now = float(np.max(np.abs(self.nonlinear_residual(u))))
            if res_now / max(init_res, 1e-30) < 1e-9 and du < tol * max(1.0, float(np.max(np.abs(u)))):
                return u, init_res, res_now, it, du_last
        fres = float(np.max(np.abs(self.nonlinear_residual(u))))
        return u, init_res, fres, max_iters, du_last

    def extract_aself(self, u, lo=0.85, hi=1.15):
        """Extract |grad u_self| in both the wide ring [0.85, 1.15] and thin shell [0.96, 1.04]."""
        dx, dz = self.dx, self.dz
        dR = np.zeros_like(u)
        dR[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2.0 * dx)
        dR[0, :] = (u[1, :] - u[0, :]) / (2.0 * dx)
        dR[-1, :] = (self.bc_vals[-1, :] - u[-2, :]) / (1.5 * dx)

        dZ = np.zeros_like(u)
        dZ[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2.0 * dz)
        dZ[:, 0] = (u[:, 1] - self.bc_vals[:, 0]) / (1.5 * dz)
        dZ[:, -1] = (self.bc_vals[:, -1] - u[:, -2]) / (1.5 * dz)

        mag = np.sqrt(dR ** 2 + dZ ** 2)
        ring_wide = (self.x > lo) & (self.x < hi)
        ring_thin = (self.x > 0.96) & (self.x < 1.04)
        rms_wide = float(np.sqrt(np.mean(mag[ring_wide] ** 2)))
        rms_thin = float(np.sqrt(np.mean(mag[ring_thin] ** 2)))
        return rms_wide, rms_thin


if __name__ == "__main__":
    t0 = time.time()
    kappa = 20.0
    xmax = 25.0
    results = {}

    print("==========================================================================")
    print("STAGE 0 — DIAGNOSTIC REPRODUCTION OF UNNORMALISED PLUMMER SOURCE (1/sqrt(3))")
    print("==========================================================================")
    g_bug = ValidatedCylindricalPDESolver(200, 400, xmax, kappa, plume_factor=1.0)
    g_bug.bc_vals = 2.0 * math.sqrt(kappa) * np.log(np.maximum(g_bug.x, 1e-6) / xmax)
    u0_bug = math.sqrt(kappa + 1.0) * 0.7 / np.sqrt(1.0 + g_bug.x ** 2)
    u_b, ir_b, fr_b, it_b, _ = g_bug.solve(u0_bug, max_iters=60, damp=0.5)
    aw_b, at_b = g_bug.extract_aself(u_b)
    target_nl = math.sqrt(kappa * (1.0 / (2.0 ** 1.5)))  # 2.659148
    print(f"[Unnormalised source diagnostic] a_wide={aw_b:.4f} | target={target_nl:.4f} | ratio={aw_b/target_nl:.4f} (= 1/sqrt(3) = {1/math.sqrt(3):.4f})")

    print("\n==========================================================================")
    print("STAGE 1 — NEWTONIAN POISSON MACHINE CHECK (-nabla^2 u = 3*kappa*(1+x^2)^(-5/2))")
    print("==========================================================================")
    target_newt = kappa / (2.0 ** 1.5)
    g_newt = ValidatedCylindricalPDESolver(200, 400, xmax, kappa, newtonian=True, plume_factor=3.0)
    g_newt.bc_vals = kappa / np.sqrt(1.0 + g_newt.x ** 2)
    u_n, ir_n, fr_n, _, _ = g_newt.solve(np.zeros_like(g_newt.x))
    aw_n, at_n = g_newt.extract_aself(u_n)
    print(f"[Newtonian 200x400] a_wide={aw_n:.5f} (ratio={aw_n/target_newt:.5f}) | a_x1={at_n:.5f} (ratio={at_n/target_newt:.5f}) | final_res={fr_n:.2e}")
    results["newtonian_200x400"] = {"target": target_newt, "a_wide": aw_n, "ratio_wide": aw_n / target_newt, "a_thin": at_n, "ratio_thin": at_n / target_newt, "final_res": fr_n}

    print("\n==========================================================================")
    print("STAGE 2 — ISOLATED NONLINEAR PDE v0 CONVERGENCE SUITE (200x400 -> 400x800)")
    print("Analytical Gauss target at x=1: sqrt(20 / 2^1.5) = 2.659148 | Tolerance: 1.00 +- 0.02")
    print("==========================================================================")
    base_wide, base_thin = None, None
    for nr, nz, tag in [(200, 400, "malha_base_200x400"), (400, 800, "malha_2x_400x800")]:
        ts = time.time()
        g = ValidatedCylindricalPDESolver(nr, nz, xmax, kappa, ub_z=0.0, newtonian=False, plume_factor=3.0)
        g.bc_vals = -math.sqrt(kappa) * np.log(np.maximum(g.x, 1e-6) / xmax)
        u0 = -math.sqrt(kappa) * np.log(np.sqrt(1.0 + g.x ** 2) / xmax)
        u_sol, ires, fres, n_it, du_f = g.solve(u0, max_iters=80, damp=0.7)
        aw, at = g.extract_aself(u_sol)
        drop_orders = math.log10(max(ires, 1e-30) / max(fres, 1e-30))
        dt = time.time() - ts
        print(f"[{tag}] ({dt:.1f}s, {n_it} iters):")
        print(f"  Wide ring [0.85, 1.15] : a_self = {aw:.5f} | target = {target_nl:.5f} | RATIO = {aw/target_nl:.5f}")
        print(f"  Thin shell [0.96, 1.04]: a_self = {at:.5f} | target = {target_nl:.5f} | RATIO = {at/target_nl:.5f}")
        print(f"  Residual {ires:.3e} -> {fres:.3e} | Drop = {drop_orders:.2f} orders (criterion >= 6.0)")
        results[tag] = {
            "a_wide": aw,
            "ratio_wide": aw / target_nl,
            "a_thin": at,
            "ratio_thin": at / target_nl,
            "init_res": ires,
            "final_res": fres,
            "drop_orders": drop_orders,
            "iters": n_it,
            "seconds": dt,
        }
        if base_wide is None:
            base_wide, base_thin = aw, at
        else:
            diff_w = abs(aw - base_wide) / base_wide * 100.0
            diff_t = abs(at - base_thin) / base_thin * 100.0
            print(f"  2x mesh drift: wide_ring = {diff_w:.4f}% | thin_shell = {diff_t:.4f}% (criterion < 2.0%)")
            results["mesh_refinement_pct"] = {"wide_pct": diff_w, "thin_pct": diff_t}

    passed = (
        abs(results["malha_base_200x400"]["ratio_wide"] - 1.0) <= 0.02
        and abs(results["malha_2x_400x800"]["ratio_wide"] - 1.0) <= 0.02
        and results["malha_base_200x400"]["drop_orders"] >= 6.0
        and results["malha_2x_400x800"]["drop_orders"] >= 6.0
        and results["mesh_refinement_pct"]["wide_pct"] < 2.0
    )
    results["OFFICIAL_VERDICT"] = "PASSED_ALL_CONVERGENCE_CRITERIA" if passed else "FAILED"
    print("\n==========================================================================")
    print(f"VERIFICATION VERDICT: {results['OFFICIAL_VERDICT']}")
    print("==========================================================================")

    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_path = os.path.join(pkg_root, "data", "pde_v0_validation_results.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"Saved to: {out_path}")
