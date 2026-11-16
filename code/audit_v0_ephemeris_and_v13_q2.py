#!/usr/bin/env python3
"""
Pre-Gaia DR4 Solar System Ephemeris Audit of Unscreened Additive v0 and
External-Field Tidal Quadrupole / Sky-Projected Anisotropy Audit of Screened v1.3
=================================================================================
Author: Jaqueline Martins (ORCID: 0009-0009-7732-5915)
Licence: MIT

Reproduces all Pre-DR4 Version 1.1.0 Amendment results:
  1. Positive Control: Symplectic Velocity-Verlet differential orbital integration
     reproduces General Relativity secular perihelion advance of Mercury (+42.98"/cy).
  2. Falsification of Unscreened Additive v0 (a_c = 1.9e-11 m s^-2):
     Exact Gauss perturbation quadrature + independent symplectic orbital integration
     for Mercury, Lunar Laser Ranging (LLR), LAGEOS-2, and GPS Block III.
  3. Exponentially Screened v1.3 (mu(x) = 1 - exp(-x), Milgrom 1999; a0 = 1.20e-10 m s^-2):
     - Local strong-field exponential suppression exp(-g_N/a0) = 0 in IEEE-754 float64;
     - Non-local Milky Way external-field phantom-source tidal quadrupole Q_2 at Saturn
       compared against the 2-sigma Cassini radio-tracking bound (Hees et al. 2014);
     - Exact 3D PDE sky-projected wide-binary quadrupole Delta_gamma in 10-30 kAU.
"""

import json
import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pde_v13_poisson_aqual_solver import PoissonAQUALPDESolver, solve_ub_from_gN

G_SI = 6.67430e-11
C_LIGHT = 299792458.0
M_SUN = 1.98847e30
M_EARTH = 5.97217e24
M_MOON = 7.3459e22
AU_M = 1.495978707e11
KAU_M = 1.0e3 * AU_M
ARCSEC_PER_RAD = 206264.80624709636
SEC_PER_JULIAN_YEAR = 365.25 * 86400.0
SEC_PER_CENTURY = 100.0 * SEC_PER_JULIAN_YEAR

A_C_V0 = 1.9e-11
A0_CANON = 1.20e-10


def gauss_secular_precession_v0(mu_gm: float, a_m: float, e: float, ac: float = A_C_V0):
    """
    Exact secular apsidal precession rate d(varpi)/dt for an inward radial perturbation
    delta_g_r(r) = sqrt(ac * mu_gm) / r (with outward radial force R_pert = -delta_g_r(r))
    evaluated via Gauss's planetary equation over true anomaly f in [0, 2*pi]:
      Delta_varpi_per_orbit = (sqrt(1 - e^2) / (e * a * n^2)) * int_0^{2*pi} R_pert(f) cos(f) (dt/df) df.
    """
    n_mean = math.sqrt(mu_gm / (a_m ** 3))
    period_s = 2.0 * math.pi / n_mean

    n_quad = 20000
    f = np.linspace(0.0, 2.0 * math.pi, n_quad, endpoint=False)
    df = 2.0 * math.pi / n_quad
    one_minus_e2 = 1.0 - e * e
    r_f = a_m * one_minus_e2 / (1.0 + e * np.cos(f))
    dt_df = (r_f ** 2) / (n_mean * (a_m ** 2) * math.sqrt(one_minus_e2))
    # Outward radial perturbing acceleration R_pert = -sqrt(ac * mu_gm) / r
    R_pert = -math.sqrt(ac * mu_gm) / r_f

    # Exact Gauss planetary equation: d(varpi)/dt = - (sqrt(1 - e^2) / (n * a * e)) * R_outward * cos(f)
    dvarpi_per_orbit_rad = (
        -(math.sqrt(one_minus_e2) / (e * a_m * n_mean))
        * float(np.sum(R_pert * np.cos(f) * dt_df) * df)
    )
    arcsec_per_yr = dvarpi_per_orbit_rad * ARCSEC_PER_RAD * (SEC_PER_JULIAN_YEAR / period_s)
    arcsec_per_cy = arcsec_per_yr * 100.0
    return {
        "period_days": period_s / 86400.0,
        "rad_per_orbit": dvarpi_per_orbit_rad,
        "arcsec_per_yr": arcsec_per_yr,
        "arcsec_per_century": arcsec_per_cy,
    }


def verlet_differential_precession(mu_gm: float, a_m: float, e: float, mode: str, steps_per_orbit: int = 80000, n_orbits: int = 25):
    """
    Symplectic Velocity-Verlet differential integration (with vs without perturbation)
    measuring Laplace-Runge-Lenz vector orientation to eliminate numerical phase drift.
    """
    period_s = 2.0 * math.pi * math.sqrt((a_m ** 3) / mu_gm)
    dt = period_s / steps_per_orbit
    total_steps = steps_per_orbit * n_orbits

    def integrate_run(include_pert: bool) -> float:
        x = a_m * (1.0 - e)
        y = 0.0
        vx = 0.0
        vy = math.sqrt(mu_gm * (1.0 + e) / (a_m * (1.0 - e)))

        def accel(px: float, py: float, pvx: float, pvy: float):
            r2 = px * px + py * py
            r = math.sqrt(r2)
            ax = -mu_gm * px / (r2 * r)
            ay = -mu_gm * py / (r2 * r)
            if include_pert:
                if mode == "v0":
                    dg = math.sqrt(A_C_V0 * mu_gm) / r
                    ax -= dg * (px / r)
                    ay -= dg * (py / r)
                elif mode == "GR":
                    # Exact 1PN Schwarzschild radial perturbation producing +6*pi*GM/(a*(1-e^2)*c^2)
                    h_z = px * pvy - py * pvx
                    dg_gr = 3.0 * mu_gm * (h_z * h_z) / ((C_LIGHT ** 2) * (r ** 4))
                    ax -= dg_gr * (px / r)
                    ay -= dg_gr * (py / r)
            return ax, ay

        ax, ay = accel(x, y, vx, vy)
        for _ in range(total_steps):
            vx_half = vx + 0.5 * dt * ax
            vy_half = vy + 0.5 * dt * ay
            x_new = x + dt * vx_half
            y_new = y + dt * vy_half
            ax_new, ay_new = accel(x_new, y_new, vx_half, vy_half)
            vx = vx_half + 0.5 * dt * ax_new
            vy = vy_half + 0.5 * dt * ay_new
            x, y, ax, ay = x_new, y_new, ax_new, ay_new

        h_z = x * vy - y * vx
        r = math.sqrt(x * x + y * y)
        ex = (vy * h_z) / mu_gm - x / r
        ey = (-vx * h_z) / mu_gm - y / r
        return math.atan2(ey, ex)

    ang_pert = integrate_run(True)
    ang_kepler = integrate_run(False)
    diff_rad = (ang_pert - ang_kepler + math.pi) % (2.0 * math.pi) - math.pi
    rad_per_orbit = diff_rad / n_orbits
    arcsec_per_cy = rad_per_orbit * ARCSEC_PER_RAD * (SEC_PER_CENTURY / period_s)
    return {"rad_per_orbit": rad_per_orbit, "arcsec_per_century": arcsec_per_cy}


def compute_v13_saturn_quadrupole_q2():
    """
    Computes the external-field-induced tidal quadrupole Q_2 at the Sun for v1.3
    via explicit 2D spherical integration of the QUMOND phantom density perturbation
    around r_M = sqrt(G M_sun / a0) = 7,030 AU.
    In QUMOND, nu(y) satisfies y * nu(y) = x(y) where (1 - exp(-x)) * x = y.
    Inside r << r_M, the external Milky Way field e_N = g_N,MW / a0 induces an axisymmetric
    quadrupole potential delta_Phi(r, theta) = -(1/2) * Q_2 * r^2 * (cos^2(theta) - 1/3).
    """
    r_m = math.sqrt(G_SI * M_SUN / A0_CANON)

    def x_of_y(y_arr: np.ndarray) -> np.ndarray:
        # Newton iteration on (1 - exp(-x))*x = y (exact to 1e-14 in 12 iterations)
        x = np.where(y_arr > 1.0, y_arr, np.sqrt(np.maximum(y_arr, 1e-30)))
        for _ in range(18):
            ex = np.exp(-np.minimum(x, 500.0))
            f = (1.0 - ex) * x - y_arr
            fp = (1.0 - ex) + x * ex
            x = np.maximum(x - f / np.maximum(fp, 1e-30), 1e-30)
        return x

    def nu_minus_1(y_arr: np.ndarray) -> np.ndarray:
        x = x_of_y(y_arr)
        return (x / np.maximum(y_arr, 1e-30)) - 1.0

    # Dimensionless QUMOND quadrupole coefficient q(e_N) defined by:
    # Q_2 = -(3 * a0 / (2 * r_M)) * int_0^infty du / u^2 int_{-1}^1 d(cos th) P_2(cos th) [...]
    # Evaluated directly via Poisson potential second derivative at the origin from div(g_phantom):
    # Phi_p(r) = -(1 / 4*pi) int d^3r' [ (nu(|y|) - 1) y(r') - (nu_e - 1) e_vec ] . (r' - r) / |r' - r|^3
    # Differentiating twice along z vs x at r = 0 gives the exact trace-free quadrupole Q_2:
    # Q_2 = -(a0 / r_M) * (3/4) * int_0^infty du / u^2 int_{-1}^1 dc [ ((nu - 1)*y_z - (nu_e - 1)*e_N)*(3 c^2 - 1) ... ]
    # Using Gauss's theorem / Desmond, Hees & Famaey (2024, MNRAS 530, 1781, Eq. 6):
    # Q_2 = (3 * G * M_sun / (2 * r_M^3)) * int_0^infty v^2 dv int_{-1}^1 dxi (nu(|y|) - nu_e) * ...
    # Let us evaluate Phi_p(r) numerically at small sphere r_test = 10 AU (well inside r_M = 7030 AU)
    # on a high-resolution (u, cos_theta) grid to avoid any sign or prefactor ambiguity!
    out = {}
    u_grid = np.geomspace(1.0e-3, 120.0, 900)
    du = np.gradient(u_grid)
    nc = 360
    c_grid = np.linspace(-1.0 + 1.0 / nc, 1.0 - 1.0 / nc, nc)
    dc = 2.0 / nc
    UU, CC = np.meshgrid(u_grid, c_grid, indexing="ij")
    SS = np.sqrt(np.maximum(0.0, 1.0 - CC * CC))

    for label, e_N in [("Solar_baryonic_gN_1p2e-10", 1.00), ("Solar_effective_eN_1p58", 1.9e-10 / A0_CANON)]:
        # Dimensionless Newtonian field y_vec = (-u^-2 * sin th, e_N - u^-2 * cos th)
        y_R = -(1.0 / (UU ** 2)) * SS
        y_Z = e_N - (1.0 / (UU ** 2)) * CC
        y_mag = np.sqrt(y_R ** 2 + y_Z ** 2)
        nm1 = nu_minus_1(y_mag)
        nm1_e = float(nu_minus_1(np.array([e_N]))[0])

        # Phantom extra field delta_g_p / a0 = (nu(|y|) - 1)*y_vec - (nu(e_N) - 1)*e_N*z_hat
        #Wait: in QUMOND, -div(g_tot) = -div(nu(|g_N|/a0) g_N) = 4*pi*G*(rho + rho_p),
        # so rho_p = -(1 / 4*pi*G) div( (nu - 1) g_N ).
        # Integrating by parts with boundary term vanishing after subtracting constant ambient (nu_e - 1)*g_Ne:
        # Phi_p(x) = -(1 / 4*pi) int d^3r' [ dF(r') . (x - r') / |x - r'|^3 ],
        # where dF(r') = (nu(|g_N|) - 1)*g_N(r') - (nu_e - 1)*g_Ne (pointing in direction of g = -grad Phi).
        # Let us double-check sign with a pure Newtonian scale factor nu - 1 = const = C > 0:
        # For point mass at origin, g_N = -G M r' / r'^3, so dF = -C G M r' / r'^3.
        # At x = 0, Phi_p -> -C G M / r < 0; specifically Phi_p(x) - Phi_p(0):
        # -(1 / 4*pi) int d^3r' dF(r') . [ (x - r')/|x - r'|^3 + r'/r'^3 ].
        dF_R = nm1 * y_R
        dF_Z = nm1 * y_Z - nm1_e * e_N

        # Exact tidal tensor T_ij = -d^2 Phi_p / (dx_i dx_j) at origin x = 0:
        # Differentiating K_k(x, r') = (x_k - r'_k)/|x - r'|^3 twice wrt x_z at x = 0:
        # For axisymmetric field around Z, trace-free tidal parameter Q_2 (Hees et al. 2014 convention:
        # delta_Phi = -(1/2) Q_2 r^2 (cos^2 th - 1/3) => d^2(delta_Phi)/dZ^2 - d^2(delta_Phi)/dX^2 = -Q_2):
        # Using exact multipole kernel at origin:
        # d^2/dZ^2 [(x - r')/|x - r'|^3]|_{x=0} - d^2/dX^2 [...]|_{x=0}:
        # Let n = r'/r' = (S cos phi, S sin phi, C). Azimuthally averaged kernel on dF_R and dF_Z:
        # Start from scalar potential of dipole density P = -dF / (4*pi):
        # Phi_p(x) = (1 / 4*pi) int d^3r' dF(r') . grad_{r'} (1 / |x - r'|).
        # Expand 1/|x - r'| for |x| << r': 1/r' + (x.r')/r'^3 + (1/2) r^2 P_2(cos th_x) * 2 P_2(cos th') / r'^3 ...
        # Wait: 1/|x - r'| = sum_l (r^l / r'^{l+1}) P_l(x_hat . r'_hat).
        # By addition theorem, azimuthally symmetric part for l = 2 is:
        # (r^2 / r'^3) P_2(cos th_x) P_2(cos th').
        # Note P_2(cos th_x) = (3/2)(cos^2 th_x - 1/3), so r^2 P_2(cos th_x) = (3/2) r^2 (cos^2 th_x - 1/3).
        # Comparing with delta_Phi = -(1/2) Q_2 r^2 (cos^2 th_x - 1/3) = -(1/3) Q_2 r^2 P_2(cos th_x):
        # Coefficient of r^2 P_2(cos th_x) in Phi_p(x) is:
        # C_2 = (1 / 4*pi) int d^3r' dF(r') . grad_{r'} [ P_2(cos th') / r'^3 ]!
        # Let's compute grad_{r'} [ P_2(c) / r'^3 ] in spherical/cylindrical coordinates:
        # Let f(r', th') = P_2(cos th') / r'^3 = (3 c^2 - 1) / (2 r'^3).
        # Radial derivative: df/dr' = -3 P_2(c) / r'^4 = -(3/2)(3 c^2 - 1) / r'^4.
        # Polar derivative: (1/r') df/dth' = (1/r'^4) * (3 c)(-sin th') = -3 c s / r'^4.
        # Converting dF = (dF_R, dF_Z) to spherical (dF_r, dF_th):
        # dF_r = dF_R * S + dF_Z * C
        # dF_th = dF_R * C - dF_Z * S
        # Therefore dF . grad(P_2 / r'^3) = [ -3 P_2(c) * dF_r - 3 c S * dF_th ] / r'^4!
        # Volume element d^3r' = 2*pi * r'^2 dr' dc. With r' = r_M * u and physical dF_phys = a0 * dF:
        # C_2 = (a0 / (2 * r_M)) * int_0^infty (du / u^2) int_{-1}^1 dc [ -3 P_2(c) dF_r - 3 c S dF_th ].
        # And since delta_Phi = C_2 * r^2 P_2(cos th_x) = -(1/3) Q_2 * r^2 P_2(cos th_x), we have:
        # Q_2 = -3 * C_2!
        dF_r = dF_R * SS + dF_Z * CC
        dF_th = dF_R * CC - dF_Z * SS
        P2 = 0.5 * (3.0 * CC * CC - 1.0)
        integrand = -3.0 * P2 * dF_r - 3.0 * CC * SS * dF_th
        integral_val = float(np.sum((integrand / (UU ** 2)) * du[:, None] * dc))
        C_2 = (A0_CANON / (2.0 * r_m)) * integral_val
        Q_2 = -3.0 * C_2
        out[label] = {
            "e_N": e_N,
            "r_M_AU": r_m / AU_M,
            "Q2_s_minus_2": Q_2,
            "abs_Q2_over_Cassini_2sigma": abs(Q_2) / 5.8e-27,
        }
    return out


def compute_v13_wide_binary_quadrupole():
    """
    Computes the sky-projected wide-binary quadrupole Delta_gamma = gamma_perp - gamma_par
    in 10-30 kAU (s = 15, 20, 25 kAU) for v1.3 using the identical 200,000-ray Monte Carlo
    sky-projection protocol as code/pde_v0_efe_dwarfs_binaries.py.
    """
    out = {}
    rng = np.random.default_rng(42)
    n_mc = 200000
    z_gc = rng.uniform(-1.0, 1.0, n_mc)
    beta = rng.uniform(0.0, 2.0 * math.pi, n_mc)
    y_perp = np.sqrt(1.0 - z_gc ** 2) * np.cos(beta)
    th_sky = np.arctan2(np.abs(y_perp), np.abs(z_gc))

    for env_lbl, ub_wb in [
        ("Solar_gtot_1p9e-10", 1.9e-10 / A0_CANON),
        ("Solar_gN_bar_1p2e-10", solve_ub_from_gN(1.2e-10)),
    ]:
        rows = {}
        for s_kau in [1.5, 15.0, 20.0, 25.0]:
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
            gam_sph = float(np.sum((gtot_r[sh] / gN_r[sh]) * ws) / np.sum(ws))

            cos_th = np.abs(sv.ZZ[sh]) / sv.x[sh]
            order = np.argsort(cos_th)
            gam_mc = np.interp(np.abs(z_gc), cos_th[order], (gtot_r[sh] / gN_r[sh])[order])
            gam_par = float(np.mean(gam_mc[th_sky < math.pi / 4.0]))
            gam_perp = float(np.mean(gam_mc[th_sky >= math.pi / 4.0]))
            rows[str(s_kau)] = {
                "gamma_sph": gam_sph,
                "gamma_sky_par": gam_par,
                "gamma_sky_perp": gam_perp,
                "Delta_gamma_perp_minus_par": gam_perp - gam_par,
            }
        g_ref = rows["1.5"]["gamma_sph"]
        mean_10_30_rel = float(np.mean([rows[k]["gamma_sph"] / g_ref for k in ("15.0", "20.0", "25.0")]))
        mean_10_30_quad = float(np.mean([rows[k]["Delta_gamma_perp_minus_par"] / g_ref for k in ("15.0", "20.0", "25.0")]))
        out[env_lbl] = {
            "bins": rows,
            "gamma_rel_10_30_mean": mean_10_30_rel,
            "Delta_gamma_rel_10_30_mean": mean_10_30_quad,
        }
    return out


if __name__ == "__main__":
    mu_sun = G_SI * M_SUN
    mu_earth_moon = G_SI * (M_EARTH + M_MOON)
    mu_earth = G_SI * M_EARTH

    gr_mercury = verlet_differential_precession(mu_sun, 0.38709893 * AU_M, 0.20563069, "GR", steps_per_orbit=30000, n_orbits=10)
    v0_mercury_gauss = gauss_secular_precession_v0(mu_sun, 0.38709893 * AU_M, 0.20563069)
    v0_mercury_verlet = verlet_differential_precession(mu_sun, 0.38709893 * AU_M, 0.20563069, "v0", steps_per_orbit=30000, n_orbits=10)

    v0_moon_gauss = gauss_secular_precession_v0(mu_earth_moon, 3.84399e8, 0.05490)
    v0_lageos2_gauss = gauss_secular_precession_v0(mu_earth, 1.2163e7, 0.0135)

    a_gps = 2.656e7
    gN_gps = mu_earth / (a_gps ** 2)
    frac_gps = math.sqrt(A_C_V0 / gN_gps)
    v_gps = math.sqrt(mu_earth / a_gps)
    # Uncompensated radial offset delta_g/g in fixed orbital radius shifts period/mean motion by (1/2)*(delta_g/g)
    gps_drift_m_per_day = 0.5 * frac_gps * v_gps * 86400.0

    q2_saturn = compute_v13_saturn_quadrupole_q2()
    wb_v13_quad = compute_v13_wide_binary_quadrupole()

    results = {
        "GR_Mercury_positive_control": gr_mercury,
        "v0_ephemeris_falsification": {
            "Mercury": {
                "gauss_arcsec_per_century": v0_mercury_gauss["arcsec_per_century"],
                "verlet_arcsec_per_century": v0_mercury_verlet["arcsec_per_century"],
                "obs_bound_arcsec_per_century": 0.6,
                "exclusion_factor": abs(v0_mercury_gauss["arcsec_per_century"]) / 0.6,
            },
            "Moon_LLR": {
                "gauss_arcsec_per_century": v0_moon_gauss["arcsec_per_century"],
                "gauss_arcsec_per_yr": v0_moon_gauss["arcsec_per_yr"],
            },
            "LAGEOS2_SLR": {
                "gauss_arcsec_per_yr": v0_lageos2_gauss["arcsec_per_yr"],
                "obs_bound_arcsec_per_yr": 0.01,
                "exclusion_factor": abs(v0_lageos2_gauss["arcsec_per_yr"]) / 0.01,
            },
            "GPS_Block_III": {
                "semi_major_axis_m": a_gps,
                "delta_g_over_gN": frac_gps,
                "along_track_drift_m_per_day": gps_drift_m_per_day,
            },
        },
        "v13_local_screening_float64": {
            "Saturn_gN_over_a0": (mu_sun / ((9.5826 * AU_M) ** 2)) / A0_CANON,
            "exp_minus_x_Saturn": math.exp(-min(750.0, (mu_sun / ((9.5826 * AU_M) ** 2)) / A0_CANON)),
        },
        "v13_Saturn_QUMOND_quadrupole_Q2": q2_saturn,
        "v13_wide_binary_quadrupole": wb_v13_quad,
    }

    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_path = os.path.join(pkg_root, "data", "audit_v0_ephemeris_and_v13_q2_results.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)

    print("GR Mercury control (+42.98\"/cy):", f"{gr_mercury['arcsec_per_century']:+.2f}\"/cy")
    print("v0 Mercury Gauss vs Verlet:", f"{v0_mercury_gauss['arcsec_per_century']:.1f}\"/cy vs {v0_mercury_verlet['arcsec_per_century']:.1f}\"/cy")
    print("v0 Moon LLR:", f"{v0_moon_gauss['arcsec_per_century']:.1f}\"/cy")
    print("v0 LAGEOS-2:", f"{v0_lageos2_gauss['arcsec_per_yr']:.1f}\"/yr")
    print("v0 GPS drift:", f"{gps_drift_m_per_day:.1f} m/day (frac={frac_gps:.3e})")
    for k, v in q2_saturn.items():
        print(f"v1.3 Saturn Q2 ({k}): Q2 = {v['Q2_s_minus_2']:.3e} s^-2 ({v['abs_Q2_over_Cassini_2sigma']:.2f}x Cassini 2sigma)")
    for k, v in wb_v13_quad.items():
        print(f"v1.3 Wide Binary ({k}): gamma_rel={v['gamma_rel_10_30_mean']:.4f}, Delta_gamma={v['Delta_gamma_rel_10_30_mean']:+.4f}")
