# coding=utf-8
"""
##########################################################################################
#    Program:       DK-RD2 ThermoGravity Framework
#    Main script:   DK-RD2_ThermoGravity_Framework.py
#    Author:        Gabriel Martín del Campo Flores
#    Contact:       gabemdelc@gmail.com
#    Created:       11/Feb/2025
#    Revised for:   May/2026
#    License:       MIT License
##########################################################################################
#
#        What Is Gravity?
#        From Emergent Mass to Thermodynamic–Geometric Dynamics
#        in the DK-RD2 Framework
#
#    Description:
#    This script implements the DK-RD2 ThermoGravity framework as presented in the
#    current paper. Its primary goal is to test the falsifiable consequences of
#    gravity as an emergent thermodynamic–geometric phenomenon.
#
#    In this model, gravity arises from the thermodynamic organization and evolution
#    of energy, encoded through a dynamic relativistic coupling:
#
#        Gab_max(T, v) ≈ G0 · [1 + (v²/c²)(T0/T)]
#
#    while the observable large-scale gravitational interaction is represented by:
#
#        Gab_eff(T, v) ≈ G0 · [1 + α_DK · (v²/c²)(T0/T)]
#
#    where α_DK is an effective geometric projection factor calibrated directly
#    from DESI BAO and cosmic chronometer observations.
#
#    Scientific focus of this version:
#      • Emergent gravity from thermodynamic energy organization
#      • Dynamic gravitational coupling Gab(T, v)
#      • Effective geometric projection factor α_DK
#      • DESI BAO + CC calibration of α_DK and r_d
#      • Statistical χ² constraint on α_DK
#      • Diagnostic reconstruction of w(z) without fundamental dark energy
#      • Thermodynamic stabilization and gravitational confinement
#      • Effective-gravity interpretation of compact halo structures
#      • Observational falsifiability against ΛCDM
#
#    Code structure:
#      - DK_RD2_Core.py
#          Core physics:
#            • Gab(T, v)
#            • Expansion dynamics
#            • Distance relations
#            • BAO observables
#            • Statistical fitting tools
#
#      - DK_RD2_CLASS.py
#          Optional CMB / geometric consistency tools.
#
#      - DK-RD2_ThermoGravity_Framework.py
#          Main driver:
#            • Loads observational datasets
#            • Computes DK-RD2 predictions
#            • Generates figures for this paper
#            • Exports statistical tables and diagnostics
#
#    Core computational outputs:
#      ✓ Gab_max(T, v) thermodynamic–relativistic coupling
#      ✓ Gab_eff(T, v) projected observable coupling
#      ✓ DESI-compatible BAO expansion diagnostics
#      ✓ Calibrated α_DK and sound horizon r_d
#      ✓ χ² profiling of α_DK
#      ✓ Effective equation of state w(z)
#      ✓ Thermodynamic activity diagnostics
#      ✓ Compact-halo effective-gravity interpretation
#      ✓ Falsifiability comparison (DK-RD2 vs ΛCDM)
#
#    Figures generated for this paper:
#
#      1 — Figure 01:
#          DESI BAO calibration of α_DK and r_d.
#
#      2 — Figure 02:
#          Statistical constraint on the DK projection factor.
#
#      3 — Figure 03:
#          Thermodynamic–relativistic coupling structure.
#
#      4 — Figure 04:
#          Diagnostic reconstruction of the effective equation of state w(z).
#
#      5 — Figure 05:
#          Thermodynamic activity versus projected gravitational coupling.
#
#      6 — Figure 06:
#          Effective-gravity interpretation of compact halo structures.
#
#    Outputs:
#      - Figures exported through generate_evidence("image", figure_number)
#      - Tables exported through generate_evidence("table", figure_number)
#      - Statistical summaries saved as *_stats.csv
#      - Output paths defined in DK_RD2_Core.py
#
#    Falsifiability:
#      DK-RD2 predicts:
#        • w(z) ≠ -1 (non-constant evolution)
#        • evolving thermodynamic gravitational coupling
#        • statistically constrained α_DK
#        • observable deviations from ΛCDM expansion history
#        • thermodynamic signatures in compact gravitational systems
#
#      If observations confirm:
#        • w(z) = -1 within uncertainties
#        • α_DK → 0
#        • no statistically significant deviations from ΛCDM
#
#      → DK-RD2 is ruled out in its present form.
#
#    Scientific interpretation:
#      In DK-RD2, gravity is not treated as a fundamental interaction,
#      but as an emergent thermodynamic–geometric manifestation of
#      evolving energy organization across spacetime.
#
#    Motto:
#      Gravity is not assumed — it emerges.
#      GabE = mc²  —  Luludns = ∞Ψ
##########################################################################################
"""
from DK_RD2_Core import * # DK-RD2 Core Utilities – Constants, Functions, and Relativistic Dynamic Gravitational Engine
##########################################################################################


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

Core_git_gabe = "https://github.com/gabemdelc/DK-RD2-What-Is-Gravity"
Core_zenodo = "Zenodo DOI: https://doi.org/10.5281/zenodo.20078175" # What is Gravity

Core_autor_text = (
    f"{Core_author}. "
    "Reproducible from public DK-RD2 code.\n"
    f"GitHub: {Core_git_gabe} | Zenodo: {Core_zenodo}"
)

# Global linestyle & color for use in all Figures to models comparison
DK_RD2_color      = "blue"     # DK-RD2 main curves
LCDM_color        = "orange"

DK_LIGHT_color    = "#7fe7ff"
LCDM_LIGHT_color  = "#ffd280"

DATA_color        = "white"
ERROR_color       = "0.75"

SECONDARY_color   = "#b57cff"

GRID_color        = "0.22"
TEXT_color        = "0.92"

DK_RD2_linestyle="dashdot"
LCDM_linestyle="--"

# Global DK-RD2 calibrated projection factor cache
alpha_dk_best_global = 0.0
rd_dk_best_global = 0.0
alpha_dk_source_global = "not computed"

def get_alpha_DK(
    desi_bao_mean_path: str,
    desi_bao_cov_path: str,
):
    """
    Return cached alpha_DK. If not available, compute it from DESI BAO + CC.
    """

    global alpha_dk_best_global
    global rd_dk_best_global
    global alpha_dk_source_global

    if alpha_dk_best_global > 0.0:
        return alpha_dk_best_global, rd_dk_best_global

    alpha_dk_best_global, rd_dk_best_global = compute_alpha_DK_from_DESI(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_dk_source_global = "DESI BAO + CC calibration"

    return alpha_dk_best_global, rd_dk_best_global

def generate_figure01(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    hz_data_path: str | None = None,
    Core_H0: float | None = None,
    Core_Omega_m_LCDM: float | None = None,
    Omega_L_LCDM: float | None = None,
    Omega_L_DK: float | None = None,
    rd_planck: float | None = None,
    rd_bounds: tuple[float, float] | None = None,
    alpha_bounds: tuple[float, float] | None = None,
    v_model=None,
    T_model=None,
):
    """
    # ============================================================
    # Figure 01 — DESI Expansion and BAO Calibration of α_DK
    # ============================================================
    #
    # Purpose
    # -------
    # Calibrate and visualize the DK-RD2 effective gravitational
    # projection factor α_DK using DESI BAO + cosmic chronometer data.
    #
    # Notes
    # -----
    # α_DK and r_d are obtained through get_alpha_DK().
    # If the global cached values are already available, they are reused.
    # If not, get_alpha_DK() computes them from DESI BAO + CC calibration.
    #
    # Returns
    # -------
    # file_fig, file_table, file_stats, alpha_dk_best
    # ============================================================
    """

    from scipy.optimize import minimize_scalar

    # ------------------------------------------------------------
    # Resolve defaults from DK_RD2_Core.py constants only.
    # ------------------------------------------------------------
    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if Core_Omega_m_LCDM is None:
        Core_Omega_m_LCDM = float(Core_OMEGA_M_LCDM)

    if Omega_L_LCDM is None:
        Omega_L_LCDM = float(Core_OMEGA_L_LCDM)

    if rd_planck is None:
        rd_planck = float(Core_rd_planck_mpc)

    if rd_bounds is None:
        rd_bounds = tuple(Core_rd_fit_bounds_mpc)

    if alpha_bounds is None:
        alpha_bounds = tuple(Core_alpha_DK_bounds)

    file_fig = generate_evidence("image", 1)
    file_table = generate_evidence("table", 1)
    file_stats = file_table.replace(".csv", "_stats.csv")

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure01 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Global/cache alpha_DK access
    # ------------------------------------------------------------
    alpha_DK, rd_DK = get_alpha_DK(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_dk_best = float(alpha_DK)
    rd_dk_best = float(rd_DK)

    alpha_footer = (
        rf"$\alpha_{{DK}}={alpha_dk_best:.4f}$, "
        rf"$r_d={rd_dk_best:.2f}\,\mathrm{{Mpc}}$ "
        r"(DESI BAO + CC calibration)"
    )

    # ------------------------------------------------------------
    # Load DESI BAO compressed data
    # ------------------------------------------------------------
    z_bao, bao_obs, bao_type, bao_cov, bao_cov_inv = load_desi_gaussian_bao(
        desi_bao_mean_path,
        desi_bao_cov_path,
    )

    sigma_bao = np.sqrt(np.diag(bao_cov))

    # ------------------------------------------------------------
    # DK calibrated-coupling background
    # ------------------------------------------------------------
    def E_DK_alpha(z_in, alpha_DK):
        z_arr = np.asarray(z_in, dtype=float)

        E_raw = E_Relativistic(
            z_arr,
            Core_Omega_m=None,
            Omega_L_value=None,
            v_model=v_model,
            T_model=T_model,
        )

        E2_alpha = 1.0 + float(alpha_DK) * (E_raw**2 - 1.0)
        return np.sqrt(np.clip(E2_alpha, 1e-300, None))

    def H_DK_alpha(z_in, alpha_DK):
        return float(Core_H0) * E_DK_alpha(z_in, alpha_DK)

    # ------------------------------------------------------------
    # LambdaCDM reference and LambdaCDM r_d fit
    # ------------------------------------------------------------
    # ------------------------------------------------------------
    # ΛCDM reference uses fixed standard cosmological parameters
    # defined globally in DK_RD2_Core.py.
    #
    # These parameters are NOT adjusted during the DK-RD2 fit.
    # DK-RD2 calibration only optimizes:
    #   - α_DK
    #   - r_d
    # ------------------------------------------------------------
    lcdm_fit = fit_rd_for_bao(
        z_bao,
        bao_obs,
        bao_type,
        bao_cov_inv,
        E_LCDM,
        Core_H0,
        rd_bounds=rd_bounds,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    bao_lcdm_planck = bao_distance_vector_over_rs(
        z_bao,
        bao_type,
        E_LCDM,
        Core_H0,
        rd_planck,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    chi2_lcdm_planck = chi2_gaussian_bao(
        bao_lcdm_planck,
        bao_obs,
        bao_cov_inv,
    )

    # ------------------------------------------------------------
    # DESI DR1 cosmic chronometer point
    # ------------------------------------------------------------
    z_cc = float(Core_DESI_CC_z)
    H_cc = float(Core_DESI_CC_H)
    H_cc_err = float(Core_DESI_CC_H_err)

    H_lcdm_cc = float(H_LCDM(
        z_cc,
        Core_H0,
        Core_Omega_m_LCDM,
        Omega_L_LCDM,
    ))

    chi2_cc_lcdm = chi2_cosmic_chronometer_point(
        H_lcdm_cc,
        H_obs=H_cc,
        sigma_H=H_cc_err,
    )

    # ------------------------------------------------------------
    # DK-RD2 model vector using cached/calibrated alpha_DK and r_d
    # ------------------------------------------------------------
    def dk_model_vector(rd_mpc, alpha_DK):
        return bao_distance_vector_over_rs(
            z_bao,
            bao_type,
            lambda zz: E_DK_alpha(zz, alpha_DK),
            Core_H0,
            rd_mpc,
        )

    bao_dk_best = dk_model_vector(rd_dk_best, alpha_dk_best)
    dk_residuals = bao_dk_best - bao_obs
    dk_pulls = dk_residuals / sigma_bao

    chi2_dk_bao = chi2_gaussian_bao(
        bao_dk_best,
        bao_obs,
        bao_cov_inv,
    )

    H_dk_cc = float(H_DK_alpha(z_cc, alpha_dk_best))

    chi2_cc_dk = chi2_cosmic_chronometer_point(
        H_dk_cc,
        H_obs=H_cc,
        sigma_H=H_cc_err,
    )

    # ------------------------------------------------------------
    # DK-RD2 diagnostic tests
    # ------------------------------------------------------------
    def dk_chi2_fixed(rd_mpc, alpha_DK):
        bao_model = dk_model_vector(rd_mpc, alpha_DK)

        chi2_bao = chi2_gaussian_bao(
            bao_model,
            bao_obs,
            bao_cov_inv,
        )

        H_model_cc = float(H_DK_alpha(z_cc, alpha_DK))

        chi2_cc = chi2_cosmic_chronometer_point(
            H_model_cc,
            H_obs=H_cc,
            sigma_H=H_cc_err,
        )

        return {
            "rd_mpc": float(rd_mpc),
            "alpha_DK": float(alpha_DK),
            "chi2_BAO": float(chi2_bao),
            "chi2_CC": float(chi2_cc),
            "chi2_total": float(chi2_bao + chi2_cc),
            "model_vector": bao_model,
        }

    # Test 1: alpha = 0, rd free
    def chi2_alpha0_rd_only(rd_mpc):
        return dk_chi2_fixed(rd_mpc, 0.0)["chi2_total"]

    res_alpha0 = minimize_scalar(
        chi2_alpha0_rd_only,
        bounds=rd_bounds,
        method="bounded",
    )

    test_alpha0 = dk_chi2_fixed(res_alpha0.x, 0.0)

    # Test 2: alpha = 1, rd free
    def chi2_alpha1_rd_only(rd_mpc):
        return dk_chi2_fixed(rd_mpc, 1.0)["chi2_total"]

    res_alpha1 = minimize_scalar(
        chi2_alpha1_rd_only,
        bounds=rd_bounds,
        method="bounded",
    )

    test_alpha1 = dk_chi2_fixed(res_alpha1.x, 1.0)

    # Test 3: rd fixed to reference, alpha free
    def chi2_rd_fixed_alpha_only(alpha_DK):
        return dk_chi2_fixed(rd_planck, alpha_DK)["chi2_total"]

    res_rd_fixed = minimize_scalar(
        chi2_rd_fixed_alpha_only,
        bounds=alpha_bounds,
        method="bounded",
    )

    test_rd_fixed = dk_chi2_fixed(rd_planck, res_rd_fixed.x)

    # ------------------------------------------------------------
    # H(z) curves
    # ------------------------------------------------------------
    z_plot = np.linspace(0.001, 2.5, 900)

    H_lcdm = H_LCDM(
        z_plot,
        Core_H0,
        Core_Omega_m_LCDM,
        Omega_L_LCDM,
    )

    H_dk = H_DK_alpha(z_plot, alpha_dk_best)

    # ------------------------------------------------------------
    # BAO visual curves
    # ------------------------------------------------------------
    z_bao_curve = np.linspace(0.05, 2.5, 600)

    DM_lcdm_curve = np.array([
        comoving_distance(
            zi,
            E_LCDM,
            Core_H0,
            Core_c_km_s=Core_c_km_s,
            Core_Omega_m=Core_Omega_m_LCDM,
            Core_Omega_L=Omega_L_LCDM,
        ) / lcdm_fit["best_rd"]
        for zi in z_bao_curve
    ])

    DH_lcdm_curve = np.array([
        (
            Core_c_km_s / H_LCDM(
                zi,
                Core_H0,
                Core_Omega_m_LCDM,
                Omega_L_LCDM,
            )
        ) / lcdm_fit["best_rd"]
        for zi in z_bao_curve
    ])

    DM_dk_curve = np.array([
        comoving_distance(
            zi,
            lambda zz: E_DK_alpha(zz, alpha_dk_best),
            Core_H0,
            Core_c_km_s=Core_c_km_s,
        ) / rd_dk_best
        for zi in z_bao_curve
    ])

    DH_dk_curve = np.array([
        (
            Core_c_km_s / H_DK_alpha(zi, alpha_dk_best)
        ) / rd_dk_best
        for zi in z_bao_curve
    ])

    # ------------------------------------------------------------
    # Save evidence table
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "z": z_bao,
        "observable": bao_type,
        "DESI_DR2_value": bao_obs,
        "sigma_diag": sigma_bao,
        "LCDM_reference_rd": bao_lcdm_planck,
        "LCDM_fit_rd": lcdm_fit["model"],
        "DKRD2_alpha_fit_rd": bao_dk_best,
        "LCDM_fit_residual": lcdm_fit["residuals"],
        "DKRD2_alpha_fit_residual": dk_residuals,
        "LCDM_fit_pull": lcdm_fit["pulls"],
        "DKRD2_alpha_fit_pull": dk_pulls,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="DESI_BAO_DR2_ALPHA_DK",
        figure_id=1,
        strict=False,
        index=False,
        meta={
            "rd_reference_mpc": rd_planck,
            "rd_lcdm_best_mpc": lcdm_fit["best_rd"],
            "rd_dkrd2_best_mpc": rd_dk_best,
            "alpha_dkrd2_best": alpha_dk_best,
            "chi2_lcdm_reference_bao": chi2_lcdm_planck,
            "chi2_lcdm_fit_bao": lcdm_fit["best_chi2"],
            "chi2_dkrd2_fit_bao": chi2_dk_bao,
            "chi2_dkrd2_fit_cc": chi2_cc_dk,
            "chi2_dkrd2_fit_total": chi2_dk_bao + chi2_cc_dk,
            "H_cc_z": z_cc,
            "H_cc": H_cc,
            "H_cc_err": H_cc_err,
        },
    )

    # ------------------------------------------------------------
    # Plot Figure 01
    # ------------------------------------------------------------
    fig, (ax_h, ax_bao) = plt.subplots(
        1,
        2,
        figsize=(15, 6.5),
        dpi=130,
    )

    # Panel A — H(z)
    ax_h.plot(
        z_plot,
        H_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.1,
        label=r"$\Lambda$CDM reference",
    )

    ax_h.plot(
        z_plot,
        H_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.3,
        label=rf"DK-RD2 calibrated coupling: $\alpha_{{\rm DK}}={alpha_dk_best:.4f}$",
    )

    ax_h.errorbar(
        [z_cc],
        [H_cc],
        yerr=[H_cc_err],
        fmt="o",
        color="black",
        capsize=4,
        label=rf"DESI DR1 CC: $H({z_cc:.2f})={H_cc:.2f}\pm{H_cc_err:.2f}$",
    )

    ax_h.text(
        0.03,
        0.96,
        r"DK-RD2 calibration test:" "\n"
        r"$E_{\rm DK,\alpha}^2(z)=1+\alpha_{\rm DK}[E_{\rm DK,raw}^2(z)-1]$" "\n"
        r"$\alpha_{\rm DK}$ and rd are jointly obtained from BAO + CC calibration.",
        transform=ax_h.transAxes,
        fontsize=8,
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="0.5",
            alpha=0.85,
        ),
    )

    ax_h.set_xlabel("Redshift z")
    ax_h.set_ylabel(r"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
    ax_h.set_title(r"Panel A — Expansion Rate")
    ax_h.grid(alpha=0.3)
    ax_h.legend(fontsize=8, loc="best")

    # Panel B — BAO observables
    ax_bao.plot(
        z_bao_curve,
        DM_lcdm_curve,
        color=LCDM_color,
        linestyle="-",
        linewidth=1.8,
        label=rf"$\Lambda$CDM $D_M/r_d$, best $r_d={lcdm_fit['best_rd']:.2f}$ Mpc",
    )

    ax_bao.plot(
        z_bao_curve,
        DH_lcdm_curve,
        color=LCDM_color,
        linestyle="--",
        linewidth=1.8,
        label=r"$\Lambda$CDM $D_H/r_d$",
    )

    ax_bao.plot(
        z_bao_curve,
        DM_dk_curve,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.0,
        label=rf"DK-RD2 $D_M/r_d$, $r_d={rd_dk_best:.2f}$ Mpc",
    )

    ax_bao.plot(
        z_bao_curve,
        DH_dk_curve,
        color=DK_RD2_color,
        linestyle="--",
        linewidth=2.0,
        label=r"DK-RD2 $D_H/r_d$",
    )

    mask_dm = bao_type == "DM_over_rs"
    mask_dh = bao_type == "DH_over_rs"
    mask_dv = bao_type == "DV_over_rs"

    ax_bao.errorbar(
        z_bao[mask_dm],
        bao_obs[mask_dm],
        yerr=sigma_bao[mask_dm],
        fmt="s",
        color="black",
        capsize=3,
        label=r"DESI DR2 $D_M/r_d$",
    )

    ax_bao.errorbar(
        z_bao[mask_dh],
        bao_obs[mask_dh],
        yerr=sigma_bao[mask_dh],
        fmt="^",
        color="0.35",
        capsize=3,
        label=r"DESI DR2 $D_H/r_d$",
    )

    ax_bao.errorbar(
        z_bao[mask_dv],
        bao_obs[mask_dv],
        yerr=sigma_bao[mask_dv],
        fmt="o",
        color="0.15",
        capsize=3,
        label=r"DESI DR2 $D_V/r_d$",
    )

    delta_chi2 = (chi2_dk_bao + chi2_cc_dk) - (lcdm_fit["best_chi2"] + chi2_cc_lcdm)

    fit_text = (
        r"DK-RD2 best fit:" "\n"
        rf"$\chi^2_{{\rm BAO,DK}} = {chi2_dk_bao:.2f}$" "\n"
        rf"$\chi^2_{{\rm CC,DK}} = {chi2_cc_dk:.3f}$" "\n"
        rf"$\chi^2_{{\rm total,DK}} = {(chi2_dk_bao + chi2_cc_dk):.2f}$" "\n"
        "\n"
        rf"$\chi^2_{{\Lambda CDM}} = {(lcdm_fit['best_chi2'] + chi2_cc_lcdm):.2f}$" "\n"
        rf"$\Delta \chi^2_{{\rm DK-\Lambda CDM}} = {delta_chi2:.2f}$"
    )

    ax_bao.text(
        0.70,
        0.77,
        fit_text,
        transform=ax_bao.transAxes,
        fontsize=9.0,
        va="top",
        ha="left",
        linespacing=1.25,
        bbox=dict(
            boxstyle="round,pad=0.55",
            facecolor="white",
            edgecolor="black",
            linewidth=1.1,
            alpha=0.92,
        ),
    )

    ax_bao.set_xlabel("Redshift z")
    ax_bao.set_ylabel(r"BAO compressed observable")
    ax_bao.set_title(r"Panel B — DESI DR2 Gaussian BAO")
    ax_bao.grid(alpha=0.3)
    ax_bao.legend(fontsize=7, loc="best")

    fig.suptitle(
        rf"Figure 01 — DESI Expansion and BAO Calibration of $\alpha_{{\rm DK}}$"
        "\n"
        rf"$\alpha_{{\rm DK}} = {alpha_dk_best:.4f},\quad r_d = {rd_dk_best:.2f}\,\mathrm{{Mpc}}$",
        fontsize=14,
    )
    """
    fig.text(
        0.5,
        0.085,
        alpha_footer,
        ha="center",
        fontsize=9,
        color=DK_RD2_color,
    )
    """
    fig.text(
        0.5,
        0.043,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=6,
        color=DK_RD2_color,
    )

    plt.tight_layout(rect=(0, 0.06, 1, 1.0))
    plt.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()

    # ------------------------------------------------------------
    # Stats CSV
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "ΛCDM reference-calibrated",
            "rd_mpc": float(rd_planck),
            "alpha_DK": np.nan,
            "chi2_BAO": float(chi2_lcdm_planck),
            "chi2_CC": float(chi2_cc_lcdm),
            "chi2_total": float(chi2_lcdm_planck + chi2_cc_lcdm),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "reference_rd_fixed",
        },
        {
            "model": "ΛCDM BAO-calibrated",
            "rd_mpc": float(lcdm_fit["best_rd"]),
            "alpha_DK": np.nan,
            "chi2_BAO": float(lcdm_fit["best_chi2"]),
            "chi2_CC": float(chi2_cc_lcdm),
            "chi2_total": float(lcdm_fit["best_chi2"] + chi2_cc_lcdm),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "rd_free",
        },
        {
            "model": "DK-RD2 calibrated coupling",
            "rd_mpc": float(rd_dk_best),
            "alpha_DK": float(alpha_dk_best),
            "chi2_BAO": float(chi2_dk_bao),
            "chi2_CC": float(chi2_cc_dk),
            "chi2_total": float(chi2_dk_bao + chi2_cc_dk),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "rd_plus_alpha_free",
        },
        {
            "model": "DK-RD2 diagnostic alpha=0",
            "rd_mpc": float(test_alpha0["rd_mpc"]),
            "alpha_DK": float(test_alpha0["alpha_DK"]),
            "chi2_BAO": float(test_alpha0["chi2_BAO"]),
            "chi2_CC": float(test_alpha0["chi2_CC"]),
            "chi2_total": float(test_alpha0["chi2_total"]),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "alpha_fixed_0_rd_free",
        },
        {
            "model": "DK-RD2 diagnostic alpha=1",
            "rd_mpc": float(test_alpha1["rd_mpc"]),
            "alpha_DK": float(test_alpha1["alpha_DK"]),
            "chi2_BAO": float(test_alpha1["chi2_BAO"]),
            "chi2_CC": float(test_alpha1["chi2_CC"]),
            "chi2_total": float(test_alpha1["chi2_total"]),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "alpha_fixed_1_rd_free",
        },
        {
            "model": "DK-RD2 diagnostic rd fixed",
            "rd_mpc": float(test_rd_fixed["rd_mpc"]),
            "alpha_DK": float(test_rd_fixed["alpha_DK"]),
            "chi2_BAO": float(test_rd_fixed["chi2_BAO"]),
            "chi2_CC": float(test_rd_fixed["chi2_CC"]),
            "chi2_total": float(test_rd_fixed["chi2_total"]),
            "N_BAO": int(len(bao_obs)),
            "N_CC": 1,
            "fit_mode": "rd_fixed_alpha_free",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=1,
        fit_mode="DESI_BAO_PLUS_CC_ALPHA_DK",
        index=False,
    )

    return file_fig, file_table, file_stats, alpha_dk_best

def generate_figure02(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
    *,
    n_alpha: int = 220,
    alpha_min: float | None = None,
    alpha_max: float | None = None,
    rd_bounds: tuple[float, float] | None = None,
    Core_H0: float | None = None,
    Omega_L_DK: float | None = None,
    v_model=None,
    T_model=None,
):
    """
    # ============================================================
    # Figure 02 — Statistical Constraint on the DK Projection Factor
    # ============================================================
    #
    # Purpose
    # -------
    # Demonstrate that alpha_DK is not arbitrary by profiling the
    # BAO+CC chi-square as a function of alpha_DK.
    #
    # Method
    # ------
    # For each alpha_DK:
    #   1. Build E_DK,alpha(z)
    #   2. Profile the best r_d by minimizing chi2_BAO+CC
    #   3. Compute Delta chi2(alpha_DK)
    #
    # The minimum emerges directly from the DESI BAO + CC likelihood
    # after profiling r_d for each alpha_DK value.
    #
    # Notes
    # -----
    # - Core_H0, rd_bounds and alpha bounds are resolved from DK_RD2_Core.py.
    # - DK-RD2 uses Omega_L_value=None in E_Relativistic, i.e. no explicit Lambda.
    # - The CC term is a single low-redshift chronometer anchor.
    # ============================================================
    """

    from scipy.optimize import minimize_scalar

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure02 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Resolve defaults from DK_RD2_Core.py constants only.
    # ------------------------------------------------------------
    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if rd_bounds is None:
        rd_bounds = tuple(Core_rd_fit_bounds_mpc)

    if alpha_min is None:
        alpha_min = float(Core_alpha_DK_bounds[0])

    if alpha_max is None:
        alpha_max = float(Core_alpha_DK_bounds[1])

    # DK-RD2: no explicit Lambda in this figure.
    # Kept only for API compatibility; the actual call below uses None.
    Omega_L_DK = None

    # ------------------------------------------------------------
    # Reference alpha_DK from global/cache calibration
    # ------------------------------------------------------------
    alpha_DK_best_global, rd_DK_best_global = get_alpha_DK(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_DK_best_global = float(alpha_DK_best_global)
    rd_DK_best_global = float(rd_DK_best_global)

    file_fig = generate_evidence("image", 2)
    file_table = generate_evidence("table", 2)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Load DESI BAO compressed data
    # ------------------------------------------------------------
    z_bao, bao_obs, bao_type, bao_cov, bao_cov_inv = load_desi_gaussian_bao(
        desi_bao_mean_path,
        desi_bao_cov_path,
    )

    # ------------------------------------------------------------
    # Cosmic chronometer anchor
    # ------------------------------------------------------------
    z_cc = float(Core_DESI_CC_z)
    H_cc = float(Core_DESI_CC_H)
    H_cc_err = float(Core_DESI_CC_H_err)

    # ------------------------------------------------------------
    # DK calibrated-alpha background
    # ------------------------------------------------------------
    def E_DK_alpha(z_in, alpha_DK):
        z_arr = np.asarray(z_in, dtype=float)

        E_raw = E_Relativistic(
            z_arr,
            Core_Omega_m=None,
            Omega_L_value=None,
            v_model=v_model,
            T_model=T_model,
        )

        E2_alpha = 1.0 + float(alpha_DK) * (E_raw**2 - 1.0)
        return np.sqrt(np.clip(E2_alpha, 1e-300, None))

    def H_DK_alpha(z_in, alpha_DK):
        return float(Core_H0) * E_DK_alpha(z_in, alpha_DK)

    def model_vector_for_alpha_rd(alpha_DK, rd_mpc):
        return bao_distance_vector_over_rs(
            z_bao,
            bao_type,
            lambda zz: E_DK_alpha(zz, alpha_DK),
            Core_H0,
            rd_mpc,
        )

    def chi2_total_for_alpha_rd(alpha_DK, rd_mpc):
        bao_model = model_vector_for_alpha_rd(alpha_DK, rd_mpc)

        chi2_bao = chi2_gaussian_bao(
            bao_model,
            bao_obs,
            bao_cov_inv,
        )

        H_model_cc = float(H_DK_alpha(z_cc, alpha_DK))

        chi2_cc = chi2_cosmic_chronometer_point(
            H_model_cc,
            H_obs=H_cc,
            sigma_H=H_cc_err,
        )

        return float(chi2_bao + chi2_cc), float(chi2_bao), float(chi2_cc)

    # ------------------------------------------------------------
    # Profile r_d for each alpha_DK
    # ------------------------------------------------------------
    alpha_grid = np.linspace(float(alpha_min), float(alpha_max), int(n_alpha))
    rows = []

    for alpha_val in alpha_grid:

        def objective_rd(rd_mpc):
            chi2_total, _, _ = chi2_total_for_alpha_rd(alpha_val, rd_mpc)
            return chi2_total

        res_rd = minimize_scalar(
            objective_rd,
            bounds=rd_bounds,
            method="bounded",
        )

        rd_best = float(res_rd.x)
        chi2_total, chi2_bao, chi2_cc = chi2_total_for_alpha_rd(
            alpha_val,
            rd_best,
        )

        rows.append({
            "alpha_DK": float(alpha_val),
            "rd_profiled_mpc": rd_best,
            "chi2_BAO": chi2_bao,
            "chi2_CC": chi2_cc,
            "chi2_total": chi2_total,
            "profile_success": bool(res_rd.success),
        })

    profile_df = pd.DataFrame(rows)

    # ------------------------------------------------------------
    # Locate profiled minimum
    # ------------------------------------------------------------
    idx_min = int(profile_df["chi2_total"].idxmin())

    alpha_best_profile = float(profile_df.loc[idx_min, "alpha_DK"])
    rd_best_profile = float(profile_df.loc[idx_min, "rd_profiled_mpc"])
    chi2_min = float(profile_df.loc[idx_min, "chi2_total"])
    chi2_bao_min = float(profile_df.loc[idx_min, "chi2_BAO"])
    chi2_cc_min = float(profile_df.loc[idx_min, "chi2_CC"])

    profile_df["delta_chi2"] = profile_df["chi2_total"] - chi2_min

    # Log-scale safe plotting column.
    # The true minimum remains delta_chi2 = 0 in the table.
    delta_floor = 1.0e-2
    profile_df["delta_chi2_plot"] = np.clip(
        profile_df["delta_chi2"].to_numpy(dtype=float),
        delta_floor,
        None,
    )

    # 1D approximate confidence thresholds
    delta_1sigma = 1.0
    delta_2sigma = 4.0
    delta_3sigma = 9.0

    def interval_from_delta(delta_value):
        mask = profile_df["delta_chi2"].values <= float(delta_value)
        if not np.any(mask):
            return np.nan, np.nan

        alpha_vals = profile_df["alpha_DK"].values[mask]
        return float(np.min(alpha_vals)), float(np.max(alpha_vals))

    alpha_1sig_low, alpha_1sig_high = interval_from_delta(delta_1sigma)
    alpha_2sig_low, alpha_2sig_high = interval_from_delta(delta_2sigma)
    alpha_3sig_low, alpha_3sig_high = interval_from_delta(delta_3sigma)

    # ------------------------------------------------------------
    # Save table
    # ------------------------------------------------------------
    dkrd2_to_csv(
        profile_df,
        file_table,
        table_kind="FIG02_ALPHA_DK_CHI2_PROFILE",
        figure_id=2,
        strict=False,
        index=False,
        meta={
            "alpha_DK_best_profile": alpha_best_profile,
            "rd_best_profile_mpc": rd_best_profile,
            "chi2_min": chi2_min,
            "alpha_DK_global": alpha_DK_best_global,
            "rd_DK_global_mpc": rd_DK_best_global,
            "alpha_range": f"{alpha_min} to {alpha_max}",
            "rd_bounds_mpc": rd_bounds,
            "profiled_parameter": "rd_mpc",
            "alpha_source": "Profiled DESI BAO + CC chi2 curve",
            "cc_anchor_note": "Single low-redshift cosmic chronometer anchor.",
            "dk_background": "Omega_L_value=None; no explicit Lambda term.",
        },
    )

    # ------------------------------------------------------------
    # Save stats
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "DK-RD2 alpha_DK profiled constraint",
            "alpha_DK_best_profile": alpha_best_profile,
            "alpha_DK_global_cached": alpha_DK_best_global,
            "rd_best_profile_mpc": rd_best_profile,
            "rd_global_cached_mpc": rd_DK_best_global,
            "chi2_min_total": chi2_min,
            "chi2_min_BAO": chi2_bao_min,
            "chi2_min_CC": chi2_cc_min,
            "alpha_1sigma_low_delta_chi2_1": alpha_1sig_low,
            "alpha_1sigma_high_delta_chi2_1": alpha_1sig_high,
            "alpha_2sigma_low_delta_chi2_4": alpha_2sig_low,
            "alpha_2sigma_high_delta_chi2_4": alpha_2sig_high,
            "alpha_3sigma_low_delta_chi2_9": alpha_3sig_low,
            "alpha_3sigma_high_delta_chi2_9": alpha_3sig_high,
            "n_alpha_grid": int(n_alpha),
            "rd_bounds_low_mpc": float(rd_bounds[0]),
            "rd_bounds_high_mpc": float(rd_bounds[1]),
            "cc_anchor_z": z_cc,
            "cc_anchor_H": H_cc,
            "cc_anchor_sigma_H": H_cc_err,
            "fit_mode": "profiled_rd_for_each_alpha",
        }
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=2,
        fit_mode="ALPHA_DK_CHI2_PROFILE_RD_PROFILED",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, (ax_chi, ax_rd) = plt.subplots(
        2,
        1,
        figsize=(10.8, 8.4),
        dpi=130,
        sharex=True,
        gridspec_kw={"height_ratios": [2.55, 1.25], "hspace": 0.10},
    )

    # ------------------------------------------------------------
    # Panel A — chi2 profile
    # ------------------------------------------------------------

    # Confidence band shading: draw widest first.
    if np.isfinite(alpha_3sig_low) and np.isfinite(alpha_3sig_high):
        ax_chi.axvspan(
            alpha_3sig_low,
            alpha_3sig_high,
            color=DK_RD2_color,
            alpha=0.045,
            label=None,
        )

    if np.isfinite(alpha_2sig_low) and np.isfinite(alpha_2sig_high):
        ax_chi.axvspan(
            alpha_2sig_low,
            alpha_2sig_high,
            color=DK_RD2_color,
            alpha=0.075,
            label=None,
        )

    if np.isfinite(alpha_1sig_low) and np.isfinite(alpha_1sig_high):
        ax_chi.axvspan(
            alpha_1sig_low,
            alpha_1sig_high,
            color=DK_RD2_color,
            alpha=0.12,
            label=None,
        )

    ax_chi.plot(
        profile_df["alpha_DK"],
        profile_df["delta_chi2_plot"],
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.4,
        label=r"Profiled $\Delta\chi^2(\alpha_{\rm DK})$",
    )

    ax_chi.axhline(
        delta_1sigma,
        color="black",
        linestyle=":",
        linewidth=1.1,
        alpha=0.85,
        label=r"$1\sigma$ equivalent: $\Delta\chi^2=1$",
    )

    ax_chi.axhline(
        delta_2sigma,
        color="0.35",
        linestyle="--",
        linewidth=1.0,
        alpha=0.75,
        label=r"$2\sigma$ equivalent: $\Delta\chi^2=4$",
    )

    ax_chi.axhline(
        delta_3sigma,
        color="0.55",
        linestyle="-.",
        linewidth=1.0,
        alpha=0.70,
        label=r"$3\sigma$ equivalent: $\Delta\chi^2=9$",
    )

    ax_chi.axvline(
        alpha_best_profile,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=1.6,
        alpha=0.95,
    )

    ax_chi.scatter(
        [alpha_best_profile],
        [delta_floor],
        color=DK_RD2_color,
        s=45,
        zorder=5,
        label=None,
    )

    ax_chi.set_yscale("log")
    ax_chi.set_ylim(delta_floor * 0.8, profile_df["delta_chi2_plot"].max() * 1.35)

    ax_chi.set_ylabel(r"$\Delta\chi^2(\alpha_{\rm DK})$")
    ax_chi.set_title(
        rf"Figure 02 — Statistical Constraint on the DK Projection Factor"
        "\n"
        rf"Best profile: $\alpha_{{\rm DK}}={alpha_best_profile:.4f}$, "
        rf"$r_d={rd_best_profile:.2f}\,\mathrm{{Mpc}}$, "
        rf"$\chi^2_{{min}}={chi2_min:.2f}$",
        fontsize=14,
    )

    ax_chi.grid(alpha=0.30, which="both")
    ax_chi.legend(fontsize=8, loc="upper right")

    # Result box
    result_text = (
        r"Profiled constraint:" "\n"
        rf"$\alpha_{{\rm DK}}={alpha_best_profile:.4f}$" "\n"
        rf"$r_d={rd_best_profile:.2f}\,\mathrm{{Mpc}}$" "\n"
        rf"$\chi^2_{{\rm BAO}}={chi2_bao_min:.2f}$" "\n"
        rf"$\chi^2_{{\rm CC}}={chi2_cc_min:.4f}$" "\n"
        rf"$\chi^2_{{\rm total}}={chi2_min:.2f}$"
    )

    ax_chi.text(
        0.055,
        0.73,
        result_text,
        transform=ax_chi.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.20",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.90,
        ),
    )

    method_text = (
        r"For each $\alpha_{\rm DK}$, $r_d$ is profiled by minimizing"
        "\n"
        r"$\chi^2_{\rm BAO+CC}(\alpha_{\rm DK},r_d)$."
        "\n"
        r"The minimum is not imposed; it emerges from the likelihood."
        "\n"
        r"The CC term is a single low-$z$ normalization anchor."
    )

    ax_chi.text(
        0.50,
        0.055,
        method_text,
        transform=ax_chi.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.4,
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.18",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.75,
            alpha=0.78,
        ),
    )

    # ------------------------------------------------------------
    # Panel B — profiled rd
    # ------------------------------------------------------------
    ax_rd.plot(
        profile_df["alpha_DK"],
        profile_df["rd_profiled_mpc"],
        color="black",
        linestyle="-",
        linewidth=1.8,
        label=r"Profiled $r_d(\alpha_{\rm DK})$",
    )

    ax_rd.axvline(
        alpha_best_profile,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=1.4,
        alpha=0.90,
    )

    ax_rd.axhline(
        rd_best_profile,
        color="0.35",
        linestyle="--",
        linewidth=1.2,
        alpha=0.80,
    )

    # Optional reference bounds from core constants, not hardcoded.
    ax_rd.set_ylim(
        float(rd_bounds[0]) - 2.0,
        float(rd_bounds[1]) + 2.0,
    )

    ax_rd.set_xlabel(r"Projection factor $\alpha_{\rm DK}$", labelpad=8)
    ax_rd.set_ylabel(r"Profiled $r_d$ [Mpc]")
    ax_rd.grid(alpha=0.30)
    ax_rd.legend(fontsize=8, loc="best")

    # ------------------------------------------------------------
    # Footer
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.030,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=5.4,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.88,
        bottom=0.16,
        hspace=0.14,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure03(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
):
    """
    # ============================================================
    # Figure 03 — DK-RD2 Coupling Structure
    # ============================================================
    #
    # Purpose
    # -------
    # Visualize the distinction between:
    #   1. The maximal thermo-relativistic gravitational coupling.
    #   2. The effective projected coupling observed cosmologically.
    #
    # alpha_DK is NOT hardcoded. It is obtained from get_alpha_DK(),
    # which reuses the global DESI BAO + CC calibration or computes it
    # from the same likelihood if it has not been cached.
    #
    # Definitions
    # -----------
    # Maximal coupling:
    #     Gab_max/G0 = 1 + (v²/c²)(T0/T)
    #
    # Effective projected coupling:
    #     Gab_eff/G0 = 1 + alpha_DK (v²/c²)(T0/T)
    #
    # Normalization
    # -------------
    # The curves are shown in the relativistic upper-bound regime:
    #     beta = v/c = 1
    #
    # Therefore, as T -> T0:
    #     Gab_max -> 2G0
    # while the observable excess coupling is reduced by alpha_DK.
    # ============================================================
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure03 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Get calibrated alpha_DK from global/cache system.
    # No hardcoded alpha_DK is used here.
    # ------------------------------------------------------------
    alpha_DK, rd_DK = get_alpha_DK(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_DK = float(alpha_DK)
    rd_DK = float(rd_DK)

    file_fig = generate_evidence("image", 3)
    file_table = generate_evidence("table", 3)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Temperature domain.
    # Core_TCMB_K is defined in DK_RD2_Core.py.
    # The upper value 1e8 K is only a plotting range for hot
    # astrophysical environments; it is not fitted.
    # ------------------------------------------------------------
    T_vals = np.logspace(
        np.log10(Core_TCMB_K),
        np.log10(Core_T_plot_max),
        700,
    )

    # ------------------------------------------------------------
    # Relativistic normalization.
    # beta = v/c = 1 shows the maximal coupling envelope.
    # Core_c_light is defined in DK_RD2_Core.py.
    # ------------------------------------------------------------
    beta = 1.0
    v_rel = beta * Core_c_light

    # ------------------------------------------------------------
    # Coupling structure.
    # Gab(), Core_G0 and Core_c_light are defined in DK_RD2_Core.py.
    # ------------------------------------------------------------
    Gab_max = Gab(T_vals, v_rel)
    Gab_max_over_G0 = Gab_max / Core_G0

    thermo_rel_excess = Gab_max_over_G0 - 1.0

    Gab_eff_over_G0 = 1.0 + alpha_DK * thermo_rel_excess
    Gab_eff = Core_G0 * Gab_eff_over_G0

    projected_excess_over_G0 = Gab_eff_over_G0 - 1.0
    unprojected_excess_over_G0 = Gab_max_over_G0 - Gab_eff_over_G0

    projection_ratio_total = Gab_eff_over_G0 / Gab_max_over_G0

    projection_ratio_excess = np.where(
        thermo_rel_excess > 0.0,
        projected_excess_over_G0 / thermo_rel_excess,
        np.nan,
    )

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    df = pd.DataFrame({
        "Temperature_K": T_vals,
        "beta_v_over_c": np.full_like(T_vals, beta),
        "alpha_DK": np.full_like(T_vals, alpha_DK),
        "rd_DK_Mpc": np.full_like(T_vals, rd_DK),

        "Gab_max_over_G0": Gab_max_over_G0,
        "Gab_eff_over_G0": Gab_eff_over_G0,

        "thermo_rel_excess_over_G0": thermo_rel_excess,
        "projected_excess_over_G0": projected_excess_over_G0,
        "unprojected_excess_over_G0": unprojected_excess_over_G0,

        "projection_ratio_total": projection_ratio_total,
        "projection_ratio_excess": projection_ratio_excess,

        "Gab_max_SI": Gab_max,
        "Gab_eff_SI": Gab_eff,
    })

    dkrd2_to_csv(
        df,
        file_table,
        table_kind="FIG03_COUPLING_STRUCTURE",
        figure_id=3,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_DK,
            "rd_DK_Mpc": rd_DK,
            "velocity_reference": "beta = v/c = 1",
            "temperature_range": f"{Core_TCMB_K} K to 1e8 K",
            "T0_source": "Core_TCMB_K from DK_RD2_Core.py",
            "G0_source": "Core_G0 from DK_RD2_Core.py",
            "c_source": "Core_c_light from DK_RD2_Core.py",
            "alpha_source": "DESI BAO + CC calibration",
            "interpretation": (
                "Projected excess coupling, not projected total coupling. "
                "Gab_max approaches 2G0 in the beta=1, T -> T0 limit."
            ),
        },
    )

    # ------------------------------------------------------------
    # Plot.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10.5, 6.8), dpi=130)

    ax.plot(
        T_vals,
        Gab_max_over_G0,
        color="red",
        linestyle=LCDM_linestyle,
        linewidth=2.3,
        label=r"$G_{ab}^{max}/G_0$  $(\beta=1,\ \alpha_{\rm DK}=1)$",
    )

    ax.plot(
        T_vals,
        Gab_eff_over_G0,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.8,
        label=rf"$G_{{ab}}^{{eff}}/G_0$  $(\alpha_{{\rm DK}}={alpha_DK:.4f})$",
    )

    ax.fill_between(
        T_vals,
        Gab_eff_over_G0,
        Gab_max_over_G0,
        alpha=0.25,
        label=r"Unprojected excess coupling",
    )

    ax.axhline(
        1.0,
        color="black",
        linewidth=1.1,
        alpha=0.85,
        label=r"Newtonian limit $G_0$",
    )

    ax.axhline(
        2.0,
        color="red",
        linestyle=":",
        linewidth=1.0,
        alpha=0.65,
        label=r"Relativistic upper limit $G_{ab}^{max}\to 2G_0$",
    )

    ax.axvline(
        Core_TCMB_K,
        color="black",
        linestyle=":",
        linewidth=1.3,
        alpha=0.9,
        label=rf"$T_0={Core_TCMB_K:.4f}\,\mathrm{{K}}$",
    )

    ax.set_ylim(0.995, 2.05)

    ax.set_xscale("log")
    ax.set_xlabel("Temperature T [K]")
    ax.set_ylabel(r"Normalized gravitational coupling $G_{ab}/G_0$")
    ax.set_title(
        rf"Figure 03 — DK-RD2 Coupling Structure"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_DK:.4f}$ from DESI BAO + CC calibration",
        fontsize=14,
    )

    # ------------------------------------------------------------
    # Secondary axis — projected excess fraction.
    # This ratio is the projected fraction of the EXCESS coupling:
    # (Gab_eff/G0 - 1) / (Gab_max/G0 - 1) = alpha_DK
    # ------------------------------------------------------------
    ax_frac = ax.twinx()

    ax_frac.plot(
        T_vals,
        projection_ratio_excess,
        color="black",
        linestyle=":",
        alpha=0.70,
        linewidth = 1.5,
        label="Projected excess fraction",
    )

    ax_frac.set_ylim(0.0, 1.05)
    ax_frac.set_ylabel(r"Projected excess fraction")
    ax_frac.tick_params(axis="y")

    ax.grid(True, which="both", linestyle="--", alpha=0.30)

    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax_frac.get_legend_handles_labels()

    ax.legend(
        lines_1 + lines_2,
        labels_1 + labels_2,
        fontsize=7.7,
        loc="upper right",
    )

    # ------------------------------------------------------------
    # Equation box.
    # ------------------------------------------------------------
    eq_text = (
        r"$G_{ab}^{max}=G_0\left[1+\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        r"$G_{ab}^{eff}=G_0\left[1+\alpha_{\rm DK}\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        rf"$\beta=v/c=1,\quad \alpha_{{\rm DK}}={alpha_DK:.4f}$"
    )

    ax.text(
        0.08,
        0.96,
        eq_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        bbox=dict(
            boxstyle="round,pad=0.42",
            facecolor="white",
            edgecolor="black",
            linewidth=1.0,
            alpha=0.90,
        ),
    )

    # ------------------------------------------------------------
    # Interpretation box.
    # ------------------------------------------------------------
    insight_text = (
        "Projection interpretation:\n"
        rf"Only {alpha_DK * 100:.2f}% of the coupling excess\n"
        "is projected into observable gravity.\n"
        r"In the $\beta=1,\ T\to T_0$ limit:" "\n"
        r"$G_{ab}^{max}\to 2G_0$."
    )

    ax.text(
        0.31,
        0.44,
        insight_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        linespacing=1.18,
        bbox=dict(
            boxstyle="round,pad=0.38",
            facecolor="white",
            edgecolor="0.45",
            linewidth=1.0,
            alpha=0.88,
        ),
    )

    # ------------------------------------------------------------
    # Footer.
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.025,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=5.8,
        color=DK_RD2_color,
    )

    plt.tight_layout(rect=(0, 0.06, 1, 0.92))
    plt.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    # ------------------------------------------------------------
    # Stats CSV.
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "DK-RD2 coupling structure",
            "alpha_DK": float(alpha_DK),
            "rd_mpc": float(rd_DK),
            "T_min_K": float(T_vals.min()),
            "T_max_K": float(T_vals.max()),
            "beta_v_over_c": float(beta),
            "Gab_max_over_G0_at_T0_beta1": float(Gab_max_over_G0[0]),
            "Gab_eff_over_G0_at_T0_beta1": float(Gab_eff_over_G0[0]),
            "projected_excess_fraction": float(alpha_DK),
            "fit_source": "DESI BAO + CC calibration",
            "fit_mode": "COUPLING_STRUCTURE_PROJECTED_EXCESS",
        }
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=3,
        fit_mode="COUPLING_STRUCTURE_PROJECTED_EXCESS",
        index=False,
    )

    return file_fig, file_table, file_stats

def generate_figure04(
        desi_bao_mean_path: str | None = None,
        desi_bao_cov_path: str | None = None,
        *,
        z_min: float = 0.001,
        z_max: float = 2.30,
        n_z: int = 900,
        Core_H0: float | None = None,
        Core_Omega_m_LCDM: float | None = None,
        Omega_L_LCDM: float | None = None,
        Omega_L_DK: float | None = None,
        v_model=None,
        T_model=None,
):
    """
    # ============================================================
    # Figure 04 — Explicit Diagnostic Reconstruction of w(z)
    # ============================================================
    #
    # Purpose
    # -------
    # Provide a direct reconstruction of the effective dark-energy-like
    # equation of state w(z) inferred when the DK-RD2 calibrated expansion
    # history is interpreted through a ΛCDM matter reference.
    #
    # Important
    # ---------
    # This is a diagnostic reconstruction. It does not represent a
    # fundamental dark energy component in DK-RD2.
    #
    # Definitions
    # -----------
    # Ω_DE_eff(z) = E²(z) - Ω_m,ref (1+z)³
    #
    # w(z) = -1 + (1+z)/3 * d ln Ω_DE_eff(z) / dz
    #
    # with:
    # Ω_m,ref = Ω_m,ΛCDM = 0.315
    #
    # α_DK is not hardcoded. It is obtained from get_alpha_DK().
    # ============================================================
    """

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure04 requires desi_bao_mean_path and desi_bao_cov_path.")

    if Core_H0 is None:
        Core_H0 = float(Core_Hubble_H0)

    if Core_Omega_m_LCDM is None:
        Core_Omega_m_LCDM = float(Core_OMEGA_M_LCDM)

    if Omega_L_LCDM is None:
        Omega_L_LCDM = float(Core_OMEGA_L_LCDM)

    # ------------------------------------------------------------
    # Calibrated alpha_DK
    # ------------------------------------------------------------
    alpha_DK, rd_DK = get_alpha_DK(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_DK = float(alpha_DK)
    rd_DK = float(rd_DK)

    file_fig = generate_evidence("image", 4)
    file_table = generate_evidence("table", 4)
    file_stats = file_table.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------
    # Redshift grid
    # ------------------------------------------------------------
    z = np.linspace(float(z_min), float(z_max), int(n_z))

    # ------------------------------------------------------------
    # ΛCDM reference
    # ------------------------------------------------------------
    E_lcdm = E_LCDM(
        z,
        Core_Omega_m=Core_Omega_m_LCDM,
        Core_Omega_L=Omega_L_LCDM,
    )

    Omega_DE_lcdm = Omega_DE_eff_from_E(
        z,
        E_lcdm,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    w_lcdm = -1.0 * np.ones_like(z)

    # ------------------------------------------------------------
    # DK-RD2 calibrated expansion history
    # ------------------------------------------------------------
    # Core DK-RD2 uses the maximal thermodynamic–relativistic coupling:
    #
    #     Gab_max(T, v) = G0 * [1 + (v²/c²)(T0/T)]
    #
    # This produces the raw DK-RD2 expansion envelope E_DK_raw(z).
    # The paper, however, compares the observable large-scale projection:
    #
    #     Gab_eff(T, v) = G0 * [1 + alpha_DK * (v²/c²)(T0/T)]
    #
    # Therefore, the observable expansion curve is obtained by projecting
    # the raw DK-RD2 correction through alpha_DK:
    #
    #     E_DK_eff²(z) = 1 + alpha_DK * [E_DK_raw²(z) - 1]
    #
    # This keeps the core formula and the paper formula fully consistent:
    # Gab_max is the complete coupling envelope, while Gab_eff is the
    # DESI-calibrated projected coupling observed at cosmological scales.
    # ------------------------------------------------------------
    E_dk_raw = E_Relativistic(
        z,
        Core_Omega_m=None,
        Omega_L_value=Omega_L_DK,
        v_model=v_model,
        T_model=T_model,
    )

    E2_dk_eff = 1.0 + alpha_DK * (E_dk_raw ** 2 - 1.0)
    E_dk_eff = np.sqrt(np.clip(E2_dk_eff, 1e-300, None))


    # ------------------------------------------------------------
    # Explicit diagnostic reconstruction
    # ------------------------------------------------------------
    Omega_DE_dk = Omega_DE_eff_from_E(
        z,
        E_dk_eff,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    w_dk = w_eff_from_E(
        z,
        E_dk_eff,
        Omega_m_ref=Core_Omega_m_LCDM,
    )

    delta_w = w_dk - w_lcdm

    finite_w = np.isfinite(w_dk)
    finite_omega = np.isfinite(Omega_DE_dk)

    # ------------------------------------------------------------
    # Helper
    # ------------------------------------------------------------
    def interp_safe(x_grid, y_grid, x_value):
        mask = np.isfinite(x_grid) & np.isfinite(y_grid)
        if np.count_nonzero(mask) < 2:
            return np.nan
        return float(np.interp(x_value, x_grid[mask], y_grid[mask]))

    w0 = interp_safe(z, w_dk, 0.01)
    w05 = interp_safe(z, w_dk, 0.5)
    w1 = interp_safe(z, w_dk, 1.0)
    w2 = interp_safe(z, w_dk, 2.0)

    omega0 = interp_safe(z, Omega_DE_dk, 0.01)
    omega05 = interp_safe(z, Omega_DE_dk, 0.5)
    omega1 = interp_safe(z, Omega_DE_dk, 1.0)
    omega2 = interp_safe(z, Omega_DE_dk, 2.0)

    # ------------------------------------------------------------
    # Save table
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "z": z,
        "alpha_DK": np.full_like(z, alpha_DK),
        "rd_DK_Mpc": np.full_like(z, rd_DK),
        "Omega_m_ref_LCDM": np.full_like(z, Core_Omega_m_LCDM),

        "E_LCDM": E_lcdm,
        "E_DK_raw_from_Gab_max": E_dk_raw,
        "E_DK_eff_projected_by_alpha_DK": E_dk_eff,

        "Omega_DE_eff_LCDM": Omega_DE_lcdm,
        "Omega_DE_eff_DKRD2_diagnostic": Omega_DE_dk,

        "w_LCDM": w_lcdm,
        "w_DKRD2_reconstructed": w_dk,
        "delta_w_DK_minus_LCDM": delta_w,
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG04_EXPLICIT_W_RECONSTRUCTION",
        figure_id=4,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_DK,
            "rd_DK_Mpc": rd_DK,
            "Omega_m_ref_LCDM": Core_Omega_m_LCDM,
            "z_range": f"{z_min} to {z_max}",
            "alpha_source": "DESI BAO + CC calibration",
            "coupling_interpretation": (
                "Gab_max is the full thermodynamic-relativistic coupling used by the core; "
                "Gab_eff applies alpha_DK as the observable large-scale geometric projection."
            ),
            "expansion_projection": (
                "E_DK_eff^2 = 1 + alpha_DK * (E_DK_raw^2 - 1)."
            ),
            "interpretation": (
                "Diagnostic w(z) reconstruction using ΛCDM matter reference; "
                "not a fundamental dark energy component in DK-RD2."
            ),
        },
    )

    # ------------------------------------------------------------
    # Save stats
    # ------------------------------------------------------------
    stats = pd.DataFrame([
        {
            "model": "ΛCDM reference",
            "alpha_DK": np.nan,
            "rd_mpc": np.nan,
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "w_z_0p01": -1.0,
            "w_z_0p5": -1.0,
            "w_z_1": -1.0,
            "w_z_2": -1.0,
            "Omega_DE_z_0p01": float(Omega_L_LCDM),
            "Omega_DE_z_0p5": float(Omega_L_LCDM),
            "Omega_DE_z_1": float(Omega_L_LCDM),
            "Omega_DE_z_2": float(Omega_L_LCDM),
            "w_min": -1.0,
            "w_max": -1.0,
            "mean_delta_w_vs_minus1": 0.0,
            "N_finite": int(len(z)),
            "fit_mode": "reference_constant",
        },
        {
            "model": "DK-RD2 explicit diagnostic w(z)",
            "alpha_DK": float(alpha_DK),
            "rd_mpc": float(rd_DK),
            "Omega_m_ref": float(Core_Omega_m_LCDM),
            "w_z_0p01": w0,
            "w_z_0p5": w05,
            "w_z_1": w1,
            "w_z_2": w2,
            "Omega_DE_z_0p01": omega0,
            "Omega_DE_z_0p5": omega05,
            "Omega_DE_z_1": omega1,
            "Omega_DE_z_2": omega2,
            "w_min": float(np.nanmin(w_dk)) if np.any(finite_w) else np.nan,
            "w_max": float(np.nanmax(w_dk)) if np.any(finite_w) else np.nan,
            "mean_delta_w_vs_minus1": float(np.nanmean(w_dk + 1.0)) if np.any(finite_w) else np.nan,
            "Omega_DE_min": float(np.nanmin(Omega_DE_dk)) if np.any(finite_omega) else np.nan,
            "Omega_DE_max": float(np.nanmax(Omega_DE_dk)) if np.any(finite_omega) else np.nan,
            "N_finite": int(np.count_nonzero(finite_w)),
            "fit_mode": "diagnostic_LCDM_reference",
        },
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=4,
        fit_mode="EXPLICIT_W_RECONSTRUCTION",
        index=False,
    )

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    fig, (ax_omega, ax_w, ax_delta) = plt.subplots(
        3,
        1,
        figsize=(11.2, 10.2),
        dpi=130,
        sharex=True,
        gridspec_kw={"height_ratios": [1.25, 2.4, 1.15], "hspace": 0.10},
    )

    # ------------------------------------------------------------
    # Panel A — Omega_DE_eff(z)
    # ------------------------------------------------------------
    ax_omega.plot(
        z,
        Omega_DE_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=1.8,
        label=r"$\Lambda$CDM: $\Omega_{\Lambda}=\mathrm{const.}$",
    )

    ax_omega.plot(
        z,
        Omega_DE_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"DK-RD2 diagnostic $\Omega_{\rm DE}^{eff}(z)$",
    )

    ax_omega.axhline(
        0.0,
        color="black",
        linewidth=0.9,
        alpha=0.45,
    )

    ax_omega.set_ylabel(r"$\Omega_{\rm DE}^{eff}(z)$")
    ax_omega.set_title(
        rf"Figure 04 — Explicit Diagnostic Reconstruction of $w(z)$"
        "\n"
        rf"$\alpha_{{\rm DK}}={alpha_DK:.4f}$, "
        rf"$r_d={rd_DK:.2f}\,\mathrm{{Mpc}}$ "
        rf"(DESI BAO + CC calibrated)",
        fontsize=14,
    )
    ax_omega.grid(alpha=0.30)
    ax_omega.legend(fontsize=8, loc="best")

    # ------------------------------------------------------------
    # Panel B — Explicit reconstructed w(z)
    # ------------------------------------------------------------
    ax_w.plot(
        z,
        w_lcdm,
        color=LCDM_color,
        linestyle=LCDM_linestyle,
        linewidth=2.0,
        label=r"$\Lambda$CDM: $w=-1$",
    )

    ax_w.plot(
        z,
        w_dk,
        color=DK_RD2_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.5,
        label=(
            rf"DK-RD2 reconstructed $w(z)$ "
            rf"$(\Omega_{{m,\rm ref}}={Core_Omega_m_LCDM:.3f})$"
        ),
    )

    ax_w.axhline(
        -1.0,
        color="black",
        linewidth=1.0,
        alpha=0.35,
    )

    ax_w.set_ylabel(r"Reconstructed $w(z)$")
    ax_w.grid(alpha=0.30)
    ax_w.legend(fontsize=8.3, loc="best")

    eq_text = (
        r"$G_{ab}^{max}=G_0\left[1+\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        r"$G_{ab}^{eff}=G_0\left[1+\alpha_{\rm DK}\frac{v^2}{c^2}\frac{T_0}{T}\right]$" "\n"
        r"$E_{\rm DK,eff}^2=1+\alpha_{\rm DK}\left(E_{\rm DK,raw}^2-1\right)$" "\n"
        r"$\Omega_{\rm DE}^{eff}=E_{\rm DK,eff}^2-\Omega_{m,\rm ref}(1+z)^3$" "\n"
        r"$w(z)=-1+\frac{1+z}{3}\frac{d\ln\Omega_{\rm DE}^{eff}}{dz}$"
    )

    ax_w.text(
        0.03,
        0.87,
        eq_text,
        transform=ax_w.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.38",
            facecolor="white",
            edgecolor="black",
            linewidth=0.9,
            alpha=0.90,
        ),
    )

    result_text = (
        r"DK-RD2 reconstructed:" "\n"
        rf"$w(0.01)={w0:.3f}$" "\n"
        rf"$w(0.5)={w05:.3f}$" "\n"
        rf"$w(1.0)={w1:.3f}$"
    )

    ax_w.text(
        0.75,
        0.4,
        result_text,
        transform=ax_w.transAxes,
        ha="left",
        va="top",
        fontsize=8.7,
        linespacing=1.25,
        bbox=dict(
            boxstyle="round,pad=0.38",
            facecolor="white",
            edgecolor="0.45",
            linewidth=0.9,
            alpha=0.88,
        ),
    )

    note_text = (
        "Diagnostic reconstruction:\n"
        rf"$\Omega_{{m,\rm ref}}={Core_Omega_m_LCDM:.3f}$ is the ΛCDM matter reference."
        "\nThis curve is not a fundamental dark-energy component in DK-RD2."
    )

    ax_w.text(
        0.03,
        0.07,
        note_text,
        transform=ax_w.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.5,
        linespacing=1.18,
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.85,
            alpha=0.88,
        ),
    )

    # ------------------------------------------------------------
    # Panel C — delta w
    # ------------------------------------------------------------
    ax_delta.plot(
        z,
        delta_w,
        color=DK_RD2_color,
        linestyle="-",
        linewidth=2.0,
        label=r"$\Delta w(z)=w_{\rm DK}(z)+1$",
    )

    ax_delta.axhline(
        0.0,
        color="black",
        linewidth=1.0,
        alpha=0.70,
    )

    ax_delta.set_xlabel("Redshift z")
    ax_delta.set_ylabel(r"$\Delta w$")
    ax_delta.grid(alpha=0.30)
    ax_delta.legend(fontsize=8, loc="best")

    highz_note = (
        "The apparent divergence at high redshift does not correspond "
        "to a physical singularity in DK-RD2.\n"
        "It arises because the diagnostic ΛCDM-like effective "
        "dark-energy density approaches zero, making the reconstructed "
        "w(z) parametrization ill-defined."
    )

    ax_w.text(
        0.03,
        0.2,
        highz_note,
        transform=ax_w.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.6,
        linespacing=1.15,
        bbox=dict(
            boxstyle="round,pad=0.30",
            facecolor="white",
            edgecolor="0.55",
            linewidth=0.8,
            alpha=0.88,
        ),
    )

    # ------------------------------------------------------------
    # Footer
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.018,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        fontsize=5.7,
        color=DK_RD2_color,
    )

    fig.subplots_adjust(
        top=0.91,
        bottom=0.12,
        hspace=0.18,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def generate_figure05(
    desi_bao_mean_path: str | None = None,
    desi_bao_cov_path: str | None = None,
):
    """
    # ============================================================
    # Figure 05 — Thermodynamic Activity vs Projected Coupling
    # ============================================================
    #
    # Purpose
    # -------
    # Illustrate the DK-RD2 distinction between:
    #   1. High energy density
    #   2. High thermodynamic activity
    #   3. Gravitational confinement
    #
    # The key message is that high energy density does not necessarily
    # imply high thermodynamic activity. A saturation regime can be highly
    # confined while having suppressed accessible degrees of freedom.
    #
    # Panels
    # ------
    # A) Conceptual environment table.
    # B) Temperature-dependent maximal and projected coupling.
    # C) Projection identity residual:
    #
    #      [(G_eff/G0 - 1)/(G_max/G0 - 1)] - alpha_DK
    #
    # Notes
    # -----
    # - alpha_DK is obtained from DESI BAO + CC calibration.
    # - Core constants are resolved from DK_RD2_Core.py.
    # - Figure/table/stats IDs are consistently set to 5.
    # ============================================================
    """

    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.patheffects as pe

    if desi_bao_mean_path is None or desi_bao_cov_path is None:
        raise ValueError("Figure05 requires desi_bao_mean_path and desi_bao_cov_path.")

    # ------------------------------------------------------------
    # Calibrated DK projection factor.
    # No hardcoded alpha_DK is used.
    # ------------------------------------------------------------
    alpha_DK, rd_DK = get_alpha_DK(
        desi_bao_mean_path=desi_bao_mean_path,
        desi_bao_cov_path=desi_bao_cov_path,
    )

    alpha_DK = float(alpha_DK)
    rd_DK = float(rd_DK)

    file_fig = generate_evidence("image", 5)
    table_file = generate_evidence("table", 5)
    file_stats = table_file.replace(".csv", "_stats.csv")

    # ------------------------------------------------------------------
    # 1. Representative environments
    # ------------------------------------------------------------------
    # These temperatures are representative physical regimes, not fitted
    # parameters. Core_TCMB_K and Core_c_light come from DK_RD2_Core.py.
    T_flare_max = 6.0e7
    T_core = 1.5e7
    T_corona = 1.0e6
    T_surface_avg = 5778.0
    T_cmb = Core_TCMB_K

    environments = [
        {
            "Environment": "Stellar core",
            "Temperature_K": T_core,
            "Energy_density": "High",
            "Thermodynamic_activity": "High",
            "Confinement": "Partial",
            "Accessible_DOF": "Many",
            "Interpretation": "Dense but dynamically active",
            "Representative_beta": 1.3e-3,
        },
        {
            "Environment": "Photosphere",
            "Temperature_K": T_surface_avg,
            "Energy_density": "Moderate",
            "Thermodynamic_activity": "High",
            "Confinement": "Weak",
            "Accessible_DOF": "Many",
            "Interpretation": "Radiative redistribution",
            "Representative_beta": 1.0,
        },
        {
            "Environment": "Corona",
            "Temperature_K": T_corona,
            "Energy_density": "Low",
            "Thermodynamic_activity": "High",
            "Confinement": "Weak",
            "Accessible_DOF": "Many",
            "Interpretation": "Low density, high particle activity",
            "Representative_beta": 1.3e-3,
        },
        {
            "Environment": "Interstellar / CMB",
            "Temperature_K": T_cmb,
            "Energy_density": "Very low",
            "Thermodynamic_activity": "Low",
            "Confinement": "Weak",
            "Accessible_DOF": "Available",
            "Interpretation": "Cosmological thermal floor",
            "Representative_beta": 1.0,
        },
        {
            "Environment": "Saturation regime",
            "Temperature_K": T_cmb,
            "Energy_density": "Extreme",
            "Thermodynamic_activity": "Suppressed",
            "Confinement": "Maximal",
            "Accessible_DOF": "Minimized",
            "Interpretation": "Thermodynamic freezing",
            "Representative_beta": 0.0,
        },
    ]

    rows = []

    for env in environments:
        beta = float(env["Representative_beta"])
        v = beta * Core_c_light

        Gab_max_ratio = Gab(env["Temperature_K"], v) / Core_G0
        dynamic_term = Gab_max_ratio - 1.0
        Gab_eff_ratio = 1.0 + alpha_DK * dynamic_term

        if abs(dynamic_term) > 1e-300:
            projection_fraction = (Gab_eff_ratio - 1.0) / dynamic_term
            projection_residual = projection_fraction - alpha_DK
        else:
            projection_fraction = np.nan
            projection_residual = np.nan

        rows.append({
            "Environment": env["Environment"],
            "Temperature_K": env["Temperature_K"],
            "Energy_density": env["Energy_density"],
            "Thermodynamic_activity": env["Thermodynamic_activity"],
            "Confinement": env["Confinement"],
            "Accessible_DOF": env["Accessible_DOF"],
            "Representative_beta": beta,
            "alpha_DK": alpha_DK,
            "Gab_max_over_G0": float(Gab_max_ratio),
            "Gab_eff_over_G0": float(Gab_eff_ratio),
            "Unprojected_component_over_G0": float(Gab_max_ratio - Gab_eff_ratio),
            "Projection_fraction_dynamic": (
                float(projection_fraction) if np.isfinite(projection_fraction) else np.nan
            ),
            "Projection_residual": (
                float(projection_residual) if np.isfinite(projection_residual) else np.nan
            ),
            "Interpretation": env["Interpretation"],
        })

    df = pd.DataFrame(rows)

    # ------------------------------------------------------------
    # Save evidence table.
    # ------------------------------------------------------------
    dkrd2_to_csv(
        df,
        table_file,
        table_kind="FIG05_STELLAR_ACTIVITY_PROJECTED_COUPLING",
        figure_id=5,
        strict=False,
        index=False,
        meta={
            "alpha_DK": alpha_DK,
            "rd_DK_Mpc": rd_DK,
            "alpha_source": "DESI BAO + CC calibration",
            "interpretation": (
                "Representative environments illustrating thermodynamic activity "
                "versus gravitational confinement."
            ),
        },
    )

    # ------------------------------------------------------------------
    # 2. Temperature grid
    # ------------------------------------------------------------------
    T_min_plot = T_cmb
    T_max_plot = T_flare_max

    T_vals = np.logspace(
        np.log10(T_min_plot),
        np.log10(T_max_plot),
        400,
    )

    logT_vals = np.log10(T_vals)
    logT_surf = np.log10(T_surface_avg)
    logT_cmb = np.log10(T_cmb)

    cold_frac = np.zeros_like(T_vals)

    mask_hot = T_vals >= T_surface_avg
    cold_frac[mask_hot] = 0.0

    mask_cold = T_vals < T_surface_avg
    frac_linear = (logT_surf - logT_vals[mask_cold]) / (logT_surf - logT_cmb)
    frac_linear = np.clip(frac_linear, 0.0, 1.0)

    # Softer exponent improves visibility of the transition region.
    gamma = 45
    cold_frac[mask_cold] = frac_linear ** gamma
    cold_frac = np.clip(cold_frac, 0.0, 1.0)

    bg_data = cold_frac[:, np.newaxis]

    # ------------------------------------------------------------------
    # 3. Figure layout
    # ------------------------------------------------------------------
    fig = plt.figure(
        "Figure 05: Thermodynamic activity vs projected gravitational coupling",
        figsize=(12.5, 13.2),
    )

    gs = gridspec.GridSpec(
        nrows=3,
        ncols=1,
        height_ratios=[1.7, 2.55, 1.05],
        hspace=0.52,
    )

    # ------------------------------------------------------------------
    # Panel A — conceptual table
    # ------------------------------------------------------------------
    ax_table = fig.add_subplot(gs[0])
    ax_table.axis("off")

    col_labels = [
        "Environment",
        "T [K]",
        "Energy\ndensity",
        "Thermodynamic\nactivity",
        "Confinement",
        "Accessible\nDOF",
        r"$G_{ab}^{max}/G_0$",
        r"$G_{ab}^{eff}/G_0$",
        "Physical\ninterpretation",
    ]

    cell_rows = []
    for _, r in df.iterrows():
        cell_rows.append([
            r["Environment"],
            f"{r['Temperature_K']:.3e}",
            r["Energy_density"],
            r["Thermodynamic_activity"],
            r["Confinement"],
            r["Accessible_DOF"],
            f"{r['Gab_max_over_G0']:.7f}",
            f"{r['Gab_eff_over_G0']:.7f}",
            r["Interpretation"],
        ])

    table = ax_table.table(
        cellText=cell_rows,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )

    # Smaller font and taller cells prevent text overflow.
    table.auto_set_font_size(False)
    table.set_fontsize(5.9)
    table.scale(1.02, 1.62)

    # Wider last column for the physical interpretation text.
    col_widths = [0.135, 0.085, 0.095, 0.125, 0.095, 0.095, 0.10, 0.10, 0.17]
    n_rows = len(cell_rows)

    for col, width in enumerate(col_widths):
        for row in range(n_rows + 1):
            table[(row, col)].set_width(width)

    env_cmap = LinearSegmentedColormap.from_list(
        "dk_activity_map",
        [
            (0.00, "#b00020"),
            (0.45, "#ff7f0e"),
            (0.70, "#ffd21f"),
            (1.00, "#050505"),
        ],
    )

    row_colors = {
        "Stellar core": "#b00020",
        "Photosphere": "#ff7f0e",
        "Corona": "#ffd21f",
        "Interstellar / CMB": "#050505",
        "Saturation regime": "#6b6b6b",
    }

    for row_idx, row_data in enumerate(cell_rows, start=1):
        env_name = row_data[0]
        bg = row_colors.get(env_name, "white")

        for col_idx in range(len(col_labels)):
            cell = table[(row_idx, col_idx)]
            cell.set_facecolor(bg)

            if env_name == "Corona":
                cell.get_text().set_color("black")
            else:
                cell.get_text().set_color("white")

    for col_idx in range(len(col_labels)):
        header = table[(0, col_idx)]
        header.set_facecolor("black")
        header.get_text().set_color("white")

    ax_table.set_title(
        rf"Figure 05 — Thermodynamic Activity vs Projected Gravitational Coupling "
        rf"($\alpha_{{\rm DK}}={alpha_DK:.4f}$)",
        fontsize=12.4,
        pad=9,
    )

    formula_text = (
        r"$G_{ab}^{eff}(T,v)=G_0\left[1+\alpha_{\rm DK}"
        r"\left(\frac{v}{c}\right)^2\left(\frac{T_0}{T}\right)\right]$"
    )

    ax_table.text(
        0.5,
        0.855,
        formula_text,
        ha="center",
        va="bottom",
        fontsize=12.0,
        transform=ax_table.transAxes,
    )

    ax_table.text(
        0.5,
        0.02,
        "Conceptual comparison supporting the distinction between thermodynamic activity and gravitational confinement.\n"
        "The effective coupling uses the DESI-calibrated projection factor; high energy density does not necessarily imply high thermodynamic activity.",
        ha="center",
        va="top",
        fontsize=7.6,
        transform=ax_table.transAxes,
    )

    # ------------------------------------------------------------------
    # Panel B — thermal background + coupling curves
    # ------------------------------------------------------------------
    ax_bg = fig.add_subplot(gs[1])

    x_min, x_max = 0.8, 2.2

    im = ax_bg.imshow(
        bg_data,
        extent=(x_min, x_max, T_min_plot, T_max_plot),
        aspect="auto",
        origin="lower",
        cmap=env_cmap,
        vmin=0.0,
        vmax=1.0,
    )

    ax_bg.set_yscale("log")
    ax_bg.set_xlim(x_min, x_max)
    ax_bg.set_xticks([])
    ax_bg.set_ylim(T_min_plot, T_max_plot)
    ax_bg.invert_yaxis()

    ax = ax_bg.twiny()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(T_min_plot, T_max_plot)
    ax.invert_yaxis()

    trans_x = ax.get_xaxis_transform()

    Gab_curve_photons_max = Gab(T_vals, Core_c_light) / Core_G_const
    Gab_curve_photons_eff = 1.0 + alpha_DK * (Gab_curve_photons_max - 1.0)

    idx_eff = np.argmin(np.abs(Gab_curve_photons_eff - (1.0 + alpha_DK)))
    Gab_real_eff = Gab_curve_photons_eff[idx_eff] * Core_G_const

    G0_mantissa, G0_exponent = sci_notation(Core_G_const)
    gab_eff_mantissa, gab_eff_exponent = sci_notation(Gab_real_eff)

    common_exp = min(G0_exponent, gab_eff_exponent)
    G0_mant_common = G0_mantissa * 10 ** (G0_exponent - common_exp)
    Gab_eff_mant_common = gab_eff_mantissa * 10 ** (gab_eff_exponent - common_exp)

    ax.plot(
        [1.0, 1.0],
        [T_min_plot, T_max_plot],
        color="black",
        linewidth=1.8,
        label=rf"$G_0 = {G0_mant_common:.2f}\times 10^{{{common_exp}}}"
              rf"\,\mathrm{{m^3\,kg^{{-1}}\,s^{{-2}}}}$",
    )

    ax.plot(
        [1.0 + alpha_DK, 1.0 + alpha_DK],
        [T_min_plot, T_max_plot],
        color=DK_RD2_color,
        linewidth=2.8,
        alpha=0.90,
        label=(
            rf"$G_{{ab}}^{{eff}}(T_0,\beta=1) "
            rf"\approx {Gab_eff_mant_common:.2f}\times 10^{{{common_exp}}}"
            rf"\,\mathrm{{m^3\,kg^{{-1}}\,s^{{-2}}}}$"
        ),
    )

    ax.plot(
        [2.0, 2.0],
        [T_min_plot, T_max_plot],
        color="red",
        linewidth=1.8,
        alpha=0.70,
        linestyle="--",
        label=r"$G_{ab}^{max}(T_0,\beta=1)=2G_0$",
    )

    species_plot = [
        {"name": "Photons", "beta": 1.0, "linestyle": "-"},
        {"name": "Neutrinos", "beta": 0.999999, "linestyle": "--"},
        {"name": "Protons", "beta": 1.3e-3, "linestyle": ":"},
    ]

    species_colors = {
        "Photons": "#00c8ff",
        "Neutrinos": "#ff9f1a",
        "Protons": "#32cd32",
    }

    projection_residuals = {}

    for sp in species_plot:
        beta = float(sp["beta"])
        v = beta * Core_c_light

        Gab_curve_max = Gab(T_vals, v) / Core_G0
        Gab_curve_eff = 1.0 + alpha_DK * (Gab_curve_max - 1.0)

        dynamic_term = Gab_curve_max - 1.0
        projection_fraction = np.where(
            np.abs(dynamic_term) > 1e-300,
            (Gab_curve_eff - 1.0) / dynamic_term,
            np.nan,
        )

        projection_residual = projection_fraction - alpha_DK
        projection_residuals[sp["name"]] = projection_residual

        ax.plot(
            Gab_curve_max,
            T_vals,
            linestyle=sp["linestyle"],
            linewidth=1.1,
            alpha=0.38,
            color=species_colors[sp["name"]],
            label=f"{sp['name']} max (β={beta:.3g})",
        )

        ax.plot(
            Gab_curve_eff,
            T_vals,
            linestyle=sp["linestyle"],
            linewidth=1.7,
            alpha=0.95,
            color=species_colors[sp["name"]],
            label=f"{sp['name']} eff",
        )

    ax.set_xlabel(r"$G_{ab}(T, v) / G_0$")
    ax.set_title(
        r"Panel B — Temperature-dependent maximal and projected gravitational coupling",
        fontsize=12.5,
        pad=10,
    )

    ax.legend(
        fontsize=6.2,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.34),
        framealpha=0.92,
    )

    ax.grid(True, axis="x", linestyle="--", linewidth=0.5, alpha=0.4)

    ax_bg.set_ylabel(
        f"Temperature T [K]\n"
        f"(bottom: stellar core ≈1.5×10⁷ K, top: CMB ≈{Core_TCMB_K:.4f} K)",
        fontsize=8.2,
    )

    cbar = fig.colorbar(im, ax=ax_bg, pad=0.02)

    cbar.ax.yaxis.set_ticks_position("right")
    cbar.ax.yaxis.set_label_position("left")

    cbar.set_label(
        "Thermodynamic regime: active (red) → constrained (black)",
        color="white",
        fontsize=7.3,
        labelpad=7,
    )

    cbar.set_ticks([0.0, 0.4, 0.7, 1.0])
    cbar.set_ticklabels([
        "Core (~10⁷ K)\nActive",
        "Corona (~10⁶ K)\nDynamic",
        "Surface (~10³ K)\nStable",
        "CMB (~10⁰ K)\nLow-energy",
    ])

    cbar.ax.tick_params(labelsize=6.5)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.get_yticklabels(), color="white")
    cbar.ax.yaxis.label.set_color("black")

    # Bottom labels. The right one is shifted left/down to avoid overlap.
    ax.text(
        1.0,
        -0.008,
        rf"$G_0 = {G0_mant_common:.2f}\times10^{{{common_exp}}}$",
        transform=trans_x,
        ha="center",
        va="top",
        fontsize=8.6,
    )

    ax.text(
        1.0 + alpha_DK,
        -0.008,
        rf"$G_{{ab}}^{{eff}} = {Gab_eff_mant_common:.2f}\times10^{{{common_exp}}}$",
        transform=trans_x,
        ha="center",
        va="top",
        fontsize=8.6,
        color=DK_RD2_color,
    )

    ax.text(
        1.985,
        -0.014,
        r"$G_{ab}^{max}=2G_0$",
        transform=trans_x,
        ha="center",
        va="top",
        fontsize=8.0,
        color="red",
    )

    # ------------------------------------------------------------------
    # Panel C — projection residual
    # ------------------------------------------------------------------
    ax_res = fig.add_subplot(gs[2])

    # Different alpha, z-order and markers help distinguish nearly
    # identical residual curves without changing the computed data.
    residual_style_map = {
        "Photons": dict(color="#00c8ff", linewidth=1.8, alpha=0.95, zorder=3),
        "Neutrinos": dict(color="#ff9f1a", linewidth=1.6, alpha=0.80, zorder=2),
        "Protons": dict(color="#32cd32", linewidth=1.4, alpha=0.75, zorder=1),
    }

    # Visual-only offsets are applied only in the plot so the nearly
    # identical residual curves can be distinguished. The CSV data remains unchanged.
    residual_visual_offsets = {
        "Photons": 2.5e-4,
        "Neutrinos": 0.0,
        "Protons": -2.5e-4,
    }

    marker_map = {
        "Photons": "o",
        "Neutrinos": "s",
        "Protons": "^",
    }

    for sp in species_plot:
        name = sp["name"]
        cfg = residual_style_map[name]

        y_raw = projection_residuals[name]
        y_plot = y_raw + residual_visual_offsets[name]

        ax_res.plot(
            T_vals,
            y_plot,
            linestyle=sp["linestyle"],
            linewidth=cfg["linewidth"],
            alpha=cfg["alpha"],
            color=cfg["color"],
            zorder=cfg["zorder"],
            label=f"{name} visual offset",
        )

        ax_res.scatter(
            T_vals[::16],
            y_plot[::16],
            s=11,
            marker=marker_map[name],
            alpha=0.75,
            color=cfg["color"],
            zorder=cfg["zorder"] + 3,
        )

    ax_res.text(
        0.02,
        0.90,
        "Small visual offsets applied only for readability.",
        transform=ax_res.transAxes,
        fontsize=6.8,
        ha="left",
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="0.6",
            alpha=0.85,
        ),
    )

    ax_res.axhline(
        0.0,
        color="black",
        linewidth=0.9,
        alpha=0.75,
    )

    ax_res.set_xscale("log")
    ax_res.set_xlabel("Temperature T [K]", fontsize=8.5)
    ax_res.set_ylabel(r"Projection residual", fontsize=8.5)
    ax_res.set_title(
        r"Panel C — Projection identity residual: "
        r"$\left[(G_{ab}^{eff}/G_0-1)/(G_{ab}^{max}/G_0-1)\right]-\alpha_{\rm DK}$",
        fontsize=9.5,
        pad=7,
    )

    ax_res.tick_params(axis="both", labelsize=7.2)
    ax_res.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.35)
    ax_res.legend(fontsize=6.8, loc="best")

    all_res = np.concatenate([
        np.asarray(v, dtype=float)[np.isfinite(v)]
        for v in projection_residuals.values()
    ])

    if all_res.size > 0:
        res_abs_max = float(np.nanmax(np.abs(all_res)))
        ylim = max(6e-4, res_abs_max * 5.0)
        ax_res.set_ylim(-ylim, ylim)

    # ------------------------------------------------------------------
    # Caption + footer
    # ------------------------------------------------------------------
    plt.figtext(
        0.5,
        0.060,
        rf"This figure demonstrates that high energy density does not necessarily imply high thermodynamic activity. "
        rf"The projected coupling uses $\alpha_{{\rm DK}}={alpha_DK:.4f}$ calibrated from DESI BAO + CC.",
        ha="center",
        va="center",
        fontsize=7.2,
    )

    plt.figtext(
        0.5,
        0.022,
        f"Figure: {file_fig} | Source: {Core_autor_text}",
        ha="center",
        va="center",
        fontsize=5.5,
        color=DK_RD2_color,
    )

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------
    residual_abs_max = float(np.nanmax(np.abs(all_res))) if all_res.size > 0 else np.nan

    stats = pd.DataFrame([
        {
            "model": "DK-RD2 stellar thermodynamic coupling",
            "alpha_DK": float(alpha_DK),
            "rd_mpc": float(rd_DK),
            "T_min_K": float(T_vals.min()),
            "T_max_K": float(T_vals.max()),
            "projection_residual_abs_max": residual_abs_max,
            "fit_source": "DESI BAO + CC calibration",
        }
    ])

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=5,
        fit_mode="STELLAR_ACTIVITY_PROJECTED_COUPLING_WITH_RESIDUAL",
        index=False,
    )

    fig.subplots_adjust(
        top=0.93,
        bottom=0.11,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, table_file, file_stats

def generate_figure06(
    sidm_notebook_path: str = "data.nb",
    *,
    n_r: int = 650,
):
    """
    # ============================================================
    # Figure 06 — Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2
    # ============================================================
    #
    # Purpose
    # -------
    # Reproduce the published compact-perturber density comparison from the
    # Mathematica notebook data.nb and add a DK-RD2 phenomenological
    # reinterpretation:
    #
    #   Observed compactness may be represented either as dense SIDM
    #   core-collapsed halos or as an effective gravitational amplification
    #   of a CDM-like baseline.
    #
    # Important
    # ---------
    # This figure does NOT claim that SIDM is wrong.
    # It shows that DK-RD2 offers an alternative effective-gravity
    # interpretation of the same compact-perturber phenomenology.
    # ============================================================
    """

    import re
    from pathlib import Path
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    file_fig = generate_evidence("image", 6)
    file_table = generate_evidence("table", 6)
    file_stats = file_table.replace(".csv", "_stats.csv")

    nb_path = Path(sidm_notebook_path)
    if not nb_path.exists():
        raise FileNotFoundError(
            f"Figure06 requires the Zenodo Mathematica notebook data.nb. "
            f"File not found: {sidm_notebook_path}"
        )

    nb_text = nb_path.read_text(encoding="utf-8", errors="ignore")

    # ------------------------------------------------------------
    # Official Figure 06 local color aliases
    # ------------------------------------------------------------
    dk_color = DK_RD2_color
    lcdm_color = LCDM_color
    dk_light = DK_LIGHT_color
    lcdm_light = LCDM_LIGHT_color
    data_color = DATA_color
    error_color = ERROR_color
    secondary_color = SECONDARY_color

    cdm_color = "black"
    sidm30_color = lcdm_color
    sidm50_color = lcdm_light
    sidm100_color = secondary_color

    b1938_color = "#4a4a4a"
    fornax_color = "#8a8a8a"
    gd1_color = "0.25"
    gd1_band_color = ERROR_color

    # ------------------------------------------------------------
    # Helper: extract Mathematica RowBox point arrays from data.nb.
    # ------------------------------------------------------------
    def _extract_nb_rowbox_pairs(var_name: str) -> np.ndarray:
        start = nb_text.find(var_name)
        if start < 0:
            raise ValueError(f"Variable {var_name} not found in {sidm_notebook_path}")

        candidate_vars = [
            "DMrhodataCDM",
            "DMrhodataS30",
            "DMrhodataS50",
            "DMrhodataS100",
            "LogLogPlot",
            "ListLogLogPlot",
        ]

        next_positions = [
            nb_text.find(v, start + 1)
            for v in candidate_vars
            if nb_text.find(v, start + 1) > start
        ]

        end = min(next_positions) if next_positions else start + 20000
        block = nb_text[start:end]

        pair_pattern = re.compile(
            r'RowBox\[\{"\{",\s*RowBox\[\{"([0-9.]+)"\s*,\s*",",\s*'
            r'(?:RowBox\[\{"([0-9.]+)"\s*,\s*"\*"\s*,\s*'
            r'RowBox\[\{"10"\s*,\s*"\^"\s*,\s*"(-?\d+)"\}\]\}\]|"([0-9.]+)")'
            r'\}\]\s*,\s*"\}"\}\]',
            re.S,
        )

        pairs = []
        for match in pair_pattern.finditer(block):
            x_val = float(match.group(1))

            if match.group(4) is not None:
                y_val = float(match.group(4))
            else:
                y_val = float(match.group(2)) * 10.0 ** int(match.group(3))

            pairs.append((x_val, y_val))

        if len(pairs) < 5:
            raise ValueError(
                f"Could not extract enough data points for {var_name}. "
                f"Check data.nb format."
            )

        arr = np.array(pairs, dtype=float)
        arr = arr[np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])]
        arr = arr[arr[:, 1] > 0.0]
        arr = arr[np.argsort(arr[:, 0])]
        return arr

    # ------------------------------------------------------------
    # 1. Extract published CDM / SIDM simulation curves
    # ------------------------------------------------------------
    cdm_data = _extract_nb_rowbox_pairs("DMrhodataCDM")
    sidm30_data = _extract_nb_rowbox_pairs("DMrhodataS30")
    sidm50_data = _extract_nb_rowbox_pairs("DMrhodataS50")
    sidm100_data = _extract_nb_rowbox_pairs("DMrhodataS100")

    # ------------------------------------------------------------
    # 2. Observational analytic profiles
    # ------------------------------------------------------------
    r = np.logspace(np.log10(0.0029), np.log10(0.50), int(n_r))

    rho0_pj = 4.3e7
    rt_pj = 0.149
    rho_powell = rho0_pj * rt_pj**4 / (r**2 * (r**2 + rt_pj**2))

    M_fornax = 1.0e6
    a_fornax = 0.020
    rho_fornax = (M_fornax / (2.0 * np.pi)) * a_fornax / (
        r * (r + a_fornax) ** 3
    )

    gd1_enc_params = np.array([
        [1.0000e-3, 8.7947e7],
        [1.2248e-2, 1.5250e8],
        [2.8573e-2, 9.6397e7],
        [3.1435e-2, 1.0663e7],
        [1.5736e-2, 4.6062e5],
        [1.0000e-3, 1.7580e5],
    ], dtype=float)

    gd1_profiles = []
    for a_i, M_i in gd1_enc_params:
        rho_i = (M_i / (2.0 * np.pi)) * a_i / (r * (r + a_i) ** 3)
        gd1_profiles.append(rho_i)

    gd1_profiles = np.array(gd1_profiles)
    rho_gd1_low = np.nanmin(gd1_profiles, axis=0)
    rho_gd1_high = np.nanmax(gd1_profiles, axis=0)
    rho_gd1_rep = np.sqrt(rho_gd1_low * rho_gd1_high)

    # ------------------------------------------------------------
    # 3. Interpolate CDM baseline from notebook data
    # ------------------------------------------------------------
    cdm_interp = interp1d(
        np.log(cdm_data[:, 0]),
        np.log(cdm_data[:, 1]),
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    rho_cdm_base = np.exp(cdm_interp(np.log(r)))
    rho_cdm_base = np.clip(rho_cdm_base, 1e-300, None)

    # ------------------------------------------------------------
    # 4. Required effective gravitational amplification
    # ------------------------------------------------------------
    mu_powell_required = rho_powell / rho_cdm_base
    mu_fornax_required = rho_fornax / rho_cdm_base
    mu_gd1_rep_required = rho_gd1_rep / rho_cdm_base
    mu_gd1_low_required = rho_gd1_low / rho_cdm_base
    mu_gd1_high_required = rho_gd1_high / rho_cdm_base

    def mu_DK_profile(r_in, A0, r_c, n):
        r_in = np.asarray(r_in, dtype=float)
        return 1.0 + A0 / (1.0 + (r_in / r_c) ** n)

    def fit_mu_profile(target_mu, label, r_min=0.003, r_max=0.12):
        fit_mask = (
            (r >= r_min)
            & (r <= r_max)
            & np.isfinite(target_mu)
            & (target_mu > 0.0)
        )

        if np.count_nonzero(fit_mask) < 10:
            raise RuntimeError(f"Not enough valid points to fit {label}")

        def objective(theta):
            A0, r_c, n = theta

            if A0 <= 0.0 or r_c <= 0.0 or n <= 0.0:
                return 1.0e99

            mu_model = mu_DK_profile(r[fit_mask], A0, r_c, n)

            if not np.all(np.isfinite(mu_model)):
                return 1.0e99

            residual = np.log10(mu_model) - np.log10(target_mu[fit_mask])
            return float(np.mean(residual**2))

        res_fit = minimize(
            objective,
            x0=np.array([20.0, 0.010, 2.0]),
            bounds=[
                (0.01, 5000.0),
                (0.001, 0.250),
                (0.20, 10.0),
            ],
            method="L-BFGS-B",
        )

        if not res_fit.success:
            raise RuntimeError(f"DK-RD2 fit failed for {label}: {res_fit.message}")

        A0_best, rc_best, n_best = [float(x) for x in res_fit.x]
        mu_fit = mu_DK_profile(r, A0_best, rc_best, n_best)
        rho_eff = rho_cdm_base * mu_fit

        return {
            "label": label,
            "A0": A0_best,
            "rc": rc_best,
            "n": n_best,
            "mu": mu_fit,
            "rho_eff": rho_eff,
            "log10_mse": float(objective(res_fit.x)),
            "r_min": r_min,
            "r_max": r_max,
        }

    fit_powell = fit_mu_profile(
        mu_powell_required,
        "B1938+666",
        r_min=0.003,
        r_max=0.12,
    )

    fit_fornax = fit_mu_profile(
        mu_fornax_required,
        "Fornax 6",
        r_min=0.003,
        r_max=0.08,
    )

    fit_gd1 = fit_mu_profile(
        mu_gd1_rep_required,
        "GD-1 representative",
        r_min=0.003,
        r_max=0.12,
    )

    # ------------------------------------------------------------
    # 5. Save evidence table
    # ------------------------------------------------------------
    table = pd.DataFrame({
        "r_kpc": r,

        "rho_CDM_interpolated_from_data_nb_Msun_kpc3": rho_cdm_base,

        "rho_Powell_B1938_PseudoJaffe_Msun_kpc3": rho_powell,
        "rho_Fornax6_Hernquist_Msun_kpc3": rho_fornax,
        "rho_GD1_low_Msun_kpc3": rho_gd1_low,
        "rho_GD1_high_Msun_kpc3": rho_gd1_high,
        "rho_GD1_representative_geometric_mean_Msun_kpc3": rho_gd1_rep,

        "mu_required_Powell_Geff_over_G0": mu_powell_required,
        "mu_required_Fornax_Geff_over_G0": mu_fornax_required,
        "mu_required_GD1_low_Geff_over_G0": mu_gd1_low_required,
        "mu_required_GD1_high_Geff_over_G0": mu_gd1_high_required,
        "mu_required_GD1_representative_Geff_over_G0": mu_gd1_rep_required,

        "mu_DKRD2_B1938_Geff_over_G0": fit_powell["mu"],
        "mu_DKRD2_Fornax_Geff_over_G0": fit_fornax["mu"],
        "mu_DKRD2_GD1_rep_Geff_over_G0": fit_gd1["mu"],

        "rho_DKRD2_B1938_effective_Msun_kpc3": fit_powell["rho_eff"],
        "rho_DKRD2_Fornax_effective_Msun_kpc3": fit_fornax["rho_eff"],
        "rho_DKRD2_GD1_rep_effective_Msun_kpc3": fit_gd1["rho_eff"],
    })

    dkrd2_to_csv(
        table,
        file_table,
        table_kind="FIG06_CORE_COLLAPSED_SIDM_DKRD2_INTERPRETATION",
        figure_id=6,
        strict=False,
        index=False,
        meta={
            "notebook_source": str(nb_path),
            "interpretation": (
                "Three separate DK-RD2 phenomenological effective-gravity "
                "profiles are fitted to B1938+666, Fornax 6, and a GD-1 "
                "representative curve. This preserves the SIDM interpretation "
                "while showing an alternative effective gravitational reading."
            ),
            "mu_DK_form": "1 + A0 / (1 + (r/r_c)^n)",
            "G0_source": "Core_G0 from DK_RD2_Core.py",
        },
    )

    stats_rows = []
    for fit in [fit_powell, fit_fornax, fit_gd1]:
        stats_rows.append({
            "model": f"DK-RD2 {fit['label']} compact amplification",
            "A0_inner_excess": fit["A0"],
            "mu_inner_limit_Geff_over_G0": 1.0 + fit["A0"],
            "r_c_kpc": fit["rc"],
            "r_c_pc": fit["rc"] * 1000.0,
            "transition_index_n": fit["n"],
            "log10_space_mse": fit["log10_mse"],
            "fit_radial_window_kpc": f"{fit['r_min']} to {fit['r_max']}",
            "fit_mode": "separate_phenomenological_effective_G_amplification",
        })

    stats = pd.DataFrame(stats_rows)

    dkrd2_stats_to_csv(
        stats,
        file_stats,
        figure_id=6,
        fit_mode="SEPARATE_PHENOMENOLOGICAL_EFFECTIVE_G_AMPLIFICATION",
        index=False,
    )

    # ------------------------------------------------------------
    # 6. Plot Figure 06
    # ------------------------------------------------------------
    fig = plt.figure(figsize=(13.5, 8.9), dpi=130)

    gs = gridspec.GridSpec(
        3,
        1,
        height_ratios=[1.55, 1.0, 1.55],
        hspace=0.18,
    )

    ax_den = fig.add_subplot(gs[0])
    ax_mu = fig.add_subplot(gs[1], sharex=ax_den)
    ax_eff = fig.add_subplot(gs[2], sharex=ax_den)

    # ------------------------------------------------------------
    # Panel A — published SIDM / observed profiles
    # ------------------------------------------------------------
    ax_den.loglog(
        cdm_data[:, 0],
        cdm_data[:, 1],
        color=cdm_color,
        linestyle="--",
        linewidth=1.7,
        label="CDM baseline",
    )

    ax_den.loglog(
        sidm30_data[:, 0],
        sidm30_data[:, 1],
        color=sidm30_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM30",
    )

    ax_den.loglog(
        sidm50_data[:, 0],
        sidm50_data[:, 1],
        color=sidm50_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM50",
    )

    ax_den.loglog(
        sidm100_data[:, 0],
        sidm100_data[:, 1],
        color=sidm100_color,
        linestyle="--",
        linewidth=1.45,
        label="SIDM100",
    )

    ax_den.loglog(
        r,
        rho_powell,
        color=b1938_color,
        linewidth=1.8,
        label="B1938+666 target",
    )

    ax_den.loglog(
        r,
        rho_fornax,
        color=fornax_color,
        linewidth=1.8,
        label="Fornax 6 target",
    )

    ax_den.fill_between(
        r,
        rho_gd1_low,
        rho_gd1_high,
        color=gd1_band_color,
        alpha=0.32,
        label="GD-1 envelope",
    )

    ax_den.loglog(
        r,
        rho_gd1_rep,
        color=gd1_color,
        linestyle="-.",
        linewidth=1.5,
        alpha=0.70,
        label="GD-1 representative",
    )

    ax_den.set_xlim(0.003, 0.50)
    ax_den.set_ylim(1e7, 2e11)
    ax_den.set_ylabel(
        r"$\rho(r)$ [$M_\odot\,{\rm kpc}^{-3}$]",
        fontsize=6.0,
    )

    ax_den.set_title(
        "Figure 06 — Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2",
        fontsize=11.8,
        pad=12,
    )

    ax_den.text(
        0.5,
        1.02,
        "Published SIDM profiles are preserved; DK-RD2 is tested as an effective gravitational amplification.",
        transform=ax_den.transAxes,
        ha="center",
        va="bottom",
        fontsize=7.2,
    )

    ax_den.grid(alpha=0.23, which="both")
    ax_den.legend(
        fontsize=6.4,
        loc="upper right",
        ncols=4,
        framealpha=0.92,
    )

    # ------------------------------------------------------------
    # Panel B — required and fitted mu profiles
    # ------------------------------------------------------------
    ax_mu.semilogx(
        r,
        mu_powell_required,
        color=b1938_color,
        linewidth=1.25,
        alpha=0.62,
        label=r"Required $\mu$: B1938+666",
    )

    ax_mu.semilogx(
        r,
        mu_fornax_required,
        color=fornax_color,
        linewidth=1.25,
        alpha=0.62,
        label=r"Required $\mu$: Fornax 6",
    )

    ax_mu.fill_between(
        r,
        mu_gd1_low_required,
        mu_gd1_high_required,
        color=gd1_band_color,
        alpha=0.28,
        label=r"Required $\mu$: GD-1 envelope",
    )

    ax_mu.semilogx(
        r,
        mu_gd1_rep_required,
        color=gd1_color,
        linestyle=":",
        linewidth=1.2,
        alpha=0.70,
        label=r"Required $\mu$: GD-1 representative",
    )

    ax_mu.semilogx(
        r,
        fit_powell["mu"],
        color=dk_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 fit: B1938+666",
    )

    ax_mu.semilogx(
        r,
        fit_fornax["mu"],
        color=dk_light,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 fit: Fornax 6",
    )

    ax_mu.semilogx(
        r,
        fit_gd1["mu"],
        color=secondary_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.15,
        label="DK-RD2 fit: GD-1 rep.",
    )

    ax_mu.axhline(
        1.0,
        color="black",
        linestyle=":",
        linewidth=1.0,
        alpha=0.80,
    )

    ax_mu.set_yscale("log")
    ax_mu.set_ylim(5e-3, 5e4)
    ax_mu.set_ylabel(
        r"$\mu_{\rm eff}(r)=G_{\rm eff}/G_0$",
        fontsize=6.0,
    )
    ax_mu.grid(alpha=0.23, which="both")
    ax_mu.legend(
        fontsize=5.4,
        loc="upper right",
        ncols=3,
        framealpha=0.92,
    )

    fit_text = (
        "Separate DK-RD2 compact fits:\n"
        rf"B1938: $\mu_0={1+fit_powell['A0']:.1f}$, "
        rf"$r_c={fit_powell['rc']*1000:.2f}$ pc, "
        rf"$n={fit_powell['n']:.2f}$\n"
        rf"Fornax: $\mu_0={1+fit_fornax['A0']:.1f}$, "
        rf"$r_c={fit_fornax['rc']*1000:.2f}$ pc, "
        rf"$n={fit_fornax['n']:.2f}$\n"
        rf"GD-1 rep.: $\mu_0={1+fit_gd1['A0']:.1f}$, "
        rf"$r_c={fit_gd1['rc']*1000:.2f}$ pc, "
        rf"$n={fit_gd1['n']:.2f}$"
    )

    ax_mu.text(
        0.025,
        0.93,
        fit_text,
        transform=ax_mu.transAxes,
        fontsize=5.4,
        va="top",
        ha="left",
        linespacing=1.05,
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.45",
            alpha=0.92,
        ),
    )

    # ------------------------------------------------------------
    # Panel C — DK effective-density proxies
    # ------------------------------------------------------------
    ax_eff.loglog(
        cdm_data[:, 0],
        cdm_data[:, 1],
        color=cdm_color,
        linestyle="--",
        linewidth=1.5,
        label="CDM baseline",
    )

    ax_eff.loglog(
        r,
        rho_powell,
        color=b1938_color,
        linewidth=1.2,
        alpha=0.55,
        label="B1938+666 target",
    )

    ax_eff.loglog(
        r,
        rho_fornax,
        color=fornax_color,
        linewidth=1.2,
        alpha=0.55,
        label="Fornax 6 target",
    )

    ax_eff.fill_between(
        r,
        rho_gd1_low,
        rho_gd1_high,
        color=gd1_band_color,
        alpha=0.28,
        label="GD-1 envelope",
    )

    ax_eff.loglog(
        r,
        fit_powell["rho_eff"],
        color=dk_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: B1938",
    )

    ax_eff.loglog(
        r,
        fit_fornax["rho_eff"],
        color=dk_light,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: Fornax",
    )

    ax_eff.loglog(
        r,
        fit_gd1["rho_eff"],
        color=secondary_color,
        linestyle=DK_RD2_linestyle,
        linewidth=2.2,
        label=r"CDM $\times\mu_{\rm DK}$: GD-1 rep.",
    )

    ax_eff.set_xlim(0.003, 0.50)
    ax_eff.set_ylim(1e7, 2e11)
    ax_eff.set_xlabel(
        "Radius r [kpc]",
        fontsize=5.6,
        labelpad=2,
    )
    ax_eff.set_ylabel(
        r"Effective density proxy [$M_\odot\,{\rm kpc}^{-3}$]",
        fontsize=6.0,
    )
    ax_eff.grid(alpha=0.23, which="both")
    ax_eff.legend(
        fontsize=5.4,
        loc="upper right",
        ncols=3,
        framealpha=0.92,
    )

    ax_eff.text(
        0.025,
        0.065,
        "Interpretation: the SIDM core-collapse result is preserved.\n"
        "DK-RD2 provides an equivalent effective-gravity description using system-dependent compact amplification.",
        transform=ax_eff.transAxes,
        fontsize=5.4,
        va="bottom",
        ha="left",
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="0.55",
            alpha=0.90,
        ),
    )

    for ax in [ax_den, ax_mu, ax_eff]:
        ax.tick_params(
            axis="both",
            labelsize=6.0,
        )

    # ------------------------------------------------------------
    # Footer
    # ------------------------------------------------------------
    fig.text(
        0.5,
        0.012,
        f"Figure: {file_fig} | Source: {Core_autor_text} |\n"
        f"External comparison data: Zenodo Mathematica notebook from http://dx.doi.org/10.5281/zenodo.19116269",
        ha="center",
        fontsize=5.2,
        color=dk_color,
    )

    fig.subplots_adjust(
        top=0.91,
        bottom=0.125,
        left=0.105,
        right=0.985,
    )

    fig.savefig(file_fig, bbox_inches="tight", dpi=300)
    plt.show()
    plt.close(fig)

    return file_fig, file_table, file_stats

def set_plot_text_defaults():
    """
    Ensure consistent, English-only plot styling across the project.
    """
    import matplotlib as mpl
    mpl.rcParams["axes.titlesize"] = 12
    mpl.rcParams["axes.labelsize"] = 11
    mpl.rcParams["legend.fontsize"] = 9
    mpl.rcParams["figure.titlesize"] = 13
    mpl.rcParams["font.family"] = "DejaVu Sans"  # or another available family
    # All labels/titles you set elsewhere should be English strings.


if __name__ == "__main__":
    import os

    Core_out_dir_path = "evidence"
    os.makedirs(Core_out_dir_path, exist_ok=True)

    set_plot_text_defaults()

    # ==============================================================
    # DK-RD2 ThermoGravity Framework — Figure selector
    #
    # Examples:
    #   exec_figs = ""        # run all paper figures
    #   exec_figs = "all"     # run all paper figures
    #   exec_figs = "1,2"     # run Figure 01 and Figure 02 only
    #   exec_figs = "1-4"     # run Figure 01..04
    #
    # Figure IDs for this paper:
    #   1 — Figure 01: Stellar thermodynamic activity vs gravitational confinement
    #   2 — Figure 02: Effective equation of state w(z): DK-RD2 vs ΛCDM
    #   3 — Figure 03: Type Ia SN distance modulus μ(z): DK-RD2 vs ΛCDM
    #   4 — Figure 04: Degravification / BOAT-like transient diagnostic
    #   5 — Figure 05: Falsifiability summary table: DK-RD2 vs ΛCDM
    #
    # Optional / inherited diagnostics:
    #   8 — Optional: Linear growth observable fσ8(z)
    #   9 — Optional: CLASS/CMB geometric consistency, if retained as appendix
    #
    # Removed from this paper main sequence:
    #   old Fig.04 H(z)        -> already covered in previous DK-RD2 paper
    #   old Fig.05 Lensing SIS -> already covered in previous DK-RD2 paper
    #   old SPHEREx proxy      -> not required here
    # ==============================================================

    ascii_art = r"""
 _____  _  __     _____  _____  ___    _____ _                               ____                 _ _
|  __ \| |/ /    |  __ \|  __ \|__ \  |_   _| |                             / __ \               (_) |
| |  | | ' /_____| |__) | |  | |  ) |   | | | |__   ___ _ __ _ __ ___   ___| |  | |_ __ __ ___   _| |_ _   _
| |  | |  <______|  _  /| |  | | / /    | | | '_ \ / _ \ '__| '_ ` _ \ / _ \ |  | | '__/ _` \ \ / / | __| | | |
| |__| | . \     | | \ \| |__| |/ /_    | | | | | |  __/ |  | | | | | | (_) | |__| | | | (_| |\ V /| | |_| |_| |
|_____/|_|\_\    |_|  \_\_____/|____|   |_| |_| |_|\___|_|  |_| |_| |_|\___/ \____/|_|  \__,_| \_/ |_|\__|\__, |
                                                                                                           __/ |
                                                                                                          |___/
"""
    print(ascii_art)
    print("What Is Gravity? From Emergent Mass to Thermodynamic–Geometric Dynamics")
    print("DK-RD2 ThermoGravity Framework")
    print(Core_autor_text)

    exec_figs = "all"

    print(
        "====================== Figure IDs for this paper ======================\n"
        "0  - Run all paper figures\n"
        "\n"
        "1  — Figure 01: DESI BAO calibration of α_DK and r_d\n"
        "2  — Figure 02: Statistical constraint on the DK projection factor\n"
        "3  — Figure 03: Thermodynamic–relativistic coupling structure\n"
        "4  — Figure 04: Diagnostic reconstruction of the effective equation of state w(z)\n"
        "5  — Figure 05: Thermodynamic activity vs projected gravitational coupling\n"
        "6  — Figure 06: Effective-gravity interpretation of compact halo structures\n"
        "\n"
        "X  — Exit or Ctrl+C"
    )
    print("selector examples:")
    print('  exec_figs = "" or "0" # run all paper figures')
    print('  exec_figs = "all"     # run all paper figures')
    print('  exec_figs = "1,2"     # run Figure 01 and Figure 02 only')
    print('  exec_figs = "1-5"     # run all main paper figures\n')

    print(f"Input selector value: {exec_figs!r}")
    user_spec = input(
        "Type a new selection (e.g. 1,2 or 1-5; 0 = all; ENTER = keep current): "
    ).strip()

    if user_spec == "":
        selection_spec = exec_figs
    elif user_spec == "0":
        selection_spec = "all"
    else:
        selection_spec = user_spec

    ALL_FIGS = {1, 2, 3, 4, 5, 6}
    OPTIONAL_FIGS = {} #{8, 9}
    AVAILABLE_FIGS = ALL_FIGS.union(OPTIONAL_FIGS)

    def _parse_exec_list(spec: str | None, available: set[int]) -> set[int]:
        """
        Parse execution specification into a set of figure numbers.

        Supported:
            None / "" / "all"  -> all main paper figures
            "1,2,5"            -> explicit list
            "1-5,8"            -> ranges plus optional diagnostics
        """
        if spec is None:
            return set(ALL_FIGS)

        s = str(spec).strip()

        if s == "" or s.lower() == "all":
            return set(ALL_FIGS)

        out: set[int] = set()
        parts = [p.strip() for p in s.split(",") if p.strip()]

        for p in parts:
            if "-" in p:
                a, b = [x.strip() for x in p.split("-", 1)]
                if a.isdigit() and b.isdigit():
                    lo, hi = int(a), int(b)
                    if lo > hi:
                        lo, hi = hi, lo
                    out.update(range(lo, hi + 1))
            else:
                if p.isdigit():
                    out.add(int(p))

        return out.intersection(available)

    RUN = _parse_exec_list(selection_spec, AVAILABLE_FIGS)

    print(ascii_art)
    print(f"▶ Exec selector = {selection_spec!r} → running: {sorted(RUN)}")
    if selection_spec == "all":
        print ("Generating Déjà vu figures...patience, the matrix... had been changed.")

    path = "data/DESI/bao_data/desi_bao_dr2/"
    desi_bao_mean_path = path + "desi_gaussian_bao_ALL_GCcomb_mean.txt"
    desi_bao_cov_path = path + "desi_gaussian_bao_ALL_GCcomb_cov.txt"

    if 1 in RUN:
        # ============================================================
        # Figure 01 — DESI Expansion and BAO Calibration of α_DK
        # ============================================================
        print("Figure01  DESI Expansion and BAO Calibration of α_DK")
        print("Generating figure... hang tight, this may take a while.")
        fig1_png, fig1_csv, fig1_stats, alpha_dk_best = generate_figure01(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 01 saved: {fig1_png}")
        print(f"✔ Table saved:     {fig1_csv}")

    # ============================================================
    # Figure 02 — Statistical Constraint on the DK Projection Factor
    # ============================================================
    if 2 in RUN:
        print("Figure 02: Statistical Constraint on the DK Projection Factor")
        print("Generating figure... hang tight, this may take a while.")
        print("... or more than a while...")
        fig2_img, fig2_table, fig2_stats = generate_figure02(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 02 saved: {fig2_img}")
        print(f"✔ Table saved:     {fig2_table}")
        print(f"✔ Stats saved:     {fig2_stats}")

    if 3 in RUN:
        print("Figure 03: DK-RD2 Coupling Structure")
        print("Generating figure... hang tight, this may take a while.")
        fig3_img, fig3_table, fig3_stats = generate_figure03(desi_bao_mean_path,desi_bao_cov_path)
        print(f"✔ Figure 03 saved: {fig3_img}")
        print(f"✔ Table saved:     {fig3_table}")
        print(f"✔ Stats saved:     {fig3_stats}")

    # ============================================================
    # Figure 04 — Effective Equation of State w_eff(z)
    # ============================================================
    if 4 in RUN:
        print("Figure 04: Explicit Diagnostic Reconstruction of w(z)")
        print("Generating figure... hang tight, this may take a while.")
        fig4_img, fig4_table, fig4_stats = generate_figure04(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 04 saved: {fig4_img}")
        print(f"✔ Table saved:     {fig4_table}")
        print(f"✔ Stats saved:     {fig4_stats}")

    # ------------------------------------------------------------------
    # Figure 05 — Stellar thermodynamic activity
    # ------------------------------------------------------------------
    if 5 in RUN:
        print("Figure 05: Stellar thermodynamic activity")
        print("Generating figure... hang tight, this may take a while.")
        print("A little more than a while...")
        fig5_img, fig5_table, fig5_stats = generate_figure05(
            desi_bao_mean_path=desi_bao_mean_path,
            desi_bao_cov_path=desi_bao_cov_path,
        )
        print(f"✔ Figure 05 saved: {fig5_img}")
        print(f"✔ Table saved:     {fig5_table}")
        print(f"✔ Stats saved:     {fig5_stats}")

    # ============================================================
    # Figure 06 — Effective Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2
    # ============================================================
    sidm_notebook_path = "data/data.nb" # from http://dx.doi.org/10.5281/zenodo.19116269
    if 6 in RUN:
        print("Figure 06: Effective Gravitational Interpretation of Core-Collapsed SIDM Halos in DK-RD2")
        print("External comparison data: Zenodo Mathematica notebook data.nb")
        print("data.nb download from http://dx.doi.org/10.5281/zenodo.19116269")
        print("Generating figure... hang tight, this may take a while.")
        fig6_img, fig6_table, fig6_stats = generate_figure06(
            sidm_notebook_path=sidm_notebook_path,
        )
        print(f"✔ Figure 06 saved: {fig6_img}")
        print(f"✔ Table saved:     {fig6_table}")
        print(f"✔ Stats saved:     {fig6_stats}")

    print("\nDone. DK-RD2 ThermoGravity figure generation completed.")
    banner = r"""
     _____  _  __     _____  _____ ___                                          
    |  __ \| |/ /    |  __ \|  __ \__ \                                         
    | |  | | ' /_____| |__) | |  | | ) |                                        
    | |  | |  <______|  _  /| |  | |/ /                                         
    | |__| | . \     | | \ \| |__| / /_                                         
    |_____/|_|\_\    |_|  \_\_____/____|      _____                 _ _         
    |__   __| |                              / ____|               (_) |        
       | |  | |__   ___ _ __ _ __ ___   ___ | |  __ _ __ __ ___   ___| |_ _   _ 
       | |  | '_ \ / _ \ '__| '_ ` _ \ / _ \| | |_ | '__/ _` \ \ / / | __| | | |
       | |  | | | |  __/ |  | | | | | | (_) | |__| | | | (_| |\ V /| | |_| |_| |
       |_|  |_|_|_|\___|_|  |_| |_| |_|\___/ \_____|_|  \__,_| \_/ |_|\__|\__, |
              |  ____|                                         | |         __/ |
              | |__ _ __ __ _ _ __ ___   _____      _____  _ __| | __     |___/ 
              |  __| '__/ _` | '_ ` _ \ / _ \ \ /\ / / _ \| '__| |/ /           
              | |  | | | (_| | | | | | |  __/\ V  V / (_) | |  |   <            
              |_|  |_|  \__,_|_| |_| |_|\___| \_/\_/ \___/|_|  |_|\_\           
    """
    print(banner)
    print("— GabE=mc² & Luludns -> ∞Ψ")
    print(Core_autor_text)
    print(f"Output Evidence files in: '{Core_out_dir_path}/'")
