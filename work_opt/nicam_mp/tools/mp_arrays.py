"""Metadata for the NICAM NSW6 cloud-microphysics mini-app dataset.

Single source of truth for:
  * ARRAYS  -- per-array role / units / one-line physical meaning
  * SCALARS -- the scalar parameters the MP step needs but that are NOT dumped
               as arrays (they come from the NICAM namelists / mod_cnst /
               mod_runconf setup and are constant for this configuration)
  * PROVENANCE -- where the dumped state came from

Consumed by dump_to_npy.py.  Editing this file changes both manifests.
"""

# ---------------------------------------------------------------------------
PROVENANCE = {
    "model": "NICAM (Fortran reference, OpenACC tree) -- kksref dump harness "
             "nicam/NICAM/nhm/driver/prg_kksref.f90, Stages C+E+F",
    "harness_binary": "nicam/NICAM/bin_cpu_local/kksref (CPU build, "
                      "nvfortran 26.3, -Mbyteswapio -Mfprelaxed -Mfpapprox)",
    "run_dir": "run/mp-dataset-realstep (config: "
               "run/mp-dataset-realstep/nhm_driver.cnf)",
    "resolution": "gl07rl00Az78: glevel=7, rlevel=0, 78 vertical layers "
                  "(+2 boundary halo levels), 10 diamond regions on 1 MPI rank",
    "state_chain": "day-10 restart -> ONE full dynamics_step (harness Stage E) "
                   "-> mp_driver (Stage F). This is exactly the model's own "
                   "sequence for production step 28801: the prg_driver.f90 "
                   "time loop is 'call dynamics_step; call physics_step' "
                   "(:588-589), physics updates the prognostic state directly, "
                   "and TIME_advance comes after physics -- so the dumped "
                   "mp_driver inputs are the state the production run hands "
                   "its step-28801 MP call.",
    "initial_state": "run/phase1/fortran-10d/restart_all_GL07RL00Az78.rgn*.nc "
                     "-- 10 model days of production full-physics spin-up "
                     "(28800 steps x 30 s), restart written at model time "
                     "864000 s",
    "model_date": "2018-01-11 00:00:00 UTC (TIME_CTIME = 6.3682848e10 s, "
                  "start 2018-01-01 + 10 days; TIME_advance runs after "
                  "physics_step, so step 28801's dynamics AND MP both see "
                  "this CTIME)",
    "mp_scheme": "MP_TYPE='NSW6' (Tomita 2008 six-category single-moment bulk "
                 "microphysics), RAIN_TYPE='COLD', opt_2moment_water=.false.",
    "physics_context": "The harness ran under the PRODUCTION configuration "
                       "(RD=MSTRNX, SF=LOUIS, TB=MYNN2, GWD=MCFARLANE, "
                       "LAND=MATSIRO, OCEAN=MIXEDLAYER, CP=NONE), so runconf/"
                       "diagvar were set up exactly as in the run that wrote "
                       "the restart. Only mp_driver itself was CALLED (the "
                       "harness replicates the physics_step clip block rather "
                       "than executing mod_physics); the MP chain reads none "
                       "of the other *_TYPE flags (verified: no reference in "
                       "mod_mp_driver/mod_mp_nsw6/mod_precip_transport/"
                       "mod_satadjust), so the mini-app may run with them set "
                       "to 'NONE'.",
    "column_locality": "The whole MP chain is column-local: no horizontal "
                       "operators, no halo exchange, no MPI, no NetCDF. Every "
                       "array is indexed (column, level[, tracer], region).",
}

# Dimension names used in the manifests -------------------------------------
DIMS = {
    "gall": "ADM_gall_in = 16641 = 129^2 -- horizontal columns of one region "
            "as handed to the physics (the 'in' grid: interior 128x128 plus "
            "one shared row and column on the i=1 / j=1 edges). Columns are "
            "completely independent for MP.",
    "k":    "ADM_kall = 80 -- vertical levels, k=1 and k=80 are the boundary "
            "(halo) levels; the active atmosphere is k=2..79 "
            "(ADM_kmin=2, ADM_kmax=79). k increases upward.",
    "v":    "TRC_VMAX = 6 -- tracer index: 1=qv 2=qc 3=qr 4=qi 5=qs 6=qg.",
    "l":    "ADM_lall = 10 -- diamond regions on this rank (l is ALWAYS the "
            "last/slowest dimension, so region 0 is a contiguous prefix).",
    "phase": "precipitation phase: 1 = liquid (rain), 2 = solid "
             "(snow+graupel+ice).",
    "sfc":  "single surface level (ADM_KNONE = 1).",
    "param": "mpF_params only: 1 = dt, 2 = ct.",
}

# ---------------------------------------------------------------------------
# role: metric | input | input-derived | scalar | output
# (name, role, units, description)
_A = [
    # --- geometry / vertical metric: constant in time, IN to all 3 entry points
    ("mpF_z",       "metric", "m",
     "Height above mean sea level at full (layer-centre) levels, GRD_vz(:,:,:,GRD_Z)."),
    ("mpF_zh",      "metric", "m",
     "Height above mean sea level at half (interface) levels, GRD_vz(:,:,:,GRD_ZH)."),
    ("mpF_phi",     "metric", "m2/s2",
     "Geopotential g*z at full levels (VMTR_PHI); used for the potential-energy budget."),
    ("mpF_gam2",    "metric", "1",
     "gamma^2 = (r/a)^2 vertical metric factor at full levels (= 1/VMTR_RGAM^2)."),
    ("mpF_gam2h",   "metric", "1",
     "gamma^2 at half levels (VMTR_GAM2H)."),
    ("mpF_gsgam2",  "metric", "1",
     "sqrt(G) * gamma^2 at full levels (VMTR_GSGAM2): the Jacobian that converts "
     "rho <-> rhog."),
    ("mpF_gsgam2h", "metric", "1",
     "sqrt(G) * gamma^2 at half levels (VMTR_GSGAM2H)."),
    ("mpF_ix",      "metric", "1",
     "x-component of the local eastward unit vector at the cell centre (GMTR_p IX)."),
    ("mpF_iy",      "metric", "1",
     "y-component of the local eastward unit vector (GMTR_p IY)."),
    ("mpF_iz",      "metric", "1",
     "z-component of the local eastward unit vector (GMTR_p IZ)."),
    ("mpF_jx",      "metric", "1",
     "x-component of the local northward unit vector (GMTR_p JX)."),
    ("mpF_jy",      "metric", "1",
     "y-component of the local northward unit vector (GMTR_p JY)."),
    ("mpF_jz",      "metric", "1",
     "z-component of the local northward unit vector (GMTR_p JZ)."),

    # --- aerosol / CCN inputs
    ("mpF_unccn",       "input", "1/m3",
     "Cloud condensation nuclei (aerosol) number concentration, diagvar I_UNCCN. "
     "Uniform 1.0e8 in this dataset."),
    ("mpF_aop_dnovrn",  "input", "1",
     "factor_n_aie: aerosol indirect-effect multiplier on Nc, diagvar I_AOP_DNOVRN. "
     "Uniform 1.0 in this dataset (so NSW6 runs with Nc = Nc_def = 50 /cm3)."),

    # --- scalars
    ("mpF_params",  "scalar", "s",
     "Two doubles: [1] dt = TIME_DTL, the physics/MP time step [s]; "
     "[2] ct = TIME_CTIME, absolute model time [s] on the NICAM calendar."),

    # --- ENTRY POINT 1: mp_driver -- INPUTS
    ("mpF_drv_rhog_in",   "input", "kg/m3",
     "rho * sqrt(G) * gamma^2: metric-weighted air density (prognostic)."),
    ("mpF_drv_rhogvx_in", "input", "kg/m2/s",
     "rhog * vx: metric-weighted momentum, Cartesian x."),
    ("mpF_drv_rhogvy_in", "input", "kg/m2/s",
     "rhog * vy: metric-weighted momentum, Cartesian y."),
    ("mpF_drv_rhogvz_in", "input", "kg/m2/s",
     "rhog * vz: metric-weighted momentum, Cartesian z."),
    ("mpF_drv_rhogw_in",  "input", "kg/m2/s",
     "rhog * w at HALF levels: metric-weighted vertical momentum."),
    ("mpF_drv_rhoge_in",  "input", "J/m3",
     "rhog * e_internal: metric-weighted internal energy."),
    ("mpF_drv_rhogq_in",  "input", "kg/m3",
     "rhog * q(v): metric-weighted tracer mass, v = qv,qc,qr,qi,qs,qg."),
    ("mpF_drv_rho_in",    "input", "kg/m3",  "Air density (diagnostic)."),
    ("mpF_drv_pre_in",    "input", "Pa",     "Air pressure (diagnostic)."),
    ("mpF_drv_tem_in",    "input", "K",      "Air temperature (diagnostic)."),
    ("mpF_drv_vx_in",     "input", "m/s",    "Wind, Cartesian x component."),
    ("mpF_drv_vy_in",     "input", "m/s",    "Wind, Cartesian y component."),
    ("mpF_drv_vz_in",     "input", "m/s",    "Wind, Cartesian z component."),
    ("mpF_drv_w_in",      "input", "m/s",    "Vertical velocity at HALF levels."),
    ("mpF_drv_q_in",      "input", "kg/kg",
     "Tracer mass mixing ratios, v = qv,qc,qr,qi,qs,qg."),

    # --- ENTRY POINT 1: mp_driver -- OUTPUTS
    ("mpF_drv_rhog_out",   "output", "kg/m3",   "rhog after one full mp_driver call."),
    ("mpF_drv_rhogvx_out", "output", "kg/m2/s", "rhogvx after mp_driver."),
    ("mpF_drv_rhogvy_out", "output", "kg/m2/s", "rhogvy after mp_driver."),
    ("mpF_drv_rhogvz_out", "output", "kg/m2/s", "rhogvz after mp_driver."),
    ("mpF_drv_rhogw_out",  "output", "kg/m2/s", "rhogw after mp_driver (half level)."),
    ("mpF_drv_rhoge_out",  "output", "J/m3",    "rhoge after mp_driver."),
    ("mpF_drv_rhogq_out",  "output", "kg/m3",   "rhogq after mp_driver."),
    ("mpF_drv_rho_out",    "output", "kg/m3",   "rho after mp_driver."),
    ("mpF_drv_pre_out",    "output", "Pa",      "pre after mp_driver."),
    ("mpF_drv_tem_out",    "output", "K",       "tem after mp_driver."),
    ("mpF_drv_vx_out",     "output", "m/s",     "vx after mp_driver."),
    ("mpF_drv_vy_out",     "output", "m/s",     "vy after mp_driver."),
    ("mpF_drv_vz_out",     "output", "m/s",     "vz after mp_driver."),
    ("mpF_drv_w_out",      "output", "m/s",     "w after mp_driver (half level)."),
    ("mpF_drv_q_out",      "output", "kg/kg",   "q after mp_driver."),
    ("mpF_drv_precip",        "output", "kg/m2/s",
     "Surface precipitation flux; index 1 = liquid (rain), 2 = solid (snow+graupel+ice)."),
    ("mpF_drv_precip_lh",     "output", "W/m2",
     "Latent-heat flux carried out of the column by precipitation."),
    ("mpF_drv_precip_rhoein", "output", "W/m2",
     "Internal-energy flux carried out by precipitation."),
    ("mpF_drv_precip_rhophi", "output", "W/m2",
     "Potential-energy flux carried out by precipitation."),
    ("mpF_drv_precip_rhokin", "output", "W/m2",
     "Kinetic-energy flux carried out by precipitation."),
    ("mpF_drv_precip_energy", "output", "W/m2",
     "Total energy flux carried out by precipitation (sum of the four above)."),
    ("mpF_drv_cldfrac",       "output", "1",
     "Diagnosed cloud fraction GDCFRC in [0,1] (input was all-zero: CP_TYPE='NONE')."),
    ("mpF_drv_qhyd",          "output", "kg/kg",
     "Diagnosed cloud-water content GDCLW handed to radiation."),

    # --- ENTRY POINT 2: precip_transport_new -- EXTRA INPUTS
    # (the state inputs are the mpF_drv_*_out arrays above)
    ("mpF_pt_Vt",     "input", "m/s",
     "Terminal fall velocity per tracer from mp_terminal_velocity (argument WPREC). "
     "Nonzero only for v=3,4,5,6 (qr,qi,qs,qg); zeroed on the k=1 and k=80 boundary "
     "levels. Dumped as an INPUT -- a mini-app reads it, it does not have to "
     "recompute it."),
    ("mpF_pt_qd_in",  "input-derived", "kg/kg",
     "Dry-air mass ratio qd = 1 - sum_v q(v), from THRMDYN_qd on the post-mp_driver state."),
    ("mpF_pt_rgs",    "input-derived", "1",
     "gam2/gsgam2 = 1/sqrt(G) at full levels."),
    ("mpF_pt_rgsh",   "input-derived", "1",
     "gam2h/gsgam2h = 1/sqrt(G) at half levels."),

    # --- ENTRY POINT 2 -- OUTPUTS
    ("mpF_pt_rhog_out",   "output", "kg/m3",   "rhog after precip_transport_new."),
    ("mpF_pt_rhogvx_out", "output", "kg/m2/s", "rhogvx after precip_transport_new."),
    ("mpF_pt_rhogvy_out", "output", "kg/m2/s", "rhogvy after precip_transport_new."),
    ("mpF_pt_rhogvz_out", "output", "kg/m2/s", "rhogvz after precip_transport_new."),
    ("mpF_pt_rhogw_out",  "output", "kg/m2/s", "rhogw after precip_transport_new."),
    ("mpF_pt_rhoge_out",  "output", "J/m3",    "rhoge after precip_transport_new."),
    ("mpF_pt_rhogq_out",  "output", "kg/m3",   "rhogq after precip_transport_new."),
    ("mpF_pt_rho_out",    "output", "kg/m3",   "rho after precip_transport_new."),
    ("mpF_pt_tem_out",    "output", "K",       "tem after precip_transport_new."),
    ("mpF_pt_pre_out",    "output", "Pa",      "pre after precip_transport_new."),
    ("mpF_pt_q_out",      "output", "kg/kg",   "q after precip_transport_new."),
    ("mpF_pt_qd_out",     "output", "kg/kg",   "qd after precip_transport_new."),
    ("mpF_pt_precip",        "output", "kg/m2/s",
     "Surface precipitation flux from this single call; 1 = liquid, 2 = solid."),
    ("mpF_pt_precip_rhoe",   "output", "W/m2",
     "Internal-energy flux removed by the falling hydrometeors."),
    ("mpF_pt_precip_lh",     "output", "W/m2",
     "Latent-heat flux removed by the falling hydrometeors."),
    ("mpF_pt_precip_rhophi", "output", "W/m2",
     "Potential-energy flux removed by the falling hydrometeors."),
    ("mpF_pt_precip_rhokin", "output", "W/m2",
     "Kinetic-energy flux removed by the falling hydrometeors."),
    ("mpF_pt_precip_flux",   "output", "kg/m2/s",
     "Vertical profile of the precipitation flux; index 1 = liquid, 2 = solid."),

    # --- ENTRY POINT 3: SATURATION_adjustment -- OUTPUTS
    # (inputs are the mpF_pt_*_out state + mpF_gsgam2 + mpF_pt_qd_out)
    ("mpF_sat_rhoge_out", "output", "J/m3",  "rhoge after SATURATION_adjustment."),
    ("mpF_sat_rhogq_out", "output", "kg/m3", "rhogq after SATURATION_adjustment."),
    ("mpF_sat_tem_out",   "output", "K",     "tem after SATURATION_adjustment."),
    ("mpF_sat_q_out",     "output", "kg/kg", "q after SATURATION_adjustment."),
]

ARRAYS = {name: {"role": role, "units": units, "desc": desc}
          for (name, role, units, desc) in _A}

# ---------------------------------------------------------------------------
# Scalars the MP step needs that are NOT in any .bin array.
# Values are those actually in force in the harness run (read back from
# run/mp-dataset/msg.pe000000 and mod_runconf.f90 / mod_cnst.f90).
SCALARS = {
    "time": {
        "dt": {"value": 30.0, "units": "s",
               "meaning": "TIME_DTL, the MP/physics time step (also mpF_params[1])"},
        "ct": {"value": 6.3682848e10, "units": "s",
               "meaning": "TIME_CTIME, absolute model time = 2018-01-11 00:00:00 "
                          "(also mpF_params[2])"},
        "MP_DIV_NUM": {"value": 1, "units": "1",
                       "meaning": "number of MP sub-steps inside one mp_driver call; "
                                  "with 1 the sub-step dt equals dt"},
    },
    "dimensions": {
        "ADM_gall_in": {"value": 16641, "units": "1", "meaning": "columns per region (129^2)"},
        "ADM_kall":    {"value": 80,    "units": "1", "meaning": "vertical levels incl. 2 halo levels"},
        "ADM_kmin":    {"value": 2,     "units": "1", "meaning": "lowest active level (1-based)"},
        "ADM_kmax":    {"value": 79,    "units": "1", "meaning": "highest active level (1-based)"},
        "ADM_lall":    {"value": 10,    "units": "1", "meaning": "regions on this rank"},
        "TRC_VMAX":    {"value": 6,     "units": "1", "meaning": "number of tracers"},
    },
    "tracer_indices": {
        "_note": "1-based Fortran indices into the v dimension of q / rhogq / Vt. "
                 "RAIN_TYPE='COLD' -> six water tracers, no aerosol/chemistry tracers.",
        "I_QV": 1, "I_QC": 2, "I_QR": 3, "I_QI": 4, "I_QS": 5, "I_QG": 6,
        "NQW_STR": 1, "NQW_END": 6, "NQW_MAX": 6,
    },
    "heat_capacities": {
        "_note": "EIN_TYPE='SIMPLE' -> CVW(v)=CVdry for all six water tracers, "
                 "CPW(I_QV)=CPdry and CPW(v)=CVdry for the condensed species. "
                 "Indexed NQW_STR..NQW_END = 1..6.",
        "CVW": [717.6400000000001] * 6,
        "CPW": [1004.69, 717.6400000000001, 717.6400000000001,
                717.6400000000001, 717.6400000000001, 717.6400000000001],
        "CVW_units": "J/kg/K", "CPW_units": "J/kg/K",
    },
    "cnst": {
        "_note": "mod_cnst values read by the MP chain (CNSTPARAM not specified -> "
                 "NICAM defaults).",
        "CNST_PI":     {"value": 3.141592653589793,     "units": "1"},
        "CNST_EGRAV":  {"value": 9.806649999999999,     "units": "m/s2",   "meaning": "gravity"},
        "CNST_RAIR":   {"value": 287.05,                "units": "J/kg/K", "meaning": "Rdry"},
        "CNST_RVAP":   {"value": 461.51,                "units": "J/kg/K", "meaning": "Rvap"},
        "CNST_CP":     {"value": 1004.69,               "units": "J/kg/K", "meaning": "CPdry"},
        "CNST_CV":     {"value": 717.6400000000001,     "units": "J/kg/K", "meaning": "CVdry"},
        "CNST_CPV":    {"value": 1846.06,               "units": "J/kg/K", "meaning": "CPvap"},
        "CNST_CVV":    {"value": 1384.55,               "units": "J/kg/K", "meaning": "CVvap"},
        "CNST_CL":     {"value": 4217.8,                "units": "J/kg/K", "meaning": "specific heat, liquid water"},
        "CNST_CI":     {"value": 2106.0,                "units": "J/kg/K", "meaning": "specific heat, ice"},
        "CNST_EPSV":   {"value": 0.6219800221013629,    "units": "1",      "meaning": "Rdry/Rvap"},
        "CNST_EPSVT":  {"value": 0.6077686814143877,    "units": "1",      "meaning": "1/EPSV - 1"},
        "CNST_LH0":    {"value": 2500840.0,             "units": "J/kg",   "meaning": "LHV0, vaporization at 0 C"},
        "CNST_LHS0":   {"value": 2834500.0,             "units": "J/kg",   "meaning": "sublimation at 0 C"},
        "CNST_LHF0":   {"value": 333660.0,              "units": "J/kg",   "meaning": "fusion at 0 C"},
        "CNST_LH00":   {"value": 3148680.781,           "units": "J/kg",   "meaning": "vaporization at 0 K"},
        "CNST_LHS00":  {"value": 2905502.611,           "units": "J/kg",   "meaning": "sublimation at 0 K"},
        "CNST_LHF00":  {"value": -243178.17,            "units": "J/kg",   "meaning": "fusion at 0 K"},
        "CNST_PSAT0":  {"value": 610.7,                 "units": "Pa",     "meaning": "saturation vapour pressure at 0 C"},
        "CNST_DWATR":  {"value": 1000.0,                "units": "kg/m3",  "meaning": "density of water"},
        "CNST_DICE":   {"value": 916.8,                 "units": "kg/m3",  "meaning": "density of ice"},
        "CNST_PRE00":  {"value": 100000.0,              "units": "Pa",     "meaning": "reference pressure"},
        "CNST_TEM00":  {"value": 273.15,                "units": "K",      "meaning": "reference temperature (0 C)"},
        "CNST_TMELT":  {"value": 273.15,                "units": "K",      "meaning": "melting point"},
        "CNST_UNDEF":  {"value": -9.9999e30,            "units": "1",      "meaning": "undefined value"},
        "CNST_EPS_ZERO": {"value": 2.220446049250313e-16, "units": "1"},
        "LHV": {"value": 2500840.0, "units": "J/kg",
                "meaning": "mod_runconf LHV (= CNST_LHV0, EIN_TYPE='SIMPLE')"},
        "LHS": {"value": 2834500.0, "units": "J/kg", "meaning": "mod_runconf LHS (= CNST_LHS0)"},
        "LHF": {"value": 333660.0,  "units": "J/kg", "meaning": "mod_runconf LHF (= CNST_LHF0)"},
    },
    "runconf": {
        "_note": "&RUNCONFPARAM switches. Only the first five change the MP "
                 "code path; the *_TYPE flags below them are the DUMP RUN's "
                 "(production) values, recorded for provenance -- the MP chain "
                 "reads none of them, so the mini-app runs with them = 'NONE'.",
        "MP_TYPE": "NSW6", "RAIN_TYPE": "COLD", "EIN_TYPE": "SIMPLE",
        "PRCIP_TRN_ECORRECT": "KIN2EIN", "opt_2moment_water": False,
        "CP_TYPE": "NONE", "TB_TYPE": "MYNN2", "RD_TYPE": "MSTRNX",
        "SF_TYPE": "LOUIS", "GWD_TYPE": "MCFARLANE", "CH_TYPE": "NONE",
        "AE_TYPE": "NONE", "LAND_TYPE": "MATSIRO", "OCEAN_TYPE": "MIXEDLAYER",
    },
    "nm_mp_driver_init": {
        "_note": "&NM_MP_DRIVER_INIT (defaults).",
        "OPT_RADIUS_EXPLICIT": False, "OPT_VOLUME_EXPLICIT": True,
        "CONST_RE_CLD": 8.0e-6, "CONST_RE_RAIN": 5.0e-4, "CONST_RE_ICE": 4.0e-5,
        "CONST_RE_SNOW": 5.0e-4, "CONST_RE_GRAUPEL": 5.0e-4, "CONST_RE_HAIL": 5.0e-4,
        "TSICE": 273.15, "TWICE": 258.15,
    },
    "nm_mp_nsw6": {
        "_note": "&NM_MP_NSW6 (defaults) -- the NSW6 tuning constants.",
        "N0R": 8.0e6, "N0S": 3.0e6, "N0G": 4.0e6,
        "RHO_G": 400.0, "RHO_S": 100.0, "C_D": 0.6, "CR": 130.0, "CS": 4.84,
        "EIW": 1.0, "ERW": 1.0, "ESW": 1.0, "EGW": 1.0, "ERI": 1.0, "ESI": 1.0,
        "EGI": 0.1, "ESR": 1.0, "EGR": 1.0, "EGS": 1.0,
        "GAMMA_SACR": 0.025, "GAMMA_GACS": 0.09, "MI": 4.19e-13, "NC_DEF": 50.0,
        "PSAUT_BETA0": 6.0e-3, "GAMMA_SAUT": 0.06, "QICRT_SAUT": 0.0,
        "BETA_GAUT": 1.0e-3, "GAMMA_GAUT": 0.09, "QSCRT_GAUT": 6.0e-4,
        "F1R": 0.78, "F2R": 0.27, "F1S": 0.65, "F2S": 0.39, "F1G": 0.78, "F2G": 0.27,
        "A_FRZ": 0.66, "B_FRZ": 100.0,
        "PRECIP_SCHEME_TYPE": "Flux-Semilag_new",
        "PRECIP_TRANSPORT_TYPE": "4WATER",
        "INC_PGAUT": False, "OPT_EXPLICIT_ICEGEN": False, "OPT_INDIR": False,
        "CNST_V_TERM_QI": -9.9999e30, "QR_AUT_ACC_TYPE": "Default",
        "ROH_FLAG": True, "OPT_SCLW": False, "CRIT_SCLW_W": 1.0e-2,
        "OPT_BETA_SAUT": False, "DS_CUTOFF": 1.5e-4, "DG_CUTOFF": 1.5e-4,
        "OPT_VT_NDW6": False,
    },
    "saturation_adjustment": {
        "_note": "Entry point 3 is called exactly as mod_mp_nsw6 does it.",
        "SATURATION_setrange": [100.0, 90.0],
        "ice_adjust": False,
    },
    "precip_transport_flags": {
        "_note": "precipitating_flag(v) passed to precip_transport_new "
                 "(entry point 2), 1-based tracer index.",
        "precipitating_flag": [False, False, True, True, True, True],
    },
}
