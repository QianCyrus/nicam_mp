!-------------------------------------------------------------------------------
!> Module run configuration
!!
!! @par Description
!!          This module is for managing run configuration
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_runconf
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: runconf_setup
  public :: runconf_set_trc_name

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: RUN_TYPE = ''

  !---< Component Selector >---

  !--- Dynamics
  integer,                public :: NON_HYDRO_ALPHA    = 1 ! Nonhydrostatic/hydrostatic flag
  integer,                public :: DYN_DIV_NUM        = 1
  character(len=H_SHORT), public :: TRC_ADV_TYPE       = 'MIURA2004'
  character(len=H_SHORT), public :: NDIFF_LOCATION     = 'IN_LARGE_STEP2'
  character(len=H_SHORT), public :: TRC_ADV_LOCATION   = 'IN_DYN_DIV_LOOP' ! or ‘OUT_DYN_DIV_LOOP’
  logical,                public :: TRC_APPLY_LIMITER  = .true.
  character(len=H_SHORT), public :: TRC_NEG_FIX        = 'OFF'
  logical,                public :: FLAG_NUDGING       = .false.
  logical,                public :: FLAG_LBC           = .false. ! [ADD] 20210322 Junya D-NICAM
  real(RP),               public :: CORIOLIS_PARAM     = 0.0_RP
  logical,                public :: first_layer_remedy = .true.
  logical,                public :: reduce_comm_tracer = .false.

  !--- Physics
  logical,                public :: run_phys_step      = .true.
  character(len=H_SHORT), public :: RAIN_TYPE          = 'AUTO'
  logical,                public :: opt_2moment_water  = .false.
  character(len=H_SHORT), public :: EIN_TYPE           = 'EXACT'
  character(len=H_SHORT), public :: PRCIP_TRN_ECORRECT = 'KIN2EIN'

  character(len=H_SHORT), public :: CP_TYPE            = 'NONE'
  character(len=H_SHORT), public :: SC_TYPE            = 'NONE'
  character(len=H_SHORT), public :: MC_TYPE            = 'NONE'
  character(len=H_SHORT), public :: MP_TYPE            = 'NONE'
  character(len=H_SHORT), public :: RD_TYPE            = 'NONE'
  character(len=H_SHORT), public :: SF_TYPE            = 'LOUIS'
  character(len=H_SHORT), public :: ROUGHNESS_SEA_TYPE = 'YQW'
  character(len=H_SHORT), public :: OCEAN_TYPE         = 'NONE'
  character(len=H_SHORT), public :: RIVER_TYPE         = 'NONE'
  character(len=H_SHORT), public :: LAND_TYPE          = 'NONE'
  character(len=H_SHORT), public :: LAKE_TYPE          = 'NONE'
  character(len=H_SHORT), public :: TB_TYPE            = 'NONE'
  character(len=H_SHORT), public :: EV_TYPE            = 'NONE' ! eddy viscosity
  character(len=H_SHORT), public :: AE_TYPE            = 'NONE'
  character(len=H_SHORT), public :: CH_TYPE            = 'NONE'
  character(len=H_SHORT), public :: GWD_TYPE(2)        = ['NONE', 'NONE']

  character(len=H_SHORT), public :: FR_TYPE            = 'NONE'
  character(len=H_SHORT), public :: AF_TYPE            = 'NONE'
  character(len=H_SHORT), public :: RAD_TYPE           = 'ISCCP'

  character(len=H_SHORT), public :: LLOUT_TYPE         = 'NONE'
  character(len=H_SHORT), public :: OUT_FILE_TYPE      = 'DEFAULT'
  character(len=H_SHORT), public :: NDG_TYPE           = 'EXTDATA'
  character(len=H_SHORT), public :: AHE_TYPE           = 'NONE' ! [add] 2023.6 for AHE (jh230052)

  integer,                public :: MP_DIV_NUM  = 1
  integer,                public :: SFC_DIV_NUM = 1
  integer,                public :: TB_DIV_NUM  = 1

  !---< thermodynamical setting >---

  real(RP), public, allocatable :: CVW(:)
  real(RP), public, allocatable :: CPW(:)

  real(RP), public              :: LHV
  real(RP), public              :: LHF
  real(RP), public              :: LHS

  real(RP), public, allocatable :: CVW_bin(:) ! for SBM [add] Y.Sato 20171227
  real(RP), public, allocatable :: CPW_bin(:) ! for SBM [add] Y.Sato 20171227

  !---< tracer ID setting >---

  integer,  public            :: PRG_vmax        ! total number of prognostic variables
  integer,  public, parameter :: PRG_vmax0  = 6

  integer,  public, parameter :: I_RHOG     =  1 ! Density x G^1/2
  integer,  public, parameter :: I_RHOGVX   =  2 ! Density x G^1/2 x Horizontal velocity (X-direction)
  integer,  public, parameter :: I_RHOGVY   =  3 ! Density x G^1/2 x Horizontal velocity (Y-direction)
  integer,  public, parameter :: I_RHOGVZ   =  4 ! Density x G^1/2 x Horizontal velocity (Z-direction)
  integer,  public, parameter :: I_RHOGW    =  5 ! Density x G^1/2 x Vertical   velocity
  integer,  public, parameter :: I_RHOGE    =  6 ! Density x G^1/2 x Energy
  integer,  public, parameter :: I_RHOGQstr =  7 ! tracers
  integer,  public            :: I_RHOGQend = -1 !

  character(len=H_SHORT), public :: PRG_name(PRG_vmax0) ! short name  of prognostic variables
  character(len=H_MID),   public :: PRG_desc(PRG_vmax0) ! description of prognostic variables
  character(len=H_SHORT), public :: PRG_unit(PRG_vmax0) ! unit        of prognostic variables

  data PRG_name / 'rhog', 'rhogvx', 'rhogvy', 'rhogvz', 'rhogw', 'rhoge' /

  data PRG_desc / 'Density * G^1/2',                     &
                  'Density * G^1/2 * H-Velocity(X_DIR)', &
                  'Density * G^1/2 * H-Velocity(Y_DIR)', &
                  'Density * G^1/2 * H-Velocity(Z_DIR)', &
                  'Density * G^1/2 * V-Velocity',        &
                  'Density * G^1/2 * Internal Energy'    /

  data PRG_unit / 'kg/m3', 'kg/m2/s', 'kg/m2/s', 'kg/m2/s', 'kg/m2/s', 'J/m3' /

  integer,  public            :: DIAG_vmax       ! total number of diagnostic variables
  integer,  public, parameter :: DIAG_vmax0 = 6

  integer,  public, parameter :: I_pre      =  1 ! Pressure
  integer,  public, parameter :: I_tem      =  2 ! Temperature
  integer,  public, parameter :: I_vx       =  3 ! Horizontal velocity (X-direction)
  integer,  public, parameter :: I_vy       =  4 ! Horizontal velocity (Y-direction)
  integer,  public, parameter :: I_vz       =  5 ! Horizontal velocity (Z-direction)
  integer,  public, parameter :: I_w        =  6 ! Vertical   velocity
  integer,  public, parameter :: I_qstr     =  7 ! tracers
  integer,  public            :: I_qend     = -1 !

  character(len=H_SHORT), public :: DIAG_name(DIAG_vmax0) ! short name  of diagnostic variables
  character(len=H_MID),   public :: DIAG_desc(DIAG_vmax0) ! description of diagnostic variables
  character(len=H_SHORT), public :: DIAG_unit(DIAG_vmax0) ! unit        of diagnostic variables

  data DIAG_name / 'pre', 'tem', 'vx', 'vy', 'vz', 'w' /

  data DIAG_desc / 'Pressure',          &
                   'Temperature',       &
                   'H-Velocity(X_DIR)', &
                   'H-Velocity(Y_DIR)', &
                   'H-Velocity(Z_DIR)', &
                   'V-Velocity'         /

  data DIAG_unit / 'Pa', 'K', 'm/s', 'm/s', 'm/s', 'm/s' /

  integer,  public            :: TRC_vmax   =  0 ! total number of tracers

  character(len=H_SHORT), public, allocatable :: TRC_name(:) ! short name  of tracers
  character(len=H_MID),   public, allocatable :: WLABEL  (:) ! description of tracers
  character(len=H_SHORT), public, allocatable :: TRC_unit(:) ! unit        of tracers
  character(len=H_SHORT), public, allocatable :: TRC_type(:) ! type        of tracers

  integer,  public            :: NQW_MAX    =  0 ! subtotal number of water mass tracers
  integer,  public            :: NQW_STR    = -1 ! start index of water mass tracers
  integer,  public            :: NQW_END    = -1 ! end   index of water mass tracers
  integer,  public            :: I_QV       = -1 ! Water vapor
  integer,  public            :: I_QC       = -1 ! Cloud water
  integer,  public            :: I_QR       = -1 ! Rain
  integer,  public            :: I_QI       = -1 ! Ice
  integer,  public            :: I_QS       = -1 ! Snow
  integer,  public            :: I_QG       = -1 ! Graupel
  integer,  public            :: I_QH       = -1 ! Hail

  integer,  public, parameter :: HYDRO_MAX  = 7

  !--- double moment cloud
  integer,  public            :: NNW_MAX    =  0 ! subtotal number of water number tracers
  integer,  public            :: NNW_STR    = -1 ! start index of water number tracers
  integer,  public            :: NNW_END    = -1 ! end   index of water number tracers
  integer,  public            :: I_NC       = -1 ! Cloud water (number)
  integer,  public            :: I_NR       = -1 ! Rain        (number)
  integer,  public            :: I_NI       = -1 ! Ice         (number)
  integer,  public            :: I_NS       = -1 ! Snow        (number)
  integer,  public            :: I_NG       = -1 ! Graupel     (number)
  integer,  public            :: I_NH       = -1 ! Hail        (number)

  logical,  public, allocatable :: flag_diagnose_number(:)

  !--- spectral bin cloud [ADD] Y.Sato 20171227 [SBM]
  integer,  public            :: NUM_BIN    = 33 ! subtotal number of water mass tracers
  integer,  public            :: BIN_MAX    =  0 ! subtotal number of SBM tracers
  integer,  public            :: BIN_STR    = -1 ! start index of SBM tracers
  integer,  public            :: BIN_END    = -1 ! end   index of SBM tracers

  !--- PDF CLOUD SCHEME
  integer,  public            :: NPW_MAX    =  0 ! subtotal number of water number tracers
  integer,  public            :: NPW_STR    = -1 ! start index of water number tracers
  integer,  public            :: NPW_END    = -1 ! end   index of water number tracers
  integer,  public            :: I_MU2       = -1 ! Cloud water (number)
  integer,  public            :: I_MU3       = -1 ! Rain        (number)
  logical,  public            :: opt_PDF_CLOUD    = .false.

  !--- water isotope tracer
  character(len=H_SHORT), private :: ISOTOPE = 'OFF'
  character(len=H_SHORT), private :: NOFRAC  = 'OFF'
  logical,  public            :: opt_water_isotope    = .false.
  logical,  public            :: opt_disable_wisofrac = .false.

  integer,  public            :: ISO_MAX    =  0 ! subtotal number of isotpe tracers
  integer,  public            :: ISO_STR    = -1 ! start index     of isotpe tracers
  integer,  public            :: ISO_STR2   = -1 ! start index     of isotpe tracers
  integer,  public            :: ISO_END    = -1 ! end   index     of isotpe tracers

  integer,  public            :: ISO1_MAX   =  0 ! subtotal number of isotpe tracers (O18)
  integer,  public            :: ISO1_STR   = -1 ! start index     of isotpe tracers (O18)
  integer,  public            :: ISO1_END   = -1 ! end   index     of isotpe tracers (O18)
  integer,  public            :: ISO1_QV    = -1 ! Water vapor
  integer,  public            :: ISO1_QC    = -1 ! Cloud water
  integer,  public            :: ISO1_QR    = -1 ! Rain
  integer,  public            :: ISO1_QI    = -1 ! Ice
  integer,  public            :: ISO1_QS    = -1 ! Snow
  integer,  public            :: ISO1_QG    = -1 ! Graupel
  integer,  public            :: ISO1_QH    = -1 ! Hail

  integer,  public            :: ISO2_MAX   =  0 ! subtotal number of isotpe tracers (HD)
  integer,  public            :: ISO2_STR   = -1 ! start index     of isotpe tracers (HD)
  integer,  public            :: ISO2_END   = -1 ! end   index     of isotpe tracers (HD)
  integer,  public            :: ISO2_QV    = -1 ! Water vapor
  integer,  public            :: ISO2_QC    = -1 ! Cloud water
  integer,  public            :: ISO2_QR    = -1 ! Rain
  integer,  public            :: ISO2_QI    = -1 ! Ice
  integer,  public            :: ISO2_QS    = -1 ! Snow
  integer,  public            :: ISO2_QG    = -1 ! Graupel
  integer,  public            :: ISO2_QH    = -1 ! Hail

  !--- turbulent tracer
  integer,  public            :: NTB_MAX    =  0 ! subtotal number of turbulent tracers
  integer,  public            :: NTB_STR    = -1 ! start index of turbulent tracers
  integer,  public            :: NTB_END    = -1 ! end   index of turbulent tracers
  integer,  public            :: I_TKE      = -1 ! turbulence kinetic energy
  integer,  public            :: I_QKEp     = -1
  integer,  public            :: I_TSQp     = -1
  integer,  public            :: I_QSQp     = -1
  integer,  public            :: I_COVp     = -1
  integer,  public            :: I_TKEsgs   = -1

  !--- aerosol tracer
  logical,  public            :: opt_aerosol_forcing     = .false.
  logical,  public            :: opt_aerosolcomp_forcing = .false. ! [ADD] 22/05/27 D.Goto
  logical,  public            :: opt_offline_aerosol     = .false.
  logical,  public            :: opt_carb_on             = .true.
  logical,  public            :: opt_dust_on             = .true.
  logical,  public            :: opt_salt_on             = .true.
  logical,  public            :: opt_sulf_on             = .true.
  logical,  public            :: opt_dust_6bin           = .false. ! [ADD] 20/07/22 D.Goto

  integer,  public            :: NQA_MAX    =  0 ! subtotal number of aerosol mass tracers
  integer,  public            :: NQA_STR    = -1 ! start index     of aerosol mass tracers
  integer,  public            :: NQA_END    = -1 ! end   index     of aerosol mass tracers
  integer,  public            :: NDU_MAX    =  0 ! subtotal number of aerosol mass tracers (soil dust)
  integer,  public            :: NQDU_STR   = -1 ! start index     of aerosol mass tracers (soil dust)
  integer,  public            :: NQDU_END   = -1 ! end   index     of aerosol mass tracers (soil dust)
  integer,  public            :: NCB_MAX    =  0 ! subtotal number of aerosol mass tracers (carbonaceous)
  integer,  public            :: NQCB_STR   = -1 ! start index     of aerosol mass tracers (carbonaceous)
  integer,  public            :: NQCB_END   = -1 ! end   index     of aerosol mass tracers (carbonaceous)
  integer,  public            :: NSU_MAX    =  0 ! subtotal number of aerosol mass tracers (sulfate)
  integer,  public            :: NQSU_STR   = -1 ! start index     of aerosol mass tracers (sulfate)
  integer,  public            :: NQSU_END   = -1 ! end   index     of aerosol mass tracers (sulfate)
  integer,  public            :: NQSU_SO4   = -1 ! SO4(2-)
  integer,  public            :: NQSU_SO2   = -1 ! SO2(gas and liquid)
  integer,  public            :: NQSU_DMS   = -1 ! DMS
  integer,  public            :: NSA_MAX    =  0 ! subtotal number of aerosol mass tracers (sea salt)
  integer,  public            :: NQSA_STR   = -1 ! start index     of aerosol mass tracers (sea salt)
  integer,  public            :: NQSA_END   = -1 ! end   index     of aerosol mass tracers (sea salt)

  logical,  public            :: opt_incloud_aerosol = .false.

  integer,  public            :: NQAIN_MAX  =  0 ! subtotal number of in-cloud aerosol mass tracers
  integer,  public            :: NQAIN_STR  = -1 ! start index     of in-cloud aerosol mass tracers
  integer,  public            :: NQAIN_END  = -1 ! end   index     of in-cloud aerosol mass tracers
  integer,  public            :: NQDUIN_MAX =  0 ! subtotal number of in-cloud aerosol mass tracers (soil dust)
  integer,  public            :: NQDUIN_STR = -1 ! start index     of in-cloud aerosol mass tracers (soil dust)
  integer,  public            :: NQDUIN_END = -1 ! end   index     of in-cloud aerosol mass tracers (soil dust)
  integer,  public            :: NQCBIN_MAX =  0 ! subtotal number of in-cloud aerosol mass tracers (carbonaceous)
  integer,  public            :: NQCBIN_STR = -1 ! start index     of in-cloud aerosol mass tracers (carbonaceous)
  integer,  public            :: NQCBIN_END = -1 ! end   index     of in-cloud aerosol mass tracers (carbonaceous)
  integer,  public            :: NQSUIN_MAX =  0 ! subtotal number of in-cloud aerosol mass tracers (sulfate)
  integer,  public            :: NQSUIN_STR = -1 ! start index     of in-cloud aerosol mass tracers (sulfate)
  integer,  public            :: NQSUIN_END = -1 ! end   index     of in-cloud aerosol mass tracers (sulfate)
  integer,  public            :: NQSO4IN    = -1 ! in-cloud SO4(2-)
  integer,  public            :: NQSAIN_MAX =  0 ! subtotal number of in-cloud aerosol mass tracers (sea salt)
  integer,  public            :: NQSAIN_STR = -1 ! start index     of in-cloud aerosol mass tracers (sea salt)
  integer,  public            :: NQSAIN_END = -1 ! end   index     of in-cloud aerosol mass tracers (sea salt)

  !--- chemical tracer
  logical,  public            :: opt_online_ozone = .true.

  integer,  public            :: NCHEM_MAX  =  0 ! subtotal number of chemical (or general purpose) tracers
  integer,  public            :: NCHEM_STR  = -1 ! start index     of chemical (or general purpose) tracers
  integer,  public            :: NCHEM_END  = -1 ! end   index     of chemical (or general purpose) tracers

  logical,  public            :: opt_tracer_initall  = .false. ! Set initial data for tracers?
                                                               ! overwrite opt_age_tracer_init, opt_radon_init, &
                                                               !           opt_syn_tracer_init, opt_passive_init

  logical,  public            :: opt_age_tracer      = .false. ! Add age of air tracer?
  logical,  public            :: opt_age_tracer_init = .false. ! Set initial data?
  integer,  public            :: AGE_MAX    =  3
  integer,  public            :: AGE_STR    = -1
  integer,  public            :: AGE_END    = -1

  logical,  public            :: opt_syn_tracer      = .false. ! Add synthetic tracer?
  logical,  public            :: opt_syn_tracer_init = .false. ! Set initial data?
  integer,  public            :: SYN_MAX    =  6
  integer,  public            :: SYN_STR    = -1
  integer,  public            :: SYN_END    = -1

  logical,  public            :: opt_passive         = .true.  ! Add general passive tracer?
  logical,  public            :: opt_passive_init    = .false. ! Set initial data?
  integer,  public            :: PASSIVE_MAX = -1
  integer,  public            :: NQC_STR     = -1
  integer,  public            :: NQC_END     = -1

  logical,  public            :: opt_radon           = .false. ! Add Radon (simple)?
  logical,  public            :: opt_radon_init      = .false. ! Set initial data?
  integer,  public            :: RDN_MAX    =  1
  integer,  public            :: RDN_STR    = -1
  integer,  public            :: RDN_END    = -1

  real(RP), public            :: RDCAY = 2.1E-6_RP ! decay rate of Radon-222 [1/s]

  !--- TM tracer
  character(len=H_SHORT), public :: TM_RUN = 'OFF'

  integer,  public            :: TM_TRC_vmax = 0
  integer,  public            :: TM_STR      = -999
  integer,  public            :: TM_END      = -999

  !---< other setting >---

  !--- cloud tracer numbers for radiation-microphysics coupling
  integer,  public, parameter :: Nparticle_mp  = 6
  integer,  public, parameter :: I_mp_QC       = 1
  integer,  public, parameter :: I_mp_QR       = 2
  integer,  public, parameter :: I_mp_QI       = 3
  integer,  public, parameter :: I_mp_QS       = 4
  integer,  public, parameter :: I_mp_QG       = 5
  integer,  public, parameter :: I_mp_QH       = 6

  character(len=H_SHORT), public :: particle_name_mp(Nparticle_mp)

  data particle_name_mp / 'Qc', 'Qr', 'Qi', 'Qs', 'Qg', 'Qh' /

  !--- No. of band for rad.
  integer,  public, parameter :: NRBND         = 3
  integer,  public, parameter :: NRBND_VIS     = 1
  integer,  public, parameter :: NRBND_NIR     = 2
  integer,  public, parameter :: NRBND_IR      = 3

  !--- direct/diffuse
  integer,  public, parameter :: NRDIR         = 2
  integer,  public, parameter :: NRDIR_DIRECT  = 1
  integer,  public, parameter :: NRDIR_DIFFUSE = 2

  !--- No. of cloud radiative forcing category (clear sky,all sky)
  integer,  public, parameter :: NCRF          = 2

  !--- No. of aerosol category
  integer,  public            :: KAPCL         = 7

  !--- for ISCCP
  integer,  public, parameter :: NTAU_ISCCP    = 7
  integer,  public, parameter :: NPRES_ISCCP   = 7

  !--- roughness  parameter
  integer,  public, parameter :: NTYPE_Z0      = 3
  integer,  public, parameter :: N_Z0M         = 1
  integer,  public, parameter :: N_Z0H         = 2
  integer,  public, parameter :: N_Z0E         = 3

  !--- solar insolation and cloud setting
  logical,  public            :: opt_rad_partialcloud   = .false.
  logical,  public            :: opt_rad_clearsky       = .false.
  logical,  public            :: opt_rad_global_tropics = .false.
  real(SP), public            :: opt_rad_fix_lon ! [rad]
  real(SP), public            :: opt_rad_fix_lat ! [rad]

  !--- ideal surface flux setting
  real(SP), public            :: opt_tb_fix_sh   ! [W/m2]
  real(SP), public            :: opt_tb_fix_lh   ! [W/m2]

  !--- for precipitation transport
  real(RP), public            :: W_VADV1D_MAX =   1.E0_RP
  real(RP), public            :: W_VADV1D_MIN = -40.E0_RP

  !--- VISIT ! [Add] 2024/8/3 A.Noda
  logical,  public            :: OPT_LECO             = .false.
  logical,  public            :: CONSTANT_DAILYMEAN   = .false. ! fix dailymean data?
  logical,  public            :: input_leco_dailymean = .false.
  logical,  public            :: LECO_ACROSS_UTC      = .false.
  logical,  public            :: LECO_DAILYMEAN_DONE  = .false.
  real(RP), public            :: tune_visit_rsfd      = 1.0_RP
  real(RP), public            :: tune_visit_clcv      = 1.0_RP
  real(RP), public            :: tune_visit_ta        = 1.0_RP
  real(RP), public            :: tune_visit_tsoil     = 1.0_RP
  real(RP), public            :: tune_visit_ts        = 1.0_RP
  real(RP), public            :: tune_visit_qa        = 1.0_RP
  real(RP), public            :: tune_visit_fprcc     = 1.0_RP
  real(RP), public            :: tune_visit_fprcl     = 1.0_RP
  real(RP), public            :: tune_visit_fsnwc     = 1.0_RP
  real(RP), public            :: tune_visit_fsnwl     = 1.0_RP
  real(RP), public            :: tune_visit_ua        = 1.0_RP
  real(RP), public            :: tune_visit_va        = 1.0_RP
  real(RP), public            :: tune_visit_pa        = 1.0_RP
  real(RP), public            :: tune_visit_bco2      = 1.0_RP
  real(RP), public            :: tune_visit_fndep     = 1.0_RP
  real(RP), public            :: tune_visit_wsoil     = 1.0_RP
  real(RP), public            :: tune_visit_rof       = 1.0_RP
  real(RP), public            :: tune_visit_wfps      = 1.0_RP
  logical,public,save :: OPT_LECO_VOC_SUBDAILY=.false.  ! 2025/6/17 A.Noda
  character(len=H_SHORT), public,save :: OPT_LECO_RSFD_TYPE='PAR' ! 2025/6/17 A.Noda
  character(len=H_SHORT), public,save :: OPT_LECO_PAR_TYPE = 'PAR_BND' ! 2025/6/17 A.Noda

  !--- for MIROC
  character(len=H_SHORT), public  :: MIROC_PHY_TYPE = '' ! MIROC6/7 (6 or 7 )SP/NSP (ON/OFF sprinters)
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine RUNCONF_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_cnst, only: &
       CNST_UNDEF4
    implicit none

    namelist / RUNCONFPARAM / &
       RUN_TYPE,                &
       NON_HYDRO_ALPHA,         &
       NDIFF_LOCATION,          &
       run_phys_step,           &
       RAIN_TYPE,               &
       RAD_TYPE ,               &
       EIN_TYPE,                &
       CORIOLIS_PARAM,          &
       first_layer_remedy,      &
       opt_rad_partialcloud,    &
       opt_rad_clearsky,        &
       opt_rad_global_tropics,  &
       opt_rad_fix_lon,         &
       opt_rad_fix_lat,         &
       opt_tb_fix_sh,           &
       opt_tb_fix_lh,           &
       TRC_ADV_TYPE,            &
       TRC_ADV_LOCATION,        &
       TRC_APPLY_LIMITER,       &
       TRC_NEG_FIX,             &
       OUT_FILE_TYPE,           &
       FLAG_NUDGING,            &
       RD_TYPE,                 &
       AF_TYPE,                 &
       FR_TYPE,                 &
       TB_TYPE,                 &
       EV_TYPE,                 &
       SF_TYPE,                 &
       MC_TYPE,                 &
       MP_TYPE,                 &
       CP_TYPE,                 &
       SC_TYPE,                 &
       GWD_TYPE,                &
       LAND_TYPE,               &
       OCEAN_TYPE,              &
       LAKE_TYPE,               &
       RIVER_TYPE,              &
       LLOUT_TYPE,              &
       NDG_TYPE,                &
       AHE_TYPE,                &
       ROUGHNESS_SEA_TYPE,      &
       PRCIP_TRN_ECORRECT,      &
       NUM_BIN,                 &
       opt_2moment_water,       &
       opt_PDF_CLOUD,           &
       opt_incloud_aerosol,     &
       AE_TYPE,                 &
       KAPCL,                   &
       opt_aerosol_forcing,     &
       opt_aerosolcomp_forcing, & ! [Add] 22/05/27 D.Goto
       opt_offline_aerosol,     &
       opt_carb_on,             &
       opt_dust_on,             &
       opt_salt_on,             &
       opt_sulf_on,             &
       opt_dust_6bin,           & ! [Add] 20/07/22 D.Goto
       CH_TYPE,                 &
       opt_online_ozone,        &
       opt_tracer_initall,      &
       opt_age_tracer,          &
       opt_age_tracer_init,     &
       AGE_MAX,                 &
       opt_syn_tracer,          &
       opt_syn_tracer_init,     &
       SYN_MAX,                 &
       opt_passive,             &
       opt_passive_init,        &
       PASSIVE_MAX,             &
       opt_radon,               &
       opt_radon_init,          &
       RDN_MAX,                 &
       RDCAY,                   &
       MP_DIV_NUM,              &
       TB_DIV_NUM,              &
       SFC_DIV_NUM,             &
       DYN_DIV_NUM,             &
       W_VADV1D_MAX,            &
       W_VADV1D_MIN,            &
       ISOTOPE,                 &
       NOFRAC,                  &
       opt_water_isotope,       &
       opt_disable_wisofrac,    &
       TM_RUN,                  &
       TM_TRC_vmax,             &
       OPT_LECO,                & ! [Add] 2024/8/3 A.Noda
       CONSTANT_DAILYMEAN,      & ! [Add] 2024/8/3 A.Noda
       input_leco_dailymean,    & ! [Add] 2024/8/3 A.Noda
       OPT_LECO_VOC_SUBDAILY,& ! 2025/6/17 A.Noda
       OPT_LECO_RSFD_TYPE   ,& !  ! 2025/6/17 A.Noda
       OPT_LECO_PAR_TYPE    ,& !  ! 2025/6/17 A.Noda
       tune_visit_rsfd,         & ! [Add] 2024/8/3 A.Noda
       tune_visit_clcv,         & ! [Add] 2024/8/3 A.Noda
       tune_visit_ta,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_tsoil,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_ts,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_qa,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_fprcc,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_fprcl,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_fsnwc,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_fsnwl,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_ua,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_va,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_pa,           & ! [Add] 2024/8/3 A.Noda
       tune_visit_bco2,         & ! [Add] 2024/8/3 A.Noda
       tune_visit_fndep,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_wsoil,        & ! [Add] 2024/8/3 A.Noda
       tune_visit_rof,          & ! [Add] 2024/8/3 A.Noda
       tune_visit_wfps,         & ! [Add] 2024/8/3 A.Noda
       MIROC_PHY_TYPE,          &
       reduce_comm_tracer

    integer  :: ierr
    !---------------------------------------------------------------------------

    opt_rad_fix_lon = CNST_UNDEF4
    opt_rad_fix_lat = CNST_UNDEF4
    opt_tb_fix_sh   = CNST_UNDEF4
    opt_tb_fix_lh   = CNST_UNDEF4

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[runconf]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RUNCONFPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** RUNCONFPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist RUNCONFPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist RUNCONFPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=RUNCONFPARAM)

    call RUNCONF_component_setup

    call RUNCONF_tracer_setup

    call RUNCONF_thermodyn_setup

    return
  end subroutine RUNCONF_setup

  !-----------------------------------------------------------------------------
  !> component check
  subroutine RUNCONF_component_setup
    !--- MINIAPP: mod_coupler (Jcup) is not reachable from the microphysics and
    !--- drags in the whole Coupler tree; CPL_jcup_active is forced .false. ---
!MINIAPP    use mod_coupler, only: &
!MINIAPP       CPL_jcup_active, &
!MINIAPP       CPL_exists_coco, &
!MINIAPP       CPL_exists_io,   &
!MINIAPP       CPL_exists_pair
    use mod_adm, only: &
       ADM_HGRID_SYSTEM, & ! [ADD] 20210322 Junya D-NICAM
       ADM_proc_stop
    implicit none
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '====== Component Settings ======'

    if ( TRC_APPLY_LIMITER ) then
       write(IO_FID_LOG,*) '--- Trace Advection : Run with Thuburn Limiter in MIURA2004 Advection'
    else
       write(IO_FID_LOG,*) '--- Trace Advection : Without Thuburn Limiter in MIURA2004 Advection'
    endif

    if ( RAIN_TYPE == 'AUTO' ) then
       write(IO_FID_LOG,*) '*** Rain Type is automatically set from MP_TYPE.'
       select case(MP_TYPE)
       case('NONE')
          RAIN_TYPE = 'DRY'
       case('LSC')
          RAIN_TYPE = 'CLOUD_PARAM'
       case('LSC_MIROC6')
          RAIN_TYPE = 'CLOUD_PARAM_I'
       case('KESSLER','G98')
          RAIN_TYPE = 'WARM'
       case('LIN','WSM6','NSW5','NSW6','NDW6','S10')
          RAIN_TYPE = 'COLD'
       case default
          write(*,*) 'xxx Not appropriate type. MP_TYPE = ', trim(MP_TYPE)
          call ADM_proc_stop
       end select
    endif
    write(IO_FID_LOG,*) '--- Rain Type : ', trim(RAIN_TYPE)

    if ( ISOTOPE == 'ON' ) then
       write(IO_FID_LOG,*) '--- Water Isotope : ON'
       opt_water_isotope    = .true. ! backward compatibility
    endif

    if ( NOFRAC == 'ON' ) then
       opt_disable_wisofrac = .true. ! backward compatibility
    endif

    if ( mod(TB_DIV_NUM,SFC_DIV_NUM) /= 0 ) then
       TB_DIV_NUM = SFC_DIV_NUM
       write(IO_FID_LOG,*) '*** FIX to TB_DIV_NUM = ', TB_DIV_NUM
    endif

    if ( opt_rad_partialcloud ) then
       if (      TB_TYPE == 'MYNN2'      &
            .OR. TB_TYPE == 'MYNN2.5'    &
            .OR. TB_TYPE == 'MYNN3'      &
            .OR. TB_TYPE == 'MYNN2_R1'   &
            .OR. TB_TYPE == 'MYNN2.5_R1' &
            .OR. TB_TYPE == 'MYNN3_R1'   ) then

          write(IO_FID_LOG,*) '--- Parctial Cloud in TB : Subgrid-scale cloud fraction calculated by MYNN'

       else
          write(IO_FID_LOG,*) 'xxx opt_rad_partialcloud option is not allowed for this TB_TYPE=', trim(TB_TYPE)
          call ADM_proc_stop
       endif

       if ( opt_rad_clearsky ) then
          write(IO_FID_LOG,*) 'xxx opt_rad_partialcloud option is not allowed when opt_rad_clearsky is on.'
          call ADM_proc_stop
       endif
    endif

    if ( CH_TYPE == 'PASSIVE' ) then ! Yamashita [Passive tracer]
       if ( opt_tracer_initall ) then ! force true all
          opt_age_tracer_init = .true.
          opt_syn_tracer_init = .true.
          opt_passive_init    = .true.
          opt_radon_init      = .true.
       endif

       if ( opt_age_tracer .AND. AGE_MAX < 1 ) then
             write(IO_FID_LOG,*) 'xxx AGE of Air tracer is ON but AGE_MAX < 1! stop.'
             call ADM_proc_stop
       endif

       if ( opt_syn_tracer .AND. SYN_MAX < 1 ) then
             write(IO_FID_LOG,*) 'xxx Synthetic tracer is ON but SYN_MAX < 1! stop.'
             call ADM_proc_stop
       endif

       if ( opt_passive .AND. PASSIVE_MAX < 1 ) then
             write(IO_FID_LOG,*) 'xxx Passive tracer is ON but PASSIVE_MAX < 1! stop.'
             call ADM_proc_stop
       endif

       if ( opt_radon .AND. RDN_MAX < 1 ) then
             write(IO_FID_LOG,*) 'xxx Radon tracer is ON but RDN_MAX < 1! stop.'
             call ADM_proc_stop
       endif
    else
       opt_age_tracer = .false.
       opt_syn_tracer = .false.
       opt_passive    = .false.
       opt_radon      = .false.
    endif

    if ( opt_age_tracer ) then
       write(IO_FID_LOG,*) '--- AGE of Air tracer is ON'
    endif
    if ( opt_syn_tracer ) then
       write(IO_FID_LOG,*) '--- Synthetic tracer is ON'
    endif
    if ( opt_passive ) then
       write(IO_FID_LOG,*) '--- Passive tracer is ON'
    endif
    if ( opt_radon ) then
       write(IO_FID_LOG,*) '--- Radon tracer is ON'
    endif

!MINIAPP    if ( CPL_jcup_active ) then
!MINIAPP       if ( CPL_exists_coco ) then
!MINIAPP          OCEAN_TYPE = 'COCO'
!MINIAPP          write(IO_FID_LOG,*) '--- Overwrite OCEAN_TYPE => COCO'
!MINIAPP       elseif( OCEAN_TYPE == 'COCO' ) then
!MINIAPP          write(*,*) 'xxx COCO is not alive. STOP!'
!MINIAPP          call ADM_proc_stop
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( CPL_exists_io ) then
!MINIAPP          LLOUT_TYPE = 'CoupledIO'
!MINIAPP          write(IO_FID_LOG,*) '--- Overwrite LLOUT_TYPE => CoupledIO'
!MINIAPP       elseif( LLOUT_TYPE == 'CoupledIO' ) then
!MINIAPP          write(*,*) 'xxx CoupledIO is not alive. STOP!'
!MINIAPP          call ADM_proc_stop
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( CPL_exists_pair ) then
!MINIAPP          ! do  nothing
!MINIAPP       elseif( NDG_TYPE == 'CPLDATA' ) then
!MINIAPP          write(*,*) 'xxx NICAM pair is not alive. STOP!'
!MINIAPP          call ADM_proc_stop
!MINIAPP       endif
!MINIAPP    endif

    if ( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! [ADD] 20210322 Junya D-NICAM
       write(IO_FID_LOG,*) '--- Lateral Boundary Treatment is ON for D-NICAM'
       FLAG_LBC = .true.
    else
       FLAG_LBC = .false.
    endif

    return
  end subroutine RUNCONF_component_setup

  !-----------------------------------------------------------------------------
  !> tracer setup
  subroutine RUNCONF_tracer_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_ch_includes, only: &
       CH_TRC_vmax, &
       CH_TRC_name, &
       CH_TRC_desc
    implicit none

    character(len=1) :: is1
    character(len=2) :: is
    character(len=4) :: is4

    integer  :: tmp_vmax
    integer  :: v
    !---------------------------------------------------------------------------

    !--- counting tracer
    TRC_vmax = 0

    !###########################################################################
    ! Tracer index management
    !###########################################################################

    !--- Mass tracer for Water
    if    ( RAIN_TYPE == 'DRY' ) then
       NQW_MAX = 1
       I_QV    = TRC_vmax + 1
    elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
       NQW_MAX = 2
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
    elseif ( RAIN_TYPE == 'CLOUD_PARAM_I' ) then
      NQW_MAX = 3
      I_QV    = TRC_vmax + 1
      I_QC    = TRC_vmax + 2
      I_QI    = TRC_vmax + 3
    elseif( RAIN_TYPE == 'WARM' ) then
       NQW_MAX = 3
       I_QV    = TRC_vmax + 1
       I_QC    = TRC_vmax + 2
       I_QR    = TRC_vmax + 3
    elseif( RAIN_TYPE == 'COLD' ) then
      NQW_MAX = 6
      I_QV    = TRC_vmax + 1
      I_QC    = TRC_vmax + 2
      I_QR    = TRC_vmax + 3
      I_QI    = TRC_vmax + 4
      I_QS    = TRC_vmax + 5
      I_QG    = TRC_vmax + 6
    elseif ( MIROC_PHY_TYPE(1:6) == 'MIROC6' ) then
      NQW_MAX = 3
      I_QV    = TRC_vmax + 1
      I_QC    = TRC_vmax + 2
      I_QI    = TRC_vmax + 3
    elseif ( MIROC_PHY_TYPE(1:6) == 'MIROC7' ) then
      NQW_MAX = 5
      I_QV    = TRC_vmax + 1
      I_QC    = TRC_vmax + 2
      I_QI    = TRC_vmax + 3
      ! NOTE MIROC7 is used I_QR=4, N_QR=5, I_QS=6, N_QS=7
      ! but NICAM is assumed sequentially water tracers
      I_QR    = TRC_vmax + 4
      I_QS    = TRC_vmax + 5
    else
       write(*,         *) 'xxx You must set RAIN_TYPE to DRY,CLOUD_PARAM,WARM,COLD or AUTO. STOP.'
       write(IO_FID_LOG,*) 'xxx You must set RAIN_TYPE to DRY,CLOUD_PARAM,WARM,COLD or AUTO. STOP.'
       call ADM_proc_stop
    endif
    NQW_STR  = TRC_vmax + min(1,NQW_MAX)
    NQW_END  = TRC_vmax + NQW_MAX
    TRC_vmax = TRC_vmax + NQW_MAX

    !--- Number tracer for Water
    if ( opt_2moment_water ) then
       if    ( RAIN_TYPE == 'DRY' ) then
          NNW_MAX = 0
       elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
          NNW_MAX = 1
          I_NC    = TRC_vmax + 1
       elseif( RAIN_TYPE == 'CLOUD_PARAM_I' ) then
          NNW_MAX = 2
          I_NC    = TRC_vmax + 1
          I_NI    = TRC_vmax + 2
       elseif( RAIN_TYPE == 'WARM' ) then
          NNW_MAX = 2
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
       elseif( RAIN_TYPE == 'COLD' ) then
          NNW_MAX = 5
          I_NC    = TRC_vmax + 1
          I_NR    = TRC_vmax + 2
          I_NI    = TRC_vmax + 3
          I_NS    = TRC_vmax + 4
          I_NG    = TRC_vmax + 5
       endif
    endif
    if ( MIROC_PHY_TYPE(1:6) == 'MIROC7' ) then
      ! NOTE MIROC7's NC/NI is last part of tracer
      NNW_MAX = 2 !
      I_NR = TRC_vmax + 1
      I_NS = TRC_vmax + 2
    end if

    NNW_STR  = TRC_vmax + min(1,NNW_MAX)
    NNW_END  = TRC_vmax + NNW_MAX
    TRC_vmax = TRC_vmax + NNW_MAX

    !--- PDF cloud scheme
    if ( opt_PDF_CLOUD .or. MIROC_PHY_TYPE(1:5) == 'MIROC' .or. MC_TYPE == 'HPC' ) then
      NPW_MAX = 2
      I_MU2 = TRC_vmax + 1
      I_MU3 = TRC_vmax + 2
      opt_PDF_CLOUD = .true.
    endif
    NPW_STR = TRC_vmax + min(1,NPW_MAX)
    NPW_END = TRC_vmax + NPW_MAX
    TRC_vmax = TRC_vmax + NPW_MAX

    !--- Tracer for turbulence
    if    ( TB_TYPE == 'MY2.5' ) then
       NTB_MAX = 1
       I_TKE   = TRC_vmax + 1
    elseif( TB_TYPE == 'MYNN2.5' .OR. TB_TYPE == 'MYNN2.5_R1' ) then
       NTB_MAX = 1
       I_QKEp  = TRC_vmax + 1
    elseif( TB_TYPE == 'MYNN3' .OR. TB_TYPE == 'MYNN3_R1' ) then
       NTB_MAX = 4
       I_QKEp  = TRC_vmax + 1
       I_TSQp  = TRC_vmax + 2
       I_QSQp  = TRC_vmax + 3
       I_COVp  = TRC_vmax + 4
    elseif( MIROC_PHY_TYPE(1:5) == 'MIROC' ) then
       NTB_MAX = 1
       I_QKEp  = TRC_vmax + 1
    endif

    if ( EV_TYPE == 'SMG15' ) then
       NTB_MAX  = NTB_MAX  + 1
       I_TKEsgs = TRC_vmax + NTB_MAX
    endif
    NTB_STR = TRC_vmax + min(1,NTB_MAX)
    NTB_END = TRC_vmax + NTB_MAX
    TRC_vmax = TRC_vmax + NTB_MAX

    !--- Tracer for aerosol
    if ( AE_TYPE /= 'NONE' ) then

       tmp_vmax = TRC_vmax

       ! soil dust
       if ( opt_dust_on ) then
          if (opt_dust_6bin) then  !! [ADD] 20/07/22 D.Goto
             NDU_MAX = 6
          else                    !! NICAM19 standard version
             NDU_MAX = 10
          endif
       endif
       NQDU_STR = tmp_vmax + min(1,NDU_MAX)
       NQDU_END = tmp_vmax + NDU_MAX
       tmp_vmax = tmp_vmax + NDU_MAX

       ! carbonaceous
       if ( opt_carb_on ) then
          NCB_MAX = 7
       endif
       NQCB_STR = tmp_vmax + min(1,NCB_MAX)
       NQCB_END = tmp_vmax + NCB_MAX
       tmp_vmax = tmp_vmax + NCB_MAX

       ! sulfate
       if ( opt_sulf_on ) then
          NSU_MAX = 3

          NQSU_SO4 = tmp_vmax + 1
          NQSU_SO2 = tmp_vmax + 2
          NQSU_DMS = tmp_vmax + 3
       endif
       NQSU_STR = tmp_vmax + min(1,NSU_MAX)
       NQSU_END = tmp_vmax + NSU_MAX
       tmp_vmax = tmp_vmax + NSU_MAX

       ! sea salt
       if ( opt_salt_on ) then
          NSA_MAX = 4
       endif
       NQSA_STR = tmp_vmax + min(1,NSA_MAX)
       NQSA_END = tmp_vmax + NSA_MAX
       tmp_vmax = tmp_vmax + NSA_MAX

       NQA_MAX = NDU_MAX + NCB_MAX + NSU_MAX + NSA_MAX

       if ( opt_incloud_aerosol ) then
          ! soil dust
          if ( opt_dust_on ) then
             NQDUIN_MAX = 10
          endif
          NQDUIN_STR = tmp_vmax + min(1,NQDUIN_MAX)
          NQDUIN_END = tmp_vmax + NQDUIN_MAX
          tmp_vmax   = tmp_vmax + NQDUIN_MAX

          ! carbonaceous
          if ( opt_carb_on ) then
             NQCBIN_MAX = 7
          endif
          NQCBIN_STR = tmp_vmax + min(1,NQCBIN_MAX)
          NQCBIN_END = tmp_vmax + NQCBIN_MAX
          tmp_vmax   = tmp_vmax + NQCBIN_MAX

          ! sea salt
          if ( opt_salt_on ) then
             NQSAIN_MAX = 4
          endif
          NQSAIN_STR = tmp_vmax + min(1,NQSAIN_MAX)
          NQSAIN_END = tmp_vmax + NQSAIN_MAX
          tmp_vmax   = tmp_vmax + NQSAIN_MAX

          ! sulfate
          if ( opt_sulf_on ) then
             NQSUIN_MAX = 1

             NQSO4IN = tmp_vmax + 1
          endif
          NQSUIN_STR = tmp_vmax + min(1,NQSUIN_MAX)
          NQSUIN_END = tmp_vmax + NQSUIN_MAX
          tmp_vmax   = tmp_vmax + NQSUIN_MAX

          NQAIN_MAX = NQDUIN_MAX + NQCBIN_MAX + NQSUIN_MAX + NQSAIN_MAX

       endif
    endif
    ! MIROC PHYSICS
    if ( MIROC_PHY_TYPE(7:8) == 'SP' ) then ! MIROC6SP/7SP
      tmp_vmax = TRC_vmax
      ! DUST
      NDU_MAX = 6
      NQDU_STR = tmp_vmax + min(1,NDU_MAX)
      NQDU_END = tmp_vmax + NDU_MAX
      tmp_vmax = tmp_vmax + NDU_MAX
      ! sea SALT
      NSA_MAX = 4
      NQSA_STR = tmp_vmax + min(1,NSA_MAX)
      NQSA_END = tmp_vmax + NSA_MAX
      tmp_vmax = tmp_vmax + NSA_MAX
      ! CARBON
      NCB_MAX = 6
      NQCB_STR = tmp_vmax + min(1,NCB_MAX)
      NQCB_END = tmp_vmax + NCB_MAX
      tmp_vmax = tmp_vmax + NCB_MAX
      ! SULFATE
      NSU_MAX = 3
      NQSU_STR = tmp_vmax + min(1,NSU_MAX)
      NQSU_END = tmp_vmax + NSU_MAX
      tmp_vmax = tmp_vmax + NSU_MAX

      NQA_MAX = NDU_MAX + NCB_MAX + NSU_MAX + NSA_MAX
    end if
    NQA_STR   = TRC_vmax + min(1,NQA_MAX)
    NQA_END   = TRC_vmax + NQA_MAX
    TRC_vmax  = TRC_vmax + NQA_MAX
    NQAIN_STR = TRC_vmax + min(1,NQAIN_MAX)
    NQAIN_END = TRC_vmax + NQAIN_MAX
    TRC_vmax  = TRC_vmax + NQAIN_MAX
    ! MIROC PHYSICS
    ! note: cloud number concentration
    if ( MIROC_PHY_TYPE(7:8) == 'SP' ) then ! MIROC6SP
      NNW_MAX = 2 ! not used but setted
      I_NC    = TRC_vmax + 1
      I_NI    = TRC_vmax + 2
      NNW_STR = TRC_VMAX + min(1,NNW_MAX)
      NNW_END = TRC_VMAX + NNW_MAX
      TRC_vmax = TRC_vmax + NNW_MAX
    end if

    !--- Tracer for water isotope
    if ( opt_water_isotope ) then
       if    ( RAIN_TYPE == 'DRY' ) then
          ISO_MAX = 2
          ISO1_QV = TRC_vmax + 1
          ISO2_QV = TRC_vmax + 2
       elseif( RAIN_TYPE == 'CLOUD_PARAM' ) then
          ISO_MAX = 4
          ISO1_QV = TRC_vmax + 1
          ISO1_QC = TRC_vmax + 2
          ISO2_QV = TRC_vmax + 3
          ISO2_QC = TRC_vmax + 4
       elseif( RAIN_TYPE == 'WARM' ) then
          ISO_MAX = 6
          ISO1_QV = TRC_vmax + 1
          ISO1_QC = TRC_vmax + 2
          ISO1_QR = TRC_vmax + 3
          ISO2_QV = TRC_vmax + 4
          ISO2_QC = TRC_vmax + 5
          ISO2_QR = TRC_vmax + 6
       elseif( RAIN_TYPE == 'COLD' ) then
          ISO_MAX = 12
          ISO1_QV = TRC_vmax + 1
          ISO1_QC = TRC_vmax + 2
          ISO1_QR = TRC_vmax + 3
          ISO1_QI = TRC_vmax + 4
          ISO1_QS = TRC_vmax + 5
          ISO1_QG = TRC_vmax + 6
          ISO2_QV = TRC_vmax + 7
          ISO2_QC = TRC_vmax + 8
          ISO2_QR = TRC_vmax + 9
          ISO2_QI = TRC_vmax + 10
          ISO2_QS = TRC_vmax + 11
          ISO2_QG = TRC_vmax + 12
       endif
    endif
    ISO_STR  = TRC_vmax + min(1,ISO_MAX)
    ISO_STR2 = ISO2_QV
    ISO_END  = TRC_vmax + ISO_MAX
    TRC_vmax = TRC_vmax + ISO_MAX

    !--- SBM Tracer [ADD] Y.Sato 20171227
    if ( MP_TYPE == 'S10' ) then
       if    ( RAIN_TYPE == 'WARM' ) then
          BIN_MAX = NUM_BIN
       elseif( RAIN_TYPE == 'COLD' ) then
          BIN_MAX = NUM_BIN * 7
       else
          write(*         ,*) 'xxx RAIN_TYPE should be WARM or COLD for MP_TYPE=S10! STOP.'
          write(IO_FID_LOG,*) 'xxx RAIN_TYPE should be WARM or COLD for MP_TYPE=S10! STOP.'
          call ADM_proc_stop
       endif
    endif
    BIN_STR  = TRC_vmax + min(1,BIN_MAX)
    BIN_END  = TRC_vmax + BIN_MAX
    TRC_vmax = TRC_vmax + BIN_MAX

    !--- TM Tracer
    if ( TM_RUN == 'ON' ) then ! Niwa [TM]
       TM_STR   = TRC_vmax + 1
       TM_END   = TRC_vmax + TM_TRC_vmax
       TRC_vmax = TRC_vmax + TM_TRC_vmax
    endif

    !--- Tracer for chemistry

    if ( CH_TYPE == 'CHASER_TROP' ) then
       NCHEM_MAX = CH_TRC_vmax
       NCHEM_STR = TRC_vmax + 1
       NCHEM_END = TRC_vmax + NCHEM_MAX
       TRC_vmax  = TRC_vmax + NCHEM_MAX
    elseif ( CH_TYPE == 'CHAPMAN' ) then
       NCHEM_MAX = 1
       NCHEM_STR = TRC_vmax + 1
       NCHEM_END = TRC_vmax + NCHEM_MAX
       TRC_vmax  = TRC_vmax + NCHEM_MAX
    endif

    if ( opt_age_tracer ) then
       AGE_STR  = TRC_vmax + 1
       AGE_END  = TRC_vmax + AGE_MAX
       TRC_vmax = TRC_vmax + AGE_MAX
    endif

    if ( opt_syn_tracer ) then
       SYN_STR  = TRC_vmax + 1
       SYN_END  = TRC_vmax + SYN_MAX
       TRC_vmax = TRC_vmax + SYN_MAX
    endif

    if ( opt_passive_init ) then
       NQC_STR  = TRC_vmax + 1
       NQC_END  = TRC_vmax + PASSIVE_MAX
       TRC_vmax = TRC_vmax + PASSIVE_MAX
    endif

    if ( opt_radon ) then
       RDN_STR  = TRC_vmax + 1
       RDN_END  = TRC_vmax + RDN_MAX
       TRC_vmax = TRC_vmax + RDN_MAX
    endif

    !###########################################################################
    ! Labeling
    !###########################################################################

    allocate( TRC_name(TRC_vmax) )
    allocate( WLABEL  (TRC_vmax) )
    allocate( TRC_unit(TRC_vmax) )
    allocate( TRC_type(TRC_vmax) ) ! Niwa [TM]

    ! default value
    do v = 1, TRC_vmax
       write(is4,'(I4.4)') v
       TRC_name(v) = 'fix_me_name_' //is4
       WLABEL  (v) = 'fix_me_label_'//is4
       TRC_unit(v) = '1'
       TRC_type(v) = ''
    enddo

    do v = 1, TRC_vmax
       if    ( v == I_QV ) then
          TRC_name(v) = 'qv'
          WLABEL  (v) = 'VAPOR'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QC ) then
          TRC_name(v) = 'qc'
          WLABEL  (v) = 'CLOUD'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QR ) then
          TRC_name(v) = 'qr'
          WLABEL  (v) = 'RAIN'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QI ) then
          TRC_name(v) = 'qi'
          WLABEL  (v) = 'ICE'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QS ) then
          TRC_name(v) = 'qs'
          WLABEL  (v) = 'SNOW'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QG ) then
          TRC_name(v) = 'qg'
          WLABEL  (v) = 'GRAUPEL'
          TRC_unit(v) = 'kg/kg'
       elseif( v == I_QH ) then
          TRC_name(v) = 'qh'
          WLABEL  (v) = 'HAIL'
          TRC_unit(v) = 'kg/kg'

       elseif( v == I_NC ) then
          TRC_name(v) = 'nc'
          WLABEL  (v) = 'CLOUD_NUM'
          TRC_unit(v) = '1/kg'
       elseif( v == I_NR ) then
          TRC_name(v) = 'nr'
          WLABEL  (v) = 'RAIN_NUM'
          TRC_unit(v) = '1/kg'
       elseif( v == I_NI ) then
          TRC_name(v) = 'ni'
          WLABEL  (v) = 'ICE_NUM'
          TRC_unit(v) = '1/kg'
       elseif( v == I_NS ) then
          TRC_name(v) = 'ns'
          WLABEL  (v) = 'SNOW_NUM'
          TRC_unit(v) = '1/kg'
       elseif( v == I_NG ) then
          TRC_name(v) = 'ng'
          WLABEL  (v) = 'GRAUPEL_NUM'
          TRC_unit(v) = '1/kg'
       elseif( v == I_NH ) then
          TRC_name(v) = 'nh'
          WLABEL  (v) = 'HAIL_NUM'
          TRC_unit(v) = '1/kg'

       elseif( v == I_MU2 ) then
            TRC_name(v) = 'mu2'
            WLABEL  (v) = 'PDFv'
            TRC_unit(v) = '1'
       elseif( v == I_MU3 ) then
            TRC_name(v) = 'mu3'
            WLABEL  (v) = 'PDFs'
            TRC_unit(v) = '1'

       elseif( v == I_TKE ) then
          TRC_name(v) = 'tke'
          WLABEL  (v) = 'TKE'
          TRC_unit(v) = 'm2/s2'
       elseif( v == I_QKEp ) then
          TRC_name(v) = 'qkep'
          WLABEL  (v) = 'QKEp'
          TRC_unit(v) = 'm2/s2'
       elseif( v == I_TSQp ) then
          TRC_name(v) = 'tsqp'
          WLABEL  (v) = 'TSQp'
          TRC_unit(v) = 'K2'
       elseif( v == I_QSQp ) then
          TRC_name(v) = 'qsqp'
          WLABEL  (v) = 'QSQp'
          TRC_unit(v) = '1'
       elseif( v == I_COVp ) then
          TRC_name(v) = 'covp'
          WLABEL  (v) = 'COVp'
          TRC_unit(v) = 'K'
       elseif( v == I_TKEsgs ) then
          TRC_name(v) = 'tke_sgs'
          WLABEL  (v) = 'TKE'
          TRC_unit(v) = 'm2/s2'

       elseif( NQDU_STR <= v .AND. v <= NQDU_END ) then
          write(is,'(I2.2)') v - NQDU_STR + 1
          TRC_name(v) = 'dust'//is
          WLABEL  (v) = 'DUST'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQCB_STR <= v .AND. v <= NQCB_END ) then
          write(is,'(I2.2)') v - NQCB_STR + 1
          TRC_name(v) = 'carb'//is
          WLABEL  (v) = 'CARBON'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQSU_STR <= v .AND. v <= NQSU_END ) then
          write(is,'(I2.2)') v - NQSU_STR + 1
          TRC_name(v) = 'sulf'//is
          WLABEL  (v) = 'SULFATE'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQSA_STR <= v .AND. v <= NQSA_END ) then
          write(is,'(I2.2)') v - NQSA_STR + 1
          TRC_name(v) = 'salt'//is
          WLABEL  (v) = 'SEASALT'//is
          TRC_unit(v) = 'kg/kg'

       elseif( NQDUIN_STR <= v .AND. v <= NQDUIN_END ) then
          write(is,'(I2.2)') v - NQDUIN_STR + 1
          TRC_name(v) = 'qcdust'//is
          WLABEL  (v) = 'DUST_INCLOUD'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQCBIN_STR <= v .AND. v <= NQCBIN_END ) then
          write(is,'(I2.2)') v - NQCBIN_STR + 1
          TRC_name(v) = 'qccb'//is
          WLABEL  (v) = 'CARBON_INCLOUD'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQSUIN_STR <= v .AND. v <= NQSUIN_END ) then
          write(is,'(I2.2)') v - NQSUIN_STR + 1
          TRC_name(v) = 'qcso4'//is
          WLABEL  (v) = 'SO4_INCLOUD'//is
          TRC_unit(v) = 'kg/kg'
       elseif( NQSAIN_STR <= v .AND. v <= NQSAIN_END ) then
          write(is,'(I2.2)') v - NQSAIN_STR + 1
          TRC_name(v) = 'qcsalt'//is
          WLABEL  (v) = 'SEASALT_INCLOUD'//is
          TRC_unit(v) = 'kg/kg'

       elseif( v == ISO1_QV ) then
          TRC_name(v) = 'qv_o18'
          WLABEL  (v) = 'O18_VAPOR'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO1_QC ) then
          TRC_name(v) = 'qc_o18'
          WLABEL  (v) = 'O18_CLOUD'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO1_QR ) then
          TRC_name(v) = 'qr_o18'
          WLABEL  (v) = 'O18_RAIN'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO1_QI ) then
          TRC_name(v) = 'qi_o18'
          WLABEL  (v) = 'O18_ICE'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO1_QS ) then
          TRC_name(v) = 'qs_o18'
          WLABEL  (v) = 'O18_SNOW'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO1_QG ) then
          TRC_name(v) = 'qg_o18'
          WLABEL  (v) = 'O18_GRAUPEL'
          TRC_unit(v) = 'kg/kg'

       elseif( v == ISO2_QV ) then
          TRC_name(v) = 'qv_hdo'
          WLABEL  (v) = 'HDO_VAPOR'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO2_QC ) then
          TRC_name(v) = 'qc_hdo'
          WLABEL  (v) = 'HDO_CLOUD'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO2_QR ) then
          TRC_name(v) = 'qr_hdo'
          WLABEL  (v) = 'HDO_RAIN'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO2_QI ) then
          TRC_name(v) = 'qi_hdo'
          WLABEL  (v) = 'HDO_ICE'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO2_QS ) then
          TRC_name(v) = 'qs_hdo'
          WLABEL  (v) = 'HDO_SNOW'
          TRC_unit(v) = 'kg/kg'
       elseif( v == ISO2_QG ) then
          TRC_name(v) = 'qg_hdo'
          WLABEL  (v) = 'HDO_GRAUPEL'
          TRC_unit(v) = 'kg/kg'

       elseif( BIN_STR <= v .AND. v <= BIN_END ) then ! Y.Sato [SBM]
          if    ( BIN_STR + NUM_BIN * 0 <= v .AND. v < BIN_STR + NUM_BIN * 1 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 0
             TRC_name(v) = 'sdf_liqu'//trim(is)
             WLABEL  (v) = 'SDF_Liq'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 1 <= v .AND. v < BIN_STR + NUM_BIN * 2 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 1
             TRC_name(v) = 'sdf_colm'//trim(is)
             WLABEL  (v) = 'SDF_Colmn'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 2 <= v .AND. v < BIN_STR + NUM_BIN * 3 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 2
             TRC_name(v) = 'sdf_plat'//trim(is)
             WLABEL  (v) = 'SDF_Plate'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 3 <= v .AND. v < BIN_STR + NUM_BIN * 4 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 3
             TRC_name(v) = 'sdf_dend'//trim(is)
             WLABEL  (v) = 'SDF_Dendrite'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 4 <= v .AND. v < BIN_STR + NUM_BIN * 5 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 4
             TRC_name(v) = 'sdf_snow'//trim(is)
             WLABEL  (v) = 'SDF_Snow'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 5 <= v .AND. v < BIN_STR + NUM_BIN * 6 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 5
             TRC_name(v) = 'sdf_grau'//trim(is)
             WLABEL  (v) = 'SDF_Graupel'//trim(is)
             TRC_unit(v) = 'kg/kg'
          elseif( BIN_STR + NUM_BIN * 6 <= v .AND. v < BIN_STR + NUM_BIN * 7 ) then
             write(is,'(I0)') v - BIN_STR + 1 - NUM_BIN * 6
             TRC_name(v) = 'sdf_hail'//trim(is)
             WLABEL  (v) = 'SDF_Hail'//trim(is)
             TRC_unit(v) = 'kg/kg'
          endif

       elseif( TM_STR <= v .AND. v <= TM_END ) then ! Niwa [TM]
          write(is,'(I2.2)') v - TM_STR + 1
          TRC_name(v) = 'TM_tmp'//is
          WLABEL  (v) = 'TM_tmp'//is
          TRC_unit(v) = 'kg/kg'
          TRC_type(v) = 'DEFAULT'

       elseif( AGE_STR <= v .AND. v <= AGE_END ) then ! Yamashita [Passive Tracer]
          write(is1,'(I1.1)') v - AGE_STR + 1
          TRC_name(v) = 'aoa'//is1
          WLABEL  (v) = 'age of air type'//is1
       elseif( SYN_STR <= v .AND. v <= SYN_END ) then ! Yamashita [Passive tracer]
          write(is,'(I2.2)') v - SYN_STR + 1
          TRC_name(v) = 'syn'//is
          WLABEL  (v) = 'synthetic tracer'//is
       elseif( NQC_STR <= v .AND. v <= NQC_END ) then ! Yamashita [Passive tracer]
          write(is,'(I2.2)') v - NQC_STR + 1
          TRC_name(v) = 'passive'//is
          WLABEL  (v) = 'passive tracer'//is
          TRC_unit(v) = 'kg/kg'
       elseif( RDN_STR <= v .AND. v <= RDN_END ) then ! Yamashita [Passive tracer]
          write(is1,'(I1.1)') v - RDN_STR + 1
          TRC_name(v) = 'radon'//is1
          WLABEL  (v) = 'radon'//is1
          TRC_unit(v) = 'Bq/kg'
       elseif( NCHEM_STR <= v .AND. v <= NCHEM_END ) then ! Yamashita [CHASER] 17/09/28
          TRC_name(v) = CH_TRC_name(v-NCHEM_STR+1)
          WLABEL  (v) = CH_TRC_desc(v-NCHEM_STR+1)
          TRC_unit(v) = 'kg/kg'
       endif
    enddo

    if ( TM_RUN == 'ON' ) then
       OUT_FILE_TYPE = 'obs' ! Force OUT_FILE_TYPE='obs'
       do v = TM_STR, TM_END
          TRC_type(v) = 'TM'
       enddo
    endif

    allocate( flag_diagnose_number(TRC_vmax) )
    flag_diagnose_number(:) = .false.

    PRG_vmax   = PRG_vmax0  + TRC_vmax
    I_RHOGQend = PRG_vmax

    DIAG_vmax  = DIAG_vmax0 + TRC_vmax
    I_qend     = DIAG_vmax

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*)                            '===== Prognostic Tracers ====='
    write(IO_FID_LOG,'(1x,A,I3)')                  '--- Number of Tracer Variables    = ', TRC_vmax
    do v = 1, TRC_vmax
       write(IO_FID_LOG,'(1x,A,I0.3,A,A16,A,A32)') '--- [', v, '] = ', trim(TRC_name(v)), ' : ', trim(WLABEL(v))
    enddo
    write(IO_FID_LOG,'(1x,A,I3,A,I0.3,A,I0.3,A)')  '--- Thermodynamic (water) Tracers = ', NQW_MAX, ' [', NQW_STR, '-', NQW_END, ']'

    return
  end subroutine RUNCONF_tracer_setup

  !-----------------------------------------------------------------------------
  ! Niwa [TM]
  subroutine runconf_set_trc_name( idx, name, label )
    implicit none

    integer,          intent(in) :: idx
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: label
    !---------------------------------------------------------------------------

    TRC_name(idx) = trim(name)
    WLABEL  (idx) = trim(label)
    !TRC_unit(idx) = trim(unit)

    write(IO_FID_LOG,*) '*** Labeling of Extra Tracers'
    write(IO_FID_LOG,'(1x,A,I4,A,A16,A,A,A)') '|ID=', idx, ':', TRC_name(idx), ':', WLABEL(idx),'|'

    return
  end subroutine runconf_set_trc_name

  !-----------------------------------------------------------------------------
  !> thermodynamic setup
  subroutine RUNCONF_thermodyn_setup
    use mod_cnst, only: &
       CPdry => CNST_CPdry, &
       CPvap => CNST_CPvap, &
       CVdry => CNST_CVdry, &
       CVvap => CNST_CVvap, &
       CL    => CNST_CL,    &
       CI    => CNST_CI,    &
       CNST_LHV,            &
       CNST_LHS,            &
       CNST_LHF,            &
       CNST_LHV00,          &
       CNST_LHS00,          &
       CNST_LHF00,          &
       CNST_LHV0,           &
       CNST_LHS0,           &
       CNST_LHF0
    implicit none

    integer  :: v
    !---------------------------------------------------------------------------
    ! 'SIMPLE': standard approximation CVD * T
    ! 'EXACT': exact formulation
    !         -> if warm rain
    !            qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
    !         -> if cold rain
    !            qd*CVD*T + qv*CVV*T + (qc+qr)*CPL*T
    !            + (qi+qs)*CPI*T

    !--- Heat capacity for thermodynamics
    allocate( CVW(NQW_STR:NQW_END) )
    allocate( CPW(NQW_STR:NQW_END) )

    !--- Heat capacity for SBM tracers [ADD] Y.Sato 20171227
    allocate(CVW_bin(BIN_STR:BIN_END))
    allocate(CPW_bin(BIN_STR:BIN_END))

    if ( EIN_TYPE == 'SIMPLE' ) then

       LHV = CNST_LHV0
       LHS = CNST_LHS0
       LHF = CNST_LHF0

       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVdry
             CPW(v) = CPdry
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QR ) then ! rain
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QI ) then ! ice
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QS ) then ! snow
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CVdry
             CPW(v) = CVdry
          elseif( v == I_QH ) then ! hail
             CVW(v) = CVdry
             CPW(v) = CVdry
          endif
       enddo

       if ( MP_TYPE == 'S10' ) then ! Y.Sato [SBM]
          if ( RAIN_TYPE == 'WARM' ) then
             do v = BIN_STR, BIN_END ! liquid water
                CVW_bin(v) = CVdry
                CPW_bin(v) = CVdry
             enddo
          elseif( RAIN_TYPE == 'COLD' ) then
             do v = BIN_STR, BIN_STR+NUM_BIN-1 ! liquid water
                CVW_bin(v) = CVdry
                CPW_bin(v) = CVdry
             enddo
             do v = BIN_STR+NUM_BIN, BIN_END ! solid water
                CVW_bin(v) = CVdry
                CPW_bin(v) = CVdry
             enddo
          endif
       endif

    elseif( EIN_TYPE == 'SIMPLE2' ) then

       LHV = CNST_LHV0
       LHS = CNST_LHS0
       LHF = CNST_LHF0

       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVvap
             CPW(v) = CPvap
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QR ) then ! rain
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QI ) then ! ice
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QS ) then ! snow
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CPvap
             CPW(v) = CPvap
          elseif( v == I_QH ) then ! hail
             CVW(v) = CPvap
             CPW(v) = CPvap
          endif
       enddo

       if ( MP_TYPE == 'S10' ) then ! Y.Sato [SBM]
          if ( RAIN_TYPE == 'WARM' ) then
             do v = BIN_STR, BIN_END ! liquid water
                CVW_bin(v) = CPvap
                CPW_bin(v) = CPvap
             enddo
          elseif( RAIN_TYPE == 'COLD' ) then
             do v = BIN_STR, BIN_STR+NUM_BIN-1 ! liquid water
                CVW_bin(v) = CPvap
                CPW_bin(v) = CPvap
             enddo
             do v = BIN_STR+NUM_BIN, BIN_END ! solid water
                CVW_bin(v) = CPvap
                CPW_bin(v) = CPvap
             enddo
          endif
       endif

    elseif( EIN_TYPE == 'EXACT' ) then

       LHV = CNST_LHV00
       LHS = CNST_LHS00
       LHF = CNST_LHF00

       do v = NQW_STR, NQW_END
          if    ( v == I_QV ) then ! vapor
             CVW(v) = CVvap
             CPW(v) = CPvap
          elseif( v == I_QC ) then ! cloud
             CVW(v) = CL
             CPW(v) = CL
          elseif( v == I_QR ) then ! rain
             CVW(v) = CL
             CPW(v) = CL
          elseif( v == I_QI ) then ! ice
             CVW(v) = CI
             CPW(v) = CI
          elseif( v == I_QS ) then ! snow
             CVW(v) = CI
             CPW(v) = CI
          elseif( v == I_QG ) then ! graupel
             CVW(v) = CI
             CPW(v) = CI
          elseif( v == I_QH ) then ! hail
             CVW(v) = CI
             CPW(v) = CI
          endif
       enddo

       if ( MP_TYPE == 'S10' ) then ! Y.Sato [SBM]
          if ( RAIN_TYPE == 'WARM' ) then
             do v = BIN_STR, BIN_END ! liquid water
                CVW_bin(v) = CL
                CPW_bin(v) = CL
             enddo
          elseif( RAIN_TYPE == 'COLD' ) then
             do v = BIN_STR, BIN_STR+NUM_BIN-1 ! liquid water
                CVW_bin(v) = CL
                CPW_bin(v) = CL
             enddo
             do v = BIN_STR+NUM_BIN, BIN_END ! solid water
                CVW_bin(v) = CI
                CPW_bin(v) = CI
             enddo
          endif
       endif

    endif

    CNST_LHV = LHV
    CNST_LHS = LHS
    CNST_LHF = LHF

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*)              '===== Thermodynamic parameters ====='
    write(IO_FID_LOG,*)              '--- Internal Energy Calc. Type : ', trim(EIN_TYPE)
    write(IO_FID_LOG,'(1x,A,F13.3)') '--- Latent Heat of vaporizaion : ', LHV
    write(IO_FID_LOG,'(1x,A,F13.3)') '--- Latent Heat of sublimation : ', LHS
    write(IO_FID_LOG,'(1x,A,F13.3)') '--- Latent Heat of fusion      : ', LHF

    write(IO_FID_LOG,*)              '--- Specific Heat Capacity at constant volume'
    do v = NQW_STR, NQW_END
       write(IO_FID_LOG,'(1x,A,A16,A,F10.3)') '--- ', trim(TRC_name(v)), ' = ', CVW(v)
    enddo

    write(IO_FID_LOG,*)              '--- Specific Heat Capacity at constant pressure'
    do v = NQW_STR, NQW_END
       write(IO_FID_LOG,'(1x,A,A16,A,F10.3)') '--- ', trim(TRC_name(v)), ' = ', CPW(v)
    enddo

    return
  end subroutine RUNCONF_thermodyn_setup

end module mod_runconf
