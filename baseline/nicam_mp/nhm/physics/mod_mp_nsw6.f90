!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          6-category one-moment bulk cloud microphysics NSW6
!!          Reference: Tomita(2008)
!!
!!          Reference : Berry(1968)                : Proc.Phys.Soc.,B66,pp.688-694
!!                      Hsie et al.(1980)          : J.Appl.Meteor.,19,pp.950-977
!!                      Grabowski(1998)            : J.Atmos.Sci.,55,pp.3283-3298
!!                      Hon and Lim(2006)          : J.Korean Meteor.Soc.,42,pp.129-151
!!                      Lin et al.(1983)           : J.Appl.Meteor.,22,pp.1065-1092
!!                      Prupacher and Klett(1978)  : Micophysics of Clouds and
!!                                                   Precipitation.
!!                                                   Kluwer Academic Publishers
!!                      Ruttledge and Hobbs(1983)  : J.Atmos.Sci.,40,pp.1185-1206
!!                      Ruttledge and Hobbs(1984)  : J.Atmos.Sci.,40,pp.2949-2977
!!                      Heymsfield and Doner(1990) : J.Atmos.Sci.,??,pp.1865-1877
!!                      Khairoutdinov and Kogan (2000)
!!                      Roh and Satoh(2014)        : J.Atmos.Sci.,71,pp.2654-2673
!!                      Aerosol Indirect Effect:
!!                      Suzuki et al. (2004)       : J.Atmos.Sci.,61,pp.179-194
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_mp_nsw6
  !-----------------------------------------------------------------------------
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !         0.00   06/12/30   add this module
  !         0.00   07/02/07   M.Satoh, rm multiple use of mod_adm in mp_nsw6
  !         0.01   07/07/24   K.Suzuki implementing SPRINTARS
  !                           1. Add OPT_INDIR for aerosol indirect effect
  !                           2. Add UNCCN to argument for indirect effect
  !                           3. Add precip_flux to argument for use in SPRINTARS
  !                07/08/21   H. Tomita : dconc, ddisp -> dconc_def, ddisp_def
  !                                       for vectorization on ES.
  !                07/11/13   H.Tomita : Add an option INC_PGAUT.
  !                07/11/22   T.Mitsui: Effective Radius is calculated in rd_driver => mp
  !                                    add calculation of rceff, rctop, rwtop, t
  !                07/11/27   H.Tomita : Bug fix (suggested by Mr.Matsuda)
  !                                      Change default
  !                                      1. INC_PGAUT = .TRUE. -> .FALSE.
  !                                      2. PRECIP_TRANSPORT_TYPE = '4WATER'
  !                                         ->  '3WATER'
  !                07/11/30   T.Mitsui: vectorization on ES.
  !                07/12/21   H.Tomita: Bug fix of Pgaut
  !                08/04/22   H.Tomita: Add option OPT_EXPLICIT_ICEGEN
  !                09/04/14   T.Mitsui: add procedures mp_nsw6_terminal_velocity
  !                                                and mp_nsw6_effective_radius
  !                09/04/15   T.Mitsui: avoid negative mass in mp_nsw6_effective_radius
  !                09-07-10   H.Tomita : precip_transport is used.
  !                09-09-1X   S.Iga : Merge '09/01/19 [Add] M.Satoh' (output conversion rate)
  !                09-10-13   S.Iga : Commented out the part 'OPT_EXPLICIT_ICEGEN'
  !                                   for vectorization.
  !                09-11-15   S.Iga : Bug Fix (my bug at 09-10-13)
  !                10-05-11   M.Satoh : add l_region in history_in
  !                10-10-07   A.T.Noda : CONST_v_term_qi and Psaut_beta0
  !                11-08-01   A.T.Noda : avoid an argument of exponential to be 0
  !                11-08-03   A.T.Noda : use 'cbrt' instead of SMALL
  !                11-08-06   A.T.Noda : add SMALL2 and SMALL3, which are used for squart root
  !                                      and cubic root computations, respectively
  !                12-03-12   teraim : Use preprocessor for OPT_EXPICIT_ICEGEN, OPT_INDIR and OPT_OUTPUT_XXX.
  !                                    Use array to reduce register spill.
  !                12-05-18   H.Yashiro : Speedup Flux-Semi-Lagrangian & Simple Upstream-Euler
  !                12-06-06   H.Yashiro : bugfix
  !                12-08-01   T.Miyakawa : enabled specification of cs and cr ('Miura parameter')
  !                                        from namelist
  !                13-04-03   H.Yashiro : bugfix
  !                14-10-11   W. Roh    : add flag Roh and Satoh 2014
  !                15-02-26   T.Seiki    : add KK(2000) scheme
  !                15-10-21   H.Yashiro : bugfix
  !                16-12-05   H.Yashiro : Full refactoring
  !                                       1)
  !                18-02-08   T.Seiki   : supercooled liquid cloud ver. (implimented by A.Noda on 20-10-22)
  !                20-08-19   A.Noda    : Add parametarization of b_aut based on NDW6 in GL10
  !                                       (OPT_beta_saut), suggested by T.Seiki
  !                20-11-17   T.Seiki   : bug-Fix
  !                21-10-08   M.Tanoue  : Add stable water isotopes
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision, only: &
     SP, DP, RP => RP_MP, RP_G => RP
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: mp_nsw6_init
  public :: mp_nsw6
  public :: mp_nsw6_terminal_velocity
  public :: mp_nsw6_effective_radius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: Bergeron_param

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: dens00 = 1.28_RP   !< standard density [kg/m3]

  ! Parameter for Marshall-Palmer distribution
  real(RP), private            :: N0r    = 8.E+6_RP  !< intercept parameter for rain    [1/m4]
  real(RP), private            :: N0s    = 3.E+6_RP  !< intercept parameter for snow    [1/m4]
  real(RP), private            :: N0g    = 4.E+6_RP  !< intercept parameter for graupel [1/m4]

  real(RP), private            :: rho_s  = 100.0_RP  !< density of snow    [kg/m3]
  real(RP), private            :: rho_g  = 400.0_RP  !< density of graupel [kg/m3]
                                                     !   graupel : 400
                                                     !   hail    : 917

  real(RP), private            :: C_d    =   0.6_RP  !< drag coefficient for graupel
  real(RP), private            :: Cr     = 130.0_RP
  real(RP), private            :: Cs     =  4.84_RP

  ! Empirical parameter
  real(RP), private            :: Ar, As, Ag
  real(RP), private            :: Br, Bs, Bg
  real(RP), private            :: Cg
  real(RP), private            :: Dr, Ds, Dg

  ! GAMMA function
  real(RP), private            :: GAM, GAM_2, GAM_3

  real(RP), private            :: GAM_1br, GAM_2br, GAM_3br
  real(RP), private            :: GAM_3dr
  real(RP), private            :: GAM_6dr
  real(RP), private            :: GAM_1brdr
  real(RP), private            :: GAM_5dr_h

  real(RP), private            :: GAM_1bs, GAM_2bs, GAM_3bs
  real(RP), private            :: GAM_3ds
  real(RP), private            :: GAM_1bsds
  real(RP), private            :: GAM_5ds_h

  real(RP), private            :: GAM_1bg, GAM_3dg
  real(RP), private            :: GAM_1bgdg
  real(RP), private            :: GAM_5dg_h

  !---< Khairoutdinov and Kogan (2000) >---
  real(RP), private            :: sw_kk2000       = 0.0_RP     !< switch for k-k scheme

  !---< Roh and Satoh (2014) >---
  logical,  private            :: Roh_flag        = .true.     !< enable setting by Roh and Satoh (2014)?
  real(RP), private            :: sw_roh2014      = 0.0_RP     !< switch for Roh scheme
  real(RP), private            :: ln10                         !< log(10)
  real(RP), private            :: coef_a(10) = (/ 5.065339_RP, -0.062659_RP, -3.032362_RP, 0.029469_RP, -0.000285_RP, &
                                                  0.31255_RP,   0.000204_RP,  0.003199_RP, 0.0_RP,      -0.015952_RP  /)
  real(RP), private            :: coef_b(10) = (/ 0.476221_RP, -0.015896_RP,  0.165977_RP, 0.007468_RP, -0.000141_RP, &
                                                  0.060366_RP,  0.000079_RP,  0.000594_RP, 0.0_RP,      -0.003577_RP  /)

  !---< Seiki and Roh (2020) >---
  logical,  private            :: OPT_SCLW        = .false.    !< option to represent Super-Cooled Liquid Water
  real(RP), private            :: CRIT_SCLW_W     = 0.01_RP    !< critical value to compute SCLW 2021/3/12 A.Noda
  real(RP), private            :: Ds_cutoff       = 150.E-6_RP !< cutoff diameter [m] for riming process
  real(RP), private            :: Dg_cutoff       = 150.E-6_RP !< cutoff diameter [m] for riming process

  real(RP), private            :: sw_sr2020       = 0.0_RP     !< switch for SCLW scheme
  real(RP), private            :: sw_phillips2007 = 0.0_RP     !< heterogenous ice nucleation scheme Phillips et al. (2007), JAS
  real(RP), private            :: sw_beta_saut    = 0.0_RP

  ! Accretion parameter
  real(RP), private            :: Eiw        = 1.0_RP      !< collection efficiency of cloud ice for cloud water
  real(RP), private            :: Erw        = 1.0_RP      !< collection efficiency of rain    for cloud water
  real(RP), private            :: Esw        = 1.0_RP      !< collection efficiency of snow    for cloud water
  real(RP), private            :: Egw        = 1.0_RP      !< collection efficiency of graupel for cloud water
  real(RP), private            :: Eri        = 1.0_RP      !< collection efficiency of rain    for cloud ice
  real(RP), private            :: Esi        = 1.0_RP      !< collection efficiency of snow    for cloud ice
  real(RP), private            :: Egi        = 0.1_RP      !< collection efficiency of graupel for cloud ice
  real(RP), private            :: Esr        = 1.0_RP      !< collection efficiency of snow    for rain
  real(RP), private            :: Egr        = 1.0_RP      !< collection efficiency of graupel for rain
  real(RP), private            :: Egs        = 1.0_RP      !< collection efficiency of graupel for snow
  real(RP), private            :: gamma_sacr = 25.E-3_RP   !< effect of low temperature for Esi
  real(RP), private            :: gamma_gacs = 90.E-3_RP   !< effect of low temperature for Egs
  real(RP), private            :: mi         = 4.19E-13_RP !< mass of one cloud ice crystal [kg]

  ! Auto-conversion parameter
  real(RP), private, parameter :: Nc_lnd     = 2000.0_RP   !< number concentration of cloud water (land)  [1/cc]
  real(RP), private, parameter :: Nc_ocn     =   50.0_RP   !< number concentration of cloud water (ocean) [1/cc]
  real(RP), private            :: Nc_def                   !< number concentration of cloud water         [1/cc]

  real(RP), private            :: beta_saut  =  1.E-3_RP   !< auto-conversion factor beta  for ice
  real(RP), private            :: gamma_saut = 25.E-3_RP   !< auto-conversion factor gamma for ice
  real(RP), private            :: beta_gaut  =  1.E-3_RP   !< auto-conversion factor beta  for snow
  real(RP), private            :: gamma_gaut = 90.E-3_RP   !< auto-conversion factor gamma for snow
  real(RP), private            :: qicrt_saut =  0.0_RP     !< mixing ratio threshold for Psaut [kg/kg]
  real(RP), private            :: qscrt_gaut =  6.E-4_RP   !< mixing ratio threshold for Pgaut [kg/kg]

  ! Evaporation, Sublimation parameter
  real(RP), private, parameter :: Ka0        = 2.428E-2_RP !< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
  real(RP), private, parameter :: dKa_dT     =  7.47E-5_RP !< Coefficient of Ka depending on temperature      [J/m/s/K/K]
  real(RP), private, parameter :: Kd0        = 2.222E-5_RP !< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
  real(RP), private, parameter :: dKd_dT     =  1.37E-7_RP !< Coefficient of Dw depending on temperature                 [m2/s/K]
  real(RP), private, parameter :: nu0        = 1.718E-5_RP !< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
  real(RP), private, parameter :: dnu_dT     =  5.28E-8_RP !< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

  real(RP), private            :: f1r        = 0.78_RP     !< ventilation factor 1 for rain
  real(RP), private            :: f2r        = 0.27_RP     !< ventilation factor 2 for rain
  real(RP), private            :: f1s        = 0.65_RP     !< ventilation factor 1 for snow
  real(RP), private            :: f2s        = 0.39_RP     !< ventilation factor 2 for snow
  real(RP), private            :: f1g        = 0.78_RP     !< ventilation factor 1 for graupel
  real(RP), private            :: f2g        = 0.27_RP     !< ventilation factor 2 for graupel

  ! Freezing parameter
  real(RP), private            :: A_frz      = 0.66_RP     !< freezing factor [/K]
  real(RP), private            :: B_frz      = 100.0_RP    !< freezing factor [/m3/s]

  ! Bergeron process parameter
  real(RP), private            :: mi40       = 2.46E-10_RP !< mass              of a 40 micron ice crystal [kg]
  real(RP), private            :: mi50       = 4.80E-10_RP !< mass              of a 50 micron ice crystal [kg]
  real(RP), private            :: vti50      = 1.0_RP      !< terminal velocity of a 50 micron ice crystal [m/s]
  real(RP), private            :: Ri50       = 5.E-5_RP    !< radius            of a 50 micron ice crystal [m]

  ! Explicit ice generation
  logical,  private            :: OPT_EXPLICIT_ICEGEN = .false. !< enable explicit ice generation?
  real(RP), private            :: sw_expice  = 0.0_RP      !< switch for explicit ice generation
  !real(RP), private, parameter :: Nc_ihtr    = 300.0_RP    !< cloud number concentration for heterogeneous ice nucleation [1/cc] [Del] 2020-11-17 T.Seiki, delete to avoid inconsistency
  real(RP), private, parameter :: Di_max     = 500.E-6_RP
  real(RP), private, parameter :: Di_a       = 11.9_RP

  !  aerosol indirect effect
  logical, private             :: OPT_INDIR  = .false.   ! enable aerosol indirect effect?

  integer,  private, parameter :: wk_nmax   = 49
  integer,  private, parameter :: I_dqv_dt  =  1 !
  integer,  private, parameter :: I_dqc_dt  =  2 !
  integer,  private, parameter :: I_dqr_dt  =  3 !
  integer,  private, parameter :: I_dqi_dt  =  4 !
  integer,  private, parameter :: I_dqs_dt  =  5 !
  integer,  private, parameter :: I_dqg_dt  =  6 !
  integer,  private, parameter :: I_delta1  =  7 ! separation switch for r->s,g
  integer,  private, parameter :: I_delta2  =  8 ! separation switch for s->g
  integer,  private, parameter :: I_spsati  =  9 ! separation switch for ice sublimation
  integer,  private, parameter :: I_iceflg  = 10 ! separation switch for T > 0
  integer,  private, parameter :: I_RLMDr   = 11
  integer,  private, parameter :: I_RLMDs   = 12
  integer,  private, parameter :: I_RLMDg   = 13
  integer,  private, parameter :: I_Piacr   = 14 ! r->s,g
  integer,  private, parameter :: I_Psacr   = 15 ! r->s,g
  integer,  private, parameter :: I_Praci   = 16 ! i->s,g
  integer,  private, parameter :: I_Pigen   = 17 ! v->i
  integer,  private, parameter :: I_Pidep   = 18 ! v->i
  integer,  private, parameter :: I_Psdep   = 19 ! v->s
  integer,  private, parameter :: I_Pgdep   = 20 ! v->g
  integer,  private, parameter :: I_Praut   = 21 ! c->r
  integer,  private, parameter :: I_Pracw   = 22 ! c->r
  integer,  private, parameter :: I_Pihom   = 23 ! c->i
  integer,  private, parameter :: I_Pihtr   = 24 ! c->i
  integer,  private, parameter :: I_Psacw   = 25 ! c->s
  integer,  private, parameter :: I_Psfw    = 26 ! c->s
  integer,  private, parameter :: I_Pgacw   = 27 ! c->g
  integer,  private, parameter :: I_Prevp   = 28 ! r->v
  integer,  private, parameter :: I_Piacr_s = 29 ! r->s
  integer,  private, parameter :: I_Psacr_s = 30 ! r->s
  integer,  private, parameter :: I_Piacr_g = 31 ! r->g
  integer,  private, parameter :: I_Psacr_g = 32 ! r->g
  integer,  private, parameter :: I_Pgacr   = 33 ! r->g
  integer,  private, parameter :: I_Pgfrz   = 34 ! r->g
  integer,  private, parameter :: I_Pisub   = 35 ! i->v
  integer,  private, parameter :: I_Pimlt   = 36 ! i->c
  integer,  private, parameter :: I_Psaut   = 37 ! i->s
  integer,  private, parameter :: I_Praci_s = 38 ! i->s
  integer,  private, parameter :: I_Psaci   = 39 ! i->s
  integer,  private, parameter :: I_Psfi    = 40 ! i->s
  integer,  private, parameter :: I_Praci_g = 41 ! i->g
  integer,  private, parameter :: I_Pgaci   = 42 ! i->g
  integer,  private, parameter :: I_Pssub   = 43 ! s->v
  integer,  private, parameter :: I_Psmlt   = 44 ! s->r
  integer,  private, parameter :: I_Pgaut   = 45 ! s->g
  integer,  private, parameter :: I_Pracs   = 46 ! s->g
  integer,  private, parameter :: I_Pgacs   = 47 ! s->g
  integer,  private, parameter :: I_Pgsub   = 48 ! g->v
  integer,  private, parameter :: I_Pgmlt   = 49 ! g->r

  character(len=H_SHORT), private :: w_name(wk_nmax)

  data w_name / 'dqv_dt ', &
                'dqc_dt ', &
                'dqr_dt ', &
                'dqi_dt ', &
                'dqs_dt ', &
                'dqg_dt ', &
                'delta1 ', &
                'delta2 ', &
                'spsati ', &
                'iceflg ', &
                'RLMDr  ', &
                'RLMDs  ', &
                'RLMDg  ', &
                'Piacr  ', &
                'Psacr  ', &
                'Praci  ', &
                'Pigen  ', &
                'Pidep  ', &
                'Psdep  ', &
                'Pgdep  ', &
                'Praut  ', &
                'Pracw  ', &
                'Pihom  ', &
                'Pihtr  ', &
                'Psacw  ', &
                'Psfw   ', &
                'Pgacw  ', &
                'Prevp  ', &
                'Piacr_s', &
                'Psacr_s', &
                'Piacr_g', &
                'Psacr_g', &
                'Pgacr  ', &
                'Pgfrz  ', &
                'Pisub  ', &
                'Pimlt  ', &
                'Psaut  ', &
                'Praci_s', &
                'Psaci  ', &
                'Psfi   ', &
                'Praci_g', &
                'Pgaci  ', &
                'Pssub  ', &
                'Psmlt  ', &
                'Pgaut  ', &
                'Pracs  ', &
                'Pgacs  ', &
                'Pgsub  ', &
                'Pgmlt  '  /

  character(len=H_SHORT), private :: precip_transport_type = '4WATER'
  character(len=H_SHORT), private :: precip_scheme_type    = 'Flux-Semilag_new'
                                                           ! 'Default' = SL_UPSTREAM
                                                           ! 'Flux-Semilag_new', same as default but fast
                                                           ! 'Upwind-Euler'
  real(RP),               private :: constVti              ! constant terminal velocity for ice [m/s], upward is positive
  real(RP),               private :: sw_constVti           = 0.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine mp_nsw6_init
    use mod_misc, only: &
       MISC_gammafunc
    use mod_adm, only: &
       ADM_proc_stop
    use mod_cnst, only: &
       UNDEF => CNST_UNDEF, &
       PI    => CNST_PI,    &
       GRAV  => CNST_GRAV,  &
       rho_w => CNST_DWATR
    implicit none

    logical                :: INC_PGAUT       = .false.
    real(RP)               :: PSAUT_BETA0
    real(RP)               :: CNST_V_TERM_QI              !< force constant terminal velocity for ice
    character(len=H_SHORT) :: qr_aut_acc_type = 'Default' !<
                                              ! 'Default' !< Berry-type scheme
                                              ! 'KK2000'  !< Khairoutdinov and Kogan (2000)
    logical                :: opt_beta_saut   = .false.   !< use parameterization for beta_saut
    logical                :: opt_vt_ndw6     = .false.   !< change terminal velocity coefficient

    namelist / nm_mp_nsw6 / &
         N0r,                   &
         N0s,                   &
         N0g,                   &
         rho_g,                 &
         rho_s,                 &
         C_d,                   &
         Cr,                    &
         Cs,                    &
         Eiw,                   &
         Erw,                   &
         Esw,                   &
         Egw,                   &
         Eri,                   &
         Esi,                   &
         Egi,                   &
         Esr,                   &
         Egr,                   &
         Egs,                   &
         gamma_sacr,            &
         gamma_gacs,            &
         mi,                    &
         Nc_def,                &
         PSAUT_BETA0,           &
         gamma_saut,            &
         qicrt_saut,            &
         beta_gaut,             &
         gamma_gaut,            &
         qscrt_gaut,            &
         f1r,                   &
         f2r,                   &
         f1s,                   &
         f2s,                   &
         f1g,                   &
         f2g,                   &
         A_frz,                 &
         B_frz,                 &
         precip_scheme_type,    &
         precip_transport_type, &
         INC_PGAUT,             &
         OPT_EXPLICIT_ICEGEN,   &
         OPT_INDIR,             &
         CNST_V_TERM_QI,        &
         qr_aut_acc_type,       & ! [Add] T.Seiki 2015/02/26
         Roh_flag,              & ! [add] WS. Roh 20141011
         OPT_SCLW,              & ! [Add] T.Seiki 2018/02/08
         CRIT_SCLW_W,           & ! [Add] A.Noda  2021/03/12
         opt_beta_saut,         & ! [Add] A.Noda  2020/08/19
         Ds_cutoff,             & ! [Add] T.Seiki 2018/02/08
         Dg_cutoff,             & ! [Add] T.Seiki 2018/02/08
         opt_vt_ndw6              ! [Add] T.Seiki 2018/02/28

    integer :: ierr
    !---------------------------------------------------------------------------

    Nc_def         = Nc_ocn
    CNST_V_TERM_QI = UNDEF
    PSAUT_BETA0    = beta_saut

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[MP NSW6]/Category[nhm physics]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_nsw6,iostat=ierr)
    if (ierr<0) then
       write(IO_FID_LOG,*) '*** nm_mp_nsw6 is not specified. use default.'
    elseif(ierr>0) then
       write(*,         *) 'xxx Not appropriate names in namelist nm_mp_nsw6. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist nm_mp_nsw6. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=nm_mp_nsw6)

    beta_saut = PSAUT_BETA0

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Calculation flag of sedimentation:'
    write(IO_FID_LOG,*) '*** QV => NO'
    write(IO_FID_LOG,*) '*** QC => NO'
    write(IO_FID_LOG,*) '*** QR => YES'
    if    ( precip_transport_type == '3WATER' ) then
       write(IO_FID_LOG,*) '*** QI => NO'
    elseif( precip_transport_type == '4WATER' ) then
       write(IO_FID_LOG,*) '*** QI => YES'
    endif
    write(IO_FID_LOG,*) '*** QS => YES'
    write(IO_FID_LOG,*) '*** QG => YES'

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Precipitation(sedimentation) scheme:'
    if    ( precip_scheme_type == 'Upwind-Euler' ) then
       write(IO_FID_LOG,*) '*** => Upwind-Euler'
    elseif( precip_scheme_type == 'Flux-Semilag_new' ) then
       write(IO_FID_LOG,*) '*** => Flux-Semilag_new'
    else
       write(IO_FID_LOG,*) '*** => Default(Flux-Semilag)'
    endif

    !--- empirical coefficients A, B, C, D
    Ar = PI * rho_w / 6.0_RP
    As = PI * rho_s / 6.0_RP
    Ag = PI * rho_g / 6.0_RP

    Br = 3.0_RP
    Bs = 3.0_RP
    Bg = 3.0_RP

    Cg = sqrt( ( 4.0_RP * rho_g * GRAV ) / ( 3.0_RP * dens00 * C_d ) )

    Dr = 0.50_RP
    Ds = 0.25_RP
    Dg = 0.50_RP

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Autoconversion & Accretion scheme for QC->Qr:'
    if    ( qr_aut_acc_type == 'Default' ) then
       write(IO_FID_LOG,*) '*** => Berry(1968) : default'
       sw_kk2000 = 0.0_RP
    elseif( qr_aut_acc_type == 'KK2000' ) then
       write(IO_FID_LOG,*) '*** => Khairoutdinov and Kogan(2000)'
       sw_kk2000 = 1.0_RP
    else
       write(*,         *) 'xxx Not appropriate qr_aut_acc_type. STOP', trim(qr_aut_acc_type)
       write(IO_FID_LOG,*) 'xxx Not appropriate qr_aut_acc_type. STOP', trim(qr_aut_acc_type)
       call ADM_proc_stop
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use setting of Roh and Satoh (2014)?:'
    if ( Roh_flag ) then ! overwrite parameters
       write(IO_FID_LOG,*) '*** => Yes'
       OPT_EXPLICIT_ICEGEN = .true.

       sw_roh2014 = 1.0_RP
       N0g        = 4.E+8_RP
       As         = 0.069_RP
       Bs         = 2.0_RP
       Esi        = 0.25_RP
       Egi        = 0.0_RP
       Egs        = 0.0_RP
    else
       write(IO_FID_LOG,*) '*** => No : default'
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use setting of Seiki and Roh (2020)?:'
    if ( OPT_SCLW ) then ! [Add] T.Seiki 2018/02/28
       write(IO_FID_LOG,*) '*** => Yes'
       sw_sr2020       = 1.0_RP ! opt_pgsub_off=T & opt_satadj_qc=F, opt_double_count is always F
       sw_phillips2007 = 1.0_RP
       !Ds_cutoff = 150.E-6_RP : default value, changeable with namelist
       !Dg_cutoff = 150.E-6_RP : default value, changeable with namelist
    else
       write(IO_FID_LOG,*) '*** => No : default'
       Ds_cutoff = 0.0_RP ! force disable => opt_cutoff_riming=F
       Dg_cutoff = 0.0_RP ! force disable => opt_cutoff_riming=F
    endif

    if ( opt_beta_saut ) then
       sw_beta_saut = 1.0_RP
    else
       sw_beta_saut = 0.0_RP
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use explicit ice generation scheme?:'
    if ( OPT_EXPLICIT_ICEGEN ) then
       write(IO_FID_LOG,*) '*** => Yes'
       sw_expice = 1.0_RP
    else
       write(IO_FID_LOG,*) '*** => No : default'
       sw_expice = 0.0_RP
    endif

    if ( .NOT. INC_PGAUT ) then
       beta_gaut = 0.0_RP
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use constant terminal velocity for ice?'
    if ( CNST_V_TERM_QI /= UNDEF ) then
       sw_constVti = 1.0_RP
       constVti    = -abs(CNST_V_TERM_QI)
    else
       sw_constVti = 0.0_RP
       constVti    = 0.0_RP ! for debug mode to avoid undefined
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use setting of NDW6 terminal velocities?'
    if ( opt_vt_ndw6 ) then ! [Add] T.Seiki 2018/02/28
       write(IO_FID_LOG,*) '*** => Yes'

       Cs = 257.5_RP
       Cg = 831.3_RP
       Ds = 0.808_RP
       Dg = 0.925_RP
    else
       write(IO_FID_LOG,*) '*** => No : default'
    endif

    GAM       = 1.0_RP ! =0!
    GAM_2     = 1.0_RP ! =1!
    GAM_3     = 2.0_RP ! =2!

    GAM_1br   = MISC_gammafunc( 1.0_RP + Br ) ! = 4!
    GAM_2br   = MISC_gammafunc( 2.0_RP + Br ) ! = 5!
    GAM_3br   = MISC_gammafunc( 3.0_RP + Br ) ! = 6!
    GAM_3dr   = MISC_gammafunc( 3.0_RP + Dr )
    GAM_6dr   = MISC_gammafunc( 6.0_RP + Dr )
    GAM_1brdr = MISC_gammafunc( 1.0_RP + Br + Dr )
    GAM_5dr_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Dr) )

    GAM_1bs   = MISC_gammafunc( 1.0_RP + Bs ) ! = 4!
    GAM_2bs   = MISC_gammafunc( 2.0_RP + Bs ) ! = 5!
    GAM_3bs   = MISC_gammafunc( 3.0_RP + Bs ) ! = 6!
    GAM_3ds   = MISC_gammafunc( 3.0_RP + Ds )
    GAM_1bsds = MISC_gammafunc( 1.0_RP + Bs + Ds )
    GAM_5ds_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Ds) )

    GAM_1bg   = MISC_gammafunc( 1.0_RP + Bg ) ! = 4!
    GAM_3dg   = MISC_gammafunc( 3.0_RP + Dg )
    GAM_1bgdg = MISC_gammafunc( 1.0_RP + Bg + Dg)
    GAM_5dg_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Dg) )

    ln10 = log(10.0_RP)

    return
  end subroutine mp_nsw6_init

  !-----------------------------------------------------------------------------
  subroutine mp_nsw6( &
       ijdim,          &
       l_region,       &
       rhog,           &
       rhogvx,         &
       rhogvy,         &
       rhogvz,         &
       rhogw,          &
       rhoge,          &
       rhogq,          &
       vx,             &
       vy,             &
       vz,             &
       w,              &
       UNCCN,          &
       factor_n_aie,   & ! [Add] 2017/06/13 T.Seiki
       rho,            &
       tem,            &
       pre,            &
       q,              &
       qd,             &
       precip,         &
       precip_ISO1,    & ! [Add] M.Tanoue 20200706
       precip_ISO2,    & ! [Add] M.Tanoue 20200706
       precip_rhoe,    &
       precip_lh_heat, &
       precip_rhophi,  &
       precip_rhokin,  &
       precip_flux,    &
       evap_liq,       & ! CHASER: 2018/11/08 Yamashita
       evap_ice,       & ! CHASER: 2018/11/08 Yamashita
       rceff,          &
       rctop,          &
       rwtop,          &
       tctop,          &
       rceff_cld,      &
       rctop_cld,      &
       rwtop_cld,      &
       tctop_cld,      &
       gsgam2,         &
       gsgam2h,        &
       gam2,           &
       gam2h,          &
       ix,             &
       iy,             &
       iz,             &
       jx,             &
       jy,             &
       jz,             &
       z,              &
       zh,             &
       dt              )
    use mod_misc, only: &
       MISC_gammafunc
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_UNDEF,          &
       CNST_EPS,            &
       CNST_PI,             &
       CNST_Rvap,           &
       CNST_CL,             &
       CNST_LHV0,           &
       CNST_LHS0,           &
       CNST_LHF0,           &
       rho_w => CNST_DWATR, &
       CNST_TEM00,          &
       CNST_Pstd
    use mod_runconf, only: &
       LHV,               &
       LHF,               &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       I_QR,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       opt_water_isotope, & ! [Add] M.Tanoue 20200821
       ISO_STR,           & ! [Add] M.Tanoue 20200905
       ISO_END,           & ! [Add] M.Tanoue 20200905
       ISO1_QV,           & ! [Add] M.Tanoue 20200821
       ISO2_QV,           & ! [Add] M.Tanoue 20200821
       ISO1_QC,           & ! [Add] M.Tanoue 20200821
       ISO2_QC,           & ! [Add] M.Tanoue 20200821
       ISO1_QR,           & ! [Add] M.Tanoue 20200821
       ISO2_QR,           & ! [Add] M.Tanoue 20200821
       ISO1_QI,           & ! [Add] M.Tanoue 20200821
       ISO2_QI,           & ! [Add] M.Tanoue 20200821
       ISO1_QS,           & ! [Add] M.Tanoue 20200821
       ISO2_QS,           & ! [Add] M.Tanoue 20200821
       ISO1_QG,           & ! [Add] M.Tanoue 20200821
       ISO2_QG,           & ! [Add] M.Tanoue 20200821
       CVW
    use mod_isotope, only: &
       NISO,         &
       DIFF_H18O,    &
       DIFF_HDO,     &
       WISO_eqfract, &
       ALP0
    use mod_thrmdyn, only: &
       THRMDYN_cv, &
       THRMDYN_qd, &
       thrmdyn_wait_on => acc_wait_on
    use mod_satadjust, only: &
       SATURATION_psat_liq,  &
       SATURATION_psat_ice,  &
       SATURATION_setrange,  &
       SATURATION_adjustment,&
       satadjustment_wait_on => acc_wait_on
    use mod_precip_transport, only: &
       precip_transport,      &
       precip_transport_new,  &
       precip_transport_euler
!--- MINIAPP: mod_history is the NetCDF/FIO history writer; the mini-app has
!--- no output stream and validates the mp_driver arguments directly. ---
#ifndef MINIAPP
    use mod_history, only: &
       history_in
#endif
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: l_region
    real(RP), intent(inout) :: rhog          (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz        (ijdim,kdim)
    real(RP), intent(inout) :: rhogw         (ijdim,kdim)
    real(RP), intent(inout) :: rhoge         (ijdim,kdim)
    real(RP), intent(inout) :: rhogq         (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: vx            (ijdim,kdim)
    real(RP), intent(in)    :: vy            (ijdim,kdim)
    real(RP), intent(in)    :: vz            (ijdim,kdim)
    real(RP), intent(in)    :: w             (ijdim,kdim)
    real(RP), intent(in)    :: UNCCN         (ijdim,kdim)
    real(RP), intent(in)    :: factor_n_aie  (ijdim)      ! [Add] 2017/06/13 T.Seiki
    real(RP), intent(inout) :: rho           (ijdim,kdim)
    real(RP), intent(inout) :: tem           (ijdim,kdim)
    real(RP), intent(inout) :: pre           (ijdim,kdim)
    real(RP), intent(inout) :: q             (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd            (ijdim,kdim)
    real(RP), intent(out)   :: precip        (ijdim,2)
    real(RP), intent(out)   :: precip_ISO1   (ijdim,2)    ! [Add] M.Tanoue 20200807
    real(RP), intent(out)   :: precip_ISO2   (ijdim,2)    ! [Add] M.Tanoue 20200807
    real(RP), intent(out)   :: precip_rhoe   (ijdim)
    real(RP), intent(out)   :: precip_lh_heat(ijdim)
    real(RP), intent(out)   :: precip_rhophi (ijdim)
    real(RP), intent(out)   :: precip_rhokin (ijdim)
    real(RP), intent(out)   :: precip_flux   (ijdim,kdim,2)
    real(RP), intent(out)   :: evap_liq      (ijdim,kdim) ! 2018/11/08 Yamashita
    real(RP), intent(out)   :: evap_ice      (ijdim,kdim) ! 2018/11/08 Yamashita
    real(RP), intent(out)   :: rceff         (ijdim,kdim)
    real(RP), intent(out)   :: rctop         (ijdim,1)
    real(RP), intent(out)   :: rwtop         (ijdim,1)
    real(RP), intent(out)   :: tctop         (ijdim,1)
    real(RP), intent(out)   :: rceff_cld     (ijdim,kdim)
    real(RP), intent(out)   :: rctop_cld     (ijdim,1)
    real(RP), intent(out)   :: rwtop_cld     (ijdim,1)
    real(RP), intent(out)   :: tctop_cld     (ijdim,1)
    real(RP), intent(in)    :: gsgam2        (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h       (ijdim,kdim)
    real(RP), intent(in)    :: gam2          (ijdim,kdim)
    real(RP), intent(in)    :: gam2h         (ijdim,kdim)
    real(RP), intent(in)    :: ix            (ijdim)
    real(RP), intent(in)    :: iy            (ijdim)
    real(RP), intent(in)    :: iz            (ijdim)
    real(RP), intent(in)    :: jx            (ijdim)
    real(RP), intent(in)    :: jy            (ijdim)
    real(RP), intent(in)    :: jz            (ijdim)
    real(RP), intent(in)    :: z             (ijdim,kdim)
    real(RP), intent(in)    :: zh            (ijdim,kdim)
    real(DP), intent(in)    :: dt

    ! working
    real(RP) :: drhogqv(ijdim,kdim)
    real(RP) :: drhogqc(ijdim,kdim)
    real(RP) :: drhogqi(ijdim,kdim)
    real(RP) :: drhogqr(ijdim,kdim)
    real(RP) :: drhogqs(ijdim,kdim)
    real(RP) :: drhogqg(ijdim,kdim)
    real(RP) :: psatl  (ijdim,kdim)
    real(RP) :: psati  (ijdim,kdim)
    real(RP) :: Nc     (ijdim,kdim)              !< Number concentration of cloud water [1/cc]

    real(RP) :: dens                             !< density
    real(RP) :: temp                             !< T [K]
    real(RP) :: qv                               !< mixing ratio of water vapor  [kg/kg]
    real(RP) :: qc                               !< mixing ratio of liquid water [kg/kg]
    real(RP) :: qr                               !< mixing ratio of rain         [kg/kg]
    real(RP) :: qi                               !< mixing ratio of ice water    [kg/kg]
    real(RP) :: qs                               !< mixing ratio of snow         [kg/kg]
    real(RP) :: qg                               !< mixing ratio of graupel      [kg/kg]
    real(RP) :: qv_t                             !< tendency     of water vapor  [kg/kg/s]
    real(RP) :: qc_t                             !< tendency     of liquid water [kg/kg/s]
    real(RP) :: qr_t                             !< tendency     of rain         [kg/kg/s]
    real(RP) :: qi_t                             !< tendency     of ice water    [kg/kg/s]
    real(RP) :: qs_t                             !< tendency     of snow         [kg/kg/s]
    real(RP) :: qg_t                             !< tendency     of graupel      [kg/kg/s]
    real(RP) :: qsatl                            !< saturated water vapor for liquid water [kg/kg]
    real(RP) :: qsati                            !< saturated water vapor for ice water    [kg/kg]
    real(RP) :: Sliq                             !< saturated ratio S for liquid water [0-1]
    real(RP) :: Sice                             !< saturated ratio S for ice water    [0-1]
    real(RP) :: Rdens                            !< 1 / density
    real(RP) :: rho_fact                         !< density factor
    real(RP) :: temc                             !< T - T0 [K]

    real(RP) :: RLMDr, RLMDr_2, RLMDr_3
    real(RP) :: RLMDs, RLMDs_2, RLMDs_3
    real(RP) :: RLMDg, RLMDg_2, RLMDg_3
    real(RP) :: RLMDr_1br, RLMDr_2br, RLMDr_3br
    real(RP) :: RLMDs_1bs, RLMDs_2bs, RLMDs_3bs
    real(RP) :: RLMDr_dr, RLMDr_3dr, RLMDr_5dr
    real(RP) :: RLMDs_ds, RLMDs_3ds, RLMDs_5ds
    real(RP) :: RLMDg_dg, RLMDg_3dg, RLMDg_5dg
    real(RP) :: RLMDr_7
    real(RP) :: RLMDr_6dr

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, Xs2
    real(RP) :: MOMs_0, MOMs_1, MOMs_2
    real(RP) :: MOMs_0bs, MOMs_1bs, MOMs_2bs
    real(RP) :: MOMs_2ds, MOMs_5ds_h, RMOMs_Vt
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: Vti, Vtr, Vts, Vtg               !< terminal velocity [m/s], upward is positive
    real(RP) :: Esi_mod, Egs_mod                 !< modified accretion efficiency
    real(RP) :: rhoqc                            !< rho * qc
    real(RP) :: Pracw_orig,  Pracw_kk            !< accretion       term by orig  & k-k scheme
    real(RP) :: Praut_berry, Praut_kk            !< auto-conversion term by berry & k-k scheme
    real(RP) :: Dc                               !< relative variance
    real(RP) :: beta_saut_mod
    real(RP) :: betai, betas                     !< sticky parameter for auto-conversion
    real(RP) :: Ka                               !< thermal diffusion coefficient of air
    real(RP) :: Kd                               !< diffusion coefficient of water vapor in air
    real(RP) :: Nu                               !< kinematic viscosity of air
    real(RP) :: Glv, Giv, Gil                    !< thermodynamic function
    real(RP) :: ventr, vents, ventg              !< ventilation factor
    real(RP) :: net, fac, fac_sw
    real(RP) :: zerosw, tmp

    !---< Seiki and Roh (2020) >---
    real(RP) :: Sice_eq                          !< equivalent saturated ratio S for ice water[0-1]
    ! mean diameter used for Cutoff size of riming process
    ! p.150-152 in Wang [2002], Ice Microdynamics, Academic Press
    real(RP) :: Ds_ave, Dg_ave
    real(RP) :: cutoff_s, cutoff_g               !< 1 or 0
    real(RP) :: Pigen_sclw
    real(RP) :: XNi_def, XNi_p07
    real(RP) :: flag_stop_nucleation_qs          !< 1 or 0
    real(RP) :: flag_stop_nucleation_qg          !< 1 or 0
    real(RP) :: flag_shallow_stratiform          !< 1 or 0
    real(RP) :: flags                            !< 1 or 0
    real(RP) :: sw1, sw2

    !---< Bergeron process >---
    real(RP) :: sw_bergeron                      !< if 0C<T<30C, sw=1
    real(RP) :: a1 (ijdim,kdim)                  !<
    real(RP) :: a2 (ijdim,kdim)                  !<
    real(RP) :: ma2(ijdim,kdim)                  !< 1-a2
    real(RP) :: dt1                              !< time during which the an ice particle of 40um grows to 50um
    real(RP) :: Ni50                             !< number concentration of ice particle of 50um

    !---< Explicit ice generation >---
    real(RP) :: sw, rhoqi, XNi, XMi, Di, Ni0, Qi0

    !---< Effective radius >---
    real(RP), parameter :: r2_min = 1.E-10_RP
    real(RP), parameter :: q_min  = 1.E-5_RP
    real(RP)            :: xf_qc                 !< mean mass   of qc [kg]
    real(RP)            :: rf_qc                 !< mean radius of qc [m]
    real(RP)            :: r2_qc                 !< r^2 moment  of qc
    real(RP)            :: r2_qr                 !< r^2 moment  of qr
    real(RP)            :: r3_qc                 !< r^3 moment  of qc
    real(RP)            :: r3_qr                 !< r^3 moment  of qr
    real(RP)            :: dgamma_a, GAM_dgam23, GAM_dgam
    real(RP)            :: coef_dgam, coef_xf, OneThird

    !---< Precipitation >---
    logical  :: preciptation_flag(nqmax)

    real(RP) :: Vt      (ijdim,kdim,nqmax)
    real(RP) :: cva     (ijdim,kdim)
    real(RP) :: rgs     (ijdim,kdim)
    real(RP) :: rgsh    (ijdim,kdim)

    real(RP) :: wk      (wk_nmax)
!    real(RP) :: ml_wk   (ijdim,kdim,wk_nmax) tentatively remove due to the stack overflow
    real(RP) :: ml_Pconv(ijdim,kdim)
    real(RP) :: ml_Pconw(ijdim,kdim)
    real(RP) :: ml_Pconi(ijdim,kdim)

    !---< ISOTOPE >--- [Add] M.Tanoue 20201030
    real(RP) :: qvi  (NISO)                      !< mixing ratio of water vapor isotope  [kg/kg]
    real(RP) :: qci  (NISO)                      !< mixing ratio of liquid water isotope [kg/kg]
    real(RP) :: qri  (NISO)                      !< mixing ratio of rain isotope         [kg/kg]
    real(RP) :: qii  (NISO)                      !< mixing ratio of ice water isotope    [kg/kg]
    real(RP) :: qsi  (NISO)                      !< mixing ratio of snow isotope         [kg/kg]
    real(RP) :: qgi  (NISO)                      !< mixing ratio of graupel isotope      [kg/kg]
    real(RP) :: qvi_t(NISO)                      !< tendency     of water vapor isotope  [kg/kg/s]
    real(RP) :: qci_t(NISO)                      !< tendency     of liquid water isotope [kg/kg/s]
    real(RP) :: qri_t(NISO)                      !< tendency     of rain  isotope        [kg/kg/s]
    real(RP) :: qii_t(NISO)                      !< tendency     of ice water isotope    [kg/kg/s]
    real(RP) :: qsi_t(NISO)                      !< tendency     of snow isotope         [kg/kg/s]
    real(RP) :: qgi_t(NISO)                      !< tendency     of graupel isotope      [kg/kg/s]

    real(RP) :: drhogqvi(ijdim,kdim,NISO)
    real(RP) :: drhogqci(ijdim,kdim,NISO)
    real(RP) :: drhogqii(ijdim,kdim,NISO)
    real(RP) :: drhogqri(ijdim,kdim,NISO)
    real(RP) :: drhogqsi(ijdim,kdim,NISO)
    real(RP) :: drhogqgi(ijdim,kdim,NISO)
    real(RP) :: wk_iso  (wk_nmax,NISO)
    integer  :: iso

    real(RP) :: fidep, fsdep, fgdep
    real(RP) :: ALP(NISO), ALPI(NISO), ALPKL(NISO), ALPKS(NISO)
    real(RP) :: Sleff, Sieff
    real(RP) :: Tau
    real(RP) :: qr_end, qri_end(NISO)
    real(RP) :: Rl_end(NISO), Rl_ini(NISO), Rl_inf(NISO)
    real(RP) :: tmphi(NISO)
    real(RP) :: RDD(NISO)
    real(RP) :: Rqv(NISO), Rqc(NISO), Rqr(NISO)
    real(RP) :: Rqi(NISO), Rqs(NISO), Rqg(NISO)
    integer  :: I_QVI(NISO), I_QCI(NISO), I_QRI(NISO)
    integer  :: I_QII(NISO), I_QSI(NISO), I_QGI(NISO)

    ! SIMD parameter
    integer, parameter :: simdlen = 512
    integer            :: blk, vec, veclen

    real(RP) :: UNDEF, EPS, PI, Rvap, CL, LHV0, LHS0, LHF0, TEM00, PRE00
    real(RP) :: Nc_def0, Ds_cutoff0, Dg_cutoff0, CRIT_SCLW_W0
    integer  :: ij, k, nq, ip
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(qd,precip,precip_ISO1,precip_ISO2,precip_rhoe,precip_lh_heat) &
    !$acc pcreate(precip_rhophi,precip_rhokin,precip_flux,evap_liq,evap_ice) &
    !$acc pcreate(rceff,rctop,rwtop,tctop,rceff_cld,rctop_cld,rwtop_cld,tctop_cld) &
    !$acc pcreate(drhogqv,drhogqc,drhogqi,drhogqr,drhogqs,drhogqg) &
    !$acc pcreate(psatl,qsatl,psati,qsati,Nc,a1,a2,ma2) &
    !$acc pcreate(Vt,cva,rgs,rgsh) &
    !$acc pcreate(wk,ml_Pconv,ml_Pconw,ml_Pconi) &
    !$acc pcopyin(vx,vy,vz,w,UNCCN,factor_n_aie) &
    !$acc pcopyin(gsgam2,gsgam2h,gam2,gam2h,ix,iy,iz,jx,jy,jz,z,zh) &
    !$acc pcopyin(coef_a,coef_b,CVW,coef_at,coef_bt) &
    !$acc pcopy(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq) &
    !$acc pcopy(rho,tem,pre,q)

    call PROF_rapstart('____MP_NSW6',3)
    thrmdyn_wait_on = .false.
    satadjustment_wait_on = .false.

    UNDEF = CNST_UNDEF
    EPS   = CNST_EPS
    PI    = CNST_PI
    Rvap  = CNST_Rvap
    CL    = CNST_CL
    LHV0  = CNST_LHV0
    LHS0  = CNST_LHS0
    LHF0  = CNST_LHF0
    TEM00 = CNST_TEM00
    PRE00 = CNST_Pstd

    Nc_def0      = Nc_def
    Ds_cutoff0   = Ds_cutoff
    Dg_cutoff0   = Dg_cutoff
    CRIT_SCLW_W0 = CRIT_SCLW_W

    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20201211
       RDD  (1) = DIFF_H18O
       RDD  (2) = DIFF_HDO
       I_QVI(1) = ISO1_QV
       I_QVI(2) = ISO2_QV
       I_QCI(1) = ISO1_QC
       I_QCI(2) = ISO2_QC
       I_QRI(1) = ISO1_QR
       I_QRI(2) = ISO2_QR
       I_QII(1) = ISO1_QI
       I_QII(2) = ISO2_QI
       I_QSI(1) = ISO1_QS
       I_QSI(2) = ISO2_QS
       I_QGI(1) = ISO1_QG
       I_QGI(2) = ISO2_QG
    endif

#ifdef MORETIMER
call PROF_rapstart('____MP_NSW6_02z')
#endif

    call negative_filter( ijdim,         & ! [IN]
                          rhog  (:,:),   & ! [INOUT]
                          rhoge (:,:),   & ! [INOUT]
                          rhogq (:,:,:), & ! [INOUT]
                          rho   (:,:),   & ! [INOUT]
                          tem   (:,:),   & ! [IN]
                          pre   (:,:),   & ! [INOUT]
                          q     (:,:,:), & ! [INOUT]
                          gsgam2(:,:)    ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_02z')
call PROF_rapstart('____MP_NSW6_03')
#endif

    if ( OPT_INDIR ) then ! aerosol indirect effect
!OCL XFILL, NOPREFETCH
       !$acc kernels async(1)pcopy(Nc) pcopyin(UNCCN,factor_n_aie)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,Nc,UNCCN,Nc_def0,factor_n_aie)
       do k  = 1, kdim
       do ij = 1, ijdim
          Nc(ij,k) = max( UNCCN(ij,k)*1.E-6_RP*factor_n_aie(ij), Nc_def0 ) ! [#/m3]->[#/cc]
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    else
!OCL XFILL, NOPREFETCH
       !$acc kernels async(1) pcopy(Nc) pcopyin(factor_n_aie)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,Nc,Nc_def0,factor_n_aie)
       do k  = 1, kdim
       do ij = 1, ijdim
          Nc(ij,k) = Nc_def0 * factor_n_aie(ij)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_03')
call PROF_rapstart('____MP_NSW6_04z')
#endif

    ! saturation water contensts
    call SATURATION_psat_liq( ijdim, kdim, tem(:,:), psatl(:,:) )
    call SATURATION_psat_ice( ijdim, kdim, tem(:,:), psati(:,:) )

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_04z')
call PROF_rapstart('____MP_NSW6_05z')
#endif

    ! Bergeron process parameters
    call Bergeron_param( ijdim,    & ! [IN]
                         tem(:,:), & ! [IN]
                         a1 (:,:), & ! [OUT]
                         a2 (:,:), & ! [OUT]
                         ma2(:,:)  ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_05z')
#endif

    call PROF_rapstart('____MP_NSW6_CORE',3)

    !$acc kernels async(1) pcopy(drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg,Vt,evap_liq,evap_ice) &
    !$acc pcopyin(rhog,rho,tem,pre,q,Nc,psatl,psati,a1,a2,ma2,coef_a,coef_b,coef_at,coef_bt)
    !$acc loop independent
    !$omp parallel do default(none), &
    !$omp private(blk,vec,veclen)                                                                             &
    !$omp private(ij,k,ip,dens,temp,qv,qc,qr,qi,qs,qg,qsatl,qsati,Sliq,Sice,Rdens,rho_fact,temc,sw_bergeron)  &
    !$omp private(RLMDr,RLMDr_dr,RLMDr_2,RLMDr_3,RLMDr_7,RLMDr_1br,RLMDr_2br,RLMDr_3br,RLMDr_3dr,             &
    !$omp         RLMDr_5dr,RLMDr_6dr,RLMDs,RLMDs_ds,RLMDs_2,RLMDs_3,RLMDs_1bs,RLMDs_2bs,RLMDs_3bs,           &
    !$omp         RLMDs_3ds,RLMDs_5ds,RLMDg,RLMDg_dg,RLMDg_2,RLMDg_3,RLMDg_3dg,RLMDg_5dg,                     &
    !$omp         MOMs_0,MOMs_1,MOMs_2,MOMs_0bs,MOMs_1bs,MOMs_2bs,MOMs_2ds,MOMs_5ds_h,RMOMs_Vt,               &
    !$omp         coef_bt,coef_at,Xs2,tems,loga_,b_,nm,zerosw,sw1,sw2,tmp)                                    &
    !$omp private(Vti,Vtr,Vts,Vtg,Esi_mod,Egs_mod,Pracw_orig,Pracw_kk,rhoqc,Dc,Praut_berry,Praut_kk,          &
    !$omp         beta_saut_mod,betai,betas,Ka,Kd,NU,Glv,Giv,Gil,ventr,vents,ventg,dt1,Ni50,sw,rhoqi,         &
    !$omp         XNi,XMi,Di,Ni0,Qi0,net,fac,fac_sw,qv_t,qc_t,qr_t,qi_t,qs_t,qg_t,wk)                         &
    !$omp private(Ds_ave,Dg_ave,cutoff_s,cutoff_g,Pigen_sclw,Sice_eq,XNi_def,XNi_p07,                         &
    !$omp         flag_shallow_stratiform,flag_stop_nucleation_qs,flag_stop_nucleation_qg,flags)              &
!   !$omp private(ml_wk)                                                                                      &
#ifdef NSW6_WISO
    !$omp private(qvi,qci,qri,qii,qsi,qgi,qvi_t,qci_t,qri_t,qii_t,qsi_t,qgi_t,iso,                            &
    !$omp         Rqv,Rqc,Rqr,Rqi,Rqs,Rqg,ALP,ALPI,ALPKL,ALPKS,Sleff,Sieff,                                   &
    !$omp         wk_iso,qr_end,qri_end,Rl_end,tmphi,Tau,Rl_inf,Rl_ini,fidep,fsdep,fgdep)                     &
    !$omp shared(drhogqvi,drhogqci,drhogqri,drhogqii,drhogqsi,drhogqgi,                                       &
    !$omp        RDD,I_QVI,I_QCI,I_QRI,I_QII,I_QSI,I_QGI,opt_water_isotope)                                   &
#endif
    !$omp shared(ijdim,kmin,kmax,rho,tem,pre,q,psatl,psati,Nc,dt,sw_expice,sw_roh2014,sw_kk2000,sw_constVti,  &
    !$omp        drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg,rhog,Vt,evap_liq,evap_ice,constVti,          &
    !$omp        UNDEF,EPS,PI,Rvap,CL,LHV0,LHS0,LHF0,TEM00,PRE00,I_QV,I_QC,I_QR,I_QI,I_QS,I_QG)               &
    !$omp shared(GAM,GAM_2,GAM_3,GAM_1br,GAM_2br,GAM_3br,GAM_3dr,GAM_6dr,GAM_1brdr,GAM_5dr_h,GAM_1bs,GAM_2bs, &
    !$omp        GAM_3bs,GAM_3ds,GAM_1bsds,GAM_5ds_h,GAM_1bg,GAM_3dg,GAM_1bgdg,GAM_5dg_h,coef_a,coef_b,ln10,  &
    !$omp        N0r,N0s,N0g,Ar,As,Ag,Br,Bs,Bg,Cr,Cs,Cg,Dr,Ds,Dg,Eiw,Erw,Esw,Egw,Eri,Esi,Egi,Esr,Egr,Egs,     &
    !$omp        gamma_sacr,gamma_gacs,mi,beta_saut,gamma_saut,beta_gaut,gamma_gaut,qicrt_saut,qscrt_gaut,    &
    !$omp        f1r,f2r,f1s,f2s,f1g,f2g,A_frz,B_frz,mi40,mi50,vti50,Ri50,a1,a2,ma2)                          &
    !$omp shared(sw_sr2020,sw_phillips2007,sw_beta_saut,Ds_cutoff0,Dg_cutoff0,w,CRIT_SCLW_W0)
    do k  = kmin, kmax
    !$acc loop independent
    do blk = 1, ijdim, simdlen ! cache blocking
       veclen = min( ijdim-blk+1, simdlen )

       !ocl loop_fission_target(ls)
       !$acc loop independent private(wk,coef_at,coef_bt)
       do vec = 1, veclen
       ij = blk + vec - 1

       dens     = rho(ij,k)
       temp     = tem(ij,k)
       qv       = max( q(ij,k,I_QV), 0.0_RP )
       qc       = max( q(ij,k,I_QC), 0.0_RP )
       qr       = max( q(ij,k,I_QR), 0.0_RP )
       qi       = max( q(ij,k,I_QI), 0.0_RP )
       qs       = max( q(ij,k,I_QS), 0.0_RP )
       qg       = max( q(ij,k,I_QG), 0.0_RP )

       ! saturation ratio S
       qsatl    = psatl(ij,k) / ( dens * Rvap * temp )
       qsati    = psati(ij,k) / ( dens * Rvap * temp )
       Sliq     = qv / max( qsatl, EPS )
       Sice     = qv / max( qsati, EPS )

       Rdens    = 1.0_RP / dens
       rho_fact = sqrt( dens00 * Rdens )
       temc     = temp - TEM00

       wk(I_delta1) = ( 0.5_RP + sign(0.5_RP, qr - 1.E-4_RP ) )

       wk(I_delta2) = ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qr ) ) &
                    * ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qs ) )

       wk(I_spsati) = 0.5_RP + sign(0.5_RP, Sice - 1.0_RP )

       wk(I_iceflg) = 0.5_RP - sign( 0.5_RP, temc ) ! 0: warm, 1: ice

       wk(I_dqv_dt) = qv / dt
       wk(I_dqc_dt) = qc / dt
       wk(I_dqr_dt) = qr / dt
       wk(I_dqi_dt) = qi / dt
       wk(I_dqs_dt) = qs / dt
       wk(I_dqg_dt) = qg / dt

       sw_bergeron = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
                   * ( 0.5_RP - sign(0.5_RP, temc           ) ) &
                   * ( 1.0_RP - sw_expice                     )

       ! slope parameter lambda (Rain)
       zerosw = 0.5_RP - sign(0.5_RP, qr - 1.E-12_RP )
       RLMDr  = sqrt(sqrt( dens * qr / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDr_dr  = sqrt( RLMDr )       ! **Dr
       RLMDr_2   = RLMDr**2
       RLMDr_3   = RLMDr**3
       RLMDr_7   = RLMDr**7
       RLMDr_1br = RLMDr**4 ! (1+Br)
       RLMDr_2br = RLMDr**5 ! (2+Br)
       RLMDr_3br = RLMDr**6 ! (3+Br)
       RLMDr_3dr = RLMDr**3 * RLMDr_dr
       RLMDr_5dr = RLMDr**5 * RLMDr_dr
       RLMDr_6dr = RLMDr**6 * RLMDr_dr

       ! slope parameter lambda (Snow)
       zerosw = 0.5_RP - sign(0.5_RP, dens * qs - 1.E-12_RP )
       RLMDs  = sqrt(sqrt( dens * qs / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
       RLMDs_2   = RLMDs**2
       RLMDs_3   = RLMDs**3
       RLMDs_1bs = RLMDs**4 ! (1+Bs)
       RLMDs_2bs = RLMDs**5 ! (2+Bs)
       RLMDs_3bs = RLMDs**6 ! (3+Bs)
       RLMDs_3ds = RLMDs**3 * RLMDs_ds
       RLMDs_5ds = RLMDs**5 * RLMDs_ds

       MOMs_0     = N0s * GAM       * RLMDs           ! Ns * 0th moment
       MOMs_1     = N0s * GAM_2     * RLMDs_2         ! Ns * 1st moment
       MOMs_2     = N0s * GAM_3     * RLMDs_3         ! Ns * 2nd moment
       MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
       MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs       ! Ns * 1+bs
       MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs       ! Ns * 2+bs
       MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds       ! Ns * 2+ds
       MOMs_5ds_h = N0s * GAM_5ds_h * sqrt(RLMDs_5ds) ! Ns * (5+ds)/2
       RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds

       !---< modification by Roh and Satoh (2014) >---
       ! bimodal size distribution of snow
       Xs2 = dens * qs / As

       tems       = min( -0.1_RP, temc )
       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
       coef_at(4) = coef_a(10)
       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
       coef_bt(4) = coef_b(10)
       ! 0th moment
       loga_  = coef_at(1)
          b_  = coef_bt(1)
       MOMs_0 = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_0
       ! 1st moment
       nm = 1.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_1 = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_1
       ! 2nd moment
       MOMs_2 = (        sw_roh2014 ) * Xs2 &
              + ( 1.0_RP-sw_roh2014 ) * MOMs_2
       ! 0 + Bs(=2) moment
       nm = 2.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_0bs = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_0bs
       ! 1 + Bs(=2) moment
       nm = 3.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_1bs = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_1bs
       ! 2 + Bs(=2) moment
       nm = 4.0_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_2bs = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_2bs
       ! 2 + Ds(=0.25) moment
       nm = 2.25_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_2ds = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
                + ( 1.0_RP-sw_roh2014 ) * MOMs_2ds
       ! ( 3 + Ds(=0.25) ) / 2  moment
       nm = 1.625_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_5ds_h = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ &
                  + ( 1.0_RP-sw_roh2014 ) * MOMs_5ds_h
       ! Bs(=2) + Ds(=0.25) moment
       nm = 2.25_RP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       RMOMs_Vt = (        sw_roh2014 ) * exp(ln10*loga_) * Xs2**b_ / ( MOMs_0bs + zerosw ) &
                + ( 1.0_RP-sw_roh2014 ) * RMOMs_Vt

       ! slope parameter lambda (Graupel)
       zerosw = 0.5_RP - sign(0.5_RP, qg - 1.E-12_RP )
       RLMDg  = sqrt(sqrt( dens * qg / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

       RLMDg_dg  = sqrt( RLMDg )       ! **Dg
       RLMDg_2   = RLMDg**2
       RLMDg_3   = RLMDg**3
       RLMDg_3dg = RLMDg**3 * RLMDg_dg
       RLMDg_5dg = RLMDg**5 * RLMDg_dg

       wk(I_RLMDr) = RLMDr
       wk(I_RLMDs) = RLMDs
       wk(I_RLMDg) = RLMDg

       !---< terminal velocity >---
       Vti = ( 1.0_RP-sw_constVti ) * (-3.29_RP) * (dens*qi)**0.16_RP &
           + (        sw_constVti ) * (constVti)
       Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       Vts = -Cs * rho_fact * RMOMs_Vt
       Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg

       ! [Add] 2018/02/08 T.Seiki
       Ds_ave = MOMs_1 / max( MOMs_0, EPS )
       Dg_ave = RLMDg

       cutoff_s = 0.5_RP + sign(0.5_RP, Ds_ave - Ds_cutoff0 ) ! if Ds_ave < Ds_cutoff, sw=0, default is always 1
       cutoff_g = 0.5_RP + sign(0.5_RP, Dg_ave - Dg_cutoff0 ) ! if Dg_ave < Dg_cutoff, sw=0, default is always 1

       !---< Nucleation >---
       ! [Pigen] ice nucleation
       Ni0 = max( exp(-0.1_RP*temc), 1.0_RP ) * 1000.0_RP
       Qi0 = 4.92E-11_RP * Ni0**1.33_RP * Rdens

       Pigen_sclw = max( min( Qi0-qi, qv-qsati ), 0.0_RP ) / dt

       wk(I_Pigen) = ( 1.0_RP - sw_sr2020 ) * Pigen_sclw

       !---< Accretion >---
       Esi_mod = min( Esi, Esi * exp( gamma_sacr * temc ) )
       Egs_mod = min( Egs, Egs * exp( gamma_gacs * temc ) )

       ! [Pracw] accretion rate of cloud water by rain
       Pracw_orig = qc * 0.25_RP * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       Pracw_kk   = 67.0_RP * (qc*qr)**1.15_RP ! eq.(33) in KK(2000)

       ! switch orig / k-k scheme
       wk(I_Pracw) = ( 1.0_RP-sw_kk2000 ) * Pracw_orig &
                   + (        sw_kk2000 ) * Pracw_kk

       ! [Psacw] accretion rate of cloud water by snow
       wk(I_Psacw) = qc * 0.25_RP * PI * Esw       * Cs * MOMs_2ds            * rho_fact * cutoff_s ! [Mod] 2018/02/08 T.Seiki

       ! [Pgacw] accretion rate of cloud water by graupel
       wk(I_Pgacw) = qc * 0.25_RP * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact * cutoff_g ! [Mod] 2018/02/08 T.Seiki

       ! [Praci] accretion rate of cloud ice by rain
       wk(I_Praci) = qi * 0.25_RP * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

       ! [Psaci] accretion rate of cloud ice by snow
       wk(I_Psaci) = qi * 0.25_RP * PI * Esi_mod   * Cs * MOMs_2ds            * rho_fact

       ! [Pgaci] accretion rate of cloud ice by grupel
       wk(I_Pgaci) = qi * 0.25_RP * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

       ! [Piacr] accretion rate of rain by cloud ice
       wk(I_Piacr) = qi * Ar / Mi * 0.25_RP * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact

       ! [Psacr] accretion rate of rain by snow
       wk(I_Psacr) = Ar * 0.25_RP * PI * Rdens * Esr * N0r       * abs(Vtr-Vts) &
                   * (          GAM_1br * RLMDr_1br * MOMs_2          &
                     + 2.0_RP * GAM_2br * RLMDr_2br * MOMs_1          &
                     +          GAM_3br * RLMDr_3br * MOMs_0          )

       ! [Pgacr] accretion rate of rain by graupel
       wk(I_Pgacr) = Ar * 0.25_RP * PI * Rdens * Egr * N0g * N0r * abs(Vtg-Vtr) &
                   * (          GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3 &
                     + 2.0_RP * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2 &
                     +          GAM_3br * RLMDr_3br * GAM   * RLMDg   )

       ! [Pracs] accretion rate of snow by rain
       wk(I_Pracs) = As * 0.25_RP * PI * Rdens * Esr       *  N0r * abs(Vtr-Vts) &
                   * (          MOMs_0bs            * GAM_3 * RLMDr_3 &
                     + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDr_2 &
                     +          MOMs_2bs            * GAM   * RLMDr   )

       ! [Pgacs] accretion rate of snow by graupel
       wk(I_Pgacs) = As * 0.25_RP * PI * Rdens * Egs_mod   * N0g * abs(Vtg-Vts) &
                   * (          MOMs_0bs            * GAM_3 * RLMDg_3 &
                     + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDg_2 &
                     +          MOMs_2bs            * GAM   * RLMDg   )

       !---< Autoconversion >---

       ! [Praut] auto-conversion rate from cloud water to rain
       rhoqc = dens * qc * 1000.0_RP ! [g/m3]
       Dc    = 0.146_RP - 5.964E-2_RP * log( Nc(ij,k) / 2000.0_RP )
       Praut_berry = Rdens * 1.67E-5_RP * rhoqc * rhoqc / ( 5.0_RP + 3.66E-2_RP * Nc(ij,k) / ( Dc * rhoqc + EPS ) )

       Praut_kk    = 1350.0_RP * qc**2.47_RP * Nc(ij,k)**(-1.79_RP) ! eq.(29) in KK(2000)

       ! switch berry / k-k scheme
       wk(I_Praut) = ( 1.0_RP-sw_kk2000 ) * Praut_berry &
                   + (        sw_kk2000 ) * Praut_kk

       ! [Psaut] auto-conversion rate from cloud ice to snow

       ! [ADD] 2020/08/19 A.Noda b_aut based on NDW6 in GL10, suggested by T.Seiki
       beta_saut_mod = ( 1.0_RP-sw_beta_saut ) * beta_saut &
                     + (        sw_beta_saut ) * 5.301_RP * (dens*qi)**0.9087_RP

       betai = min( beta_saut_mod, beta_saut_mod * exp( gamma_saut * temc ) )
       wk(I_Psaut) = max( betai*(qi-qicrt_saut), 0.0_RP )

       ! [Pgaut] auto-conversion rate from snow to graupel
       !betas = min( beta_gaut, beta_gaut * exp( gamma_gaut * temc ) )
       !wk(I_Pgaut) = max( betas*(qs-qscrt_gaut), 0.0_RP ) [mod] 2018/02/08 T.Seiki
       wk(I_Pgaut) = 0.0_RP

       !---< Evaporation, Sublimation, Melting, and Freezing >---

       Ka  = ( Ka0 + dKa_dT * temc )
       Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
       NU  = ( nu0 + dnu_dT * temc ) * Rdens

       Glv = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsatl) )
       Giv = 1.0_RP / ( LHS0/(Ka*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsati) )
       Gil = ( Ka * temc ) / LHF0

       ! [Prevp] evaporation rate of rain
       ventr = f1r * GAM_2 * RLMDr_2 + f2r * sqrt( Cr * rho_fact / NU * RLMDr_5dr ) * GAM_5dr_h

       wk(I_Prevp) = 2.0_RP * PI * Rdens * N0r * ( 1.0_RP-min(Sliq,1.0_RP) ) * Glv * ventr

       ! [Pidep,Pisub] deposition/sublimation rate for ice
       rhoqi = max(dens*qi,EPS)

       XNi_def = min( max( 5.38E+7_RP * rhoqi**0.75_RP, 1.E+3_RP ), 1.E+6_RP )

       ! [Add] 2018/02/08 T.Seiki, Phillips et al. (2007), JAS
       Sice_eq = min( 3.0_RP, (qv+qi) / max(qsati,EPS) )

       sw1 = ( 0.5_RP - sign(0.5_RP, temc + 30.0_RP ) ) ! if T < -30C, sw1=1
       sw2 = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
           * ( 0.5_RP - sign(0.5_RP, temc           ) ) ! if -30C < T < 0C, sw2=1

       XNi_p07 = (        sw1 ) * 1.E+3_RP * exp( 12.96_RP*(Sice_eq-1.1_RP) * 0.3_RP   )           & ! if        temc < -30C
               + (        sw2 ) * 1.E+3_RP * exp( 12.96_RP*(Sice_eq-1.0_RP) - 0.639_RP ) * 0.06_RP & ! if -30C < temc <   0C
               + ( 1.0_RP-sw1 ) * ( 1.0_RP-sw2 ) * 1.E+3_RP                                          ! if   0C < temc

       XNi = ( 1.0_RP-sw_phillips2007 ) * XNi_def &
           + (        sw_phillips2007 ) * XNi_p07
       XMi = rhoqi / XNi
       Di  = min( Di_a * sqrt(XMi), Di_max )

       tmp = ( 1.0_RP - sw_sr2020 ) * 4.0_RP * Di * XNi * Rdens * ( Sice-1.0_RP ) * Giv ! [mod] 2018/02/08 T.Seiki

       wk(I_Pidep) = (        wk(I_spsati) ) * ( tmp) ! Sice > 1
       wk(I_Pisub) = ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

       ! [Pihom] homogenious freezing at T < -40C
       sw = ( 0.5_RP - sign(0.5_RP, temc + 40.0_RP ) ) ! if T < -40C, sw=1

       wk(I_Pihom) = sw * qc / dt

       ! [Pihtr] heteroginous freezing at -40C < T < 0C
       sw = ( 0.5_RP + sign(0.5_RP, temc + 40.0_RP ) ) &
          * ( 0.5_RP - sign(0.5_RP, temc           ) ) ! if -40C < T < 0C, sw=1

       wk(I_Pihtr) = sw * ( dens / rho_w * qc**2 / ( Nc(ij,k) * 1.E+6_RP ) ) & ! [Mod] 2020-11-17 T.Seiki, delete Nc_ihtr to avoid inconsistency
                   * B_frz * min( exp(-A_frz*temc)-1.0_RP, 1.E+18_RP )

       ! [Pimlt] ice melting at T > 0C
       sw = ( 0.5_RP + sign(0.5_RP, temc           ) ) ! if T > 0C, sw=1

       wk(I_Pimlt) = sw * qi / dt

       ! [Psdep,Pssub] deposition/sublimation rate for snow
       vents = f1s * MOMs_1          + f2s * sqrt( Cs * rho_fact / NU             ) * MOMs_5ds_h

       tmp = 2.0_RP * PI * Rdens *       ( Sice-1.0_RP ) * Giv * vents

       flag_shallow_stratiform = ( 0.5_RP - sign(0.5_RP, w(ij,k) - CRIT_SCLW_W0 ) ) & ! [mod] 2021/03/01 A.Noda
                               * ( 0.5_RP + sign(0.5_RP, qc - 1.E-6_RP          ) )   ! [Add] 2018/02/08 T.Seiki

       flag_stop_nucleation_qs = 0.5_RP - sign(0.5_RP, qs - 1.E-6_RP ) ! [Add] 2018/02/08 T.Seiki

       flags = 1.0_RP - ( 1.0_RP - flag_stop_nucleation_qs ) * ( 1.0_RP - flag_shallow_stratiform ) ! flag_stop_nucleation_qs .OR. flag_shallow_stratiform

       wk(I_Psdep) = ( 1.0_RP - sw_sr2020 * flags ) * (        wk(I_spsati) ) * ( tmp) ! Sice > 1
       wk(I_Pssub) =                                  ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

       ! [Psmlt] melting rate of snow
       wk(I_Psmlt) = 2.0_RP * PI * Rdens *       Gil * vents &
                   + CL * temc / LHF0 * ( wk(I_Psacw) + wk(I_Psacr) )
       wk(I_Psmlt) = max( wk(I_Psmlt), 0.0_RP )

       ! [Pgdep/pgsub] deposition/sublimation rate for graupel
       ventg = f1g * GAM_2 * RLMDg_2 + f2g * sqrt( Cg * rho_fact / NU * RLMDg_5dg ) * GAM_5dg_h

       tmp = 2.0_RP * PI * Rdens * N0g * ( Sice-1.0_RP ) * Giv * ventg

       flag_stop_nucleation_qg = 0.5_RP - sign(0.5_RP, qg - 1.E-6_RP ) ! [Add] 2018/02/08 T.Seiki

       flags = 1.0_RP - ( 1.0_RP - flag_stop_nucleation_qg ) * ( 1.0_RP - flag_shallow_stratiform ) ! flag_stop_nucleation_qg .OR. flag_shallow_stratiform

       wk(I_Pgdep) = ( 1.0_RP - sw_sr2020 * flags ) * (        wk(I_spsati) ) * ( tmp) ! Sice > 1
       wk(I_Pgsub) =                                  ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

       ! [Pgmlt] melting rate of graupel
       wk(I_Pgmlt) = 2.0_RP * PI * Rdens * N0g * Gil * ventg &
                   + CL * temc / LHF0 * ( wk(I_Pgacw) + wk(I_Pgacr) )
       wk(I_Pgmlt) = max( wk(I_Pgmlt), 0.0_RP )

       ! [Pgfrz] freezing rate of graupel
       wk(I_Pgfrz) = 120.0_RP * PI * Ar * Rdens * N0r * B_frz * min( exp(-A_frz*temc)-1.0_RP, 1.E+18_RP ) * RLMDr_7

       ! [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
       dt1  = ( mi50**ma2(ij,k) - mi40**ma2(ij,k) ) / ( a1(ij,k) * ma2(ij,k) )
       Ni50 = qi * dt / ( mi50 * dt1 )

       wk(I_Psfw) = ( 1.0_RP - sw_sr2020 ) * Ni50 * ( a1(ij,k) * mi50**a2(ij,k)                &
                                                    + PI * Eiw * dens * qc * Ri50*Ri50 * vti50 )
       wk(I_Psfi) = ( 1.0_RP - sw_sr2020 ) * qi / dt1

#ifdef NSW6_WISO
       if ( opt_water_isotope ) then ! [Add] M.Tanoue 20201211
          !<= [Add] M.Tanoue 20200215
          if ( qv > 0.0_RP ) then
             do iso = 1, NISO
                qvi(iso) = max( q(ij,k,I_QVI(iso)), 0.0_RP )
                Rqv(iso) = qvi(iso) / qv
                wk_iso(I_dqv_dt, iso) = qvi(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qvi(iso) = 0.0_RP
                Rqv(iso) = 0.0_RP
                wk_iso(I_dqv_dt, iso) = 0.0_RP
             enddo
          endif
          if ( qc > 0.0_RP ) then
             do iso = 1, NISO
                qci(iso) = max( q(ij,k,I_QCI(iso)), 0.0_RP )
                Rqc(iso) = qci(iso) / qc
                wk_iso(I_dqc_dt, iso) = qci(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qci(iso) = 0.0_RP
                Rqc(iso) = 0.0_RP
                wk_iso(I_dqc_dt, iso) = 0.0_RP
             enddo
          endif
          if ( qr > 0.0_RP ) then
             do iso = 1, NISO
                qri(iso) = max( q(ij,k,I_QRI(iso)), 0.0_RP )
                Rqr(iso) = qri(iso) / qr
                wk_iso(I_dqr_dt, iso) = qri(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qri(iso) = 0.0_RP
                Rqr(iso) = 0.0_RP
                wk_iso(I_dqr_dt, iso) = 0.0_RP
             enddo
          endif
          if ( qi > 0.0_RP ) then
             do iso = 1, NISO
                qii(iso) = max( q(ij,k,I_QII(iso)), 0.0_RP )
                Rqi(iso) = qii(iso) / qi
                wk_iso(I_dqi_dt, iso) = qii(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qii(iso) = 0.0_RP
                Rqi(iso) = 0.0_RP
                wk_iso(I_dqi_dt, iso) = 0.0_RP
             enddo
          endif
          if ( qs > 0.0_RP ) then
             do iso = 1, NISO
                qsi(iso) = max( q(ij,k,I_QSI(iso)), 0.0_RP )
                Rqs(iso) = qsi(iso) / qs
                wk_iso(I_dqs_dt, iso) = qsi(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qsi(iso) = 0.0_RP
                Rqs(iso) = 0.0_RP
                wk_iso(I_dqs_dt, iso) = 0.0_RP
             enddo
          endif
          if ( qg > 0.0_RP ) then
             do iso = 1, NISO
                qgi(iso) = max( q(ij,k,I_QGI(iso)), 0.0_RP )
                Rqg(iso) = qgi(iso) / qg
                wk_iso(I_dqg_dt, iso) = qgi(iso) / dt
             enddo
          else
             do iso = 1, NISO
                qgi(iso) = 0.0_RP
                Rqg(iso) = 0.0_RP
                wk_iso(I_dqg_dt, iso) = 0.0_RP
             enddo
          endif
          !=> [Add] M.Tanoue 20200215

          !---< Phase change with fractionation >---

          call WISO_eqfract( temp, real(ALP0,RP_G), ALP, ALPI )
          Sleff = (1.0_RP - Glv * (1.0_RP - Sliq ** (-1.0_RP)) ) ** (-1.0_RP)
          Sieff = (1.0_RP - Giv * (1.0_RP - Sice ** (-1.0_RP)) ) ** (-1.0_RP)
          do iso = 1, NISO
             ALPKL(iso) = ALP (iso) / ( 1.0_RP + (ALP (iso) * RDD(iso) - 1.0_RP ) * (1.0_RP-Sleff**(-1.0_RP)) )
             ALPKS(iso) = ALPI(iso) / ( 1.0_RP + (ALPI(iso) * RDD(iso) - 1.0_RP ) * (1.0_RP-Sieff**(-1.0_RP)) )
          enddo

          ! raindrop evaporation
          if ( qr > 0.0_RP ) then
             Tau = 2.0_RP * PI * Rdens * N0r * ventr / qr
             qr_end = qr - wk(I_Prevp) * dt
             do iso = 1, NISO
                tmphi  (iso) = min(Sliq, 1.0_RP) + Glv * (ALP(iso)*RDD(iso)-1.0_RP) * (min(Sliq, 1.0_RP)-1.0_RP)
                Rl_ini (iso) = qri(iso) / qr
                Rl_inf (iso) = qvi(iso) / qv * ALP(iso) * min(Sliq, 1.0_RP) / tmphi(iso)
                Rl_end (iso) = Rl_inf(iso) + (Rl_ini(iso) - Rl_inf(iso)) * exp(-dt*Tau/ALP(iso)*tmphi(iso))
                qri_end(iso) = max( Rl_end(iso)*qr_end, 0.0_RP )
                wk_ISO(I_Prevp, iso) = max( ( qri(iso)-qri_end(iso) ) / dt, 0.0_RP )
             enddo
          else
             wk_ISO(I_Prevp, 1) = 0.0_RP
             wk_ISO(I_Prevp, 2) = 0.0_RP
          endif

          !<= [Add] M.Tanoue 20220218
          if ( OPT_SCLW ) then
             do iso = 1, NISO
                wk_ISO(I_Pigen, iso) = 0.0_RP
                wk_ISO(I_Psfw , iso) = 0.0_RP
             enddo
          else
             do iso = 1, NISO
                ! Ice generation
                !wk_ISO(I_Pigen, iso) = max(wk(I_Pigen), 0.0_RP) * ALPKS(iso) * Rqv(iso)
                wk_ISO(I_Pigen, iso) = max(wk(I_Pigen), 0.0_RP) * Rqv(iso)
                ! Bergeron process
                wk_ISO(I_Psfw , iso) = Ni50 * Rqc(iso) * ( ALPKS(iso)/ALPKL(iso) * a1(ij,k) * mi50**a2(ij,k) &
                    + PI * Eiw * dens * qc * Ri50*Ri50 * vti50 )
             enddo
          endif
          !=> [Add] M.Tanoue 20220218>

          if ( qv > 0.0_RP ) then
             fidep = wk(I_Pidep)*dt/qv
             fsdep = wk(I_Psdep)*dt/qv
             fgdep = wk(I_Pgdep)*dt/qv
             do iso = 1, NISO
                wk_ISO(I_Pidep, iso) = qvi(iso)/dt * (1.0_RP - (1.0_RP - fidep)**ALPKS(iso))
                wk_ISO(I_Psdep, iso) = qvi(iso)/dt * (1.0_RP - (1.0_RP - fsdep)**ALPKS(iso))
                wk_ISO(I_Pgdep, iso) = qvi(iso)/dt * (1.0_RP - (1.0_RP - fgdep)**ALPKS(iso))
             enddo
          else
             do iso = 1, NISO
                wk_ISO(I_Pidep, iso) = 0.0_RP
                wk_ISO(I_Psdep, iso) = 0.0_RP
                wk_ISO(I_Pgdep, iso) = 0.0_RP
             enddo
          endif
          !---< Phase change without fractionation >---

          do iso = 1, NISO
             wk_ISO(I_Pracw, iso) = wk(I_Pracw) * Rqc(iso)
             wk_ISO(I_Psacw, iso) = wk(I_Psacw) * Rqc(iso)
             wk_ISO(I_Pgacw, iso) = wk(I_Pgacw) * Rqc(iso)
             wk_ISO(I_Praut, iso) = wk(I_Praut) * Rqc(iso)
             wk_ISO(I_Pihom, iso) = wk(I_Pihom) * Rqc(iso)
             wk_ISO(I_Pihtr, iso) = wk(I_Pihtr) * Rqc(iso)
             wk_ISO(I_Praci, iso) = wk(I_Praci) * Rqi(iso)
             wk_ISO(I_Psaci, iso) = wk(I_Psaci) * Rqi(iso)
             wk_ISO(I_Pgaci, iso) = wk(I_Pgaci) * Rqi(iso)
             wk_ISO(I_Psaut, iso) = wk(I_Psaut) * Rqi(iso)
             wk_ISO(I_Pisub, iso) = wk(I_Pisub) * Rqi(iso)
             wk_ISO(I_Pimlt, iso) = wk(I_Pimlt) * Rqi(iso)
             wk_ISO(I_Psfi,  iso) = wk(I_Psfi ) * Rqi(iso)
             wk_ISO(I_Piacr, iso) = wk(I_Piacr) * Rqr(iso)
             wk_ISO(I_Psacr, iso) = wk(I_Psacr) * Rqr(iso)
             wk_ISO(I_Pgacr, iso) = wk(I_Pgacr) * Rqr(iso)
             wk_ISO(I_Pgfrz, iso) = wk(I_Pgfrz) * Rqr(iso)
             wk_ISO(I_Pracs, iso) = wk(I_Pracs) * Rqs(iso)
             wk_ISO(I_Pgacs, iso) = wk(I_Pgacs) * Rqs(iso)
             wk_ISO(I_Pgaut, iso) = wk(I_Pgaut) * Rqs(iso)
             wk_ISO(I_Pssub, iso) = wk(I_Pssub) * Rqs(iso)
             wk_ISO(I_Psmlt, iso) = wk(I_Psmlt) * Rqs(iso)
             wk_ISO(I_Pgsub, iso) = wk(I_Pgsub) * Rqg(iso)
             wk_ISO(I_Pgmlt, iso) = wk(I_Pgmlt) * Rqg(iso)
          enddo
       endif
#endif

       !---< limiter >---
       wk(I_Pigen) = min( wk(I_Pigen), wk(I_dqv_dt) ) * (        wk(I_iceflg) ) * sw_expice
       wk(I_Pidep) = min( wk(I_Pidep), wk(I_dqv_dt) ) * (        wk(I_iceflg) ) * sw_expice
       wk(I_Psdep) = min( wk(I_Psdep), wk(I_dqv_dt) ) * (        wk(I_iceflg) )
       wk(I_Pgdep) = min( wk(I_Pgdep), wk(I_dqv_dt) ) * (        wk(I_iceflg) )

       wk(I_Pracw) = wk(I_Pracw)                           &
                   + wk(I_Psacw) * ( 1.0_RP-wk(I_iceflg) ) & ! c->r by s
                   + wk(I_Pgacw) * ( 1.0_RP-wk(I_iceflg) )   ! c->r by g

       wk(I_Praut) = min( wk(I_Praut), wk(I_dqc_dt) )
       wk(I_Pracw) = min( wk(I_Pracw), wk(I_dqc_dt) )
       wk(I_Pihom) = min( wk(I_Pihom), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_expice
       wk(I_Pihtr) = min( wk(I_Pihtr), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_expice
       wk(I_Psacw) = min( wk(I_Psacw), wk(I_dqc_dt) ) * (        wk(I_iceflg) )
       wk(I_Psfw ) = min( wk(I_Psfw ), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_bergeron
       wk(I_Pgacw) = min( wk(I_Pgacw), wk(I_dqc_dt) ) * (        wk(I_iceflg) )

       wk(I_Prevp) = min( wk(I_Prevp), wk(I_dqr_dt) )
       wk(I_Piacr) = min( wk(I_Piacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
       wk(I_Psacr) = min( wk(I_Psacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
       wk(I_Pgacr) = min( wk(I_Pgacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
       wk(I_Pgfrz) = min( wk(I_Pgfrz), wk(I_dqr_dt) ) * (        wk(I_iceflg) )

       wk(I_Pisub) = min( wk(I_Pisub), wk(I_dqi_dt) ) * (        wk(I_iceflg) ) * sw_expice
       wk(I_Pimlt) = min( wk(I_Pimlt), wk(I_dqi_dt) ) * ( 1.0_RP-wk(I_iceflg) ) * sw_expice
       wk(I_Psaut) = min( wk(I_Psaut), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
       wk(I_Praci) = min( wk(I_Praci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
       wk(I_Psaci) = min( wk(I_Psaci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
       wk(I_Psfi ) = min( wk(I_Psfi ), wk(I_dqi_dt) ) * (        wk(I_iceflg) ) * sw_bergeron
       wk(I_Pgaci) = min( wk(I_Pgaci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )

       wk(I_Pssub) = min( wk(I_Pssub), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
       wk(I_Psmlt) = min( wk(I_Psmlt), wk(I_dqs_dt) ) * ( 1.0_RP-wk(I_iceflg) )
       wk(I_Pgaut) = min( wk(I_Pgaut), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
       wk(I_Pracs) = min( wk(I_Pracs), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
       wk(I_Pgacs) = min( wk(I_Pgacs), wk(I_dqs_dt) )

       wk(I_Pgsub) = min( wk(I_Pgsub), wk(I_dqg_dt) ) * (        wk(I_iceflg) )
       wk(I_Pgmlt) = min( wk(I_Pgmlt), wk(I_dqg_dt) ) * ( 1.0_RP-wk(I_iceflg) )

       wk(I_Piacr_s) = ( 1.0_RP - wk(I_delta1) ) * wk(I_Piacr)
       wk(I_Piacr_g) = (          wk(I_delta1) ) * wk(I_Piacr)
       wk(I_Praci_s) = ( 1.0_RP - wk(I_delta1) ) * wk(I_Praci)
       wk(I_Praci_g) = (          wk(I_delta1) ) * wk(I_Praci)
       wk(I_Psacr_s) = (          wk(I_delta2) ) * wk(I_Psacr)
       wk(I_Psacr_g) = ( 1.0_RP - wk(I_delta2) ) * wk(I_Psacr)
       wk(I_Pracs  ) = ( 1.0_RP - wk(I_delta2) ) * wk(I_Pracs)

#ifdef NSW6_WISO
       if ( opt_water_isotope ) then ! [Add] M.Tanoue 20201211
          do iso = 1, NISO
             wk_ISO(I_Pigen,iso) = min( wk_ISO(I_Pigen,iso), wk_ISO(I_dqv_dt,iso) ) * (        wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Pidep,iso) = min( wk_ISO(I_Pidep,iso), wk_ISO(I_dqv_dt,iso) ) * (        wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Psdep,iso) = min( wk_ISO(I_Psdep,iso), wk_ISO(I_dqv_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pgdep,iso) = min( wk_ISO(I_Pgdep,iso), wk_ISO(I_dqv_dt,iso) ) * (        wk(I_iceflg) )

             wk_ISO(I_Pracw,iso) = wk_ISO(I_Pracw,iso)                &
                      + wk_ISO(I_Psacw,iso) * ( 1.0_RP-wk(I_iceflg) ) & ! c->r by s
                      + wk_ISO(I_Pgacw,iso) * ( 1.0_RP-wk(I_iceflg) )   ! c->r by g

             wk_ISO(I_Praut,iso) = min( wk_ISO(I_Praut,iso), wk_ISO(I_dqc_dt,iso) )
             wk_ISO(I_Pracw,iso) = min( wk_ISO(I_Pracw,iso), wk_ISO(I_dqc_dt,iso) )
             wk_ISO(I_Pihom,iso) = min( wk_ISO(I_Pihom,iso), wk_ISO(I_dqc_dt,iso) ) * (        wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Pihtr,iso) = min( wk_ISO(I_Pihtr,iso), wk_ISO(I_dqc_dt,iso) ) * (        wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Psacw,iso) = min( wk_ISO(I_Psacw,iso), wk_ISO(I_dqc_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Psfw ,iso) = min( wk_ISO(I_Psfw ,iso), wk_ISO(I_dqc_dt,iso) ) * (        wk(I_iceflg) ) * sw_bergeron
             wk_ISO(I_Pgacw,iso) = min( wk_ISO(I_Pgacw,iso), wk_ISO(I_dqc_dt,iso) ) * (        wk(I_iceflg) )

             wk_ISO(I_Prevp,iso) = min( wk_ISO(I_Prevp,iso), wk_ISO(I_dqr_dt,iso) )
             wk_ISO(I_Piacr,iso) = min( wk_ISO(I_Piacr,iso), wk_ISO(I_dqr_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Psacr,iso) = min( wk_ISO(I_Psacr,iso), wk_ISO(I_dqr_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pgacr,iso) = min( wk_ISO(I_Pgacr,iso), wk_ISO(I_dqr_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pgfrz,iso) = min( wk_ISO(I_Pgfrz,iso), wk_ISO(I_dqr_dt,iso) ) * (        wk(I_iceflg) )

             wk_ISO(I_Pisub,iso) = min( wk_ISO(I_Pisub,iso), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Pimlt,iso) = min( wk_ISO(I_Pimlt,iso), wk_ISO(I_dqi_dt,iso) ) * ( 1.0_RP-wk(I_iceflg) ) * sw_expice
             wk_ISO(I_Psaut,iso) = min( wk_ISO(I_Psaut,iso), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Praci,iso) = min( wk_ISO(I_Praci,iso), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Psaci,iso) = min( wk_ISO(I_Psaci,iso), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Psfi,iso)  = min( wk_ISO(I_Psfi,iso ), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) ) * sw_bergeron
             wk_ISO(I_Pgaci,iso) = min( wk_ISO(I_Pgaci,iso), wk_ISO(I_dqi_dt,iso) ) * (        wk(I_iceflg) )

             wk_ISO(I_Pssub,iso) = min( wk_ISO(I_Pssub,iso), wk_ISO(I_dqs_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Psmlt,iso) = min( wk_ISO(I_Psmlt,iso), wk_ISO(I_dqs_dt,iso) ) * ( 1.0_RP-wk(I_iceflg) )
             wk_ISO(I_Pgaut,iso) = min( wk_ISO(I_Pgaut,iso), wk_ISO(I_dqs_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pracs,iso) = min( wk_ISO(I_Pracs,iso), wk_ISO(I_dqs_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pgacs,iso) = min( wk_ISO(I_Pgacs,iso), wk_ISO(I_dqs_dt,iso) )

             wk_ISO(I_Pgsub,iso) = min( wk_ISO(I_Pgsub,iso), wk_ISO(I_dqg_dt,iso) ) * (        wk(I_iceflg) )
             wk_ISO(I_Pgmlt,iso) = min( wk_ISO(I_Pgmlt,iso), wk_ISO(I_dqg_dt,iso) ) * ( 1.0_RP-wk(I_iceflg) )

             wk_ISO(I_Piacr_s,iso) = ( 1.0_RP - wk(I_delta1) ) * wk_ISO(I_Piacr,iso)
             wk_ISO(I_Piacr_g,iso) = (          wk(I_delta1) ) * wk_ISO(I_Piacr,iso)
             wk_ISO(I_Praci_s,iso) = ( 1.0_RP - wk(I_delta1) ) * wk_ISO(I_Praci,iso)
             wk_ISO(I_Praci_g,iso) = (          wk(I_delta1) ) * wk_ISO(I_Praci,iso)
             wk_ISO(I_Psacr_s,iso) = (          wk(I_delta2) ) * wk_ISO(I_Psacr,iso)
             wk_ISO(I_Psacr_g,iso) = ( 1.0_RP - wk(I_delta2) ) * wk_ISO(I_Psacr,iso)
             wk_ISO(I_Pracs,iso)   = ( 1.0_RP - wk(I_delta2) ) * wk_ISO(I_Pracs,iso)
          enddo
       endif
#endif

       ! [QC]
       net = &
           + wk(I_Pimlt  ) & ! [prod] i->c
           - wk(I_Praut  ) & ! [loss] c->r
           - wk(I_Pracw  ) & ! [loss] c->r
           - wk(I_Pihom  ) & ! [loss] c->i
           - wk(I_Pihtr  ) & ! [loss] c->i
           - wk(I_Psacw  ) & ! [loss] c->s
           - wk(I_Psfw   ) & ! [loss] c->s
           - wk(I_Pgacw  )   ! [loss] c->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqc_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Pimlt  ) = wk(I_Pimlt  ) * fac
       wk(I_Praut  ) = wk(I_Praut  ) * fac
       wk(I_Pracw  ) = wk(I_Pracw  ) * fac
       wk(I_Pihom  ) = wk(I_Pihom  ) * fac
       wk(I_Pihtr  ) = wk(I_Pihtr  ) * fac
       wk(I_Psacw  ) = wk(I_Psacw  ) * fac
       wk(I_Psfw   ) = wk(I_Psfw   ) * fac
       wk(I_Pgacw  ) = wk(I_Pgacw  ) * fac

       ! [QI]
       net = &
           + wk(I_Pigen  ) & ! [prod] v->i
           + wk(I_Pidep  ) & ! [prod] v->i
           + wk(I_Pihom  ) & ! [prod] c->i
           + wk(I_Pihtr  ) & ! [prod] c->i
           - wk(I_Pisub  ) & ! [loss] i->v
           - wk(I_Pimlt  ) & ! [loss] i->c
           - wk(I_Psaut  ) & ! [loss] i->s
           - wk(I_Praci_s) & ! [loss] i->s
           - wk(I_Psaci  ) & ! [loss] i->s
           - wk(I_Psfi   ) & ! [loss] i->s
           - wk(I_Praci_g) & ! [loss] i->g
           - wk(I_Pgaci  )   ! [loss] i->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqi_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Pigen  ) = wk(I_Pigen  ) * fac
       wk(I_Pidep  ) = wk(I_Pidep  ) * fac
       wk(I_Pihom  ) = wk(I_Pihom  ) * fac
       wk(I_Pihtr  ) = wk(I_Pihtr  ) * fac
       wk(I_Pisub  ) = wk(I_Pisub  ) * fac
       wk(I_Pimlt  ) = wk(I_Pimlt  ) * fac
       wk(I_Psaut  ) = wk(I_Psaut  ) * fac
       wk(I_Praci_s) = wk(I_Praci_s) * fac
       wk(I_Psaci  ) = wk(I_Psaci  ) * fac
       wk(I_Psfi   ) = wk(I_Psfi   ) * fac
       wk(I_Praci_g) = wk(I_Praci_g) * fac
       wk(I_Pgaci  ) = wk(I_Pgaci  ) * fac

       ! [QR]
       net = &
           + wk(I_Praut  ) & ! [prod] c->r
           + wk(I_Pracw  ) & ! [prod] c->r
           + wk(I_Psmlt  ) & ! [prod] s->r
           + wk(I_Pgmlt  ) & ! [prod] g->r
           - wk(I_Prevp  ) & ! [loss] r->v
           - wk(I_Piacr_s) & ! [loss] r->s
           - wk(I_Psacr_s) & ! [loss] r->s
           - wk(I_Piacr_g) & ! [loss] r->g
           - wk(I_Psacr_g) & ! [loss] r->g
           - wk(I_Pgacr  ) & ! [loss] r->g
           - wk(I_Pgfrz  )   ! [loss] r->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqr_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Praut  ) = wk(I_Praut  ) * fac
       wk(I_Pracw  ) = wk(I_Pracw  ) * fac
       wk(I_Psmlt  ) = wk(I_Psmlt  ) * fac
       wk(I_Pgmlt  ) = wk(I_Pgmlt  ) * fac
       wk(I_Prevp  ) = wk(I_Prevp  ) * fac
       wk(I_Piacr_s) = wk(I_Piacr_s) * fac
       wk(I_Psacr_s) = wk(I_Psacr_s) * fac
       wk(I_Piacr_g) = wk(I_Piacr_g) * fac
       wk(I_Psacr_g) = wk(I_Psacr_g) * fac
       wk(I_Pgacr  ) = wk(I_Pgacr  ) * fac
       wk(I_Pgfrz  ) = wk(I_Pgfrz  ) * fac

       ! [QV]
       net = &
           + wk(I_Prevp  ) & ! [prod] r->v
           + wk(I_Pisub  ) & ! [prod] i->v
           + wk(I_Pssub  ) & ! [prod] s->v
           + wk(I_Pgsub  ) & ! [prod] g->v
           - wk(I_Pigen  ) & ! [loss] v->i
           - wk(I_Pidep  ) & ! [loss] v->i
           - wk(I_Psdep  ) & ! [loss] v->s
           - wk(I_Pgdep  )   ! [loss] v->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqv_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Prevp  ) = wk(I_Prevp  ) * fac
       wk(I_Pisub  ) = wk(I_Pisub  ) * fac
       wk(I_Pssub  ) = wk(I_Pssub  ) * fac
       wk(I_Pgsub  ) = wk(I_Pgsub  ) * fac
       wk(I_Pigen  ) = wk(I_Pigen  ) * fac
       wk(I_Pidep  ) = wk(I_Pidep  ) * fac
       wk(I_Psdep  ) = wk(I_Psdep  ) * fac
       wk(I_Pgdep  ) = wk(I_Pgdep  ) * fac

       ! [QS]
       net = &
           + wk(I_Psdep  ) & ! [prod] v->s
           + wk(I_Psacw  ) & ! [prod] c->s
           + wk(I_Psfw   ) & ! [prod] c->s
           + wk(I_Piacr_s) & ! [prod] r->s
           + wk(I_Psacr_s) & ! [prod] r->s
           + wk(I_Psaut  ) & ! [prod] i->s
           + wk(I_Praci_s) & ! [prod] i->s
           + wk(I_Psaci  ) & ! [prod] i->s
           + wk(I_Psfi   ) & ! [prod] i->s
           - wk(I_Pssub  ) & ! [loss] s->v
           - wk(I_Psmlt  ) & ! [loss] s->r
           - wk(I_Pgaut  ) & ! [loss] s->g
           - wk(I_Pracs  ) & ! [loss] s->g
           - wk(I_Pgacs  )   ! [loss] s->g

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqs_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Psdep  ) = wk(I_Psdep  ) * fac
       wk(I_Psacw  ) = wk(I_Psacw  ) * fac
       wk(I_Psfw   ) = wk(I_Psfw   ) * fac
       wk(I_Piacr_s) = wk(I_Piacr_s) * fac
       wk(I_Psacr_s) = wk(I_Psacr_s) * fac
       wk(I_Psaut  ) = wk(I_Psaut  ) * fac
       wk(I_Praci_s) = wk(I_Praci_s) * fac
       wk(I_Psaci  ) = wk(I_Psaci  ) * fac
       wk(I_Psfi   ) = wk(I_Psfi   ) * fac
       wk(I_Pssub  ) = wk(I_Pssub  ) * fac
       wk(I_Psmlt  ) = wk(I_Psmlt  ) * fac
       wk(I_Pgaut  ) = wk(I_Pgaut  ) * fac
       wk(I_Pracs  ) = wk(I_Pracs  ) * fac
       wk(I_Pgacs  ) = wk(I_Pgacs  ) * fac

       ! [QG]
       net = &
           + wk(I_Pgdep  ) & ! [prod] v->g
           + wk(I_Pgacw  ) & ! [prod] c->g
           + wk(I_Piacr_g) & ! [prod] r->g
           + wk(I_Psacr_g) & ! [prod] r->g
           + wk(I_Pgacr  ) & ! [prod] r->g
           + wk(I_Pgfrz  ) & ! [prod] r->g
           + wk(I_Praci_g) & ! [prod] i->g
           + wk(I_Pgaci  ) & ! [prod] i->g
           + wk(I_Pgaut  ) & ! [prod] s->g
           + wk(I_Pracs  ) & ! [prod] s->g
           + wk(I_Pgacs  ) & ! [prod] s->g
           - wk(I_Pgsub  ) & ! [loss] g->v
           - wk(I_Pgmlt  )   ! [loss] g->r

       fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
       fac    = (          fac_sw ) &
              + ( 1.0_RP - fac_sw ) * min( -wk(I_dqg_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

       wk(I_Pgdep  ) = wk(I_Pgdep  ) * fac
       wk(I_Pgacw  ) = wk(I_Pgacw  ) * fac
       wk(I_Piacr_g) = wk(I_Piacr_g) * fac
       wk(I_Psacr_g) = wk(I_Psacr_g) * fac
       wk(I_Pgacr  ) = wk(I_Pgacr  ) * fac
       wk(I_Pgfrz  ) = wk(I_Pgfrz  ) * fac
       wk(I_Praci_g) = wk(I_Praci_g) * fac
       wk(I_Pgaci  ) = wk(I_Pgaci  ) * fac
       wk(I_Pgaut  ) = wk(I_Pgaut  ) * fac
       wk(I_Pracs  ) = wk(I_Pracs  ) * fac
       wk(I_Pgacs  ) = wk(I_Pgacs  ) * fac
       wk(I_Pgsub  ) = wk(I_Pgsub  ) * fac
       wk(I_Pgmlt  ) = wk(I_Pgmlt  ) * fac

       qc_t = + wk(I_Pimlt  ) & ! [prod] i->c
              - wk(I_Praut  ) & ! [loss] c->r
              - wk(I_Pracw  ) & ! [loss] c->r
              - wk(I_Pihom  ) & ! [loss] c->i
              - wk(I_Pihtr  ) & ! [loss] c->i
              - wk(I_Psacw  ) & ! [loss] c->s
              - wk(I_Psfw   ) & ! [loss] c->s
              - wk(I_Pgacw  )   ! [loss] c->g

       qr_t = + wk(I_Praut  ) & ! [prod] c->r
              + wk(I_Pracw  ) & ! [prod] c->r
              + wk(I_Psmlt  ) & ! [prod] s->r
              + wk(I_Pgmlt  ) & ! [prod] g->r
              - wk(I_Prevp  ) & ! [loss] r->v
              - wk(I_Piacr_s) & ! [loss] r->s
              - wk(I_Psacr_s) & ! [loss] r->s
              - wk(I_Piacr_g) & ! [loss] r->g
              - wk(I_Psacr_g) & ! [loss] r->g
              - wk(I_Pgacr  ) & ! [loss] r->g
              - wk(I_Pgfrz  )   ! [loss] r->g

       qi_t = + wk(I_Pigen  ) & ! [prod] v->i
              + wk(I_Pidep  ) & ! [prod] v->i
              + wk(I_Pihom  ) & ! [prod] c->i
              + wk(I_Pihtr  ) & ! [prod] c->i
              - wk(I_Pisub  ) & ! [loss] i->v
              - wk(I_Pimlt  ) & ! [loss] i->c
              - wk(I_Psaut  ) & ! [loss] i->s
              - wk(I_Praci_s) & ! [loss] i->s
              - wk(I_Psaci  ) & ! [loss] i->s
              - wk(I_Psfi   ) & ! [loss] i->s
              - wk(I_Praci_g) & ! [loss] i->g
              - wk(I_Pgaci  )   ! [loss] i->g

       qs_t = + wk(I_Psdep  ) & ! [prod] v->s
              + wk(I_Psacw  ) & ! [prod] c->s
              + wk(I_Psfw   ) & ! [prod] c->s
              + wk(I_Piacr_s) & ! [prod] r->s
              + wk(I_Psacr_s) & ! [prod] r->s
              + wk(I_Psaut  ) & ! [prod] i->s
              + wk(I_Praci_s) & ! [prod] i->s
              + wk(I_Psaci  ) & ! [prod] i->s
              + wk(I_Psfi   ) & ! [prod] i->s
              - wk(I_Pssub  ) & ! [loss] s->v
              - wk(I_Psmlt  ) & ! [loss] s->r
              - wk(I_Pgaut  ) & ! [loss] s->g
              - wk(I_Pracs  ) & ! [loss] s->g
              - wk(I_Pgacs  )   ! [loss] s->g

       qg_t = + wk(I_Pgdep  ) & ! [prod] v->g
              + wk(I_Pgacw  ) & ! [prod] c->g
              + wk(I_Piacr_g) & ! [prod] r->g
              + wk(I_Psacr_g) & ! [prod] r->g
              + wk(I_Pgacr  ) & ! [prod] r->g
              + wk(I_Pgfrz  ) & ! [prod] r->g
              + wk(I_Praci_g) & ! [prod] i->g
              + wk(I_Pgaci  ) & ! [prod] i->g
              + wk(I_Pgaut  ) & ! [prod] s->g
              + wk(I_Pracs  ) & ! [prod] s->g
              + wk(I_Pgacs  ) & ! [prod] s->g
              - wk(I_Pgsub  ) & ! [loss] g->v
              - wk(I_Pgmlt  )   ! [loss] g->r

       qc_t = max( qc_t, -wk(I_dqc_dt) )
       qr_t = max( qr_t, -wk(I_dqr_dt) )
       qi_t = max( qi_t, -wk(I_dqi_dt) )
       qs_t = max( qs_t, -wk(I_dqs_dt) )
       qg_t = max( qg_t, -wk(I_dqg_dt) )

       qv_t = - ( qc_t &
                + qr_t &
                + qi_t &
                + qs_t &
                + qg_t )

       drhogqv(ij,k) = rhog(ij,k) * qv_t * dt
       drhogqc(ij,k) = rhog(ij,k) * qc_t * dt
       drhogqr(ij,k) = rhog(ij,k) * qr_t * dt
       drhogqi(ij,k) = rhog(ij,k) * qi_t * dt
       drhogqs(ij,k) = rhog(ij,k) * qs_t * dt
       drhogqg(ij,k) = rhog(ij,k) * qg_t * dt

       Vt(ij,k,I_QR) = Vtr
       Vt(ij,k,I_QI) = Vti
       Vt(ij,k,I_QS) = Vts
       Vt(ij,k,I_QG) = Vtg

       evap_liq(ij,k) = wk(I_Prevp)                             ! 2018/11/08 Yamashita
       evap_ice(ij,k) = wk(I_Pisub) + wk(I_Pssub) + wk(I_Pgsub) ! 2018/11/08 Yamashita

#ifdef NSW6_WISO
       if ( opt_water_isotope ) then ! [Add] M.Tanoue 20201211
          do iso = 1, NISO
             ! [QC]
             net = &
              + wk_ISO(I_Pimlt, iso  ) & ! [prod] i->c
              - wk_ISO(I_Praut, iso  ) & ! [loss] c->r
              - wk_ISO(I_Pracw, iso  ) & ! [loss] c->r
              - wk_ISO(I_Pihom, iso  ) & ! [loss] c->i
              - wk_ISO(I_Pihtr, iso  ) & ! [loss] c->i
              - wk_ISO(I_Psacw, iso  ) & ! [loss] c->s
              - wk_ISO(I_Psfw , iso  ) & ! [loss] c->s
              - wk_ISO(I_Pgacw, iso  )   ! [loss] c->g

             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqc_dt,iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Pimlt,iso) = wk_ISO(I_Pimlt,iso) * fac
             wk_ISO(I_Praut,iso) = wk_ISO(I_Praut,iso) * fac
             wk_ISO(I_Pracw,iso) = wk_ISO(I_Pracw,iso) * fac
             wk_ISO(I_Pihom,iso) = wk_ISO(I_Pihom,iso) * fac
             wk_ISO(I_Pihtr,iso) = wk_ISO(I_Pihtr,iso) * fac
             wk_ISO(I_Psacw,iso) = wk_ISO(I_Psacw,iso) * fac
             wk_ISO(I_Psfw,iso ) = wk_ISO(I_Psfw,iso ) * fac
             wk_ISO(I_Pgacw,iso) = wk_ISO(I_Pgacw,iso) * fac

             ! [QI]
             net = &
              + wk_ISO(I_Pigen  ,iso) & ! [prod] v->i
              + wk_ISO(I_Pidep  ,iso) & ! [prod] v->i
              + wk_ISO(I_Pihom  ,iso) & ! [prod] c->i
              + wk_ISO(I_Pihtr  ,iso) & ! [prod] c->i
              - wk_ISO(I_Pisub  ,iso) & ! [loss] i->v
              - wk_ISO(I_Pimlt  ,iso) & ! [loss] i->c
              - wk_ISO(I_Psaut  ,iso) & ! [loss] i->s
              - wk_ISO(I_Praci_s,iso) & ! [loss] i->s
              - wk_ISO(I_Psaci  ,iso) & ! [loss] i->s
              - wk_ISO(I_Psfi   ,iso) & ! [loss] i->s
              - wk_ISO(I_Praci_g,iso) & ! [loss] i->g
              - wk_ISO(I_Pgaci  ,iso)   ! [loss] i->g
             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac   = (          fac_sw ) &
                    + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqi_dt, iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Pigen,iso  ) = wk_ISO(I_Pigen,iso  ) * fac
             wk_ISO(I_Pidep,iso  ) = wk_ISO(I_Pidep,iso  ) * fac
             wk_ISO(I_Pihom,iso  ) = wk_ISO(I_Pihom,iso  ) * fac
             wk_ISO(I_Pihtr,iso  ) = wk_ISO(I_Pihtr,iso  ) * fac
             wk_ISO(I_Pisub,iso  ) = wk_ISO(I_Pisub,iso  ) * fac
             wk_ISO(I_Pimlt,iso  ) = wk_ISO(I_Pimlt,iso  ) * fac
             wk_ISO(I_Psaut,iso  ) = wk_ISO(I_Psaut,iso  ) * fac
             wk_ISO(I_Praci_s,iso) = wk_ISO(I_Praci_s,iso) * fac
             wk_ISO(I_Psaci,iso  ) = wk_ISO(I_Psaci,iso  ) * fac
             wk_ISO(I_Psfi,iso   ) = wk_ISO(I_Psfi,iso   ) * fac
             wk_ISO(I_Praci_g,iso) = wk_ISO(I_Praci_g,iso) * fac
             wk_ISO(I_Pgaci,iso  ) = wk_ISO(I_Pgaci,iso  ) * fac

             ! [QR]
             net = &
              + wk_ISO(I_Praut  ,iso) & ! [prod] c->r
              + wk_ISO(I_Pracw  ,iso) & ! [prod] c->r
              + wk_ISO(I_Psmlt  ,iso) & ! [prod] s->r
              + wk_ISO(I_Pgmlt  ,iso) & ! [prod] g->r
              - wk_ISO(I_Prevp  ,iso) & ! [loss] r->v
              - wk_ISO(I_Piacr_s,iso) & ! [loss] r->s
              - wk_ISO(I_Psacr_s,iso) & ! [loss] r->s
              - wk_ISO(I_Piacr_g,iso) & ! [loss] r->g
              - wk_ISO(I_Psacr_g,iso) & ! [loss] r->g
              - wk_ISO(I_Pgacr  ,iso) & ! [loss] r->g
              - wk_ISO(I_Pgfrz  ,iso)   ! [loss] r->g
             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac   = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqr_dt,iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Praut,iso  ) = wk_ISO(I_Praut,iso  ) * fac
             wk_ISO(I_Pracw,iso  ) = wk_ISO(I_Pracw,iso  ) * fac
             wk_ISO(I_Psmlt,iso  ) = wk_ISO(I_Psmlt,iso  ) * fac
             wk_ISO(I_Pgmlt,iso  ) = wk_ISO(I_Pgmlt,iso  ) * fac
             wk_ISO(I_Prevp,iso  ) = wk_ISO(I_Prevp,iso  ) * fac
             wk_ISO(I_Piacr_s,iso) = wk_ISO(I_Piacr_s,iso) * fac
             wk_ISO(I_Psacr_s,iso) = wk_ISO(I_Psacr_s,iso) * fac
             wk_ISO(I_Piacr_g,iso) = wk_ISO(I_Piacr_g,iso) * fac
             wk_ISO(I_Psacr_g,iso) = wk_ISO(I_Psacr_g,iso) * fac
             wk_ISO(I_Pgacr,iso  ) = wk_ISO(I_Pgacr,iso  ) * fac
             wk_ISO(I_Pgfrz,iso  ) = wk_ISO(I_Pgfrz,iso  ) * fac

             ! [QV]
             net = &
              + wk_ISO(I_Prevp  ,iso) & ! [prod] r->v
              + wk_ISO(I_Pisub  ,iso) & ! [prod] i->v
              + wk_ISO(I_Pssub  ,iso) & ! [prod] s->v
              + wk_ISO(I_Pgsub  ,iso) & ! [prod] g->v
              - wk_ISO(I_Pigen  ,iso) & ! [loss] v->i
              - wk_ISO(I_Pidep  ,iso) & ! [loss] v->i
              - wk_ISO(I_Psdep  ,iso) & ! [loss] v->s
              - wk_ISO(I_Pgdep  ,iso)   ! [loss] v->g
             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac    = (          fac_sw ) &
                     + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqv_dt,iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Prevp,iso  ) = wk_ISO(I_Prevp,iso  ) * fac
             wk_ISO(I_Pisub,iso  ) = wk_ISO(I_Pisub,iso  ) * fac
             wk_ISO(I_Pssub,iso  ) = wk_ISO(I_Pssub,iso  ) * fac
             wk_ISO(I_Pgsub,iso  ) = wk_ISO(I_Pgsub,iso  ) * fac
             wk_ISO(I_Pigen,iso  ) = wk_ISO(I_Pigen,iso  ) * fac
             wk_ISO(I_Pidep,iso  ) = wk_ISO(I_Pidep,iso  ) * fac
             wk_ISO(I_Psdep,iso  ) = wk_ISO(I_Psdep,iso  ) * fac
             wk_ISO(I_Pgdep,iso  ) = wk_ISO(I_Pgdep,iso  ) * fac

             ! [QS]
             net = &
              + wk_ISO(I_Psdep  ,iso) & ! [prod] v->s
              + wk_ISO(I_Psacw  ,iso) & ! [prod] c->s
              + wk_ISO(I_Psfw   ,iso) & ! [prod] c->s
              + wk_ISO(I_Piacr_s,iso) & ! [prod] r->s
              + wk_ISO(I_Psacr_s,iso) & ! [prod] r->s
              + wk_ISO(I_Psaut  ,iso) & ! [prod] i->s
              + wk_ISO(I_Praci_s,iso) & ! [prod] i->s
              + wk_ISO(I_Psaci  ,iso) & ! [prod] i->s
              + wk_ISO(I_Psfi   ,iso) & ! [prod] i->s
              - wk_ISO(I_Pssub  ,iso) & ! [loss] s->v
              - wk_ISO(I_Psmlt  ,iso) & ! [loss] s->r
              - wk_ISO(I_Pgaut  ,iso) & ! [loss] s->g
              - wk_ISO(I_Pracs  ,iso) & ! [loss] s->g
              - wk_ISO(I_Pgacs  ,iso)   ! [loss] s->g
             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac   = (          fac_sw ) &
                    + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqs_dt, iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Psdep,iso  ) = wk_ISO(I_Psdep,iso  ) * fac
             wk_ISO(I_Psacw,iso  ) = wk_ISO(I_Psacw,iso  ) * fac
             wk_ISO(I_Psfw,iso   ) = wk_ISO(I_Psfw,iso   ) * fac
             wk_ISO(I_Piacr_s,iso) = wk_ISO(I_Piacr_s,iso) * fac
             wk_ISO(I_Psacr_s,iso) = wk_ISO(I_Psacr_s,iso) * fac
             wk_ISO(I_Psaut,iso  ) = wk_ISO(I_Psaut,iso  ) * fac
             wk_ISO(I_Praci_s,iso) = wk_ISO(I_Praci_s,iso) * fac
             wk_ISO(I_Psaci,iso  ) = wk_ISO(I_Psaci,iso  ) * fac
             wk_ISO(I_Psfi,iso   ) = wk_ISO(I_Psfi,iso   ) * fac
             wk_ISO(I_Pssub,iso  ) = wk_ISO(I_Pssub,iso  ) * fac
             wk_ISO(I_Psmlt,iso  ) = wk_ISO(I_Psmlt,iso  ) * fac
             wk_ISO(I_Pgaut,iso  ) = wk_ISO(I_Pgaut,iso  ) * fac
             wk_ISO(I_Pracs,iso  ) = wk_ISO(I_Pracs,iso  ) * fac
             wk_ISO(I_Pgacs,iso  ) = wk_ISO(I_Pgacs,iso  ) * fac

             ! [QG]
             net = &
              + wk_ISO(I_Pgdep  ,iso) & ! [prod] v->g
              + wk_ISO(I_Pgacw  ,iso) & ! [prod] c->g
              + wk_ISO(I_Piacr_g,iso) & ! [prod] r->g
              + wk_ISO(I_Psacr_g,iso) & ! [prod] r->g
              + wk_ISO(I_Pgacr  ,iso) & ! [prod] r->g
              + wk_ISO(I_Pgfrz  ,iso) & ! [prod] r->g
              + wk_ISO(I_Praci_g,iso) & ! [prod] i->g
              + wk_ISO(I_Pgaci  ,iso) & ! [prod] i->g
              + wk_ISO(I_Pgaut  ,iso) & ! [prod] s->g
              + wk_ISO(I_Pracs  ,iso) & ! [prod] s->g
              + wk_ISO(I_Pgacs  ,iso) & ! [prod] s->g
              - wk_ISO(I_Pgsub  ,iso) & ! [loss] g->v
              - wk_ISO(I_Pgmlt  ,iso)   ! [loss] g->r
             fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
             fac   = (          fac_sw ) &
                    + ( 1.0_RP - fac_sw ) * min( -wk_ISO(I_dqg_dt, iso)/(net-fac_sw), 1.0_RP ) ! loss limiter
             wk_ISO(I_Pgdep,iso  ) = wk_ISO(I_Pgdep,iso  ) * fac
             wk_ISO(I_Pgacw,iso  ) = wk_ISO(I_Pgacw,iso  ) * fac
             wk_ISO(I_Piacr_g,iso) = wk_ISO(I_Piacr_g,iso) * fac
             wk_ISO(I_Psacr_g,iso) = wk_ISO(I_Psacr_g,iso) * fac
             wk_ISO(I_Pgacr,iso  ) = wk_ISO(I_Pgacr,iso  ) * fac
             wk_ISO(I_Pgfrz,iso  ) = wk_ISO(I_Pgfrz,iso  ) * fac
             wk_ISO(I_Praci_g,iso) = wk_ISO(I_Praci_g,iso) * fac
             wk_ISO(I_Pgaci,iso  ) = wk_ISO(I_Pgaci,iso  ) * fac
             wk_ISO(I_Pgaut,iso  ) = wk_ISO(I_Pgaut,iso  ) * fac
             wk_ISO(I_Pracs,iso  ) = wk_ISO(I_Pracs,iso  ) * fac
             wk_ISO(I_Pgacs,iso  ) = wk_ISO(I_Pgacs,iso  ) * fac
             wk_ISO(I_Pgsub,iso  ) = wk_ISO(I_Pgsub,iso  ) * fac
             wk_ISO(I_Pgmlt,iso  ) = wk_ISO(I_Pgmlt,iso  ) * fac

             qci_t(iso) =                  &
                   + wk_ISO(I_Pimlt,iso  ) & ! [prod] i->c
                   - wk_ISO(I_Praut,iso  ) & ! [loss] c->r
                   - wk_ISO(I_Pracw,iso  ) & ! [loss] c->r
                   - wk_ISO(I_Pihom,iso  ) & ! [loss] c->i
                   - wk_ISO(I_Pihtr,iso  ) & ! [loss] c->i
                   - wk_ISO(I_Psacw,iso  ) & ! [loss] c->s
                   - wk_ISO(I_Psfw,iso   ) & ! [loss] c->s
                   - wk_ISO(I_Pgacw,iso  )   ! [loss] c->g

             qri_t(iso) =                  &
                   + wk_ISO(I_Praut,iso  ) & ! [prod] c->r
                   + wk_ISO(I_Pracw,iso  ) & ! [prod] c->r
                   + wk_ISO(I_Psmlt,iso  ) & ! [prod] s->r
                   + wk_ISO(I_Pgmlt,iso  ) & ! [prod] g->r
                   - wk_ISO(I_Prevp,iso  ) & ! [loss] r->v
                   - wk_ISO(I_Piacr_s,iso) & ! [loss] r->s
                   - wk_ISO(I_Psacr_s,iso) & ! [loss] r->s
                   - wk_ISO(I_Piacr_g,iso) & ! [loss] r->g
                   - wk_ISO(I_Psacr_g,iso) & ! [loss] r->g
                   - wk_ISO(I_Pgacr,iso  ) & ! [loss] r->g
                   - wk_ISO(I_Pgfrz,iso  )   ! [loss] r->g

             qii_t(iso) =                  &
                   + wk_ISO(I_Pigen,iso  ) & ! [prod] v->i
                   + wk_ISO(I_Pidep,iso  ) & ! [prod] v->i
                   + wk_ISO(I_Pihom,iso  ) & ! [prod] c->i
                   + wk_ISO(I_Pihtr,iso  ) & ! [prod] c->i
                   - wk_ISO(I_Pisub,iso  ) & ! [loss] i->v
                   - wk_ISO(I_Pimlt,iso  ) & ! [loss] i->c
                   - wk_ISO(I_Psaut,iso  ) & ! [loss] i->s
                   - wk_ISO(I_Praci_s,iso) & ! [loss] i->s
                   - wk_ISO(I_Psaci,iso  ) & ! [loss] i->s
                   - wk_ISO(I_Psfi,iso   ) & ! [loss] i->s
                   - wk_ISO(I_Praci_g,iso) & ! [loss] i->g
                   - wk_ISO(I_Pgaci,iso  )   ! [loss] i->g

             qsi_t(iso) =                  &
                   + wk_ISO(I_Psdep,iso  ) & ! [prod] v->s
                   + wk_ISO(I_Psacw,iso  ) & ! [prod] c->s
                   + wk_ISO(I_Psfw,iso   ) & ! [prod] c->s
                   + wk_ISO(I_Piacr_s,iso) & ! [prod] r->s
                   + wk_ISO(I_Psacr_s,iso) & ! [prod] r->s
                   + wk_ISO(I_Psaut,iso  ) & ! [prod] i->s
                   + wk_ISO(I_Praci_s,iso) & ! [prod] i->s
                   + wk_ISO(I_Psaci,iso  ) & ! [prod] i->s
                   + wk_ISO(I_Psfi,iso   ) & ! [prod] i->s
                   - wk_ISO(I_Pssub,iso  ) & ! [loss] s->v
                   - wk_ISO(I_Psmlt,iso  ) & ! [loss] s->r
                   - wk_ISO(I_Pgaut,iso  ) & ! [loss] s->g
                   - wk_ISO(I_Pracs,iso  ) & ! [loss] s->g
                   - wk_ISO(I_Pgacs,iso  )   ! [loss] s->g

             qgi_t(iso) =                  &
                   + wk_ISO(I_Pgdep,iso  ) & ! [prod] v->g
                   + wk_ISO(I_Pgacw,iso  ) & ! [prod] c->g
                   + wk_ISO(I_Piacr_g,iso) & ! [prod] r->g
                   + wk_ISO(I_Psacr_g,iso) & ! [prod] r->g
                   + wk_ISO(I_Pgacr,iso  ) & ! [prod] r->g
                   + wk_ISO(I_Pgfrz,iso  ) & ! [prod] r->g
                   + wk_ISO(I_Praci_g,iso) & ! [prod] i->g
                   + wk_ISO(I_Pgaci,iso  ) & ! [prod] i->g
                   + wk_ISO(I_Pgaut,iso  ) & ! [prod] s->g
                   + wk_ISO(I_Pracs,iso  ) & ! [prod] s->g
                   + wk_ISO(I_Pgacs,iso  ) & ! [prod] s->g
                   - wk_ISO(I_Pgsub,iso  ) & ! [loss] g->v
                   - wk_ISO(I_Pgmlt,iso  )   ! [loss] g->r

             qci_t(iso) = max( qci_t(iso), -wk_ISO(I_dqc_dt,iso) )
             qri_t(iso) = max( qri_t(iso), -wk_ISO(I_dqr_dt,iso) )
             qii_t(iso) = max( qii_t(iso), -wk_ISO(I_dqi_dt,iso) )
             qsi_t(iso) = max( qsi_t(iso), -wk_ISO(I_dqs_dt,iso) )
             qgi_t(iso) = max( qgi_t(iso), -wk_ISO(I_dqg_dt,iso) )

             qvi_t(iso) = - ( qci_t(iso) &
                   + qri_t(iso) &
                   + qii_t(iso) &
                   + qsi_t(iso) &
                   + qgi_t(iso) )

             drhogqvi(ij,k,iso) = rhog(ij,k) * qvi_t(iso) * dt
             drhogqci(ij,k,iso) = rhog(ij,k) * qci_t(iso) * dt
             drhogqri(ij,k,iso) = rhog(ij,k) * qri_t(iso) * dt
             drhogqii(ij,k,iso) = rhog(ij,k) * qii_t(iso) * dt
             drhogqsi(ij,k,iso) = rhog(ij,k) * qsi_t(iso) * dt
             drhogqgi(ij,k,iso) = rhog(ij,k) * qgi_t(iso) * dt

             Vt(ij,k,I_QRI(iso)) = Vtr
             Vt(ij,k,I_QII(iso)) = Vti
             Vt(ij,k,I_QSI(iso)) = Vts
             Vt(ij,k,I_QGI(iso)) = Vtg
          enddo
       endif
#endif

       !--- Individual tendency for each process
!        do ip = 1, wk_nmax
!           ml_wk(ij,k,ip) = wk(ip)
!        enddo

       enddo ! vector length loop

    enddo ! cache blocking
    enddo ! k loop
    !$omp end parallel do
    !$acc end kernels

    call PROF_rapend  ('____MP_NSW6_CORE',3)

#ifdef MORETIMER
call PROF_rapstart('____MP_NSW6_06')
#endif

    !--- update rhogq
    !$acc kernels async(1) pcopy(rhogq) pcopyin(drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhogq,drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg, &
    !$omp        I_QV,I_QC,I_QR,I_QI,I_QS,I_QG )
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) + drhogqv(ij,k)
       rhogq(ij,k,I_QC) = rhogq(ij,k,I_QC) + drhogqc(ij,k)
       rhogq(ij,k,I_QR) = rhogq(ij,k,I_QR) + drhogqr(ij,k)
       rhogq(ij,k,I_QI) = rhogq(ij,k,I_QI) + drhogqi(ij,k)
       rhogq(ij,k,I_QS) = rhogq(ij,k,I_QS) + drhogqs(ij,k)
       rhogq(ij,k,I_QG) = rhogq(ij,k,I_QG) + drhogqg(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef NSW6_WISO
    if ( opt_water_isotope )then ! [Add] M.Tanoue 20201211
       !$acc kernels async(1) pcopy(rhogq) pcopyin(drhogqvi,drhogqci,drhogqri,drhogqii,drhogqsi,drhogqgi)
       !$omp parallel do default(none),private(ij,k,iso), &
       !$omp shared(ijdim,kmin,kmax,rhogq,drhogqvi,drhogqci,drhogqri,drhogqii,drhogqsi,drhogqgi, &
       !$omp        I_QV,I_QC,I_QR,I_QI,I_QS,I_QG,I_QVI,I_QCI,I_QRI,I_QII,I_QSI,I_QGI            )
       do k   = kmin, kmax
       do ij  = 1, ijdim
          if ( rhogq(ij,k,I_QV) == 0.0_DP ) then
             rhogq(ij,k,I_QVI(1)) = 0.0_DP
             rhogq(ij,k,I_QVI(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QVI(iso)) = rhogq(ij,k,I_QVI(iso)) + drhogqvi(ij,k,iso)
                rhogq(ij,k,I_QVI(iso)) = max( rhogq(ij,k,I_QVI(iso)), 0.0_DP )
             enddo
          endif
          if ( rhogq(ij,k,I_QC) == 0.0_DP ) then
             rhogq(ij,k,I_QCI(1)) = 0.0_DP
             rhogq(ij,k,I_QCI(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QCI(iso)) = rhogq(ij,k,I_QCI(iso)) + drhogqci(ij,k,iso)
                rhogq(ij,k,I_QCI(iso)) = max( rhogq(ij,k,I_QCI(iso)), 0.0_DP )
             enddo
          endif
          if ( rhogq(ij,k,I_QR) == 0.0_DP ) then
             rhogq(ij,k,I_QRI(1)) = 0.0_DP
             rhogq(ij,k,I_QRI(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QRI(iso)) = rhogq(ij,k,I_QRI(iso)) + drhogqri(ij,k,iso)
                rhogq(ij,k,I_QRI(iso)) = max( rhogq(ij,k,I_QRI(iso)), 0.0_DP )
             enddo
          endif
          if ( rhogq(ij,k,I_QI) == 0.0_DP ) then
             rhogq(ij,k,I_QII(1)) = 0.0_DP
             rhogq(ij,k,I_QII(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QII(iso)) = rhogq(ij,k,I_QII(iso)) + drhogqii(ij,k,iso)
                rhogq(ij,k,I_QII(iso)) = max( rhogq(ij,k,I_QII(iso)), 0.0_DP )
             enddo
          endif
          if ( rhogq(ij,k,I_QS) == 0.0_DP ) then
             rhogq(ij,k,I_QSI(1)) = 0.0_DP
             rhogq(ij,k,I_QSI(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QSI(iso)) = rhogq(ij,k,I_QSI(iso)) + drhogqsi(ij,k,iso)
                rhogq(ij,k,I_QSI(iso)) = max( rhogq(ij,k,I_QSI(iso)), 0.0_DP )
             enddo
          endif
          if ( rhogq(ij,k,I_QG) == 0.0_DP ) then
             rhogq(ij,k,I_QGI(1)) = 0.0_DP
             rhogq(ij,k,I_QGI(2)) = 0.0_DP
          else
             do iso = 1, NISO
                rhogq(ij,k,I_QGI(iso)) = rhogq(ij,k,I_QGI(iso)) + drhogqgi(ij,k,iso)
                rhogq(ij,k,I_QGI(iso)) = max( rhogq(ij,k,I_QGI(iso)), 0.0_DP )
             enddo
          endif
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif
#endif

    !--- update rhoge
    !$acc kernels async(1) pcopy(rhoge) pcopyin(drhogqv,drhogqi,drhogqs,drhogqg)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhoge,drhogqv,drhogqi,drhogqs,drhogqg,LHV,LHF)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhoge(ij,k) = rhoge(ij,k) - LHV * drhogqv(ij,k) &
                                 + LHF * drhogqi(ij,k) &
                                 + LHF * drhogqs(ij,k) &
                                 + LHF * drhogqg(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !--- Effective Radius of Liquid Water

    ! dgamma_a is parameter of Gamma-Dist.
    ! 'Berry and Reinhardt(1974-a) eq.(16-a,b)'
    ! 0.38 * {var1/2x} =~ {var1/2 r}
    !         dgamma_a = 1.0_RP /{var x}
    Dc         = 0.146_RP - 5.964E-2_RP * log( Nc_def / 2000.0_RP )
    dgamma_a   = 0.1444_RP / Dc**2
    GAM_dgam   = MISC_gammafunc( dgamma_a )
    GAM_dgam23 = MISC_gammafunc( dgamma_a + 2.0_RP/3.0_RP )
    coef_dgam  = GAM_dgam23 / GAM_dgam * dgamma_a**(-2.0_RP/3.0_RP)
    coef_xf    = 3.0_RP / 4.0_RP / PI / rho_w
    OneThird   = 1.0_DP / 3.0_DP

    !$acc kernels async(1) pcopy(rctop,tctop,rctop_cld,tctop_cld)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,rctop,tctop,rctop_cld,tctop_cld,UNDEF)
    do ij = 1, ijdim
       rctop    (ij,1) = UNDEF
       tctop    (ij,1) = UNDEF
       rctop_cld(ij,1) = UNDEF
       tctop_cld(ij,1) = UNDEF
    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc kernels loop seq pcopy(rceff,rctop,tctop,rceff_cld,rctop_cld,tctop_cld) async(1) &
    !$acc pcopyin(rho,tem,q,Nc)
    !$omp parallel default(none), &
    !$omp private(ij,k,dens,temp,qc,qr,rhoqc,zerosw,RLMDr,sw,xf_qc,rf_qc,r2_qc,r2_qr,r3_qc,r3_qr) &
    !$omp shared(ijdim,kmin,kmax,rceff,rctop,tctop,rceff_cld,rctop_cld,tctop_cld, &
    !$omp        rho,tem,q,Nc,coef_dgam,coef_xf,Ar,OneThird,N0r,GAM_1br,GAM_3,UNDEF,I_QC,I_QR)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do k  = kmin, kmax
       !$acc loop independent
       !$omp do
       do ij = 1, ijdim
          dens  = rho(ij,k)
          temp  = tem(ij,k)
          qc    = max( q(ij,k,I_QC), 0.0_RP )
          qr    = max( q(ij,k,I_QR), 0.0_RP )
          rhoqc = dens * qc * 1000.0_RP ! [g/m3]

          ! slope parameter lambda (Rain)
          zerosw = 0.5_RP - sign(0.5_RP, qr - 1.E-12_RP )
          RLMDr  = sqrt(sqrt( dens * qr / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

          xf_qc = rhoqc / Nc(ij,k) * 1.E-9_RP                    ! mean mass   of qc [kg]
          rf_qc = ( coef_xf * xf_qc + 1.E-16 )**OneThird         ! mean radius of qc [m]
          r2_qc = coef_dgam * rf_qc**2 * ( Nc(ij,k) * 1.E+6_RP ) ! r^2 moment  of qc
          r2_qr = 0.25_RP * N0r * GAM_3 * RLMDr**3               ! r^2 moment  of qr
          r3_qc = coef_xf * dens * qc                            ! r^3 moment  of qc
          r3_qr = coef_xf * dens * qr                            ! r^3 moment  of qr

          zerosw = 0.5_RP - sign(0.5_RP, (r2_qc+r2_qr)-r2_min )
          rceff(ij,k) = ( r3_qc + r3_qr ) / ( r2_qc + r2_qr + zerosw ) * ( 1.0_RP - zerosw )

          sw = ( 0.5_RP + sign(0.5_RP, (qc+qr)-q_min ) ) ! if qc+qr > 0.01[g/kg], sw=1

          rctop(ij,1) = (        sw ) * rceff(ij,k) &
                      + ( 1.0_RP-sw ) * rctop(ij,1)
          tctop(ij,1) = (        sw ) * temp &
                      + ( 1.0_RP-sw ) * tctop(ij,1)

          zerosw = 0.5_RP - sign(0.5_RP, r2_qc-r2_min )
          rceff_cld(ij,k) = r3_qc / ( r2_qc + zerosw ) * ( 1.0_RP - zerosw )

          sw = ( 0.5_RP + sign(0.5_RP, qc-q_min ) ) ! if qc > 0.01[g/kg], sw=1

          rctop_cld(ij,1) = (        sw ) * rceff_cld(ij,k) &
                          + ( 1.0_RP-sw ) * rctop_cld(ij,1)
          tctop_cld(ij,1) = (        sw ) * temp &
                          + ( 1.0_RP-sw ) * tctop_cld(ij,1)
       enddo
       !$omp end do

    enddo
    !$omp end parallel
    !$acc end kernels

    !$acc kernels async(1) pcopy(rceff,rceff_cld,rwtop,rwtop_cld) pcopyin(tctop,rctop,tctop_cld,rctop_cld)
    !$omp parallel do default(none),private(ij,sw), &
    !$omp shared(ijdim,kmin,kmax,rceff,tctop,rctop,rwtop,rceff_cld,tctop_cld,rctop_cld,rwtop_cld,TEM00,UNDEF)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do ij = 1, ijdim
       rceff    (ij,kmin-1) = 0.0_RP
       rceff_cld(ij,kmin-1) = 0.0_RP
       rceff    (ij,kmax+1) = 0.0_RP
       rceff_cld(ij,kmax+1) = 0.0_RP

       sw = ( 0.5_RP + sign(0.5_RP, tctop(ij,1)-TEM00 ) ) ! if T > 0C, sw=1

       rwtop(ij,1) = (        sw ) * rctop(ij,1) &
                   + ( 1.0_RP-sw ) * UNDEF

       sw = ( 0.5_RP + sign(0.5_RP, tctop_cld(ij,1)-TEM00 ) ) ! if T > 0C, sw=1

       rwtop_cld(ij,1) = (        sw ) * rctop_cld(ij,1) &
                       + ( 1.0_RP-sw ) * UNDEF
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_06')
call PROF_rapstart('____MP_NSW6_07')
#endif

    !---< preciptation >---

!OCL XFILL
    !$acc kernels async(1) pcopy(Vt)
    !$omp parallel do default(none),private(ij,nq), &
    !$omp shared(nqmax,ijdim,kmin,kmax,Vt)
!OCL NOPREFETCH
    do nq = 1, nqmax
    do ij = 1, ijdim
       Vt(ij,kmin-1,nq) = 0.0_RP
       Vt(ij,kmax+1,nq) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !--- update mass concentration
!OCL XFILL
    !$acc kernels async(1) pcopy(q) pcopyin(rhog,rhogq)
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,q,rhogq,rhog)
!OCL NOPREFETCH
    do nq = NQW_STR, NQW_END
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20200905
!OCL XFILL
       !$acc kernels async(1) pcopy(q) pcopyin(rhog,rhogq)
       !$omp parallel do default(none),private(ij,k,nq), &
       !$omp shared(ijdim,kdim,ISO_STR,ISO_END,q,rhogq,rhog)
!OCL NOPREFETCH
       do nq = ISO_STR, ISO_END
       do k  = 1, kdim
       do ij = 1, ijdim
          q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
       enddo
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )
    call THRMDYN_cv( ijdim, kdim, qd(:,:), q(:,:,:), cva(:,:) )

!OCL XFILL
    !$acc kernels async(1) pcopy(tem,rgs,rgsh) pcopyin(rhog,rhoge,cva,gam2,gam2h,gsgam2,gsgam2h)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,tem,rhoge,rhog,cva,rgs,rgsh,gam2,gam2h,gsgam2,gsgam2h)
!OCL NOPREFETCH
    do k  = 1, kdim
    do ij = 1, ijdim
       tem (ij,k) = rhoge(ij,k) / ( rhog(ij,k) * cva(ij,k) )
       rgs (ij,k) = gam2 (ij,k) / gsgam2 (ij,k)
       rgsh(ij,k) = gam2h(ij,k) / gsgam2h(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_07')
call PROF_rapstart('____MP_NSW6_08z')
#endif

    preciptation_flag(:)    = .false.
    preciptation_flag(I_QR) = .true.
    preciptation_flag(I_QI) = .true.
    preciptation_flag(I_QS) = .true.
    preciptation_flag(I_QG) = .true.

    if ( precip_transport_type == '3WATER' ) then
       preciptation_flag(I_QI) = .false.
    endif

    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20200908
       preciptation_flag(ISO1_QR) = .true.
       preciptation_flag(ISO1_QI) = .true.
       preciptation_flag(ISO1_QS) = .true.
       preciptation_flag(ISO1_QG) = .true.
       preciptation_flag(ISO2_QR) = .true.
       preciptation_flag(ISO2_QI) = .true.
       preciptation_flag(ISO2_QS) = .true.
       preciptation_flag(ISO2_QG) = .true.

       if ( precip_transport_type == '3WATER' ) then
          preciptation_flag(ISO1_QI) = .false.
          preciptation_flag(ISO2_QI) = .false.
       endif
    endif

    if ( precip_scheme_type == 'Upwind-Euler' ) then

       call precip_transport_euler( ijdim,             &
                                    rhog,              &
                                    rhogvx,            &
                                    rhogvy,            &
                                    rhogvz,            &
                                    rhogw,             &
                                    rhoge,             &
                                    rhogq,             &
                                    rho,               &
                                    tem,               &
                                    pre,               &
                                    vx,                &
                                    vy,                &
                                    vz,                &
                                    w,                 &
                                    q,                 &
                                    qd,                &
                                    z,                 &
                                    zh,                &
                                    Vt,                &
                                    preciptation_flag, &
                                    precip,            &
                                    precip_rhoe,       &
                                    precip_lh_heat,    &
                                    precip_rhophi,     &
                                    precip_rhokin,     &
                                    precip_flux,       &
                                    gsgam2,            &
                                    gsgam2h,           &
                                    rgs,               &
                                    rgsh,              &
                                    ix,                &
                                    iy,                &
                                    iz,                &
                                    jx,                &
                                    jy,                &
                                    jz,                &
                                    dt                 )

    elseif( precip_scheme_type == 'Flux-Semilag_new' ) then

       call precip_transport_new( ijdim,             &
                                  rhog,              &
                                  rhogvx,            &
                                  rhogvy,            &
                                  rhogvz,            &
                                  rhogw,             &
                                  rhoge,             &
                                  rhogq,             &
                                  rho,               &
                                  tem,               &
                                  pre,               &
                                  vx,                &
                                  vy,                &
                                  vz,                &
                                  w,                 &
                                  q,                 &
                                  qd,                &
                                  z,                 &
                                  zh,                &
                                  Vt,                &
                                  preciptation_flag, &
                                  precip,            &
                                  precip_ISO1,       & ! [Add] M.Tanoue 20200908
                                  precip_ISO2,       & ! [Add] M.Tanoue 20200908
                                  precip_rhoe,       &
                                  precip_lh_heat,    &
                                  precip_rhophi,     &
                                  precip_rhokin,     &
                                  precip_flux,       &
                                  gsgam2,            &
                                  gsgam2h,           &
                                  rgs,               &
                                  rgsh,              &
                                  ix,                &
                                  iy,                &
                                  iz,                &
                                  jx,                &
                                  jy,                &
                                  jz,                &
                                  dt                 )

    else

       call precip_transport( ijdim,             &
                              rhog,              &
                              rhogvx,            &
                              rhogvy,            &
                              rhogvz,            &
                              rhogw,             &
                              rhoge,             &
                              rhogq,             &
                              rho,               &
                              tem,               &
                              pre,               &
                              vx,                &
                              vy,                &
                              vz,                &
                              w,                 &
                              q,                 &
                              qd,                &
                              z,                 &
                              zh,                &
                              Vt,                &
                              preciptation_flag, &
                              precip,            &
                              precip_rhoe,       &
                              precip_lh_heat,    &
                              precip_rhophi,     &
                              precip_rhokin,     &
                              precip_flux,       &
                              gsgam2,            &
                              gsgam2h,           &
                              rgs,               &
                              rgsh,              &
                              ix,                &
                              iy,                &
                              iz,                &
                              jx,                &
                              jy,                &
                              jz,                &
                              dt                 )

    endif
    ! for OpenACC
    thrmdyn_wait_on = .false.
    satadjustment_wait_on = .false.
#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_08z')
call PROF_rapstart('____MP_NSW6_09')
#endif

!OCL XFILL
    !$acc kernels async(1) pcopy(ml_Pconv,ml_Pconw,ml_Pconi) pcopyin(q)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,ml_Pconv,ml_Pconw,ml_Pconi,q,I_QV,I_QC,I_QI)
    do k  = 1, kdim
    do ij = 1, ijdim
       ml_Pconv(ij,k) = q(ij,k,I_QV)
       ml_Pconw(ij,k) = q(ij,k,I_QC)
       ml_Pconi(ij,k) = q(ij,k,I_QI)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_09')
call PROF_rapstart('____MP_NSW6_10z')
#endif

    if ( OPT_EXPLICIT_ICEGEN ) then

       if ( .NOT. OPT_SCLW ) then ! [Add] 2018/02/08 T.Seiki add option not to use saturation adjustment
          call SATURATION_setrange( 100.0_RP, 90.0_RP )

          call SATURATION_adjustment( ijdim,             & ! [IN]
                                      kdim,              & ! [IN]
                                      rhog  (:,:),       & ! [IN]
                                      rhoge (:,:),       & ! [INOUT]
                                      rhogq (:,:,:),     & ! [INOUT]
                                      tem   (:,:),       & ! [INOUT]
                                      q     (:,:,:),     & ! [INOUT]
                                      qd    (:,:),       & ! [IN]
                                      gsgam2(:,:),       & ! [IN]
                                      ice_adjust=.false. ) ! [IN]

          if ( opt_water_isotope ) then ! [Add] M.Tanoue 20201210
             call isotope_condensation_liq( ijdim,            &  !--[IN]
                                            kdim,             &  !--[IN]
                                            nqmax,            &  !--[IN]
                                            rho      (:,:),   &  !--[IN]
                                            rhog     (:,:),   &  !--[IN]
                                            rhogq    (:,:,:), &  !--[INOUT]
                                            tem      (:,:),   &  !--[IN]
                                            pre      (:,:),   &  !--[IN]
                                            q        (:,:,:)  )  !--[INOUT]
          endif

       else
          ! condensation(evaporation) of qc and qi is sponstaneously solved together.
          call explicit_condensation( ijdim,        & ! [IN]
                                      kmin,         & ! [IN]
                                      kmax,         & ! [IN]
                                      kdim,         & ! [IN]
                                      nqmax,        & ! [IN]
                                      rhog (:,:),   & ! [IN]
                                      rhoge(:,:),   & ! [INOUT]
                                      rhogq(:,:,:), & ! [INOUT] rhogqv,qc,qi are updated
                                      rho  (:,:),   & ! [IN]
                                      tem  (:,:),   & ! [INOUT]
                                      pre  (:,:),   & ! [IN]
                                      q    (:,:,:), & ! [INOUT]
                                      Nc   (:,:),   & ! [IN]
                                      w    (:,:),   & ! [IN]
                                      dt            ) ! [IN]
       endif

    else

       call SATURATION_setrange( 273.16_RP, 233.16_RP )

       call SATURATION_adjustment( ijdim,             & ! [IN]
                                   kdim,              & ! [IN]
                                   rhog  (:,:),       & ! [IN]
                                   rhoge (:,:),       & ! [INOUT]
                                   rhogq (:,:,:),     & ! [INOUT]
                                   tem   (:,:),       & ! [INOUT]
                                   q     (:,:,:),     & ! [INOUT]
                                   qd    (:,:),       & ! [IN]
                                   gsgam2(:,:),       & ! [IN]
                                   ice_adjust=.true.  ) ! [IN]

    endif

    ! for OpenACC
    thrmdyn_wait_on = .false.
    satadjustment_wait_on = .false.
#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_10z')
call PROF_rapstart('____MP_NSW6_11')
#endif

    !$acc kernels async(1) pcopy(ml_Pconv,ml_Pconw,ml_Pconi) pcopyin(q)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,ml_Pconv,ml_Pconw,ml_Pconi,q,I_QV,I_QC,I_QI)
    do k  = 1, kdim
    do ij = 1, ijdim
       ml_Pconv(ij,k) = q(ij,k,I_QV) - ml_Pconv(ij,k)
       ml_Pconw(ij,k) = q(ij,k,I_QC) - ml_Pconw(ij,k)
       ml_Pconi(ij,k) = q(ij,k,I_QI) - ml_Pconi(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    ! Individual tendency output
!     do ip = 1, wk_nmax
!        call history_in( 'ml_'//trim(w_name(ip)), ml_wk(:,:,ip) )
!     enddo
!--- MINIAPP: history_in removed together with mod_history. ---
#ifndef MINIAPP
    call history_in( 'ml_Pconv', ml_Pconv(:,:) )
    call history_in( 'ml_Pconw', ml_Pconw(:,:) )
    call history_in( 'ml_Pconi', ml_Pconi(:,:) )
#endif

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_11')
call PROF_rapstart('____MP_NSW6_12z')
#endif

    call negative_filter( ijdim,         & ! [IN]
                          rhog  (:,:),   & ! [INOUT]
                          rhoge (:,:),   & ! [INOUT]
                          rhogq (:,:,:), & ! [INOUT]
                          rho   (:,:),   & ! [INOUT]
                          tem   (:,:),   & ! [IN]
                          pre   (:,:),   & ! [INOUT]
                          q     (:,:,:), & ! [INOUT]
                          gsgam2(:,:)    ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('____MP_NSW6_12z')
#endif
    thrmdyn_wait_on = .true.
    satadjustment_wait_on = .true.
    !$acc wait(1)
    call PROF_rapend  ('____MP_NSW6',3)

    !$acc end data
    return
  end subroutine mp_nsw6

  !-----------------------------------------------------------------------------
  ![Add] 2018/02/08 T.Seiki
  subroutine explicit_condensation( &
       ijdim, &
       kmin,  &
       kmax,  &
       kdim,  &
       nqmax, &
       rhog,  &
       rhoge, &
       rhogq, &
       rho,   &
       tem,   &
       pre,   &
       q,     &
       Nc,    &
       w,     &
       dt     )
    use mod_cnst, only: &
       CNST_UNDEF, &
       CNST_EPS,   &
       CNST_PI,    &
       CNST_GRAV,  &
       CNST_Rvap,  &
       CNST_LHV0,  &
       CNST_LHS0,  &
       CNST_LHF0,  &
       CNST_TEM00, &
       CNST_Pstd
    use mod_runconf, only: &
       I_QV, I_QC, I_QI,          &
       opt_water_isotope,         &
       ISO1_QV, ISO1_QC, ISO1_QI, &
       ISO2_QV, ISO2_QC, ISO2_QI, &
       LHV, LHF, LHS
    use mod_thrmdyn, only: &
       thrmdyn_qd, &
       thrmdyn_cv, &
       thrmdyn_cp
    use mod_satadjust, only: &
       SATURATION_psat_liq,  &
       SATURATION_psat_ice,  &
       moist_dqsw_dtem_rho,  &
       moist_dqsw_dtem_dpre, &
       moist_dqsi_dtem_rho,  &
       moist_dqsi_dtem_dpre
    use mod_isotope, only: &
       NISO,         &
       DIFF_H18O,    &
       DIFF_HDO,     &
       WISO_eqfract, &
       ALP0
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kmin
    integer,  intent(in)    :: kmax
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: nqmax
    real(RP), intent(in)    :: rhog (ijdim,kdim)
    real(RP), intent(inout) :: rhoge(ijdim,kdim)
    real(RP), intent(inout) :: rhogq(ijdim,kdim,nqmax)
    real(RP), intent(in)    :: rho  (ijdim,kdim)
    real(RP), intent(inout) :: tem  (ijdim,kdim)
    real(RP), intent(in)    :: pre  (ijdim,kdim)
    real(RP), intent(inout) :: q    (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: Nc   (ijdim,kdim)
    real(RP), intent(in)    :: w    (ijdim,kdim)
    real(DP), intent(in)    :: dt

    real(RP) :: dqswdtem_rho(ijdim,kdim)         !< (dqsw/dtem)_rho
    real(RP) :: dqswdtem_pre(ijdim,kdim)         !< (dqsw/dtem)_pre
    real(RP) :: dqswdpre_tem(ijdim,kdim)         !< (dqsw/dpre)_tem
    real(RP) :: dqsidtem_rho(ijdim,kdim)         !< (dqsi/dtem)_rho
    real(RP) :: dqsidtem_pre(ijdim,kdim)         !< (dqsi/dtem)_pre
    real(RP) :: dqsidpre_tem(ijdim,kdim)         !< (dqsi/dpre)_tem
    real(RP) :: qd          (ijdim,kdim)
    real(RP) :: cpa         (ijdim,kdim)         !< specific heat of atmosphere at constant pressure
    real(RP) :: cva         (ijdim,kdim)         !< specific heat of atmosphere at constant volume
    real(RP) :: psatl       (ijdim,kdim)
    real(RP) :: psati       (ijdim,kdim)

    real(RP) :: drhogqc     (ijdim,kdim)
    real(RP) :: drhogqi     (ijdim,kdim)

    real(RP) :: dens                             !< density
    real(RP) :: temp                             !< T [K]
    real(RP) :: qv                               !< mixing ratio of water vapor  [kg/kg]
    real(RP) :: qc                               !< mixing ratio of liquid water [kg/kg]
    real(RP) :: qi                               !< mixing ratio of ice water    [kg/kg]
    real(RP) :: qv_t                             !< tendency     of water vapor  [kg/kg/s]
    real(RP) :: qc_t                             !< tendency     of liquid water [kg/kg/s]
    real(RP) :: qi_t                             !< tendency     of ice water    [kg/kg/s]
    real(RP) :: qsatl                            !< saturated water vapor for liquid water [kg/kg]
    real(RP) :: qsati                            !< saturated water vapor for ice water    [kg/kg]
    real(RP) :: Sliq                             !< saturated ratio S for liquid water [0-1]
    real(RP) :: Sice                             !< saturated ratio S for ice water    [0-1]
    real(RP) :: Rdens                            !< 1 / density
    real(RP) :: temc                             !< T - T0 [K]
    real(RP) :: rhoqc                            !< rho * qc
    real(RP) :: rhoqi                            !< rho * qi
    real(RP) :: dc_ave                           ! mean volume diameter of cloud water [m]
    real(RP) :: Ka                               !< thermal diffusion coefficient of air
    real(RP) :: Kd                               !< diffusion coefficient of water vapor in air
    real(RP) :: Glv, Giv                         !< thermodynamic function

    !---< Seiki and Roh (2020) >---
    real(RP) :: Sice_eq                          !< equivalent saturated ratio S for ice water[0-1]
    real(RP) :: XNi_def, XNi_p07
    real(RP) :: XNi, XMi, Di
    real(RP) :: sw1, sw2

    real(RP) :: r_tau_cnd_c, r_tau_dep_i
    real(RP) :: aliqliq, asolliq
    real(RP) :: aliqsol, asolsol
    real(RP) :: Pdynliq, Pradliq
    real(RP) :: Pdynsol, Pradsol
    real(RP) :: Acnd, Adep
    real(RP) :: r_tau_cnd, r_tau_dep
    real(RP) :: sw_cnd, sw_dep

    ! cloud water condensation & cloud ice deposition
    real(RP) :: plccnd, plccnd1, plccnd2         !< condensation rate of cloud water [kg/m3/s]
    real(RP) :: plidep, plidep1, plidep2         !< condensation rate of cloud water [kg/m3/s]
    real(RP) :: tau_cnd, tau_dep
    real(RP) :: coef_a_cnd, coef_b_cnd
    real(RP) :: coef_a_dep, coef_b_dep

    ! limiter
    real(RP) :: qci_t
    real(RP) :: qv_t_min, qv_t_max
    real(RP) :: fac

    ! WISO
    real(RP) :: RDD  (NISO)
    real(RP) :: ALP  (NISO)
    real(RP) :: ALPI (NISO)
    real(RP) :: ALPKL(NISO)
    real(RP) :: ALPKS(NISO)
    real(RP) :: Rqv  (NISO)
    real(RP) :: Rqc  (NISO)
    real(RP) :: Rqi  (NISO)
    real(RP) :: drhogqc1, drhogqc2
    real(RP) :: drhogqi1, drhogqi2
    real(RP) :: rhogqv, rhogqv1, rhogqv2
    real(RP) :: rhogqc, rhogqc1, rhogqc2
    real(RP) :: rhogqi, rhogqi1, rhogqi2
    real(RP) :: Sleff, Sieff
    real(RP) :: fidep, ff
    integer  :: iso

    real(RP) :: UNDEF, EPS, PI, GRAV, Rvap, LHV0, LHS0, LHF0, TEM00, PRE00
    real(RP) :: dt_RP
    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(dqswdtem_rho,dqswdtem_pre,dqswdpre_tem,dqsidtem_rho,dqsidtem_pre,dqsidpre_tem) &
    !$acc pcreate(qd,cpa,cva,psatl,psati,drhogqc,drhogqi) &
    !$acc pcreate(RDD,ALP,ALPI,ALPKL,ALPKS,Rqv,Rqc,Rqi) &
    !$acc pcopyin(rhog,rho,pre,Nc,w) &
    !$acc pcopy(rhoge,rhogq,tem,q)

    call PROF_rapstart('_____MP_NSW6_EXCOND',3)

    UNDEF = CNST_UNDEF
    EPS   = CNST_EPS
    PI    = CNST_PI
    GRAV  = CNST_GRAV
    Rvap  = CNST_Rvap
    LHV0  = CNST_LHV0
    LHS0  = CNST_LHS0
    LHF0  = CNST_LHF0
    TEM00 = CNST_TEM00
    PRE00 = CNST_Pstd

    dt_RP = real(dt,kind=RP)

    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20220318
       RDD(1) = DIFF_H18O
       RDD(2) = DIFF_HDO
    endif

    call SATURATION_psat_liq ( ijdim, kdim, tem(:,:), psatl(:,:) )
    call moist_dqsw_dtem_rho ( ijdim, kdim, tem(:,:), rho(:,:), dqswdtem_rho(:,:) )
    call moist_dqsw_dtem_dpre( ijdim, kdim, tem(:,:), pre(:,:), dqswdtem_pre(:,:), dqswdpre_tem(:,:) )

    call SATURATION_psat_ice ( ijdim, kdim, tem(:,:), psati(:,:) )
    call moist_dqsi_dtem_rho ( ijdim, kdim, tem(:,:), rho(:,:), dqsidtem_rho(:,:) )
    call moist_dqsi_dtem_dpre( ijdim, kdim, tem(:,:), pre(:,:), dqsidtem_pre(:,:), dqsidpre_tem(:,:) )

    call thrmdyn_qd( ijdim, kdim, q(:,:,:), qd(:,:) )
    call thrmdyn_cv( ijdim, kdim, qd(:,:), q(:,:,:), cva(:,:) )
    call thrmdyn_cp( ijdim, kdim, qd(:,:), q(:,:,:), cpa(:,:) )

    !$acc kernels pcopy(drhogqc,drhogqi) pcopyin(rhog,rho,tem,pre,q,w,Nc,qd,cpa,cva,psatl,psati) &
    !$acc         pcopyin(dqswdtem_rho,dqswdtem_pre,dqswdpre_tem,dqsidtem_rho,dqsidtem_pre,dqsidpre_tem)
    !$acc loop independent
    !$omp parallel do default(none), &
    !$omp private(ij,k,dens,temp,qv,qc,qi,qc_t,qi_t,qsatl,qsati,Sliq,Sice,Rdens,temc,    &
    !$omp         rhoqc,rhoqi,dc_ave,Ka,Kd,Glv,Giv,Sice_eq,XNi_def,XNi_p07,XNi,XMi,Di,sw1,sw2, &
    !$omp         r_tau_cnd_c,r_tau_dep_i,aliqliq,asolliq,aliqsol,asolsol,Pdynliq,       &
    !$omp         Pradliq,Pdynsol,Pradsol,Acnd,Adep,r_tau_cnd,r_tau_dep,sw_cnd,sw_dep,   &
    !$omp         plccnd,plccnd1,plccnd2,plidep,plidep1,plidep2,tau_cnd,tau_dep,         &
    !$omp         coef_a_cnd,coef_b_cnd,coef_a_dep,coef_b_dep,qci_t,qv_t_min,qv_t_max,fac) &
    !$omp shared(ijdim,kmin,kmax,drhogqc,drhogqi,rhog,rho,tem,pre,q,w,Nc,                &
    !$omp        psatl,psati,dqswdtem_rho,dqswdtem_pre,dqswdpre_tem,dqsidtem_rho,        &
    !$omp        dqsidtem_pre,dqsidpre_tem,qd,cva,cpa,                                   &
    !$omp        I_QV,I_QC,I_QI,EPS,PI,GRAV,Rvap,LHV0,LHS0,LHF0,TEM00,PRE00,LHV,LHF,LHS, &
    !$omp        sw_phillips2007,dt_RP)
    do k  = kmin, kmax
    !ocl loop_fission_target(ls)
    do ij = 1, ijdim

       dens     = rho(ij,k)
       temp     = tem(ij,k)
       qv       = max( q(ij,k,I_QV), 0.0_RP )
       qc       = max( q(ij,k,I_QC), 0.0_RP )
       qi       = max( q(ij,k,I_QI), 0.0_RP )

       ! saturation ratio S
       qsatl    = psatl(ij,k) / ( dens * Rvap * temp )
       qsati    = psati(ij,k) / ( dens * Rvap * temp )
       Sliq     = qv / max( qsatl, EPS )
       Sice     = qv / max( qsati, EPS )

       Rdens    = 1.0_RP / dens
       temc     = temp - TEM00
       rhoqc    = max(dens*qc,EPS)
       rhoqi    = max(dens*qi,EPS)
       dc_ave   = 0.124_RP * ( rhoqc / (Nc(ij,k)*1.E+6_RP) )**0.33333_RP

       Ka  = ( Ka0 + dKa_dT * temc )
       Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
       Glv = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsatl) )
       Giv = 1.0_RP / ( LHS0/(Ka*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsati) )

       XNi_def = min( max( 5.38E+7_RP * rhoqi**0.75_RP, 1.E+3_RP ), 1.E+6_RP )

       ! [Add] 2018/02/08 T.Seiki, Phillips et al. (2007), JAS
       Sice_eq = min( 3.0_RP, (qv+qi) / max(qsati,EPS) )

       sw1 = ( 0.5_RP - sign(0.5_RP, temc + 30.0_RP ) ) ! if T < -30C, sw1=1
       sw2 = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
           * ( 0.5_RP - sign(0.5_RP, temc           ) ) ! if -30C < T < 0C, sw2=1

       XNi_p07 = (        sw1 ) * 1.E+3_RP * exp( 12.96_RP*(Sice_eq-1.1_RP) * 0.3_RP   )           & ! if        temc < -30C
               + (        sw2 ) * 1.E+3_RP * exp( 12.96_RP*(Sice_eq-1.0_RP) - 0.639_RP ) * 0.06_RP & ! if -30C < temc <   0C
               + ( 1.0_RP-sw1 ) * ( 1.0_RP-sw2 ) * 1.E+3_RP                                          ! if   0C < temc

       XNi = ( 1.0_RP-sw_phillips2007 ) * XNi_def &
           + (        sw_phillips2007 ) * XNi_p07
       XMi = rhoqi / XNi
       Di  = min( Di_a * sqrt(XMi), Di_max )

       r_tau_cnd_c = ( 2.0_RP * PI * dc_ave * Nc(ij,k)*1.E+6_RP * Glv ) / ( dens * qsatl )
       r_tau_dep_i = ( 4.0_RP      * Di     * XNi               * Giv ) / ( dens * qsati )

       aliqliq     = 1.0_RP + LHV / cva(ij,k) * dqswdtem_rho(ij,k)
       asolliq     = 1.0_RP + LHS / cva(ij,k) * dqswdtem_rho(ij,k)
       aliqsol     = 1.0_RP + LHV / cva(ij,k) * dqsidtem_rho(ij,k)
       asolsol     = 1.0_RP + LHS / cva(ij,k) * dqsidtem_rho(ij,k)

       Pdynliq     = w(ij,k) * GRAV * ( dqswdtem_pre(ij,k) / cpa(ij,k) + dens * dqswdpre_tem(ij,k) )
       Pdynsol     = w(ij,k) * GRAV * ( dqsidtem_pre(ij,k) / cpa(ij,k) + dens * dqsidpre_tem(ij,k) )
       Pradliq     = 0.0_RP
       Pradsol     = 0.0_RP

       Acnd        = Pdynliq + Pradliq - r_tau_dep_i * max( qsatl-qsati, 0.0_RP )
       Adep        = Pdynsol + Pradsol + r_tau_cnd_c * max( qsatl-qsati, 0.0_RP )

       r_tau_cnd   = aliqliq * r_tau_cnd_c + asolliq * r_tau_dep_i
       r_tau_dep   = aliqsol * r_tau_cnd_c + asolsol * r_tau_dep_i

       sw_cnd = 0.5_RP + sign(0.5_RP, r_tau_cnd - EPS ) ! if r_tau_cnd > EPS, sw=1
       sw_dep = 0.5_RP + sign(0.5_RP, r_tau_dep - EPS ) ! if r_tau_dep > EPS, sw=1

       !--- cloud water condensation

       ! (case1) condensation/evaporation is negligible
       plccnd1    = 2.0_RP * PI * Rdens * dc_ave * Nc(ij,k) * 1.E+6_RP * (Sliq-1.0_RP) * Glv
       ! (case2) condensation is calculated assuming time variability of Sliq in a time-step using a relaxation equation.
       tau_cnd    = 1.0_RP / max(r_tau_cnd,EPS)
       coef_a_cnd = Acnd * tau_cnd
       coef_b_cnd = tau_cnd/dt_RP * ( qv - qsatl - Acnd*tau_cnd ) * ( exp(-dt_RP*r_tau_cnd)-1.0_RP )
       plccnd2    = coef_a_cnd * r_tau_cnd_c - coef_b_cnd * r_tau_cnd_c

       plccnd = ( 1.0_RP-sw_cnd ) * plccnd1 &
              + (        sw_cnd ) * plccnd2

       ! filter for condensation/evaporation of qc
       if    ( plccnd >= 0.0_RP .AND. qv >= qsatl ) then ! always saturated
          qc_t = min( plccnd, (qv-qsatl)/dt_RP )
       elseif( plccnd <= 0.0_RP .AND. qv <= qsatl ) then ! always evaporation
          qc_t = max( plccnd, (qv-qsatl)/dt_RP, -qc/dt_RP )
       elseif( plccnd <= 0.0_RP .AND. qv >= qsatl ) then ! partially unsaturated during time-step
          qc_t = max( plccnd, -qc/dt_RP )
       else ! partial condensation is not accepted here
          qc_t = 0.0_RP
       endif

       !--- cloud ice deposition

       ! (case1) condensation/evaporation is negligible
       plidep1    = 4.0_RP * Di * XNi * Rdens * (Sice-1.0_RP) * Giv
       ! (case2) condensation is calculated assuming time variability of Sliq in a time-step using a relaxation equation.
       tau_dep    = 1.0_RP / max(r_tau_dep,EPS)
       coef_a_dep = Adep * tau_dep
       coef_b_dep = tau_dep/dt_RP * ( qv - qsati - Adep*tau_dep ) * ( exp(-dt_RP*r_tau_dep)-1.0_RP )
       plidep2    = coef_a_dep * r_tau_dep_i - coef_b_dep * r_tau_dep_i

       plidep = ( 1.0_RP-sw_dep ) * plidep1 &
              + (        sw_dep ) * plidep2

       ! filter for deposition/sublimation of qi
       if    ( plidep >= 0.0_RP .AND. qv >= qsati ) then ! always saturated
          qi_t = min( plidep, (qv-qsati)/dt_RP )
       elseif( plidep <= 0.0_RP .AND. qv <= qsati ) then ! always evaporation
          qi_t = max( plidep, (qv-qsati)/dt_RP, -qi/dt_RP )
       elseif( plidep <= 0.0_RP .AND. qv >= qsati ) then ! partially unsaturated during time-step
          qi_t = max( plidep, -qi/dt_RP )
       else ! partial condensation is not accepted here
          qi_t = 0.0_RP
       endif

       ! filter not to overestimate "cnd+dep" or "evp+sub"
       qci_t = qc_t + qi_t

       qv_t_min = max( qv - min(qsatl,qsati), 0.0_RP ) / dt_RP
       qv_t_max = min( qv - max(qsatl,qsati), 0.0_RP ) / dt_RP

       if    ( qci_t >  EPS ) then
          fac = min( qci_t, qv_t_min ) / qci_t
       elseif( qci_t < -EPS ) then
          fac = max( qci_t, qv_t_max ) / qci_t
       else
          fac = 1.0_RP
       endif

       drhogqc(ij,k) = rhog(ij,k) * qc_t * dt_RP * fac
       drhogqi(ij,k) = rhog(ij,k) * qi_t * dt_RP * fac
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20220318
!       !$acc kernels pcopy(rhogq,q) pcopyin(rhog,rho,tem,pre,psatl,psati) &
!       !$acc         pcopyin(RDD,ALP,ALPI,drhogqc,drhogqi)
       !$omp parallel do default(none), &
       !$omp private(ij,k,dens,temp,qv,rhogqv,rhogqc,rhogqi,rhogqv1,rhogqv2,rhogqc1,rhogqc2, &
       !$omp         rhogqi1,rhogqi2,qsatl,qsati,Sliq,Sice,Ka,Kd,Glv,Giv,Sleff,Sieff,        &
       !$omp         ALPKL,ALPKS,Rqv,Rqc,Rqi,fidep,ff,drhogqc1,drhogqc2,drhogqi1,drhogqi2    ) &
       !$omp shared(ijdim,kmin,kmax,rhogq,q,rhog,rho,tem,pre,psatl,psati,temc,drhogqc,drhogqi, &
       !$omp        I_QV,I_QC,I_QI,ISO1_QV,ISO1_QC,ISO1_QI,ISO2_QV,ISO2_QC,ISO2_QI,            &
       !$omp        RDD,ALP,ALPI,EPS,Rvap,LHV0,LHS0,TEM00,PRE00,ALP0                           )
       do k  = kmin, kmax
       do ij = 1, ijdim

          dens     = rho(ij,k)
          temp     = tem(ij,k)
          qv       = max( q    (ij,k,I_QV)   , 0.0_RP )
          rhogqv   = max( rhogq(ij,k,I_QV)   , 0.0_RP )
          rhogqc   = max( rhogq(ij,k,I_QC)   , 0.0_RP )
          rhogqi   = max( rhogq(ij,k,I_QI)   , 0.0_RP )
          rhogqv1  =      rhogq(ij,k,ISO1_QV)
          rhogqv2  =      rhogq(ij,k,ISO2_QV)
          rhogqc1  =      rhogq(ij,k,ISO1_QC)
          rhogqc2  =      rhogq(ij,k,ISO2_QC)
          rhogqi1  =      rhogq(ij,k,ISO1_QI)
          rhogqi2  =      rhogq(ij,k,ISO2_QI)

          ! saturation ratio S
          qsatl    = psatl(ij,k) / ( dens * Rvap * temp )
          qsati    = psati(ij,k) / ( dens * Rvap * temp )
          Sliq     = qv / max( qsatl, EPS )
          Sice     = qv / max( qsati, EPS )

          temc     = temp - TEM00

          Ka       = ( Ka0 + dKa_dT * temc )
          Kd       = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
          Glv      = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsatl) )
          Giv      = 1.0_RP / ( LHS0/(Ka*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*qsati) )

          Sleff    = 1.0_RP / ( 1.0_RP - Glv * ( 1.0_RP - 1.0_RP / Sliq ) )
          Sieff    = 1.0_RP / ( 1.0_RP - Giv * ( 1.0_RP - 1.0_RP / Sice ) )

          call WISO_eqfract( temp, real(ALP0,kind=RP_G), ALP, ALPI )

          do iso = 1, NISO
             ALPKL(iso) = ALP (iso) / ( 1.0_RP + ( ALP (iso) * RDD(iso) - 1.0_RP ) * ( 1.0_RP - 1.0_RP / Sleff ) )
             ALPKS(iso) = ALPI(iso) / ( 1.0_RP + ( ALPI(iso) * RDD(iso) - 1.0_RP ) * ( 1.0_RP - 1.0_RP / Sieff ) )
          enddo

          if ( rhogqv >  0.0_RP ) then
             Rqv(1) =          rhogqv1       / rhogqv
             Rqv(2) =          rhogqv2       / rhogqv
             fidep  =          drhogqi(ij,k) / rhogqv
             ff     = 1.0_RP - drhogqc(ij,k) / rhogqv
          else
             Rqv(1) = 0.0_RP
             Rqv(2) = 0.0_RP
             fidep  = 0.0_RP
             ff     = 0.0_RP
          endif

          if ( rhogqc > 0.0_RP ) then
             Rqc(1) = rhogqc1 / rhogqc
             Rqc(2) = rhogqc2 / rhogqc
          else
             Rqc(1) = 0.0_RP
             Rqc(2) = 0.0_RP
          endif

          if ( rhogqi > 0.0_RP ) then
             Rqi(1) = rhogqi1 / rhogqi
             Rqi(2) = rhogqi2 / rhogqi
          else
             Rqi(1) = 0.0_RP
             Rqi(2) = 0.0_RP
          endif

          ! Cloud water evaporation/condensation
          ! *** Under developing (2022/04/28) ***
          if ( drhogqc(ij,k) >= 0.0_RP ) then
             ! Condensation (Saturation adjustment)
             ! *** Currently, SATADJ produces a large d-excess values of precipitation over cold regions. ***
             !drhogqci(1) = ( rhogq(ij,k,ISO1_QV)+rhogq(ij,k,ISO1_QC) ) * ALP(1)*(rhogqc+drhogqc(ij,k)) &
             !            / ( rhogqv-drhogqc(ij,k) + ALP(1)*(rhogqc+drhogqc(ij,k)) ) - rhogq(ij,k,ISO1_QC)
             !drhogqci(2) = ( rhogq(ij,k,ISO2_QV)+rhogq(ij,k,ISO2_QC) ) * ALP(2)*(rhogqc+drhogqc(ij,k)) &
             !            / ( rhogqv-drhogqc(ij,k) + ALP(2)*(rhogqc+drhogqc(ij,k)) ) - rhogq(ij,k,ISO2_QC)

             ! Condensation (Closed system)
             drhogqc1 = rhogqv * ( 1.0_DP-ff ) * ( ff * ( 1.0_DP - 1.0_DP/ALPKL(1) ) + 1.0_DP/ALPKL(1) )
             drhogqc2 = rhogqv * ( 1.0_DP-ff ) * ( ff * ( 1.0_DP - 1.0_DP/ALPKL(1) ) + 1.0_DP/ALPKL(1) )
          else
             ! Evaporation (Saturation adjustment)
             ! *** Currently, SATADJ produces a large d-excess values of precipitation over cold regions. ***
             !drhogqci(1) = ( rhogq(ij,k,ISO1_QV)+rhogq(ij,k,ISO1_QC) ) * ALP(1)*(rhogqc+drhogqc(ij,k)) &
             !            / ( rhogqv-drhogqc(ij,k) + ALP(1)*(rhogqc+drhogqc(ij,k)) ) - rhogq(ij,k,ISO1_QC)
             !drhogqci(2) = ( rhogq(ij,k,ISO2_QV)+rhogq(ij,k,ISO2_QC) ) * ALP(2)*(rhogqc+drhogqc(ij,k)) &
             !            / ( rhogqv-drhogqc(ij,k) + ALP(2)*(rhogqc+drhogqc(ij,k)) ) - rhogq(ij,k,ISO2_QC)

             ! Evaporation (Closed system)
             if ( rhogqc + drhogqc(ij,k) > EPS ) then  !! avoid numerically instability
                drhogqc1 = rhogqv * ( 1.0_DP-ff ) * ( ff * ( 1.0_DP - 1.0_DP/ALPKL(1) ) + 1.0_DP/ALPKL(1) )
                drhogqc2 = rhogqv * ( 1.0_DP-ff ) * ( ff * ( 1.0_DP - 1.0_DP/ALPKL(1) ) + 1.0_DP/ALPKL(1) )
             else
                ! Occur numerically instability under rhogqc < EPS (about -30 order).
                !drhogqci(1) = drhogqc(ij,k)*Rqv(1)*ALPKL(1)
                !drhogqci(2) = drhogqc(ij,k)*Rqv(2)*ALPKL(2)
                ! (Temporally) Saturation adjustment
                drhogqc1 = ( rhogqv1 + rhogqc1 ) * ALP(1) * ( rhogqc + drhogqc(ij,k) ) &
                         / ( rhogqv - drhogqc(ij,k) + ALP(1) * ( rhogqc + drhogqc(ij,k) ) ) - rhogqc1
                drhogqc2 = ( rhogqv2 + rhogqc2 ) * ALP(2) * ( rhogqc + drhogqc(ij,k) ) &
                         / ( rhogqv - drhogqc(ij,k) + ALP(2) * ( rhogqc + drhogqc(ij,k) ) ) - rhogqc2
             endif
          endif

          ! Cloud ice deposition/sublimation
          if ( drhogqi(ij,k) >= 0.0_RP ) then
             ! Deposition (with isotope fractionation)
             drhogqi1 = rhogqv1 * ( 1.0_RP - ( 1.0_RP - fidep )**ALPKS(1) )
             drhogqi2 = rhogqv2 * ( 1.0_RP - ( 1.0_RP - fidep )**ALPKS(2) )
          else
             ! Sublimation (without isotope fractionation)
             drhogqi1 = drhogqi(ij,k) * Rqi(1)
             drhogqi2 = drhogqi(ij,k) * Rqi(2)
          endif

          rhogq(ij,k,ISO1_QV) = max( rhogqv1 - drhogqc1 -drhogqi1, 0.0_RP )
          rhogq(ij,k,ISO1_QC) = max( rhogqc1 + drhogqc1          , 0.0_RP )
          rhogq(ij,k,ISO1_QI) = max( rhogqi1 + drhogqi1          , 0.0_RP )
          q    (ij,k,ISO1_QV) = rhogq(ij,k,ISO1_QV) / rhog(ij,k)
          q    (ij,k,ISO1_QC) = rhogq(ij,k,ISO1_QC) / rhog(ij,k)
          q    (ij,k,ISO1_QI) = rhogq(ij,k,ISO1_QI) / rhog(ij,k)

          rhogq(ij,k,ISO2_QV) = max( rhogqv2 - drhogqc2 -drhogqi2, 0.0_RP )
          rhogq(ij,k,ISO2_QC) = max( rhogqc2 + drhogqc2          , 0.0_RP )
          rhogq(ij,k,ISO2_QI) = max( rhogqi2 + drhogqi2          , 0.0_RP )
          q    (ij,k,ISO2_QV) = rhogq(ij,k,ISO2_QV) / rhog(ij,k)
          q    (ij,k,ISO2_QC) = rhogq(ij,k,ISO2_QC) / rhog(ij,k)
          q    (ij,k,ISO2_QI) = rhogq(ij,k,ISO2_QI) / rhog(ij,k)
       enddo
       enddo
       !$omp end parallel do
!       !$acc end kernels
    endif

    !$acc kernels pcopy(rhogq,rhoge,q,tem) pcopyin(rhog,drhogqc,drhogqi,cva)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhogq,rhoge,q,tem,rhog,drhogqc,drhogqi,cva,I_QV,I_QC,I_QI,LHV,LHF)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) - drhogqc(ij,k) - drhogqi(ij,k)
       rhogq(ij,k,I_QC) = rhogq(ij,k,I_QC) + drhogqc(ij,k)
       rhogq(ij,k,I_QI) = rhogq(ij,k,I_QI) + drhogqi(ij,k)
       rhoge(ij,k)      = rhoge(ij,k)      + drhogqc(ij,k) * LHV &
                                           + drhogqi(ij,k) * LHF

       q    (ij,k,I_QV) = rhogq(ij,k,I_QV) / rhog(ij,k)
       q    (ij,k,I_QC) = rhogq(ij,k,I_QC) / rhog(ij,k)
       q    (ij,k,I_QI) = rhogq(ij,k,I_QI) / rhog(ij,k)
       tem  (ij,k)      = rhoge(ij,k)      / rhog(ij,k) / cva(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    call PROF_rapend  ('_____MP_NSW6_EXCOND',3)

    !$acc end data
    return
  end subroutine explicit_condensation

  !-----------------------------------------------------------------------------
!  subroutine isotope_condensation_pre( &
!     ijdim,  &  !--[IN]
!     kdim,   &  !--[IN]
!     nqmax,  &  !--[IN]
!     rho,    &  !--[IN]
!     tem,    &  !--[IN]
!     pre,    &  !--[IN]
!     qsatl,  &  !--[IN]
!     q,      &  !--[IN]
!     qe_old  )  !--[OUT]
!    use mod_adm, &
!       only: ADM_LOG_FID
!    use mod_runconf, only: &
!         I_QV,    &
!         I_QC,    &
!         I_QI,    &
!         ISO1_QV, &
!         ISO2_QV, &
!         ISO1_QC, &
!         ISO2_QC, &
!         ISO1_QI, &
!         ISO2_QI
!    use mod_cnst, only: &
!         CNST_Pstd,           &
!         CNST_Rvap,           &
!         CNST_EPS,            &
!         CNST_LHV0,           &
!         CNST_TEM00
!    use mod_isotope, only: &
!       NISO,         &
!       DIFF_H18O,    &
!       DIFF_HDO,     &
!       WISO_eqfract, &
!       ALP0
!    implicit none
!    !
!    integer, intent(in)     :: ijdim
!    integer, intent(in)     :: kdim
!    integer, intent(in)     :: nqmax
!    !
!    real(RP), intent(in)    :: rho    (ijdim,kdim)
!    real(RP), intent(in)    :: tem    (ijdim,kdim)
!    real(RP), intent(in)    :: pre    (ijdim,kdim)
!    real(RP), intent(in)    :: qsatl  (ijdim,kdim)
!    real(RP), intent(in)    :: q      (ijdim,kdim,nqmax)
!    real(RP), intent(out)   :: qe_old (ijdim,kdim,NISO)
!    !
!    !
!    real(RP) :: temp, temc, rtem
!    real(RP) :: Ka, Kd, Glv, Giv, Sliq
!    real(RP) :: PRE00, EPS, Rvap, LHV0, TEM00
!    !
!    real(RP) :: ALP(NISO), ALPI(NISO), ALPS(NISO), AEVP(NISO)
!    real(RP) :: Sleff, ALPKL(NISO)
!    !
!    integer  :: ij, k
!    !
!    PRE00 = CNST_Pstd
!    EPS   = CNST_EPS
!    Rvap  = CNST_Rvap
!    LHV0  = CNST_LHV0
!    TEM00 = CNST_TEM00
!    PRE00 = CNST_Pstd
!
!    do k  = 1, kdim
!       do ij = 1, ijdim
!          temp        = tem(ij,k)
!          temc        = temp - TEM00
!          rtem        = 1.0_RP / max( temp, 10.0_RP )
!          Ka          = ( Ka0 + dKa_dT * temc )
!          Kd          = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
!          Glv         = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*rho(ij,k)*qsatl(ij,k)) )
!          Sliq        = q(ij,k,I_QV) / max( qsatl(ij,k), EPS )
!          call WISO_eqfract( temp, ALP0, ALP, ALPI )
!          Sleff = (1.0_RP - Glv * (1.0_RP - Sliq ** (-1.0_RP)) ) ** (-1.0_RP)
!          ALPKL(1) = ALP(1) / ( 1.0_RP + ( ALP(1) *DIFF_H18O-1.0_RP ) * (1.0_RP-Sleff**(-1.0_RP) ))
!          ALPKL(2) = ALP(2) / ( 1.0_RP + ( ALP(2) *DIFF_HDO -1.0_RP ) * (1.0_RP-Sleff**(-1.0_RP) ))
!          qe_old(ij,k,1) = q(ij,k,I_QV) + ALPKL(1)*q(ij,k,I_QC)
!          qe_old(ij,k,2) = q(ij,k,I_QV) + ALPKL(2)*q(ij,k,I_QC)
!       enddo
!    enddo
!
!  end subroutine isotope_condensation_pre

  !-----------------------------------------------------------------------------
  subroutine isotope_condensation_liq( &
       ijdim, &
       kdim,  &
       nqmax, &
       rho,   &
       rhog,  &
       rhogq, &
       tem,   &
       pre,   &
       q      )
    use mod_runconf, only: &
       I_QV,    &
       I_QC,    &
       ISO1_QV, &
       ISO2_QV, &
       ISO1_QC, &
       ISO2_QC
    use mod_isotope, only: &
       NISO,         &
       WISO_eqfract, &
       ALP0
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: nqmax
    real(RP), intent(in)    :: rho  (ijdim,kdim)
    real(RP), intent(in)    :: rhog (ijdim,kdim)
    real(RP), intent(inout) :: rhogq(ijdim,kdim,nqmax)
    real(RP), intent(in)    :: tem  (ijdim,kdim)
    real(RP), intent(in)    :: pre  (ijdim,kdim)
    real(RP), intent(inout) :: q    (ijdim,kdim,nqmax)

    real(RP) :: ALP (NISO)
    real(RP) :: ALPI(NISO)
    real(RP) :: qvci(NISO)
    real(RP) :: temp, qv, qc

    integer  :: ij, k
    !---------------------------------------------------------------------------

    do k  = 1, kdim
    do ij = 1, ijdim
       temp    = tem(ij,k)
       qv      = q(ij,k,I_QV)
       qc      = q(ij,k,I_QC)
       qvci(1) = q(ij,k,ISO1_QV) + q(ij,k,ISO1_QC)
       qvci(2) = q(ij,k,ISO2_QV) + q(ij,k,ISO2_QC)

       if ( qc > 0.0_RP ) then
          call WISO_eqfract( temp, real(ALP0,RP_G), ALP, ALPI )
          q(ij,k,ISO1_QC) = qvci(1) * ALP(1)*qc / ( qv + ALP(1)*qc )
          q(ij,k,ISO2_QC) = qvci(2) * ALP(2)*qc / ( qv + ALP(2)*qc )
          q(ij,k,ISO1_QV) = qvci(1) - q(ij,k,ISO1_QC)
          q(ij,k,ISO2_QV) = qvci(2) - q(ij,k,ISO2_QC)
       else
          q(ij,k,ISO1_QC) = 0.0_RP
          q(ij,k,ISO2_QC) = 0.0_RP
          q(ij,k,ISO1_QV) = qvci(1)
          q(ij,k,ISO2_QV) = qvci(2)
       endif
       rhogq(ij,k,ISO1_QV) = q(ij,k,ISO1_QV) * rhog(ij,k)
       rhogq(ij,k,ISO2_QV) = q(ij,k,ISO2_QV) * rhog(ij,k)
       rhogq(ij,k,ISO1_QC) = q(ij,k,ISO1_QC) * rhog(ij,k)
       rhogq(ij,k,ISO2_QC) = q(ij,k,ISO2_QC) * rhog(ij,k)
    enddo
    enddo

    return
  end subroutine isotope_condensation_liq

  !-----------------------------------------------------------------------------
!  subroutine isotope_condensation_all( &
!     ijdim,  &  !--[IN]
!     kdim,   &  !--[IN]
!     nqmax,  &  !--[IN]
!     rho,    &  !--[IN]
!     rhog,   &  !--[IN]
!     rhogq,  &  !--[INOUT]
!     tem,    &  !--[IN]
!     pre,    &  !--[IN]
!     qe_old, &  !--[IN]
!     qv_old, &  !--[IN]
!     qc_old, &  !--[IN]
!     qi_old, &  !--[IN]
!     q       )  !--[INOUT]
!    use mod_adm, &
!       only: ADM_LOG_FID
!    use mod_runconf, only: &
!         I_QV,    &
!         I_QC,    &
!         I_QI,    &
!         ISO1_QV, &
!         ISO2_QV, &
!         ISO1_QC, &
!         ISO2_QC, &
!         ISO1_QI, &
!         ISO2_QI
!    use mod_cnst, only: &
!         CNST_Pstd,           &
!         CNST_Rvap,           &
!         CNST_EPS,            &
!         CNST_LHV0,           &
!         CNST_LHS0,           &
!         CNST_TEM00
!    use mod_satadjust, only: &
!       SATURATION_psat_liq,  &
!       SATURATION_psat_ice
!    use mod_isotope, only: &
!       NISO,         &
!       DIFF_H18O,    &
!       DIFF_HDO,     &
!       WISO_eqfract, &
!       ALP0
!    implicit none
!    !
!    integer, intent(in)     :: ijdim
!    integer, intent(in)     :: kdim
!    integer, intent(in)     :: nqmax
!    !
!    real(RP), intent(in)    :: rho    (ijdim,kdim)
!    real(RP), intent(in)    :: rhog   (ijdim,kdim)
!    real(RP), intent(inout) :: rhogq  (ijdim,kdim,nqmax)
!    real(RP), intent(in)    :: tem    (ijdim,kdim)
!    real(RP), intent(in)    :: pre    (ijdim,kdim)
!    real(RP), intent(in)    :: qe_old (ijdim,kdim,NISO)
!    real(RP), intent(in)    :: qv_old (ijdim,kdim)
!    real(RP), intent(in)    :: qc_old (ijdim,kdim)
!    real(RP), intent(in)    :: qi_old (ijdim,kdim)
!    real(RP), intent(inout) :: q      (ijdim,kdim,nqmax)
!    !
!    real(RP) :: qsatl (ijdim,kdim)
!    real(RP) :: qsati (ijdim,kdim)
!    real(RP) :: psatl (ijdim,kdim)
!    real(RP) :: psati (ijdim,kdim)
!    !
!    real(RP) :: temp, temc, rtem
!    real(RP) :: Ka, Kd, Glv, Giv, Sliq, Sice
!    real(RP) :: PRE00, EPS, Rvap, LHV0, LHS0, TEM00
!    !
!    real(RP) :: ALP(NISO), ALPI(NISO), ALPS(NISO), AEVP(NISO)
!    real(RP) :: Sleff, Sieff, ALPKL(NISO), ALPKS(NISO)
!    real(RP) :: qv, qc, qi, qc_t, qi_t
!    real(RP) :: qvi(NISO), qci(NISO), qii(NISO)
!    real(RP) :: qii_t(NISO)
!    real(RP) :: qe(NISO), qvci(NISO), qsumi(NISO)
!    !
!    integer  :: ij, k, iso
!    !
!    PRE00 = CNST_Pstd
!    EPS   = CNST_EPS
!    Rvap  = CNST_Rvap
!    LHV0  = CNST_LHV0
!    LHS0  = CNST_LHS0
!    TEM00 = CNST_TEM00
!    PRE00 = CNST_Pstd
!
!    call SATURATION_psat_liq( ijdim, kdim, tem(:,:), psatl(:,:) )
!    call SATURATION_psat_ice( ijdim, kdim, tem(:,:), psati(:,:) )
!
!    do k  = 1, kdim
!       do ij = 1, ijdim
!          qv = max(q(ij,k,I_QV), 0.0_RP)
!          qc = max(q(ij,k,I_QC), 0.0_RP)
!          qi = max(q(ij,k,I_QI), 0.0_RP)
!
!          qc_t   = qc - qc_old(ij,k)
!          qi_t   = qi - qi_old(ij,k)
!
!          !if ( abs(qc_t) == 0.0_RP .and. abs(qi_t) == 0.0_RP ) then
!          !   !! do not update for the case of qc_t=0 and qi_t=0
!          !   cycle
!          !endif
!
!          temp = tem(ij,k)
!          temc = temp - TEM00
!          rtem = 1.0_RP / max( temp, 10.0_RP )
!          Ka   = ( Ka0 + dKa_dT * temc )
!          Kd   = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
!          qsatl(ij,k) = psatl(ij,k) / ( rho(ij,k) * Rvap * tem(ij,k) )
!          qsati(ij,k) = psati(ij,k) / ( rho(ij,k) * Rvap * tem(ij,k) )
!
!          Glv  = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*rho(ij,k)*qsatl(ij,k)) )
!          Giv  = 1.0_RP / ( LHS0/(Ka*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*rho(ij,k)*qsati(ij,k)) )
!          Sliq = qv / max( qsatl(ij,k), EPS )
!          Sice = qv / max( qsati(ij,k), EPS )
!
!          qvi(1) = max(q(ij,k,ISO1_QV), 0.0_RP)
!          qvi(2) = max(q(ij,k,ISO2_QV), 0.0_RP)
!          qci(1) = max(q(ij,k,ISO1_QC), 0.0_RP)
!          qci(2) = max(q(ij,k,ISO2_QC), 0.0_RP)
!          qii(1) = max(q(ij,k,ISO1_QI), 0.0_RP)
!          qii(2) = max(q(ij,k,ISO2_QI), 0.0_RP)
!
!          call WISO_eqfract( temp, ALP0, ALP, ALPI )
!          Sleff = (1.0_RP - Glv * (1.0_RP - Sliq ** (-1.0_RP)) ) ** (-1.0_RP)
!          Sieff = (1.0_RP - Giv * (1.0_RP - Sice ** (-1.0_RP)) ) ** (-1.0_RP)
!          ALPKL(1) = ALP(1)  / ( 1.0_RP + ( ALP(1) *DIFF_H18O-1.0_RP ) * (1.0_RP-Sleff**(-1.0_RP) ))
!          ALPKL(2) = ALP(2)  / ( 1.0_RP + ( ALP(2) *DIFF_HDO -1.0_RP ) * (1.0_RP-Sleff**(-1.0_RP) ))
!          ALPKS(1) = ALPI(1) / ( 1.0_RP + ( ALPI(1)*DIFF_H18O-1.0_RP ) * (1.0_RP-Sieff**(-1.0_RP) ))
!          ALPKS(2) = ALPI(2) / ( 1.0_RP + ( ALPI(2)*DIFF_HDO -1.0_RP ) * (1.0_RP-Sieff**(-1.0_RP) ))
!
!          do iso = 1, NISO
!             qsumi(iso) = qvi(iso) + qci(iso) + qvi(iso)
!             qe   (iso) = qv + ALPKL(iso)*qc
!          enddo
!
!          if ( qv+qc == 0.0_RP ) then
!             do iso = 1, NISO
!                qvi(iso) = 0.0_RP
!                qci(iso) = 0.0_RP
!                qii(iso) = qsumi(iso)
!             enddo
!          else
!             if ( qi_old(ij,k) > 0.0_RP ) then
!                if ( qi_t >= 0.0_RP ) then
!                   !! fractionation (from Blossy et al., 2010)
!                   do iso = 1, NISO
!                      qii_t(iso) = qii(iso) * ( (qi/qi_old(ij,k)) ** (ALPKS(iso)) -1.0_RP )
!                   enddo
!                else
!                   !! no fractionation
!                   do iso = 1, NISO
!                      qii_t(iso) = qii(iso) * (  qi/qi_old(ij,k) -1.0_RP )
!                   enddo
!                endif
!                !! calculate qvi and qci from mass balance
!                do iso = 1, NISO
!                   qvci (iso) = qvi(iso) + qci(iso) - qii_t(iso)
!                   qvi  (iso) = qv * qvci(iso) / qe  (iso)
!                   qci  (iso) = max( qvci(iso) - qvi (iso), 0.0_RP )
!                   qii  (iso) = max( qii (iso) + qii_t(iso), 0.0_RP )
!                enddo
!             else
!                !! in the case of qi_old=0.0
!                !! qii_t is calculated by qsumi-qvci
!                do iso = 1, NISO
!                   qvci (iso) = (qe(iso)+ALPKS(iso)*qi_t) * (qvi(iso)+qci(iso)) / qe_old(ij,k,iso)
!                   qvi  (iso) = qv * qvci (iso) / qe  (iso)
!                   qci  (iso) = max( qvci (iso) - qvi (iso), 0.0_RP )
!                   qii  (iso) = qsumi(iso) - qvci(iso)
!                enddo
!             endif
!          endif
!
!          ! Tracer update
!          q(ij,k,ISO1_QV) = qvi(1)
!          q(ij,k,ISO2_QV) = qvi(2)
!          q(ij,k,ISO1_QC) = qci(1)
!          q(ij,k,ISO2_QC) = qci(2)
!          q(ij,k,ISO1_QI) = qii(1)
!          q(ij,k,ISO2_QI) = qii(2)
!
!          rhogq(ij,k,ISO1_QV) = q(ij,k,ISO1_QV) * rhog(ij,k)
!          rhogq(ij,k,ISO2_QV) = q(ij,k,ISO2_QV) * rhog(ij,k)
!          rhogq(ij,k,ISO1_QC) = q(ij,k,ISO1_QC) * rhog(ij,k)
!          rhogq(ij,k,ISO2_QC) = q(ij,k,ISO2_QC) * rhog(ij,k)
!          rhogq(ij,k,ISO1_QI) = q(ij,k,ISO1_QI) * rhog(ij,k)
!          rhogq(ij,k,ISO2_QI) = q(ij,k,ISO2_QI) * rhog(ij,k)
!
!    enddo
!    enddo
!
!  end subroutine isotope_condensation_all

  !-----------------------------------------------------------------------------
  subroutine negative_filter( &
       ijdim, &
       rhog,  &
       rhoge, &
       rhogq, &
       rho,   &
       tem,   &
       pre,   &
       q,     &
       gsgam2 )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       opt_water_isotope, & ! [Add] M.Tanoue 20200907
       ISO_STR,           & ! [Add] M.Tanoue 20200907
       ISO_STR2,          & ! [Add] M.Tanoue 20200907
       ISO_END,           & ! [Add] M.Tanoue 20200907
       ISO1_QV,           & ! [Add] M.Tanoue 20200907
       ISO2_QV              ! [Add] M.Tanoue 20200907
    use mod_thrmdyn, only: &
       THRMDYN_cv, &
       THRMDYN_qd
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(inout) :: rhog  (ijdim,kdim)
    real(RP), intent(inout) :: rhoge (ijdim,kdim)
    real(RP), intent(inout) :: rhogq (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: pre   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: gsgam2(ijdim,kdim)

    real(RP) :: qd   (ijdim,kdim)
    real(RP) :: cva  (ijdim,kdim)
    real(RP) :: diffq(ijdim,kdim)
    real(RP) :: Rdry, Rvap

    real(RP) :: diffqi1(ijdim,kdim) ! [Add] M.Tanoue 20200907
    real(RP) :: diffqi2(ijdim,kdim) ! [Add] M.Tanoue 20200907

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(qd,cva,diffq) &
    !$acc pcopyin(tem,gsgam2)

    call PROF_rapstart('____MP_NSW6_NEGFILTER',3)

    Rdry = CNST_Rdry
    Rvap = CNST_Rvap

    !$omp parallel default(none),private(ij,k,nq,diffq), &
    !$omp shared(ijdim,kmin,kmax,NQW_STR,NQW_END,rhogq,rhog,q,rho,gsgam2,I_QV)

    !$acc kernels async(1) pcopy(diffq)
    !$omp do
    do k  = kmin, kmax
    do ij = 1, ijdim
       diffq(ij,k) = 0.0_RP
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    do nq = NQW_STR+1, NQW_END
       !$acc kernels async(1) pcopy(rhogq,diffq)
       !$omp do
       do k  = kmin, kmax
       do ij = 1, ijdim
          diffq(ij,k) = diffq(ij,k) + rhogq(ij,k,nq)
          ! remove negative value of hydrometeors (mass)
          rhogq(ij,k,nq) = max( rhogq(ij,k,nq), 0.0_RP )
          ! difference between before and after correction
          diffq(ij,k) = diffq(ij,k) - rhogq(ij,k,nq)
       enddo
       enddo
       !$omp end do
       !$acc end kernels
    enddo

    ! Compensate for the lack of hydrometeors by the water vapor

    !$acc kernels async(1) pcopy(rhogq) pcopyin(diffq)
    !$omp do
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) + diffq(ij,k)
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    !$acc kernels async(1) pcopy(rhog,rhogq,rho) pcopyin(diffq,gsgam2)
    !$omp do
    do k  = kmin, kmax
    do ij = 1, ijdim
       diffq(ij,k) = rhogq(ij,k,I_QV)
       ! remove negative value of water vapor (mass)
       rhogq(ij,k,I_QV) = max( rhogq(ij,k,I_QV), 0.0_RP )
       ! difference between before and after correction
       diffq(ij,k) = diffq(ij,k) - rhogq(ij,k,I_QV)

       ! Apply correction to total density
       rhog(ij,k) = rhog(ij,k) * ( 1.0_RP - diffq(ij,k) ) ! diffq is negative
       rho (ij,k) = rhog(ij,k) / gsgam2(ij,k)
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    !$acc kernels async(1) pcopy(q) pcopyin(rhogq,rhog)
    !$omp do
    do nq = NQW_STR, NQW_END
    do k  = kmin, kmax
    do ij = 1, ijdim
       !--- update mass concentration
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    !$omp end parallel

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )
    call THRMDYN_cv( ijdim, kdim, qd(:,:), q(:,:,:), cva(:,:) )

    !$acc kernels async(1) pcopy(rhoge,pre) pcopyin(rhog,rho,tem,qd,q,cva)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhoge,pre,rho,tem,rhog,qd,q,cva,Rdry,Rvap,I_QV)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhoge(ij,k) = tem(ij,k) * rhog(ij,k) * cva(ij,k)
       pre  (ij,k) = rho(ij,k) * ( qd(ij,k) * Rdry + q(ij,k,I_QV) * Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels


    ! isotope correction
    if ( opt_water_isotope ) then ! [Add] M.Tanoue 20200907
       do k  = kmin, kmax
       do ij = 1, ijdim
          diffqi1(ij,k) = 0.0_RP
          diffqi2(ij,k) = 0.0_RP
       enddo
       enddo

       do nq = ISO_STR+1, ISO_STR2-1
          do k  = kmin, kmax
          do ij = 1, ijdim
             diffqi1(ij,k)    = diffqi1(ij,k) + rhogq(ij,k,nq  )
             diffqi2(ij,k)    = diffqi2(ij,k) + rhogq(ij,k,nq+6)
             rhogq(ij,k,nq  ) = max( rhogq(ij,k,nq  ), 0.0_RP )
             rhogq(ij,k,nq+6) = max( rhogq(ij,k,nq+6), 0.0_RP )
             diffqi1(ij,k)    = diffqi1(ij,k) - rhogq(ij,k,nq  )
             diffqi2(ij,k)    = diffqi2(ij,k) - rhogq(ij,k,nq+6)
          enddo
          enddo
       enddo

       do k  = kmin, kmax
       do ij = 1, ijdim
          rhogq(ij,k,ISO1_QV) = rhogq(ij,k,ISO1_QV) + diffqi1(ij,k)
          rhogq(ij,k,ISO2_QV) = rhogq(ij,k,ISO2_QV) + diffqi2(ij,k)
       enddo
       enddo

       do nq = ISO_STR, ISO_END
       do k  = kmin, kmax
       do ij = 1, ijdim
          !--- update mass concentration
          q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
       enddo
       enddo
       enddo
    endif

    !$acc wait(1)
    call PROF_rapend  ('____MP_NSW6_NEGFILTER',3)

    !$acc end data
    return
  end subroutine negative_filter

  !-----------------------------------------------------------------------------
  subroutine Bergeron_param( &
       ijdim, &
       tem,   &
       a1,    &
       a2,    &
       ma2    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(out) :: a1 (ijdim,kdim)
    real(RP), intent(out) :: a2 (ijdim,kdim)
    real(RP), intent(out) :: ma2(ijdim,kdim)

    real(RP) :: a1_tab(32)
    real(RP) :: a2_tab(32)

    data a1_tab / 0.0001E-7_RP, 0.7939E-7_RP, 0.7841E-6_RP, 0.3369E-5_RP, 0.4336E-5_RP, &
                  0.5285E-5_RP, 0.3728E-5_RP, 0.1852E-5_RP, 0.2991E-6_RP, 0.4248E-6_RP, &
                  0.7434E-6_RP, 0.1812E-5_RP, 0.4394E-5_RP, 0.9145E-5_RP, 0.1725E-4_RP, &
                  0.3348E-4_RP, 0.1725E-4_RP, 0.9175E-5_RP, 0.4412E-5_RP, 0.2252E-5_RP, &
                  0.9115E-6_RP, 0.4876E-6_RP, 0.3473E-6_RP, 0.4758E-6_RP, 0.6306E-6_RP, &
                  0.8573E-6_RP, 0.7868E-6_RP, 0.7192E-6_RP, 0.6513E-6_RP, 0.5956E-6_RP, &
                  0.5333E-6_RP, 0.4834E-6_RP  /

    data a2_tab / 0.0100_RP, 0.4006_RP, 0.4831_RP, 0.5320_RP, 0.5307_RP, &
                  0.5319_RP, 0.5249_RP, 0.4888_RP, 0.3849_RP, 0.4047_RP, &
                  0.4318_RP, 0.4771_RP, 0.5183_RP, 0.5463_RP, 0.5651_RP, &
                  0.5813_RP, 0.5655_RP, 0.5478_RP, 0.5203_RP, 0.4906_RP, &
                  0.4447_RP, 0.4126_RP, 0.3960_RP, 0.4149_RP, 0.4320_RP, &
                  0.4506_RP, 0.4483_RP, 0.4460_RP, 0.4433_RP, 0.4413_RP, &
                  0.4382_RP, 0.4361_RP  /

    real(RP) :: temc
    integer  :: itemc
    real(RP) :: fact
    real(RP) :: TEM00

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data &
    !!$acc pcopyin(a1_tab,a2_tab)

    TEM00 = CNST_TEM00

    !$acc kernels async(1) pcopy(a1,a2,ma2) pcopyin(tem,a1_tab,a2_tab)
    !$omp parallel do default(none),private(ij,k,temc,itemc,fact), &
    !$omp shared(ijdim,kmin,kmax,a1,a2,ma2,tem,a1_tab,a2_tab,TEM00)
    do k  = kmin, kmax
    do ij = 1, ijdim
       temc  = min( max( tem(ij,k)-TEM00, -30.99_RP ), 0.0_RP ) ! 0C <= T  <  31C
       itemc = int( -temc ) + 1                                 ! 1  <= iT <= 31
       fact  = - ( temc + real(itemc-1,kind=RP) )

       a1 (ij,k) = ( 1.0_RP-fact ) * a1_tab(itemc  ) &
                 + (        fact ) * a1_tab(itemc+1)

       a2 (ij,k) = ( 1.0_RP-fact ) * a2_tab(itemc  ) &
                 + (        fact ) * a2_tab(itemc+1)

       ma2(ij,k) = 1.0_RP - a2(ij,k)

       a1 (ij,k) = a1(ij,k) * 1.E-3_RP**ma2(ij,k) ! [g->kg]
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !!$acc end data
    return
  end subroutine Bergeron_param

  !-----------------------------------------------------------------------------
  subroutine mp_nsw6_terminal_velocity( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       w,     &
       qc,    &
       qr,    &
       qi,    &
       qs,    &
       qg,    &
       vt_qc, &
       vt_qr, &
       vt_qi, &
       vt_qs, &
       vt_qg  )
    use mod_adm, only: &
       kmin => adm_kmin, &
       kmax => adm_kmax
    use mod_cnst, only: &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: rho  (ijdim,kdim)   !< density
    real(DP), intent(in)  :: tem  (ijdim,kdim)   !< temperature
    real(DP), intent(in)  :: w    (ijdim,kdim)   !< vertical velocity [m/s]    ![Add] 2016/06/04 W.Roh
    real(DP), intent(in)  :: qc   (ijdim,kdim)   !< mixing ratio of cloud
    real(DP), intent(in)  :: qr   (ijdim,kdim)   !< mixing ratio of rain
    real(DP), intent(in)  :: qi   (ijdim,kdim)   !< mixing ratio of ice
    real(DP), intent(in)  :: qs   (ijdim,kdim)   !< mixing ratio of snow
    real(DP), intent(in)  :: qg   (ijdim,kdim)   !< mixing ratio of graupel
    real(DP), intent(out) :: vt_qc(ijdim,kdim)   !< terminal velocity of cloud
    real(DP), intent(out) :: vt_qr(ijdim,kdim)   !< terminal velocity of rain
    real(DP), intent(out) :: vt_qi(ijdim,kdim)   !< terminal velocity of ice
    real(DP), intent(out) :: vt_qs(ijdim,kdim)   !< terminal velocity of snow
    real(DP), intent(out) :: vt_qg(ijdim,kdim)   !< terminal velocity of graupel

    real(DP) :: dens                             !< density
    real(DP) :: rho_fact                         !< density factor

    real(DP) :: RLMDr, RLMDr_dr
    real(DP) :: RLMDs, RLMDs_ds, RLMDs_1bs
    real(DP) :: RLMDg, RLMDg_dg

    !---< Roh and Satoh (2014) >---
    real(DP) :: tems, Xs2
    real(DP) :: MOMs_0bs, RMOMs_Vt
    real(DP) :: coef_at(4), coef_bt(4)
    real(DP) :: loga_, b_, nm

    real(DP) :: zerosw
    real(DP) :: EPS, TEM00

    integer  :: ij,k
    !---------------------------------------------------------------------------

    EPS   = 1.E-16_DP
    TEM00 = CNST_TEM00

    !$omp parallel do default(none), &
    !$omp private(ij,k,dens,rho_fact,RLMDr,RLMDr_dr,RLMDs,RLMDs_ds,RLMDs_1bs,RLMDg,RLMDg_dg, &
    !$omp         tems,Xs2,MOMs_0bs,RMOMs_Vt,coef_at,coef_bt,loga_,b_,nm,zerosw),            &
    !$omp shared(ijdim,kmin,kmax,vt_qc,vt_qr,vt_qi,vt_qs,vt_qg,rho,tem,qr,qi,qs,qg,          &
    !$omp        GAM_1br,GAM_1brdr,GAM_1bs,GAM_1bsds,GAM_1bg,GAM_1bgdg,coef_a,coef_b,        &
    !$omp        N0r,N0s,N0g,Ar,As,Ag,Cr,Cs,Cg,sw_roh2014,EPS,TEM00)
    do k  = kmin, kmax
    do ij = 1, ijdim

       dens     = rho(ij,k)
       rho_fact = sqrt( dens00 / dens )

       ! slope parameter lambda (Rain)
       zerosw = 0.5_DP - sign(0.5_DP, qr(ij,k) - 1.E-12_DP )
       RLMDr  = sqrt(sqrt( dens * qr(ij,k) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_DP - zerosw )
       RLMDr_dr  = sqrt( RLMDr )       ! **Dr

       ! slope parameter lambda (Snow)
       zerosw = 0.5_DP - sign(0.5_DP, qs(ij,k) - 1.E-12_DP )
       RLMDs  = sqrt(sqrt( dens * qs(ij,k) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_DP - zerosw )

       RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
       RLMDs_1bs = RLMDs**4 ! (1+Bs)

       MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
       RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds

       ! bimodal size distribution of snow
       zerosw = 0.5_DP - sign(0.5_DP, dens * qs(ij,k) - 1.E-12_DP )
       Xs2    = dens * qs(ij,k) / As * ( 1.0_DP - zerosw )

       tems       = min( -0.1_DP, tem(ij,k)-TEM00 )
       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
       coef_at(4) = coef_a(10)
       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
       coef_bt(4) = coef_b(10)
       ! 0 + Bs(=2) moment
       nm = 2.0_DP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
       MOMs_0bs = (        sw_roh2014 ) * 10.0_DP**loga_ * (Xs2+zerosw)**b_ * ( 1.0_DP - zerosw ) &
                + ( 1.0_DP-sw_roh2014 ) * MOMs_0bs
       ! Bs(=2) + Ds(=0.25) moment
       nm = 2.25_DP
       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )

       zerosw = 0.5_DP - sign(0.5_DP, abs(MOMs_0bs) - EPS )
       RMOMs_Vt = (        sw_roh2014 ) * 10.0_DP**loga_ * Xs2**b_ / ( MOMs_0bs + zerosw ) &
                + ( 1.0_DP-sw_roh2014 ) * RMOMs_Vt

       ! slope parameter lambda (Graupel)
       zerosw = 0.5_DP - sign(0.5_DP, qg(ij,k) - 1.E-12_DP )
       RLMDg  = sqrt(sqrt( dens * qg(ij,k) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_DP - zerosw )
       RLMDg_dg  = sqrt( RLMDg )       ! **Dg

       !---< terminal velocity >---

       vt_qc(ij,k) = 0.0_DP
       zerosw = 0.5_DP + sign(0.5_DP, qi(ij,k) - 1.E-8_DP )
       vt_qi(ij,k) = -3.29_DP * ( dens * qi(ij,k) * zerosw )**0.16_DP
       vt_qr(ij,k) = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
       vt_qs(ij,k) = -Cs * rho_fact * RMOMs_Vt
       vt_qg(ij,k) = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg

    enddo
    enddo

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,kmax,vt_qc,vt_qr,vt_qi,vt_qs,vt_qg)
    do ij = 1, ijdim
       vt_qc(ij,kmin-1) = 0.0_DP
       vt_qr(ij,kmin-1) = 0.0_DP
       vt_qi(ij,kmin-1) = 0.0_DP
       vt_qs(ij,kmin-1) = 0.0_DP
       vt_qg(ij,kmin-1) = 0.0_DP

       vt_qc(ij,kmax+1) = 0.0_DP
       vt_qr(ij,kmax+1) = 0.0_DP
       vt_qi(ij,kmax+1) = 0.0_DP
       vt_qs(ij,kmax+1) = 0.0_DP
       vt_qg(ij,kmax+1) = 0.0_DP
    enddo
    !$omp end parallel do

    return
  end subroutine mp_nsw6_terminal_velocity

  !-----------------------------------------------------------------------------
  subroutine mp_nsw6_effective_radius( &
       ijdim, kdim, kmin, kmax, &
       rho, ccn,                &
       factor_n_aie,            &
       qc,  qr,  qi,  qs,  qg,  &
       rec, rer, rei, res, reg  )
    use mod_cnst, only: &
       PI    => CNST_PI,   &
       rho_w => CNST_DWATR
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(DP), intent(in)  :: rho         (ijdim,kdim) ! air density [kg/m3]
    real(DP), intent(in)  :: ccn         (ijdim,kdim) ! cloud condensation nuclei [1/m3]
    real(DP), intent(in)  :: factor_n_aie(ijdim)      ! ant. indirect effect for NSW6 [add] C.Kodama 2017.11.27
    real(DP), intent(in)  :: qc          (ijdim,kdim)
    real(DP), intent(in)  :: qr          (ijdim,kdim)
    real(DP), intent(in)  :: qi          (ijdim,kdim)
    real(DP), intent(in)  :: qs          (ijdim,kdim)
    real(DP), intent(in)  :: qg          (ijdim,kdim)
    real(DP), intent(out) :: rec         (ijdim,kdim)
    real(DP), intent(out) :: rer         (ijdim,kdim)
    real(DP), intent(out) :: rei         (ijdim,kdim)
    real(DP), intent(out) :: res         (ijdim,kdim)
    real(DP), intent(out) :: reg         (ijdim,kdim)

    ! [Add] 2016/03/29 T.Seiki
    real(DP), parameter :: rho_i = 916.7_DP
    real(DP) :: nc
    real(DP) :: Xni
    real(DP) :: RLMDr, LMDs, LMDg
    real(DP) :: factor, OneThird, zerosw

    real(DP), parameter:: EPS = 1.E-16_DP

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(rec,rer,rei,res,reg) &
    !$acc pcopyin(rho,ccn,factor_n_aie,qc,qr,qi,qs,qg)

    OneThird = 1.0_DP / 3.0_DP

    ! [Add] 2016/03/29 T.Seiki
    factor = 3.0_DP / ( 4.0_DP * PI * rho_w )

    if ( OPT_INDIR ) then
       !$omp parallel do default(none),private(ij,k,nc), &
       !$omp shared(ijdim,kdim,rec,ccn,factor_n_aie,rho,qc,Nc_def,factor,OneThird)
       !$acc kernels pcopy(rec) pcopyin(ccn,factor_n_aie,rho,qc)
       do k  = 1, kdim
       do ij = 1, ijdim
!          nc = ccn(ij,k) ! [/m3]
          nc = max( ccn(ij,k)*factor_n_aie(ij), Nc_def*1.E+6_DP ) ! [mod] C.Kodama 2017.11.27
          rec(ij,k) = 1.1_DP * ( factor * rho(ij,k) * max( qc(ij,k), EPS ) / nc )**OneThird
          rec(ij,k) = min( 1.E-3_DP, max( 1.E-6_DP, rec(ij,k) ) ) ! limiter
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    else
       nc = Nc_def * 1.E+6_DP ! [/m3]

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,rec,rho,qc,Nc_def,factor,nc,OneThird)
       !$acc kernels pcopy(rec) pcopyin(rho,qc)
       do k  = 1, kdim
       do ij = 1, ijdim
!          nc = Nc_def * factor_n_aie(ij) * 1.E+6_DP ! [mod] C.Kodama 2017.11.27
          rec(ij,k) = 1.1_DP * ( factor * rho(ij,k) * max( qc(ij,k), EPS ) / nc )**OneThird
          rec(ij,k) = min( 1.E-3_DP, max( 1.E-6_DP, rec(ij,k) ) ) ! limiter
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    factor = 3.0_DP / ( 4.0_DP * PI * rho_i )

    if ( OPT_EXPLICIT_ICEGEN ) then
       !$omp parallel do default(none),private(ij,k,Xni), &
       !$omp shared(ijdim,kdim,rei,rho,qi,factor,OneThird)
       !$acc kernels pcopy(rei) pcopyin(rho,qi)
       do k  = 1, kdim
       do ij = 1, ijdim
          Xni = 5.38E+7_DP * ( rho(ij,k) * max( qi(ij,k), EPS ) )**0.75_DP
          Xni = min( 1.E+6_DP, max( 1.E+3_DP, Xni ) ) ! limiter

          rei(ij,k) = ( factor * rho(ij,k) * max( qi(ij,k), EPS ) / Xni )**OneThird
          rei(ij,k) = min( 1.E-3_DP, max( 1.E-6_DP, rei(ij,k) ) ) ! limiter
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    else
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,rei)
       !$acc kernels pcopy(rei)
       do k  = 1, kdim
       do ij = 1, ijdim
          rei(ij,k) = 40.E-6_DP
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    factor = 3.0_DP / ( 2.0_DP * PI * rho_i )

    !$omp parallel do default(none),private(ij,k,zerosw,RLMDr,LMDg), &
    !$omp shared(ijdim,kdim,rer,reg,rho,qr,qg,Ar,N0r,GAM_1br,Ag,N0g,GAM_1bg,factor)
    !$acc kernels pcopy(rer,reg) pcopyin(rho,qr,qg)
    do k  = 1, kdim
    do ij = 1, ijdim
       ! slope parameter lambda (Rain)
       zerosw = 0.5_DP - sign(0.5_DP, qr(ij,k) - 1.E-12_DP )
       RLMDr  = sqrt(sqrt( rho(ij,k) * qr(ij,k) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_DP - zerosw )
       ! slope parameter lambda (Graupel)
       zerosw = 0.5_DP - sign(0.5_DP, qg(ij,k) - 1.E-12_DP )
       LMDg   = sqrt(sqrt( Ag * N0g * GAM_1bg / ( rho(ij,k) * qg(ij,k) + zerosw ) )) * ( 1.0_DP - zerosw )

       rer(ij,k) = 1.5_DP * RLMDr
       reg(ij,k) = factor * rho(ij,k) * qg(ij,k) / N0g * LMDg**3
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    ! default value
    if ( Roh_flag ) then ! [add] WS.Roh 20141011

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,res,As)
       !$acc kernels pcopy(res)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! ratio of (pi/4) to 0.45 is a correction to circle with maximum dimension.
          ! 0.45 is derived from A=0.45*D**2 (approximateion of 0.2285*D**2 from Mitchell, 1996, JAS)
          res(ij,k) = ( 3.0_DP * As / ( PI * rho_i ) ) * ( PI * 0.25_DP / 0.45_DP ) ! = 125um when as=0.069
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

    else
       factor = 3.0_DP / ( 2.0_DP * PI * rho_i )

       !$omp parallel do default(none),private(ij,k,zerosw,LMDs), &
       !$omp shared(ijdim,kdim,res,rho,qs,As,N0s,GAM_1bs,factor)
       !$acc kernels pcopy(res) pcopyin(rho,qs)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! slope parameter lambda (Graupel)
          zerosw = 0.5_DP - sign(0.5_DP, qs(ij,k) - 1.E-12_DP )
          LMDs   = sqrt(sqrt( As * N0s * GAM_1bs / ( rho(ij,k) * qs(ij,k) + zerosw ) )) * ( 1.0_DP - zerosw )

          res(ij,k) = factor * rho(ij,k) * qs(ij,k) / N0s * LMDs**3
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

    endif

    !$acc end data
    return
  end subroutine mp_nsw6_effective_radius

  !-----------------------------------------------------------------------------
  !--- MINIAPP: ml2sl computes a column integral with GMTR_area and
  !--- VMTR_VOLUME (mod_gmtr + mod_vmtr + mod_oprt + mod_comm).  It has no
  !--- caller anywhere in the NSW6 chain.  Preserved, commented out. ---
!MINIAPP  subroutine ml2sl( &
!MINIAPP       ijdim, &
!MINIAPP       l,     &
!MINIAPP       rho,   &
!MINIAPP       v3d,   &
!MINIAPP       v2d    )
!MINIAPP    use mod_adm, only: &
!MINIAPP       k0   => ADM_KNONE, &
!MINIAPP       kdim => ADM_kall,  &
!MINIAPP       ADM_kmin,          &
!MINIAPP       ADM_kmax
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_area
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_VOLUME
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer,  intent(in)  :: ijdim
!MINIAPP    integer,  intent(in)  :: l
!MINIAPP    real(RP), intent(in)  :: rho(ijdim,kdim)
!MINIAPP    real(RP), intent(in)  :: v3d(ijdim,kdim)
!MINIAPP    real(RP), intent(out) :: v2d(ijdim,k0)
!MINIAPP
!MINIAPP    integer :: k
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    v2d(:,k0) = 0.0_RP
!MINIAPP    do k = ADM_kmin, ADM_kmax
!MINIAPP       v2d(:,k0) = v2d(:,k0) + rho(:,k) * v3d(:,k) * VMTR_VOLUME(:,k,l) / GMTR_area(:,l)
!MINIAPP    enddo
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine ml2sl

end module mod_mp_nsw6
