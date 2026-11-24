!-------------------------------------------------------------------------------
!> Module mod_mp_opt_bridge
!!
!! iso_c_binding bridge between NICAM mp_driver and the ultra-fused IEEE FP64
!! CUDA C++ / CPU operator suite (ANZIAM 2027).
!!
!! Runtime environment controls:
!!   NICAM_OPT_MODE        : 0 = original baseline (default if unset = 0)
!!                           1 = ENG           (engineering-only fused coalesced)
!!                           2 = ENG+MATH      (+ Roh dyadic + Laguerre-HH4 + Fubini)
!!                           3 = ENG+MATH+MONO (+ monotone backward-envelope departure)
!!   NICAM_COMPUTE_RCEFF   : 1 = include coalesced rceff column kernel (default 1)
!!                           0 = omit unused miniapp rceff diagnostic
!!   NICAM_COLLECT_COUNTERS: 1 = collect mathematical work counters (default 0)
!<
!-------------------------------------------------------------------------------
module mod_mp_opt_bridge
  use iso_c_binding
#ifdef _OPENACC
  use openacc
#endif
  use mod_precision, only: DP, RP => RP_MP, RP_PHY
  use mod_stdio
  implicit none
  private

  public :: NSW6Params_C
  public :: NSW6Counters_C
  public :: mp_opt_get_mode
  public :: mp_opt_reload_env
  public :: mp_opt_reset_counters
  public :: mp_opt_report_counters
  public :: mp_opt_driver_dispatch

  type, bind(C) :: NSW6Params_C
     integer(c_int) :: ijdim
     integer(c_int) :: kdim
     integer(c_int) :: kmin0
     integer(c_int) :: kmax0
     integer(c_int) :: nqmax
     integer(c_int) :: i_qv, i_qc, i_qr, i_qi, i_qs, i_qg

     integer(c_int) :: opt_mode
     integer(c_int) :: ein_is_simple
     integer(c_int) :: prcip_kin2ein
     integer(c_int) :: precip_4water
     integer(c_int) :: opt_indir
     integer(c_int) :: opt_sclw
     integer(c_int) :: opt_explicit_icegen
     integer(c_int) :: compute_rceff
     integer(c_int) :: is_top_tem
     integer(c_int) :: is_top_epl
     integer(c_int) :: is_btm_tem
     integer(c_int) :: is_btm_epl
     integer(c_int) :: mp_div_num
     integer(c_int) :: collect_counters

     real(c_double) :: dt
     real(c_double) :: PI
     real(c_double) :: EPS
     real(c_double) :: UNDEF
     real(c_double) :: GRAV
     real(c_double) :: Rdry
     real(c_double) :: Rvap
     real(c_double) :: CVdry
     real(c_double) :: CL
     real(c_double) :: LHV
     real(c_double) :: LHS
     real(c_double) :: LHF
     real(c_double) :: LHV0
     real(c_double) :: LHS0
     real(c_double) :: LHF0
     real(c_double) :: PSAT0
     real(c_double) :: TEM00
     real(c_double) :: PRE00
     real(c_double) :: rho_w
     real(c_double) :: CPovR_liq
     real(c_double) :: CVovR_liq
     real(c_double) :: LovR_liq
     real(c_double) :: CPovR_ice
     real(c_double) :: CVovR_ice
     real(c_double) :: LovR_ice
     real(c_double) :: CVW(6)
     real(c_double) :: W_VADV1D_MIN
     real(c_double) :: dtemp_criteria

     real(c_double) :: dens00
     real(c_double) :: N0r, N0s, N0g
     real(c_double) :: Ar, As, Ag
     real(c_double) :: Br, Bs, Bg
     real(c_double) :: Cr, Cs, Cg
     real(c_double) :: Dr, Ds, Dg
     real(c_double) :: GAM, GAM_2, GAM_3
     real(c_double) :: GAM_1br, GAM_2br, GAM_3br, GAM_3dr, GAM_6dr, GAM_1brdr, GAM_5dr_h
     real(c_double) :: GAM_1bs, GAM_2bs, GAM_3bs, GAM_3ds, GAM_1bsds, GAM_5ds_h
     real(c_double) :: GAM_1bg, GAM_3dg, GAM_1bgdg, GAM_5dg_h

     real(c_double) :: Eiw, Erw, Esw, Egw, Eri, Esi, Egi, Esr, Egr, Egs
     real(c_double) :: gamma_sacr, gamma_gacs, mi
     real(c_double) :: Nc_def
     real(c_double) :: beta_saut, gamma_saut, beta_gaut, gamma_gaut, qicrt_saut, qscrt_gaut
     real(c_double) :: Ka0, dKa_dT, Kd0, dKd_dT, nu0, dnu_dT
     real(c_double) :: f1r, f2r, f1s, f2s, f1g, f2g
     real(c_double) :: A_frz, B_frz
     real(c_double) :: mi40, mi50, vti50, Ri50
     real(c_double) :: Di_max, Di_a
     real(c_double) :: sw_kk2000, sw_roh2014, sw_expice, sw_sr2020
     real(c_double) :: sw_phillips2007, sw_beta_saut, sw_constVti, constVti
     real(c_double) :: Ds_cutoff, Dg_cutoff, CRIT_SCLW_W

     real(c_double) :: coef_a(10)
     real(c_double) :: coef_b(10)
     real(c_double) :: ln10
     real(c_double) :: a1_tab(32)
     real(c_double) :: a2_tab(32)

     real(c_double) :: rceff_coef_dgam
     real(c_double) :: rceff_coef_xf
  end type NSW6Params_C

  type, bind(C) :: NSW6Counters_C
     integer(c_long_long) :: total_cells
     integer(c_long_long) :: dry_clear_cells
     integer(c_long_long) :: snow_active_cells
     integer(c_long_long) :: sat_active_cells
     integer(c_long_long) :: sat_total_iters
     integer(c_long_long) :: n_exp_calls
     integer(c_long_long) :: n_pow_calls
     integer(c_long_long) :: n_log_calls
     integer(c_long_long) :: dry_precip_cols
     integer(c_long_long) :: precip_layers_touched
     integer(c_long_long) :: precip_crossings_full
     integer(c_long_long) :: precip_crossings_half
     integer(c_long_long) :: precip_top_overshoots
     integer(c_long_long) :: neg_q_before_postfilter
  end type NSW6Counters_C

  interface
     subroutine nsw6_opt_reset_counters(stream_ptr) bind(C, name="nsw6_opt_reset_counters")
       import :: c_ptr
       type(c_ptr), value :: stream_ptr
     end subroutine nsw6_opt_reset_counters

     subroutine nsw6_opt_get_counters(out_cnt, stream_ptr) bind(C, name="nsw6_opt_get_counters")
       import :: c_ptr, NSW6Counters_C
       type(NSW6Counters_C), intent(out) :: out_cnt
       type(c_ptr), value :: stream_ptr
     end subroutine nsw6_opt_get_counters

     subroutine nsw6_opt_driver_launch( &
          P, &
          rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
          rho, pre, tem, vx, vy, vz, w, q, &
          precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy, &
          gdcfrc, gdclw, &
          unccn, factor_n_aie, z, zh, &
          ix, iy, iz, jx, jy, jz, &
          phi, gam2, gam2h, gsgam2, gsgam2h, &
          GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact, &
          stream_ptr) bind(C, name="nsw6_opt_driver_launch")
       import :: c_ptr, c_double, NSW6Params_C
       type(NSW6Params_C), intent(in) :: P
       real(c_double), intent(inout) :: rhog(*), rhogvx(*), rhogvy(*), rhogvz(*), rhogw(*), rhoge(*), rhogq(*)
       real(c_double), intent(inout) :: rho(*), pre(*), tem(*), vx(*), vy(*), vz(*), w(*), q(*)
       real(c_double), intent(out)   :: precip(*), precip_lh(*), precip_rhoein(*), precip_rhophi(*), precip_rhokin(*), precip_energy(*)
       real(c_double), intent(inout) :: gdcfrc(*)
       real(c_double), intent(out)   :: gdclw(*)
       real(c_double), intent(in)    :: unccn(*), factor_n_aie(*), z(*), zh(*)
       real(c_double), intent(in)    :: ix(*), iy(*), iz(*), jx(*), jy(*), jz(*)
       real(c_double), intent(in)    :: phi(*), gam2(*), gam2h(*), gsgam2(*), gsgam2h(*)
       real(c_double), intent(in)    :: GRD_rdgz(*), GRD_rdgzh(*), GRD_afact(*), GRD_bfact(*), GRD_cfact(*), GRD_dfact(*)
       type(c_ptr), value :: stream_ptr
     end subroutine nsw6_opt_driver_launch
  end interface

  logical, save :: env_initialized = .false.
  integer, save :: cached_opt_mode = 0
  integer, save :: cached_compute_rceff = 1
  integer, save :: cached_collect_counters = 0

contains

  subroutine init_opt_env()
    character(len=32) :: val
    integer :: stat, ios, ival
    if (env_initialized) return
    env_initialized = .true.

    cached_opt_mode = 0
    call get_environment_variable('NICAM_OPT_MODE', val, status=stat)
    if (stat == 0 .and. len_trim(val) > 0) then
       read(val, *, iostat=ios) ival
       if (ios == 0) cached_opt_mode = ival
    endif

    cached_compute_rceff = 1
    call get_environment_variable('NICAM_COMPUTE_RCEFF', val, status=stat)
    if (stat == 0 .and. len_trim(val) > 0) then
       read(val, *, iostat=ios) ival
       if (ios == 0) cached_compute_rceff = ival
    endif

    cached_collect_counters = 0
    call get_environment_variable('NICAM_COLLECT_COUNTERS', val, status=stat)
    if (stat == 0 .and. len_trim(val) > 0) then
       read(val, *, iostat=ios) ival
       if (ios == 0) cached_collect_counters = ival
    endif
  end subroutine init_opt_env

  integer function mp_opt_get_mode()
    call init_opt_env()
    mp_opt_get_mode = cached_opt_mode
  end function mp_opt_get_mode

  subroutine mp_opt_reload_env()
    env_initialized = .false.
    call init_opt_env()
  end subroutine mp_opt_reload_env

  subroutine mp_opt_reset_counters()
    type(c_ptr) :: strm
    call init_opt_env()
    if (cached_opt_mode <= 0) return
#ifdef _OPENACC
    strm = acc_get_cuda_stream(1)
#else
    strm = c_null_ptr
#endif
    call nsw6_opt_reset_counters(strm)
  end subroutine mp_opt_reset_counters

  subroutine mp_opt_report_counters(label)
    character(len=*), intent(in) :: label
    type(NSW6Counters_C) :: cnt
    type(c_ptr) :: strm
    call init_opt_env()
    if (cached_opt_mode <= 0 .or. cached_collect_counters == 0) return
#ifdef _OPENACC
    strm = acc_get_cuda_stream(1)
#else
    strm = c_null_ptr
#endif
    call nsw6_opt_get_counters(cnt, strm)
    write(*,'(1x,A,A,I2,A)') '====== MATH WORK COUNTERS [', trim(label), cached_opt_mode, '] ======'
    write(*,'(1x,A,I14)') '--- total interior cells        : ', cnt%total_cells
    write(*,'(1x,A,I14)') '--- dry clear-sky skipped cells : ', cnt%dry_clear_cells
    write(*,'(1x,A,I14)') '--- snow-active cells (Xs2>0)   : ', cnt%snow_active_cells
    write(*,'(1x,A,I14)') '--- sat-adjust active cells     : ', cnt%sat_active_cells
    write(*,'(1x,A,I14)') '--- sat-adjust total iterations : ', cnt%sat_total_iters
    write(*,'(1x,A,I14)') '--- exp() evaluations           : ', cnt%n_exp_calls
    write(*,'(1x,A,I14)') '--- pow() evaluations           : ', cnt%n_pow_calls
    write(*,'(1x,A,I14)') '--- log() evaluations           : ', cnt%n_log_calls
    write(*,'(1x,A,I14)') '--- dry precip columns skipped  : ', cnt%dry_precip_cols
    write(*,'(1x,A,I14)') '--- precip layers touched       : ', cnt%precip_layers_touched
    write(*,'(1x,A,I14)') '--- trajectory crossings (full) : ', cnt%precip_crossings_full
    write(*,'(1x,A,I14)') '--- trajectory crossings (half) : ', cnt%precip_crossings_half
    write(*,'(1x,A,I14)') '--- top overshoots              : ', cnt%precip_top_overshoots
    write(*,'(1x,A,I14)') '--- neg-q cells pre-postfilter  : ', cnt%neg_q_before_postfilter
  end subroutine mp_opt_report_counters

  subroutine mp_opt_driver_dispatch( &
       ijdim, kdim, nqmax, dt, &
       rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
       rho, pre, tem, vx, vy, vz, w, q, &
       precip, precip_ISO1, precip_ISO2, &
       precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy, &
       GDCFRC, GDCLW, &
       unccn, factor_n_aie, z, zh, &
       ix, iy, iz, jx, jy, jz, &
       phi, gam2, gam2h, gsgam2, gsgam2h, &
       handled )
    use mod_adm, only: kmin => ADM_kmin, kmax => ADM_kmax
    use mod_grd, only: GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact
    use mod_cnst, only: &
       CNST_PI, CNST_EPS, CNST_UNDEF, CNST_GRAV, CNST_Rdry, CNST_Rvap, CNST_CVdry, &
       CNST_CL, CNST_LHV0, CNST_LHS0, CNST_LHF0, CNST_PSAT0, CNST_TEM00, CNST_Pstd, CNST_DWATR
    use mod_runconf, only: &
       EIN_TYPE, PRCIP_TRN_ECORRECT, MP_DIV_NUM, opt_water_isotope, &
       I_QV, I_QC, I_QR, I_QI, I_QS, I_QG, LHV, LHS, LHF, CVW, W_VADV1D_MIN
    use mod_satadjust, only: &
       CPovR_liq, CVovR_liq, LovR_liq, CPovR_ice, CVovR_ice, LovR_ice
    use mod_bndcnd, only: BNDCND_get_thermo_flags
    use mod_mp_nsw6, only: mp_nsw6_fill_opt_params
    implicit none

    integer,      intent(in)    :: ijdim, kdim, nqmax
    real(DP),     intent(in)    :: dt
    real(RP_PHY), target, intent(inout) :: rhog(ijdim,kdim), rhogvx(ijdim,kdim), rhogvy(ijdim,kdim)
    real(RP_PHY), target, intent(inout) :: rhogvz(ijdim,kdim), rhogw(ijdim,kdim), rhoge(ijdim,kdim)
    real(RP_PHY), target, intent(inout) :: rhogq(ijdim,kdim,nqmax)
    real(RP_PHY), target, intent(inout) :: rho(ijdim,kdim), pre(ijdim,kdim), tem(ijdim,kdim)
    real(RP_PHY), target, intent(inout) :: vx(ijdim,kdim), vy(ijdim,kdim), vz(ijdim,kdim), w(ijdim,kdim)
    real(RP_PHY), target, intent(inout) :: q(ijdim,kdim,nqmax)
    real(RP_PHY), target, intent(out)   :: precip(ijdim,2), precip_ISO1(ijdim,2), precip_ISO2(ijdim,2)
    real(RP_PHY), target, intent(out)   :: precip_lh(ijdim), precip_rhoein(ijdim), precip_rhophi(ijdim)
    real(RP_PHY), target, intent(out)   :: precip_rhokin(ijdim), precip_energy(ijdim)
    real(RP_PHY), target, intent(inout) :: GDCFRC(ijdim,kdim)
    real(RP_PHY), target, intent(out)   :: GDCLW(ijdim,kdim)
    real(RP_PHY), target, intent(in)    :: unccn(ijdim,kdim), factor_n_aie(ijdim)
    real(RP_PHY), target, intent(in)    :: z(ijdim,kdim), zh(ijdim,kdim)
    real(RP_PHY), target, intent(in)    :: ix(ijdim), iy(ijdim), iz(ijdim)
    real(RP_PHY), target, intent(in)    :: jx(ijdim), jy(ijdim), jz(ijdim)
    real(RP_PHY), target, intent(in)    :: phi(ijdim,kdim), gam2(ijdim,kdim), gam2h(ijdim,kdim)
    real(RP_PHY), target, intent(in)    :: gsgam2(ijdim,kdim), gsgam2h(ijdim,kdim)
    logical,      intent(out)   :: handled

    type(NSW6Params_C) :: P
    type(c_ptr) :: strm
    logical :: ok_nsw6, b_top_tem, b_top_epl, b_btm_tem, b_btm_epl
    integer :: ij

    handled = .false.
    call init_opt_env()
    if (cached_opt_mode <= 0) return

    ! Guard: require double precision, standard 6-tracer NSW6, no isotopes, MP_DIV_NUM==1, KIN2EIN
    if (RP /= DP .or. RP_PHY /= DP) return
    if (opt_water_isotope) return
    if (nqmax /= 6 .or. MP_DIV_NUM /= 1) return
    if (PRCIP_TRN_ECORRECT /= 'KIN2EIN') return

    P%ijdim = ijdim
    P%kdim  = kdim
    P%kmin0 = kmin - 1
    P%kmax0 = kmax - 1
    P%nqmax = nqmax
    P%i_qv  = I_QV - 1
    P%i_qc  = I_QC - 1
    P%i_qr  = I_QR - 1
    P%i_qi  = I_QI - 1
    P%i_qs  = I_QS - 1
    P%i_qg  = I_QG - 1

    P%opt_mode         = cached_opt_mode
    P%ein_is_simple    = merge(1, 0, EIN_TYPE == 'SIMPLE')
    P%prcip_kin2ein    = 1
    P%compute_rceff    = cached_compute_rceff
    P%mp_div_num       = MP_DIV_NUM
    P%collect_counters = cached_collect_counters

    call BNDCND_get_thermo_flags(b_top_tem, b_top_epl, b_btm_tem, b_btm_epl)
    P%is_top_tem = merge(1, 0, b_top_tem)
    P%is_top_epl = merge(1, 0, b_top_epl)
    P%is_btm_tem = merge(1, 0, b_btm_tem)
    P%is_btm_epl = merge(1, 0, b_btm_epl)

    P%dt     = dt / real(MP_DIV_NUM, kind=DP)
    P%PI     = CNST_PI
    P%EPS    = CNST_EPS
    P%UNDEF  = CNST_UNDEF
    P%GRAV   = CNST_GRAV
    P%Rdry   = CNST_Rdry
    P%Rvap   = CNST_Rvap
    P%CVdry  = CNST_CVdry
    P%CL     = CNST_CL
    P%LHV    = LHV
    P%LHS    = LHS
    P%LHF    = LHF
    P%LHV0   = CNST_LHV0
    P%LHS0   = CNST_LHS0
    P%LHF0   = CNST_LHF0
    P%PSAT0  = CNST_PSAT0
    P%TEM00  = CNST_TEM00
    P%PRE00  = CNST_Pstd
    P%rho_w  = CNST_DWATR
    P%CPovR_liq = CPovR_liq
    P%CVovR_liq = CVovR_liq
    P%LovR_liq  = LovR_liq
    P%CPovR_ice = CPovR_ice
    P%CVovR_ice = CVovR_ice
    P%LovR_ice  = LovR_ice
    P%CVW(1:6)  = CVW(1:6)
    P%W_VADV1D_MIN   = W_VADV1D_MIN
    P%dtemp_criteria = 10.0_DP**(-(precision(1.0_DP)+1)/2)

    call mp_nsw6_fill_opt_params( &
         P%precip_4water, P%opt_indir, P%opt_sclw, P%opt_explicit_icegen, &
         P%dens00, P%N0r, P%N0s, P%N0g, P%Ar, P%As, P%Ag, P%Br, P%Bs, P%Bg, &
         P%Cr, P%Cs, P%Cg, P%Dr, P%Ds, P%Dg, P%GAM, P%GAM_2, P%GAM_3, &
         P%GAM_1br, P%GAM_2br, P%GAM_3br, P%GAM_3dr, P%GAM_6dr, P%GAM_1brdr, P%GAM_5dr_h, &
         P%GAM_1bs, P%GAM_2bs, P%GAM_3bs, P%GAM_3ds, P%GAM_1bsds, P%GAM_5ds_h, &
         P%GAM_1bg, P%GAM_3dg, P%GAM_1bgdg, P%GAM_5dg_h, &
         P%Eiw, P%Erw, P%Esw, P%Egw, P%Eri, P%Esi, P%Egi, P%Esr, P%Egr, P%Egs, &
         P%gamma_sacr, P%gamma_gacs, P%mi, P%Nc_def, &
         P%beta_saut, P%gamma_saut, P%beta_gaut, P%gamma_gaut, P%qicrt_saut, P%qscrt_gaut, &
         P%Ka0, P%dKa_dT, P%Kd0, P%dKd_dT, P%nu0, P%dnu_dT, &
         P%f1r, P%f2r, P%f1s, P%f2s, P%f1g, P%f2g, P%A_frz, P%B_frz, &
         P%mi40, P%mi50, P%vti50, P%Ri50, P%Di_max, P%Di_a, &
         P%sw_kk2000, P%sw_roh2014, P%sw_expice, P%sw_sr2020, &
         P%sw_phillips2007, P%sw_beta_saut, P%sw_constVti, P%constVti, &
         P%Ds_cutoff, P%Dg_cutoff, P%CRIT_SCLW_W, &
         P%coef_a, P%coef_b, P%ln10, P%a1_tab, P%a2_tab, &
         P%rceff_coef_dgam, P%rceff_coef_xf, ok_nsw6 )

    if (.not. ok_nsw6) return

    ! Launch directly on OpenACC device pointers (zero copies, zero pcreate, zero waits!)
#ifdef _OPENACC
    strm = acc_get_cuda_stream(1)
    !$acc host_data use_device( &
    !$acc   rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
    !$acc   rho, pre, tem, vx, vy, vz, w, q, &
    !$acc   precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy, &
    !$acc   GDCFRC, GDCLW, unccn, factor_n_aie, z, zh, &
    !$acc   ix, iy, iz, jx, jy, jz, phi, gam2, gam2h, gsgam2, gsgam2h, &
    !$acc   GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact )
    call nsw6_opt_driver_launch( &
         P, &
         rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
         rho, pre, tem, vx, vy, vz, w, q, &
         precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy, &
         GDCFRC, GDCLW, &
         unccn, factor_n_aie, z, zh, &
         ix, iy, iz, jx, jy, jz, &
         phi, gam2, gam2h, gsgam2, gsgam2h, &
         GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact, &
         strm )
    !$acc end host_data
#else
    strm = c_null_ptr
    do ij = 1, ijdim
       precip_ISO1(ij,1) = 0.0_DP
       precip_ISO1(ij,2) = 0.0_DP
       precip_ISO2(ij,1) = 0.0_DP
       precip_ISO2(ij,2) = 0.0_DP
    enddo
    call nsw6_opt_driver_launch( &
         P, &
         rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
         rho, pre, tem, vx, vy, vz, w, q, &
         precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy, &
         GDCFRC, GDCLW, &
         unccn, factor_n_aie, z, zh, &
         ix, iy, iz, jx, jy, jz, &
         phi, gam2, gam2h, gsgam2, gsgam2h, &
         GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact, &
         strm )
#endif

    handled = .true.
  end subroutine mp_opt_driver_dispatch

end module mod_mp_opt_bridge
