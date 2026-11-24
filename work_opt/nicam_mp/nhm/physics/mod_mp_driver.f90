!-------------------------------------------------------------------------------
!
!+  Microphysics driver module
!
!-------------------------------------------------------------------------------
module mod_mp_driver
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module sets up and load the microphysics process.
  !
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17  Imported from Satoh's nonhydro model.
  !      0.03      04-05-31  Add lin scheme.
  !                          Modify the interface of sub[mp_driver].
  !                05-10-07  M.Satoh: Add g98w5 scheme.
  !                05-11-16  M.Satoh: zero set: GDCLW, GDCFRC
  !                06-01-26  H.Tomita: Add some mp modules.
  !                07-07-24  K.Suzuki: implementing SPRINTARS
  !                          1. Add unccn to argument for aerosol indirect effect
  !                          2. Add precip_flux to argument for use in SPRINTARS
  !                07-11-22  T.Mitsui: Effective Radius is calculated in rd_driver => mp_driver
  !                08-02-19  T.Mitsui: Add argument rceff_cld and l_region
  !                08-03-06  H.Tomita: Add WSM3 scheme.
  !                08-04-12  T.Mitsui: Add NDW6 scheme.
  !                08-04-14  H.Tomita: Add WSM6 scheme.
  !                08-04-22  H.Tomita: Saturation adjustment process is moved
  !                                    to mod_mp_nsw6 for NSW6 scheme.
  !                08-04-23  H.Tomita: Add the option for number of calculation
  !                                    of microphsics.
  !                08-05-30  T.Mitsui: Mod rceff => re_liquid, Add re_solid
  !                                    and del rceff_tmp
  !                09-04-14  T.Mitsui: add subroutine mp_terminal_velocity,
  !                                    mp_diag_volume and mp_effective_radius,
  !                                    add ndw6_1m for 1-Moment ver. of NDW6
  !                09-04-15  T.Mitsui: Add option to skip effective_radius and diag_volume
  !                09-08-05  T.Mitsui: cfra_cum is adm_knone
  !                09-08-18  T.Mitsui: add argument frhoge_rad and qke in mp_driver and
  !                                    there remains bug in the change of 09-08-05
  !                09-09-1X   S.Iga : Merge '09/01/19 [Add] M.Satoh'
  !                09-09-18   T.Mitsui: recovery of 09-09-03 T.Mitsui
  !                11-08-29   A.Noda: bugfix (precip_flux=0 for LSC)
  !                12-02-01   T.Seiki: add argument of tracers in precipitation for output
  !                12-07-23  H.Yashiro : [add] Water Isotope by K.Yoshimura
  !                12-10-18  T.Seiki   : recovery of precip_flux(11-08-29) for SPRINTARS
  !                16-03-29  T.Seiki   : add w in subroutine mp_effective_radius
  !                16-05-13  T.Seiki   : fix cloud volume used for radiation (NSW6)
  !                17-12-27  Y.Sato    : add spectral bin microphysics (S10)
  !                20-08-28  M.Tanoue  : add cloud water tracers (LSC)
  !                20-10-30  M.Tanoue  : add isotope tracers (NSW6)
  !      ----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision, only: &
     SP, DP, RP => RP_MP, RP_PHY
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: mp_init
  public :: mp_driver
  public :: mp_terminal_velocity
!MINIAPP  public :: mp_diag_volume
!MINIAPP  public :: mp_effective_radius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: MP_TYPE = 'NONE'

  logical,  private :: opt_radius_explicit = .false.
  logical,  private :: opt_volume_explicit = .true.

  real(RP), private :: const_re_cld        =   8.E-6_RP
  real(RP), private :: const_re_rain       = 500.E-6_RP
  real(RP), private :: const_re_ice        =  40.E-6_RP
  real(RP), private :: const_re_snow       = 500.E-6_RP
  real(RP), private :: const_re_graupel    = 500.E-6_RP
  real(RP), private :: const_re_hail       = 500.E-6_RP

  real(RP), private :: TSICE               = 273.15_RP ! temperature of ice/liq = 0 (for tentative use)
  real(RP), private :: TWICE               = 258.15_RP ! temperature of ice/liq = 1 (for tentative use)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine mp_init( MP_TYPE_in )
    use mod_adm, only: &
       ADM_proc_stop
!--- MINIAPP: only MP_TYPE='NSW6' is compiled into the mini-app; the other
!--- schemes and the cloudvar/radvar diagnostic containers (which pull in
!--- mod_comm) are guarded out with -DMINIAPP. ---
#ifndef MINIAPP
    use mod_mp_lsc, only: &
       mp_lsc_init
    use mod_mp_kessler, only: &
       mp_kessler_init
    use mod_mp_g98, only: &
       mp_g98_init
    use mod_mp_lin, only: &
       mp_lin_init
    use mod_mp_wsm6, only: &
       mp_wsm6_init
    use mod_mp_nsw5, only: &
       mp_nsw5_init
#endif
    use mod_mp_nsw6, only: &
       mp_nsw6_init
#ifndef MINIAPP
    use mod_mp_ndw6, only: &
       mp_ndw6_init
    use mod_mp_s10, only: &
       mp_s10_init
    use mod_mp_lsc_miroc6, only: &
       mp_lsc_miroc6_init
    use mod_cloudvar, only: &
       cloudvar_setup
#endif
    implicit none

    character(len=*), intent(in) :: MP_TYPE_in

    namelist / nm_mp_driver_init / &
       opt_radius_explicit, &
       opt_volume_explicit, &
       const_re_cld,        &
       const_re_rain,       &
       const_re_ice,        &
       const_re_snow,       &
       const_re_graupel,    &
       const_re_hail,       &
       TSICE,               &
       TWICE

    integer :: ierr
    !---------------------------------------------------------------------------

    MP_TYPE = trim(MP_TYPE_in)

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[Cloud Microphysics Driver]/Category[nhm physics]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_driver_init,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** nm_mp_driver_init is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,         *) 'xxx Not appropriate names in namelist nm_mp_driver_init. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist nm_mp_driver_init. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=nm_mp_driver_init)

    write(IO_FID_LOG,*)

    select case(MP_TYPE)
!--- MINIAPP: all non-NSW6 branches are guarded out (see above). ---
#ifndef MINIAPP
    case('NONE')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = Not used'
    case('LSC')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = Large Scale Condensation'
       call mp_lsc_init
       call cloudvar_setup
    case('LSC_MIROC6')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = Large Scale Condensation (MIROC6)'
       call mp_lsc_miroc6_init
       call cloudvar_setup
    case('KESSLER')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = KESSLER'
       call mp_kessler_init
       call cloudvar_setup
    case('G98')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = G98'
       call mp_g98_init
       call cloudvar_setup
    case('LIN')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = LIN'
       call mp_lin_init
       call cloudvar_setup
    case('WSM6')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = WSM6'
       call mp_wsm6_init
       call cloudvar_setup
    case('NSW5')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = NSW5'
       call mp_nsw5_init
       call cloudvar_setup
#endif
    case('NSW6')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = NSW6'
       call mp_nsw6_init
!--- MINIAPP: cloudvar_setup allocates the MP->CH/AE diagnostic container and
!--- needs mod_comm; the mini-app validates none of those fields. ---
#ifndef MINIAPP
       call cloudvar_setup
#endif
#ifndef MINIAPP
    case('NDW6')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = NDW6'
       call mp_ndw6_init
       call cloudvar_setup
    case('S10')
       write(IO_FID_LOG,*) '--- Cloud Microphysics Type = S10'
       call mp_s10_init
       call cloudvar_setup
#endif
    case default
       write(*,*) 'xxx Not appropriate Cloud Microphysics Type. MP_TYPE = ', trim(MP_TYPE)
       call ADM_proc_stop
    end select

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A,L1)') '--- Calculate optical volume of hydrometeors?   : ', opt_volume_explicit
    write(IO_FID_LOG,'(1x,A,L1)') '--- Calculate effective radius of hydrometeors? : ', opt_radius_explicit
    if ( .NOT. opt_radius_explicit ) then
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of liquid cloud = ', const_re_cld
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of rain         = ', const_re_rain
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of ice cloud    = ', const_re_ice
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of snow         = ', const_re_snow
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of graupel      = ', const_re_graupel
       write(IO_FID_LOG,'(1x,A,ES15.8)') '--- Constant Re of hail         = ', const_re_hail
    endif

    return
  end subroutine mp_init

  !-----------------------------------------------------------------------------
  subroutine mp_driver( &
       ijdim,         &
       kdim,          &
       l_region,      &
       RP_PHY_rhog,          &
       RP_PHY_rhogvx,        &
       RP_PHY_rhogvy,        &
       RP_PHY_rhogvz,        &
       RP_PHY_rhogw,         &
       RP_PHY_rhoge,         &
       RP_PHY_rhogq,         &
       RP_PHY_rho,           &
       RP_PHY_pre,           &
       RP_PHY_tem,           &
       RP_PHY_vx,            &
       RP_PHY_vy,            &
       RP_PHY_vz,            &
       RP_PHY_w,             &
       RP_PHY_q,             &
       RP_PHY_precip,        &
       RP_PHY_precip_ISO1,   &
       RP_PHY_precip_ISO2,   &
       RP_PHY_precip_lh,     &
       RP_PHY_precip_rhoein, &
       RP_PHY_precip_rhophi, &
       RP_PHY_precip_rhokin, &
       RP_PHY_precip_energy, &
       RP_PHY_GDCFRC,        &
       RP_PHY_GDCLW,         &
       RP_PHY_CBMFX,         &
       RP_PHY_unccn,         &
       RP_PHY_factor_n_aie,  &
       RP_PHY_rhogqke,       &
       RP_PHY_z,             &
       RP_PHY_zh,            &
       RP_PHY_ix,            &
       RP_PHY_iy,            &
       RP_PHY_iz,            &
       RP_PHY_jx,            &
       RP_PHY_jy,            &
       RP_PHY_jz,            &
       RP_PHY_phi,           &
       RP_PHY_gam2,          &
       RP_PHY_gam2h,         &
       RP_PHY_gsgam2,        &
       RP_PHY_gsgam2h,       &
       dt,            &
       ct             )
    use mod_cnst, only: &
       CNST_UNDEF
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       MP_DIV_NUM,        &
       NCRF,              &
       NQW_STR,           &
       NQW_END,           &
       I_QC,              &
       I_QR
    use mod_satadjust, only: &
       SATURATION_setrange, &
       SATURATION_adjustment
    use mod_cnvvar, only: &
       CNVVAR_prg2diag_in
!--- MINIAPP: only NSW6 is compiled in; the other schemes, the cloudvar/radvar
!--- containers (mod_comm), mod_rd_driver and mod_history are guarded out. ---
#ifndef MINIAPP
    use mod_mp_lsc, only: &
       mp_lsc
    use mod_mp_lsc_miroc6, only: &
       mp_lsc_miroc6
    use mod_mp_kessler, only: &
       mp_kessler
    use mod_mp_g98, only: &
       mp_g98
    use mod_mp_lin, only: &
       mp_lin
    use mod_mp_wsm6, only: &
       mp_wsm6
    use mod_mp_nsw5, only: &
       mp_nsw5
#endif
    use mod_mp_nsw6, only: &
       mp_nsw6
    use mod_mp_opt_bridge, only: &
       mp_opt_get_mode, &
       mp_opt_driver_dispatch
#ifndef MINIAPP
    use mod_mp_ndw6, only: &
       mp_ndw6
    use mod_mp_s10, only: &
       mp_s10
    use mod_cloudvar, only: &
       cloudvar_setr1_in, &
       cloudv,            &
       wait_on_cloudvar => acc_wait_on,  &
       I_RFXL,            &
       I_SFXL,            &
       I_Re_qr,           &
       I_Vt_qr,           &
       I_RWTDL,           &
       I_IWTDL,           &
       I_ISEDL,           &
       I_REVP,            &
       I_REVS,            &
       I_GDPFRC,          &
       I_FLIQ,            &
       I_ESWAGG
    use mod_radvar, only: &
       radvar_getr2_in, &
       radvar_setr1_in, &
       radv,            &
       wait_on_radvar => acc_wait_on, &
       I_RFLXSD,        &
       I_RFLXSU,        &
       I_RFLXLD,        &
       I_RFLXLU,        &
       I_RCEFF,         &
       I_RCEFF_SOLID,   &
       I_RCEFF_CLD,     &
       I_RCTOP,         &
       I_RWTOP,         &
       I_TCTOP
    use mod_rd_driver, only: &
       rd_heatingrate
    use mod_history, only: &
       history_in
#endif
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: l_region
    real(RP_PHY), intent(inout) :: RP_PHY_rhog         (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhogvx       (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhogvy       (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhogvz       (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhogw        (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhoge        (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_rhogq        (ijdim,kdim,nqmax)
    real(RP_PHY), intent(inout) :: RP_PHY_rho          (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_pre          (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_tem          (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_vx           (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_vy           (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_vz           (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_w            (ijdim,kdim)
    real(RP_PHY), intent(inout) :: RP_PHY_q            (ijdim,kdim,nqmax)
    real(RP_PHY), intent(out)   :: RP_PHY_precip       (ijdim,2)          ! surface precipitation flux (rain/snow)      [kg/m2/s]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_ISO1  (ijdim,2)          ! surface precipitation flux (rain/snow,d18O) [kg/m2/s]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_ISO2  (ijdim,2)          ! surface precipitation flux (rain/snow,dD)   [kg/m2/s]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_lh    (ijdim)            ! latent heat      of precipitation [W/m2]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_rhoein(ijdim)            ! internal  energy of precipitation [W/m2]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_rhophi(ijdim)            ! potential energy of precipitation [W/m2]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_rhokin(ijdim)            ! kinetic   energy of precipitation [W/m2]
    real(RP_PHY), intent(out)   :: RP_PHY_precip_energy(ijdim)            ! total     energy of precipitation [W/m2]
    real(RP_PHY), intent(inout) :: RP_PHY_GDCFRC       (ijdim,kdim)       ! cloud fraction [0-1]
    real(RP_PHY), intent(out)   :: RP_PHY_GDCLW        (ijdim,kdim)       ! cloud water    [kg/kg]
    real(RP_PHY), intent(in)    :: RP_PHY_CBMFX        (ijdim,kdim)       ! cloud-base mass flux [kg/m2/s]
    real(RP_PHY), intent(in)    :: RP_PHY_unccn        (ijdim,kdim)       ! aerosol number (1/m3)
    real(RP_PHY), intent(in)    :: RP_PHY_factor_n_aie (ijdim)            ! factor of aerosol indirect effet
    real(RP_PHY), intent(in)    :: RP_PHY_rhogqke      (ijdim,kdim)       ! rhog * TKE*2
    real(RP_PHY), intent(in)    :: RP_PHY_z            (ijdim,kdim)       ! altitude (full level) [m]
    real(RP_PHY), intent(in)    :: RP_PHY_zh           (ijdim,kdim)       ! altitude (half level) [m]
    real(RP_PHY), intent(in)    :: RP_PHY_ix           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_iy           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_iz           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_jx           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_jy           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_jz           (ijdim)
    real(RP_PHY), intent(in)    :: RP_PHY_phi          (ijdim,kdim)
    real(RP_PHY), intent(in)    :: RP_PHY_gam2         (ijdim,kdim)
    real(RP_PHY), intent(in)    :: RP_PHY_gam2h        (ijdim,kdim)
    real(RP_PHY), intent(in)    :: RP_PHY_gsgam2       (ijdim,kdim)
    real(RP_PHY), intent(in)    :: RP_PHY_gsgam2h      (ijdim,kdim)
    real(DP), intent(in)    :: dt
    real(DP), intent(in)    :: ct

    real(RP) :: qd        (ijdim,kdim)      ! dry air [kg/kg]
    real(RP) :: q_tot     (ijdim,kdim)

    real(RP) :: Re_liq    (ijdim,kdim)      ! [MP->RD]    Effective Radius of cloud+rain     [m]
    real(RP) :: Re_ice    (ijdim,kdim)      ! [MP->RD]    Effective Radius of ice+snow       [m]
    real(RP) :: Re_qc     (ijdim,kdim)      ! [MP->RD]    Effective Radius of cloud          [m]
    real(RP) :: rctop     (ijdim,1)         ! [MP->RD]    Effective Radius of Cloud Top      [m]
    real(RP) :: rwtop     (ijdim,1)         ! [MP->RD]    Effective Radius of Warm-Cloud Top [m]
    real(RP) :: tctop     (ijdim,1)         ! [MP->RD]    Cloud Top Temperature              [K]
    real(RP) :: rctop_cld (ijdim,1)         ! [MP->RD]    Effective Radius of Cloud Top      [m]
    real(RP) :: rwtop_cld (ijdim,1)         ! [MP->RD]    Effective Radius of Warm-Cloud To  [m]
    real(RP) :: tctop_cld (ijdim,1)         ! [MP->RD]    Cloud Top Temperature              [K]
    real(RP) :: GPREC     (ijdim,kdim,2)    ! [MP->AE,CH] precipitation flux                 [kg/m2/s]
    real(RP) :: Re_qr     (ijdim,kdim)      ! [MP->AE]    effective radius  of rain          [m]
    real(RP) :: Vt_qr0    (ijdim,kdim)      ! [MP->AE]    terminal velocity of rain          [m/s]
    real(RP) :: RWTDL     (ijdim,kdim)      ! [MP->CH]    rain tendency                      [kg/kg/s]
    real(RP) :: IWTDL     (ijdim,kdim)      ! [MP->CH]    ice  tendency                      [kg/kg/s]
    real(RP) :: ISEDL     (ijdim,kdim)      ! [MP->CH]    ice  sedimentation                 [kg/kg/s]
    real(RP) :: REVP      (ijdim,kdim)      ! [MP->CH]    rain evaporation rate              [m2/kg]
    real(RP) :: REVS      (ijdim,kdim)      ! [MP->CH]    snow evaporation rate              [m2/kg]
    real(RP) :: GDPFRC    (ijdim,kdim)      ! [MP->CH]    cloud fraction before prec.        [0-1]
    real(RP) :: FLIQ      (ijdim,kdim)      ! [MP->MP]    cloud liquid water fraction        [0-1]
    real(RP) :: ESWAGG    (ijdim,kdim)      ! [MP->MP]    riming collection efficiency

    ! for NDW6 only
    real(RP) :: RFLXSD    (ijdim,kdim,NCRF) ! shortwave downward   [W/m2]
    real(RP) :: RFLXSU    (ijdim,kdim,NCRF) ! shortwave upward     [W/m2]
    real(RP) :: RFLXLD    (ijdim,kdim,NCRF) ! longwave  downward   [W/m2]
    real(RP) :: RFLXLU    (ijdim,kdim,NCRF) ! longwave  upward     [W/m2]
    real(RP) :: frhoge_af (ijdim,kdim)
    real(RP) :: frhogqv_af(ijdim,kdim)
    real(RP) :: frhoge_rad(ijdim,kdim)      ! energy tendency by radiation
    real(RP) :: qke       (ijdim,kdim)

    ! fractional step
    real(RP) :: precip_div       (ijdim,2)
    real(RP) :: precip_ISO1_div  (ijdim,2)
    real(RP) :: precip_ISO2_div  (ijdim,2)
    real(RP) :: precip_lh_div    (ijdim)
    real(RP) :: precip_rhoein_div(ijdim)
    real(RP) :: precip_rhophi_div(ijdim)
    real(RP) :: precip_rhokin_div(ijdim)
    real(RP) :: GPREC_div        (ijdim,kdim,2)
    real(DP) :: dt_mp       ! fractional dt
    real(DP) :: fraction_mp ! 1 / MP_DIV_NUM
    real(DP) :: ct_div      ! current time in fractional step

    real(RP), parameter :: SMALL = 1.E-10_RP

    real(RP) :: rhog         (ijdim,kdim)
    real(RP) :: rhogvx       (ijdim,kdim)
    real(RP) :: rhogvy       (ijdim,kdim)
    real(RP) :: rhogvz       (ijdim,kdim)
    real(RP) :: rhogw        (ijdim,kdim)
    real(RP) :: rhoge        (ijdim,kdim)
    real(RP) :: rhogq        (ijdim,kdim,nqmax)
    real(RP) :: rho          (ijdim,kdim)
    real(RP) :: pre          (ijdim,kdim)
    real(RP) :: tem          (ijdim,kdim)
    real(RP) :: vx           (ijdim,kdim)
    real(RP) :: vy           (ijdim,kdim)
    real(RP) :: vz           (ijdim,kdim)
    real(RP) :: w            (ijdim,kdim)
    real(RP) :: q            (ijdim,kdim,nqmax)
    real(RP) :: precip       (ijdim,2)          ! surface precipitation flux (rain/snow)      [kg/m2/s]
    real(RP) :: precip_ISO1  (ijdim,2)          ! surface precipitation flux (rain/snow,d18O) [kg/m2/s]
    real(RP) :: precip_ISO2  (ijdim,2)          ! surface precipitation flux (rain/snow,dD)   [kg/m2/s]
    real(RP) :: precip_lh    (ijdim)            ! latent heat      of precipitation [W/m2]
    real(RP) :: precip_rhoein(ijdim)            ! internal  energy of precipitation [W/m2]
    real(RP) :: precip_rhophi(ijdim)            ! potential energy of precipitation [W/m2]
    real(RP) :: precip_rhokin(ijdim)            ! kinetic   energy of precipitation [W/m2]
    real(RP) :: precip_energy(ijdim)            ! total     energy of precipitation [W/m2]
    real(RP) :: GDCFRC       (ijdim,kdim)       ! cloud fraction [0-1]
    real(RP) :: GDCLW        (ijdim,kdim)       ! cloud water    [kg/kg]
    real(RP) :: CBMFX        (ijdim,kdim)       ! cloud-base mass flux [kg/m2/s]
    real(RP) :: unccn        (ijdim,kdim)       ! aerosol number (1/m3)
    real(RP) :: factor_n_aie (ijdim)            ! factor of aerosol indirect effet
    real(RP) :: rhogqke      (ijdim,kdim)       ! rhog * TKE*2
    real(RP) :: z            (ijdim,kdim)       ! altitude (full level) [m]
    real(RP) :: zh           (ijdim,kdim)       ! altitude (half level) [m]
    real(RP) :: ix           (ijdim)
    real(RP) :: iy           (ijdim)
    real(RP) :: iz           (ijdim)
    real(RP) :: jx           (ijdim)
    real(RP) :: jy           (ijdim)
    real(RP) :: jz           (ijdim)
    real(RP) :: phi          (ijdim,kdim)
    real(RP) :: gam2         (ijdim,kdim)
    real(RP) :: gam2h        (ijdim,kdim)
    real(RP) :: gsgam2       (ijdim,kdim)
    real(RP) :: gsgam2h      (ijdim,kdim)
    integer  :: ij, k, nq, m
    logical  :: opt_handled
    !---------------------------------------------------------------------------
    if ( MP_TYPE == 'NSW6' .and. mp_opt_get_mode() > 0 ) then
       call mp_opt_driver_dispatch( &
            ijdim, kdim, nqmax, dt, &
            RP_PHY_rhog, RP_PHY_rhogvx, RP_PHY_rhogvy, RP_PHY_rhogvz, RP_PHY_rhogw, RP_PHY_rhoge, RP_PHY_rhogq, &
            RP_PHY_rho, RP_PHY_pre, RP_PHY_tem, RP_PHY_vx, RP_PHY_vy, RP_PHY_vz, RP_PHY_w, RP_PHY_q, &
            RP_PHY_precip, RP_PHY_precip_ISO1, RP_PHY_precip_ISO2, &
            RP_PHY_precip_lh, RP_PHY_precip_rhoein, RP_PHY_precip_rhophi, RP_PHY_precip_rhokin, RP_PHY_precip_energy, &
            RP_PHY_GDCFRC, RP_PHY_GDCLW, &
            RP_PHY_unccn, RP_PHY_factor_n_aie, RP_PHY_z, RP_PHY_zh, &
            RP_PHY_ix, RP_PHY_iy, RP_PHY_iz, RP_PHY_jx, RP_PHY_jy, RP_PHY_jz, &
            RP_PHY_phi, RP_PHY_gam2, RP_PHY_gam2h, RP_PHY_gsgam2, RP_PHY_gsgam2h, &
            opt_handled )
       if ( opt_handled ) return
    endif

    !$acc data &
    !$acc pcreate(precip,precip_ISO1,precip_ISO2,precip_lh) &
    !$acc pcreate(precip_rhoein,precip_rhophi,precip_rhokin,precip_energy,GDCFRC,GDCLW) &
    !$acc pcreate(RP_PHY_precip,RP_PHY_precip_ISO1,RP_PHY_precip_ISO2,RP_PHY_precip_lh) &
    !$acc pcreate(RP_PHY_precip_rhoein,RP_PHY_precip_rhophi,RP_PHY_precip_rhokin,RP_PHY_precip_energy,RP_PHY_GDCLW) &
    !$acc pcreate(qd,q_tot) &
    !$acc pcreate(Re_liq,Re_ice,Re_qc,rctop,rwtop,tctop,rctop_cld,rwtop_cld,tctop_cld) &
    !$acc pcreate(GPREC,Re_qr,Vt_qr0,RWTDL,IWTDL,ISEDL,REVP,REVS,GDPFRC,FLIQ,ESWAGG) &
    !$acc pcreate(precip_div,precip_ISO1_div,precip_ISO2_div,precip_lh_div) &
    !$acc pcreate(precip_rhoein_div,precip_rhophi_div,precip_rhokin_div,GPREC_div) &
    !$acc pcreate(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq,rho,pre,tem,vx,vy,vz,w,q ) &
    !$acc pcreate(CBMFX,unccn,factor_n_aie,rhogqke,z,zh,ix,iy,iz,jx,jy,jz,phi,gam2,gam2h,gsgam2,gsgam2h) &
    !$acc pcopyin(RP_PHY_CBMFX,RP_PHY_unccn,RP_PHY_factor_n_aie,RP_PHY_rhogqke,RP_PHY_z,RP_PHY_zh,RP_PHY_ix,RP_PHY_iy,RP_PHY_iz,RP_PHY_jx,RP_PHY_jy,RP_PHY_jz,RP_PHY_phi,RP_PHY_gam2,RP_PHY_gam2h,RP_PHY_gsgam2,RP_PHY_gsgam2h) &
!--- MINIAPP: cloudv / radv belong to the removed diagnostic containers. ---
#ifndef MINIAPP
    !$acc pcopyin(RP_PHY_GDCFRC) &
    !$acc pcopy(cloudv,radv)
#else
    !$acc pcopyin(RP_PHY_GDCFRC)
#endif

    ! single/double precision conversion adapter
    !$acc kernels async(1)
    do k = 1, kdim
    do ij = 1, ijdim
       rhog   (ij,k) = RP_PHY_rhog   (ij,k)
       rhogvx (ij,k) = RP_PHY_rhogvx (ij,k)
       rhogvy (ij,k) = RP_PHY_rhogvy (ij,k)
       rhogvz (ij,k) = RP_PHY_rhogvz (ij,k)
       rhogw  (ij,k) = RP_PHY_rhogw  (ij,k)
       rhoge  (ij,k) = RP_PHY_rhoge  (ij,k)
       rho    (ij,k) = RP_PHY_rho    (ij,k)
       pre    (ij,k) = RP_PHY_pre    (ij,k)
       tem    (ij,k) = RP_PHY_tem    (ij,k)
       vx     (ij,k) = RP_PHY_vx     (ij,k)
       vy     (ij,k) = RP_PHY_vy     (ij,k)
       vz     (ij,k) = RP_PHY_vz     (ij,k)
       w      (ij,k) = RP_PHY_w      (ij,k)
       CBMFX  (ij,k) = RP_PHY_CBMFX  (ij,k)
       unccn  (ij,k) = RP_PHY_unccn  (ij,k)
       rhogqke(ij,k) = RP_PHY_rhogqke(ij,k)
       z      (ij,k) = RP_PHY_z      (ij,k)
       zh     (ij,k) = RP_PHY_zh     (ij,k)
       phi    (ij,k) = RP_PHY_phi    (ij,k)
       gam2   (ij,k) = RP_PHY_gam2   (ij,k)
       gam2h  (ij,k) = RP_PHY_gam2h  (ij,k)
       gsgam2 (ij,k) = RP_PHY_gsgam2 (ij,k)
       gsgam2h(ij,k) = RP_PHY_gsgam2h(ij,k)
       GDCFRC (ij,k) = RP_PHY_GDCFRC(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1)
    do ij = 1, ijdim
      ix (ij) = RP_PHY_ix (ij)
      iy (ij) = RP_PHY_iy (ij)
      iz (ij) = RP_PHY_iz (ij)
      jx (ij) = RP_PHY_jx (ij)
      jy (ij) = RP_PHY_jy (ij)
      jz (ij) = RP_PHY_jz (ij)
      factor_n_aie(ij) = RP_PHY_factor_n_aie(ij)
    enddo
    !$acc end kernels

    !$acc kernels async(1)
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
      rhogq (ij,k,nq)  = RP_PHY_rhogq (ij,k,nq)
      q     (ij,k,nq)  = RP_PHY_q     (ij,k,nq)
    enddo
    enddo
    enddo
    !$acc end kernels

    fraction_mp = 1.0_RP / real(MP_DIV_NUM,kind=RP)
    dt_mp       = dt * fraction_mp

    !$omp parallel default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,precip,precip_ISO1,precip_ISO2,           &
    !$omp        precip_lh,precip_rhoein,precip_rhophi,precip_rhokin, &
    !$omp        GPREC,RWTDL,IWTDL,ISEDL,REVS,REVP,GDPFRC             )
!OCL XFILL
    !$acc kernels pcopy(precip,precip_ISO1,precip_ISO2) &
    !$acc         pcopy(precip_lh,precip_rhoein,precip_rhophi,precip_rhokin) async(1)
    !$omp do
    do ij = 1, ijdim
       precip       (ij,1) = 0.0_RP
       precip       (ij,2) = 0.0_RP
       precip_ISO1  (ij,1) = 0.0_RP
       precip_ISO1  (ij,2) = 0.0_RP
       precip_ISO2  (ij,1) = 0.0_RP
       precip_ISO2  (ij,2) = 0.0_RP
       precip_lh    (ij)   = 0.0_RP
       precip_rhoein(ij)   = 0.0_RP
       precip_rhophi(ij)   = 0.0_RP
       precip_rhokin(ij)   = 0.0_RP
    enddo
    !$omp end do nowait
    !$acc end kernels
!OCL XFILL
    !$acc kernels pcopy(GPREC,RWTDL,IWTDL,ISEDL,REVS,REVP,GDPFRC) async(1)
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       GPREC (ij,k,1) = 0.0_RP
       GPREC (ij,k,2) = 0.0_RP
       RWTDL (ij,k)   = 0.0_RP
       IWTDL (ij,k)   = 0.0_RP
       ISEDL (ij,k)   = 0.0_RP
       REVS  (ij,k)   = 0.0_RP
       REVP  (ij,k)   = 0.0_RP
       GDPFRC(ij,k)   = 0.0_RP
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    !$omp end parallel

    do m = 1, MP_DIV_NUM
       !$acc wait(1)
       !$omp parallel default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,precip_div,precip_ISO1_div,precip_ISO2_div,GPREC_div,   &
       !$omp        precip_lh_div,precip_rhoein_div,precip_rhophi_div,precip_rhokin_div)
!OCL XFILL
       !$acc kernels pcopy(precip_div,precip_ISO1_div,precip_ISO2_div) &
       !$acc         pcopy(precip_lh_div,precip_rhoein_div,precip_rhophi_div,precip_rhokin_div)
       !$omp do
       do ij = 1, ijdim
          precip_div       (ij,1) = 0.0_RP
          precip_div       (ij,2) = 0.0_RP
          precip_ISO1_div  (ij,1) = 0.0_RP
          precip_ISO1_div  (ij,2) = 0.0_RP
          precip_ISO2_div  (ij,1) = 0.0_RP
          precip_ISO2_div  (ij,2) = 0.0_RP
          precip_lh_div    (ij)   = 0.0_RP
          precip_rhoein_div(ij)   = 0.0_RP
          precip_rhophi_div(ij)   = 0.0_RP
          precip_rhokin_div(ij)   = 0.0_RP
       enddo
       !$omp end do nowait
       !$acc end kernels
!OCL XFILL
       !$acc kernels pcopy(GPREC_div)
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          GPREC_div(ij,k,1) = 0.0_RP
          GPREC_div(ij,k,2) = 0.0_RP
       enddo
       enddo
       !$omp end do
       !$acc end kernels

       !$omp end parallel

       if ( MP_TYPE == 'NONE' ) then

          if ( m == MP_DIV_NUM ) then
             !$acc kernels async(1)
             !$omp parallel workshare
             GDCFRC   (:,:) = 0.0_RP
             Re_liq   (:,:) = 10.E-6_RP ! dummy
             Re_ice   (:,:) = 20.E-6_RP ! dummy
             Re_qc    (:,:) = 10.E-6_RP ! dummy
             Re_qr    (:,:) = CNST_UNDEF
             Vt_qr0   (:,:) = CNST_UNDEF
             rctop    (:,:) = 10.E-6_RP ! dummy
             rwtop    (:,:) = 10.E-6_RP ! dummy
             tctop    (:,:) = CNST_UNDEF
             rctop_cld(:,:) = 10.E-6_RP ! dummy
             rwtop_cld(:,:) = 10.E-6_RP ! dummy
             tctop_cld(:,:) = CNST_UNDEF
             !$omp end parallel workshare
             !$acc end kernels
          endif

       elseif( MP_TYPE == 'LSC' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_lsc( ijdim,                    & ! [IN]
                       rhog             (:,:),   & ! [INOUT]
                       rhogvx           (:,:),   & ! [INOUT]
                       rhogvy           (:,:),   & ! [INOUT]
                       rhogvz           (:,:),   & ! [INOUT]
                       rhogw            (:,:),   & ! [INOUT]
                       rhoge            (:,:),   & ! [INOUT]
                       rhogq            (:,:,:), & ! [INOUT]
                       precip_div       (:,:),   & ! [OUT]
                       precip_ISO1_div  (:,:),   & ! [OUT]
                       precip_ISO2_div  (:,:),   & ! [OUT]
                       precip_rhoein_div(:),     & ! [OUT]
                       GDCFRC           (:,:),   & ! [OUT]
                       GDCLW            (:,:),   & ! [OUT]
                       GPREC_div        (:,:,:), & ! [OUT]
                       Re_liq           (:,:),   & ! [OUT]
                       Re_ice           (:,:),   & ! [OUT]
                       rctop            (:,:),   & ! [OUT]
                       rwtop            (:,:),   & ! [OUT]
                       tctop            (:,:),   & ! [OUT]
                       RWTDL            (:,:),   & ! [OUT]
                       IWTDL            (:,:),   & ! [OUT]
                       ISEDL            (:,:),   & ! [OUT]
                       REVP             (:,:),   & ! [OUT]
                       REVS             (:,:),   & ! [OUT]
                       GDPFRC           (:,:),   & ! [OUT]
                       CBMFX            (:,:),   & ! [IN]
                       rho              (:,:),   & ! [IN]
                       tem              (:,:),   & ! [IN]
                       pre              (:,:),   & ! [IN]
                       gsgam2           (:,:),   & ! [IN]
                       gsgam2h          (:,:),   & ! [IN]
                       q                (:,:,:), & ! [IN]
                       unccn            (:,:),   & ! [IN]
                       z                (:,:),   & ! [IN]
                       zh               (:,:),   & ! [IN]
                       dt_mp                     ) ! [IN]
#endif

          if ( m == MP_DIV_NUM ) then
             !$acc kernels async(1)
             !$omp parallel workshare
             Re_qc    (:,:) = Re_liq(:,:) ! only qc is prognostic
             Re_qr    (:,:) = CNST_UNDEF
             Vt_qr0   (:,:) = CNST_UNDEF
             rctop_cld(:,:) = rctop(:,:)
             rwtop_cld(:,:) = rwtop(:,:)
             tctop_cld(:,:) = tctop(:,:)
             !$omp end parallel workshare
             !$acc end kernels
          endif
       elseif( MP_TYPE == 'LSC_MIROC6' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_lsc_miroc6( ijdim,                    & ! [IN]
                              l_region,                 & ! [IN]
                              rhog             (:,:),   & ! [INOUT]
                              rhogvx           (:,:),   & ! [INOUT]
                              rhogvy           (:,:),   & ! [INOUT]
                              rhogvz           (:,:),   & ! [INOUT]
                              rhogw            (:,:),   & ! [INOUT]
                              rhoge            (:,:),   & ! [INOUT]
                              rhogq            (:,:,:), & ! [INOUT]
                              precip_div       (:,:),   & ! [OUT]
                              precip_rhoein_div(:),     & ! [OUT]
                              GDCFRC           (:,:),   & ! [INOUT]
                              GDCLW            (:,:),   & ! [OUT]
                              GPREC_div        (:,:,:), & ! [OUT]
!                              Re_liq           (:,:),   & ! [OUT]
!                              Re_ice           (:,:),   & ! [OUT]
!                              rctop            (:,:),   & ! [OUT]
!                              rwtop            (:,:),   & ! [OUT]
!                              tctop            (:,:),   & ! [OUT]
                              RWTDL            (:,:),   & ! [OUT]
                              IWTDL            (:,:),   & ! [OUT]
                              ISEDL            (:,:),   & ! [OUT]
                              REVP             (:,:),   & ! [OUT]
                              REVS             (:,:),   & ! [OUT]
                              GDPFRC           (:,:),   & ! [OUT]
                              FLIQ             (:,:),   & ! [OUT]
                              ESWAGG           (:,:),   & ! [OUT]
                              rho              (:,:),   & ! [IN]
                              tem              (:,:),   & ! [IN]
                              pre              (:,:),   & ! [IN]
                              gsgam2           (:,:),   & ! [IN]
                              gsgam2h          (:,:),   & ! [IN]
                              q                (:,:,:), & ! [IN]
                              unccn            (:,:),   & ! [IN]
                              z                (:,:),   & ! [IN]
                              zh               (:,:),   & ! [IN]
                              dt_mp                     ) ! [IN]
#endif

          if ( m == MP_DIV_NUM ) then
             !$acc kernels async(1)
             !$omp parallel workshare
             Re_liq   (:,:) = 10.E-6_RP ! dummy
             Re_ice   (:,:) = 20.E-6_RP ! dummy
             Re_qc    (:,:) = 10.E-6_RP ! dummy
             Re_qr    (:,:) = CNST_UNDEF
             Vt_qr0   (:,:) = CNST_UNDEF
             rctop_cld(:,:) = 10.E-6_RP ! dummy
             rwtop_cld(:,:) = 10.E-6_RP ! dummy
             tctop_cld(:,:) = CNST_UNDEF
             !$omp end parallel workshare
             !$acc end kernels
          endif
       elseif( MP_TYPE == 'KESSLER' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_kessler( ijdim,             & ! [IN]
                           rhog,              & ! [INOUT]
                           rhogvx,            & ! [INOUT]
                           rhogvy,            & ! [INOUT]
                           rhogvz,            & ! [INOUT]
                           rhogw,             & ! [INOUT]
                           rhoge,             & ! [INOUT]
                           rhogq,             & ! [INOUT]
                           vx,                & ! [IN]
                           vy,                & ! [IN]
                           vz,                & ! [IN]
                           w,                 & ! [IN]
                           rho,               & ! [INOUT]
                           tem,               & ! [INOUT]
                           pre,               & ! [INOUT]
                           q,                 & ! [INOUT]
                           qd,                & ! [OUT]
                           precip(:,1),       & ! [OUT]
                           precip_rhoein_div, & ! [OUT]
                           precip_lh_div,     & ! [OUT]
                           precip_rhophi_div, & ! [OUT]
                           precip_rhokin_div, & ! [OUT]
                           GPREC_div,         & ! [OUT]
                           gsgam2,            & ! [IN]
                           gsgam2h,           & ! [IN]
                           gam2,              & ! [IN]
                           gam2h,             & ! [IN]
                           ix,                & ! [IN]
                           iy,                & ! [IN]
                           iz,                & ! [IN]
                           jx,                & ! [IN]
                           jy,                & ! [IN]
                           jz,                & ! [IN]
                           z,                 & ! [IN]
                           zh,                & ! [IN]
                           dt_mp              ) ! [IN]
#endif

          call SATURATION_setrange( 10.0_RP, 5.0_RP ) ! [IN]

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

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_liq   (:,:) = 10.E-6_RP ! dummy
             Re_ice   (:,:) = 20.E-6_RP ! dummy
             Re_qc    (:,:) = 10.E-6_RP ! dummy
             Re_qr    (:,:) = CNST_UNDEF
             Vt_qr0   (:,:) = CNST_UNDEF
             rctop    (:,:) = 10.E-6_RP ! dummy
             rwtop    (:,:) = 10.E-6_RP ! dummy
             rctop_cld(:,:) = 10.E-6_RP ! dummy
             rwtop_cld(:,:) = 10.E-6_RP ! dummy
             !$omp end parallel workshare

             do k  = 1, kdim
             do ij = 1, ijdim
                if ( rho(ij,k)*(q(ij,k,I_QC)+q(ij,k,I_QR)) > 1.E-5_RP ) then
                   tctop    (ij,1) = tem(ij,k)
                   tctop_cld(ij,1) = tem(ij,k)
                endif
             enddo
             enddo
          endif

       elseif ( MP_TYPE == 'G98' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_g98( ijdim,             & ! [IN]
                       rhog,              & ! [INOUT]
                       rhogvx,            & ! [INOUT]
                       rhogvy,            & ! [INOUT]
                       rhogvz,            & ! [INOUT]
                       rhogw,             & ! [INOUT]
                       rhoge,             & ! [INOUT]
                       rhogq,             & ! [INOUT]
                       vx,                & ! [IN]
                       vy,                & ! [IN]
                       vz,                & ! [IN]
                       w,                 & ! [IN]
                       unccn,             & ! [IN]
                       rho,               & ! [INOUT]
                       tem,               & ! [INOUT]
                       pre,               & ! [INOUT]
                       q,                 & ! [INOUT]
                       qd,                & ! [OUT]
                       precip(:,1),       & ! [OUT]
                       precip_rhoein_div, & ! [OUT]
                       precip_lh_div,     & ! [OUT]
                       precip_rhophi_div, & ! [OUT]
                       precip_rhokin_div, & ! [OUT]
                       GPREC_div,         & ! [OUT]
                       Re_liq,            & ! [OUT]
                       rctop,             & ! [OUT]
                       rwtop,             & ! [OUT]
                       tctop,             & ! [OUT]
                       Re_qc,             & ! [OUT]
                       rctop_cld,         & ! [OUT]
                       rwtop_cld,         & ! [OUT]
                       tctop_cld,         & ! [OUT]
                       gsgam2,            & ! [IN]
                       gsgam2h,           & ! [IN]
                       gam2,              & ! [IN]
                       gam2h,             & ! [IN]
                       ix,                & ! [IN]
                       iy,                & ! [IN]
                       iz,                & ! [IN]
                       jx,                & ! [IN]
                       jy,                & ! [IN]
                       jz,                & ! [IN]
                       z,                 & ! [IN]
                       zh,                & ! [IN]
                       dt_mp              ) ! [IN]
#endif

          call SATURATION_setrange( 273.15_RP, 258.15_RP ) ! [IN]

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

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_ice(:,:) = 20.E-6_RP ! dummy
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
          endif

       elseif ( MP_TYPE == 'LIN' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_lin( ijdim,             & ! [IN]
                       rhog,              & ! [INOUT]
                       rhogvx,            & ! [INOUT]
                       rhogvy,            & ! [INOUT]
                       rhogvz,            & ! [INOUT]
                       rhogw,             & ! [INOUT]
                       rhoge,             & ! [INOUT]
                       rhogq,             & ! [INOUT]
                       vx,                & ! [IN]
                       vy,                & ! [IN]
                       vz,                & ! [IN]
                       w,                 & ! [IN]
                       unccn,             & ! [IN]
                       rho,               & ! [INOUT]
                       tem,               & ! [INOUT]
                       pre,               & ! [INOUT]
                       q,                 & ! [INOUT]
                       qd,                & ! [OUT]
                       precip_div,        & ! [OUT]
                       precip_rhoein_div, & ! [OUT]
                       precip_lh_div,     & ! [OUT]
                       precip_rhophi_div, & ! [OUT]
                       precip_rhokin_div, & ! [OUT]
                       GPREC_div,         & ! [OUT]
                       Re_liq,            & ! [OUT]
                       rctop,             & ! [OUT]
                       rwtop,             & ! [OUT]
                       tctop,             & ! [OUT]
                       Re_qc,             & ! [OUT]
                       rctop_cld,         & ! [OUT]
                       rwtop_cld,         & ! [OUT]
                       tctop_cld,         & ! [OUT]
                       gsgam2,            & ! [IN]
                       gsgam2h,           & ! [IN]
                       gam2,              & ! [IN]
                       gam2h,             & ! [IN]
                       ix,                & ! [IN]
                       iy,                & ! [IN]
                       iz,                & ! [IN]
                       jx,                & ! [IN]
                       jy,                & ! [IN]
                       jz,                & ! [IN]
                       z,                 & ! [IN]
                       zh,                & ! [IN]
                       dt_mp              ) ! [IN]
#endif

          ! saturation adjustment is already done in Lin scheme.

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_ice(:,:) = 20.E-6_RP ! dummy
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
          endif

       elseif ( MP_TYPE == 'WSM6' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_wsm6( ijdim,             & ! [IN]
                        rhog,              & ! [INOUT]
                        rhogvx,            & ! [INOUT]
                        rhogvy,            & ! [INOUT]
                        rhogvz,            & ! [INOUT]
                        rhogw,             & ! [INOUT]
                        rhoge,             & ! [INOUT]
                        rhogq,             & ! [INOUT]
                        vx,                & ! [IN]
                        vy,                & ! [IN]
                        vz,                & ! [IN]
                        w,                 & ! [IN]
                        unccn,             & ! [IN]
                        rho,               & ! [INOUT]
                        tem,               & ! [INOUT]
                        pre,               & ! [INOUT]
                        q,                 & ! [INOUT]
                        qd,                & ! [OUT]
                        precip_div,        & ! [OUT]
                        precip_rhoein_div, & ! [OUT]
                        precip_lh_div,     & ! [OUT]
                        precip_rhophi_div, & ! [OUT]
                        precip_rhokin_div, & ! [OUT]
                        GPREC_div,         & ! [OUT]
                        Re_liq,            & ! [OUT]
                        rctop,             & ! [OUT]
                        rwtop,             & ! [OUT]
                        tctop,             & ! [OUT]
                        Re_qc,             & ! [OUT]
                        rctop_cld,         & ! [OUT]
                        rwtop_cld,         & ! [OUT]
                        tctop_cld,         & ! [OUT]
                        gsgam2,            & ! [IN]
                        gsgam2h,           & ! [IN]
                        gam2,              & ! [IN]
                        gam2h,             & ! [IN]
                        ix,                & ! [IN]
                        iy,                & ! [IN]
                        iz,                & ! [IN]
                        jx,                & ! [IN]
                        jy,                & ! [IN]
                        jz,                & ! [IN]
                        z,                 & ! [IN]
                        zh,                & ! [IN]
                        dt_mp              ) ! [IN]
#endif

          call SATURATION_setrange( 100.0_RP, 90.0_RP ) ! [IN]

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

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_ice(:,:) = 20.E-6_RP ! dummy
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
          endif

       elseif ( MP_TYPE == 'NSW5' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_nsw5( ijdim,             & ! [IN]
                        rhog,              & ! [INOUT]
                        rhogvx,            & ! [INOUT]
                        rhogvy,            & ! [INOUT]
                        rhogvz,            & ! [INOUT]
                        rhogw,             & ! [INOUT]
                        rhoge,             & ! [INOUT]
                        rhogq,             & ! [INOUT]
                        vx,                & ! [IN]
                        vy,                & ! [IN]
                        vz,                & ! [IN]
                        w,                 & ! [IN]
                        unccn,             & ! [IN]
                        rho,               & ! [INOUT]
                        tem,               & ! [INOUT]
                        pre,               & ! [INOUT]
                        q,                 & ! [INOUT]
                        qd,                & ! [OUT]
                        precip_div,        & ! [OUT]
                        precip_rhoein_div, & ! [OUT]
                        precip_lh_div,     & ! [OUT]
                        precip_rhophi_div, & ! [OUT]
                        precip_rhokin_div, & ! [OUT]
                        GPREC_div,         & ! [OUT]
                        Re_liq,            & ! [OUT]
                        rctop,             & ! [OUT]
                        rwtop,             & ! [OUT]
                        tctop,             & ! [OUT]
                        Re_qc,             & ! [OUT]
                        rctop_cld,         & ! [OUT]
                        rwtop_cld,         & ! [OUT]
                        tctop_cld,         & ! [OUT]
                        gsgam2,            & ! [IN]
                        gsgam2h,           & ! [IN]
                        gam2,              & ! [IN]
                        gam2h,             & ! [IN]
                        ix,                & ! [IN]
                        iy,                & ! [IN]
                        iz,                & ! [IN]
                        jx,                & ! [IN]
                        jy,                & ! [IN]
                        jz,                & ! [IN]
                        z,                 & ! [IN]
                        zh,                & ! [IN]
                        dt_mp              ) ! [IN]
#endif

          call SATURATION_setrange( 273.16_RP, 233.16_RP ) ! [IN]

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

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_ice(:,:) = 20.E-6_RP ! dummy
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
          endif

       elseif ( MP_TYPE == 'NSW6' ) then

          call mp_nsw6( ijdim,             & ! [IN]
                        l_region,          & ! [IN]
                        rhog,              & ! [INOUT]
                        rhogvx,            & ! [INOUT]
                        rhogvy,            & ! [INOUT]
                        rhogvz,            & ! [INOUT]
                        rhogw,             & ! [INOUT]
                        rhoge,             & ! [INOUT]
                        rhogq,             & ! [INOUT]
                        vx,                & ! [IN]
                        vy,                & ! [IN]
                        vz,                & ! [IN]
                        w,                 & ! [IN]
                        unccn,             & ! [IN]
                        factor_n_aie,      & ! [IN]
                        rho,               & ! [INOUT]
                        tem,               & ! [INOUT]
                        pre,               & ! [INOUT]
                        q,                 & ! [INOUT]
                        qd,                & ! [OUT]
                        precip_div,        & ! [OUT]
                        precip_ISO1_div,   & ! [OUT]
                        precip_ISO2_div,   & ! [OUT]
                        precip_rhoein_div, & ! [OUT]
                        precip_lh_div,     & ! [OUT]
                        precip_rhophi_div, & ! [OUT]
                        precip_rhokin_div, & ! [OUT]
                        GPREC_div,         & ! [OUT]
                        REVP,              & ! [OUT]
                        REVS,              & ! [OUT]
                        Re_liq,            & ! [OUT]
                        rctop,             & ! [OUT]
                        rwtop,             & ! [OUT]
                        tctop,             & ! [OUT]
                        Re_qc,             & ! [OUT]
                        rctop_cld,         & ! [OUT]
                        rwtop_cld,         & ! [OUT]
                        tctop_cld,         & ! [OUT]
                        gsgam2,            & ! [IN]
                        gsgam2h,           & ! [IN]
                        gam2,              & ! [IN]
                        gam2h,             & ! [IN]
                        ix,                & ! [IN]
                        iy,                & ! [IN]
                        iz,                & ! [IN]
                        jx,                & ! [IN]
                        jy,                & ! [IN]
                        jz,                & ! [IN]
                        z,                 & ! [IN]
                        zh,                & ! [IN]
                        dt_mp              ) ! [IN]

          if ( m == MP_DIV_NUM ) then
             !$acc kernels async(1)
             !$omp parallel workshare
             Re_ice(:,:) = 20.E-6_RP ! dummy
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
             !$acc end kernels
          endif

       elseif( MP_TYPE == 'NDW6' ) then

!--- MINIAPP: NDW6 branch (radvar + rd_driver) removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          ct_div = ct + dt_mp * real(m,kind=DP)

          if ( m == 1 ) then
             call radvar_getr2_in( RFLXSD(:,:,:), vid=I_RFLXSD, klev1=kdim, klev2=NCRF, l_region=l_region ) ! [OUT]
             call radvar_getr2_in( RFLXSU(:,:,:), vid=I_RFLXSU, klev1=kdim, klev2=NCRF, l_region=l_region ) ! [OUT]
             call radvar_getr2_in( RFLXLD(:,:,:), vid=I_RFLXLD, klev1=kdim, klev2=NCRF, l_region=l_region ) ! [OUT]
             call radvar_getr2_in( RFLXLU(:,:,:), vid=I_RFLXLU, klev1=kdim, klev2=NCRF, l_region=l_region ) ! [OUT]

             ! For unresolved supersaturation ( Ghan et al.,1997, GRL, Morrison etal.,2005, JAS)
             call rd_heatingrate( ijdim,       & ! [IN] : number of horizontal grid
                                  RFLXSD    (:,:,:), & ! [IN] : downward short wave(half lev)
                                  RFLXSU    (:,:,:), & ! [IN] : upward shot wave(half lev)
                                  RFLXLD    (:,:,:), & ! [IN] : upward long wave(half lev)
                                  RFLXLU    (:,:,:), & ! [IN] : upward long wave(half lev)
                                  frhoge_rad(:,:),   & ! [OUT]: energy tendency
                                  gam2h     (:,:)    ) ! [IN] : metric coefficient

             !$omp parallel workshare
             frhoge_af (:,:) = 0.0_RP
             frhogqv_af(:,:) = 0.0_RP
             qke       (:,:) = rhogqke(:,:) / rhog(:,:)
             !$omp end parallel workshare
          endif

          call mp_ndw6( ijdim,             & ! [IN]
                        l_region,          & ! [IN]
                        m, MP_DIV_NUM,     & ! [IN]
                        rhog,              & ! [INOUT]
                        rhogvx,            & ! [INOUT]
                        rhogvy,            & ! [INOUT]
                        rhogvz,            & ! [INOUT]
                        rhogw,             & ! [INOUT]
                        rhoge,             & ! [INOUT]
                        rhogq,             & ! [INOUT]
                        vx,                & ! [IN]
                        vy,                & ! [IN]
                        vz,                & ! [IN]
                        w,                 & ! [IN]
                        unccn,             & ! [IN]
                        rho,               & ! [INOUT]
                        tem,               & ! [INOUT]
                        pre,               & ! [INOUT]
                        q,                 & ! [INOUT]
                        qd,                & ! [OUT]
                        precip_div,        & ! [OUT]
                        precip_rhoein_div, & ! [OUT]
                        precip_lh_div,     & ! [OUT]
                        precip_rhophi_div, & ! [OUT]
                        precip_rhokin_div, & ! [OUT]
                        GPREC_div,         & ! [OUT]
                        Re_liq,            & ! [OUT]
                        Re_ice,            & ! [OUT]
                        rctop,             & ! [OUT]
                        rwtop,             & ! [OUT]
                        tctop,             & ! [OUT]
                        Re_qc,             & ! [OUT]
                        Re_qr,             & ! [OUT]
                        rctop_cld,         & ! [OUT]
                        rwtop_cld,         & ! [OUT]
                        tctop_cld,         & ! [OUT]
                        Vt_qr0,            & ! [OUT]
                        frhoge_af,         & ! [IN] energy tendency by additional forcing
                        frhogqv_af,        & ! [IN] vapor  tendency by additional forcing
                        frhoge_rad,        & ! [IN] energy tendency by radiation
                        qke,               & ! [IN] 2*TKE
                        gsgam2,            & ! [IN]
                        gsgam2h,           & ! [IN]
                        gam2,              & ! [IN]
                        gam2h,             & ! [IN]
                        ix,                & ! [IN]
                        iy,                & ! [IN]
                        iz,                & ! [IN]
                        jx,                & ! [IN]
                        jy,                & ! [IN]
                        jz,                & ! [IN]
                        z,                 & ! [IN]
                        zh,                & ! [IN]
                        dt_mp,             & ! [IN]
                        ct_div             ) ! [IN] 09/04/14 [Add] T.Mitsui
#endif

       elseif ( MP_TYPE == 'S10' ) then

!--- MINIAPP: non-NSW6 scheme call removed with -DMINIAPP. ---
#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP) && !defined(MINIAPP)
          call mp_s10( ijdim,             & ! [IN]
                       l_region,          & ! [IN]
                       rhog,              & ! [INOUT]
                       rhogvx,            & ! [INOUT]
                       rhogvy,            & ! [INOUT]
                       rhogvz,            & ! [INOUT]
                       rhogw,             & ! [INOUT]
                       rhoge,             & ! [INOUT]
                       rhogq,             & ! [INOUT]
                       vx,                & ! [IN]
                       vy,                & ! [IN]
                       vz,                & ! [IN]
                       w,                 & ! [IN]
                       unccn,             & ! [IN]
                       rho,               & ! [INOUT]
                       tem,               & ! [INOUT]
                       pre,               & ! [INOUT]
                       q,                 & ! [INOUT]
                       qd,                & ! [OUT]
                       precip_div,        & ! [OUT]
                       precip_rhoein_div, & ! [OUT]
                       precip_lh_div,     & ! [OUT]
                       precip_rhophi_div, & ! [OUT]
                       precip_rhokin_div, & ! [OUT]
                       GPREC_div,         & ! [OUT]
                       Re_liq,            & ! [OUT]
                       Re_ice,            & ! [OUT]
                       rctop,             & ! [OUT]
                       rwtop,             & ! [OUT]
                       tctop,             & ! [OUT]
                       Re_qc,             & ! [OUT]
                       rctop_cld,         & ! [OUT]
                       rwtop_cld,         & ! [OUT]
                       tctop_cld,         & ! [OUT]
                       gsgam2,            & ! [IN]
                       gsgam2h,           & ! [IN]
                       gam2,              & ! [IN]
                       gam2h,             & ! [IN]
                       ix,                & ! [IN]
                       iy,                & ! [IN]
                       iz,                & ! [IN]
                       jx,                & ! [IN]
                       jy,                & ! [IN]
                       jz,                & ! [IN]
                       z,                 & ! [IN]
                       zh,                & ! [IN]
                       dt_mp              ) ! [IN]
#endif

          if ( m == MP_DIV_NUM ) then
             !$omp parallel workshare
             Re_qr (:,:) = CNST_UNDEF
             Vt_qr0(:,:) = CNST_UNDEF
             !$omp end parallel workshare
          endif

       endif

       !$omp parallel default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,fraction_mp,                                              &
       !$omp        precip_div,precip_ISO1_div,precip_ISO2_div,GPREC_div,                &
       !$omp        precip_lh_div,precip_rhoein_div,precip_rhophi_div,precip_rhokin_div, &
       !$omp        precip,precip_ISO1,precip_ISO2,GPREC,                                &
       !$omp        precip_lh,precip_rhoein,precip_rhophi,precip_rhokin                  )

       !$acc kernels pcopy(precip_div,precip_ISO1_div,precip_ISO2_div) &
       !$acc         pcopy(precip_lh_div,precip_rhoein_div,precip_rhophi_div,precip_rhokin_div) &
       !$acc         pcopy(precip,precip_ISO1,precip_ISO2) &
       !$acc         pcopy(precip_lh,precip_rhoein,precip_rhophi,precip_rhokin) async(1)
       !$omp do
       do ij = 1, ijdim
          precip       (ij,1) = precip       (ij,1) + precip_div       (ij,1) * fraction_mp
          precip       (ij,2) = precip       (ij,2) + precip_div       (ij,2) * fraction_mp
          precip_ISO1  (ij,1) = precip_ISO1  (ij,1) + precip_ISO1_div  (ij,1) * fraction_mp
          precip_ISO1  (ij,2) = precip_ISO1  (ij,2) + precip_ISO1_div  (ij,2) * fraction_mp
          precip_ISO2  (ij,1) = precip_ISO2  (ij,1) + precip_ISO2_div  (ij,1) * fraction_mp
          precip_ISO2  (ij,2) = precip_ISO2  (ij,2) + precip_ISO2_div  (ij,2) * fraction_mp
          precip_lh    (ij)   = precip_lh    (ij)   + precip_lh_div    (ij)   * fraction_mp
          precip_rhoein(ij)   = precip_rhoein(ij)   + precip_rhoein_div(ij)   * fraction_mp
          precip_rhophi(ij)   = precip_rhophi(ij)   + precip_rhophi_div(ij)   * fraction_mp
          precip_rhokin(ij)   = precip_rhokin(ij)   + precip_rhokin_div(ij)   * fraction_mp
       enddo
       !$omp end do nowait
       !$acc end kernels

       !$acc kernels pcopy(GPREC_div,GPREC) async(1)
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          GPREC(ij,k,1) = GPREC(ij,k,1) + GPREC_div(ij,k,1) * fraction_mp
          GPREC(ij,k,2) = GPREC(ij,k,2) + GPREC_div(ij,k,2) * fraction_mp
       enddo
       enddo
       !$omp end do
       !$acc end kernels

       !$omp end parallel
    enddo



    if ( MP_TYPE /= 'NONE' ) then
!OCL XFILL
       !$acc kernels pcopy(precip_energy) pcopyin(precip_lh,precip_rhoein,precip_rhophi,precip_rhokin) async(1)
       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,precip_energy,precip_lh,precip_rhoein,precip_rhophi,precip_rhokin)
       do ij = 1, ijdim
          precip_energy(ij) = precip_lh    (ij) &
                            + precip_rhoein(ij) &
                            + precip_rhophi(ij) &
                            + precip_rhokin(ij)
       enddo
       !$omp end parallel do
       !$acc end kernels

       !$acc wait(1)
       call CNVVAR_prg2diag_in( ijdim,              & ! [IN]
                                rho    (:,:),       & ! [OUT]
                                pre    (:,:),       & ! [OUT]
                                tem    (:,:),       & ! [OUT]
                                vx     (:,:),       & ! [OUT]
                                vy     (:,:),       & ! [OUT]
                                vz     (:,:),       & ! [OUT]
                                w      (:,:),       & ! [OUT]
                                q      (:,:,:),     & ! [OUT]
                                rhog   (:,:),       & ! [IN]
                                rhogvx (:,:),       & ! [IN]
                                rhogvy (:,:),       & ! [IN]
                                rhogvz (:,:),       & ! [IN]
                                rhogw  (:,:),       & ! [IN]
                                rhoge  (:,:),       & ! [IN]
                                rhogq  (:,:,:),     & ! [IN]
                                phi    (:,:),       & ! [IN]
                                gsgam2 (:,:),       & ! [IN]
                                gsgam2h(:,:),       & ! [IN]
                                set_boundary=.true. ) ! [IN]

       ! Diagnose cloud fraction
       if ( MP_TYPE /= 'LSC' .and. MP_TYPE /= 'LSC_MIROC6' ) then
#ifdef _OPENACC
          !$acc kernels pcopy(GDCFRC) pcopyin(q_tot) async(1)
          do k  = 1, kdim
          do ij = 1, ijdim
             q_tot(ij,k) = 0.0_RP

             !$acc loop seq
             do nq = NQW_STR+1, NQW_END
                q_tot(ij,k) = q_tot(ij,k) + q(ij,k,nq)
             enddo

             if ( q_tot(ij,k) >= SMALL ) then
                GDCFRC(ij,k) = 1.0_RP
             else
                GDCFRC(ij,k) = 0.0_RP
             endif

             GDCLW(ij,k) = CNST_UNDEF
          enddo
          enddo
          !$acc end kernels
#else
          !$omp parallel default(none),private(ij,k,nq), &
          !$omp shared(ijdim,kdim,NQW_STR,NQW_END,GDCFRC,q_tot,q,GDCLW,CNST_UNDEF)
!OCL XFILL
          !$omp do
          do k  = 1, kdim
          do ij = 1, ijdim
             q_tot (ij,k) = 0.0_RP
             GDCFRC(ij,k) = 0.0_RP
             GDCLW (ij,k) = CNST_UNDEF
          enddo
          enddo
          !$omp end do

          do nq = NQW_STR+1, NQW_END
             !$omp do
             do k  = 1, kdim
             do ij = 1, ijdim
                q_tot(ij,k) = q_tot(ij,k) + q(ij,k,nq)
             enddo
             enddo
             !$omp end do
          enddo

          !$omp do
          do k  = 1, kdim
          do ij = 1, ijdim
             if( q_tot(ij,k) >= SMALL ) GDCFRC(ij,k) = 1.0_RP
          enddo
          enddo
          !$omp end do

          !$omp end parallel
#endif
       endif
!--- MINIAPP: the MP->CH/AE/RD diagnostic hand-off (cloudvar/radvar/history)
!--- is not part of the validated mp_driver output set and needs mod_comm. ---
#ifndef MINIAPP
       wait_on_cloudvar = .false.
       wait_on_radvar   = .false.

       call cloudvar_setr1_in( GPREC (:,:,1), vid=I_RFXL  , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( GPREC (:,:,2), vid=I_SFXL  , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( Re_qr (:,:)  , vid=I_Re_qr , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( Vt_qr0(:,:)  , vid=I_Vt_qr , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( RWTDL (:,:)  , vid=I_RWTDL , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( IWTDL (:,:)  , vid=I_IWTDL , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( ISEDL (:,:)  , vid=I_ISEDL , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( REVP  (:,:)  , vid=I_REVP  , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( REVS  (:,:)  , vid=I_REVS  , klevs=kdim, l_region=l_region ) ! [IN]
       call cloudvar_setr1_in( GDPFRC(:,:)  , vid=I_GDPFRC, klevs=kdim, l_region=l_region ) ! [IN]
       if ( MP_TYPE == 'LSC_MIROC6' ) then
         call cloudvar_setr1_in( FLIQ  (:,:)  , vid=I_FLIQ  , klevs=kdim, l_region=l_region ) ! [IN]
         call cloudvar_setr1_in( ESWAGG(:,:)  , vid=I_ESWAGG, klevs=kdim, l_region=l_region ) ! [IN]
       endif

       call radvar_setr1_in( Re_liq(:,:), vid=I_RCEFF      , klevs=kdim, l_region=l_region ) ! [IN]
       call radvar_setr1_in( Re_ice(:,:), vid=I_RCEFF_SOLID, klevs=kdim, l_region=l_region ) ! [IN]
       call radvar_setr1_in( Re_qc (:,:), vid=I_RCEFF_CLD  , klevs=kdim, l_region=l_region ) ! [IN]
       call radvar_setr1_in( rctop (:,:), vid=I_RCTOP      , klevs=1   , l_region=l_region ) ! [IN]
       call radvar_setr1_in( rwtop (:,:), vid=I_RWTOP      , klevs=1   , l_region=l_region ) ! [IN]
       call radvar_setr1_in( tctop (:,:), vid=I_TCTOP      , klevs=1   , l_region=l_region ) ! [IN]

       !$acc wait(1)
       wait_on_cloudvar = .true.
       wait_on_radvar   = .true.
       
       call history_in( 'ml_re_liq'   , Re_liq   (:,:) )
       call history_in( 'ml_re_sol'   , Re_ice   (:,:) )
       call history_in( 'sl_rctop_cld', rctop_cld(:,:) )
       call history_in( 'sl_rwtop_cld', rwtop_cld(:,:) )
       call history_in( 'sl_tctop_cld', tctop_cld(:,:) )
#else
       !$acc wait(1)
#endif
    endif

    ! single/double precision conversion adapter

    !$acc kernels async(1)
    do k  = 1, kdim
    do ij = 1, ijdim
       RP_PHY_rhog  (ij,k) = rhog  (ij,k)
       RP_PHY_rhogvx(ij,k) = rhogvx(ij,k)
       RP_PHY_rhogvy(ij,k) = rhogvy(ij,k)
       RP_PHY_rhogvz(ij,k) = rhogvz(ij,k)
       RP_PHY_rhogw (ij,k) = rhogw (ij,k)
       RP_PHY_rhoge (ij,k) = rhoge (ij,k)
       RP_PHY_rho   (ij,k) = rho   (ij,k)
       RP_PHY_pre   (ij,k) = pre   (ij,k)
       RP_PHY_tem   (ij,k) = tem   (ij,k)
       RP_PHY_vx    (ij,k) = vx    (ij,k)
       RP_PHY_vy    (ij,k) = vy    (ij,k)
       RP_PHY_vz    (ij,k) = vz    (ij,k)
       RP_PHY_w     (ij,k) = w     (ij,k)
       RP_PHY_GDCFRC(ij,k) = GDCFRC(ij,k)
       RP_PHY_GDCLW (ij,k) = GDCLW (ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1)
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       RP_PHY_rhogq (ij,k,nq) = rhogq(ij,k,nq)
       RP_PHY_q     (ij,k,nq) = q    (ij,k,nq)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1)
    do ij = 1, ijdim
       RP_PHY_precip       (ij,1) = precip       (ij,1)
       RP_PHY_precip       (ij,2) = precip       (ij,2)
       RP_PHY_precip_ISO1  (ij,1) = precip_ISO1  (ij,1)
       RP_PHY_precip_ISO1  (ij,2) = precip_ISO1  (ij,2)
       RP_PHY_precip_ISO2  (ij,1) = precip_ISO2  (ij,1)
       RP_PHY_precip_ISO2  (ij,2) = precip_ISO2  (ij,2)
       RP_PHY_precip_lh    (ij)   = precip_lh    (ij)
       RP_PHY_precip_rhoein(ij)   = precip_rhoein(ij)
       RP_PHY_precip_rhophi(ij)   = precip_rhophi(ij)
       RP_PHY_precip_rhokin(ij)   = precip_rhokin(ij)
       RP_PHY_precip_energy(ij)   = precip_energy(ij)
    enddo
    !$acc end kernels

    !$acc wait(1)
    !$acc end data
    return
  end subroutine mp_driver

  !-----------------------------------------------------------------------------
  subroutine mp_terminal_velocity( &
       ijdim, kdim,      & ! in
       rho, tem, pre, w, & ! in
       q,                & ! in
       vt_qc, vt_qr, vt_qi, vt_qs, vt_qg, vt_qh,& ! out
       vt_nc, vt_nr, vt_ni, vt_ns, vt_ng, vt_nh ) ! out
    use mod_adm, only: &
         kmin => adm_kmin, &
         kmax => adm_kmax
!--- MINIAPP: only the NSW6 terminal velocity is compiled in. ---
#ifndef MINIAPP
    use mod_mp_g98, only: &
         mp_g98_terminal_velocity
    use mod_mp_nsw5, only: &
         mp_nsw5_terminal_velocity
#endif
    use mod_mp_nsw6, only: &
         mp_nsw6_terminal_velocity
#ifndef MINIAPP
    use mod_mp_ndw6, only: &
         mp_ndw6_terminal_velocity, &
         mp_ndw6_diag_number ! [Add] 09/04/14 T.Mitsui
    use mod_mp_s10, only: &
         mp_s10_terminal_velocity
#endif
    use mod_runconf, only: &
         nqmax=>TRC_VMAX,  &
         MP_TYPE,          &
         opt_2moment_water, & ! [Add] 09/04/14 T.Mitsui
         I_QV,              & ! [Add] 09/04/14 T.Mitsui
         I_QC, I_QR, I_QI, I_QS, I_QG, &
         I_NC, I_NR, I_NI, I_NS, I_NG
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: rho(ijdim,kdim)   ! density
    real(DP), intent(in)  :: tem(ijdim,kdim)   ! temperature
    real(DP), intent(in)  :: pre(ijdim,kdim)   ! pressure
    real(DP), intent(in)  :: w  (ijdim,kdim)   ! vertical velocity
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! mixing ratio of hydrometeors
    real(DP), intent(out) :: vt_qc(ijdim,kdim) ! terminal velocity of cloud mass
    real(DP), intent(out) :: vt_qr(ijdim,kdim) ! terminal velocity of rain mass
    real(DP), intent(out) :: vt_qi(ijdim,kdim) ! terminal velocity of ice mass
    real(DP), intent(out) :: vt_qs(ijdim,kdim) ! terminal velocity of snow mass
    real(DP), intent(out) :: vt_qg(ijdim,kdim) ! terminal velocity of graupel mass
    real(DP), intent(out) :: vt_qh(ijdim,kdim) ! terminal velocity of hail mass
    real(DP), intent(out) :: vt_nc(ijdim,kdim) ! terminal velocity of cloud mass
    real(DP), intent(out) :: vt_nr(ijdim,kdim) ! terminal velocity of rain mass
    real(DP), intent(out) :: vt_ni(ijdim,kdim) ! terminal velocity of ice mass
    real(DP), intent(out) :: vt_ns(ijdim,kdim) ! terminal velocity of snow mass
    real(DP), intent(out) :: vt_ng(ijdim,kdim) ! terminal velocity of graupel mass
    real(DP), intent(out) :: vt_nh(ijdim,kdim) ! terminal velocity of hail mass
    ! [Add] 09/04/14 T.Mitsui
    real(DP) :: qnc(ijdim,kdim)
    real(DP) :: qnr(ijdim,kdim)
    real(DP) :: qni(ijdim,kdim)
    real(DP) :: qns(ijdim,kdim)
    real(DP) :: qng(ijdim,kdim)
    real(DP) :: rrho(ijdim,kdim)
    vt_qc(:,:)=0.0_DP
    vt_qr(:,:)=0.0_DP
    vt_qi(:,:)=0.0_DP
    vt_qs(:,:)=0.0_DP
    vt_qg(:,:)=0.0_DP
    vt_qh(:,:)=0.0_DP
    if ( trim(MP_TYPE) /= 'NDW6' ) then
       vt_nc(:,:)=0.0_DP
       vt_nr(:,:)=0.0_DP
       vt_ni(:,:)=0.0_DP
       vt_ns(:,:)=0.0_DP
       vt_ng(:,:)=0.0_DP
    endif
    vt_nh(:,:)=0.0_DP
    if ( trim(MP_TYPE)=='NONE' ) then
       ! Notice, T.Mitsui never use these mp_schemes.
       !         Please calculate for yourself who want to.
    elseif( trim(MP_TYPE)=='KESSLER') then
    elseif( trim(MP_TYPE)=='LIN') then
    elseif( trim(MP_TYPE)=='WSM6') then
    elseif( trim(MP_TYPE)=='G98') then
!--- MINIAPP: non-NSW6 terminal velocities removed with -DMINIAPP. ---
#ifndef MINIAPP
       call mp_g98_terminal_velocity( &
            ijdim, kdim, & ! in
            rho, tem,    & ! in
            q(:,:,I_QC), q(:,:,I_QR),  & ! in
            vt_qc, vt_qr, vt_qi, vt_qs ) ! out
    elseif( trim(MP_TYPE)=='NSW5') then
       call mp_nsw5_terminal_velocity( &
            ijdim, kdim, & ! in
            rho, tem,    & ! in
            q(:,:,I_QC), q(:,:,I_QR),  & ! in
            q(:,:,I_QI), q(:,:,I_QS),  & ! in
            vt_qc, vt_qr, vt_qi, vt_qs ) ! out
#endif
    elseif( trim(MP_TYPE)=='NSW6') then
       call mp_nsw6_terminal_velocity( &
            ijdim, kdim, & ! in
            rho, tem,    & ! in
            w,    & ! [Add] 2016/06/04 WS.ROH
            q(:,:,I_QC), q(:,:,I_QR),  & ! in
            q(:,:,I_QI), q(:,:,I_QS), q(:,:,I_QG),& ! in
            vt_qc, vt_qr, vt_qi, vt_qs, vt_qg     ) ! out
    elseif( trim(MP_TYPE)=='NDW6' ) then
#ifndef MINIAPP
       ! 09/04/14 T.Mitsui [Add] conditioning
       if (opt_2moment_water) then
          qnc(:,:) = q(:,:,I_NC)
          qnr(:,:) = q(:,:,I_NR)
          qni(:,:) = q(:,:,I_NI)
          qns(:,:) = q(:,:,I_NS)
          qng(:,:) = q(:,:,I_NG)
       else
          call mp_ndw6_diag_number( &
               ijdim, kdim,   & ! in
               kmin, kmax,    & ! in
               rho, tem,      & ! in
               q(:,:,I_QV),   & ! in
               q(:,:,I_QC),   & ! in
               q(:,:,I_QR),   & ! in
               q(:,:,I_QI),   & ! in
               q(:,:,I_QS),   & ! in
               q(:,:,I_QG),   & ! in
               qnc, qnr, qni, qns, qng ) ! out
          rrho(:,:) = 1.0_DP/rho(:,:)
          qnc(:,:) = qnc(:,:)*rrho(:,:)
          qnr(:,:) = qnr(:,:)*rrho(:,:)
          qni(:,:) = qni(:,:)*rrho(:,:)
          qns(:,:) = qns(:,:)*rrho(:,:)
          qng(:,:) = qng(:,:)*rrho(:,:)
       endif
       call mp_ndw6_terminal_velocity( &
            ijdim, kdim,   & ! in
            ! [Mod] 10/08/03 T.Mitsui, add argument 'pre'
!!$         rho, tem,      & ! in
            rho, tem, pre, & ! in
            q(:,:,I_QC), q(:,:,I_QR), q(:,:,I_QI), q(:,:,I_QS), q(:,:,I_QG), & ! in
            ! [Mod] 09/04/14 T.Mitsui
!!$         q(:,:,I_NC), q(:,:,I_NR), q(:,:,I_NI), q(:,:,I_NS), q(:,:,I_NG), & ! in
            qnc(:,:), qnr(:,:), qni(:,:), qns(:,:), qng(:,:), & ! in
            vt_qc, vt_qr, vt_qi, vt_qs, vt_qg, & ! out
            vt_nc, vt_nr, vt_ni, vt_ns, vt_ng  ) ! out
#endif
    elseif( trim(MP_TYPE)=='LSC') then
       vt_qr(:,:)=-5.0_DP ! given in namelist variable as VTERM
       vt_qs(:,:)=-5.0_DP ! given in namelist variable as VTERM
    elseif( trim(MP_TYPE)=='S10') then
#ifndef MINIAPP
       call mp_s10_terminal_velocity( &
            ijdim, kdim, & ! in
            q,           & ! in
            vt_qc,       & ! out
            vt_qr,       & ! out
            vt_qi,       & ! out
            vt_qs,       & ! out
            vt_qg,       & ! out
            vt_qh        ) ! out
#endif
    else
       write(*,*) &
            'Msg : Sub[mp_terminal_velocity]/Mod[mp_driver]'
       write(*,*) &
            ' *** WARNING : Not appropriate type!! CHECK!!'
    endif

    return
  end subroutine mp_terminal_velocity

  !-----------------------------------------------------------------------------
  ! Temporarily values are imported from each subroutine and set density as constant.
  ! Please prepare subroutines for yourself to synchronize parameters with cloud microphysics.
  !--- MINIAPP: mp_diag_volume / mp_effective_radius are MP->RD diagnostics.  physics_step calls them AFTER mp_driver; they are not part of the validated mp_driver output set and they pull in every other MP scheme plus mod_history / mod_cloudvar.  Preserved, commented out. ---
!MINIAPP  subroutine mp_diag_volume( &
!MINIAPP       ijdim,         &
!MINIAPP       kmin,          &
!MINIAPP       kmax,          &
!MINIAPP       kdim,          &
!MINIAPP       l_region,      &
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY)
!MINIAPP       rho,           &
!MINIAPP       tem,           &
!MINIAPP       pre,           &
!MINIAPP       q,             &
!MINIAPP       qc_cum,        &
!MINIAPP       cfrac,         &
!MINIAPP       cfrac_cum,     &
!MINIAPP       cloud_volume,  &
!MINIAPP       cumulus_volume )
!MINIAPP#else
!MINIAPP       RP_PHY_rho,           &
!MINIAPP       RP_PHY_tem,           &
!MINIAPP       RP_PHY_pre,           &
!MINIAPP       RP_PHY_q,             &
!MINIAPP       RP_PHY_qc_cum,        &
!MINIAPP       RP_PHY_cfrac,         &
!MINIAPP       RP_PHY_cfrac_cum,     &
!MINIAPP       RP_PHY_cloud_volume,  &
!MINIAPP       RP_PHY_cumulus_volume )
!MINIAPP#endif
!MINIAPP    use mod_cnst, only: &
!MINIAPP       DWATR => CNST_DWATR, &
!MINIAPP       DICE  => CNST_DICE
!MINIAPP    use mod_runconf, only: &
!MINIAPP       MP_TYPE,                   &
!MINIAPP       opt_2moment_water,         &
!MINIAPP       nqmax=>TRC_VMAX,           &
!MINIAPP       Nparticle_mp,              &
!MINIAPP       NQW_STR, NQW_END,          &
!MINIAPP       I_QV,                      &
!MINIAPP       I_QC, I_QR, I_QI,          &
!MINIAPP       I_QS, I_QG, I_QH,          &
!MINIAPP       I_NC, I_NR, I_NI,          &
!MINIAPP       I_NS, I_NG,                &
!MINIAPP       I_mp_QC, I_mp_QR, I_mp_QI, &
!MINIAPP       I_mp_QS, I_mp_QG, I_mp_QH
!MINIAPP    use mod_mp_lsc, only: &
!MINIAPP       mp_lsc_diag_volume
!MINIAPP    use mod_mp_lsc_miroc6, only: &
!MINIAPP       mp_lsc_miroc6_diag_volume
!MINIAPP    use mod_mp_g98, only: &
!MINIAPP       mp_g98_diag_volume
!MINIAPP    use mod_mp_ndw6, only: &
!MINIAPP       mp_ndw6_diag_volume, &
!MINIAPP       mp_ndw6_diag_number
!MINIAPP    use mod_mp_s10, only: &
!MINIAPP       mp_s10_diag_volume
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer,  intent(in)  :: ijdim
!MINIAPP    integer,  intent(in)  :: kmin
!MINIAPP    integer,  intent(in)  :: kmax
!MINIAPP    integer,  intent(in)  :: kdim
!MINIAPP    integer,  intent(in)  :: l_region
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY)
!MINIAPP    real(DP), intent(in)  :: q             (ijdim,kdim,nqmax)        ! [kg/kg] mixing ratio
!MINIAPP    real(DP), intent(in)  :: qc_cum        (ijdim,kdim)              ! [kg/kg] cumulus mixing ratio
!MINIAPP    real(DP), intent(in)  :: cfrac         (ijdim,kdim)              ! [0:1] cloud fraction
!MINIAPP    real(DP), intent(in)  :: cfrac_cum     (ijdim,1)                 ! [0:1] cumulus cloud fraction
!MINIAPP    real(DP), intent(in)  :: rho           (ijdim,kdim)              ! [kg/m3] air density
!MINIAPP    real(DP), intent(in)  :: tem           (ijdim,kdim)              ! [K] temperature
!MINIAPP    real(DP), intent(in)  :: pre           (ijdim,kdim)              ! [Pa] Pressure
!MINIAPP    real(DP), intent(out) :: cloud_volume  (ijdim,kdim,Nparticle_mp) ! [m3/m3]
!MINIAPP    real(DP), intent(out) :: cumulus_volume(ijdim,kdim,2)            ! [m3/m3]
!MINIAPP#else
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_q             (ijdim,kdim,nqmax)        ! [kg/kg] mixing ratio
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_qc_cum        (ijdim,kdim)              ! [kg/kg] cumulus mixing ratio
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_cfrac         (ijdim,kdim)              ! [0:1] cloud fraction
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_cfrac_cum     (ijdim,1)                 ! [0:1] cumulus cloud fraction
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_rho           (ijdim,kdim)              ! [kg/m3] air density
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_tem           (ijdim,kdim)              ! [K] temperature
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_pre           (ijdim,kdim)              ! [Pa] Pressure
!MINIAPP    real(RP_PHY), intent(out) :: RP_PHY_cloud_volume  (ijdim,kdim,Nparticle_mp) ! [m3/m3]
!MINIAPP    real(RP_PHY), intent(out) :: RP_PHY_cumulus_volume(ijdim,kdim,2)            ! [m3/m3]
!MINIAPP    real(DP) :: q             (ijdim,kdim,nqmax)        ! [kg/kg] mixing ratio
!MINIAPP    real(DP) :: qc_cum        (ijdim,kdim)              ! [kg/kg] cumulus mixing ratio
!MINIAPP    real(DP) :: cfrac         (ijdim,kdim)              ! [0:1] cloud fraction
!MINIAPP    real(DP) :: cfrac_cum     (ijdim,1)                 ! [0:1] cumulus cloud fraction
!MINIAPP    real(DP) :: rho           (ijdim,kdim)              ! [kg/m3] air density
!MINIAPP    real(DP) :: tem           (ijdim,kdim)              ! [K] temperature
!MINIAPP    real(DP) :: pre           (ijdim,kdim)              ! [Pa] Pressure
!MINIAPP    real(DP) :: cloud_volume  (ijdim,kdim,Nparticle_mp) ! [m3/m3]
!MINIAPP    real(DP) :: cumulus_volume(ijdim,kdim,2)            ! [m3/m3]
!MINIAPP#endif
!MINIAPP    real(DP) :: nc(ijdim,kdim)
!MINIAPP    real(DP) :: nr(ijdim,kdim)
!MINIAPP    real(DP) :: ni(ijdim,kdim)
!MINIAPP    real(DP) :: ns(ijdim,kdim)
!MINIAPP    real(DP) :: ng(ijdim,kdim)
!MINIAPP    real(DP) :: FLIQ
!MINIAPP
!MINIAPP    integer  :: ij, k, nq
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data &
!MINIAPP    !$acc pcreate(cloud_volume,cumulus_volume) &
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY)
!MINIAPP    !$acc pcopyin(q,qc_cum,cfrac,cfrac_cum,rho,tem,pre)
!MINIAPP#else
!MINIAPP    !$acc pcopyin(RP_PHY_q,RP_PHY_qc_cum,RP_PHY_cfrac,RP_PHY_cfrac_cum,RP_PHY_rho,RP_PHY_tem,RP_PHY_pre) &
!MINIAPP    !$acc pcopyout(RP_PHY_cloud_volume,RP_PHY_cumulus_volume) &
!MINIAPP    !$acc pcreate(q,qc_cum,cfrac,cfrac_cum,rho,tem,pre)
!MINIAPP#endif
!MINIAPP
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY)
!MINIAPP    ! do nothing
!MINIAPP#else
!MINIAPP    !$acc kernels
!MINIAPP    do k  = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       rho(ij,k)     = real(RP_PHY_rho   (ij,k), DP)
!MINIAPP       tem(ij,k)     = real(RP_PHY_tem   (ij,k), DP)
!MINIAPP       pre(ij,k)     = real(RP_PHY_pre   (ij,k), DP)
!MINIAPP       qc_cum (ij,k) = real(RP_PHY_qc_cum(ij,k), DP)
!MINIAPP       cfrac  (ij,k) = real(RP_PHY_cfrac (ij,k), DP)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels
!MINIAPP    do nq = 1, nqmax
!MINIAPP    do k  = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       q(ij,k,nq) = real(RP_PHY_q(ij,k,nq), DP)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels
!MINIAPP    do ij = 1, ijdim
!MINIAPP      cfrac_cum(ij,1) = real(RP_PHY_cfrac_cum(ij,1), DP)
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP#endif
!MINIAPP    if ( opt_volume_explicit ) then
!MINIAPP
!MINIAPP       if ( MP_TYPE == 'NDW6' ) then
!MINIAPP
!MINIAPP          if ( opt_2moment_water ) then
!MINIAPP             !$omp parallel do default(none),private(ij,k), &
!MINIAPP             !$omp shared(ijdim,kdim,nc,nr,ni,ns,ng,rho,q,I_NC,I_NR,I_NI,I_NS,I_NG)
!MINIAPP             !$acc kernels pcopy(nc,nr,ni,ns,ng) pcopyin(rho,q)
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                nc(ij,k) = rho(ij,k) * q(ij,k,I_NC)
!MINIAPP                nr(ij,k) = rho(ij,k) * q(ij,k,I_NR)
!MINIAPP                ni(ij,k) = rho(ij,k) * q(ij,k,I_NI)
!MINIAPP                ns(ij,k) = rho(ij,k) * q(ij,k,I_NS)
!MINIAPP                ng(ij,k) = rho(ij,k) * q(ij,k,I_NG)
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end parallel do
!MINIAPP             !$acc end kernels
!MINIAPP          else
!MINIAPP             call mp_ndw6_diag_number( ijdim,       & ! [IN]
!MINIAPP                                       kdim,        & ! [IN]
!MINIAPP                                       kmin,        & ! [IN]
!MINIAPP                                       kmax,        & ! [IN]
!MINIAPP                                       rho(:,:),    & ! [IN]
!MINIAPP                                       tem(:,:),    & ! [IN]
!MINIAPP                                       q(:,:,I_QV), & ! [IN]
!MINIAPP                                       q(:,:,I_QC), & ! [IN]
!MINIAPP                                       q(:,:,I_QR), & ! [IN]
!MINIAPP                                       q(:,:,I_QI), & ! [IN]
!MINIAPP                                       q(:,:,I_QS), & ! [IN]
!MINIAPP                                       q(:,:,I_QG), & ! [IN]
!MINIAPP                                       nc(:,:),     & ! [OUT]
!MINIAPP                                       nr(:,:),     & ! [OUT]
!MINIAPP                                       ni(:,:),     & ! [OUT]
!MINIAPP                                       ns(:,:),     & ! [OUT]
!MINIAPP                                       ng(:,:)      ) ! [OUT]
!MINIAPP          endif
!MINIAPP
!MINIAPP          call mp_ndw6_diag_volume( ijdim,                     & ! [IN]
!MINIAPP                                    kdim,                      & ! [IN]
!MINIAPP                                    kmin,                      & ! [IN]
!MINIAPP                                    kmax,                      & ! [IN]
!MINIAPP                                    rho         (:,:),         & ! [IN]
!MINIAPP                                    nc          (:,:),         & ! [IN]
!MINIAPP                                    nr          (:,:),         & ! [IN]
!MINIAPP                                    ni          (:,:),         & ! [IN]
!MINIAPP                                    ns          (:,:),         & ! [IN]
!MINIAPP                                    ng          (:,:),         & ! [IN]
!MINIAPP                                    q           (:,:,I_QC),    & ! [IN]
!MINIAPP                                    q           (:,:,I_QR),    & ! [IN]
!MINIAPP                                    q           (:,:,I_QI),    & ! [IN]
!MINIAPP                                    q           (:,:,I_QS),    & ! [IN]
!MINIAPP                                    q           (:,:,I_QG),    & ! [IN]
!MINIAPP                                    cloud_volume(:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                    cloud_volume(:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                    cloud_volume(:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                    cloud_volume(:,:,I_mp_QS), & ! [OUT]
!MINIAPP                                    cloud_volume(:,:,I_mp_QG)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NSW6' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume,rho,q,I_QC,I_QR,I_QI,I_QS,I_QG)
!MINIAPP          !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,I_QC) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = rho(ij,k) * q(ij,k,I_QR) / DWATR
!MINIAPP             ! [Fix] 2016/05/13 T.Seiki: here we should use optically effective volume
!MINIAPP             !                           see Fu etal.(1996),JC or Seiki etal.(2014),JGR
!MINIAPP             !cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,I_QI) / 400.0_DP ! rho_g
!MINIAPP             !cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,I_QS) / 100.0_DP ! rho_s
!MINIAPP             !cloud_volume(ij,k,I_mp_QG) = rho(ij,k) * q(ij,k,I_QG) / 400.0_DP ! rho_g
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,I_QI) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,I_QS) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = rho(ij,k) * q(ij,k,I_QG) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NSW5' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume,rho,q,I_QC,I_QR,I_QI,I_QS)
!MINIAPP          !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,I_QC) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = rho(ij,k) * q(ij,k,I_QR) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,I_QI) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,I_QS) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'WSM6' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume,rho,q,I_QC,I_QR,I_QI,I_QS,I_QG)
!MINIAPP          !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,I_QC) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = rho(ij,k) * q(ij,k,I_QR) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,I_QI) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,I_QS) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = rho(ij,k) * q(ij,k,I_QG) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'LIN' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume,rho,q,I_QC,I_QR,I_QI,I_QS,I_QG)
!MINIAPP          !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,I_QC) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = rho(ij,k) * q(ij,k,I_QR) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,I_QI) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,I_QS) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = rho(ij,k) * q(ij,k,I_QG) / DICE
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'KESSLER' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume,rho,q,I_QC)
!MINIAPP          !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,I_QC) / DWATR
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'G98' ) then
!MINIAPP
!MINIAPP          call mp_g98_diag_volume( ijdim,                     & ! [IN]
!MINIAPP                                   kdim,                      & ! [IN]
!MINIAPP                                   kmin,                      & ! [IN]
!MINIAPP                                   kmax,                      & ! [IN]
!MINIAPP                                   rho         (:,:),         & ! [IN]
!MINIAPP                                   tem         (:,:),         & ! [IN]
!MINIAPP                                   pre         (:,:),         & ! [IN]
!MINIAPP                                   q           (:,:,I_QC),    & ! [IN]
!MINIAPP                                   q           (:,:,I_QR),    & ! [IN]
!MINIAPP                                   cloud_volume(:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QS)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'LSC' ) then
!MINIAPP
!MINIAPP          call mp_lsc_diag_volume( ijdim,                     & ! [IN]
!MINIAPP                                   kdim,                      & ! [IN]
!MINIAPP                                   kmin,                      & ! [IN]
!MINIAPP                                   kmax,                      & ! [IN]
!MINIAPP                                   rho         (:,:),         & ! [IN]
!MINIAPP                                   tem         (:,:),         & ! [IN]
!MINIAPP                                   q           (:,:,I_QC),    & ! [IN]
!MINIAPP                                   cfrac       (:,:),         & ! [IN]
!MINIAPP                                   cloud_volume(:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QI)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP       elseif( MP_TYPE == 'LSC_MIROC6' ) then
!MINIAPP
!MINIAPP          call mp_lsc_miroc6_diag_volume( ijdim,                     & ! [IN]
!MINIAPP                                   kdim,                      & ! [IN]
!MINIAPP                                   kmin,                      & ! [IN]
!MINIAPP                                   kmax,                      & ! [IN]
!MINIAPP                                   rho         (:,:),         & ! [IN]
!MINIAPP                                   q           (:,:,I_QC),    & ! [IN]
!MINIAPP                                   q           (:,:,I_QI),    & ! [IN]
!MINIAPP                                   cfrac       (:,:),         & ! [IN]
!MINIAPP                                   cloud_volume(:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QI)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'S10' ) then
!MINIAPP
!MINIAPP          call mp_s10_diag_volume( ijdim,                     & ! [IN]
!MINIAPP                                   kdim,                      & ! [IN]
!MINIAPP                                   kmin,                      & ! [IN]
!MINIAPP                                   kmax,                      & ! [IN]
!MINIAPP                                   rho         (:,:),         & ! [IN]
!MINIAPP                                   q           (:,:,:),       & ! [IN]
!MINIAPP                                   cloud_volume(:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QS), & ! [OUT]
!MINIAPP                                   cloud_volume(:,:,I_mp_QG)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NONE' ) then
!MINIAPP
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,cloud_volume)
!MINIAPP          !$acc kernels pcopy(cloud_volume)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             cloud_volume(ij,k,I_mp_QC) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QI) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             cloud_volume(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       endif
!MINIAPP
!MINIAPP    else ! set default
!MINIAPP
!MINIAPP       !$omp parallel default(none),private(ij,k,nq), &
!MINIAPP       !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cloud_volume,rho,q, &
!MINIAPP       !$omp        I_QC,I_QR,I_QI,I_QS,I_QG,I_QH                  )
!MINIAPP
!MINIAPP       !$acc kernels pcopy(cloud_volume)
!MINIAPP       !$omp do
!MINIAPP       do nq = 1, Nparticle_mp
!MINIAPP       do k  = 1, kdim
!MINIAPP       do ij = 1, ijdim
!MINIAPP          cloud_volume(ij,k,nq) = 0.0_DP
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do
!MINIAPP       !$acc end kernels
!MINIAPP
!MINIAPP       !$acc kernels pcopy(cloud_volume) pcopyin(rho,q)
!MINIAPP       do nq = NQW_STR+1, NQW_END
!MINIAPP          if ( nq == I_QC ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QC) = rho(ij,k) * q(ij,k,nq) / DWATR
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          elseif( nq == I_QR ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QR) = rho(ij,k) * q(ij,k,nq) / DWATR
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          elseif( nq == I_QI ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QI) = rho(ij,k) * q(ij,k,nq) / DICE
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          elseif( nq == I_QS ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QS) = rho(ij,k) * q(ij,k,nq) / DICE
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          elseif( nq == I_QG ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QG) = rho(ij,k) * q(ij,k,nq) / DICE
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          elseif( nq == I_QH ) then
!MINIAPP             !$omp do
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                cloud_volume(ij,k,I_mp_QH) = rho(ij,k) * q(ij,k,nq) / DICE
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end do
!MINIAPP          endif
!MINIAPP       enddo
!MINIAPP       !$omp end parallel
!MINIAPP       !$acc end kernels
!MINIAPP
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! cumulus_volume with liquid/ice separation
!MINIAPP
!MINIAPP    !$omp parallel default(none),private(ij,k,FLIQ), &
!MINIAPP    !$omp shared(ijdim,kmin,kmax,cumulus_volume,tem,rho,qc_cum,cfrac_cum,TWICE,TSICE)
!MINIAPP
!MINIAPP    !$acc kernels pcopy(cumulus_volume) pcopyin(tem,rho,qc_cum,cfrac_cum)
!MINIAPP    !$omp do
!MINIAPP    do k  = kmin, kmax
!MINIAPP    do ij = 1, ijdim
!MINIAPP       FLIQ = ( tem(ij,k)-TWICE ) / ( TSICE-TWICE )
!MINIAPP       FLIQ = min( max( FLIQ, 0.0_DP ), 1.0_DP )
!MINIAPP
!MINIAPP       if ( cfrac_cum(ij,1) >= 1.E-12_DP ) then
!MINIAPP          cumulus_volume(ij,k,1) = rho(ij,k) * (        FLIQ ) * qc_cum(ij,k) / ( cfrac_cum(ij,1) * DWATR )
!MINIAPP          cumulus_volume(ij,k,2) = rho(ij,k) * ( 1.0_DP-FLIQ ) * qc_cum(ij,k) / ( cfrac_cum(ij,1) * DICE  )
!MINIAPP       else
!MINIAPP          cumulus_volume(ij,k,1) = 0.0_DP
!MINIAPP          cumulus_volume(ij,k,2) = 0.0_DP
!MINIAPP       endif
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$omp end do
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels pcopy(cumulus_volume)
!MINIAPP    !$omp do
!MINIAPP    do ij = 1, ijdim
!MINIAPP       cumulus_volume(ij,kmin-1,1) = 0.0_DP
!MINIAPP       cumulus_volume(ij,kmin-1,2) = 0.0_DP
!MINIAPP       cumulus_volume(ij,kmax+1,1) = 0.0_DP
!MINIAPP       cumulus_volume(ij,kmax+1,2) = 0.0_DP
!MINIAPP    enddo
!MINIAPP    !$omp end do
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$omp end parallel
!MINIAPP
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY)
!MINIAPP    ! do nothing
!MINIAPP#else
!MINIAPP    !$acc kernels
!MINIAPP    do k  = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       RP_PHY_cumulus_volume(ij,k,1) = real(cumulus_volume(ij,k,1), RP_PHY)
!MINIAPP       RP_PHY_cumulus_volume(ij,k,2) = real(cumulus_volume(ij,k,2), RP_PHY)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels
!MINIAPP    do nq = 1, Nparticle_mp
!MINIAPP    do k  = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       RP_PHY_cloud_volume(ij,k,nq) = real(cloud_volume(ij,k,nq), RP_PHY)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP#endif
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine mp_diag_volume

  !-----------------------------------------------------------------------------
  !--- MINIAPP: mp_diag_volume / mp_effective_radius are MP->RD diagnostics.  physics_step calls them AFTER mp_driver; they are not part of the validated mp_driver output set and they pull in every other MP scheme plus mod_history / mod_cloudvar.  Preserved, commented out. ---
!MINIAPP  subroutine mp_effective_radius( &
!MINIAPP       ijdim,        &
!MINIAPP       kmin,         &
!MINIAPP       kmax,         &
!MINIAPP       kdim,         &
!MINIAPP       l_region,     &
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP)
!MINIAPP       q,            &
!MINIAPP       GDCFRC,       &
!MINIAPP       rho,          &
!MINIAPP       tem,          &
!MINIAPP       pre,          &
!MINIAPP       ccn,          &
!MINIAPP       factor_n_aie, &
!MINIAPP       w,            &
!MINIAPP       re            )
!MINIAPP#else
!MINIAPP       RP_PHY_q,            &
!MINIAPP       RP_PHY_GDCFRC,       &
!MINIAPP       RP_PHY_rho,          &
!MINIAPP       RP_PHY_tem,          &
!MINIAPP       RP_PHY_pre,          &
!MINIAPP       RP_PHY_ccn,          &
!MINIAPP       RP_PHY_factor_n_aie, &
!MINIAPP       RP_PHY_w,            &
!MINIAPP       RP_PHY_re            )
!MINIAPP#endif
!MINIAPP    use mod_runconf, only: &
!MINIAPP       nqmax => TRC_VMAX,         &
!MINIAPP       MP_TYPE,                   &
!MINIAPP       opt_2moment_water,         &
!MINIAPP       Nparticle_mp,              &
!MINIAPP       I_QV,                      &
!MINIAPP       I_QC, I_QR, I_QI,          &
!MINIAPP       I_QS, I_QG,                &
!MINIAPP       I_NC, I_NR, I_NI,          &
!MINIAPP       I_NS, I_NG,                &
!MINIAPP       I_mp_QC, I_mp_QR, I_mp_QI, &
!MINIAPP       I_mp_QS, I_mp_QG, I_mp_QH
!MINIAPP    use mod_mp_lsc, only: &
!MINIAPP       mp_lsc_effective_radius
!MINIAPP    use mod_mp_lsc_miroc6, only: &
!MINIAPP       mp_lsc_miroc6_effective_radius
!MINIAPP    use mod_mp_g98, only: &
!MINIAPP       mp_g98_effective_radius
!MINIAPP    use mod_mp_nsw5, only: &
!MINIAPP       mp_nsw5_effective_radius
!MINIAPP    use mod_mp_nsw6, only: &
!MINIAPP       mp_nsw6_effective_radius
!MINIAPP    use mod_mp_ndw6, only: &
!MINIAPP       mp_ndw6_effective_radius, &
!MINIAPP       mp_ndw6_diag_number
!MINIAPP    use mod_mp_s10, only: &
!MINIAPP       mp_s10_effective_radius
!MINIAPP    use mod_history, only: &
!MINIAPP       history_in
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer,  intent(in)  :: ijdim
!MINIAPP    integer,  intent(in)  :: kmin
!MINIAPP    integer,  intent(in)  :: kmax
!MINIAPP    integer,  intent(in)  :: kdim
!MINIAPP    integer,  intent(in)  :: l_region
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP)
!MINIAPP    real(RP), intent(in)  :: q           (ijdim,kdim,nqmax)        ! mixing ratio [kg/kg]
!MINIAPP    real(RP), intent(in)  :: GDCFRC      (ijdim,kdim)              ! cloud fraction [0-1]
!MINIAPP    real(RP), intent(in)  :: rho         (ijdim,kdim)              ! air density [kg/m3]
!MINIAPP    real(RP), intent(in)  :: tem         (ijdim,kdim)              ! temperature [K]
!MINIAPP    real(RP), intent(in)  :: pre         (ijdim,kdim)              ! pressure [Pa]
!MINIAPP    real(RP), intent(in)  :: ccn         (ijdim,kdim)              ! cloud condensation nuclei [1/m3]
!MINIAPP    real(RP), intent(in)  :: factor_n_aie(ijdim)                   ! ant. indirect effect for NSW6 [add] 2017.11.27 C.Kodama
!MINIAPP    real(RP), intent(in)  :: w           (ijdim,kdim)              ! vertical velocity [m/s]
!MINIAPP    real(RP), intent(out) :: re          (ijdim,kdim,Nparticle_mp) ! effective radius [m]
!MINIAPP#else
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_q           (ijdim,kdim,nqmax)        ! mixing ratio [kg/kg]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_GDCFRC      (ijdim,kdim)              ! cloud fraction [0-1]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_rho         (ijdim,kdim)              ! air density [kg/m3]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_tem         (ijdim,kdim)              ! temperature [K]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_pre         (ijdim,kdim)              ! pressure [Pa]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_ccn         (ijdim,kdim)              ! cloud condensation nuclei [1/m3]
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_factor_n_aie(ijdim)                   ! ant. indirect effect for NSW6 [add] 2017.11.27 C.Kodama
!MINIAPP    real(RP_PHY), intent(in)  :: RP_PHY_w           (ijdim,kdim)              ! vertical velocity [m/s]
!MINIAPP    real(RP_PHY), intent(out) :: RP_PHY_re          (ijdim,kdim,Nparticle_mp) ! effective radius [m]
!MINIAPP
!MINIAPP    real(DP) :: q           (ijdim,kdim,nqmax)        ! mixing ratio [kg/kg]
!MINIAPP    real(DP) :: GDCFRC      (ijdim,kdim)              ! cloud fraction [0-1]
!MINIAPP    real(DP) :: rho         (ijdim,kdim)              ! air density [kg/m3]
!MINIAPP    real(DP) :: tem         (ijdim,kdim)              ! temperature [K]
!MINIAPP    real(DP) :: pre         (ijdim,kdim)              ! pressure [Pa]
!MINIAPP    real(DP) :: ccn         (ijdim,kdim)              ! cloud condensation nuclei [1/m3]
!MINIAPP    real(DP) :: factor_n_aie(ijdim)                   ! ant. indirect effect for NSW6 [add] 2017.11.27 C.Kodama
!MINIAPP    real(DP) :: w           (ijdim,kdim)              ! vertical velocity [m/s]
!MINIAPP    real(DP) :: re          (ijdim,kdim,Nparticle_mp) ! effective radius [m]
!MINIAPP#endif
!MINIAPP    ! number concentration[1/m3]
!MINIAPP    real(DP) :: nc(ijdim,kdim)
!MINIAPP    real(DP) :: nr(ijdim,kdim)
!MINIAPP    real(DP) :: ni(ijdim,kdim)
!MINIAPP    real(DP) :: ns(ijdim,kdim)
!MINIAPP    real(DP) :: ng(ijdim,kdim)
!MINIAPP
!MINIAPP    integer  :: ij, k, nq
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data &
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP)
!MINIAPP    !$acc pcreate(re) &
!MINIAPP    !$acc pcopyin(q,GDCFRC,rho,tem,pre,ccn,factor_n_aie,w) &
!MINIAPP#else
!MINIAPP    !$acc pcopyin(RP_PHY_q,RP_PHY_GDCFRC,RP_PHY_rho,RP_PHY_tem,RP_PHY_pre, &
!MINIAPP    !$acc          RP_PHY_ccn,RP_PHY_factor_n_aie,RP_PHY_w) &
!MINIAPP    !$acc pcopyout(RP_PHY_re) &
!MINIAPP    !$acc pcreate(q,GDCFRC,rho,tem,pre,ccn,factor_n_aie,w,re) &
!MINIAPP#endif
!MINIAPP    !$acc pcreate(nc,nr,ni,ns,ng)
!MINIAPP
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP)
!MINIAPP    ! do nothing
!MINIAPP#else
!MINIAPP    !$acc kernels async(1)
!MINIAPP    do nq = 1, nqmax
!MINIAPP    do k  = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       q(ij,k,nq)            = real(RP_PHY_q           (ij,k,nq), DP)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1)
!MINIAPP    do k = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       GDCFRC       (ij,k)   = real(RP_PHY_GDCFRC      (ij,k), DP)
!MINIAPP       rho          (ij,k)   = real(RP_PHY_rho         (ij,k), DP)
!MINIAPP       tem          (ij,k)   = real(RP_PHY_tem         (ij,k), DP)
!MINIAPP       pre          (ij,k)   = real(RP_PHY_pre         (ij,k), DP)
!MINIAPP       ccn          (ij,k)   = real(RP_PHY_ccn         (ij,k), DP)
!MINIAPP       w            (ij,k)   = real(RP_PHY_w           (ij,k), DP)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1)
!MINIAPP    do ij = 1, ijdim
!MINIAPP       factor_n_aie(ij)      = real(RP_PHY_factor_n_aie(ij), DP)
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP    !$acc wait(1)
!MINIAPP#endif
!MINIAPP    if ( opt_radius_explicit ) then
!MINIAPP
!MINIAPP       if ( MP_TYPE == 'NDW6' ) then
!MINIAPP          if ( opt_2moment_water ) then
!MINIAPP             !$acc kernels pcopy(nc,nr,ni,ns,ng) pcopyin(rho,q)
!MINIAPP             !$omp parallel do default(none),private(ij,k), &
!MINIAPP             !$omp shared(ijdim,kdim,nc,nr,ni,ns,ng,rho,q,I_NC,I_NR,I_NI,I_NS,I_NG)
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                nc(ij,k) = rho(ij,k) * q(ij,k,I_NC)
!MINIAPP                nr(ij,k) = rho(ij,k) * q(ij,k,I_NR)
!MINIAPP                ni(ij,k) = rho(ij,k) * q(ij,k,I_NI)
!MINIAPP                ns(ij,k) = rho(ij,k) * q(ij,k,I_NS)
!MINIAPP                ng(ij,k) = rho(ij,k) * q(ij,k,I_NG)
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end parallel do
!MINIAPP             !$acc end kernels
!MINIAPP          else
!MINIAPP             call mp_ndw6_diag_number( ijdim,       & ! [IN]
!MINIAPP                                       kdim,        & ! [IN]
!MINIAPP                                       kmin,        & ! [IN]
!MINIAPP                                       kmax,        & ! [IN]
!MINIAPP                                       rho(:,:),    & ! [IN]
!MINIAPP                                       tem(:,:),    & ! [IN]
!MINIAPP                                       q(:,:,I_QV), & ! [IN]
!MINIAPP                                       q(:,:,I_QC), & ! [IN]
!MINIAPP                                       q(:,:,I_QR), & ! [IN]
!MINIAPP                                       q(:,:,I_QI), & ! [IN]
!MINIAPP                                       q(:,:,I_QS), & ! [IN]
!MINIAPP                                       q(:,:,I_QG), & ! [IN]
!MINIAPP                                       nc(:,:),     & ! [OUT]
!MINIAPP                                       nr(:,:),     & ! [OUT]
!MINIAPP                                       ni(:,:),     & ! [OUT]
!MINIAPP                                       ns(:,:),     & ! [OUT]
!MINIAPP                                       ng(:,:)      ) ! [OUT]
!MINIAPP          endif
!MINIAPP
!MINIAPP          call mp_ndw6_effective_radius( ijdim,            & ! [IN]
!MINIAPP                                         kdim,             & ! [IN]
!MINIAPP                                         kmin,             & ! [IN]
!MINIAPP                                         kmax,             & ! [IN]
!MINIAPP                                         rho(:,:),         & ! [IN]
!MINIAPP                                         tem(:,:),         & ! [IN]
!MINIAPP                                         pre(:,:),         & ! [IN]
!MINIAPP                                         nc (:,:),         & ! [IN]
!MINIAPP                                         nr (:,:),         & ! [IN]
!MINIAPP                                         ni (:,:),         & ! [IN]
!MINIAPP                                         ns (:,:),         & ! [IN]
!MINIAPP                                         ng (:,:),         & ! [IN]
!MINIAPP                                         q  (:,:,I_QC),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QR),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QI),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QS),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QG),    & ! [IN]
!MINIAPP                                         re (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QS), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QG)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NSW6' ) then
!MINIAPP
!MINIAPP          call mp_nsw6_effective_radius( ijdim,                     & ! [IN]
!MINIAPP                                         kdim,                      & ! [IN]
!MINIAPP                                         kmin,                      & ! [IN]
!MINIAPP                                         kmax,                      & ! [IN]
!MINIAPP                                         rho         (:,:),         & ! [IN]
!MINIAPP                                         ccn         (:,:),         & ! [IN]
!MINIAPP                                         factor_n_aie(:),           & ! [IN] [add] C.Kodama 2017/11/27
!MINIAPP                                         q           (:,:,I_QC),    & ! [IN]
!MINIAPP                                         q           (:,:,I_QR),    & ! [IN]
!MINIAPP                                         q           (:,:,I_QI),    & ! [IN]
!MINIAPP                                         q           (:,:,I_QS),    & ! [IN]
!MINIAPP                                         q           (:,:,I_QG),    & ! [IN]
!MINIAPP                                         re          (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                         re          (:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                         re          (:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                         re          (:,:,I_mp_QS), & ! [OUT]
!MINIAPP                                         re          (:,:,I_mp_QG)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NSW5' ) then
!MINIAPP
!MINIAPP          call mp_nsw5_effective_radius( ijdim,            & ! [IN]
!MINIAPP                                         kdim,             & ! [IN]
!MINIAPP                                         kmin,             & ! [IN]
!MINIAPP                                         kmax,             & ! [IN]
!MINIAPP                                         rho(:,:),         & ! [IN]
!MINIAPP                                         tem(:,:),         & ! [IN]
!MINIAPP                                         pre(:,:),         & ! [IN]
!MINIAPP                                         ccn(:,:),         & ! [IN]
!MINIAPP                                         q  (:,:,I_QC),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QR),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QI),    & ! [IN]
!MINIAPP                                         q  (:,:,I_QS),    & ! [IN]
!MINIAPP                                         re (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                         re (:,:,I_mp_QS)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'WSM6' ) then
!MINIAPP
!MINIAPP          ! Tentative, constant value
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QC) =   8.E-6_DP
!MINIAPP             re(ij,k,I_mp_QR) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QI) =  20.E-6_DP
!MINIAPP             re(ij,k,I_mp_QS) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QG) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QH) =   0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'LIN' ) then
!MINIAPP
!MINIAPP          ! Tentative, constant value
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QC) =   8.E-6_DP
!MINIAPP             re(ij,k,I_mp_QR) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QI) =  20.E-6_DP
!MINIAPP             re(ij,k,I_mp_QS) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QG) = 100.E-6_DP
!MINIAPP             re(ij,k,I_mp_QH) =   0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'KESSLER' ) then
!MINIAPP
!MINIAPP          ! Tentative, constant value
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QC) = 8.E-6_DP
!MINIAPP             re(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QI) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'G98' ) then
!MINIAPP
!MINIAPP          call mp_g98_effective_radius( ijdim,            & ! [IN]
!MINIAPP                                        kdim,             & ! [IN]
!MINIAPP                                        kmin,             & ! [IN]
!MINIAPP                                        kmax,             & ! [IN]
!MINIAPP                                        rho(:,:),         & ! [IN]
!MINIAPP                                        tem(:,:),         & ! [IN]
!MINIAPP                                        pre(:,:),         & ! [IN]
!MINIAPP                                        ccn(:,:),         & ! [IN]
!MINIAPP                                        q  (:,:,I_QC),    & ! [IN]
!MINIAPP                                        q  (:,:,I_QR),    & ! [IN]
!MINIAPP                                        re (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                        re (:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                        re (:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                        re (:,:,I_mp_QS)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'LSC' ) then
!MINIAPP
!MINIAPP          call mp_lsc_effective_radius( ijdim,               & ! [IN]
!MINIAPP                                        kdim,                & ! [IN]
!MINIAPP                                        kmin,                & ! [IN]
!MINIAPP                                        kmax,                & ! [IN]
!MINIAPP                                        rho   (:,:),         & ! [IN]
!MINIAPP                                        tem   (:,:),         & ! [IN]
!MINIAPP                                        ccn   (:,:),         & ! [IN]
!MINIAPP                                        GDCFRC(:,:),         & ! [IN]
!MINIAPP                                        q     (:,:,I_QC),    & ! [IN]
!MINIAPP                                        re    (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                        re    (:,:,I_mp_QI)  ) ! [OUT]
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP       elseif( MP_TYPE == 'LSC_MIROC6' ) then
!MINIAPP          if ( opt_2moment_water ) then
!MINIAPP             !$acc kernels pcopy(nc,nr,ni,ns,ng) pcopyin(rho,q)
!MINIAPP             !$omp parallel do default(none),private(ij,k), &
!MINIAPP             !$omp shared(ijdim,kdim,nc,ni,rho,q,I_NC,I_NI )
!MINIAPP             do k  = 1, kdim
!MINIAPP             do ij = 1, ijdim
!MINIAPP                nc(ij,k) = rho(ij,k) * q(ij,k,I_NC)
!MINIAPP                ni(ij,k) = rho(ij,k) * q(ij,k,I_NI)
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP             !$omp end parallel do
!MINIAPP             !$acc end kernels
!MINIAPP          endif
!MINIAPP          call mp_lsc_miroc6_effective_radius( ijdim,                & ! [IN]
!MINIAPP                                               kdim,                 & ! [IN]
!MINIAPP                                               kmin,                 & ! [IN]
!MINIAPP                                               kmax,                 & ! [IN]
!MINIAPP                                               rho    (:,:),         & ! [IN]
!MINIAPP                                               tem    (:,:),         & ! [IN]
!MINIAPP                                               ccn    (:,:),         & ! [IN]
!MINIAPP                                               GDCFRC (:,:),         & ! [IN]
!MINIAPP                                               q      (:,:,I_QC),    & ! [IN]
!MINIAPP                                               q      (:,:,I_QI),    & ! [IN]
!MINIAPP                                               nc     (:,:),         & ! [IN]
!MINIAPP                                               ni     (:,:),         & ! [IN]
!MINIAPP                                               re     (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                               re     (:,:,I_mp_QI)  ) ! [OUT]
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP       elseif( MP_TYPE == 'S10' ) then
!MINIAPP
!MINIAPP          call mp_s10_effective_radius( ijdim,               & ! [IN]
!MINIAPP                                        kdim,                & ! [IN]
!MINIAPP                                        kmin,                & ! [IN]
!MINIAPP                                        kmax,                & ! [IN]
!MINIAPP                                        rho   (:,:),         & ! [IN]
!MINIAPP                                        tem   (:,:),         & ! [IN]
!MINIAPP                                        q     (:,:,:),       & ! [IN]
!MINIAPP                                        re    (:,:,I_mp_QC), & ! [OUT]
!MINIAPP                                        re    (:,:,I_mp_QR), & ! [OUT]
!MINIAPP                                        re    (:,:,I_mp_QI), & ! [OUT]
!MINIAPP                                        re    (:,:,I_mp_QS), & ! [OUT]
!MINIAPP                                        re    (:,:,I_mp_QG)  ) ! [OUT]
!MINIAPP
!MINIAPP       elseif( MP_TYPE == 'NONE' ) then
!MINIAPP
!MINIAPP          !$acc kernels pcopy(re)
!MINIAPP          !$omp parallel do default(none),private(ij,k), &
!MINIAPP          !$omp shared(ijdim,kdim,re)
!MINIAPP          do k  = 1, kdim
!MINIAPP          do ij = 1, ijdim
!MINIAPP             re(ij,k,I_mp_QC) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QR) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QI) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QS) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QG) = 0.0_DP
!MINIAPP             re(ij,k,I_mp_QH) = 0.0_DP
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$omp end parallel do
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP       endif
!MINIAPP
!MINIAPP    else ! set default
!MINIAPP
!MINIAPP       !$acc kernels pcopy(re)
!MINIAPP       !$omp parallel do default(none),private(ij,k), &
!MINIAPP       !$omp shared(ijdim,kdim,re,const_re_cld,const_re_rain,const_re_ice, &
!MINIAPP       !$omp        const_re_snow,const_re_graupel,const_re_hail)
!MINIAPP       do k  = 1, kdim
!MINIAPP       do ij = 1, ijdim
!MINIAPP          re(ij,k,I_mp_QC) = const_re_cld
!MINIAPP          re(ij,k,I_mp_QR) = const_re_rain
!MINIAPP          re(ij,k,I_mp_QI) = const_re_ice
!MINIAPP          re(ij,k,I_mp_QS) = const_re_snow
!MINIAPP          re(ij,k,I_mp_QG) = const_re_graupel
!MINIAPP          re(ij,k,I_mp_QH) = const_re_hail
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end parallel do
!MINIAPP       !$acc end kernels
!MINIAPP
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! [Add] 2016/03/29 T.Seiki
!MINIAPP    call history_in( 'ml_re_cld'    , re(:,:,I_mp_QC) )
!MINIAPP    call history_in( 'ml_re_rain'   , re(:,:,I_mp_QR) )
!MINIAPP    call history_in( 'ml_re_ice'    , re(:,:,I_mp_QI) )
!MINIAPP    call history_in( 'ml_re_snow'   , re(:,:,I_mp_QS) )
!MINIAPP    call history_in( 'ml_re_graupel', re(:,:,I_mp_QG) )
!MINIAPP    call history_in( 'ml_re_hail'   , re(:,:,I_mp_QH) )
!MINIAPP#if !defined(SINGLE) && !defined(SINGLE_PHY) && !defined(SINGLE_MP)
!MINIAPP    ! do nothing
!MINIAPP#else
!MINIAPP    !$acc kernels
!MINIAPP    do nq = 1, Nparticle_mp
!MINIAPP    do k = 1, kdim
!MINIAPP    do ij = 1, ijdim
!MINIAPP       RP_PHY_re(ij,k,nq) = real(re(ij,k,nq), RP_PHY)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP#endif
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine mp_effective_radius

end module mod_mp_driver
