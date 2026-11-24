!-------------------------------------------------------------------------------
!
!+  Water isotope module
!
!-------------------------------------------------------------------------------
module mod_isotope
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module defines isotope parameters for water.
  !
  !++ Current Corresponding Author : K.Yoshimura
  !
  !++ History:
  !      Version    Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00       11-04-14  [new]
  !      0.01       20-08-28  [Add] WISO_setup
  !                           [Del] OPT_NOFRAC
  !      -----------------------------------------------------------------------
  !
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
  public :: WISO_setup
  public :: WISO_fract
  public :: WISO_eqfract
  public :: WISO_eqtime

  interface WISO_fract
     module procedure WISO_fract_SP
     module procedure WISO_fract_DP
  end interface WISO_fract

  interface WISO_eqfract
     module procedure WISO_eqfract_SP
     module procedure WISO_eqfract_DP
  end interface WISO_eqfract

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public, parameter :: NISO = 2 ! number of isotope spicies

  real(RP), public :: DIFF_ISO1  = 1.0_RP
  real(RP), public :: DIFF_ISO2  = 1.0_RP
  real(RP), public :: DIFF_H18O  = 1.02849_RP ! diffusion coefficient for H18O
  real(RP), public :: DIFF_HDO   = 1.02512_RP ! diffusion coefficient for HDO
  real(RP), public :: EQFRAC_CP  = 0.45_RP    ! isotope exchange ratio for CP (AS,CHIKIRA)
  real(RP), public :: EQFRAC_MP  = 0.95_RP    ! isotope exchange ratio for MP (LSC)

  logical,  public :: OPT_NOKINETIC         = .false.
  logical,  public :: OPT_KINETIC_SICE_ORIG = .true.  ! Jouzel & Merlivat (1984)
  logical,  public :: OPT_EQTIME            = .false. ! Calculate isotope exchange ratio (LSC, CHIKIRA)
  real(RP), public :: ALP0(2,2,3)                     ! coefficients for equilibrium fractionation factor
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: DEGREE_FREE      = 0.58_RP
  real(RP), private :: JM84_SICE_a      = 1.0_RP
  real(RP), private :: JM84_SICE_b      = 0.002_RP

  logical,  private :: OPT_FRACT_LAMBER = .false. ! Use ALPI for HDO from Lamber et al. (2017)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine WISO_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_runconf, only: &
       opt_disable_wisofrac
    implicit none

    namelist / PARAM_WISO / &
       DIFF_H18O,             &
       DIFF_HDO,              &
       DEGREE_FREE,           &
       JM84_SICE_a,           &
       JM84_SICE_b,           &
       EQFRAC_CP,             &
       EQFRAC_MP,             &
       OPT_KINETIC_SICE_ORIG, &
       OPT_NOKINETIC,         &
       OPT_EQTIME,            &
       OPT_FRACT_LAMBER

    integer :: dummy1
    integer :: dummy2

    namelist / FRACT_PARAM / &
       dummy1

    namelist / FRACT_OPT / &
       dummy2

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[isotope]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_WISO,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** PARAM_WISO is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist PARAM_WISO. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist PARAM_WISO. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=PARAM_WISO)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=FRACT_PARAM,iostat=ierr)
    if( ierr > 0 ) then
       write(*         ,*) 'xxx FRACT_PARAM is obsolete parameter. Please use PARAM_WISO. STOP.'
       write(IO_FID_LOG,*) 'xxx FRACT_PARAM is obsolete parameter. Please use PARAM_WISO. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=FRACT_PARAM)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=FRACT_OPT,iostat=ierr)
    if( ierr > 0 ) then
       write(*         ,*) 'xxx FRACT_OPT is obsolete parameter. Please use PARAM_WISO. STOP.'
       write(IO_FID_LOG,*) 'xxx FRACT_OPT is obsolete parameter. Please use PARAM_WISO. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=FRACT_OPT)

    if ( opt_disable_wisofrac ) then
       DIFF_H18O   = 1.0_RP
       DIFF_HDO    = 1.0_RP
       DEGREE_FREE = 1.0_RP
       EQFRAC_CP   = 1.0_RP
       EQFRAC_MP   = 1.0_RP
       ALP0(:,:,:) = 0.0_RP
    else
       ! -- from Majoube (1970,1971)
       ALP0(1,1,:) = (/ -0.002067_RP, -0.4156_RP,  1137.0_RP /) ! ALP (1) = EXP( -0.002067 - 0.4156/TGTX + 1137./(TGTX**2) )
       ALP0(2,1,:) = (/ -0.028224_RP,  11.839_RP,     0.0_RP /) ! ALPI(1) = EXP( -0.028224 + 11.839/TGTX )

       ! -- from Mejoube (1971)
       ALP0(1,2,:) = (/  0.052612_RP, -76.248_RP, 24844.0_RP /) ! ALP (2) = EXP( 0.052612 - 76.248/TGTX + 24844./(TGTX**2))
      !ALP0(2,2,:) = (/ -0.0934_RP,       0.0_RP, 16288.0_RP /) ! ALPI(2) = EXP( -0.0934 + 16288./(TGTX**2) )
       ALP0(2,2,:) = (/ -0.0945_RP,       0.0_RP, 16289.0_RP /)

       if ( OPT_FRACT_LAMBER ) then
          ! -- from Lamber et al. (2017) PNAS
          ALP0(2,2,:) = (/ -0.0559_RP, 0.0_RP, 13525.0_RP /) ! ALPI(2) = EXP( -0.0559 + 13525./(TGTX**2) )
       endif
    endif

    DIFF_ISO1 = DIFF_H18O**DEGREE_FREE
    DIFF_ISO2 = DIFF_HDO **DEGREE_FREE

    return
  end subroutine WISO_setup

  !-----------------------------------------------------------------------------
  subroutine WISO_fract_SP( &
       ALP,              &
       ALPI,             &
       ALPS,             &
       TGTX,             &
       ALP0,             &
       kinetic_off,      &
       kinetic_sice_orig )
    use mod_cnst, only: &
       CNST_TMELT, &
       CNST_EMELT, &
       CNST_Rvap
    implicit none

    real(SP), intent(out) :: ALP (NISO)
    real(SP), intent(out) :: ALPI(NISO)
    real(SP), intent(out) :: ALPS(NISO)
    real(SP), intent(in)  :: TGTX
    real(RP), intent(in)  :: ALP0(2,NISO,3) ! coefficients for equilibrium fractionation factor
    logical,  intent(in)  :: kinetic_off
    logical,  intent(in)  :: kinetic_sice_orig

    real(SP), parameter :: TSICE = 273.15_SP
    real(SP), parameter :: TWICE = 253.15_SP
    real(SP) :: SICE
    real(SP) :: SICE0, SICE1, Z
    !---------------------------------------------------------------------------

    call WISO_eqfract( TGTX, ALP0, & ! [IN]
                       ALP, ALPI   ) ! [OUT]

    if ( kinetic_off ) then
       ALPS(1) = ALPI(1)
       ALPS(2) = ALPI(2)

       !AEVP(1) = ALP (1)
       !AEVP(2) = ALP (2)
    else
       if ( kinetic_sice_orig ) then
          SICE = JM84_SICE_a - JM84_SICE_b * min( TGTX - CNST_TMELT, 0.0_RP )
       else
          if    ( TGTX >= TSICE ) then
             SICE = 1.0_SP
          elseif( TGTX <= TWICE ) then
             SICE = 1.02_SP - 0.0038_SP * ( TGTX - CNST_TMELT )
          else
             Z     = (TGTX-TWICE) / (TSICE-TWICE)
             SICE0 = exp( -CNST_EMELT/CNST_Rvap * ( 1.0_SP/CNST_TMELT - 1.0_SP/TGTX ) )
             SICE1 = 1.02_SP - 0.0038_SP * ( TGTX - CNST_TMELT )
             SICE  = ( SICE1-SICE0 ) / log(10.0_SP+1.0_SP) * log(10.0_SP*(1.0_SP-Z)+1.0_SP) + SICE0
          endif
       endif

       !-- from Jouzel and Merlivat (1984)  1/0.9723, 1/0.9755
       ALPS(1) = ALPI(1) * SICE / ( 1.0_SP + ALPI(1) * (SICE-1.0_SP) * DIFF_H18O )
       ALPS(2) = ALPI(2) * SICE / ( 1.0_SP + ALPI(2) * (SICE-1.0_SP) * DIFF_HDO  )

       !-- from Stewart (1975) (1/0.9723)**(0.58) (1/0.9755)**(0.58)
       !AEVP(1) = 1.0164D0 * ALP(1)
       !AEVP(2) = 1.0144D0 * ALP(2)
    endif

    return
  end subroutine WISO_fract_SP

  !-----------------------------------------------------------------------------
  subroutine WISO_fract_DP( &
       ALP,              &
       ALPI,             &
       ALPS,             &
       TGTX,             &
       ALP0,             &
       kinetic_off,      &
       kinetic_sice_orig )
    use mod_cnst, only: &
       CNST_TMELT, &
       CNST_EMELT, &
       CNST_Rvap
    implicit none

    real(DP), intent(out) :: ALP (NISO)
    real(DP), intent(out) :: ALPI(NISO)
    real(DP), intent(out) :: ALPS(NISO)
    real(DP), intent(in)  :: TGTX
    real(RP), intent(in)  :: ALP0(2,NISO,3) ! coefficients for equilibrium fractionation factor
    logical,  intent(in)  :: kinetic_off
    logical,  intent(in)  :: kinetic_sice_orig

    real(DP), parameter :: TSICE = 273.15_DP
    real(DP), parameter :: TWICE = 253.15_DP
    real(DP) :: SICE
    real(DP) :: SICE0, SICE1, Z
    !---------------------------------------------------------------------------

    call WISO_eqfract( TGTX, ALP0, & ! [IN]
                       ALP, ALPI   ) ! [OUT]

    if ( kinetic_off ) then
       ALPS(1) = ALPI(1)
       ALPS(2) = ALPI(2)

       !AEVP(1) = ALP (1)
       !AEVP(2) = ALP (2)
    else
       if ( kinetic_sice_orig ) then
          SICE = JM84_SICE_a - JM84_SICE_b * min( TGTX - CNST_TMELT, 0.0_RP )
       else
          if    ( TGTX >= TSICE ) then
             SICE = 1.0_DP
          elseif( TGTX <= TWICE ) then
             SICE = 1.02_DP - 0.0038_DP * ( TGTX - CNST_TMELT )
          else
             Z     = (TGTX-TWICE) / (TSICE-TWICE)
             SICE0 = exp( -CNST_EMELT/CNST_Rvap * ( 1.0_DP/CNST_TMELT - 1.0_DP/TGTX ) )
             SICE1 = 1.02_DP - 0.0038_DP * ( TGTX - CNST_TMELT )
             SICE  = ( SICE1-SICE0 ) / log(10.0_DP+1.0_DP) * log(10.0_DP*(1.0_DP-Z)+1.0_DP) + SICE0
          endif
       endif

       !-- from Jouzel and Merlivat (1984)  1/0.9723, 1/0.9755
       ALPS(1) = ALPI(1) * SICE / ( 1.0_DP + ALPI(1) * (SICE-1.0_DP) * DIFF_H18O )
       ALPS(2) = ALPI(2) * SICE / ( 1.0_DP + ALPI(2) * (SICE-1.0_DP) * DIFF_HDO  )

       !-- from Stewart (1975) (1/0.9723)**(0.58) (1/0.9755)**(0.58)
       !AEVP(1) = 1.0164D0 * ALP(1)
       !AEVP(2) = 1.0144D0 * ALP(2)
    endif

    return
  end subroutine WISO_fract_DP

  !-----------------------------------------------------------------------------
  subroutine WISO_eqfract_SP( &
       TGTX, & ! [IN]
       ALP0, & ! [IN]
       ALP,  & ! [OUT]
       ALPI  ) ! [OUT]
    implicit none

    real(SP), intent(in)  :: TGTX           ! air temperature [K]
    real(RP), intent(in)  :: ALP0(2,NISO,3) ! coefficients for equilibrium fractionation factor
    real(SP), intent(out) :: ALP (NISO)     ! equilibrium fractionation factor for vapor-liquid phase [-]
    real(SP), intent(out) :: ALPI(NISO)     ! equilibrium fractionation factor for vapor-solid phase  [-]

    real(SP) :: TGTX2
    !---------------------------------------------------------------------------

    TGTX2 = TGTX * TGTX

    ALP (1) = exp( ALP0(1,1,1) + ALP0(1,1,2)/TGTX + ALP0(1,1,3)/TGTX2 )
    ALPI(1) = exp( ALP0(2,1,1) + ALP0(2,1,2)/TGTX + ALP0(2,1,3)/TGTX2 )
    ALP (2) = exp( ALP0(1,2,1) + ALP0(1,2,2)/TGTX + ALP0(1,2,3)/TGTX2 )
    ALPI(2) = exp( ALP0(2,2,1) + ALP0(2,2,2)/TGTX + ALP0(2,2,3)/TGTX2 )

    return
  end subroutine WISO_eqfract_SP

  !-----------------------------------------------------------------------------
  subroutine WISO_eqfract_DP( &
       TGTX, & ! [IN]
       ALP0, & ! [IN]
       ALP,  & ! [OUT]
       ALPI  ) ! [OUT]
    implicit none

    real(DP), intent(in)  :: TGTX           ! air temperature [K]
    real(RP), intent(in)  :: ALP0(2,NISO,3) ! coefficients for equilibrium fractionation factor
    real(DP), intent(out) :: ALP (NISO)     ! equilibrium fractionation factor for vapor-liquid phase [-]
    real(DP), intent(out) :: ALPI(NISO)     ! equilibrium fractionation factor for vapor-solid phase  [-]

    real(DP) :: TGTX2
    !---------------------------------------------------------------------------

    TGTX2 = TGTX * TGTX

    ALP (1) = exp( ALP0(1,1,1) + ALP0(1,1,2)/TGTX + ALP0(1,1,3)/TGTX2 )
    ALPI(1) = exp( ALP0(2,1,1) + ALP0(2,1,2)/TGTX + ALP0(2,1,3)/TGTX2 )
    ALP (2) = exp( ALP0(1,2,1) + ALP0(1,2,2)/TGTX + ALP0(1,2,3)/TGTX2 )
    ALPI(2) = exp( ALP0(2,2,1) + ALP0(2,2,2)/TGTX + ALP0(2,2,3)/TGTX2 )

    return
  end subroutine WISO_eqfract_DP

  !-----------------------------------------------------------------------------
  ! Isotope exchange ratio between falling raindrop and ambient vapor based on
  ! precipitaiton intensity. (Nasbaumer et al., 2017, JAMES)
  subroutine WISO_eqtime( &
       fe,    &
       T,     &
       pres,  &
       prec,  &
       rho,   &
       Vrain, &
       ALP,   &
       dZ     )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    real(RP), intent(out) :: fe(NISO)  ! equilibration ratio [-]
    real(RP), intent(in)  :: T         ! temperature [K]
    real(RP), intent(in)  :: pres      ! pressure [Pa]
    real(RP), intent(in)  :: prec      ! rainfall rate [kg/m2/s]
    real(RP), intent(in)  :: rho       ! density
    real(RP), intent(in)  :: Vrain     ! rain terminal velocity
    real(RP), intent(in)  :: ALP(NISO) ! equilibration factor [-]
    real(RP), intent(in)  :: dZ        ! lengh of vertical layer [m]

    real(RP) :: u, Dia(NISO), tauf, b, r, Re, Sc(NISO), x(NISO), fv(NISO), es, taue(NISO)
    !---------------------------------------------------------------------------

    ! viscosity of air [kg/m/s]
    u = ( 1.72_RP * 10.0_RP**(-5.0_RP) ) * (T/273.0_RP)**(1.5_RP) * (393.0_RP/(T+120.0_RP))

    ! diffusion of isotopic water vapor through air
    Dia(1) = ( 2.11_RP * 10.0_RP**(-5.0_RP) ) * DIFF_ISO1 * (T/273.15_RP)**(1.94_RP) * (101325.0_RP/pres)
    Dia(2) = ( 2.11_RP * 10.0_RP**(-5.0_RP) ) * DIFF_ISO2 * (T/273.15_RP)**(1.94_RP) * (101325.0_RP/pres)

    ! time scale for drop to fall through vertical layer
    tauf = dZ / Vrain

    ! raindrop radius (Williams and Gage, 2009)
    b = prec * 3600.0_RP ! conversion to [mm/h]
    r = 0.002_RP / ( 4.1_RP * b **(-0.21_RP) )

    ! Reynold number
    Re = 2.0_RP * r * rho * Vrain / u

    ! Schmidt number
    Sc(1) = u / rho / Dia(1)
    Sc(2) = u / rho / Dia(2)

    ! Ventilation factor
    x(1) = Re**(1.0_RP/2.0_RP) * Sc(1)**(1.0_RP/3.0_RP)
    x(2) = Re**(1.0_RP/2.0_RP) * Sc(2)**(1.0_RP/3.0_RP)

    if ( x(1)  >=  1.4_RP ) then
       fv(1) = 0.78_RP + 0.308 * x(1)
    else
       fv(1) = 1.0_RP + 0.108 * x(1) * x(1)
    endif

    if ( x(2)  >=  1.4_RP ) then
       fv(2) = 0.78_RP + 0.308 * x(2)
    else
       fv(2) = 1.0_RP + 0.108 * x(2) * x(2)
    endif

    ! isotopologue-specific e-folding equilibration time (unknow R and es)
    es = 6.1078_RP * 10**( 7.5 * T / (T + 237.3_RP) )
    taue(1) = ALP(1) * r*r * rho * CNST_Rvap * T / (3.0_RP * fv(1) * Dia(1) * es)
    taue(2) = ALP(2) * r*r * rho * CNST_Rvap * T / (3.0_RP * fv(2) * Dia(2) * es)

    ! equilibration ratio
    fe(1) = max( min( 1.0_RP - exp(-tauf/taue(1)), 1.0_RP ), 0.0_RP )
    fe(2) = max( min( 1.0_RP - exp(-tauf/taue(2)), 1.0_RP ), 0.0_RP )

    return
  end subroutine WISO_eqtime

end module mod_isotope
