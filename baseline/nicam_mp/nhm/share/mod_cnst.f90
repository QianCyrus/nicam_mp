!-------------------------------------------------------------------------------
!> module CONSTANT
!!
!! @par Description
!!         Physical constants
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_cnst
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
  public :: CNST_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP),   public, parameter :: CNST_PI       = 3.1415926535897932_RP     ! pi
  real(RP),   public            :: CNST_D2R                                  ! degree to radian
  real(RP),   public            :: CNST_EPS      = 1.E-16_RP                 ! small number
  real(SP),   public            :: CNST_EPS_SP   = 1.E-16_SP                 ! small number (REAL4)
  real(DP),   public            :: CNST_EPS_DP   = 1.E-16_DP                 ! small number (REAL8)
  real(RP),   public            :: CNST_EPS1     = 0.99999999999999_RP       ! small number
  real(SP),   public            :: CNST_EPS1_SP  = 0.999999999_SP            ! small number (REAL4)
  real(DP),   public            :: CNST_EPS1_DP  = 0.99999999999999_DP       ! small number (REAL8)
  real(RP),   public            :: CNST_HUGE     = 1.E+30_RP                 ! huge  number
  real(SP),   public            :: CNST_HUGE_SP  = 1.E+30_SP                 ! huge  number (REAL4)
  real(DP),   public            :: CNST_HUGE_DP  = 1.E+30_DP                 ! huge  number (REAL8)
  real(RP),   public            :: CNST_TINY     = 1.E-30_RP                 ! tiny  number

  integer(2), public, parameter :: CNST_UNDEF2   = -32767                    ! undefined value (INT2)
  real(SP),   public, parameter :: CNST_UNDEF4   = -9.9999E30_SP             ! undefined value (REAL4)
  real(DP),   public, parameter :: CNST_UNDEF8   = -9.9999E30_DP             ! undefined value (REAL8)
  real(RP),   public            :: CNST_UNDEF    = -9.9999E30_RP
  !$acc declare create(CNST_UNDEF)
  real(SP),   public, parameter :: CNST_UNDEF_LL = -9.99E34_SP               ! undefined value for lat-lon data

  ! adopted constants
  real(RP),   public            :: CNST_RADIUS   = 6371.01E+3_RP             ! radius of the planet [m]
  real(RP),   public            :: CNST_OHM      = 7.292115E-5_RP            ! angular velocity of the planet [1/s]
  real(RP),   public            :: CNST_GRAV     = 9.80665_RP                ! standard acceleration of gravity [m/s2]

  ! physical constants
  real(RP),   public, parameter :: CNST_STB      = 5.670367E-8_RP            ! Stefan-Boltzman constant [W/m2/K4]
  real(RP),   public, parameter :: CNST_KARMAN   = 0.4_RP                    ! von Karman constant
  real(RP),   public, parameter :: CNST_R        = 8.3144621_RP              ! universal gas constant [J/mol/K]

  ! dry air constants
  real(RP),   public            :: CNST_Mdry     =   28.97_RP                ! mass weight (dry air)                     [g/mol]
  real(RP),   public            :: CNST_Rdry     =  287.05_RP                ! specific gas constant (dry air)           [J/kg/K]
  real(RP),   public            :: CNST_CPdry    = 1004.69_RP                ! specific heat (dry air,constant pressure) [J/kg/K]
  real(RP),   public            :: CNST_CVdry                                ! specific heat (dry air,constant volume)   [J/kg/K]
  real(RP),   public            :: CNST_LAPS     =  6.5E-3_RP                ! lapse rate of ISA                         [K/m]
  real(RP),   public            :: CNST_LAPSdry                              ! dry adiabatic lapse rate                  [K/m]

  ! water constants
  real(RP),   public            :: CNST_Mvap     =   18.02_RP                ! mass weight (water vapor)                      [g/mol]
  real(RP),   public            :: CNST_Rvap     =  461.51_RP                ! specific gas constant (water vapor)            [J/kg/K]
  !$acc declare copyin(CNST_Rvap)
  real(RP),   public            :: CNST_CPvap    = 1846.06_RP                ! specific heat (water vapor, constant pressure) [J/kg/K]
  real(RP),   public            :: CNST_CVvap                                ! specific heat (water vapor, constant volume)   [J/kg/K]
  real(RP),   public, parameter :: CNST_CL       =  4217.8_RP                ! specific heat (liquid water)                   [J/kg/K]
  real(RP),   public, parameter :: CNST_CI       =  2106.0_RP                ! specific heat (ice)                            [J/kg/K]

  real(RP),   public            :: CNST_EPSvap                               ! Rdry / Rvap
  !$acc declare create(CNST_EPSvap)
  real(RP),   public            :: CNST_EPSTvap                              ! 1 / epsilon - 1

  real(RP),   public            :: CNST_EMELT    = 3.3366E+5_RP              ! = CNST_LHF0
  !$acc declare copyin(CNST_EMELT)
  real(RP),   public, parameter :: CNST_TMELT    = 273.15_RP                 ! Freeze point of water
  !$acc declare copyin(CNST_TMELT)
  real(RP),   public            :: CNST_TFRZS    = 271.35_RP                 ! Freeze point of sea

  real(RP),   public            :: CNST_LHV                                  ! latent heat of vaporizaion for use
  real(RP),   public            :: CNST_LHS                                  ! latent heat of sublimation for use
  real(RP),   public            :: CNST_LHF                                  ! latent heat of fusion      for use
  real(RP),   public            :: CNST_LHV0     = 2.50084E+6_RP             ! latent heat of vaporizaion at 0C [J/kg]
  !$acc declare copyin(CNST_LHV0)
  real(RP),   public            :: CNST_LHV00                                ! latent heat of vaporizaion at 0K [J/kg]
  real(RP),   public            :: CNST_LHS0     = 2.8345E+6_RP              ! latent heat of sublimation at 0C [J/kg]
  real(RP),   public            :: CNST_LHS00                                ! latent heat of sublimation at 0K [J/kg]
  real(RP),   public            :: CNST_LHF0                                 ! latent heat of fusion      at 0C [J/kg]
  real(RP),   public            :: CNST_LHF00                                ! latent heat of fusion      at 0K [J/kg]
  real(RP),   public, parameter :: CNST_PSAT0    =   610.7_RP                ! saturate pressure of water vapor at 0C [Pa]
  !$acc declare copyin(CNST_PSAT0)
  real(RP),   public, parameter :: CNST_DWATR    =  1000.0_RP                ! density of water [kg/m3]
  real(RP),   public, parameter :: CNST_DICE     =   916.8_RP                ! density of ice   [kg/m3]

  ! standards
  real(RP),   public            :: CNST_SOUND                                ! speed of sound (dry air at 0C) [m/s]

  real(RP),   public            :: CNST_Pstd     = 101325.0_RP               ! standard pressure [Pa]
  real(RP),   public            :: CNST_PRE00    = 100000.0_RP               ! pressure reference [Pa]
  real(RP),   public            :: CNST_Tstd     =   288.15_RP               ! standard temperature (15C) [K]
  real(RP),   public, parameter :: CNST_TEM00    =   273.15_RP               ! temperature reference (0C) [K]
  real(RP),   public, parameter :: CNST_PPM      =    1.E-6_RP               ! parts par million

  !%%%%% Obsolete variables for backward compatibility %%%%%
  real(DP),   public, parameter :: CNST_EPS_ZERO = 1.0D-99
  real(DP),   public, parameter :: CNST_MAX_REAL = 1.0D+99

  real(RP),   public            :: CNST_ERADIUS ! => CNST_RADIUS             ! radius of the planet [m]
  real(RP),   public            :: CNST_EGRAV   ! => CNST_GRAV               ! standard acceleration of gravity [m/s2]

  real(RP),   public            :: CNST_RAIR    ! => CNST_Rdry               ! specific gas constant (dry air)           [J/kg/K]
  real(RP),   public            :: CNST_CP      ! => CNST_CPdry              ! specific heat (dry air,constant pressure) [J/kg/K]
  real(RP),   public            :: CNST_CV      ! => CNST_CVdry              ! specific heat (dry air,constant volume)   [J/kg/K]
  real(RP),   public            :: CNST_GAMMA   !                            ! CPdry / CVdry
  real(RP),   public            :: CNST_KAPPA   !                            ! Rdry  / CPdry
  real(RP),   public            :: CNST_CPV     ! => CNST_CPvap              ! specific heat (water vapor, constant pressure) [J/kg/K]
  real(RP),   public            :: CNST_CVV     ! => CNST_CVvap              ! specific heat (water vapor, constant volume)   [J/kg/K]
  real(RP),   public            :: CNST_EPSV    ! => CNST_EPSvap             ! Rdry / Rvap
  real(RP),   public            :: CNST_EPSVT   ! => CNST_EPSTvap            ! 1 / epsilon - 1
  real(RP),   public            :: CNST_LH0     ! => CNST_LHV0               ! latent heat of vaporizaion at 0C [J/kg]
  real(RP),   public            :: CNST_LH00    ! => CNST_LHV00              ! latent heat of vaporizaion at 0K [J/kg]

  real(RP),   public            :: CNST_PRES0   ! => CNST_Pstd               ! standard pressure [Pa]
  real(RP),   public, parameter :: CNST_TQICE     =  CNST_TMELT              ! Freeze point of water
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  subroutine CNST_setup
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    real(RP) :: earth_radius                 ! Earth radius
    real(RP) :: earth_angvel                 ! Anguler velocity of the earth
    real(RP) :: small_planet_factor = 1.0_RP ! small planet factor
    real(RP) :: earth_gravity                ! Gravitational accelaration
    real(RP) :: gas_cnst                     ! Gas constant of dry air
    real(RP) :: gas_cnst_vap                 ! Gas constant of water vapour
    real(RP) :: specific_heat_pre            ! Specific heat of air( const pre )
    real(RP) :: specific_heat_pre_vap        ! Specific heat of water vapour ( const pre )
    real(RP) :: latent_heat_vap              ! latent heat of vaporization LH0 ( 0 deg )
    real(RP) :: latent_heat_sub              ! latent heat of sublimation LHS0 ( 0 deg )

    namelist / CNSTPARAM / &
       earth_radius,          &
       earth_angvel,          &
       small_planet_factor,   &
       earth_gravity,         &
       gas_cnst,              &
       gas_cnst_vap,          &
       specific_heat_pre,     &
       specific_heat_pre_vap, &
       latent_heat_vap,       &
       latent_heat_sub

    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- initialization of controled parameters
    earth_gravity         = CNST_GRAV
    earth_radius          = CNST_RADIUS
    earth_angvel          = CNST_OHM
    gas_cnst              = CNST_Rdry
    gas_cnst_vap          = CNST_Rvap
    specific_heat_pre     = CNST_CPdry
    specific_heat_pre_vap = CNST_CPvap
    latent_heat_vap       = CNST_LHV0
    latent_heat_sub       = CNST_LHS0

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[cnst]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=CNSTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** CNSTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist CNSTPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist CNSTPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=CNSTPARAM)

    CNST_GRAV   = earth_gravity
    CNST_RADIUS = earth_radius / small_planet_factor
    CNST_OHM    = earth_angvel * small_planet_factor
    CNST_Rdry   = gas_cnst
    CNST_Rvap   = gas_cnst_vap
    CNST_CPdry  = specific_heat_pre
    CNST_CPvap  = specific_heat_pre_vap
    CNST_LHV0   = latent_heat_vap
    CNST_LHS0   = latent_heat_sub

    if    ( RP == SP ) then
       CNST_UNDEF = real(CNST_UNDEF4,kind=RP)
    elseif( RP == DP ) then
       CNST_UNDEF = real(CNST_UNDEF8,kind=RP)
    else
       write(*,*) 'xxx unsupported precision: ', RP
       call ADM_proc_stop
    endif

    !CNST_PI      = 4.E0_RP * atan( 1.0_RP )
    CNST_D2R     = CNST_PI / 180.0_RP
    CNST_EPS     =          epsilon(0.0_RP)
    CNST_EPS_SP  =          epsilon(0.0_SP)
    CNST_EPS_DP  =          epsilon(0.0_DP)
    CNST_EPS1    = 1.0_RP - epsilon(0.0_RP)
    CNST_EPS1_SP = 1.0_SP - epsilon(0.0_SP)
    CNST_EPS1_DP = 1.0_DP - epsilon(0.0_DP)
    CNST_HUGE    =             huge(0.0_RP) * 1.E-4_RP
    CNST_HUGE_SP =             huge(0.0_SP) * 1.E-4_SP
    CNST_HUGE_DP =             huge(0.0_DP) * 1.E-4_DP
    CNST_TINY    =             tiny(0.0_RP)

    CNST_CVdry   = CNST_CPdry - CNST_Rdry
    CNST_LAPSdry = CNST_GRAV / CNST_CPdry

    CNST_CVvap   = CNST_CPvap - CNST_Rvap
    CNST_EPSvap  = CNST_Rdry  / CNST_Rvap
    CNST_EPSTvap = 1.0_RP     / CNST_EPSvap - 1.0_RP

    CNST_LHF0    = CNST_LHS0 - CNST_LHV0
    CNST_EMELT   = CNST_LHF0

    CNST_LHV00   = CNST_LHV0 - ( CNST_CPvap - CNST_CL ) * CNST_TEM00
    CNST_LHS00   = CNST_LHS0 - ( CNST_CPvap - CNST_CI ) * CNST_TEM00
    CNST_LHF00   = CNST_LHF0 - ( CNST_CL    - CNST_CI ) * CNST_TEM00

    CNST_SOUND = sqrt( CNST_CPdry * CNST_Rdry / ( CNST_CPdry - CNST_Rdry ) * CNST_TEM00 )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*)           '====== Precision ======'
    write(IO_FID_LOG,'(1x,A,I3)') '--- kind      (floating point value) = ', kind     (1.0_RP)
    write(IO_FID_LOG,'(1x,A,I3)') '--- precision (floating point value) = ', precision(1.0_RP)
    write(IO_FID_LOG,'(1x,A,I3)') '--- range     (floating point value) = ', range    (1.0_RP)

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '====== List of constants ======'
    write(IO_FID_LOG,*) '--- PI                                                : PI      = ', CNST_PI
    write(IO_FID_LOG,*) '--- Small number                                      : EPS     = ', CNST_EPS
    write(IO_FID_LOG,*) '--- Small number (1-EPS)                              : EPS1    = ', CNST_EPS1
    write(IO_FID_LOG,*) '--- Huge  number                                      : HUGE    = ', CNST_HUGE
    write(IO_FID_LOG,*) '--- Tiny  number                                      : TINY    = ', CNST_TINY
    write(IO_FID_LOG,*) '--- undefined number (INT2)                           : UNDEF2  = ', CNST_UNDEF2
    write(IO_FID_LOG,*) '--- undefined number (REAL,general use)               : UNDEF   = ', CNST_UNDEF
    write(IO_FID_LOG,*) '--- undefined number (REAL4)                          : UNDEF4  = ', CNST_UNDEF4
    write(IO_FID_LOG,*) '--- undefined number (REAL8)                          : UNDEF8  = ', CNST_UNDEF8

    write(IO_FID_LOG,*) '--- radius of the planet                          [m] : RADIUS  = ', CNST_RADIUS
    write(IO_FID_LOG,*) '--- angular velocity of the planet              [1/s] : OHM     = ', CNST_OHM
    write(IO_FID_LOG,*) '--- standard acceleration of gravity           [m/s2] : GRAV    = ', CNST_GRAV

    write(IO_FID_LOG,*) '--- Stefan-Boltzman constant                [W/m2/K4] : STB     = ', CNST_STB
    write(IO_FID_LOG,*) '--- von Karman constant                               : KARMAN  = ', CNST_KARMAN
    write(IO_FID_LOG,*) '--- universal gas constant                  [J/mol/K] : R       = ', CNST_R

    write(IO_FID_LOG,*) '--- mass weight (dry air)                     [g/mol] : Mdry    = ', CNST_Mdry
    write(IO_FID_LOG,*) '--- specific gas constant (dry air)          [J/kg/K] : Rdry    = ', CNST_Rdry
    write(IO_FID_LOG,*) '--- specific heat (dry air, const. pressure) [J/kg/K] : CPdry   = ', CNST_CPdry
    write(IO_FID_LOG,*) '--- specific heat (dry air, const. volume)   [J/kg/K] : Cvdry   = ', CNST_CVdry
    write(IO_FID_LOG,*) '--- lapse rate of ISA                           [K/m] : LAPS    = ', CNST_LAPS
    write(IO_FID_LOG,*) '--- dry adiabatic lapse rate                    [K/m] : LAPSdry = ', CNST_LAPSdry

    write(IO_FID_LOG,*) '--- mass weight (water vapor)                 [g/mol] : Mvap    = ', CNST_Mvap
    write(IO_FID_LOG,*) '--- specific gas constant (water vapor)      [J/kg/K] : Rvap    = ', CNST_Rvap
    write(IO_FID_LOG,*) '--- specific heat (vapor, const. pressure)   [J/kg/K] : CPvap   = ', CNST_CPvap
    write(IO_FID_LOG,*) '--- specific heat (vapor, const. volume)     [J/kg/K] : CVvap   = ', CNST_CVvap
    write(IO_FID_LOG,*) '--- specific heat (liquid water)             [J/kg/K] : CL      = ', CNST_CL
    write(IO_FID_LOG,*) '--- specific heat (ice)                      [J/kg/K] : CI      = ', CNST_CI
    write(IO_FID_LOG,*) '--- Rdry / Rvap                                       : EPSvap  = ', CNST_EPSvap
    write(IO_FID_LOG,*) '--- 1 / EPSvap - 1                                    : EPSTvap = ', CNST_EPSTvap

    write(IO_FID_LOG,*) '--- latent heat of vaporizaion at 0C           [J/kg] : LHV0    = ', CNST_LHV0
    write(IO_FID_LOG,*) '--- latent heat of sublimation at 0C           [J/kg] : LHS0    = ', CNST_LHS0
    write(IO_FID_LOG,*) '--- latent heat of fusion      at 0C           [J/kg] : LHF0    = ', CNST_LHF0
    write(IO_FID_LOG,*) '--- latent heat of vaporizaion at 0K           [J/kg] : LHV00   = ', CNST_LHV00
    write(IO_FID_LOG,*) '--- latent heat of sublimation at 0K           [J/kg] : LHS00   = ', CNST_LHS00
    write(IO_FID_LOG,*) '--- latent heat of fusion      at 0K           [J/kg] : LHF00   = ', CNST_LHF00
    write(IO_FID_LOG,*) '--- saturate pressure of water vapor at 0C       [Pa] : PSAT0   = ', CNST_PSAT0
    write(IO_FID_LOG,*) '--- density of water                          [kg/m3] : DWATR   = ', CNST_DWATR
    write(IO_FID_LOG,*) '--- density of ice                            [kg/m3] : DICE    = ', CNST_DICE

    write(IO_FID_LOG,*) '--- speed of sound (dry air at 0C)              [m/s] : SOUND   = ', CNST_SOUND
    write(IO_FID_LOG,*) '--- standard pressure                            [Pa] : Pstd    = ', CNST_Pstd
    write(IO_FID_LOG,*) '--- pressure reference                           [Pa] : PRE00   = ', CNST_PRE00
    write(IO_FID_LOG,*) '--- standard temperature (15C)                    [K] : Tstd    = ', CNST_Tstd
    write(IO_FID_LOG,*) '--- temperature reference (0C)                    [K] : TEM00   = ', CNST_TEM00


    !%%%%% Update variable on GPU for acc routine %%%%%
    !$acc update device(CNST_EPSvap)
    !$acc update device(CNST_UNDEF)

    !%%%%% Obsolete variables for backward compatibility %%%%%
    CNST_ERADIUS = CNST_RADIUS
    CNST_EGRAV   = CNST_GRAV

    CNST_RAIR  = CNST_Rdry
    CNST_CP    = CNST_CPdry
    CNST_CV    = CNST_CVdry
    CNST_GAMMA = CNST_CPdry / CNST_CVdry
    CNST_KAPPA = CNST_Rdry / CNST_CPdry
    CNST_CPV   = CNST_CPvap
    CNST_CVV   = CNST_CVvap
    CNST_EPSV  = CNST_EPSvap
    CNST_EPSVT = CNST_EPSTvap
    CNST_LH0   = CNST_LHV0
    CNST_LH00  = CNST_LHV00

    CNST_PRES0 = CNST_Pstd
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return
  end subroutine CNST_setup

end module mod_cnst
