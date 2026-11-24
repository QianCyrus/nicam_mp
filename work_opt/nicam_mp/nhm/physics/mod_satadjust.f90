!-------------------------------------------------------------------------------
!> Module saturation process
!!
!! @par Description
!!         This module is for saturation processes
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_satadjust
  !-----------------------------------------------------------------------------
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17  Imported from Satoh's nonhydro model.
  !      0.03      04-05-31  Add the ice adjustment routines.
  !                          e.g. sub[SATURATION_psat_ice].
  !                05-11-10  M.Satoh: use nqmax, initialize q
  !                07-01-19  K.Suzuki: bug fix for number of saturated point
  !                07-01-26  H.Tomita :
  !                          Rename moist_psat to SATURATION_psat_liq.
  !                          Rename moist_qsat to SATURATION_qsat_liq.
  !                          Reconstracut the saturation adjustment
  !                          including ice phase.
  !                          Add sub[SATURATION_rh] and
  !                              sub[SATURATION_setrange].
  !                07-06-05  M.Satoh : bug fix
  !                            I_QI was used in SATURATION_adjustment
  !                            regardless of RAIN_TYPE
  !                08-04-12  T.Mitsui: add subroutine omit qi in adjustment
  !                                    [SATURATION_adjustment_water]
  !                                    [moist_tsat2_water]
  !                08-04-14  H.Tomita: add option RT in [SATURATION_setrange].
  !                                    for adjustment of cloud water only.
  !                                    This is the same thing with
  !                                    [SATURATION_adjustment_water]?
  !                08-06-13  T.Mitsui: fix SATURATION_enthalpy_nolist,
  !                                    add SATURATION_dewtem(dew point temperature)
  !                09-08-18  T.Mitsui: add derivative of qsat_water and _ice
  !                                    moist_dqsw_dtem_rho and moist_dqsw_dtem_dpre
  !                                    moist_dqsi_dtem_rho and moist_dqsi_dtem_dpre
  !                10-05-11  M.Satoh: add intent(in)
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
#ifdef _OPENACC
  use openacc
#endif
  use mod_precision
  use mod_stdio
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: SATURATION_setup
  public :: SATURATION_setrange

  interface SATURATION_setrange
     module procedure SATURATION_setrange_SP
     module procedure SATURATION_setrange_DP
  end interface SATURATION_setrange

  public :: SATURATION_alpha
  public :: SATURATION_dalphadT

  public :: SATURATION_psat_all
  public :: SATURATION_psat_liq
  public :: SATURATION_psat_ice

  public :: SATURATION_qsat_liq
  public :: SATURATION_qsat_ice
  public :: SATURATION_qsat_liq2
  public :: SATURATION_qsat_ice2

  interface SATURATION_alpha
     module procedure SATURATION_alpha_0D
     module procedure SATURATION_alpha_1D
     module procedure SATURATION_alpha_2D
  end interface SATURATION_alpha

  interface SATURATION_dalphadT
     module procedure SATURATION_dalphadT_0D
     module procedure SATURATION_dalphadT_1D
     module procedure SATURATION_dalphadT_2D
  end interface SATURATION_dalphadT

  interface SATURATION_psat_all
     module procedure SATURATION_psat_all_0D
     module procedure SATURATION_psat_all_1D
     module procedure SATURATION_psat_all_2D
  end interface SATURATION_psat_all

  interface SATURATION_psat_liq
     module procedure SATURATION_psat_liq_0D
     module procedure SATURATION_psat_liq_1D_SP
     module procedure SATURATION_psat_liq_1D_DP
     module procedure SATURATION_psat_liq_2D_SP
     module procedure SATURATION_psat_liq_2D_DP
  end interface SATURATION_psat_liq

  interface SATURATION_psat_ice
     module procedure SATURATION_psat_ice_0D
     module procedure SATURATION_psat_ice_1D_SP
     module procedure SATURATION_psat_ice_1D_DP
     module procedure SATURATION_psat_ice_2D_SP
     module procedure SATURATION_psat_ice_2D_DP
  end interface SATURATION_psat_ice

  interface SATURATION_qsat_liq
     module procedure SATURATION_qsat_liq_0D
     module procedure SATURATION_qsat_liq_1D_SP
     module procedure SATURATION_qsat_liq_1D_DP
     module procedure SATURATION_qsat_liq_2D_SP
     module procedure SATURATION_qsat_liq_2D_DP
  end interface SATURATION_qsat_liq

  interface SATURATION_qsat_liq2
     module procedure SATURATION_qsat_liq2_0D
     module procedure SATURATION_qsat_liq2_1D_SP
     module procedure SATURATION_qsat_liq2_1D_DP
     module procedure SATURATION_qsat_liq2_2D_SP
     module procedure SATURATION_qsat_liq2_2D_DP
  end interface SATURATION_qsat_liq2

  interface SATURATION_qsat_ice
     module procedure SATURATION_qsat_ice_0D
     module procedure SATURATION_qsat_ice_1D_SP
     module procedure SATURATION_qsat_ice_1D_DP
     module procedure SATURATION_qsat_ice_2D_SP
     module procedure SATURATION_qsat_ice_2D_DP
  end interface SATURATION_qsat_ice

  interface SATURATION_qsat_ice2
     module procedure SATURATION_qsat_ice2_0D
     module procedure SATURATION_qsat_ice2_1D_SP
     module procedure SATURATION_qsat_ice2_1D_DP
     module procedure SATURATION_qsat_ice2_2D_SP
     module procedure SATURATION_qsat_ice2_2D_DP
  end interface SATURATION_qsat_ice2

  public :: SATURATION_rh

  interface SATURATION_rh
     module procedure SATURATION_rh_SP
     module procedure SATURATION_rh_DP
  end interface SATURATION_rh

  public :: SATURATION_dewtem

  public :: SATURATION_adjustment

  interface SATURATION_adjustment
     module procedure SATURATION_adjustment_SP
     module procedure SATURATION_adjustment_DP
  end interface SATURATION_adjustment

  public :: SATURATION_enthalpy

  public :: moist_dqsw_dtem_rho
  public :: moist_dqsi_dtem_rho
  public :: moist_dqsw_dtem_dpre
  public :: moist_dqsi_dtem_dpre


  interface moist_dqsw_dtem_rho
     module procedure moist_dqsw_dtem_rho_SP
     module procedure moist_dqsw_dtem_rho_DP
  end interface moist_dqsw_dtem_rho

  interface moist_dqsi_dtem_rho
     module procedure moist_dqsi_dtem_rho_SP
     module procedure moist_dqsi_dtem_rho_DP
  end interface moist_dqsi_dtem_rho

  interface moist_dqsw_dtem_dpre
     module procedure moist_dqsw_dtem_dpre_SP
     module procedure moist_dqsw_dtem_dpre_DP
  end interface moist_dqsw_dtem_dpre

  interface moist_dqsi_dtem_dpre
     module procedure moist_dqsi_dtem_dpre_SP
     module procedure moist_dqsi_dtem_dpre_DP
  end interface moist_dqsi_dtem_dpre

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: CPovR_liq
  real(RP), public :: CPovR_ice
  real(RP), public :: CVovR_liq
  real(RP), public :: CVovR_ice
  real(RP), public :: LovR_liq
  real(RP), public :: LovR_ice

  logical, public :: acc_wait_on = .true.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: satadjust_all_ACC
  private :: satadjust_liq_ACC
  private :: satadjust_all
  private :: satadjust_liq

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: TEM_MIN   = 10.0_RP   !< minimum temperature [K]
  real(RP), private, parameter :: DTEM_EPS0 = 1.0E-8_RP ! temperature convergence criterion

  character(len=H_SHORT), private :: ALPHA_TYPE             = 'LINEAR'

  real(RP),               private :: SATURATION_ULIMIT_TEMP = 273.15_RP     !< upper limit temperature
  real(RP),               private :: SATURATION_LLIMIT_TEMP = 233.15_RP     !< lower limit temperature

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine SATURATION_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_cnst, only: &
       CNST_Rvap,  &
       CNST_CPvap, &
       CNST_CVvap, &
       CNST_CL,    &
       CNST_CI
    use mod_runconf, only: &
       EIN_TYPE, &
       LHV,      &
       LHS
    implicit none

    NAMELIST / SATURATIONPARAM / &
       ALPHA_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[saturation]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=SATURATIONPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** SATURATIONPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist SATURATIONPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist SATURATIONPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=SATURATIONPARAM)

    if ( EIN_TYPE == 'EXACT' ) then
       CPovR_liq = ( CNST_CPvap - CNST_CL ) / CNST_Rvap
       CPovR_ice = ( CNST_CPvap - CNST_CI ) / CNST_Rvap
       CVovR_liq = ( CNST_CVvap - CNST_CL ) / CNST_Rvap
       CVovR_ice = ( CNST_CVvap - CNST_CI ) / CNST_Rvap
    elseif(      EIN_TYPE == 'SIMPLE'  &
            .OR. EIN_TYPE == 'SIMPLE2' ) then
       CPovR_liq = 0.0_RP
       CPovR_ice = 0.0_RP
       CVovR_liq = 0.0_RP
       CVovR_ice = 0.0_RP
    endif
    LovR_liq = LHV / CNST_Rvap
    LovR_ice = LHS / CNST_Rvap

    return
  end subroutine SATURATION_setup

  !-----------------------------------------------------------------------------
  subroutine SATURATION_setrange_SP( Tw, Ti )
    implicit none

    real(SP), intent(in) :: Tw
    real(SP), intent(in) :: Ti
    !---------------------------------------------------------------------------

    SATURATION_ULIMIT_TEMP = real(Tw,kind=RP)
    SATURATION_LLIMIT_TEMP = real(Ti,kind=RP)

    return
  end subroutine SATURATION_setrange_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_setrange_DP( Tw, Ti )
    implicit none

    real(DP), intent(in) :: Tw
    real(DP), intent(in) :: Ti
    !---------------------------------------------------------------------------

    SATURATION_ULIMIT_TEMP = real(Tw,kind=RP)
    SATURATION_LLIMIT_TEMP = real(Ti,kind=RP)

    return
  end subroutine SATURATION_setrange_DP

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine SATURATION_alpha_0D( &
       tem,   &
       alpha  )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    real(RP), intent(in)  :: tem   !< temperature [K]
    real(RP), intent(out) :: alpha !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       if    ( tem > SATURATION_ULIMIT_TEMP ) then
          alpha = 1.0_RP
       elseif( tem < SATURATION_LLIMIT_TEMP ) then
          alpha = 0.0_RP
       else
          fact = ( SATURATION_ULIMIT_TEMP - tem                    ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
       endif

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       fact = ( tem                    - SATURATION_LLIMIT_TEMP ) &
            / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

       alpha = min( max( fact, 0.0_RP ), 1.0_RP )

    endif

    return
  end subroutine SATURATION_alpha_0D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (1D)
  subroutine SATURATION_alpha_1D( &
       ijdim, &
       tem,   &
       alpha  )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem  (ijdim) !< temperature [K]
    real(RP), intent(out) :: alpha(ijdim) !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI

    integer :: ij
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       !$omp parallel default(none),private(ij,fact), &
       !$omp shared(ijdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)

       !$omp do
       do ij = 1, ijdim
          fact = ( SATURATION_ULIMIT_TEMP - tem(ij)                ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij) = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
       enddo
       !$omp end do

       !$omp do
       do ij = 1, ijdim
          if    ( tem(ij) > SATURATION_ULIMIT_TEMP ) then
             alpha(ij) = 1.0_RP
          elseif( tem(ij) < SATURATION_LLIMIT_TEMP ) then
             alpha(ij) = 0.0_RP
          endif
       enddo
       !$omp end do

       !$omp end parallel

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,fact), &
       !$omp shared(ijdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do ij = 1, ijdim
          fact = ( tem(ij)                - SATURATION_LLIMIT_TEMP ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij) = min( max( fact, 0.0_RP ), 1.0_RP )
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_alpha_1D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (2D)
  subroutine SATURATION_alpha_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       alpha  )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem  (ijdim,kdim) !< temperature [K]
    real(RP), intent(out) :: alpha(ijdim,kdim) !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI

    integer :: ij, k
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       !$omp parallel default(none),private(ij,k,fact), &
       !$omp shared(ijdim,kdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do k  = 1, kdim

          !$omp do
          do ij = 1, ijdim
             fact = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
                  / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             alpha(ij,k) = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
          enddo
          !$omp end do

          !$omp do
          do ij = 1, ijdim
             if    ( tem(ij,k) > SATURATION_ULIMIT_TEMP ) then
                alpha(ij,k) = 1.0_RP
             elseif( tem(ij,k) < SATURATION_LLIMIT_TEMP ) then
                alpha(ij,k) = 0.0_RP
             endif
          enddo
          !$omp end do

       enddo
       !$omp end parallel

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,fact), &
       !$omp shared(ijdim,kdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do k  = 1, kdim
       do ij = 1, ijdim
          fact = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij,k) = min( max( fact, 0.0_RP ), 1.0_RP )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_alpha_2D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(tem), 0D
  subroutine SATURATION_dalphadT_0D( &
       tem,      &
       dalpha_dT )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    real(RP), intent(in)  :: tem       !< temperature [K]
    real(RP), intent(out) :: dalpha_dT !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       ! if Tup < tem, dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem )
       ! if Tdn > tem, dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, tem - SATURATION_LLIMIT_TEMP )

       fact = ( SATURATION_ULIMIT_TEMP - tem                    ) &
            / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

       dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                 * 0.5_RP * PI * sin( fact * PI )

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       ! if Tup < tem, dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem )
       ! if Tdn > tem, dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, tem - SATURATION_LLIMIT_TEMP )

       dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

    endif

    return
  end subroutine SATURATION_dalphadT_0D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(tem), 1D
  subroutine SATURATION_dalphadT_1D( &
       ijdim,    &
       tem,      &
       dalpha_dT )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem      (ijdim) !< temperature [K]
    real(RP), intent(out) :: dalpha_dT(ijdim) !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI

    integer :: ij
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       !$omp parallel do default(none),private(ij,fact,lim1,lim2), &
       !$omp shared(ijdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij) - SATURATION_LLIMIT_TEMP )

          fact = ( SATURATION_ULIMIT_TEMP - tem(ij)                ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                        * 0.5_RP * PI * sin( fact * PI )
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,lim1,lim2), &
       !$omp shared(ijdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij) - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_dalphadT_1D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 2D
  subroutine SATURATION_dalphadT_2D( &
       ijdim,    &
       kdim,     &
       tem,      &
       dalpha_dT )
    use mod_cnst, only: &
       CNST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem      (ijdim,kdim) !< temperature [K]
    real(RP), intent(out) :: dalpha_dT(ijdim,kdim) !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI

    integer :: ij, k
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       !$omp parallel do default(none),private(ij,k,fact,lim1,lim2), &
       !$omp shared(ijdim,kdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )

          fact = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij,k) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                            * 0.5_RP * PI * sin( fact * PI )
       enddo
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,lim1,lim2), &
       !$omp shared(ijdim,kdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij,k) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_dalphadT_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (0D)
  subroutine SATURATION_psat_all_0D( &
       tem,   &
       psat   )
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    !---------------------------------------------------------------------------

    call SATURATION_alpha   ( tem, alpha )
    call SATURATION_psat_liq( tem, psatl )
    call SATURATION_psat_ice( tem, psati )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    return
  end subroutine SATURATION_psat_all_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (1D)
  subroutine SATURATION_psat_all_1D( &
       ijdim, &
       tem,   &
       psat   )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(RP) :: alpha(ijdim), psatl(ijdim), psati(ijdim)

    integer :: ij
    !---------------------------------------------------------------------------

    call SATURATION_alpha   ( ijdim, tem(:), alpha(:) )
    call SATURATION_psat_liq( ijdim, tem(:), psatl(:) )
    call SATURATION_psat_ice( ijdim, tem(:), psati(:) )

    do ij = 1, ijdim
       psat(ij) = psatl(ij) * (          alpha(ij) ) &
                + psati(ij) * ( 1.0_RP - alpha(ij) )
    enddo

    return
  end subroutine SATURATION_psat_all_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (2D)
  subroutine SATURATION_psat_all_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_PI,    &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: rtem
    real(RP) :: PI, RTEM00, PSAT0

    integer :: ij, k
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CNST_TEM00
    PSAT0  = CNST_PSAT0

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CNST_PI

       !$omp parallel do default(none),private(ij,k,alpha,rtem,psatl,psati), &
       !$omp shared(ijdim,kdim,psat,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI, &
       !$omp        CPovR_liq,LovR_liq,CPovR_ice,LovR_ice,RTEM00,PSAT0)
       do k  = 1, kdim
       do ij = 1, ijdim

          if    ( tem(ij,k) > SATURATION_ULIMIT_TEMP ) then
             alpha = 1.0_RP
          elseif( tem(ij,k) < SATURATION_LLIMIT_TEMP ) then
             alpha = 0.0_RP
          else
             alpha = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
                   / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             alpha = 0.5_RP * ( 1.0_RP + cos( alpha * PI ) )
          endif

          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )

          psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )

          psat(ij,k) = psatl * (          alpha ) &
                     + psati * ( 1.0_RP - alpha )
       enddo
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,alpha,rtem,psatl,psati), &
       !$omp shared(ijdim,kdim,psat,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP, &
       !$omp        CPovR_liq,LovR_liq,CPovR_ice,LovR_ice,RTEM00,PSAT0)
       do k  = 1, kdim
       do ij = 1, ijdim
          alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
                / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )

          psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )

          psat(ij,k) = psatl * (          alpha ) &
                     + psati * ( 1.0_RP - alpha )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_psat_all_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine SATURATION_psat_liq_0D( &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CNST_TEM00
    PSAT0  = CNST_PSAT0

    rtem = 1.0_RP / max( tem, TEM_MIN )

    psat = PSAT0 * ( tem * RTEM00 )**CPovR_liq &
                 * exp( LovR_liq * ( RTEM00 - rtem ) )

    return
  end subroutine SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_liq_1D_SP( &
       ijdim, &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(SP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(SP) :: rtem
    real(SP) :: RTEM00, PSAT0, TMIN

    integer :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_SP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

!OCL SWP
    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0,TMIN)
    do ij = 1, ijdim
       rtem = 1.0_SP / max( tem(ij), TMIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_liq_1D_SP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_liq_1D_DP( &
       ijdim, &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(DP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(DP) :: rtem
    real(DP) :: RTEM00, PSAT0, TMIN

    integer :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_DP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

!OCL SWP
    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0,TMIN)
    do ij = 1, ijdim
       rtem = 1.0_DP / max( tem(ij), TMIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_liq_1D_DP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_liq_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(SP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(SP) :: rtem
    real(SP) :: RTEM00, PSAT0, TMIN

    integer :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_SP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

!OCL SWP
    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0,TMIN)
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_SP / max( tem(ij,k), TMIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                          * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_liq_2D_SP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_liq_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(DP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(DP) :: rtem
    real(DP) :: RTEM00, PSAT0, TMIN

    integer :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_DP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

!OCL SWP
    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0,TMIN)
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_DP / max( tem(ij,k), TMIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                          * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_liq_2D_DP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine SATURATION_psat_ice_0D( &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CNST_TEM00
    PSAT0  = CNST_PSAT0

    rtem = 1.0_RP / max( tem, TEM_MIN )

    psat = PSAT0 * ( tem * RTEM00 )**CPovR_ice &
                 * exp( LovR_ice * ( RTEM00 - rtem ) )

    return
  end subroutine SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_ice_1D_SP( &
       ijdim, &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(SP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(SP) :: rtem
    real(SP) :: RTEM00, PSAT0, TMIN

    integer :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_SP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0,TMIN)
!OCL SWP
    do ij = 1, ijdim
       rtem = 1.0_SP / max( tem(ij), TMIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_ice_1D_SP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_ice_1D_DP( &
       ijdim, &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(DP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(DP) :: rtem
    real(DP) :: RTEM00, PSAT0, TMIN

    integer :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_DP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

    !$acc kernels async(1) pcopy(psat) pcopyin(tem) if(acc_is_present(psat))
    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0,TMIN)
!OCL SWP
    do ij = 1, ijdim
       rtem = 1.0_DP / max( tem(ij), TMIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_ice_1D_DP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_ice_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(SP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(SP) :: rtem
    real(SP) :: RTEM00, PSAT0, TMIN

    integer :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))
    RTEM00 = 1.0_SP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

    !$acc kernels async(1) pcopy(psat) pcopyin(tem)
    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0,TMIN)
!OCL SWP
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_SP / max( tem(ij,k), TMIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                          * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_ice_2D_SP

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_ice_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
    use mod_cnst, only: &
       CNST_TEM00, &
       CNST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(DP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(DP) :: rtem
    real(DP) :: RTEM00, PSAT0, TMIN

    integer :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem) if(acc_is_present(psat))

    RTEM00 = 1.0_DP / CNST_TEM00
    PSAT0  = CNST_PSAT0
    TMIN   = TEM_MIN

    !$acc kernels async(1) pcopy(psat) pcopyin(tem)
    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0,TMIN)
!OCL SWP
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_DP / max( tem(ij,k), TMIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                          * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine SATURATION_psat_ice_2D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_0D( &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    real(RP), intent(in)  :: pre
    real(RP), intent(in)  :: tem
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: EPSvap
    !---------------------------------------------------------------------------

    EPSvap = CNST_EPSvap

    call SATURATION_psat_liq_0D( tem, psat )

    qsat = EPSvap * psat / ( pre - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine SATURATION_qsat_liq_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_1D_SP( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: pre (ijdim)
    real(SP), intent(in)  :: tem (ijdim)
    real(SP), intent(out) :: qsat(ijdim)

    real(SP) :: psat(ijdim)
    real(SP) :: EPSvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_liq_1D_SP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_SP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq_1D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_1D_DP( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: pre (ijdim)
    real(DP), intent(in)  :: tem (ijdim)
    real(DP), intent(out) :: qsat(ijdim)

    real(DP) :: psat(ijdim)
    real(DP) :: EPSvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_liq_1D_DP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_DP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq_1D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: pre (ijdim,kdim)
    real(SP), intent(in)  :: tem (ijdim,kdim)
    real(SP), intent(out) :: qsat(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim)
    real(SP) :: EPSvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_liq_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_SP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq_2D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: pre (ijdim,kdim)
    real(DP), intent(in)  :: tem (ijdim,kdim)
    real(DP), intent(out) :: qsat(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim)
    real(DP) :: EPSvap

    integer :: ij, k
    !---------------------------------------------------------------------------
   !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_liq_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_DP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq_2D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_0D( &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    real(RP), intent(in)  :: pre
    real(RP), intent(in)  :: tem
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: EPSvap
    !---------------------------------------------------------------------------

    EPSvap = CNST_EPSvap

    call SATURATION_psat_ice_0D( tem, psat )

    qsat = EPSvap * psat / ( pre - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine SATURATION_qsat_ice_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_1D_SP( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: pre (ijdim)
    real(SP), intent(in)  :: tem (ijdim)
    real(SP), intent(out) :: qsat(ijdim)

    real(SP) :: psat(ijdim)
    real(SP) :: EPSvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_ice_1D_SP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_SP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice_1D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_1D_DP( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: pre (ijdim)
    real(DP), intent(in)  :: tem (ijdim)
    real(DP), intent(out) :: qsat(ijdim)

    real(DP) :: psat(ijdim)
    real(DP) :: EPSvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))
    EPSvap = CNST_EPSvap

    call SATURATION_psat_ice_1D_DP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_DP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice_1D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: pre (ijdim,kdim)
    real(SP), intent(in)  :: tem (ijdim,kdim)
    real(SP), intent(out) :: qsat(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim)
    real(SP) :: EPSvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_ice_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_SP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice_2D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
    use mod_cnst, only: &
       CNST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: pre (ijdim,kdim)
    real(DP), intent(in)  :: tem (ijdim,kdim)
    real(DP), intent(out) :: qsat(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim)
    real(DP) :: EPSvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,pre) pcreate(psat) if(acc_is_present(qsat))

    EPSvap = CNST_EPSvap

    call SATURATION_psat_ice_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,pre) if(acc_is_present(qsat))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_DP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice_2D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_rh_SP( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       qv,    &
       rh     )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TMELT
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: rho(ijdim,kdim)
    real(SP), intent(in)  :: tem(ijdim,kdim)
    real(SP), intent(in)  :: qv (ijdim,kdim)
    real(SP), intent(out) :: rh (ijdim,kdim)

    real(SP) :: psat_liq(ijdim,kdim)
    real(SP) :: psat_ice(ijdim,kdim)
    real(SP) :: rh_liq, rh_ice, alpha
    real(SP) :: Rvap, TMELT

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcreate(psat_liq,psat_ice) if(acc_is_present(rh))
    Rvap  = CNST_Rvap
    TMELT = CNST_TMELT

    call SATURATION_psat_liq_2D_SP( ijdim, kdim, tem(:,:), psat_liq(:,:) )
    call SATURATION_psat_ice_2D_SP( ijdim, kdim, tem(:,:), psat_ice(:,:) )

    !$acc kernels pcopy(rh) pcopyin(qv,rho,tem,psat_liq,psat_ice) if(acc_is_present(rh))
    !$omp parallel do default(none),private(ij,k,rh_liq,rh_ice,alpha), &
    !$omp shared(ijdim,kdim,rh,rho,tem,qv,psat_liq,psat_ice,Rvap,TMELT)
    do k  = 1, kdim
    do ij = 1, ijdim
       rh_liq = qv(ij,k) / psat_liq(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       rh_ice = qv(ij,k) / psat_ice(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       alpha  = 0.5_SP + sign(0.5_SP,tem(ij,k)-TMELT)

       rh(ij,k) = (        alpha ) * rh_liq &
                + ( 1.0_SP-alpha ) * rh_ice
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc end data
    return
  end subroutine SATURATION_rh_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_rh_DP( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       qv,    &
       rh     )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TMELT
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: rho(ijdim,kdim)
    real(DP), intent(in)  :: tem(ijdim,kdim)
    real(DP), intent(in)  :: qv (ijdim,kdim)
    real(DP), intent(out) :: rh (ijdim,kdim)

    real(DP) :: psat_liq(ijdim,kdim)
    real(DP) :: psat_ice(ijdim,kdim)
    real(DP) :: rh_liq, rh_ice, alpha
    real(DP) :: Rvap, TMELT

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcreate(psat_liq,psat_ice) if(acc_is_present(rh))
    Rvap  = CNST_Rvap
    TMELT = CNST_TMELT

    call SATURATION_psat_liq_2D_DP( ijdim, kdim, tem(:,:), psat_liq(:,:) )
    call SATURATION_psat_ice_2D_DP( ijdim, kdim, tem(:,:), psat_ice(:,:) )

    !$acc kernels pcopy(rh) pcopyin(qv,rho,tem,psat_liq,psat_ice) if(acc_is_present(rh))
    !$omp parallel do default(none),private(ij,k,rh_liq,rh_ice,alpha), &
    !$omp shared(ijdim,kdim,rh,rho,tem,qv,psat_liq,psat_ice,Rvap,TMELT)
    do k  = 1, kdim
    do ij = 1, ijdim
       rh_liq = qv(ij,k) / psat_liq(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       rh_ice = qv(ij,k) / psat_ice(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       alpha  = 0.5_DP + sign(0.5_DP,tem(ij,k)-TMELT)

       rh(ij,k) = (        alpha ) * rh_liq &
                + ( 1.0_DP-alpha ) * rh_ice
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc end data
    return
  end subroutine SATURATION_rh_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_dewtem( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       qv,    &
       wtem   )
    use mod_cnst, only: &
       CNST_EPSvap, &
       CNST_PSAT0,  &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim)
    real(RP), intent(in)  :: pre (ijdim,kdim)
    real(RP), intent(in)  :: qd  (ijdim,kdim)
    real(RP), intent(in)  :: qv  (ijdim,kdim)
    real(RP), intent(out) :: wtem(ijdim,kdim)       ! dew point temperature

    real(RP) :: prev(ijdim,kdim)
    logical  :: flag(ijdim,kdim)
    real(RP) :: rtem, lv, psatl, f, dfdtem
    real(RP) :: EPSvap, RTEM00, PSAT0

    real(RP), parameter :: criteria = 1.E-8_RP
    integer,  parameter :: itelim   = 20

    integer :: ite
    integer :: ij, k
    !---------------------------------------------------------------------------

    EPSvap = CNST_EPSvap
    RTEM00 = 1.0_RP / CNST_TEM00
    PSAT0  = CNST_PSAT0

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,prev,wtem,flag,tem,pre,qd,qv,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       prev(ij,k) = pre(ij,k) * qv(ij,k) / ( EPSvap * qd(ij,k) + qv(ij,k) )

       wtem(ij,k) = tem(ij,k) * 0.98_RP

       flag(ij,k) = .false.
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel default(none),private(ij,k,rtem,lv,psatl,f,dfdtem), &
    !$omp shared(ite,ijdim,kdim,wtem,flag,prev,CPovR_liq,LovR_liq,RTEM00,PSAT0)
    do ite = 1, itelim
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          if ( .NOT. flag(ij,k) ) then
             rtem   = 1.0_RP / max( wtem(ij,k), TEM_MIN )

             lv     = LovR_liq + wtem(ij,k) * CPovR_liq

             psatl  = PSAT0 * ( wtem(ij,k) * RTEM00 )**CPovR_liq &
                            * exp( LovR_liq * ( RTEM00 - rtem ) )

             f      = psatl - prev(ij,k)
             dfdtem = lv * psatl * rtem**2

             wtem(ij,k) = wtem(ij,k) - f / dfdtem

             if ( abs(f) < prev(ij,k) * criteria ) then
                flag(ij,k) = .true.
             endif
          endif
       enddo
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    return
  end subroutine SATURATION_dewtem

  !-----------------------------------------------------------------------------
  subroutine SATURATION_adjustment_SP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhoge,     &
       rhogq,     &
       tem,       &
       q,         &
       qd,        &
       gsgam2,    &
       ice_adjust )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       I_QI,              &
       LHV,               &
       LHF
    use mod_thrmdyn, only: &
       THRMDYN_cv, &
       thrmdyn_wait_on => acc_wait_on
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(SP), intent(in)    :: rhog  (ijdim,kdim)
    real(SP), intent(inout) :: rhoge (ijdim,kdim)
    real(SP), intent(inout) :: rhogq (ijdim,kdim,nqmax)
    real(SP), intent(inout) :: tem   (ijdim,kdim)
    real(SP), intent(inout) :: q     (ijdim,kdim,nqmax)
    real(SP), intent(in)    :: qd    (ijdim,kdim)
    real(SP), intent(in)    :: gsgam2(ijdim,kdim)
    logical,  intent(in)    :: ice_adjust

    real(RP) :: tem0     (ijdim,kdim)
    real(RP) :: qd0      (ijdim,kdim)
    real(RP) :: q0       (ijdim,kdim,nqmax)
    real(RP) :: ein_moist(ijdim,kdim)
    real(RP) :: qsum     (ijdim,kdim)
    real(RP) :: CVtot    (ijdim,kdim)
    real(RP) :: rho      (ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(tem0,qd0,q0,ein_moist,qsum,CVtot,rho) &
    !$acc pcopyin(rhog,qd,gsgam2)

    call PROF_rapstart('____Saturation_Adjustment',2)
    thrmdyn_wait_on  = .false.

    !$acc kernels async(1) pcopy(q0) pcopyin(q)
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,q0,q)
    do nq = NQW_STR, NQW_END
    do k  = 1, kdim
    do ij = 1, ijdim
       q0(ij,k,nq) = q(ij,k,nq)
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    ! ein_moist = U1(rho,qsum,T1) : 'unsaturated temperature'
    if ( I_QI > 0 .AND. ice_adjust ) then
       !$acc kernels async(1) pcopy(qd0,ein_moist,qsum,q0) pcopyin(rhog,rhoge,q)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,qd0,qd,ein_moist,qsum,rhoge,rhog,q0,q,LHV,LHF,I_QV,I_QC,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          qd0      (ij,k) = qd(ij,k)

          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q0(ij,k,I_QV) * LHV &
                          - q0(ij,k,I_QI) * LHF

          qsum     (ij,k) = q0(ij,k,I_QV) &
                          + q0(ij,k,I_QC) &
                          + q0(ij,k,I_QI)

          q0(ij,k,I_QV) = qsum(ij,k)
          q0(ij,k,I_QC) = 0.0_RP
          q0(ij,k,I_QI) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    else
       !$acc kernels async(1) pcopy(qd0,ein_moist,qsum,q0) pcopyin(rhog,rhoge)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,qd0,qd,ein_moist,qsum,rhoge,rhog,q0,q,LHV,I_QV,I_QC)
       do k  = 1, kdim
       do ij = 1, ijdim
          qd0      (ij,k) = qd(ij,k)

          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q0(ij,k,I_QV) * LHV

          qsum     (ij,k) = q0(ij,k,I_QV) &
                          + q0(ij,k,I_QC)

          q0(ij,k,I_QV) = qsum(ij,k)
          q0(ij,k,I_QC) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd0  (:,:),   & ! [IN]
                     q0   (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(tem0,rho) pcopyin(rhog,gsgam2,q0,ein_moist,CVtot)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,tem0,ein_moist,q0,rhog,gsgam2,CVtot,I_QV,LHV)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho (ij,k) = rhog(ij,k) / gsgam2(ij,k)
       tem0(ij,k) = ( ein_moist(ij,k) - q0(ij,k,I_QV) * LHV ) / CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( I_QI > 0 .AND. ice_adjust ) then
#ifdef _OPENACC
       call satadjust_all_ACC( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#else
       call satadjust_all    ( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#endif
    else
#ifdef _OPENACC
       call satadjust_liq_ACC( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#else
       call satadjust_liq    ( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#endif
    endif

    !$acc kernels async(1) pcopy(rhogq,q) pcopyin(rhog,q0)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhogq,rhog,q,q0,I_QV,I_QC)
    do k  = 1, kdim
    do ij = 1, ijdim
       q    (ij,k,I_QV) = q0  (ij,k,I_QV)
       q    (ij,k,I_QC) = q0  (ij,k,I_QC)

       rhogq(ij,k,I_QV) = rhog(ij,k) * q0(ij,k,I_QV)
       rhogq(ij,k,I_QC) = rhog(ij,k) * q0(ij,k,I_QC)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( I_QI > 0 .AND. ice_adjust ) then
       !$acc kernels async(1) pcopy(rhogq,q) pcopyin(rhog,q0)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,rhogq,rhog,q,q0,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          q    (ij,k,I_QI) = q0  (ij,k,I_QI)

          rhogq(ij,k,I_QI) = rhog(ij,k) * q0(ij,k,I_QI)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd0  (:,:),   & ! [IN]
                     q0   (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(rhoge,tem) pcopyin(tem0,CVtot,rhog)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhoge,rhog,tem,tem0,CVtot)
    do k  = 1, kdim
    do ij = 1, ijdim
       tem  (ij,k) = tem0(ij,k)

       rhoge(ij,k) = rhog(ij,k) * tem0(ij,k) * CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    thrmdyn_wait_on  = .true.
    if ( acc_wait_on ) then
         !$acc wait(1)
    endif
    call PROF_rapend  ('____Saturation_Adjustment',2)

    !$acc end data
    return
  end subroutine SATURATION_adjustment_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_adjustment_DP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhoge,     &
       rhogq,     &
       tem,       &
       q,         &
       qd,        &
       gsgam2,    &
       ice_adjust )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       I_QI,              &
       LHV,               &
       LHF
    use mod_thrmdyn, only: &
       THRMDYN_cv,        &
       thrmdyn_wait_on => acc_wait_on
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(DP), intent(in)    :: rhog  (ijdim,kdim)
    real(DP), intent(inout) :: rhoge (ijdim,kdim)
    real(DP), intent(inout) :: rhogq (ijdim,kdim,nqmax)
    real(DP), intent(inout) :: tem   (ijdim,kdim)
    real(DP), intent(inout) :: q     (ijdim,kdim,nqmax)
    real(DP), intent(in)    :: qd    (ijdim,kdim)
    real(DP), intent(in)    :: gsgam2(ijdim,kdim)
    logical,  intent(in)    :: ice_adjust

    real(RP) :: tem0     (ijdim,kdim)
    real(RP) :: qd0      (ijdim,kdim)
    real(RP) :: q0       (ijdim,kdim,nqmax)
    real(RP) :: ein_moist(ijdim,kdim)
    real(RP) :: qsum     (ijdim,kdim)
    real(RP) :: CVtot    (ijdim,kdim)
    real(RP) :: rho      (ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(tem0,qd0,q0,ein_moist,qsum,CVtot,rho) &
    !$acc pcopyin(rhog,qd,gsgam2)

    call PROF_rapstart('____Saturation_Adjustment',2)
    thrmdyn_wait_on  = .false.

    !$acc kernels async(1) pcopy(q0) pcopyin(q)
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,q0,q)
    do nq = NQW_STR, NQW_END
    do k  = 1, kdim
    do ij = 1, ijdim
       q0(ij,k,nq) = q(ij,k,nq)
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    ! ein_moist = U1(rho,qsum,T1) : 'unsaturated temperature'
    if ( I_QI > 0 .AND. ice_adjust ) then
       !$acc kernels async(1) pcopy(qd0,ein_moist,qsum,q0) pcopyin(rhog,rhoge,q)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,qd0,qd,ein_moist,qsum,rhoge,rhog,q0,q,LHV,LHF,I_QV,I_QC,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          qd0      (ij,k) = qd(ij,k)

          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q0(ij,k,I_QV) * LHV &
                          - q0(ij,k,I_QI) * LHF

          qsum     (ij,k) = q0(ij,k,I_QV) &
                          + q0(ij,k,I_QC) &
                          + q0(ij,k,I_QI)

          q0(ij,k,I_QV) = qsum(ij,k)
          q0(ij,k,I_QC) = 0.0_RP
          q0(ij,k,I_QI) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    else
       !$acc kernels async(1) pcopy(qd0,ein_moist,qsum,q0) pcopyin(rhog,rhoge)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,qd0,qd,ein_moist,qsum,rhoge,rhog,q0,q,LHV,I_QV,I_QC)
       do k  = 1, kdim
       do ij = 1, ijdim
          qd0      (ij,k) = qd(ij,k)

          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q0(ij,k,I_QV) * LHV

          qsum     (ij,k) = q0(ij,k,I_QV) &
                          + q0(ij,k,I_QC)

          q0(ij,k,I_QV) = qsum(ij,k)
          q0(ij,k,I_QC) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd0  (:,:),   & ! [IN]
                     q0   (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(tem0,rho) pcopyin(rhog,gsgam2,q0,ein_moist,CVtot)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,tem0,ein_moist,q0,rhog,gsgam2,CVtot,I_QV,LHV)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho (ij,k) = rhog(ij,k) / gsgam2(ij,k)
       tem0(ij,k) = ( ein_moist(ij,k) - q0(ij,k,I_QV) * LHV ) / CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( I_QI > 0 .AND. ice_adjust ) then
#ifdef _OPENACC
       call satadjust_all_ACC( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#else
       call satadjust_all    ( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#endif
    else
#ifdef _OPENACC
       call satadjust_liq_ACC( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#else
       call satadjust_liq    ( ijdim,           & ! [IN]
                               kdim,            & ! [IN]
                               rho      (:,:),  & ! [IN]
                               ein_moist(:,:),  & ! [IN]
                               qsum     (:,:),  & ! [IN]
                               tem0     (:,:),  & ! [INOUT]
                               q0       (:,:,:) ) ! [INOUT]
#endif
    endif

    !$acc kernels async(1) pcopy(rhogq,q) pcopyin(rhog,q0)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhogq,rhog,q,q0,I_QV,I_QC)
    do k  = 1, kdim
    do ij = 1, ijdim
       q    (ij,k,I_QV) = q0  (ij,k,I_QV)
       q    (ij,k,I_QC) = q0  (ij,k,I_QC)

       rhogq(ij,k,I_QV) = rhog(ij,k) * q0(ij,k,I_QV)
       rhogq(ij,k,I_QC) = rhog(ij,k) * q0(ij,k,I_QC)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( I_QI > 0 .AND. ice_adjust ) then
       !$acc kernels async(1) pcopy(rhogq,q) pcopyin(rhog,q0)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,rhogq,rhog,q,q0,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          q    (ij,k,I_QI) = q0  (ij,k,I_QI)

          rhogq(ij,k,I_QI) = rhog(ij,k) * q0(ij,k,I_QI)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd0  (:,:),   & ! [IN]
                     q0   (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(rhoge,tem) pcopyin(tem0,CVtot,rhog)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhoge,rhog,tem,tem0,CVtot)
    do k  = 1, kdim
    do ij = 1, ijdim
       tem  (ij,k) = tem0(ij,k)

       rhoge(ij,k) = rhog(ij,k) * tem0(ij,k) * CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    thrmdyn_wait_on  = .true.
    if ( acc_wait_on ) then
         !$acc wait(1)
    endif
    call PROF_rapend  ('____Saturation_Adjustment',2)

    !$acc end data
    return
  end subroutine SATURATION_adjustment_DP

  !-----------------------------------------------------------------------------
  subroutine satadjust_all_ACC( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
    use mod_cnst, only: &
       CNST_EPS,   &
       CNST_CVdry, &
       CNST_Rvap,  &
       CNST_PSAT0, &
       CNST_TEM00
    use mod_adm, only: &
       ADM_proc_stop,    &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       I_QI,              &
       CVW,               &
       LHV,               &
       LHF
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd      (ijdim,kdim)
    integer  :: exe_flag(ijdim,kdim)

    real(RP) :: rtem, alpha, psatl, psati, psat, qsatl, qsati, qsat
    real(RP) :: CVtot, Emoist_new, dtemp, lim1, lim2
    real(RP) :: dalpha_dT, dqsatl_dT, dqsati_dT, dqsat_dT, dqc_dT, dqi_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    integer, parameter :: itelim = 100

    real(RP) :: dtemp_criteria
    integer  :: counts
    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(qd,exe_flag) &
    !$acc pcopyin(rho,Emoist,qsum)

    call PROF_rapstart('_____SatAdj_CORE',3)

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CNST_EPS
    CVdry  = CNST_CVdry
    Rvap   = CNST_Rvap
    PSAT0  = CNST_PSAT0
    RTEM00 = 1.0_RP / CNST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !$acc kernels async(1) pcopy(tem,q,exe_flag) pcopyin(rho,Emoist,qsum,qd,CVW)
    do k  = kmin, kmax
    do ij = 1, ijdim
       alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
             / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
       psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
       psat  = psatl * (          alpha ) &
             + psati * ( 1.0_RP - alpha )

       qsat  = psat  / ( rho(ij,k) * Rvap * tem(ij,k) )

       exe_flag(ij,k) = 0

       if ( qsum(ij,k)-qsat > EPS ) then
          exe_flag(ij,k) = 1

          !$acc loop seq
          do ite = 1, itelim
             ! treat nvfort bug??
             if ( exe_flag(ij,k) == 0  )  cycle
             if ( tem(ij,k)*0.0_RP /= 0.0_RP ) cycle
             alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
                   / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
             alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

             rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

             psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
             psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
             psat  = psatl * (          alpha ) &
                   + psati * ( 1.0_RP - alpha )

             qsatl = psatl / ( rho(ij,k) * Rvap * tem(ij,k) )
             qsati = psati / ( rho(ij,k) * Rvap * tem(ij,k) )
             qsat  = psat  / ( rho(ij,k) * Rvap * tem(ij,k) )

             ! Separation
             q(ij,k,I_QV) = qsat
             q(ij,k,I_QC) = ( qsum(ij,k)-qsat ) * (        alpha )
             q(ij,k,I_QI) = ( qsum(ij,k)-qsat ) * ( 1.0_RP-alpha )

             CVtot = qd(ij,k) * CVdry
             !$acc loop seq
             do nq = NQW_STR, NQW_END
                CVtot = CVtot + q(ij,k,nq) * CVW(nq)
             enddo

             Emoist_new = tem(ij,k) * CVtot + q(ij,k,I_QV) * LHV - q(ij,k,I_QI) * LHF

             ! dX/dT
             lim1       = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
             lim2       = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )
             dalpha_dT  = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             dqsatl_dT  = ( LovR_liq / tem(ij,k)**2 + CVovR_liq / tem(ij,k) ) * qsatl
             dqsati_dT  = ( LovR_ice / tem(ij,k)**2 + CVovR_ice / tem(ij,k) ) * qsati
             dqsat_dT   = qsatl * dalpha_dT + dqsatl_dT * (        alpha ) &
                        - qsati * dalpha_dT + dqsati_dT * ( 1.0_RP-alpha )

             dqc_dT     =  ( qsum(ij,k)-qsat ) * dalpha_dT - dqsat_dT * (        alpha )
             dqi_dT     = -( qsum(ij,k)-qsat ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha )

             dCVtot_dT  = dqsat_dT * CVW(I_QV) &
                        + dqc_dT   * CVW(I_QC) &
                        + dqi_dT   * CVW(I_QI)

             dEmoist_dT = tem(ij,k) * dCVtot_dT &
                        + CVtot                 &
                        + dqsat_dT  * LHV       &
                        - dqi_dT    * LHF

             dtemp      = ( Emoist_new - Emoist(ij,k) ) / dEmoist_dT

             tem(ij,k)  = tem(ij,k) - dtemp

             if ( abs(dtemp) < dtemp_criteria ) then
                exe_flag(ij,k) = 0
                !exit
             endif

             !if( tem(ij,k)*0.0_RP /= 0.0_RP ) exit
          enddo
       endif
    enddo
    enddo
    !$acc end kernels

    counts = 0
    !$acc kernels async(1) pcopyin(exe_flag)
    !$acc loop reduction(+:counts)
    do k  = kmin, kmax
    !$acc loop reduction(+:counts)
    do ij = 1, ijdim
       counts = counts + exe_flag(ij,k)
    enddo
    enddo
    !$acc end kernels
    !$acc wait(1)
    if ( counts /= 0 ) then
       write(*,*) 'xxx [satadjust_all_ACC] not converged! Stop.', counts
       !$acc update self(exe_flag,rho,tem,q)
       do k  = kmin, kmax
       do ij = 1, ijdim
          if ( exe_flag(ij,k) > 0 ) then
             write(*,*) ij,k,rho(ij,k),tem(ij,k),q(ij,k,I_QV),q(ij,k,I_QC),q(ij,k,I_QI)
          endif
       enddo
       enddo

       call ADM_proc_stop
    endif

    call PROF_rapend  ('_____SatAdj_CORE',3)

    !$acc end data
    return
  end subroutine satadjust_all_ACC

  !-----------------------------------------------------------------------------
  subroutine satadjust_liq_ACC( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
    use mod_cnst, only: &
       CNST_EPS,   &
       CNST_CVdry, &
       CNST_Rvap,  &
       CNST_PSAT0, &
       CNST_TEM00
    use mod_adm, only: &
       ADM_proc_stop,    &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       CVW,               &
       LHV
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd      (ijdim,kdim)
    integer  :: exe_flag(ijdim,kdim)

    real(RP) :: rtem, psat, qsat
    real(RP) :: CVtot, Emoist_new, dtemp
    real(RP) :: dqsat_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    integer, parameter :: itelim = 100

    real(RP) :: dtemp_criteria
    integer  :: counts
    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(qd,exe_flag,counts) &
    !$acc pcopyin(rho,Emoist,qsum)

    call PROF_rapstart('_____SatAdj_CORE',3)

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CNST_EPS
    CVdry  = CNST_CVdry
    Rvap   = CNST_Rvap
    PSAT0  = CNST_PSAT0
    RTEM00 = 1.0_RP / CNST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !!$acc kernels pcopy(tem,q,exe_flag) pcopyin(rho,Emoist,qsum,qd,CVW)
    !$acc parallel loop async(1) collapse(2) private(ite,rtem,psat,qsat,CVtot,Emoist_new,dqsat_dT,dCVtot_dT,dEmoist_dT,dtemp,ite,nq)
    do k  = kmin, kmax
    do ij = 1, ijdim

       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psat = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

       qsat = psat / ( rho(ij,k) * Rvap * tem(ij,k) )

       exe_flag(ij,k) = 0

       if ( qsum(ij,k)-qsat > EPS ) then
          exe_flag(ij,k) = 1

          !$acc loop seq
          do ite = 1, itelim
            ! treat nvfort bug??
             if ( exe_flag(ij,k) == 0  )  cycle
             if ( tem(ij,k)*0.0_RP /= 0.0_RP ) cycle
             rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

             psat = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

             qsat = psat / ( rho(ij,k) * Rvap * tem(ij,k) )

             ! Separation
             q(ij,k,I_QV) = qsat
             q(ij,k,I_QC) = qsum(ij,k)-qsat

             CVtot = qd(ij,k) * CVdry
             !$acc loop seq
             do nq = NQW_STR, NQW_END
                CVtot = CVtot + q(ij,k,nq) * CVW(nq)
             enddo

             Emoist_new = tem(ij,k) * CVtot + q(ij,k,I_QV) * LHV

             ! dX/dT
             dqsat_dT   = ( LovR_liq / tem(ij,k)**2 + CVovR_liq / tem(ij,k) ) * qsat

             dCVtot_dT  = dqsat_dT * ( CVW(I_QV) - CVW(I_QC) )

             dEmoist_dT = tem(ij,k) * dCVtot_dT &
                        + CVtot                 &
                        + dqsat_dT  * LHV

             dtemp      = ( Emoist_new - Emoist(ij,k) ) / dEmoist_dT

             tem(ij,k)  = tem(ij,k) - dtemp
             ! treat nvfort bug??
             if ( abs(dtemp) < dtemp_criteria ) then
                exe_flag(ij,k) = 0
                !exit
             endif

             !if( tem(ij,k)*0.0_RP /= 0.0_RP ) exit
          enddo
       endif
    enddo
    enddo
    !!$acc end kernels
    !$acc end parallel

    counts = 0
    !$acc update device(counts) async(1)
    !$acc kernels pcopyin(exe_flag) async(1)
    !$acc loop reduction(+:counts)
    do k  = kmin, kmax
    !$acc loop reduction(+:counts)
    do ij = 1, ijdim
       counts = counts + exe_flag(ij,k)
    enddo
    enddo
    !$acc end kernels


    !$acc update self(counts) async(1)
    !$acc wait(1)
    if ( counts /= 0 ) then
       write(*,*) 'xxx [satadjust_liq_ACC] not converged! Stop.', counts
       !$acc update self(exe_flag,rho,tem,q)
       do k  = kmin, kmax
       do ij = 1, ijdim
          if ( exe_flag(ij,k) > 0 ) then
             write(*,*) ij,k,rho(ij,k),tem(ij,k),q(ij,k,I_QV),q(ij,k,I_QC)
          endif
       enddo
       enddo

       call ADM_proc_stop
    endif

    call PROF_rapend  ('_____SatAdj_CORE',3)

    !$acc end data
    return
  end subroutine satadjust_liq_ACC

  !-----------------------------------------------------------------------------
  subroutine satadjust_all( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
    use mod_cnst, only: &
       CNST_EPS,   &
       CNST_CVdry, &
       CNST_Rvap,  &
       CNST_PSAT0, &
       CNST_TEM00
    use mod_adm, only: &
       ADM_proc_stop,    &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       I_QI,              &
       CVW,               &
       LHV,               &
       LHF
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd  (ijdim,kdim)
    real(RP) :: qsat(ijdim,kdim)

    real(RP) :: rtem, psatl, psati, psat
    real(RP) :: Emoist_new, dtemp, lim1, lim2
    real(RP) :: dalpha_dT, dqsatl_dT, dqsati_dT, dqsat_dT, dqc_dT, dqi_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    real(RP) :: dtemp_criteria
    integer, parameter :: itelim  = 100

    ! SIMD parameter
    integer, parameter :: simdlen = 512
    integer            :: blk, vec, veclen
    integer            :: n, nmax(kdim)

    ! List vector
    integer  :: idx_list   (ijdim,kdim)
    real(SP) :: exe_flag   (ijdim,kdim)
    real(RP) :: qd_list    (ijdim,kdim)
    real(RP) :: qsat_list  (ijdim,kdim)
    real(RP) :: rho_list   (ijdim,kdim)
    real(RP) :: Emoist_list(ijdim,kdim)
    real(RP) :: qsum_list  (ijdim,kdim)
    real(RP) :: tem_list   (ijdim,kdim)
    real(RP) :: q_list     (ijdim,kdim,nqmax)

    real(RP) :: alpha      (ijdim,kdim)
    real(RP) :: qsatl      (ijdim,kdim)
    real(RP) :: qsati      (ijdim,kdim)
    real(RP) :: CVtot      (ijdim,kdim)

    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------

    call PROF_rapstart('_____SatAdj_CORE',3)

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CNST_EPS
    CVdry  = CNST_CVdry
    Rvap   = CNST_Rvap
    PSAT0  = CNST_PSAT0
    RTEM00 = 1.0_RP / CNST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !$omp parallel do default(none), &
    !$omp private(ij,k,nq,n,blk,ite,vec,veclen,                                               &
    !$omp         rtem,psatl,psati,psat,Emoist_new,dtemp,lim1,lim2,                           &
    !$omp         dalpha_dT,dqsatl_dT,dqsati_dT,dqsat_dT,dqc_dT,dqi_dT,dCVtot_dT,dEmoist_dT), &
    !$omp shared (ijdim,kmin,kmax,NQW_STR,NQW_END,nmax,                           &
    !$omp         tem,q,qd,qsat,rho,Emoist,qsum,                                  &
    !$omp         tem_list,q_list,qsat_list,exe_flag,idx_list,                    &
    !$omp         alpha,qsatl,qsati,CVtot,qd_list,rho_list,Emoist_list,qsum_list, &
    !$omp         CPovR_liq,CVovR_liq,LovR_liq,CPovR_ice,CVovR_ice,LovR_ice,      &
    !$omp         I_QV,I_QC,I_QI,CVW,LHV,LHF,RTEM00,PSAT0,Rvap,CVdry,EPS,         &
    !$omp         SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,dtemp_criteria    )
    do k  = kmin, kmax

       !OCL loop_fission_target(ls)
       do ij = 1, ijdim
          alpha(ij,k) = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
                      / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
          alpha(ij,k) = min( max( alpha(ij,k), 0.0_RP ), 1.0_RP )

          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
          psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
          psat  = psatl * (          alpha(ij,k) ) &
                + psati * ( 1.0_RP - alpha(ij,k) )

          qsat(ij,k) = psat  / ( rho(ij,k) * Rvap * tem(ij,k) )
       enddo

       n = 0
!OCL SIMD
       do ij = 1, ijdim
          if ( qsum(ij,k)-qsat(ij,k) > EPS ) then
             n = n + 1

             idx_list(n,k) = ij
          endif
       enddo

       nmax(k) = n

       do n = 1, nmax(k)
          ij = idx_list(n,k)

          exe_flag   (n,k) = 1.0
          qd_list    (n,k) = qd    (ij,k)
          qsat_list  (n,k) = qsat  (ij,k)
          rho_list   (n,k) = rho   (ij,k)
          Emoist_list(n,k) = Emoist(ij,k)
          qsum_list  (n,k) = qsum  (ij,k)
          tem_list   (n,k) = tem   (ij,k)
       enddo

       do nq = NQW_STR, NQW_END
       do n  = 1, nmax(k)
          ij = idx_list(n,k)

          q_list(n,k,nq) = q(ij,k,nq)
       enddo
       enddo

       do blk = 1, nmax(k), simdlen ! cache blocking
          veclen = min( nmax(k)-blk+1, simdlen )

          iteration: do ite = 1, itelim

             !OCL loop_fission_target(ls)
             do vec = 1, veclen
                n = blk + vec - 1

                alpha(n,k) = ( tem_list(n,k)          - SATURATION_LLIMIT_TEMP ) &
                                / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
                alpha(n,k) = min( max( alpha(n,k), 0.0_RP ), 1.0_RP )

                rtem = 1.0_RP / max( tem_list(n,k), TEM_MIN )

                psatl = PSAT0 * ( tem_list(n,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
                psati = PSAT0 * ( tem_list(n,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
                psat  = psatl * (          alpha(n,k) ) &
                      + psati * ( 1.0_RP - alpha(n,k) )

                qsatl    (n,k) = psatl / ( rho_list(n,k) * Rvap * tem_list(n,k) )
                qsati    (n,k) = psati / ( rho_list(n,k) * Rvap * tem_list(n,k) )
                qsat_list(n,k) = psat  / ( rho_list(n,k) * Rvap * tem_list(n,k) )

                ! Separation
                q_list(n,k,I_QV) = qsat_list(n,k)
                q_list(n,k,I_QC) = ( qsum_list(n,k)-qsat_list(n,k) ) * (        alpha(n,k) )
                q_list(n,k,I_QI) = ( qsum_list(n,k)-qsat_list(n,k) ) * ( 1.0_RP-alpha(n,k) )

                CVtot(n,k) = qd_list(n,k) * CVdry
             enddo

             do nq  = NQW_STR, NQW_END
             do vec = 1, veclen
                n = blk + vec - 1

                CVtot(n,k) = CVtot(n,k) + q_list(n,k,nq) * CVW(nq)
             enddo
             enddo

             !OCL loop_fission_target(ls)
             do vec = 1, veclen
                n = blk + vec - 1

                Emoist_new = tem_list(n,k) * CVtot(n,k) + q_list(n,k,I_QV) * LHV - q_list(n,k,I_QI) * LHF

                ! dX/dT
                lim1       = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem_list(n,k) )
                lim2       = 0.5_RP + sign( 0.5_RP, tem_list(n,k) - SATURATION_LLIMIT_TEMP )
                dalpha_dT  = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

                dqsatl_dT  = ( LovR_liq / tem_list(n,k)**2 + CVovR_liq / tem_list(n,k) ) * qsatl(n,k)
                dqsati_dT  = ( LovR_ice / tem_list(n,k)**2 + CVovR_ice / tem_list(n,k) ) * qsati(n,k)
                dqsat_dT   = qsatl(n,k) * dalpha_dT + dqsatl_dT * (        alpha(n,k) ) &
                           - qsati(n,k) * dalpha_dT + dqsati_dT * ( 1.0_RP-alpha(n,k) )

                dqc_dT     =  ( qsum_list(n,k)-qsat_list(n,k) ) * dalpha_dT - dqsat_dT * (        alpha(n,k) )
                dqi_dT     = -( qsum_list(n,k)-qsat_list(n,k) ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha(n,k) )

                dCVtot_dT  = dqsat_dT * CVW(I_QV) &
                           + dqc_dT   * CVW(I_QC) &
                           + dqi_dT   * CVW(I_QI)

                dEmoist_dT = tem_list(n,k) * dCVtot_dT &
                           + CVtot(n,k)                &
                           + dqsat_dT  * LHV           &
                           - dqi_dT    * LHF

                dtemp      = ( Emoist_new - Emoist_list(n,k) ) / dEmoist_dT

                tem_list(n,k) = tem_list(n,k) - dtemp

                exe_flag(n,k) = 0.5_RP + sign(0.5_RP, abs(dtemp)-dtemp_criteria)
             enddo

             if ( sum(exe_flag(blk:blk+veclen-1,k)) <= 0.0 ) then
                exit iteration
             endif

          enddo iteration ! iteration

       enddo ! outer block

       if ( sum(exe_flag(1:nmax(k),k)) > 0.0 ) then
          write(*,*) 'xxx [satadjust_all] not converged! Stop.', k, int(sum(exe_flag(1:nmax(k),k)))
          call ADM_proc_stop
       endif

       do n = 1, nmax(k)
          ij = idx_list(n,k)

          tem(ij,k)      = tem_list(n,k)
          q  (ij,k,I_QV) = q_list  (n,k,I_QV)
          q  (ij,k,I_QC) = q_list  (n,k,I_QC)
          q  (ij,k,I_QI) = q_list  (n,k,I_QI)
       enddo

    enddo ! k loop
    !$omp end parallel do

    call PROF_rapend  ('_____SatAdj_CORE',3)

    return
  end subroutine satadjust_all

  !-----------------------------------------------------------------------------
  subroutine satadjust_liq( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
    use mod_cnst, only: &
       CNST_EPS,   &
       CNST_CVdry, &
       CNST_Rvap,  &
       CNST_PSAT0, &
       CNST_TEM00
    use mod_adm, only: &
       ADM_proc_stop,    &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QC,              &
       CVW,               &
       LHV
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd  (ijdim,kdim)
    real(RP) :: qsat(ijdim,kdim)

    real(RP) :: rtem, psat
    real(RP) :: Emoist_new, dtemp
    real(RP) :: dqsat_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    real(RP) :: dtemp_criteria
    integer, parameter :: itelim = 100

    ! SIMD parameter
    integer, parameter :: simdlen = 512
    integer            :: blk, vec, veclen
    integer            :: n, nmax(kdim)

    ! List vector
    integer  :: idx_list   (ijdim,kdim)
    real(SP) :: exe_flag   (ijdim,kdim)
    real(RP) :: qd_list    (ijdim,kdim)
    real(RP) :: qsat_list  (ijdim,kdim)
    real(RP) :: rho_list   (ijdim,kdim)
    real(RP) :: Emoist_list(ijdim,kdim)
    real(RP) :: qsum_list  (ijdim,kdim)
    real(RP) :: tem_list   (ijdim,kdim)
    real(RP) :: q_list     (ijdim,kdim,nqmax)

    real(RP) :: CVtot      (ijdim,kdim)

    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------

    call PROF_rapstart('_____SatAdj_CORE',3)

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CNST_EPS
    CVdry  = CNST_CVdry
    Rvap   = CNST_Rvap
    PSAT0  = CNST_PSAT0
    RTEM00 = 1.0_RP / CNST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !$omp parallel do default(none), &
    !$omp private(ij,k,nq,n,blk,ite,vec,veclen,   &
    !$omp         rtem,psat,Emoist_new,dtemp,     &
    !$omp         dqsat_dT,dCVtot_dT,dEmoist_dT), &
    !$omp shared (ijdim,kmin,kmax,NQW_STR,NQW_END,nmax,          &
    !$omp         tem,q,qd,qsat,rho,Emoist,qsum,                 &
    !$omp         tem_list,q_list,qsat_list,exe_flag,idx_list,   &
    !$omp         CVtot,qd_list,rho_list,Emoist_list,qsum_list,  &
    !$omp         CPovR_liq,CVovR_liq,LovR_liq,                  &
    !$omp         I_QV,I_QC,CVW,LHV,RTEM00,PSAT0,Rvap,CVdry,EPS, &
    !$omp         dtemp_criteria                                 )
    do k  = kmin, kmax

       do ij = 1, ijdim
          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psat = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

          qsat(ij,k) = psat / ( rho(ij,k) * Rvap * tem(ij,k) )
       enddo

       n = 0
!OCL SIMD
       do ij = 1, ijdim
          if ( qsum(ij,k)-qsat(ij,k) > EPS ) then
             n = n + 1

             idx_list(n,k) = ij
          endif
       enddo
       nmax(k) = n

       do n = 1, nmax(k)
          ij = idx_list(n,k)

          exe_flag   (n,k) = 1.0
          qd_list    (n,k) = qd    (ij,k)
          qsat_list  (n,k) = qsat  (ij,k)
          rho_list   (n,k) = rho   (ij,k)
          Emoist_list(n,k) = Emoist(ij,k)
          qsum_list  (n,k) = qsum  (ij,k)
          tem_list   (n,k) = tem   (ij,k)
       enddo

       do nq = NQW_STR, NQW_END
       do n = 1, nmax(k)
          ij = idx_list(n,k)

          q_list(n,k,nq) = q(ij,k,nq)
       enddo
       enddo

       do blk = 1, nmax(k), simdlen ! cache blocking
          veclen = min( nmax(k)-blk+1, simdlen )

          iteration: do ite = 1, itelim

             !OCL loop_fission_target(ls)
             do vec = 1, veclen
                n = blk + vec - 1

                rtem = 1.0_RP / max( tem_list(n,k), TEM_MIN )

                psat = PSAT0 * ( tem_list(n,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

                qsat_list(n,k) = psat / ( rho_list(n,k) * Rvap * tem_list(n,k) )

                ! Separation
                q_list(n,k,I_QV) = qsat_list(n,k)
                q_list(n,k,I_QC) = qsum_list(n,k)-qsat_list(n,k)

                CVtot(n,k) = qd_list(n,k) * CVdry
             enddo

             do nq  = NQW_STR, NQW_END
             do vec = 1, veclen
                n = blk + vec - 1

                CVtot(n,k) = CVtot(n,k) + q_list(n,k,nq) * CVW(nq)
             enddo
             enddo

             !OCL loop_fission_target(ls)
             do vec = 1, veclen
                n = blk + vec - 1

                Emoist_new = tem_list(n,k) * CVtot(n,k) + q_list(n,k,I_QV) * LHV

                ! dX/dT
                dqsat_dT   = ( LovR_liq / tem_list(n,k)**2 + CVovR_liq / tem_list(n,k) ) * qsat_list(n,k)

                dCVtot_dT  = dqsat_dT * ( CVW(I_QV) - CVW(I_QC) )

                dEmoist_dT = tem_list(n,k) * dCVtot_dT &
                           + CVtot(n,k)                &
                           + dqsat_dT    * LHV

                dtemp      = ( Emoist_new - Emoist_list(n,k) ) / dEmoist_dT

                tem_list(n,k) = tem_list(n,k) - dtemp

                exe_flag(n,k) = 0.5_RP + sign(0.5_RP, abs(dtemp)-dtemp_criteria)
             enddo

             if ( sum(exe_flag(blk:blk+veclen-1,k)) <= 0.0 ) then
                exit iteration
             endif

          enddo iteration ! iteration

       enddo ! outer block

       if ( sum(exe_flag(1:nmax(k),k)) > 0.0 ) then
          write(*,*) 'xxx [satadjust_liq] not converged! Stop.', k, int(sum(exe_flag(1:nmax(k),k)))
          call ADM_proc_stop
       endif

       do n = 1, nmax(k)
          ij = idx_list(n,k)

          tem(ij,k)      = tem_list(n,k)
          q  (ij,k,I_QV) = q_list  (n,k,I_QV)
          q  (ij,k,I_QC) = q_list  (n,k,I_QC)
       enddo

    enddo
    !$omp end parallel do

    call PROF_rapend  ('_____SatAdj_CORE',3)

    return
  end subroutine satadjust_liq

  !-----------------------------------------------------------------------------
  subroutine SATURATION_enthalpy( &
       ijdim,       &
       kdim,        &
       tem,         &
       ent,         &
       pre,         &
       qw,          &
       adiabat_type )
    use mod_cnst, only: &
       CNST_Rdry,   &
       CNST_Rvap,   &
       CNST_CPdry,  &
       CNST_CPvap,  &
       CNST_LHV0,   &
       CNST_PSAT0,  &
       CNST_EPSvap, &
       CNST_PRE00,  &
       CNST_TEM00
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV,              &
       I_QC
    use mod_thrmdyn, only: &
       THRMDYN_ent
    implicit none

    integer,          intent(in)    :: ijdim
    integer,          intent(in)    :: kdim
    real(DP),         intent(inout) :: tem(ijdim,kdim)
    real(DP),         intent(in)    :: ent(ijdim,kdim)
    real(DP),         intent(in)    :: pre(ijdim,kdim)
    real(DP),         intent(in)    :: qw (ijdim,kdim)
    character(len=*), intent(in)    :: adiabat_type
    !                                = 'RMA' reversible: qv = qv*, qc = qw - qv*
    !                                = 'PMA' pseudo    : if qw > qv* : qv = qv*, qc=0
    !                                = 'DA'  qc = 0

    real(DP) :: qd     (ijdim,kdim)
    real(DP) :: qsat   (ijdim,kdim)
    real(DP) :: q      (ijdim,kdim,nqmax)
    real(DP) :: pred   (ijdim,kdim)
    real(DP) :: prev   (ijdim,kdim)
    real(DP) :: tem_old(ijdim,kdim)
    real(DP) :: ents   (ijdim,kdim)
    real(DP) :: dents  (ijdim,kdim)
    real(DP) :: dtem   (ijdim,kdim)
    logical  :: flag   (ijdim,kdim)

    real(DP) :: PREMIN     = 1.E-10_DP
    real(DP) :: DENTS_FACT = 1.0_DP
    real(DP) :: DENTS_MIN  = 1.E-10_DP
    real(DP) :: TEMMIN     = 1.0_DP

    integer, parameter :: itelim = 100

    integer :: ite
    integer :: ij, k
    !---------------------------------------------------------------------------

    flag(:,:) = .false.

    q(:,:,:) = 0.0_RP

    if ( adiabat_type == 'DA' ) then

       tem (:,:)      = exp( ( ent(:,:) + CNST_Rdry * log( pre(:,:) / CNST_PRE00 ) ) / CNST_CPdry ) * CNST_TEM00
       qd  (:,:)      = 1.0_RP-qw(:,:)
       q   (:,:,I_QV) =        qw(:,:)
       pred(:,:)      = pre(:,:) * CNST_EPSvap * qd(:,:) / ( CNST_EPSvap * qd(:,:) + q(:,:,I_QV) )
       pred(:,:)      = max( pred(:,:), PREMIN )
       prev(:,:)      = pre(:,:) * q(:,:,I_QV)            / ( CNST_EPSvap * qd(:,:) + q(:,:,I_QV) )
       prev(:,:)      = max( prev(:,:), PREMIN )

       tem(:,:) = exp( ( ent(:,:) &
                       + qd(:,:)      * ( CNST_Rdry * log( pred(:,:) / CNST_PRE00 )                              )&
                       + q (:,:,I_QV) * ( CNST_Rvap * log( prev(:,:) / CNST_PSAT0 ) - CNST_LHV0 / CNST_TEM00 ) ) &
                     / ( qd(:,:) * CNST_CPdry + q(:,:,I_QV) * CNST_CPvap ) ) * CNST_TEM00

    elseif( adiabat_type == 'RMA' .OR. adiabat_type == 'PMA' ) then

       tem(:,:) = exp( ( ent(:,:) + CNST_Rdry * log( pre(:,:) / CNST_PRE00 ) ) / CNST_CPdry ) * CNST_TEM00

       do ite = 1, itelim

          tem_old(:,:) = tem(:,:)

          call SATURATION_qsat_liq_2D_DP( ijdim, kdim, tem(:,:), pre(:,:), qsat(:,:) )

          q(:,:,I_QV) = min ( qw(:,:), qsat(:,:) )

          if    ( adiabat_type == 'RMA' ) then
             q (:,:,I_QC) = qw(:,:) - q(:,:,I_QV)
             qd(:,:)      = 1.0_RP - qw(:,:)
          elseif( adiabat_type == 'PMA' ) then
             q (:,:,I_QC) = 0.0_RP
             qd(:,:)      = 1.0_RP - q(:,:,I_QV)
          endif

          call THRMDYN_ent( ijdim, & ! [IN]
                            kdim,  & ! [IN]
                            ents,  & ! [OUT]
                            tem,   & ! [IN]
                            pre,   & ! [IN]
                            q,     & ! [IN]
                            qd     ) ! [IN]

          do k  = 1, kdim
          do ij = 1, ijdim
             dents(:,:) = ( CNST_CPdry * qd(:,:) + CNST_CPvap * q(:,:,I_QV) ) / tem(:,:)    &
                        + (                   CNST_LHV0  * q(:,:,I_QV) ) / tem(:,:)**2 &
                        * ( CNST_LHV0 / ( tem(:,:) * CNST_Rvap ) - 1.0_RP )
             dents(:,:) = max ( dents(:,:) * DENTS_FACT, DENTS_MIN )
          enddo
          enddo

          do k  = 1, kdim
          do ij = 1, ijdim
             if ( .NOT. flag(ij,k) ) then
                tem (ij,k) = tem(ij,k) - ( ents(ij,k) - ent(ij,k) ) / dents(ij,k)
                tem (ij,k) = max ( tem(ij,k), TEMMIN )
                dtem(ij,k) = abs( tem(ij,k) - tem_old(ij,k) )
             endif

             if ( dtem(ij,k) < DTEM_EPS0 ) then
                dtem(ij,k) = 0.0_RP
                flag(ij,k) = .true.
             endif
          enddo
          enddo

          if ( maxval( dtem(:,:) ) < DTEM_EPS0 ) exit

       enddo ! iteration
    else
       write(IO_FID_LOG,*) '### SATURATION_enthalpy: invalid adiabat_type=', adiabat_type
    endif

    return
  end subroutine SATURATION_enthalpy

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  subroutine moist_dqsw_dtem_rho_SP( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem    (ijdim,kdim)
    real(SP), intent(in)  :: rho    (ijdim,kdim)
    real(SP), intent(out) :: dqsdtem(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(SP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CNST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CNST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_SP )

    return
  end subroutine moist_dqsw_dtem_rho_SP

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  subroutine moist_dqsw_dtem_rho_DP( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem    (ijdim,kdim)
    real(DP), intent(in)  :: rho    (ijdim,kdim)
    real(DP), intent(out) :: dqsdtem(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(DP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CNST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CNST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_DP )

    return
  end subroutine moist_dqsw_dtem_rho_DP

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  subroutine moist_dqsi_dtem_rho_SP( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem    (ijdim,kdim)
    real(SP), intent(in)  :: rho    (ijdim,kdim)
    real(SP), intent(out) :: dqsdtem(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(SP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CNST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CNST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_SP )

    return
  end subroutine moist_dqsi_dtem_rho_SP

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  subroutine moist_dqsi_dtem_rho_DP( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
    use mod_cnst, only: &
       CNST_Rvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem    (ijdim,kdim)
    real(DP), intent(in)  :: rho    (ijdim,kdim)
    real(DP), intent(out) :: dqsdtem(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(DP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CNST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CNST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_DP )

    return
  end subroutine moist_dqsi_dtem_rho_DP

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsw_dtem_dpre_SP( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
    use mod_cnst, only: &
       CNST_EPSvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem    (ijdim,kdim)
    real(SP), intent(in)  :: pre    (ijdim,kdim)
    real(SP), intent(out) :: dqsdtem(ijdim,kdim)
    real(SP), intent(out) :: dqsdpre(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(SP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(SP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CNST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_SP-CNST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CNST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CNST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsw_dtem_dpre_SP

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsw_dtem_dpre_DP( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
    use mod_cnst, only: &
       CNST_EPSvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem    (ijdim,kdim)
    real(DP), intent(in)  :: pre    (ijdim,kdim)
    real(DP), intent(out) :: dqsdtem(ijdim,kdim)
    real(DP), intent(out) :: dqsdpre(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(DP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(DP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CNST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_DP-CNST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CNST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CNST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsw_dtem_dpre_DP

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsi_dtem_dpre_SP( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
    use mod_cnst, only: &
       CNST_EPSvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem    (ijdim,kdim)
    real(SP), intent(in)  :: pre    (ijdim,kdim)
    real(SP), intent(out) :: dqsdtem(ijdim,kdim)
    real(SP), intent(out) :: dqsdpre(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(SP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(SP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CNST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_SP-CNST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CNST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CNST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsi_dtem_dpre_SP

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsi_dtem_dpre_DP( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
    use mod_cnst, only: &
       CNST_EPSvap, &
       CNST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem    (ijdim,kdim)
    real(DP), intent(in)  :: pre    (ijdim,kdim)
    real(DP), intent(out) :: dqsdtem(ijdim,kdim)
    real(DP), intent(out) :: dqsdpre(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(DP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(DP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CNST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_DP-CNST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CNST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CNST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsi_dtem_dpre_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq2_0D( &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    real(RP), intent(in)  :: tem
    real(RP), intent(in)  :: rho
    real(RP), intent(in)  :: qv
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: Rvap
    !---------------------------------------------------------------------------

    Rvap = CNST_Rvap

    call SATURATION_psat_liq_0D( tem, psat )

    qsat = psat / ( rho*(1.0_RP-qv)*tem*Rvap + psat )

    return
  end subroutine SATURATION_qsat_liq2_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq2_1D_SP( &
       ijdim, &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: tem (ijdim)
    real(SP), intent(in)  :: rho (ijdim)
    real(SP), intent(in)  :: qv  (ijdim)
    real(SP), intent(out) :: qsat(ijdim)

    real(SP) :: psat(ijdim)
    real(SP) :: Rvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_liq_1D_SP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,rho,tem,psat,qv,Rvap)
    do ij = 1, ijdim
       qsat(ij) = psat(ij) / ( rho(ij)*(1.0_SP-qv(ij))*tem(ij)*Rvap + psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq2_1D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq2_1D_DP( &
       ijdim, &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: tem (ijdim)
    real(DP), intent(in)  :: rho (ijdim)
    real(DP), intent(in)  :: qv  (ijdim)
    real(DP), intent(out) :: qsat(ijdim)

    real(DP) :: psat(ijdim)
    real(DP) :: Rvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_liq_1D_DP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,rho,tem,psat,qv,Rvap)
    do ij = 1, ijdim
       qsat(ij) = psat(ij) / ( rho(ij)*(1.0_DP-qv(ij))*tem(ij)*Rvap + psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq2_1D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq2_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem (ijdim,kdim)
    real(SP), intent(in)  :: rho (ijdim,kdim)
    real(SP), intent(in)  :: qv  (ijdim,kdim)
    real(SP), intent(out) :: qsat(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim)
    real(SP) :: Rvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcreate(qsat,psat) pcopyin(tem,rho,qv)

    Rvap = CNST_Rvap

    call SATURATION_psat_liq_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,rho,tem,psat,qv,Rvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = psat(ij,k) / ( rho(ij,k)*(1.0_SP-qv(ij,k))*tem(ij,k)*Rvap + psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq2_2D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq2_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem (ijdim,kdim)
    real(DP), intent(in)  :: rho (ijdim,kdim)
    real(DP), intent(in)  :: qv  (ijdim,kdim)
    real(DP), intent(out) :: qsat(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim)
    real(DP) :: Rvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcreate(qsat,psat) pcopyin(tem,rho,qv)

    Rvap = CNST_Rvap

    call SATURATION_psat_liq_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,rho,tem,psat,qv,Rvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = psat(ij,k) / ( rho(ij,k)*(1.0_DP-qv(ij,k))*tem(ij,k)*Rvap + psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_liq2_2D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice2_0D( &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    real(RP), intent(in)  :: tem
    real(RP), intent(in)  :: rho
    real(RP), intent(in)  :: qv
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: Rvap
    !---------------------------------------------------------------------------

    Rvap = CNST_Rvap

    call SATURATION_psat_ice_0D( tem, psat )

    qsat = psat / ( rho*(1.0_RP-qv)*tem*Rvap + psat )

    return
  end subroutine SATURATION_qsat_ice2_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice2_1D_SP( &
       ijdim, &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: tem (ijdim)
    real(SP), intent(in)  :: rho (ijdim)
    real(SP), intent(in)  :: qv  (ijdim)
    real(SP), intent(out) :: qsat(ijdim)

    real(SP) :: psat(ijdim)
    real(SP) :: Rvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_ice_1D_SP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,rho,tem,psat,qv,Rvap)
    do ij = 1, ijdim
       qsat(ij) = psat(ij) / ( rho(ij)*(1.0_SP-qv(ij))*tem(ij)*Rvap + psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice2_1D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice2_1D_DP( &
       ijdim, &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: tem (ijdim)
    real(DP), intent(in)  :: rho (ijdim)
    real(DP), intent(in)  :: qv  (ijdim)
    real(DP), intent(out) :: qsat(ijdim)

    real(DP) :: psat(ijdim)
    real(DP) :: Rvap

    integer :: ij
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_ice_1D_DP( ijdim, tem(:), psat(:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,rho,tem,psat,qv,Rvap)
    do ij = 1, ijdim
       qsat(ij) = psat(ij) / ( rho(ij)*(1.0_DP-qv(ij))*tem(ij)*Rvap + psat(ij) )
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice2_1D_DP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice2_2D_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem (ijdim,kdim)
    real(SP), intent(in)  :: rho (ijdim,kdim)
    real(SP), intent(in)  :: qv  (ijdim,kdim)
    real(SP), intent(out) :: qsat(ijdim,kdim)

    real(SP) :: psat(ijdim,kdim)
    real(SP) :: Rvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_ice_2D_SP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,rho,tem,psat,qv,Rvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = psat(ij,k) / ( rho(ij,k)*(1.0_SP-qv(ij,k))*tem(ij,k)*Rvap + psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice2_2D_SP

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice2_2D_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       rho,   &
       qv,    &
       qsat   )
    use mod_cnst, only: &
       CNST_Rvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem (ijdim,kdim)
    real(DP), intent(in)  :: rho (ijdim,kdim)
    real(DP), intent(in)  :: qv  (ijdim,kdim)
    real(DP), intent(out) :: qsat(ijdim,kdim)

    real(DP) :: psat(ijdim,kdim)
    real(DP) :: Rvap

    integer :: ij, k
    !---------------------------------------------------------------------------
    !$acc data pcopyin(tem,rho,qv) pcreate(psat) if (acc_is_present(qsat))

    Rvap = CNST_Rvap

    call SATURATION_psat_ice_2D_DP( ijdim, kdim, tem(:,:), psat(:,:) )

    !$acc kernels async(1) pcopy(qsat) pcopyin(psat,tem,rho,qv)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,rho,tem,psat,qv,Rvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = psat(ij,k) / ( rho(ij,k)*(1.0_DP-qv(ij,k))*tem(ij,k)*Rvap + psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    !$acc wait(1)

    !$acc end data
    return
  end subroutine SATURATION_qsat_ice2_2D_DP

end module mod_satadjust
