!-------------------------------------------------------------------------------
!> Module boundary conditions
!!
!! @par Description
!!          This module provides the subroutines for boundary conditions
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_bndcnd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
#ifdef _OPENACC
  use openacc
#endif
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: BNDCND_setup
  public :: BNDCND_get_thermo_flags
  public :: BNDCND_all
  public :: BNDCND_thermo
  public :: BNDCND_rhovxvyvz
  public :: BNDCND_rhow
  public :: BNDCND_pre_sfc

  interface BNDCND_all
     module procedure BNDCND_all_SP
     module procedure BNDCND_all_DP
  end interface BNDCND_all

  interface BNDCND_thermo
     module procedure BNDCND_thermo_SP
     module procedure BNDCND_thermo_DP
  end interface BNDCND_thermo

  interface BNDCND_rhovxvyvz
     module procedure BNDCND_rhovxvyvz_SP
     module procedure BNDCND_rhovxvyvz_DP
  end interface BNDCND_rhovxvyvz

  interface BNDCND_rhow
     module procedure BNDCND_rhow_SP
     module procedure BNDCND_rhow_DP
  end interface BNDCND_rhow

  interface BNDCND_pre_sfc
     module procedure BNDCND_pre_sfc_SP
     module procedure BNDCND_pre_sfc_SPDP
     module procedure BNDCND_pre_sfc_DP
  end interface BNDCND_pre_sfc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: acc_wait_on = .true.
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !--- Vertical boundary condition for temperature at the top
  character(len=H_SHORT), private :: BND_TYPE_T_TOP    != 'TEM' : tem(kmax+1) = tem(kmax) (default)
                                                       != 'EPL' : lagrange extrapolation
                                                       != 'ASIS': leave as it is

  !--- Vertical boundary condition for temperature at the ground
  character(len=H_SHORT), private :: BND_TYPE_T_BOTTOM != 'TEM' : tem(kmin-1) = tem(kmin) (default)
                                                       != 'EPL' : lagrange extrapolation
                                                       != 'ASIS': leave as it is

  !--- Vertical boundary condition for momentum at the top
  character(len=H_SHORT), private :: BND_TYPE_M_TOP    != 'RIGID' : rigid surface
                                                       != 'FREE'  : free surface (default)

  !--- Vertical boundary condition for momentum at the ground
  character(len=H_SHORT), private :: BND_TYPE_M_BOTTOM != 'RIGID' : rigid surface (default)
                                                       != 'FREE'  : free surface

  logical, private :: is_top_tem   = .false.
  logical, private :: is_top_epl   = .false.
  logical, private :: is_btm_tem   = .false.
  logical, private :: is_btm_epl   = .false.
  logical, private :: is_top_rigid = .false.
  logical, private :: is_top_free  = .false.
  logical, private :: is_btm_rigid = .false.
  logical, private :: is_btm_free  = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine BNDCND_setup
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    namelist / BNDCNDPARAM / &
       BND_TYPE_T_TOP,    &
       BND_TYPE_T_BOTTOM, &
       BND_TYPE_M_TOP,    &
       BND_TYPE_M_BOTTOM

    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- set default
    BND_TYPE_T_TOP    = 'TEM'
    BND_TYPE_T_BOTTOM = 'TEM'
    BND_TYPE_M_TOP    = 'FREE'
    BND_TYPE_M_BOTTOM = 'RIGID'

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[bndcnd]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=BNDCNDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** BNDCNDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist BNDCNDPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=BNDCNDPARAM)

    write(IO_FID_LOG,*)

    if    ( BND_TYPE_T_TOP == 'TEM' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, top   ) : equal to uppermost atmosphere'
       is_top_tem = .true.
    elseif( BND_TYPE_T_TOP == 'EPL' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, top   ) : lagrange extrapolation'
       is_top_epl = .true.
    elseif( BND_TYPE_T_TOP == 'ASIS' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, top   ) : leave as it is'
       is_top_epl = .true.
    else
       write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_T_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_T_BOTTOM == 'TEM' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, bottom) : equal to lowermost atmosphere'
       is_btm_tem = .true.
    elseif( BND_TYPE_T_BOTTOM == 'EPL' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, bottom) : lagrange extrapolation'
       is_btm_epl = .true.
    elseif( BND_TYPE_T_BOTTOM == 'ASIS' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (temperature, top   ) : leave as it is'
    else
       write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_T_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_M_TOP == 'RIGID' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (momentum,    top   ) : rigid'
       is_top_rigid = .true.
    elseif( BND_TYPE_M_TOP == 'FREE' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (momentum,    top   ) : free'
       is_top_free  = .true.
    else
       write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_M_TOP. STOP.'
       call ADM_proc_stop
    endif

    if    ( BND_TYPE_M_BOTTOM == 'RIGID' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (momentum,    bottom) : rigid'
       is_btm_rigid = .true.
    elseif( BND_TYPE_M_BOTTOM == 'FREE' ) then
       write(IO_FID_LOG,*) '--- Boundary setting type (momentum,    bottom) : free'
       is_btm_free  = .true.
    else
       write(IO_FID_LOG,*) 'xxx Invalid BND_TYPE_M_BOTTOM. STOP.'
       call ADM_proc_stop
    endif

    return
  end subroutine BNDCND_setup

  !-----------------------------------------------------------------------------
  subroutine BNDCND_get_thermo_flags( top_tem, top_epl, btm_tem, btm_epl )
    implicit none
    logical, intent(out) :: top_tem, top_epl, btm_tem, btm_epl
    top_tem = is_top_tem
    top_epl = is_top_epl
    btm_tem = is_btm_tem
    btm_epl = is_btm_epl
  end subroutine BNDCND_get_thermo_flags

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for all prognostic variables.
  subroutine BNDCND_all_SP( &
       ijdim,      &
       kdim,       &
       ldim,       &
       rho,        &
       vx,         &
       vy,         &
       vz,         &
       w,          &
       ein,        &
       tem,        &
       pre,        &
       rhog,       &
       rhogvx,     &
       rhogvy,     &
       rhogvz,     &
       rhogw,      &
       rhoge,      &
       gsqrtgam2,  &
       phi,        &
       C2Wfact,    &
       C2Wfact_Gz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    use mod_cnst, only: &
       CNST_CVdry
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(SP), intent(inout) :: rho       (ijdim,kdim,ldim) ! density
    real(SP), intent(inout) :: vx        (ijdim,kdim,ldim) ! horizontal wind (x)
    real(SP), intent(inout) :: vy        (ijdim,kdim,ldim) ! horizontal wind (y)
    real(SP), intent(inout) :: vz        (ijdim,kdim,ldim) ! horizontal wind (z)
    real(SP), intent(inout) :: w         (ijdim,kdim,ldim) ! vertical   wind
    real(SP), intent(inout) :: ein       (ijdim,kdim,ldim) ! internal energy
    real(SP), intent(inout) :: tem       (ijdim,kdim,ldim) ! temperature
    real(SP), intent(inout) :: pre       (ijdim,kdim,ldim) ! pressure
    real(SP), intent(inout) :: rhog      (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvx    (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvy    (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvz    (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogw     (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhoge     (ijdim,kdim,ldim)
    real(SP), intent(in)    :: gsqrtgam2 (ijdim,kdim,ldim)
    real(SP), intent(in)    :: phi       (ijdim,kdim,ldim)   ! geopotential
    real(SP), intent(in)    :: C2Wfact   (ijdim,kdim,2,ldim)
    real(SP), intent(in)    :: C2Wfact_Gz(ijdim,kdim,6,ldim)

    integer  :: kmin, kmax
    real(SP) :: CVdry

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(gsqrtgam2,phi,C2Wfact,C2Wfact_Gz) if(acc_is_present(rho))

    kmin = ADM_kmin
    kmax = ADM_kmax

    CVdry = CNST_CVdry

    !--- Thermodynamical variables ( rho, ein, tem, pre, rhog, rhoge ), q = 0 at boundary

    do l = 1, ldim
       call BNDCND_thermo( ijdim,      & ! [IN]
                           tem(:,:,l), & ! [INOUT]
                           rho(:,:,l), & ! [INOUT]
                           pre(:,:,l), & ! [INOUT]
                           phi(:,:,l)  ) ! [IN]
    enddo

    !$acc kernels async(1) pcopy(rhog,ein,rhoge) pcopyin(rho,tem,gsqrtgam2) if(acc_is_present(rhog))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,rhog,rhoge,ein,rho,tem,gsqrtgam2,CVdry)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          rhog (ij,kmax+1,l) = rho(ij,kmax+1,l) * gsqrtgam2(ij,kmax+1,l)
          rhog (ij,kmin-1,l) = rho(ij,kmin-1,l) * gsqrtgam2(ij,kmin-1,l)

          ein  (ij,kmax+1,l) = CVdry * tem(ij,kmax+1,l)
          ein  (ij,kmin-1,l) = CVdry * tem(ij,kmin-1,l)

          rhoge(ij,kmax+1,l) = rhog(ij,kmax+1,l) * ein(ij,kmax+1,l)
          rhoge(ij,kmin-1,l) = rhog(ij,kmin-1,l) * ein(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogvx, rhogvy, rhogvz, vx, vy, vz )

    call BNDCND_rhovxvyvz( ijdim,         & ! [IN]
                           kdim,          & ! [IN]
                           ldim,          & ! [IN]
                           rhog  (:,:,:), & ! [IN]
                           rhogvx(:,:,:), & ! [INOUT]
                           rhogvy(:,:,:), & ! [INOUT]
                           rhogvz(:,:,:)  ) ! [INOUT]

    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(rhogvx,rhogvy,rhogvz,rhog) if(acc_is_present(vx))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,vx,vy,vz,rhog,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          vx(ij,kmax+1,l) = rhogvx(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vy(ij,kmax+1,l) = rhogvy(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vz(ij,kmax+1,l) = rhogvz(ij,kmax+1,l) / rhog(ij,kmax+1,l)

          vx(ij,kmin-1,l) = rhogvx(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vy(ij,kmin-1,l) = rhogvy(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vz(ij,kmin-1,l) = rhogvz(ij,kmin-1,l) / rhog(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogw, w )
    do l = 1, ldim
       call BNDCND_rhow( ijdim,              & ! [IN]
                         rhogvx    (:,:,l),  & ! [IN]
                         rhogvy    (:,:,l),  & ! [IN]
                         rhogvz    (:,:,l),  & ! [IN]
                         rhogw     (:,:,l),  & ! [INOUT]
                         C2Wfact_Gz(:,:,:,l) ) ! [IN]
    enddo

    !$acc kernels async(1) pcopy(w) pcopyin(rhog,rhogw,C2Wfact) if(acc_is_present(w))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,w,rhogw,rhog,C2Wfact)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          w(ij,kmax+1,l) = rhogw(ij,kmax+1,l) / ( C2Wfact(ij,kmax+1,1,l) * rhog(ij,kmax+1,l) &
                                                + C2Wfact(ij,kmax+1,2,l) * rhog(ij,kmax  ,l) )
          w(ij,kmin  ,l) = rhogw(ij,kmin  ,l) / ( C2Wfact(ij,kmin  ,1,l) * rhog(ij,kmin  ,l) &
                                                + C2Wfact(ij,kmin  ,2,l) * rhog(ij,kmin-1,l) )
          w(ij,kmin-1,l) = 0.0_SP
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_all_SP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for all prognostic variables.
  subroutine BNDCND_all_DP( &
       ijdim,      &
       kdim,       &
       ldim,       &
       rho,        &
       vx,         &
       vy,         &
       vz,         &
       w,          &
       ein,        &
       tem,        &
       pre,        &
       rhog,       &
       rhogvx,     &
       rhogvy,     &
       rhogvz,     &
       rhogw,      &
       rhoge,      &
       gsqrtgam2,  &
       phi,        &
       C2Wfact,    &
       C2Wfact_Gz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    use mod_cnst, only: &
       CNST_CVdry
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(DP), intent(inout) :: rho       (ijdim,kdim,ldim) ! density
    real(DP), intent(inout) :: vx        (ijdim,kdim,ldim) ! horizontal wind (x)
    real(DP), intent(inout) :: vy        (ijdim,kdim,ldim) ! horizontal wind (y)
    real(DP), intent(inout) :: vz        (ijdim,kdim,ldim) ! horizontal wind (z)
    real(DP), intent(inout) :: w         (ijdim,kdim,ldim) ! vertical   wind
    real(DP), intent(inout) :: ein       (ijdim,kdim,ldim) ! internal energy
    real(DP), intent(inout) :: tem       (ijdim,kdim,ldim) ! temperature
    real(DP), intent(inout) :: pre       (ijdim,kdim,ldim) ! pressure
    real(DP), intent(inout) :: rhog      (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvx    (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvy    (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvz    (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogw     (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhoge     (ijdim,kdim,ldim)
    real(DP), intent(in)    :: gsqrtgam2 (ijdim,kdim,ldim)
    real(DP), intent(in)    :: phi       (ijdim,kdim,ldim)   ! geopotential
    real(DP), intent(in)    :: C2Wfact   (ijdim,kdim,2,ldim)
    real(DP), intent(in)    :: C2Wfact_Gz(ijdim,kdim,6,ldim)

    integer  :: kmin, kmax
    real(DP) :: CVdry

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(gsqrtgam2,phi,C2Wfact,C2Wfact_Gz) if(acc_is_present(rho))

    kmin = ADM_kmin
    kmax = ADM_kmax

    CVdry = CNST_CVdry

    !--- Thermodynamical variables ( rho, ein, tem, pre, rhog, rhoge ), q = 0 at boundary

    do l = 1, ldim
       call BNDCND_thermo( ijdim,      & ! [IN]
                           tem(:,:,l), & ! [INOUT]
                           rho(:,:,l), & ! [INOUT]
                           pre(:,:,l), & ! [INOUT]
                           phi(:,:,l)  ) ! [IN]
    enddo

    !$acc kernels async(1) pcopy(rhog,ein,rhoge) pcopyin(rho,tem,gsqrtgam2) if(acc_is_present(rhog))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,rhog,rhoge,ein,rho,tem,gsqrtgam2,CVdry)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          rhog (ij,kmax+1,l) = rho(ij,kmax+1,l) * gsqrtgam2(ij,kmax+1,l)
          rhog (ij,kmin-1,l) = rho(ij,kmin-1,l) * gsqrtgam2(ij,kmin-1,l)

          ein  (ij,kmax+1,l) = CVdry * tem(ij,kmax+1,l)
          ein  (ij,kmin-1,l) = CVdry * tem(ij,kmin-1,l)

          rhoge(ij,kmax+1,l) = rhog(ij,kmax+1,l) * ein(ij,kmax+1,l)
          rhoge(ij,kmin-1,l) = rhog(ij,kmin-1,l) * ein(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogvx, rhogvy, rhogvz, vx, vy, vz )

    call BNDCND_rhovxvyvz( ijdim,         & ! [IN]
                           kdim,          & ! [IN]
                           ldim,          & ! [IN]
                           rhog  (:,:,:), & ! [IN]
                           rhogvx(:,:,:), & ! [INOUT]
                           rhogvy(:,:,:), & ! [INOUT]
                           rhogvz(:,:,:)  ) ! [INOUT]

    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(rhogvx,rhogvy,rhogvz,rhog) if(acc_is_present(vx))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,vx,vy,vz,rhog,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          vx(ij,kmax+1,l) = rhogvx(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vy(ij,kmax+1,l) = rhogvy(ij,kmax+1,l) / rhog(ij,kmax+1,l)
          vz(ij,kmax+1,l) = rhogvz(ij,kmax+1,l) / rhog(ij,kmax+1,l)

          vx(ij,kmin-1,l) = rhogvx(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vy(ij,kmin-1,l) = rhogvy(ij,kmin-1,l) / rhog(ij,kmin-1,l)
          vz(ij,kmin-1,l) = rhogvz(ij,kmin-1,l) / rhog(ij,kmin-1,l)
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels

    !--- Momentum ( rhogw, w )
    do l = 1, ldim
       call BNDCND_rhow( ijdim,              & ! [IN]
                         rhogvx    (:,:,l),  & ! [IN]
                         rhogvy    (:,:,l),  & ! [IN]
                         rhogvz    (:,:,l),  & ! [IN]
                         rhogw     (:,:,l),  & ! [INOUT]
                         C2Wfact_Gz(:,:,:,l) ) ! [IN]
    enddo

    !$acc kernels async(1) pcopy(w) pcopyin(rhog,rhogw,C2Wfact) if(acc_is_present(w))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,w,rhogw,rhog,C2Wfact)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
!OCL XFILL
       !$omp do
       do ij = 1, ijdim
          w(ij,kmax+1,l) = rhogw(ij,kmax+1,l) / ( C2Wfact(ij,kmax+1,1,l) * rhog(ij,kmax+1,l) &
                                                + C2Wfact(ij,kmax+1,2,l) * rhog(ij,kmax  ,l) )
          w(ij,kmin  ,l) = rhogw(ij,kmin  ,l) / ( C2Wfact(ij,kmin  ,1,l) * rhog(ij,kmin  ,l) &
                                                + C2Wfact(ij,kmin  ,2,l) * rhog(ij,kmin-1,l) )
          w(ij,kmin-1,l) = 0.0_DP
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_all_DP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for thermodynamical variables
  subroutine BNDCND_thermo_SP( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    use mod_cnst, only: &
       CNST_GRAV, &
       CNST_Rdry
    implicit none

    integer,  intent(in)    :: ijdim           ! number of horizontal grid
    real(SP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(SP), intent(inout) :: rho(ijdim,kdim) ! density
    real(SP), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(SP), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer  :: kmin, kmax
    real(SP) :: GRAV, Rdry

    integer  :: ij

    real(SP) :: z,z1,z2,z3,p1,p2,p3
    real(SP) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(phi) if(acc_is_present(rho))

    kmin = ADM_kmin
    kmax = ADM_kmax

    GRAV = CNST_GRAV
    Rdry = CNST_Rdry

    !$acc kernels async(1) pcopy(tem,pre,rho) pcopyin(phi) if(acc_is_present(rho))
    !$omp parallel do default(none),private(ij,z,z1,z2,z3), &
    !$omp shared(ijdim,kmax,kmin,is_top_tem,is_top_epl,is_btm_tem,is_btm_epl, &
    !$omp        rho,pre,tem,phi,Rdry,GRAV)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do ij = 1, ijdim

       if    ( is_top_tem ) then

          tem(ij,kmax+1) = tem(ij,kmax) ! dT/dz = 0

       elseif( is_top_epl ) then
          z  = phi(ij,kmax+1) / GRAV
          z1 = phi(ij,kmax  ) / GRAV
          z2 = phi(ij,kmax-1) / GRAV
          z3 = phi(ij,kmax-2) / GRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )
       endif

       if   ( is_btm_tem ) then

          tem(ij,kmin-1) = tem(ij,kmin) ! dT/dz = 0

       elseif( is_btm_epl ) then
          z1 = phi(ij,kmin+2) / GRAV
          z2 = phi(ij,kmin+1) / GRAV
          z3 = phi(ij,kmin  ) / GRAV
          z  = phi(ij,kmin-1) / GRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       endif

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density ( equation of state )
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( Rdry * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( Rdry * tem(ij,kmin-1) )

    enddo
    !$omp end parallel do
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_thermo_SP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for thermodynamical variables
  subroutine BNDCND_thermo_DP( &
       ijdim, &
       tem,   &
       rho,   &
       pre,   &
       phi    )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    use mod_cnst, only: &
       CNST_GRAV, &
       CNST_Rdry
    implicit none

    integer,  intent(in)    :: ijdim           ! number of horizontal grid
    real(DP), intent(inout) :: tem(ijdim,kdim) ! temperature
    real(DP), intent(inout) :: rho(ijdim,kdim) ! density
    real(DP), intent(inout) :: pre(ijdim,kdim) ! pressure
    real(DP), intent(in)    :: phi(ijdim,kdim) ! geopotential

    integer  :: kmin, kmax
    real(DP) :: GRAV, Rdry

    integer  :: ij

    real(DP) :: z,z1,z2,z3,p1,p2,p3
    real(DP) :: lag_intpl
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3))/((z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3))/((z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2))/((z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(phi) if(acc_is_present(rho))

    kmin = ADM_kmin
    kmax = ADM_kmax

    GRAV = CNST_GRAV
    Rdry = CNST_Rdry

    !$acc kernels async(1) pcopy(tem,pre,rho) pcopyin(phi) if(acc_is_present(rho))
    !$omp parallel do default(none),private(ij,z,z1,z2,z3), &
    !$omp shared(ijdim,kmax,kmin,is_top_tem,is_top_epl,is_btm_tem,is_btm_epl, &
    !$omp        rho,pre,tem,phi,Rdry,GRAV)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do ij = 1, ijdim

       if    ( is_top_tem ) then

          tem(ij,kmax+1) = tem(ij,kmax) ! dT/dz = 0

       elseif( is_top_epl ) then
          z  = phi(ij,kmax+1) / GRAV
          z1 = phi(ij,kmax  ) / GRAV
          z2 = phi(ij,kmax-1) / GRAV
          z3 = phi(ij,kmax-2) / GRAV

          tem(ij,kmax+1) = lag_intpl( z ,                 &
                                      z1, tem(ij,kmax  ), &
                                      z2, tem(ij,kmax-1), &
                                      z3, tem(ij,kmax-2)  )
       endif

       if   ( is_btm_tem ) then

          tem(ij,kmin-1) = tem(ij,kmin) ! dT/dz = 0

       elseif( is_btm_epl ) then
          z1 = phi(ij,kmin+2) / GRAV
          z2 = phi(ij,kmin+1) / GRAV
          z3 = phi(ij,kmin  ) / GRAV
          z  = phi(ij,kmin-1) / GRAV

          tem(ij,kmin-1) = lag_intpl( z,                  &
                                      z1, tem(ij,kmin+2), &
                                      z2, tem(ij,kmin+1), &
                                      z3, tem(ij,kmin  )  )
       endif

       !--- set the boundary of pressure ( hydrostatic balance )
       pre(ij,kmax+1) = pre(ij,kmax-1) - rho(ij,kmax) * ( phi(ij,kmax+1) - phi(ij,kmax-1) )
       pre(ij,kmin-1) = pre(ij,kmin+1) - rho(ij,kmin) * ( phi(ij,kmin-1) - phi(ij,kmin+1) )

       !--- set the boundary of density ( equation of state )
       rho(ij,kmax+1) = pre(ij,kmax+1) / ( Rdry * tem(ij,kmax+1) )
       rho(ij,kmin-1) = pre(ij,kmin-1) / ( Rdry * tem(ij,kmin-1) )

    enddo
    !$omp end parallel do
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_thermo_DP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for horizontal momentum
  subroutine BNDCND_rhovxvyvz_SP( &
       ijdim,  &
       kdim,   &
       ldim,   &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(SP), intent(in)    :: rhog  (ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvx(ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvy(ijdim,kdim,ldim)
    real(SP), intent(inout) :: rhogvz(ijdim,kdim,ldim)

    integer  :: kmin, kmax

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhog) if(acc_is_present(rhogvx))

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels async(1) pcopy(rhogvx,rhogvy,rhogvz) pcopyin(rhog) if(acc_is_present(rhogvx))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhog,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
       !$omp do
       do ij = 1, ijdim

          if    ( is_top_rigid ) then
             rhogvx(ij,kmax+1,l) = -rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) = -rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) = -rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          elseif( is_top_free  ) then
             rhogvx(ij,kmax+1,l) =  rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) =  rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) =  rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          endif

          if    ( is_btm_rigid ) then
             rhogvx(ij,kmin-1,l) = -rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) = -rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) = -rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          elseif( is_btm_free  ) then
             rhogvx(ij,kmin-1,l) =  rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) =  rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) =  rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          endif

       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_rhovxvyvz_SP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for horizontal momentum
  subroutine BNDCND_rhovxvyvz_DP( &
       ijdim,  &
       kdim,   &
       ldim,   &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogvz  )
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    integer,  intent(in)    :: ldim
    real(DP), intent(in)    :: rhog  (ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvx(ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvy(ijdim,kdim,ldim)
    real(DP), intent(inout) :: rhogvz(ijdim,kdim,ldim)

    integer  :: kmin, kmax

    integer  :: ij, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhog) if(acc_is_present(rhogvx))

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels async(1) pcopy(rhogvx,rhogvy,rhogvz) pcopyin(rhog) if(acc_is_present(rhogvx))
    !$omp parallel default(none),private(ij,l), &
    !$omp shared(ijdim,ldim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhog,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do l  = 1, ldim
       !$omp do
       do ij = 1, ijdim

          if    ( is_top_rigid ) then
             rhogvx(ij,kmax+1,l) = -rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) = -rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) = -rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          elseif( is_top_free  ) then
             rhogvx(ij,kmax+1,l) =  rhogvx(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvy(ij,kmax+1,l) =  rhogvy(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
             rhogvz(ij,kmax+1,l) =  rhogvz(ij,kmax,l) / rhog(ij,kmax,l) * rhog(ij,kmax+1,l)
          endif

          if    ( is_btm_rigid ) then
             rhogvx(ij,kmin-1,l) = -rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) = -rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) = -rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          elseif( is_btm_free  ) then
             rhogvx(ij,kmin-1,l) =  rhogvx(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvy(ij,kmin-1,l) =  rhogvy(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
             rhogvz(ij,kmin-1,l) =  rhogvz(ij,kmin,l) / rhog(ij,kmin,l) * rhog(ij,kmin-1,l)
          endif

       enddo
       !$omp end do
    enddo
    !$omp end parallel
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_rhovxvyvz_DP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for vertical momentum
  subroutine BNDCND_rhow_SP( &
       ijdim,     &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact_Gz )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    real(SP), intent(in)    :: rhogvx    (ijdim,kdim)
    real(SP), intent(in)    :: rhogvy    (ijdim,kdim)
    real(SP), intent(in)    :: rhogvz    (ijdim,kdim)
    real(SP), intent(inout) :: rhogw     (ijdim,kdim)
    real(SP), intent(in)    :: C2Wfact_Gz(ijdim,kdim,6)

    integer  :: kmin, kmax

    integer  :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhogvx,rhogvy,rhogvz,C2Wfact_Gz) if(acc_is_present(rhogw))

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels async(1) pcopy(rhogw) pcopyin(rhogvx,rhogvy,rhogvz,C2Wfact_Gz) if(acc_is_present(rhogw))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhogw,C2Wfact_Gz,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do ij = 1, ijdim

       if    ( is_top_rigid ) then
          rhogw(ij,kmax+1) = 0.0_SP
       elseif( is_top_free  ) then
          rhogw(ij,kmax+1) = - ( C2Wfact_Gz(ij,kmax+1,1) * rhogvx(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,2) * rhogvx(ij,kmax  ) &
                               + C2Wfact_Gz(ij,kmax+1,3) * rhogvy(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,4) * rhogvy(ij,kmax  ) &
                               + C2Wfact_Gz(ij,kmax+1,5) * rhogvz(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,6) * rhogvz(ij,kmax  ) )
       endif

       if    ( is_btm_rigid ) then
          rhogw(ij,kmin) = 0.0_SP
       elseif( is_btm_free  ) then
          rhogw(ij,kmin) = - ( C2Wfact_Gz(ij,kmin,1) * rhogvx(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,2) * rhogvx(ij,kmin-1) &
                             + C2Wfact_Gz(ij,kmin,3) * rhogvy(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,4) * rhogvy(ij,kmin-1) &
                             + C2Wfact_Gz(ij,kmin,5) * rhogvz(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,6) * rhogvz(ij,kmin-1) )
       endif

       rhogw(ij,kmin-1) = 0.0_SP

    enddo
    !$omp end parallel do
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_rhow_SP

  !-----------------------------------------------------------------------------
  !> Boundary condition setting for vertical momentum
  subroutine BNDCND_rhow_DP( &
       ijdim,     &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact_Gz )
    use mod_adm, only: &
       kdim => ADM_kall, &
       ADM_kmin,         &
       ADM_kmax
    implicit none

    integer,  intent(in)    :: ijdim
    real(DP), intent(in)    :: rhogvx    (ijdim,kdim)
    real(DP), intent(in)    :: rhogvy    (ijdim,kdim)
    real(DP), intent(in)    :: rhogvz    (ijdim,kdim)
    real(DP), intent(inout) :: rhogw     (ijdim,kdim)
    real(DP), intent(in)    :: C2Wfact_Gz(ijdim,kdim,6)

    integer  :: kmin, kmax

    integer  :: ij
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhogvx,rhogvy,rhogvz,C2Wfact_Gz) if(acc_is_present(rhogw))

    kmin = ADM_kmin
    kmax = ADM_kmax

    !$acc kernels async(1) pcopy(rhogw) pcopyin(rhogvx,rhogvy,rhogvz,C2Wfact_Gz) if(acc_is_present(rhogw))
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmax,kmin,is_top_rigid,is_top_free,is_btm_rigid,is_btm_free, &
    !$omp        rhogw,C2Wfact_Gz,rhogvx,rhogvy,rhogvz)
!OCL PREFETCH_SEQUENTIAL(SOFT)
    do ij = 1, ijdim

       if    ( is_top_rigid ) then
          rhogw(ij,kmax+1) = 0.0_DP
       elseif( is_top_free  ) then
          rhogw(ij,kmax+1) = - ( C2Wfact_Gz(ij,kmax+1,1) * rhogvx(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,2) * rhogvx(ij,kmax  ) &
                               + C2Wfact_Gz(ij,kmax+1,3) * rhogvy(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,4) * rhogvy(ij,kmax  ) &
                               + C2Wfact_Gz(ij,kmax+1,5) * rhogvz(ij,kmax+1) &
                               + C2Wfact_Gz(ij,kmax+1,6) * rhogvz(ij,kmax  ) )
       endif

       if    ( is_btm_rigid ) then
          rhogw(ij,kmin) = 0.0_DP
       elseif( is_btm_free  ) then
          rhogw(ij,kmin) = - ( C2Wfact_Gz(ij,kmin,1) * rhogvx(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,2) * rhogvx(ij,kmin-1) &
                             + C2Wfact_Gz(ij,kmin,3) * rhogvy(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,4) * rhogvy(ij,kmin-1) &
                             + C2Wfact_Gz(ij,kmin,5) * rhogvz(ij,kmin  ) &
                             + C2Wfact_Gz(ij,kmin,6) * rhogvz(ij,kmin-1) )
       endif

       rhogw(ij,kmin-1) = 0.0_DP

    enddo
    !$omp end parallel do
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_rhow_DP

  !-----------------------------------------------------------------------------
  subroutine BNDCND_pre_sfc_SP( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    use mod_cnst, only: &
       GRAV => CNST_GRAV
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: rho    (ijdim,kdim)
    real(SP), intent(in)  :: pre    (ijdim,kdim)
    real(SP), intent(in)  :: z      (ijdim,kdim)
    real(SP), intent(in)  :: z_srf  (ijdim)
    real(SP), intent(out) :: rho_srf(ijdim)
    real(SP), intent(out) :: pre_srf(ijdim)

    real(SP) :: z_ks, z_k1, z_k2, z_k3
    integer  :: ij

    real(SP) :: zz,z1,z2,z3,p1,p2,p3
    real(SP) :: lag_intpl
    lag_intpl(zz,z1,p1,z2,p2,z3,p3) = ((zz-z2)*(zz-z3))/((z1-z2)*(z1-z3))*p1 &
                                    + ((zz-z1)*(zz-z3))/((z2-z1)*(z2-z3))*p2 &
                                    + ((zz-z1)*(zz-z2))/((z3-z1)*(z3-z2))*p3
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))

    !$acc kernels async(1) pcopy(rho_srf,pre_srf) pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))
    do ij = 1, ijdim
       z_ks = z_srf(ij)
       z_k1 = z(ij,kmin  )
       z_k2 = z(ij,kmin+1)
       z_k3 = z(ij,kmin+2)

       ! surface density ( extrapolation )
       rho_srf(ij) = lag_intpl( z_ks, z_k1, rho(ij,kmin  ), &
                                      z_k2, rho(ij,kmin+1), &
                                      z_k3, rho(ij,kmin+2)  )

       ! surface pressure ( hydrostatic balance )
       pre_srf(ij) = pre(ij,kmin) &
                   + 0.5_SP * ( rho_srf(ij) + rho(ij,kmin) ) * GRAV * ( z_k1 - z_ks )
    enddo
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_pre_sfc_SP

  !-----------------------------------------------------------------------------
  subroutine BNDCND_pre_sfc_SPDP( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    use mod_cnst, only: &
       GRAV => CNST_GRAV
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: rho    (ijdim,kdim)
    real(DP), intent(in)  :: pre    (ijdim,kdim)
    real(SP), intent(in)  :: z      (ijdim,kdim)
    real(SP), intent(in)  :: z_srf  (ijdim)
    real(DP), intent(out) :: rho_srf(ijdim)
    real(DP), intent(out) :: pre_srf(ijdim)

    real(DP) :: z_ks, z_k1, z_k2, z_k3
    integer  :: ij

    real(DP) :: zz,z1,z2,z3,p1,p2,p3
    real(DP) :: lag_intpl
    lag_intpl(zz,z1,p1,z2,p2,z3,p3) = ((zz-z2)*(zz-z3))/((z1-z2)*(z1-z3))*p1 &
                                    + ((zz-z1)*(zz-z3))/((z2-z1)*(z2-z3))*p2 &
                                    + ((zz-z1)*(zz-z2))/((z3-z1)*(z3-z2))*p3
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))

    !$acc kernels async(1) pcopy(rho_srf,pre_srf) pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))
    do ij = 1, ijdim
       z_ks = real(z_srf(ij)   ,kind=DP)
       z_k1 = real(z(ij,kmin  ),kind=DP)
       z_k2 = real(z(ij,kmin+1),kind=DP)
       z_k3 = real(z(ij,kmin+2),kind=DP)

       ! surface density ( extrapolation )
       rho_srf(ij) = lag_intpl( z_ks, z_k1, rho(ij,kmin  ), &
                                      z_k2, rho(ij,kmin+1), &
                                      z_k3, rho(ij,kmin+2)  )

       ! surface pressure ( hydrostatic balance )
       pre_srf(ij) = pre(ij,kmin) &
                   + 0.5_DP * ( rho_srf(ij) + rho(ij,kmin) ) * GRAV * ( z_k1 - z_ks )
    enddo
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_pre_sfc_SPDP

  !-----------------------------------------------------------------------------
  subroutine BNDCND_pre_sfc_DP( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    use mod_cnst, only: &
       GRAV => CNST_GRAV
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: rho    (ijdim,kdim)
    real(DP), intent(in)  :: pre    (ijdim,kdim)
    real(DP), intent(in)  :: z      (ijdim,kdim)
    real(DP), intent(in)  :: z_srf  (ijdim)
    real(DP), intent(out) :: rho_srf(ijdim)
    real(DP), intent(out) :: pre_srf(ijdim)

    real(DP) :: z_ks, z_k1, z_k2, z_k3
    integer  :: ij

    real(DP) :: zz,z1,z2,z3,p1,p2,p3
    real(DP) :: lag_intpl
    lag_intpl(zz,z1,p1,z2,p2,z3,p3) = ((zz-z2)*(zz-z3))/((z1-z2)*(z1-z3))*p1 &
                                    + ((zz-z1)*(zz-z3))/((z2-z1)*(z2-z3))*p2 &
                                    + ((zz-z1)*(zz-z2))/((z3-z1)*(z3-z2))*p3
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))

    !$acc kernels async(1) pcopy(rho_srf,pre_srf) pcopyin(z_srf,z,rho,pre) if(acc_is_present(rho_srf))
    do ij = 1, ijdim
       z_ks = z_srf(ij)
       z_k1 = z(ij,kmin  )
       z_k2 = z(ij,kmin+1)
       z_k3 = z(ij,kmin+2)

       ! surface density ( extrapolation )
       rho_srf(ij) = lag_intpl( z_ks, z_k1, rho(ij,kmin  ), &
                                      z_k2, rho(ij,kmin+1), &
                                      z_k3, rho(ij,kmin+2)  )

       ! surface pressure ( hydrostatic balance )
       pre_srf(ij) = pre(ij,kmin) &
                   + 0.5_DP * ( rho_srf(ij) + rho(ij,kmin) ) * GRAV * ( z_k1 - z_ks )
    enddo
    !$acc end kernels


    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine BNDCND_pre_sfc_DP

end module mod_bndcnd
