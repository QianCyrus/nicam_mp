!-------------------------------------------------------------------------------
!
!+  vertical advection module
!
!-------------------------------------------------------------------------------
module mod_vadv1d
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module contains subroutines for 1d veritical advection.
  !
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17  Imported from Satoh's nonhydro model.
  !      0.03      04-05-31  Default : -> SL0_UPWIND
  !                09-07-10  H.Tomita : add the new version 'vadv1d_getflux'.
  !                11-11-18  C.Kodama: [mod] sl0_upwind and sl0_upwind_new
  !                          to avoid illegal memory reference.
  !                12-02-01  T.Seiki : fix undefined reference
  !                12-02-07  T.Seiki : typo
  !                12-05-18  H.Yashiro : vadv1d_getflux_new, faster version
  !                12-06-06  H.Yashiro : bugfix
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision, only: &
     SP, DP, RP => RP_MP
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
  public :: vadv1d_init
  public :: vadv1d_getflux      !--- Defalut VESRION
  public :: vadv1d_setup
  public :: vadv1d_prep
  public :: vadv1d_getflux_new
  public :: advection_v1d_flux  !--- OLD OLD VERSION
  public :: vadv1d_getflux_0s_DP
  public :: vadv1d_getflux_1s_DP
  public :: vadv1d_getflux_2s_DP

  interface vadv1d_getflux
     module procedure vadv1d_getflux_SP
     module procedure vadv1d_getflux_DP
  end interface vadv1d_getflux

  interface vadv1d_setup
     module procedure vadv1d_setup_SP
     module procedure vadv1d_setup_DP
  end interface vadv1d_setup

  interface vadv1d_prep
     module procedure vadv1d_prep_SP
     module procedure vadv1d_prep_DP
  end interface vadv1d_prep

  interface vadv1d_getflux_new
     module procedure vadv1d_getflux_new_SP
     module procedure vadv1d_getflux_new_DP
  end interface vadv1d_getflux_new

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: sl0_upwind_new     !--- Defalut VESRION
  private :: sl0_upwind
  private :: ppm
  private :: csl3_hyman
  private :: csl3_uno
  private :: csl3_cw
  private :: sub_fh_cell
  private :: sub_fh_3rd

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  character(len=H_SHORT), private :: ADVECTION_TYPE_V1D = 'SL0_UPWIND'
  !                                                        = 'UPWIND'
  !                                                        = 'SL0_UPWIND'
  !                                                        = 'PPM'
  !                                                        = 'CSL3_HYMAN'
  !                                                        = 'CSL3_UNO'
  !                                                        = 'CSL3_CW'
  logical :: flg_getflux_0s_DP = .false.
  integer, allocatable :: kcell_max_lim(:)
  integer, allocatable :: kcell_min_lim(:)

  integer, allocatable :: kcell(:,:)
  real(SP), allocatable :: dw_SP(:,:)
  real(DP), allocatable :: dw_DP(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vadv1d_init( )
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_kall
    implicit none
#ifdef _OPENACC
    allocate ( kcell(ADM_gall_in,ADM_kall) )
    allocate ( dw_SP(ADM_gall_in,ADM_kall) )
    allocate ( dw_DP(ADM_gall_in,ADM_kall) )
    !$acc enter data create(kcell,dw_SP,dw_DP)
#endif
    return
  end subroutine vadv1d_init
  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_SP( &
       ijdim, &
       kdim,  &
       kmin,  &
       kmax,  &
       dz,    &
       zh,    &
       rhof,  &
       wp,    &
       frhof, &
       dt     )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(SP), intent(in)  :: dz   (ijdim,kdim)
    real(SP), intent(in)  :: zh   (ijdim,kdim)
    real(SP), intent(in)  :: rhof (ijdim,kdim) ! p-grid
    real(SP), intent(in)  :: wp   (ijdim,kdim) ! p-grid
    real(SP), intent(out) :: frhof(ijdim,kdim) ! w-grid
    real(DP), intent(in)  :: dt

    integer  :: kcell    (ijdim,kdim) ! w-grid
    integer  :: kcell_max(kdim)       ! w-grid
    integer  :: kcell_min(kdim)       ! w-grid
    real(SP) :: wh       (ijdim,kdim) ! w-grid
    real(SP) :: zdis     (ijdim,kdim) ! w-grid
    real(SP) :: zzmax, zzmin
    integer  :: ij, k, kk
    integer  :: k2
    !---------------------------------------------------------------------------

    !--- vetical velocity at the half level
    do k = kmin+1, kmax
       wh(:,k) = ( wp(:,k-1) + wp(:,k) ) / 2.0_SP
    enddo
    !------ bottom boundary for wh
    wh(:,kmin) = wp(:,kmin)
    wh(:,kmin-1) = wp(:,kmin-1)
    !------ topboundary for wh : same as inner region
    wh(:,kmax+1) = wp(:,kmax)
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) / 2.0_SP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(:,k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) ) &
            / ( dz(:,k-1) + dz(:,k) ) * 2.0_SP ) ) / 6.0_SP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_SP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(:,k) / 2.0_SP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) / 2.0_SP
    !--- calculation of kcell
    !------  top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    ! [fix] 2012/02/01, T.Seiki, to avoid undefined reference
!!$ do k = kmin, kmax+1
    do k  = 1, kdim
    do ij = 1, ijdim
       kcell(ij,k) = k
    enddo
    enddo
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       zzmin = minval(zdis(:,k))
       kcell_min(k) = k
       kcell_max(k) = k
       if ( zzmax > 0.0_SP ) then
          do k2 = k, kmin, -1
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmax  &
                  .AND. zh(ij,k) - zzmax < zh(ij,k2+1) ) then
                kcell_min(k) = min( kcell_min(k), k2 )
             endif
          enddo
          enddo
       endif
       if ( zzmin < 0.0_SP ) then
          do k2 = k, kmax
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmin  &
                  .AND. zh(ij,k) - zzmin < zh(ij,k2+1) ) then
                kcell_max(k) = max( kcell_max(k), k2 )
             endif
          enddo
          enddo
       endif
    enddo
    !------ determine the kcell at each point.
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          kcell(:,k) = k
       else
          kcell(:,k) = 0
          do k2 = kcell_min(k), kcell_max(k)
             kcell(:,k) = max ( kcell(:,k), &
                  int ( k2 &
                  * sign ( 1.0_SP, ( zh(:,k) - zdis(:,k) ) - zh(:,k2) ) &
                  * sign ( 1.0_SP, zh(:,k2+1) - ( zh(:,k) - zdis(:,k) ) ) ) )
          enddo
       endif
    enddo
    !--- set zero value for norain, see if block in do loop
    frhof(:,:) = 0.0_SP
    !------ integration in the integer cells
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          ! frhof(:,:) = 0.0_SP
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             frhof(:,k) &
                  = frhof(:,k) &
                  + rhof(:,k2) * dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_SP ) / 2.0_SP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_SP ) / 2.0_SP &
                  + ( - 1.0_SP ) * ( sign ( 1, k2 - k ) + 1.0_SP ) / 2.0_SP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_SP ) / 2.0_SP )
             zdis(:,k) = zdis(:,k) - dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_SP ) / 2.0_SP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_SP ) / 2.0_SP &
                  + ( - 1.0_SP ) * ( sign ( 1, k2 - k ) + 1.0_SP ) / 2.0_SP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_SP ) / 2.0_SP )
          enddo
       endif
    enddo

    ! top boundary: rigid
    !  zdis(:,kmax+1) = 0.0_SP
    !  frhof(:,kmax+1) = 0.0_SP

    !--- integration in the non-integer cell ( depending on the scheme below )
    do k  = kmin, kmax
    do ij = 1, ijdim
       kk = max( min( kcell(ij,k), kmax ), kmin )

       frhof(ij,k) = frhof(ij,k) + rhof(ij,kk) * zdis(ij,k)
    enddo
    enddo

    return
  end subroutine vadv1d_getflux_SP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_DP( &
       ijdim, &
       kdim,  &
       kmin,  &
       kmax,  &
       dz,    &
       zh,    &
       rhof,  &
       wp,    &
       frhof, &
       dt     )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(DP), intent(in)  :: dz   (ijdim,kdim)
    real(DP), intent(in)  :: zh   (ijdim,kdim)
    real(DP), intent(in)  :: rhof (ijdim,kdim) ! p-grid
    real(DP), intent(in)  :: wp   (ijdim,kdim) ! p-grid
    real(DP), intent(out) :: frhof(ijdim,kdim) ! w-grid
    real(DP), intent(in)  :: dt

    integer  :: kcell    (ijdim,kdim) ! w-grid
    integer  :: kcell_max(kdim)       ! w-grid
    integer  :: kcell_min(kdim)       ! w-grid
    real(DP) :: wh       (ijdim,kdim) ! w-grid
    real(DP) :: zdis     (ijdim,kdim) ! w-grid
    real(DP) :: zzmax, zzmin
    integer  :: ij, k, kk
    integer  :: k2, k3
    !---------------------------------------------------------------------------

    if(.not. flg_getflux_0s_DP)then
      call vadv1d_getflux_0s_DP &
            ( ijdim &
            , kdim  &
            , kmin  &
            , kmax  &
            , dz    &
            , zh    &
            , dt    &
            )
    endif

    !--- vetical velocity at the half level
    do k = kmin+1, kmax
       wh(:,k) = ( wp(:,k-1) + wp(:,k) ) / 2.0_DP
    enddo
    !------ bottom boundary for wh
    wh(:,kmin) = wp(:,kmin)
    wh(:,kmin-1) = wp(:,kmin-1)
    !------ topboundary for wh : same as inner region
    wh(:,kmax+1) = wp(:,kmax)
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) / 2.0_DP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(:,k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) ) &
            / ( dz(:,k-1) + dz(:,k) ) * 2.0_DP ) ) / 6.0_DP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_DP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(:,k) / 2.0_DP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) / 2.0_DP
    !--- calculation of kcell
    !------  top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    ! [fix] 2012/02/01, T.Seiki, to avoid undefined reference
!!$ do k = kmin, kmax+1
    do k  = 1, kdim
    do ij = 1, ijdim
       kcell(ij,k) = k
    enddo
    enddo
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       zzmin = minval(zdis(:,k))
       kcell_min(k) = k
       kcell_max(k) = k
       if ( zzmax > 0.0_DP ) then
          k3 = max(kmin,kcell_min_lim(k))
          do k2 = k, k3, -1
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmax  &
                  .AND. zh(ij,k) - zzmax < zh(ij,k2+1) ) then
                kcell_min(k) = min( kcell_min(k), k2 )
             endif
          enddo
          enddo
       endif
       if ( zzmin < 0.0_DP ) then
          k3 = min(kmax,kcell_max_lim(k))
          do k2 = k, k3
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmin  &
                  .AND. zh(ij,k) - zzmin < zh(ij,k2+1) ) then
                kcell_max(k) = max( kcell_max(k), k2 )
             endif
          enddo
          enddo
       endif
    enddo
    !------ determine the kcell at each point.
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          kcell(:,k) = k
       else
          kcell(:,k) = 0
          do k2 = kcell_min(k), kcell_max(k)
             kcell(:,k) = max ( kcell(:,k), &
                  int ( k2 &
                  * sign ( 1.0_DP, ( zh(:,k) - zdis(:,k) ) - zh(:,k2) ) &
                  * sign ( 1.0_DP, zh(:,k2+1) - ( zh(:,k) - zdis(:,k) ) ) ) )
          enddo
       endif
    enddo
    !--- set zero value for norain, see if block in do loop
    frhof(:,:) = 0.0_DP
    !------ integration in the integer cells
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          ! frhof(:,:) = 0.0_DP
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             frhof(:,k) &
                  = frhof(:,k) &
                  + rhof(:,k2) * dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
             zdis(:,k) = zdis(:,k) - dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
          enddo
       endif
    enddo

    ! top boundary: rigid
    !  zdis(:,kmax+1) = 0.0_SP
    !  frhof(:,kmax+1) = 0.0_SP

    !--- integration in the non-integer cell ( depending on the scheme below )
    do k  = kmin, kmax
    do ij = 1, ijdim
       kk = max( min( kcell(ij,k), kmax ), kmin )

       frhof(ij,k) = frhof(ij,k) + rhof(ij,kk) * zdis(ij,k)
    enddo
    enddo

    return
  end subroutine vadv1d_getflux_DP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_setup_SP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       zh,        &
       kcell_lim, &
       dt         )
    use mod_runconf, only: &
       W_VADV1D_MIN
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(SP), intent(in)  :: dz       (ijdim,kdim)
    real(SP), intent(in)  :: zh       (ijdim,kdim)
    integer,  intent(out) :: kcell_lim(kdim)
    real(DP), intent(in)  :: dt

#ifndef _OPENACC
    integer  :: kcell(ijdim,kdim)
#endif
    real(SP) :: zdis, ztmp
    integer  :: vmax

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------
    !!$acc data pcreate(kcell) pcopyin(dz,zh) pcopy(kcell_lim)

    call PROF_rapstart('_____Precip_vadv1d_prep',3)

    !--- vertical velocity at the half level

!OCL XFILL, PREFETCH_STRONG
    !$acc kernels pcopy(kcell)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,kcell)
    do k  = kmin, kmax
    do ij = 1, ijdim
       kcell(ij,k) = k
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !--- Movement distance of the cell wall by the fall

    zdis = - W_VADV1D_MIN * dt

    !--- search number of source cell

    !OCL simd
    !$acc kernels pcopy(kcell) pcopyin(zh)
    !$acc loop independent
    !$omp parallel do default(none),private(ij,k,k2,ztmp), &
    !$omp shared(ijdim,kmin,kmax,kcell,zh,zdis), &
    !$omp schedule(static,1)
    do k2 = kmin, kmax
#ifdef _OPENACC
    !$acc loop independent
    do ij = 1, ijdim
    !$acc loop seq
    do k  = k2, kmax
#else
    do k  = k2, kmax
    do ij = 1, ijdim
#endif
       ztmp = zh(ij,k2) + zdis
       if (       ztmp >  zh(ij,k  ) &
            .AND. ztmp <= zh(ij,k+1) ) then
          kcell(ij,k2) = k
       endif
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef _OPENACC
    !$acc kernels pcopy(kcell_lim) pcopyin(kcell)
    do k = kmin, kmax
       vmax = 0
       !$acc loop reduction(max:vmax)
       do ij = 1, ijdim
          vmax = max(vmax,kcell(ij,k))
       enddo
       kcell_lim(k) = vmax
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(k), &
    !$omp shared(kmin,kmax,kcell_lim,kcell)
    do k = kmin, kmax
       kcell_lim(k) = maxval(kcell(:,k))
    enddo
    !$omp end parallel do
#endif

    call PROF_rapend  ('_____Precip_vadv1d_prep',3)

    !!$acc end data
    return
  end subroutine vadv1d_setup_SP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_setup_DP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       zh,        &
       kcell_lim, &
       dt         )
    use mod_runconf, only: &
       W_VADV1D_MIN
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(DP), intent(in)  :: dz       (ijdim,kdim)
    real(DP), intent(in)  :: zh       (ijdim,kdim)
    integer,  intent(out) :: kcell_lim(kdim)
    real(DP), intent(in)  :: dt

#ifndef _OPENACC
    integer  :: kcell(ijdim,kdim)
#endif
    real(DP) :: zdis, ztmp
    integer  :: vmax

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------
    !!$acc data pcreate(kcell) pcopyin(dz,zh) pcopy(kcell_lim)

    call PROF_rapstart('_____Precip_vadv1d_prep',3)

    !--- vertical velocity at the half level

!OCL XFILL, PREFETCH_STRONG
    !$acc kernels pcopy(kcell) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,kcell)
    do k  = kmin, kmax
    do ij = 1, ijdim
       kcell(ij,k) = k
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !--- Movement distance of the cell wall by the fall

    zdis = - W_VADV1D_MIN * dt

    !--- search number of source cell

    !OCL simd
    !$acc kernels pcopy(kcell) pcopyin(zh) async(1)
    !$acc loop independent
    !$omp parallel do default(none),private(ij,k,k2,ztmp), &
    !$omp shared(ijdim,kmin,kmax,kcell,zh,zdis), &
    !$omp schedule(static,1)
    do k2 = kmin, kmax
#ifdef _OPENACC
    !$acc loop independent
    do ij = 1, ijdim
    !$acc loop seq
    do k  = k2, kmax
#else
    do k  = k2, kmax
    do ij = 1, ijdim
#endif
       ztmp = zh(ij,k2) + zdis
       if (       ztmp >  zh(ij,k  ) &
            .AND. ztmp <= zh(ij,k+1) ) then
          kcell(ij,k2) = k
       endif
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef _OPENACC
    !$acc kernels pcopy(kcell_lim) pcopyin(kcell) async(1)
    do k = kmin, kmax
       vmax = 0
       !$acc loop reduction(max:vmax)
       do ij = 1, ijdim
          vmax = max(vmax,kcell(ij,k))
       enddo
       kcell_lim(k) = vmax
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(k), &
    !$omp shared(kmin,kmax,kcell_lim,kcell)
    do k = kmin, kmax
       kcell_lim(k) = maxval(kcell(:,k))
    enddo
    !$omp end parallel do
#endif

    call PROF_rapend  ('_____Precip_vadv1d_prep',3)

    !!$acc end data
    return
  end subroutine vadv1d_setup_DP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_prep_SP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       dz2,       &
       zh,        &
       w,         &
       kcell_lim, &
       zdis,      &
       kcell,     &
       kcell_max, &
       dt         )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(SP), intent(in)  :: dz       (ijdim,kdim)
    real(SP), intent(in)  :: dz2      (ijdim,kdim)
    real(SP), intent(in)  :: zh       (ijdim,kdim)
    real(SP), intent(in)  :: w        (ijdim,kdim)
    integer,  intent(in)  :: kcell_lim(kdim)
    real(SP), intent(out) :: zdis     (ijdim,kdim)
    integer,  intent(out) :: kcell    (ijdim,kdim)
    integer,  intent(out) :: kcell_max(kdim)
    real(DP), intent(in)  :: dt

#ifndef _OPENACC
    real(SP) :: dw_SP (ijdim,kdim)
#endif
    real(SP) :: ztmp
    integer  :: vmax

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------
    !!$acc data &
    !!$acc pcreate(dw) &
    !!$acc pcopyin(dz,dz2,zh,w,kcell_lim) &
    !!$acc pcopyin(zdis,kcell,kcell_max)

    call PROF_rapstart('_____Precip_vadv1d_prep',3)
#ifdef MORETIMER
call PROF_rapstart('______VADV1D_01')
#endif

    !--- vertical velocity at the half level

!OCL XFILL, PREFETCH_STRONG
    !$acc kernels pcopy(dw_SP,kcell) pcopyin(w) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,dw_SP,kcell,w)
    do k  = kmin, kmax
    do ij = 1, ijdim
       dw_SP   (ij,k) = w(ij,k+1) - w(ij,k)
       kcell(ij,k) = k
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_01')
call PROF_rapstart('______VADV1D_02')
#endif

    !--- Movement distance of the cell wall by the fall

    !$acc kernels pcopy(zdis) pcopyin(dz,w,dz2,dw_SP) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,zdis,w,dz,dw_SP,dz2,dt)
    do k  = kmin+1, kmax
    do ij = 1, ijdim
       zdis(ij,k) = - w(ij,k) * dt                                                                                          &
                    + w(ij,k) * dt**2 / 2.0_SP *     ( dw_SP(ij,k)+dw_SP(ij,k-1) )/dz2(ij,k)                                      &
                    - w(ij,k) * dt**3 / 6.0_SP * ( ( ( dw_SP(ij,k)+dw_SP(ij,k-1) )/dz2(ij,k) )**2                                 &
                                               + 2.0_SP * w(ij,k)/dz2(ij,k) * ( dw_SP(ij,k)/dz(ij,k) - dw_SP(ij,k-1)/dz(ij,k-1) ) )
       zdis(ij,k) = max( zdis(ij,k), 0.0_SP )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc kernels pcopy(zdis) pcopyin(dz,w,dw_SP) async(1)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,zdis,w,dt,dw_SP,dz)
    do ij = 1, ijdim
       zdis(ij,kmin) = - w(ij,kmin) * dt &
                       + w(ij,kmin) * dt**2 / 2.0_SP * dw_SP(ij,kmin) / dz(ij,kmin)
       zdis(ij,kmin) = max( zdis(ij,kmin), 0.0_SP )
    enddo
    !$omp end parallel do
    !$acc end kernels
#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_02')
call PROF_rapstart('______VADV1D_03')
#endif

    !--- search number of source cell

    !OCL simd
    !OCL loop_nointerchange
    !$acc kernels pcopy(kcell) pcopyin(kcell_lim,zdis,zh) async(1)
    !$omp parallel do default(none),private(ij,k,k2,ztmp), &
    !$omp shared(ijdim,kmin,kmax,kcell_lim,kcell,zh,zdis), &
    !$omp schedule(static,1)
    !$acc loop independent
    do k2 = kmin, kmax
#ifdef _OPENACC
    !$acc loop independent
    do ij = 1, ijdim
    !$acc loop seq
    do k  = k2, kcell_lim(k2)
#else
    do k  = k2, kcell_lim(k2)
    do ij = 1, ijdim
#endif
       ztmp = zh(ij,k2) + zdis(ij,k2)
       if (       ztmp >  zh(ij,k  ) &
            .AND. ztmp <= zh(ij,k+1) ) then
          kcell(ij,k2) = k
       endif
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_03')
call PROF_rapstart('______VADV1D_04')
#endif

#ifdef _OPENACC
    !$acc kernels pcopy(kcell_max) pcopyin(kcell) async(1)
    do k = kmin, kmax
       vmax = 0
       !$acc loop reduction(max:vmax)
       do ij = 1, ijdim
          vmax = max(vmax,kcell(ij,k))
       enddo
       kcell_max(k) = vmax
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(k), &
    !$omp shared(kmin,kmax,kcell_max,kcell)
    do k = kmin, kmax
       kcell_max(k) = maxval(kcell(:,k))
    enddo
    !$omp end parallel do
#endif

#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_04')
#endif

    call PROF_rapend  ('_____Precip_vadv1d_prep',3)

    !!$acc end data
    return
  end subroutine vadv1d_prep_SP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_prep_DP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       dz2,       &
       zh,        &
       w,         &
       kcell_lim, &
       zdis,      &
       kcell,     &
       kcell_max, &
       dt         )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(DP), intent(in)  :: dz       (ijdim,kdim)
    real(DP), intent(in)  :: dz2      (ijdim,kdim)
    real(DP), intent(in)  :: zh       (ijdim,kdim)
    real(DP), intent(in)  :: w        (ijdim,kdim)
    integer,  intent(in)  :: kcell_lim(kdim)
    real(DP), intent(out) :: zdis     (ijdim,kdim)
    integer,  intent(out) :: kcell    (ijdim,kdim)
    integer,  intent(out) :: kcell_max(kdim)
    real(DP), intent(in)  :: dt

#ifndef _OPENACC
    real(DP) :: dw_DP (ijdim,kdim)
#endif
    real(DP) :: ztmp
    integer  :: vmax

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------
    !!$acc data &
    !!$acc pcreate(dw_DP) &
    !!$acc pcopyin(dz,dz2,zh,w,kcell_lim) &
    !!$acc pcopyin(zdis,kcell,kcell_max)

    call PROF_rapstart('_____Precip_vadv1d_prep',3)
#ifdef MORETIMER
call PROF_rapstart('______VADV1D_01')
#endif

    !--- vertical velocity at the half level

!OCL XFILL, PREFETCH_STRONG
    !$acc kernels pcopy(dw_DP,kcell) pcopyin(w) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,dw_DP,kcell,w)
    do k  = kmin, kmax
    do ij = 1, ijdim
       dw_DP   (ij,k) = w(ij,k+1) - w(ij,k)
       kcell(ij,k) = k
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_01')
call PROF_rapstart('______VADV1D_02')
#endif

    !--- Movement distance of the cell wall by the fall

    !$acc kernels pcopy(zdis) pcopyin(dz,w,dz2,dw_DP) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,zdis,w,dz,dw_DP,dz2,dt)
    do k  = kmin+1, kmax
    do ij = 1, ijdim
       zdis(ij,k) = - w(ij,k) * dt                                                                                          &
                    + w(ij,k) * dt**2 / 2.0_DP *     ( dw_DP(ij,k)+dw_DP(ij,k-1) )/dz2(ij,k)                                      &
                    - w(ij,k) * dt**3 / 6.0_DP * ( ( ( dw_DP(ij,k)+dw_DP(ij,k-1) )/dz2(ij,k) )**2                                 &
                                               + 2.0_DP * w(ij,k)/dz2(ij,k) * ( dw_DP(ij,k)/dz(ij,k) - dw_DP(ij,k-1)/dz(ij,k-1) ) )
       zdis(ij,k) = max( zdis(ij,k), 0.0_DP )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    !$acc kernels pcopy(zdis) pcopyin(dz,w,dw_DP) async(1)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,zdis,w,dt,dw_DP,dz)
    do ij = 1, ijdim
       zdis(ij,kmin) = - w(ij,kmin) * dt &
                       + w(ij,kmin) * dt**2 / 2.0_DP * dw_DP(ij,kmin) / dz(ij,kmin)
       zdis(ij,kmin) = max( zdis(ij,kmin), 0.0_DP )
    enddo
    !$omp end parallel do
    !$acc end kernels
#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_02')
call PROF_rapstart('______VADV1D_03')
#endif

    !--- search number of source cell

    !OCL simd
    !OCL loop_nointerchange
    !$acc kernels pcopy(kcell) pcopyin(kcell_lim,zdis,zh) async(1)
    !$omp parallel do default(none),private(ij,k,k2,ztmp), &
    !$omp shared(ijdim,kmin,kmax,kcell_lim,kcell,zh,zdis), &
    !$omp schedule(static,1)
    !$acc loop independent
    do k2 = kmin, kmax
#ifdef _OPENACC
    !$acc loop independent
    do ij = 1, ijdim
    !$acc loop seq
    do k  = k2, kcell_lim(k2)
#else
    do k  = k2, kcell_lim(k2)
    do ij = 1, ijdim
#endif
       ztmp = zh(ij,k2) + zdis(ij,k2)
       if (       ztmp >  zh(ij,k  ) &
            .AND. ztmp <= zh(ij,k+1) ) then
          kcell(ij,k2) = k
       endif
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_03')
call PROF_rapstart('______VADV1D_04')
#endif

#ifdef _OPENACC
    !$acc kernels pcopy(kcell_max) pcopyin(kcell) async(1)
    do k = kmin, kmax
       vmax = 0
       !$acc loop reduction(max:vmax)
       do ij = 1, ijdim
          vmax = max(vmax,kcell(ij,k))
       enddo
       kcell_max(k) = vmax
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(k), &
    !$omp shared(kmin,kmax,kcell_max,kcell)
    do k = kmin, kmax
       kcell_max(k) = maxval(kcell(:,k))
    enddo
    !$omp end parallel do
#endif

#ifdef MORETIMER
call PROF_rapend  ('______VADV1D_04')
#endif

    call PROF_rapend  ('_____Precip_vadv1d_prep',3)

    !!$acc end data
    return
  end subroutine vadv1d_prep_DP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_new_SP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       rhoq,      &
       zdis0,     &
       kcell,     &
       kcell_max, &
       frhoq      )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(SP), intent(in)  :: dz       (ijdim,kdim)
    real(SP), intent(in)  :: rhoq     (ijdim,kdim)
    real(SP), intent(in)  :: zdis0    (ijdim,kdim) ! [bugfix] H.Yashiro 20120606
    integer,  intent(in)  :: kcell    (ijdim,kdim)
    integer,  intent(in)  :: kcell_max(kdim)
    real(SP), intent(out) :: frhoq    (ijdim,kdim)

#ifdef _OPENACC
    real(SP) :: zdis
#else
    real(SP) :: zdis(ijdim,kdim)
#endif
    real(SP) :: mask

    integer  :: ij, k, k2, kmax2
    !---------------------------------------------------------------------------
    !!$acc data &
    !!$acc pcopyin(dz,rhoq,zdis0,kcell,kcell_max)

    call PROF_rapstart('_____Precip_vadv1d_getflux',3)

    !------ integration in the integer cells

#ifdef _OPENACC
    !$acc kernels pcopy(frhoq) pcopy(dz,rhoq,kcell,kcell_max) async(1)
    do ij = 1, ijdim
       !$acc loop independent
       do k2 = 1, kdim
          frhoq(ij,k2) = 0.0_DP
          if (k2 < kmin .OR. k2 > kmax) cycle
          kmax2 = kcell_max(k2)-1

          zdis = zdis0(ij,k2)
          !$acc loop seq
          do k  = k2, kmax2
             mask = 0.5_SP * ( sign(1,kcell(ij,k2)-1-k) + 1.0_SP ) ! if k < kcell(ij,k2), mask = 1.0

             frhoq(ij,k2) = frhoq(ij,k2) - mask * dz(ij,k) * rhoq(ij,k) ! sum column mass rhoq*dz
             zdis         = zdis         - mask * dz(ij,k)              ! residual
          enddo

          ! residual (simple upwind)
          mask = 0.5_SP * ( sign(1,kcell(ij,k2)-k2) + 1.0_SP ) ! if k2 <= kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * zdis * rhoq(ij,kcell(ij,k2)) ! sum column mass rhoq*dz
       enddo
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(ij,k2), &
    !$omp shared(ijdim,kdim,frhoq)
!OCL NOPREFETCH
    do k2 = 1, kdim
    do ij = 1, ijdim
       frhoq(ij,k2) = 0.0_SP
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
    !OCL XFILL
!OCL NOPREFETCH
    do k2 = kmin, kmax
       do ij = 1, ijdim
          zdis(ij,k2) = zdis0(ij,k2)
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
!OCL NOPREFETCH
    do k2 = kmin, kmax
       kmax2 = kcell_max(k2)-1
       !OCL loop_nointerchange
       do k  = k2, kmax2
       do ij = 1, ijdim
          mask = 0.5_SP * ( sign(1,kcell(ij,k2)-1-k) + 1.0_SP ) ! if k < kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * dz(ij,k) * rhoq(ij,k) ! sum column mass rhoq*dz
          zdis (ij,k2) = zdis (ij,k2) - mask * dz(ij,k)              ! residual
       enddo
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
!OCL NOPREFETCH
    do k2 = kmin, kmax

       ! residual (simple upwind)
       do ij = 1, ijdim
          mask = 0.5_SP * ( sign(1,kcell(ij,k2)-k2) + 1.0_SP ) ! if k2 <= kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * zdis(ij,k2) * rhoq(ij,kcell(ij,k2)) ! sum column mass rhoq*dz
       enddo
    enddo
    !$omp end parallel do
#endif

    call PROF_rapend  ('_____Precip_vadv1d_getflux',3)

    !!$acc end data
    return
  end subroutine vadv1d_getflux_new_SP

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_new_DP( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       rhoq,      &
       zdis0,     &
       kcell,     &
       kcell_max, &
       frhoq      )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(DP), intent(in)  :: dz       (ijdim,kdim)
    real(DP), intent(in)  :: rhoq     (ijdim,kdim)
    real(DP), intent(in)  :: zdis0    (ijdim,kdim) ! [bugfix] H.Yashiro 20120606
    integer,  intent(in)  :: kcell    (ijdim,kdim)
    integer,  intent(in)  :: kcell_max(kdim)
    real(DP), intent(out) :: frhoq    (ijdim,kdim)

#ifdef _OPENACC
    real(DP) :: zdis
#else
    real(DP) :: zdis(ijdim,kdim)
#endif
    real(DP) :: mask

    integer  :: ij, k, k2, kmax2
    !---------------------------------------------------------------------------
    !!$acc data &
    !!$acc pcopyin(dz,rhoq,zdis0,kcell,kcell_max)

    call PROF_rapstart('_____Precip_vadv1d_getflux',3)

    !------ integration in the integer cells
#ifdef _OPENACC
    !$acc kernels pcopy(frhoq) pcopy(dz,rhoq,kcell,kcell_max) async(1)
    do ij = 1, ijdim
       !$acc loop independent
       do k2 = 1, kdim
          frhoq(ij,k2) = 0.0_DP
          if (k2 < kmin .OR. k2 > kmax) cycle
          kmax2 = kcell_max(k2)-1

          zdis = zdis0(ij,k2)
          !$acc loop seq
          do k  = k2, kmax2
             mask = 0.5_DP * ( sign(1,kcell(ij,k2)-1-k) + 1.0_DP ) ! if k < kcell(ij,k2), mask = 1.0

             frhoq(ij,k2) = frhoq(ij,k2) - mask * dz(ij,k) * rhoq(ij,k) ! sum column mass rhoq*dz
             zdis         = zdis         - mask * dz(ij,k)              ! residual
          enddo

          ! residual (simple upwind)
          mask = 0.5_DP * ( sign(1,kcell(ij,k2)-k2) + 1.0_DP ) ! if k2 <= kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * zdis * rhoq(ij,kcell(ij,k2)) ! sum column mass rhoq*dz
       enddo
    enddo
    !$acc end kernels
#else
!OCL XFILL
    !$omp parallel do default(none),private(ij,k2), &
    !$omp shared(ijdim,kdim,frhoq)
!OCL NOPREFETCH
    do k2 = 1, kdim
    do ij = 1, ijdim
       frhoq(ij,k2) = 0.0_DP
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
    !OCL XFILL
!OCL NOPREFETCH
    do k2 = kmin, kmax
       do ij = 1, ijdim
          zdis(ij,k2) = zdis0(ij,k2)
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
!OCL NOPREFETCH
    do k2 = kmin, kmax
       kmax2 = kcell_max(k2)-1
       !OCL loop_nointerchange
       do k  = k2, kmax2
       do ij = 1, ijdim
          mask = 0.5_DP * ( sign(1,kcell(ij,k2)-1-k) + 1.0_DP ) ! if k < kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * dz(ij,k) * rhoq(ij,k) ! sum column mass rhoq*dz
          zdis (ij,k2) = zdis (ij,k2) - mask * dz(ij,k)              ! residual
       enddo
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,k2,kmax2,mask), &
    !$omp shared(ijdim,kmin,kmax,zdis,zdis0,kcell_max,kcell,frhoq,rhoq,dz)
!OCL NOPREFETCH
    do k2 = kmin, kmax

       ! residual (simple upwind)
       do ij = 1, ijdim
          mask = 0.5_DP * ( sign(1,kcell(ij,k2)-k2) + 1.0_DP ) ! if k2 <= kcell(ij,k2), mask = 1.0

          frhoq(ij,k2) = frhoq(ij,k2) - mask * zdis(ij,k2) * rhoq(ij,kcell(ij,k2)) ! sum column mass rhoq*dz
       enddo
    enddo
    !$omp end parallel do
#endif

    call PROF_rapend  ('_____Precip_vadv1d_getflux',3)

    !!$acc end data
    return
  end subroutine vadv1d_getflux_new_DP

  !-----------------------------------------------------------------------------
  subroutine advection_v1d_flux( &
        ijdim,                   & !--- IN : number of horizontal grids
        rhof_DP,                    & !--- IN : rho X f X gsgam2
        wp,                      & !--- IN : verical velocity at the full level
        frhof,                   & !--- OUT : vertical flux ( X dt )
        dt,                      & !--- IN : time interval
        ADV_TYPE                 & !--- IN : (optional)
        )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       dz => GRD_dgz, &
       zh => GRD_gzh
    implicit none

    integer, intent(in) :: ijdim
    real(DP), intent(in) :: rhof_DP(1:ijdim,1:kdim)   ! p-grid
    real(DP), intent(in) :: wp(1:ijdim,1:kdim)     ! p-grid
    real(DP), intent(out) :: frhof(1:ijdim,1:kdim) ! w-grid
    real(DP), intent(in) :: dt
    character(len=*), intent(in), optional :: ADV_TYPE
    real(DP) :: wh(1:ijdim,1:kdim)          ! w-grid
    real(DP) :: rhofh_cell(1:ijdim,1:kdim)  ! w-grid
    integer :: k, k2
    integer :: kcell(1:ijdim,1:kdim) ! w-grid
    real(DP) :: zdis(1:ijdim,1:kdim)  ! w-grid
    integer :: kcell_max(1:kdim)    ! w-grid
    integer :: kcell_min(1:kdim)    ! w-grid
    real(DP) :: zzmax, zzmin
    real(DP) :: rhof(1:ijdim,1:kdim)   ! p-grid
    real(DP) :: rhofh(1:ijdim,1:kdim)
    character(len=H_SHORT) :: ADVECTION_TYPE_V1D_IN

    rhof(:,kmin:kmax) = rhof_DP(:,kmin:kmax)

    if (present(ADV_TYPE)) then
       ADVECTION_TYPE_V1D_IN = trim(ADV_TYPE)
    else
       !--- use default
       ADVECTION_TYPE_V1D_IN = trim(ADVECTION_TYPE_V1D)
    endif

    !--- vetical velocity at the half level
    do k = kmin+1, kmax
       wh(:,k) = ( wp(:,k-1) + wp(:,k) ) / 2.0_DP
    enddo
    !------ bottom boundary for wh
    wh(:,kmin) = wp(:,kmin)
    wh(:,kmin-1) = wp(:,kmin-1)
    !------ topboundary for wh : same as inner region
    wh(:,kmax+1) = wp(:,kmax)
    !--- the simplest way ( upwind )
    if ( ADVECTION_TYPE_V1D_IN == 'UPWIND' ) then
       frhof(:,kmax+1) = 0.0_DP
       frhof(:,kmin-1) = 0.0_DP
       frhof(:,kmin:kmax) = rhof(:,kmin:kmax) * wh(:,kmin:kmax) * dt
       return
    endif
    !--- for 3rd order CSL
    if ( ADVECTION_TYPE_V1D_IN(1:4) == 'CSL3' ) then
       call sub_fh_3rd ( ijdim, rhof, rhofh )
    endif
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(k-1) + dz(k) ) / 2.0_DP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(k-1) + dz(k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(k-1) ) &
            / ( dz(k-1) + dz(k) ) * 2.0_DP ) ) / 6.0_DP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_DP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(k) / 2.0_DP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(k-1) / 2.0_DP
    !--- calculation of kcell
    !------  top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    ! [fix] 2012/02/01, T.Seiki, to avoid undefined reference
!!$ do k = kmin, kmax+1
    do k = 1, kdim
       kcell(:,k) = k
    enddo
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       zzmin = minval(zdis(:,k))
       kcell_min(k) = k
       kcell_max(k) = k
       if ( zzmax > 0.0_DP ) then
          do k2 = k, kmin, -1
             if ( zh(k2) <= zh(k) - zzmax &
                  .AND. zh(k) - zzmax < zh(k2+1) ) then
                kcell_min(k) = k2
                exit
             endif
          enddo
       endif
       if ( zzmin < 0.0_DP ) then
          do k2 = k, kmax
             if ( zh(k2) <= zh(k) - zzmin &
                  .AND. zh(k) - zzmin < zh(k2+1) ) then
                kcell_max(k) = k2
                exit
             endif
          enddo
       endif
    enddo
    !------ determine the kcell at each point.
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          kcell(:,k) = k
       else
          kcell(:,k) = 0
          do k2 = kcell_min(k), kcell_max(k)
             kcell(:,k) = max ( kcell(:,k), &
                  int ( k2 &
                  * sign ( 1.0_DP, ( zh(k) - zdis(:,k) ) - zh(k2) ) &
                  * sign ( 1.0_DP, zh(k2+1) - ( zh(k) - zdis(:,k) ) ) ) )
          enddo
       endif
    enddo
    !--- set zero value for norain, see if block in do loop
    frhof(:,:) = 0.0_DP
    !------ integration in the integer cells
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          ! frhof(:,:) = 0.0_DP
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             frhof(:,k) &
                  = frhof(:,k) &
                  + rhof(:,k2) * dz(k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
             zdis(:,k) = zdis(:,k) - dz(k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
          enddo
       endif
    enddo
    ! top boundary: rigid
    !  zdis(:,kmax+1) = 0.0_DP
    !  frhof(:,kmax+1) = 0.0_DP
    !--- integration in the non-integer cell ( depending on the scheme below )
    if ( ADVECTION_TYPE_V1D_IN == 'SL0_UPWIND' ) then
       call sl0_upwind ( ijdim, rhof, kcell, rhofh_cell)
    elseif ( ADVECTION_TYPE_V1D_IN == 'PPM' ) then
       call ppm ( ijdim, rhof, rhofh_cell, kcell, zdis)
    elseif ( ADVECTION_TYPE_V1D_IN == 'CSL3_HYMAN' ) then
       call csl3_hyman ( ijdim, rhof, rhofh, rhofh_cell, kcell, zdis)
    elseif ( ADVECTION_TYPE_V1D_IN == 'CSL3_UNO' ) then
       call csl3_uno (ijdim, rhof, rhofh, rhofh_cell, kcell, zdis)
    elseif ( ADVECTION_TYPE_V1D_IN == 'CSL3_CW' ) then
       call csl3_cw (ijdim, rhof, rhofh, rhofh_cell, kcell, zdis)
    endif
    frhof(:,:) = frhof(:,:) + rhofh_cell(:,:) * zdis(:,:)
  end subroutine advection_v1d_flux
  !-----------------------------------------------------------------------------
  subroutine sl0_upwind(&
       ijdim,           & !--- IN : number of horizontal grids
       f,               & !--- IN : f
       kcell,           & !--- IN : target cell number
       fh               & !--- OUT : flux at the cell wall
       )
    !------
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz,              &
         dzh => GRD_dgzh,            &
         z => GRD_gz,                &
         zh => GRD_gzh
    !------ constant in a box
    implicit none

    integer, intent(in)    :: ijdim
    real(DP), intent(in)    :: f(1:ijdim,1:kdim)
    integer, intent(inout) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(out)   :: fh(1:ijdim,1:kdim)
    integer :: ij, k
    fh(:,:) = 0.0_DP
    ! <-- [add] C.Kodama 2011.11.18
    !           to avoid illegal memory reference
    where ( kcell(:,:) == 0 )
       kcell(:,:) = kmin
    elsewhere ( kcell(:,:) == kmax+1 )
       kcell(:,:) = kmax
    end where
    ! -->
    !<--- including the condition of [fh(:,kmax+1) = 0.0_DP]
    !<--- at the top boundary
    do k = kmin, kmax
       do ij = 1, ijdim
          fh(ij,k) = f(ij,kcell(ij,k))
       enddo
    enddo
    ! [del] C.Kodama 2011.11.18
!    do k = kmin, kmax
!       where ( kcell(:,k) == 0 )
!          fh(:,k) =  f(:,kmin)
!       end where
!       where ( kcell(:,k) == kmax+1 )
!          fh(:,k) =  f(:,kmax)
!       end where
!    enddo

    return
  end subroutine sl0_upwind


  subroutine sl0_upwind_new(&
       ijdim,           & !--- IN : number of horizontal grids
       kdim,            & !--- IN
       kmin,            & !--- IN
       kmax,            & !--- IN
       f,               & !--- IN : f
       kcell,           & !--- IN : target cell number
       fh               & !--- OUT : flux at the cell wall
       )
    !------
    !------ constant in a box
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: kmin
    integer, intent(in)    :: kmax
    real(DP), intent(in)    :: f(1:ijdim,1:kdim)
    integer, intent(inout) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(out)   :: fh(1:ijdim,1:kdim)
    integer :: ij, k
    fh(:,:) = 0.0_DP
    ! <-- [add] C.Kodama 2011.11.18
    !           to avoid illegal memory reference
    where ( kcell(:,:) == 0 )
       kcell(:,:) = kmin
    elsewhere ( kcell(:,:) == kmax+1 )
       kcell(:,:) = kmax
    end where
    ! -->
    !<--- including the condition of [fh(:,kmax+1) = 0.0_DP]
    !<--- at the top boundary
    do k = kmin, kmax
       do ij = 1, ijdim
          fh(ij,k) = f(ij,kcell(ij,k))
       enddo
    enddo
    ! [del] C.Kodama 2011.11.18
!    do k = kmin, kmax
!       where ( kcell(:,k) == 0 )
!          fh(:,k) =  f(:,kmin)
!       end where
!       where ( kcell(:,k) == kmax+1 )
!          fh(:,k) =  f(:,kmax)
!       end where
!    enddo

    return
  end subroutine sl0_upwind_new
  !-----------------------------------------------------------------------------
  subroutine ppm (&
       ijdim,     & !--- IN : number of horizontal grids
       f,         & !--- IN : f
       fh_cell,   & !--- OUT : flux at the cell wall
       kcell,     & !--- IN : target cell number
       zdis       & !--- IN : distance from kcell
       )
    !------
    !------ P. Colella and P.R. Woodward, J. Comput. Phys. 54, 174 (1984)
    !------ based on pmt.f by Xiao
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz,              &
         dzh => GRD_dgzh,            &
         z => GRD_gz,                &
         zh => GRD_gzh
    implicit none

    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh_cell(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis(1:ijdim,1:kdim)
    real(DP) :: fh(1:ijdim,1:kdim)
    real(DP) :: fr(1:ijdim,1:kdim)
    real(DP) :: fl(1:ijdim,1:kdim)
    real(DP) :: f6(1:ijdim,1:kdim)
    real(DP) :: dr, dl, ds, dn
    real(DP) :: sgn
    real(DP) :: zzz
    integer :: ij, k
    integer :: k2
    call sub_fh_3rd ( ijdim, f, fh )
    do k = kmin, kmax
       fr(:,k) = fh(:,k+1)
       fl(:,k) = fh(:,k)
    enddo
    f6(:,:) = 0.0_DP
    do k = kmin, kmax
       do ij = 1, ijdim
          dl = fl(ij,k) - f(ij,k)
          dr = fr(ij,k) - f(ij,k)
          ds = ( dr - dl ) ** 2 / 6.0_DP
          dn = ( dl ** 2 - dr ** 2 ) / 2.0_DP
          !---  if dl * dr > 0.0; then replace
          sgn = ( sign ( 1.0_DP, - ( dl * dr ) ) + 1.0_DP ) / 2.0_DP
          fr(ij,k) = fr(ij,k) * sgn + f(ij,k) * ( 1.0_DP - sgn )
          fl(ij,k) = fl(ij,k) * sgn + f(ij,k) * ( 1.0_DP - sgn )
          !---  if dn - ds > 0; then replace
          sgn = ( sign ( 1.0_DP, - ( dn - ds ) ) + 1.0_DP ) / 2.0_DP
          fl(ij,k) = fl(ij,k) * sgn &
               + ( 4.0_DP * f(ij,k) - 3.0_DP * fr(ij,k) ) * ( 1.0_DP - sgn )
          !             + ( 3.0_DP * f(ij,k) - 2.0_DP * fr(ij,k) ) * ( 1.0_DP - sgn )
          !--- if dn + ds < 0; then replace
          sgn = ( sign ( 1.0_DP, dn + ds ) + 1.0_DP ) / 2.0_DP
          fr(ij,k) = fr(ij,k) * sgn &
               + ( 4.0_DP * f(ij,k) - 3.0_DP * fl(ij,k) ) * ( 1.0_DP - sgn )
          !             + ( 3.0_DP * f(ij,k) - 2.0_DP * fl(ij,k) ) * ( 1.0_DP - sgn )
          f6(ij,k) = 6.0_DP * ( f(ij,k) - 0.5d0 * ( fl(ij,k) + fr(ij,k) ) )
       enddo
    enddo
    fh_cell(:,:) = 0.0_DP
    !<--- including the condition of [fh(:,kmax+1) = 0.0_DP]
    !<--- at the top boundary
    do k = kmin, kmax
       do ij = 1, ijdim
          k2 = kcell(ij,k)
          zzz = zdis(ij,k) / dz(k2)
          sgn = ( sign ( 1.0_DP, zzz ) + 1.0_DP ) / 2.0_DP
          fh_cell(ij,k) &
               = ( fr(ij,k2) &
               - 0.5d0 * ( ( fr(ij,k2) - fl(ij,k2) ) &
               - ( 1.0_DP - 2.0_DP / 3.0_DP * zzz ) * f6(ij,k2) ) &
               * zzz ) * sgn &
               + ( fl(ij,k2) &
               + 0.5d0 * ( ( fr(ij,k2) - fl(ij,k2) ) &
               + ( 1.0_DP - 2.0_DP / 3.0_DP * ( - zzz ) ) * f6(ij,k2) ) &
               * ( - zzz ) ) * ( 1.0_DP - sgn )
       enddo
    enddo
    do k = kmin, kmax
       where ( kcell(:,k) == 0 )
          fh_cell(:,k) =  f(:,kmin)
       end where
       where ( kcell(:,k) == kmax+1 )
          fh_cell(:,k) =  f(:,kmax)
       end where
    enddo

    return
  end subroutine ppm
  !-----------------------------------------------------------------------------
  subroutine csl3_hyman (&
       ijdim,            & !--- IN : number of horizontal grids
       f,                & !--- IN : f at the cell
       fh,               & !--- INOUT : f at the cell wall
       fh_cell,          & !--- OUT : flux at the cell wall
       kcell,            & !--- IN : target cell number
       zdis              & !--- IN : distance from kcell
       )
    !------
    !------ using Hyman's ineterpolation
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz,              &
         dzh => GRD_dgzh,            &
         z => GRD_gz,                &
         zh => GRD_gzh
    implicit none

    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(inout) :: fh(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh_cell(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis(1:ijdim,1:kdim)
    real(DP) :: h(1:ijdim,1:kdim)
    real(DP) :: deltal(1:ijdim,1:kdim)
    integer :: k
    h(:,:) = 0.0_DP
    deltal(:,:) = 0.0_DP
    do k = kmin, kmax
       h(:,k) = 3.0_DP / 2.0_DP * f(:,k) &
            - 0.25d0 * ( fh(:,k) + fh(:,k+1) )
    enddo
    !--- HYMAN APPROXIMATION ---
    do k = kmin+2, kmax-2
       deltal(:,k) &
            = ( - h(:,k+2) + 8.0_DP * h(:,k+1) &
            - 8.0_DP * h(:,k-1) + h(:,k-2) ) &
            / ( - dzh(k+2) + 7.0_DP * dzh(k+1) + 7.0_DP * dzh(k) - dzh(k-1) )
    enddo
    !! boundaries: Is this consistent with Hyman?
    do k = kmin+1, kmax-1, kmax-kmin-2
       deltal(:,k) &
            = ( dzh(k+1) ** 2 * h(:,k+1) &
            - ( dzh(k+1) ** 2 - dzh(k) ** 2 ) * h(:,k) &
            - dzh(k) ** 2 * h(:,k-1) ) &
            / dzh(k) / dzh(k+1) / ( dzh(k) + dzh(k+1) )
    enddo
    k = kmin
    deltal(:,k) = ( h(:,k+1) - h(:,k) ) / dzh(k+1)
    k = kmax
    deltal(:,k) = ( h(:,k) - h(:,k-1) ) / dzh(k)
    !---
    call sub_fh_cell (ijdim, f, fh, fh_cell, kcell, zdis, deltal)

    return
  end subroutine csl3_hyman
  !-----------------------------------------------------------------------------
  subroutine csl3_uno(&
       ijdim,         & !--- IN : number of horizontal grids
       f,             & !--- IN : f at the cell
       fh,            & !--- INOUT : f at the cell wall
       fh_cell,       & !--- OUT : flux at the cell wall
       kcell,         & !--- IN : target cell number
       zdis           & !--- IN : distance from kcell
       )
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz,              &
         dzh => GRD_dgzh,            &
         z => GRD_gz,                &
         zh => GRD_gzh
    implicit none

    !------
    !------ CSL3 with UNO
    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(inout) :: fh(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh_cell(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis(1:ijdim,1:kdim)
    real(DP) :: h(1:ijdim,1:kdim)
    real(DP) :: deltal(1:ijdim,1:kdim)
    real(DP) :: cd(1:ijdim,1:kdim)
    real(DP) :: xm(1:ijdim,1:kdim)
    integer :: ij, k
    real(DP) :: sp, sm
    real(DP) :: sgn
    h(:,:) = 0.0_DP
    deltal(:,:) = 0.0_DP
    do k = kmin, kmax
       h(:,k) = 3.0_DP / 2.0_DP * f(:,k) &
            - 0.25d0 * ( fh(:,k) + fh(:,k+1) )
    enddo
    !--- UNO reconstruction ---
    do k = kmin+1, kmax-1
       cd(:,k) = ( ( h(:,k+1) - h(:,k) ) / dzh(k+1) &
            - ( h(:,k) - h(:,k-1) ) / dzh(k) ) &
            / ( dzh(k+1) + dzh(k) )
    enddo
    do k = kmin+2, kmax-1
       do ij = 1, ijdim
          ! xm(a,b) : if |a| < |b|, xm = a; otherwise, xm = b
          ! xm(k) = xm(cd(k-1), cd(k))
          sgn = ( sign ( 1.0_DP, abs(cd(ij,k-1)) - abs(cd(ij,k)) ) &
               + 1.0_DP ) / 2.0_DP
          xm(ij,k) = cd(ij,k-1) * ( 1.0_DP - sgn ) + cd(ij,k) * sgn
          !     write(UNIT_LOGFILE,*) 'xm', cd(k-1), cd(k), xm(k), sgn
       enddo
    enddo
    do k = kmin+2, kmax-2
       do ij = 1, ijdim
          sp = ( h(ij,k+1) - h(ij,k) ) / dzh(k+1) - dzh(k+1) * xm(ij,k+1)
          sm = ( h(ij,k) - h(ij,k-1) ) / dzh(k) + dzh(k) * xm(ij,k)
          sgn = ( sign ( 1.0_DP, abs(sp) - abs(sm) ) + 1.0_DP ) / 2.0_DP
          ! deltal(k) = minmod(sp,sm)
          deltal(ij,k) = ( sp * ( 1.0_DP - sgn ) + sm * sgn ) &
               * ( sign ( 1.0_DP, sp * sm ) + 1.0_DP ) / 2.0_DP
       enddo
    enddo
    !! boundaries: Is this consistent with UNO?
    do k = kmin+1, kmax-1, kmax-kmin-2
       do ij = 1, ijdim
          sp = ( h(ij,k+1) - h(ij,k) ) / dzh(k+1)
          sm = ( h(ij,k) - h(ij,k-1) ) / dzh(k)
          sgn = ( sign ( 1.0_DP, abs(sp) - abs(sm) ) + 1.0_DP ) / 2.0_DP
          ! deltal(k) = minmod(sp,sm)
          deltal(ij,k) = ( sp * ( 1.0_DP - sgn ) + sm * sgn ) &
               * ( sign ( 1.0_DP, sp * sm ) + 1.0_DP ) / 2.0_DP
       enddo
    enddo
    k = kmin
    deltal(:,k) = ( fh(:,k+1) - fh(:,k) ) / dz(k)
    k = kmax
    deltal(:,k) = ( fh(:,k+1) - fh(:,k) ) / dz(k)
    !---
    call sub_fh_cell (ijdim, f, fh, fh_cell, kcell, zdis, deltal)

    return
  end subroutine csl3_uno
  !-----------------------------------------------------------------------------
  subroutine csl3_cw (&
       ijdim,         & !--- IN : number of horizontal grids
       f,             & !--- IN : f at the cell
       fh,            & !--- INOUT : f at the cell wall
       fh_cell,       & !--- OUT : flux at the cell wall
       kcell,         & !--- IN : target cell number
       zdis           & !--- IN : distance from kcell
       )
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz,              &
         dzh => GRD_dgzh
    implicit none

    !------
    !------ CSL3 with the same limiter as PPM.
    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(inout) :: fh(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh_cell(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis(1:ijdim,1:kdim)
    real(DP) :: h(1:ijdim,1:kdim)
    real(DP) :: deltal(1:ijdim,1:kdim)
    real(DP) :: dfi, dfp, dfn
    integer :: ij, k
    h(:,:) = 0.0_DP
    deltal(:,:) = 0.0_DP

    do k = kmin, kmax
       h(:,k) = 3.0_DP / 2.0_DP * f(:,k) &
            - 0.25d0 * ( fh(:,k) + fh(:,k+1) )
    enddo
    !--- CW reconstruction
    do k = kmin+1, kmax-1
       do ij = 1, ijdim
          dfp = h(ij,k+1) - h(ij,k)
          dfn = h(ij,k) - h(ij,k-1)
          dfi = dz(k) / ( dz(k-1) + dz(k) + dz(k+1) ) &
               * ( ( 2.0_DP * dz(k-1) + dz(k) ) / ( dz(k+1) + dz(k) ) &
               * ( h(ij,k+1) - h(ij,k) )  &
               + ( dz(k) + 2.0_DP * dz(k+1) ) / ( dz(k-1) + dz(k) ) &
               * ( h(ij,k) - h(ij,k-1) ) )
          deltal(ij,k) &
               = min( abs(dfi), 2.0_DP * abs(dfp), 2.0_DP * abs(dfn)) &
               * sign(1.0_DP, dfi) &
               * ( sign ( 1.0_DP, dfp * dfn ) + 1.0_DP ) / 2.0_DP &
               / dz(k)
       enddo
    enddo
    !! boundaries: Is this consistent with CW?
    k = kmin
    deltal(:,k) = ( h(:,k+1) - h(:,k) ) / dzh(k+1)
    k = kmax
    deltal(:,k) = ( h(:,k) - h(:,k-1) ) / dzh(k)
    !---
    call sub_fh_cell (ijdim, f, fh, fh_cell, kcell, zdis, deltal)

    return
  end subroutine csl3_cw
  !-----------------------------------------------------------------------------
  subroutine sub_fh_cell (&
       ijdim,             & !--- IN : number of horizontal grids
       f,                 & !--- IN : f at the cell
       fh,                & !--- INOUT : f at the cell wall
       fh_cell,           & !--- OUT : flux at the cell wall
       kcell,             & !--- IN : target cell number
       zdis,              & !--- IN : distance from kcell
       deltal             & !--- IN : df/dz
       )
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz
    implicit none

    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(inout) :: fh(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh_cell(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis(1:ijdim,1:kdim)
    real(DP), intent(in) :: deltal(1:ijdim,1:kdim)
    real(DP) :: fh_old(1:ijdim,1:kdim)
    real(DP) :: sgu, sdz, sdz2, sdz3
    real(DP) :: a, b, c, d
    real(DP) :: zzz
    integer :: ij, k
    integer :: k2
    fh_old(:,: ) = fh(:,:)
    fh_cell(:,:) = 0.0_DP
    do k = kmin, kmax
       do ij = 1, ijdim
          k2 = kcell(ij,k)
          sgu = sign( 1.0_DP, zdis(ij,k) )
          sdz = dz(k2)
          sdz2 = dz(k2) ** 2
          sdz3 = dz(k2) ** 3
          zzz = - zdis(ij,k)
          c =  - sgu * 6.0_DP / sdz * f(ij,k2) &
               + 6.0_DP / sdz &
               * ( fh_old(ij,k2+1) * ( sgu + 1.0_DP ) / 2.0_DP &
               - fh_old(ij,k2)  * ( - sgu + 1.0_DP ) / 2.0_DP ) &
               - 2.0_DP * deltal(ij,k2)
          a =  - 4.0_DP * deltal(ij,k2) / sdz2 &
               + 4.0_DP / sdz3 * ( fh_old(ij,k2+1) - fh_old(ij,k2) )
          b =  - 6.0_DP / sdz2 * f(ij,k2) &
               - sgu * 6.0_DP * deltal(ij,k2) / sdz &
               + 3.0_DP / sdz2 &
               * ( ( 3.0_DP * fh_old(ij,k2+1) - fh_old(ij,k2) ) &
               * ( sgu + 1.0_DP ) / 2.0_DP &
               + ( 3.0_DP * fh_old(ij,k2) - fh_old(ij,k2+1) ) &
               * ( - sgu + 1.0_DP ) / 2.0_DP )

          d = fh_old(ij,k2+1) * ( sgu + 1.0_DP ) / 2.0_DP &
               + fh_old(ij,k2) * ( - sgu + 1.0_DP ) / 2.0_DP
          fh(ij,k) = d + c * zzz + b * zzz ** 2 + a * zzz ** 3

          fh_cell(ij,k) = d + c * zzz / 2.0_DP + b * zzz ** 2 / 3.0_DP &
               + a * zzz ** 3 / 4.0_DP
       enddo
    enddo
    !  fh_cell(:,kmax+1) = 0.0_DP         ! top boundary
    !  fh_cell(:,kmax-1) = fh_cell(:,kmax-2)
    !  fh_cell(:,kmax)   = fh_cell(:,kmax-1)
    !  fh_cell(:,2)    = fh_cell(:,3)
    !  fh_cell(:,1)    = fh_cell(:,2)
    do k = kmin, kmax
       where ( kcell(:,k) == 0 )
          fh_cell(:,k) =  f(:,kmin)
       end where
       where ( kcell(:,k) == kmax+1 )
          fh_cell(:,k) =  f(:,kmax)
       end where
    enddo

    return
  end subroutine sub_fh_cell
  !-----------------------------------------------------------------------------
  subroutine sub_fh_3rd ( &
       ijdim,             & !--- IN : number of horizontal grids
       f,                 & !--- IN : f at the cell
       fh                 & !--- INOUT : f at the cell wall
       )
    use mod_adm, only: &
         kdim => ADM_kall,           &
         kmin => ADM_kmin,           &
         kmax => ADM_kmax
    use mod_grd, only: &
         dz => GRD_dgz
    implicit none

    integer, intent(in) :: ijdim
    real(DP), intent(in) :: f(1:ijdim,1:kdim)
    real(DP), intent(out) :: fh(1:ijdim,1:kdim)
    real(DP) :: deltal(1:ijdim,1:kdim)
    real(DP) :: dfi, dfp, dfn
    integer :: ij, k
    do k = kmin+1, kmax-1
       do ij = 1, ijdim
          dfp = f(ij,k+1) - f(ij,k)
          dfn = f(ij,k) - f(ij,k-1)
          dfi = dz(k) / ( dz(k-1) + dz(k) + dz(k+1) ) &
               * ( ( 2.0_DP * dz(k-1) + dz(k) ) / ( dz(k+1) + dz(k) ) &
               * ( f(ij,k+1) - f(ij,k) )  &
               + ( dz(k) + 2.0_DP * dz(k+1) ) / ( dz(k-1) + dz(k) ) &
               * ( f(ij,k) - f(ij,k-1) ) )
          deltal(ij,k) &
               = min( abs(dfi), 2.0_DP * abs(dfp), 2.0_DP * abs(dfn)) &
               * sign(1.0_DP, dfi) &
               * ( sign ( 1.0_DP, dfp * dfn ) + 1.0_DP ) / 2.0_DP
       enddo
    enddo

    fh(:,:) = 0.0_DP
    do k = kmin+1, kmax-2
       do ij = 1, ijdim
          !      afh = 7.0_DP / 12.0_DP * (f(k+1)+f(k)) - (f(k+2)-f(k-1)) / 12.0_DP
          ! or   afh = 0.5d0 * ( f(k+1) + f(k) ) &
          !          + ( deltal(k) - deltal(k+1) ) / 6.0_DP
          fh(ij,k+1) = f(ij,k) &
               + dz(k) / ( dz(k) + dz(k+1) ) * ( f(ij,k+1) - f(ij,k) ) &
               + 1.0_DP / ( dz(k-1) + dz(k) + dz(k+1) + dz(k+2) ) &
               * ( 2.0_DP * dz(k+1) * dz(k) / ( dz(k) + dz(k+1) ) &
               * ( ( dz(k-1) + dz(k) ) / ( 2.0_DP * dz(k) + dz(k+1) ) &
               - ( dz(k+2) + dz(k+1) ) / ( 2.0_DP * dz(k+1) + dz(k) ) ) &
               * ( f(ij,k+1) - f(ij,k) ) &
               - dz(k) &
               * ( dz(k-1) + dz(k) ) / ( 2.0_DP * dz(k) + dz(k+1) ) &
               * deltal(ij,k+1) &
               + dz(k+1) &
               * ( dz(k+1) + dz(k+2) ) / ( dz(k) + 2.0_DP * dz(k+1) ) &
               * deltal(ij,k) )
       enddo
    enddo
    fh(:,kmin+1) = ( f(:,kmin) + f(:,kmin+1) ) / 2.0_DP
    fh(:,kmax) = ( f(:,kmax-1) + f(:,kmax) ) / 2.0_DP
    fh(:,kmin) = f(:,kmin)
    fh(:,kmax+1) = f(:,kmax)

    return
  end subroutine sub_fh_3rd

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_0s_DP( &
       ijdim, &
       kdim,  &
       kmin,  &
       kmax,  &
       dz,    &
       zh,    &
       dt     )
    use mod_runconf, only: &
       W_VADV1D_MAX, &
       W_VADV1D_MIN
    implicit none

    integer,  intent(in) :: ijdim
    integer,  intent(in) :: kdim
    integer,  intent(in) :: kmin
    integer,  intent(in) :: kmax
    real(DP), intent(in) :: dz(1:ijdim,1:kdim)
    real(DP), intent(in) :: zh(1:ijdim,1:kdim)
    real(DP), intent(in) :: dt

    real(DP) :: zdis(1:ijdim,1:kdim)
    real(DP) :: wh(1:ijdim,1:kdim)          ! w-grid
    integer  :: ij, k
    integer  :: k2
    real(DP) :: zzmax, zzmin
    !---------------------------------------------------------------------------

    allocate( kcell_max_lim(kdim) , kcell_min_lim(kdim) )

    wh(:,:) = W_VADV1D_MIN
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) / 2.0_DP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(:,k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) ) &
            / ( dz(:,k-1) + dz(:,k) ) * 2.0_DP ) ) / 6.0_DP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_DP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(:,k) / 2.0_DP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) / 2.0_DP
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmin = minval(zdis(:,k))
       kcell_max_lim(k) = k
       do k2 = k, kmax
       do ij = 1, ijdim
          if (       zh(ij,k2) <= zh(ij,k) - zzmin  &
               .AND. zh(ij,k) - zzmin < zh(ij,k2+1) ) then
             kcell_max_lim(k) = max( kcell_max_lim(k), k2 )
          endif
       enddo
       enddo
    enddo

    wh(:,:) = W_VADV1D_MAX
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) / 2.0_DP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(:,k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) ) &
            / ( dz(:,k-1) + dz(:,k) ) * 2.0_DP ) ) / 6.0_DP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_DP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(:,k) / 2.0_DP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) / 2.0_DP
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       kcell_min_lim(k) = k
       do k2 = k, kmin, -1
       do ij = 1, ijdim
          if (       zh(ij,k2) <= zh(ij,k) - zzmax  &
               .AND. zh(ij,k) - zzmax < zh(ij,k2+1) ) then
             kcell_min_lim(k) = min( kcell_min_lim(k), k2 )
          endif
       enddo
       enddo
    enddo

    flg_getflux_0s_DP = .true.

    return
  end subroutine vadv1d_getflux_0s_DP
  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_1s_DP( &
       ijdim,            & !--- IN : number of horizontal grids
       kdim,             & !--- IN
       kmin,             & !--- IN
       kmax,             & !--- IN
       dz,               & !--- IN
       zh,               & !--- IN
       wp,               & !--- IN : verical velocity at the full level
       dt                & !--- IN : time interval
       , kcell     & !--- OUT
       , zdis      & !--- OUT
       , kcell_max & !--- OUT
       , kcell_min & !--- OUT
       )
    implicit none

    integer,  intent(in) :: ijdim
    integer,  intent(in) :: kdim
    integer,  intent(in) :: kmin
    integer,  intent(in) :: kmax
    real(DP), intent(in) :: dz(1:ijdim,1:kdim)
    real(DP), intent(in) :: zh(1:ijdim,1:kdim)
    real(DP), intent(in) :: wp(1:ijdim,1:kdim)     ! p-grid
    real(DP), intent(in) :: dt

    integer,  intent(out) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(out) :: zdis(1:ijdim,1:kdim)
    integer,  intent(out) :: kcell_max(1:kdim)
    integer,  intent(out) :: kcell_min(1:kdim)

    real(DP) :: wh(1:ijdim,1:kdim)          ! w-grid
    integer  :: ij, k
    integer  :: k2
    integer  :: k3
    real(DP) :: zzmax, zzmin

    !--- vetical velocity at the half level
    do k = kmin+1, kmax
       wh(:,k) = ( wp(:,k-1) + wp(:,k) ) / 2.0_DP
    enddo
    !------ bottom boundary for wh
    wh(:,kmin) = wp(:,kmin)
    wh(:,kmin-1) = wp(:,kmin-1)
    !------ topboundary for wh : same as inner region
    wh(:,kmax+1) = wp(:,kmax)
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) / 2.0_DP &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(:,k-1) + dz(:,k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(:,k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) ) &
            / ( dz(:,k-1) + dz(:,k) ) * 2.0_DP ) ) / 6.0_DP
    enddo
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.0_DP
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(:,k) / 2.0_DP
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(:,k-1) / 2.0_DP
    !--- calculation of kcell
    !------  top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    ! [fix] 2012/02/01, T.Seiki, to avoid undefined reference
!!$ do k = kmin, kmax+1
    do k  = 1, kdim
    do ij = 1, ijdim
       kcell(ij,k) = k
    enddo
    enddo
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       zzmin = minval(zdis(:,k))
       kcell_min(k) = k
       kcell_max(k) = k
       if ( zzmax > 0.0_DP ) then
          k3 = max(kmin,kcell_min_lim(k))
          do k2 = k, k3, -1
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmax  &
                  .AND. zh(ij,k) - zzmax < zh(ij,k2+1) ) then
                kcell_min(k) = min( kcell_min(k), k2 )
             endif
          enddo
          enddo
       endif
       if ( zzmin < 0.0_DP ) then
          k3 = min(kmax,kcell_max_lim(k))
          do k2 = k, k3
          do ij = 1, ijdim
             if (       zh(ij,k2) <= zh(ij,k) - zzmin  &
                  .AND. zh(ij,k) - zzmin < zh(ij,k2+1) ) then
                kcell_max(k) = max( kcell_max(k), k2 )
             endif
          enddo
          enddo
       endif
    enddo
    !------ determine the kcell at each point.
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          kcell(:,k) = k
       else
          kcell(:,k) = 0
          do k2 = kcell_min(k), kcell_max(k)
             kcell(:,k) = max ( kcell(:,k), &
                  int ( k2 &
                  * sign ( 1.0_DP, ( zh(:,k) - zdis(:,k) ) - zh(:,k2) ) &
                  * sign ( 1.0_DP, zh(:,k2+1) - ( zh(:,k) - zdis(:,k) ) ) ) )
          enddo
       endif
    enddo

    return
  end subroutine vadv1d_getflux_1s_DP
  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_2s_DP( &
       ijdim,            & !--- IN : number of horizontal grids
       kdim,             & !--- IN
       kmin,             & !--- IN
       kmax,             & !--- IN
       dz,               & !--- IN
       rhof,             & !--- IN : rho X f X gam2
       frhof             & !--- OUT : vertical flux ( X dt )
       , kcell     & !--- in
       , zdis_     & !--- in
       , kcell_max & !--- in
       , kcell_min & !--- in
       )
    implicit none

    integer,  intent(in) :: ijdim
    integer,  intent(in) :: kdim
    integer,  intent(in) :: kmin
    integer,  intent(in) :: kmax
    real(DP), intent(in) :: dz(1:ijdim,1:kdim)
    real(DP), intent(in) :: rhof(1:ijdim,1:kdim)   ! p-grid
    real(DP), intent(out) :: frhof(1:ijdim,1:kdim) ! w-grid

    integer,  intent(in) :: kcell(1:ijdim,1:kdim)
    real(DP), intent(in) :: zdis_(1:ijdim,1:kdim)
    integer,  intent(in) :: kcell_max(1:kdim)
    integer,  intent(in) :: kcell_min(1:kdim)

    integer  :: ij, k, kk
    integer  :: k2
    real(DP) :: zdis(1:ijdim,1:kdim)  ! w-grid

    zdis(:,:) = zdis_(:,:)
    !--- set zero value for norain, see if block in do loop
    frhof(:,:) = 0.0_DP
    !------ integration in the integer cells
    do k = kmin, kmax
       if ( kcell_min(k) == k .AND. kcell_max(k) == k ) then
          ! frhof(:,:) = 0.0_DP
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             frhof(:,k) &
                  = frhof(:,k) &
                  + rhof(:,k2) * dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
             zdis(:,k) = zdis(:,k) - dz(:,k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.0_DP ) / 2.0_DP &
                  + ( - 1.0_DP ) * ( sign ( 1, k2 - k ) + 1.0_DP ) / 2.0_DP &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.0_DP ) / 2.0_DP )
          enddo
       endif
    enddo

    ! top boundary: rigid
    !  zdis(:,kmax+1) = 0.0_SP
    !  frhof(:,kmax+1) = 0.0_SP

    !--- integration in the non-integer cell ( depending on the scheme below )
    do k  = kmin, kmax
    do ij = 1, ijdim
       kk = max( min( kcell(ij,k), kmax ), kmin )

       frhof(ij,k) = frhof(ij,k) + rhof(ij,kk) * zdis(ij,k)
    enddo
    enddo

    return
  end subroutine vadv1d_getflux_2s_DP

end module mod_vadv1d
