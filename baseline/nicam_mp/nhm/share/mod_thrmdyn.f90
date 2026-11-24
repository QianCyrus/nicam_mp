!-------------------------------------------------------------------------------
!> Module thermodyanics
!!
!! @par Description
!!          This module calculates thermodyanical variables
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_thrmdyn
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
  logical, public :: acc_wait_on = .true.
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: THRMDYN_qd
  public :: THRMDYN_cv
  public :: THRMDYN_cp
  public :: THRMDYN_rho
  public :: THRMDYN_ein
  public :: THRMDYN_th
  public :: THRMDYN_eth
  public :: THRMDYN_ent
  public :: THRMDYN_rhoein
  public :: THRMDYN_tempre
  public :: THRMDYN_phydro

  interface THRMDYN_qd
     module procedure THRMDYN_qd_ijk_SP
     module procedure THRMDYN_qd_ijk_DP
     module procedure THRMDYN_qd_ijkl_SP
     module procedure THRMDYN_qd_ijkl_DP
  end interface THRMDYN_qd

  interface THRMDYN_cv
     module procedure THRMDYN_cv_ijk_SP
     module procedure THRMDYN_cv_ijk_DP
  end interface THRMDYN_cv

  interface THRMDYN_cp
     module procedure THRMDYN_cp_ijk_SP
     module procedure THRMDYN_cp_ijk_DP
  end interface THRMDYN_cp

  interface THRMDYN_rho
     module procedure THRMDYN_rho_ijk_SP
     module procedure THRMDYN_rho_ijk_DP
  end interface THRMDYN_rho

  interface THRMDYN_ein
     module procedure THRMDYN_ein_ijk_SP
     module procedure THRMDYN_ein_ijk_DP
  end interface THRMDYN_ein

  interface THRMDYN_th
     module procedure THRMDYN_th_ijk_SP
     module procedure THRMDYN_th_ijk_DP
     module procedure THRMDYN_th_ijkl_SP
     module procedure THRMDYN_th_ijkl_DP
  end interface THRMDYN_th

  interface THRMDYN_eth
     module procedure THRMDYN_eth_ijkl_SP
     module procedure THRMDYN_eth_ijkl_DP
  end interface THRMDYN_eth

  interface THRMDYN_ent
     module procedure THRMDYN_ent_ijk_SP
     module procedure THRMDYN_ent_ijk_DP
  end interface THRMDYN_ent

  interface THRMDYN_rhoein
     module procedure THRMDYN_rhoein_ijk_SP
     module procedure THRMDYN_rhoein_ijk_DP
     module procedure THRMDYN_rhoein_ijkl_SP
     module procedure THRMDYN_rhoein_ijkl_DP
  end interface THRMDYN_rhoein

  interface THRMDYN_tempre
     module procedure THRMDYN_tempre_ijk_SP
     module procedure THRMDYN_tempre_ijk_DP
     module procedure THRMDYN_tempre_ijkl
  end interface THRMDYN_tempre

  interface THRMDYN_phydro
     module procedure THRMDYN_phydro_ijkl_SP
     module procedure THRMDYN_phydro_ijkl_DP
  end interface THRMDYN_phydro

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijk_SP( &
       ijdim, &
       kdim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

#ifdef _OPENACC
    !!$acc data pcopyin(q) if(acc_is_present(qd))

    !$acc kernels async(1) pcopy(qd) pcopyin(q) if(acc_is_present(qd))
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_SP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,qd,q)

!OCL XFILL,NOPREFETCH
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

!OCL NOPREFETCH
    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_qd_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijk_DP( &
       ijdim, &
       kdim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

#ifdef _OPENACC
    !!$acc data pcopyin(q) if(acc_is_present(qd))

    !$acc kernels async(1) pcopy(qd) pcopyin(q) if(acc_is_present(qd))
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_DP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,qd,q)

!OCL XFILL,NOPREFETCH
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

!OCL NOPREFETCH
    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_qd_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: q (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: qd(ijdim,kdim,ldim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

#ifdef _OPENACC
    !!$acc data pcopyin(q) if(acc_is_present(qd))

    !$acc kernels async(1) pcopy(qd) pcopyin(q) if(acc_is_present(qd))
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_SP

       !$acc loop seq
       do nq = NQW_STR,NQW_END
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo

    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,qd,q)

!OCL XFILL,NOPREFETCH
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_SP
    enddo
    enddo
    enddo
    !$omp end do

!OCL NOPREFETCH
    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_qd_ijkl_SP

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       q,     &
       qd     )
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: q (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: qd(ijdim,kdim,ldim)       ! dry air mass concentration [kg/kg]

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

#ifdef _OPENACC
    !!$acc data pcopyin(q) if(acc_is_present(qd))

    !$acc kernels async(1) pcopy(qd) pcopyin(q) if(acc_is_present(qd))
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_DP

       !$acc loop seq
       do nq = NQW_STR,NQW_END
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo

    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,qd,q)

!OCL XFILL,NOPREFETCH
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_DP
    enddo
    enddo
    enddo
    !$omp end do

!OCL NOPREFETCH
    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_qd_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cv_ijk_SP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cv     )
    use mod_cnst, only: &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: cv(ijdim,kdim)       ! specific heat [J/kg/K]

    real(SP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry

#ifdef _OPENACC
    !!$acc data pcopyin(qd,q,CVW) if(acc_is_present(cv))

    !$acc kernels async(1) pcopy(cv) pcopyin(qd,q,CVW) if(acc_is_present(cv))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo

    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_cv_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cv_ijk_DP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cv     )
    use mod_cnst, only: &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: cv(ijdim,kdim)       ! specific heat [J/kg/K]

    real(DP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry

#ifdef _OPENACC
    !!$acc data pcopyin(qd,q,CVW) if(acc_is_present(cv))

    !$acc kernels async(1) pcopy(cv) pcopyin(qd,q,CVW) if(acc_is_present(cv))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo

    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_cv_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cp_ijk_SP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cp     )
    use mod_cnst, only: &
       CNST_CPdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CPW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: cp(ijdim,kdim)       ! specific heat [J/kg/K]

    real(SP) :: CPdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry = CNST_CPdry

#ifdef _OPENACC
    !!$acc data pcopyin(qd,q,CPW) if(acc_is_present(cp))

    !$acc kernels async(1) pcopy(cp) pcopyin(qd,q,CPW) if(acc_is_present(cp))
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo

    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cp,qd,q,CPW,CPdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_cp_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cp_ijk_DP( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cp     )
    use mod_cnst, only: &
       CNST_CPdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CPW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: cp(ijdim,kdim)       ! specific heat [J/kg/K]

    real(DP) :: CPdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry = CNST_CPdry

#ifdef _OPENACC
    !!$acc data pcopyin(qd,q,CPW) if(acc_is_present(cp))

    !$acc kernels async(1) pcopy(cp) pcopyin(qd,q,CPW) if(acc_is_present(cp))
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo

    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,cp,qd,q,CPW,CPdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    return
  end subroutine THRMDYN_cp_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(SP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim)       ! density [kg/m3]

    real(SP) :: Rdry, Rvap

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(pre,tem,qd,q) if(acc_is_present(rho))

    Rdry = CNST_Rdry
    Rvap = CNST_Rvap

!OCL XFILL
    !$acc kernels async(1) pcopy(rho) pcopyin(pre,tem,qd,q) if(acc_is_present(rho))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,pre,tem,qd,q,Rdry,Rvap,I_QV)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_rho_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       I_QV
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(DP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim)       ! density [kg/m3]

    real(DP) :: Rdry, Rvap

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(pre,tem,qd,q) if(acc_is_present(rho))

    Rdry = CNST_Rdry
    Rvap = CNST_Rvap

!OCL XFILL
    !$acc kernels async(1) pcopy(rho) pcopyin(pre,tem,qd,q) if(acc_is_present(rho))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,pre,tem,qd,q,Rdry,Rvap,I_QV)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_rho_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate internal energy
  subroutine THRMDYN_ein_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       qd,    &
       q,     &
       ein    )
    use mod_cnst, only: &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry

#ifdef _OPENACC
    !$acc data pcreate(cv) pcopyin(tem,qd,q,CVW) if(acc_is_present(ein))

    !$acc kernels async(1) pcopy(ein,cv) pcopyin(tem,qd,q,CVW) if(acc_is_present(ein))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo

       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,ein,tem,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_ein_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate internal energy
  subroutine THRMDYN_ein_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       qd,    &
       q,     &
       ein    )
    use mod_cnst, only: &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: CVdry

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry

#ifdef _OPENACC
    !$acc data pcreate(cv) pcopyin(tem,qd,q,CVW) if(acc_is_present(ein))

    !$acc kernels async(1) pcopy(ein,cv) pcopyin(tem,qd,q,CVW) if(acc_is_present(ein))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo

       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,ein,tem,cv,qd,q,CVW,CVdry)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       enddo
       !$omp end do
    enddo

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_ein_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       th     )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(SP), intent(out) :: th (ijdim,kdim) ! potential temperature [K]

    real(SP) :: RovCP, PRE00

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem,pre) if(acc_is_present(th))

    RovCP = CNST_Rdry / CNST_CPdry
    PRE00 = CNST_PRE00

!OCL XFILL,SWP
    !$acc kernels async(1) pcopy(th) pcopyin(tem,pre) if(acc_is_present(th))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,th,tem,pre,RovCP,PRE00)
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k) = tem(ij,k) * ( PRE00 / pre(ij,k) )**RovCP
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_th_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       th     )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(DP), intent(out) :: th (ijdim,kdim) ! potential temperature [K]

    real(DP) :: RovCP, PRE00

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem,pre) if(acc_is_present(th))

    RovCP = CNST_Rdry / CNST_CPdry
    PRE00 = CNST_PRE00

!OCL XFILL,SWP
    !$acc kernels async(1) pcopy(th) pcopyin(tem,pre) if(acc_is_present(th))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,th,tem,pre,RovCP,PRE00)
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k) = tem(ij,k) * ( PRE00 / pre(ij,k) )**RovCP
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_th_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       th     )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: tem(ijdim,kdim,ldim) ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(SP), intent(out) :: th (ijdim,kdim,ldim) ! potential temperature [K]

    real(SP) :: RovCP, PRE00

    integer  :: ij, k, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem,pre) if(acc_is_present(th))

    RovCP = CNST_Rdry / CNST_CPdry
    PRE00 = CNST_PRE00

!OCL XFILL,SWP
    !$acc kernels async(1) pcopy(th) pcopyin(tem,pre) if(acc_is_present(th))
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,th,tem,pre,RovCP,PRE00), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k,l) = tem(ij,k,l) * ( PRE00 / pre(ij,k,l) )**RovCP
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_th_ijkl_SP

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       th     )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: tem(ijdim,kdim,ldim) ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(DP), intent(out) :: th (ijdim,kdim,ldim) ! potential temperature [K]

    real(DP) :: RovCP, PRE00

    integer  :: ij, k, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(tem,pre) if(acc_is_present(th))

    RovCP = CNST_Rdry / CNST_CPdry
    PRE00 = CNST_PRE00

!OCL XFILL,SWP
    !$acc kernels async(1) pcopy(th) pcopyin(tem,pre) if(acc_is_present(th))
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,th,tem,pre,RovCP,PRE00), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k,l) = tem(ij,k,l) * ( PRE00 / pre(ij,k,l) )**RovCP
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_th_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: ein(ijdim,kdim,ldim) ! internal energy [J]
    real(SP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure        [Pa]
    real(SP), intent(in)  :: rho(ijdim,kdim,ldim) ! density         [kg/m3]
    real(SP), intent(out) :: eth(ijdim,kdim,ldim) ! enthalpy

    integer  :: ij, k, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(ein,pre,rho) if(acc_is_present(eth))

!OCL XFILL
    !$acc kernels async(1) pcopy(eth) pcopyin(ein,pre,rho) if(acc_is_present(eth))
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,eth,ein,pre,rho), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k,l) = ein(ij,k,l) + pre(ij,k,l) / rho(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_eth_ijkl_SP

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: ein(ijdim,kdim,ldim) ! internal energy [J]
    real(DP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure        [Pa]
    real(DP), intent(in)  :: rho(ijdim,kdim,ldim) ! density         [kg/m3]
    real(DP), intent(out) :: eth(ijdim,kdim,ldim) ! enthalpy

    integer  :: ij, k, l
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(ein,pre,rho) if(acc_is_present(eth))

!OCL XFILL
    !$acc kernels async(1) pcopy(eth) pcopyin(ein,pre,rho) if(acc_is_present(eth))
    !$omp parallel do default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,eth,ein,pre,rho), &
    !$omp collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k,l) = ein(ij,k,l) + pre(ij,k,l) / rho(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif

    !!$acc end data
    return
  end subroutine THRMDYN_eth_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate entropy
  subroutine THRMDYN_ent_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       ent    )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_Rvap,  &
       CNST_LHV,   &
       CNST_LHF,   &
       CNST_PRE00, &
       CNST_TEM00, &
       CNST_PSAT0, &
       CNST_EPSvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)
    real(SP), intent(in)  :: pre(ijdim,kdim)
    real(SP), intent(in)  :: qd (ijdim,kdim)
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax)
    real(SP), intent(out) :: ent(ijdim,kdim)

    real(SP) :: Pdry
    real(SP) :: Pvap
    real(SP) :: LH(nqmax)
    real(SP) :: CPdry, Rdry, Rvap, TEM00, PRE00, PSAT0, EPSvap

    real(SP), parameter :: EPS = 1.E-5_SP

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry  = CNST_CPdry
    Rdry   = CNST_Rdry
    Rvap   = CNST_Rvap
    TEM00  = CNST_TEM00
    PRE00  = CNST_PRE00
    PSAT0  = CNST_PSAT0
    EPSvap = CNST_EPSvap

    do nq = NQW_STR, NQW_END
       if ( nq == I_QV ) then
          LH(nq) =  CNST_LHV / TEM00
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          LH(nq) = -CNST_LHF / TEM00
       else
          LH(nq) = 0.0_SP
       endif
    enddo

#ifdef _OPENACC
    !$acc data pcopyin(tem,pre,qd,q,CVW,LH) if(acc_is_present(ent))

    !$acc kernels async(1) pcopy(ent) pcopyin(tem,pre,qd,q) if(acc_is_present(ent))
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1) pcopy(ent) pcopyin(tem,q,CVW,LH) if(acc_is_present(ent))
    do k  = 1, kdim
    do ij = 1, ijdim
       !$acc loop seq
       do nq = NQW_STR, NQW_END
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,Pdry,Pvap), &
    !$omp shared(ijdim,kdim,nq,NQW_STR,NQW_END,ent,tem,pre,qd,q,     &
    !$omp        LH,CVW,CPdry,Rdry,Rvap,TEM00,PRE00,PSAT0,EPSvap,I_QV)

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
       enddo
       !$omp end do
    enddo

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_ent_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate entropy
  subroutine THRMDYN_ent_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       ent    )
    use mod_cnst, only: &
       CNST_Rdry,  &
       CNST_CPdry, &
       CNST_Rvap,  &
       CNST_LHV,   &
       CNST_LHF,   &
       CNST_PRE00, &
       CNST_TEM00, &
       CNST_PSAT0, &
       CNST_EPSvap
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)
    real(DP), intent(in)  :: pre(ijdim,kdim)
    real(DP), intent(in)  :: qd (ijdim,kdim)
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax)
    real(DP), intent(out) :: ent(ijdim,kdim)

    real(DP) :: Pdry
    real(DP) :: Pvap
    real(DP) :: LH(nqmax)
    real(DP) :: CPdry, Rdry, Rvap, TEM00, PRE00, PSAT0, EPSvap

    real(DP), parameter :: EPS = 1.E-10_DP

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CPdry  = CNST_CPdry
    Rdry   = CNST_Rdry
    Rvap   = CNST_Rvap
    TEM00  = CNST_TEM00
    PRE00  = CNST_PRE00
    PSAT0  = CNST_PSAT0
    EPSvap = CNST_EPSvap

    do nq = NQW_STR, NQW_END
       if ( nq == I_QV ) then
          LH(nq) =  CNST_LHV / TEM00
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          LH(nq) = -CNST_LHF / TEM00
       else
          LH(nq) = 0.0_DP
       endif
    enddo

#ifdef _OPENACC
    !$acc data pcopyin(tem,pre,qd,q,CVW,LH) if(acc_is_present(ent))

    !$acc kernels async(1) pcopy(ent) pcopyin(tem,pre,qd,q) if(acc_is_present(ent))
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1) pcopy(ent) pcopyin(tem,q,CVW,LH) if(acc_is_present(ent))
    do k  = 1, kdim
    do ij = 1, ijdim
       !$acc loop seq
       do nq = NQW_STR, NQW_END
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,Pdry,Pvap), &
    !$omp shared(ijdim,kdim,nq,NQW_STR,NQW_END,ent,tem,pre,qd,q,     &
    !$omp        LH,CVW,CPdry,Rdry,Rvap,TEM00,PRE00,PSAT0,EPSvap,I_QV)

    !$acc kernels async(1) pcopy(ent) pcopyin(tem,pre,qd,q)
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSvap*qd(ij,k) / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)    / ( EPSvap*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$omp end do
    !$acc end kernels

    do nq = NQW_STR, NQW_END
       !$acc kernels async(1) pcopy(ent) pcopyin(tem,q,CVW,LH)
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
       enddo
       !$omp end do
       !$acc end kernels
    enddo

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_ent_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijk_SP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(SP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: qd(ijdim,kdim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))

    !$acc kernels async(1) pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo

       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_rhoein_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijk_DP( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(DP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: qd(ijdim,kdim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))

    !$acc kernels async(1) pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo

       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = cv(ij,k) + qd(ij,k) * CVdry

       rho(ij,k) = pre(ij,k) / ( ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k) )
       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_rhoein_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(SP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(SP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(SP), intent(out) :: ein(ijdim,kdim,ldim)       ! internal energy [J]

    real(SP) :: cv(ijdim,kdim,ldim)
    real(SP) :: qd(ijdim,kdim,ldim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))

    !$acc kernels async(1) pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_SP
       qd(ij,k,l) = 1.0_SP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo

       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_SP
       qd(ij,k,l) = 1.0_SP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_rhoein_ijkl_SP

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(DP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(DP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(DP), intent(out) :: ein(ijdim,kdim,ldim)       ! internal energy [J]

    real(DP) :: cv(ijdim,kdim,ldim)
    real(DP) :: qd(ijdim,kdim,ldim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))

    !$acc kernels async(1) pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW) if(acc_is_present(rho))
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_DP
       qd(ij,k,l) = 1.0_DP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo

       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,rho,ein,tem,pre,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_DP
       qd(ij,k,l) = 1.0_DP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / ( ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l) )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_rhoein_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijk_SP( &
       ijdim, &
       kdim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(SP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(SP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]
    real(SP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    real(SP) :: cv(ijdim,kdim)
    real(SP) :: qd(ijdim,kdim)
    real(SP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))

    !$acc kernels async(1) pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo

       cv (ij,k) = cv (ij,k) + qd(ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv(ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL,NOPREFETCH
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_SP
       qd(ij,k) = 1.0_SP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k) = cv (ij,k) + qd (ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv (ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_tempre_ijk_SP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijk_DP( &
       ijdim, &
       kdim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(DP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(DP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]
    real(DP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    real(DP) :: cv(ijdim,kdim)
    real(DP) :: qd(ijdim,kdim)
    real(DP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))

    !$acc kernels async(1) pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo

       cv (ij,k) = cv (ij,k) + qd(ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv(ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL,NOPREFETCH
    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = 0.0_DP
       qd(ij,k) = 1.0_DP
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k) = cv (ij,k) + qd (ij,k) * CVdry
       tem(ij,k) = ein(ij,k) / cv (ij,k)
       pre(ij,k) = rho(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_tempre_ijk_DP

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    use mod_cnst, only: &
       CNST_Rdry, &
       CNST_Rvap, &
       CNST_CVdry
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QV,              &
       CVW
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(RP), intent(in)  :: ein(ijdim,kdim,ldim)       ! internal energy [J]
    real(RP), intent(in)  :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(RP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(RP), intent(out) :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]

    real(RP) :: cv(ijdim,kdim,ldim)
    real(RP) :: qd(ijdim,kdim,ldim)
    real(RP) :: CVdry, Rdry, Rvap

    integer  :: ij, k, l, nq
    !---------------------------------------------------------------------------

    CVdry = CNST_CVdry
    Rdry  = CNST_Rdry
    Rvap  = CNST_Rvap

#ifdef _OPENACC
    !$acc data pcreate(cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))

    !$acc kernels async(1) pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW) if(acc_is_present(tem))
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_RP
       qd(ij,k,l) = 1.0_RP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo

       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
       pre(ij,k,l) = rho(ij,k,l) * ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data
#else
    !$omp parallel default(none),private(ij,k,l,nq), &
    !$omp shared(ijdim,kdim,ldim,NQW_STR,NQW_END,tem,pre,ein,rho,q,cv,qd,CVW,CVdry,Rdry,Rvap,I_QV)

!OCL XFILL,NOPREFETCH
    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_RP
       qd(ij,k,l) = 1.0_RP
    enddo
    enddo
    enddo
    !$omp end do

    do nq = NQW_STR, NQW_END
       !$omp do collapse(2)
       do l  = 1, ldim
       do k  = 1, kdim
       do ij = 1, ijdim
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       enddo
       enddo
       !$omp end do
    enddo

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv (ij,k,l) = cv (ij,k,l) + qd (ij,k,l) * CVdry
       tem(ij,k,l) = ein(ij,k,l) / cv (ij,k,l)
       pre(ij,k,l) = rho(ij,k,l) * ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap ) * tem(ij,k,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel
#endif

    !$acc wait(1)

    return
  end subroutine THRMDYN_tempre_ijkl

  !-----------------------------------------------------------------------------
  !> calculate hydrostatic pressure
  subroutine THRMDYN_phydro_ijkl_DP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       z,     &
       zh,    &
       phydro )
    use mod_cnst, only: &
       CNST_EPS,  &
       CNST_GRAV
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(DP), intent(in)  :: tem   (ijdim,kdim,ldim)       ! temperature [K]
    real(DP), intent(in)  :: pre   (ijdim,kdim,ldim)       ! pressure    [Pa]
    real(DP), intent(in)  :: q     (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(DP), intent(in)  :: z     (ijdim,kdim,ldim)       ! altitude (cell center) [m]
    real(DP), intent(in)  :: zh    (ijdim,kdim,ldim)       ! altitude (cell face)   [m]
    real(DP), intent(out) :: phydro(ijdim,kdim,ldim)       ! hydrostatic pressure [Pa]

    real(DP) :: rho    (ijdim,kdim,ldim) ! density     [kg/m3]
    real(DP) :: ein    (ijdim,kdim,ldim) ! internal energy [J]
    real(DP) :: phydroh(ijdim,kdim,ldim) ! hydrostatic pressure (cell face) [Pa]
    real(DP) :: diff   (ijdim,ldim)
    real(DP) :: EPS, GRAV

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    call THRMDYN_rhoein( ijdim,         & ! [IN]
                         kdim,          & ! [IN]
                         ldim,          & ! [IN]
                         tem (:,:,:),   & ! [IN]
                         pre (:,:,:),   & ! [IN]
                         q   (:,:,:,:), & ! [IN]
                         rho (:,:,:),   & ! [OUT]
                         ein (:,:,:)    ) ! [OUT]

    EPS  = CNST_EPS
    GRAV = CNST_GRAV

    !$omp parallel default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,phydroh,phydro,rho,pre,z,zh,diff,EPS,GRAV)

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       phydroh(ij,kdim  ,l) = pre    (ij,kdim-1,l) - rho(ij,kdim-1,l) * GRAV * ( zh(ij,kdim,l) - z (ij,kdim-1,l) )
       phydroh(ij,kdim-1,l) = phydroh(ij,kdim  ,l) + rho(ij,kdim-1,l) * GRAV * ( zh(ij,kdim,l) - zh(ij,kdim-1,l) )
       phydro (ij,kdim-1,l) = pre    (ij,kdim-1,l)

       diff(ij,l) = 0.0_DP
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do k  = kdim-2, 2, -1
    do ij = 1, ijdim
       phydroh(ij,k,l) = phydroh(ij,k+1,l) + rho(ij,k,l) * GRAV * ( zh(ij,k+1,l) - zh(ij,k,l) )
       phydro (ij,k,l) = phydro (ij,k+1,l) + 0.5_DP * ( rho(ij,k+1,l) + rho(ij,k,l) ) * GRAV * ( z(ij,k+1,l) - z(ij,k,l) )

       diff(ij,l) = diff(ij,l) + ( pre(ij,k,l) - phydro(ij,k,l) ) * ( zh(ij,k+1,l) - zh(ij,k,l) )
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       phydroh(ij,1   ,l) = phydroh(ij,2     ,l) + rho(ij,1,l) * GRAV * ( zh(ij,2,l) - zh(ij,1,l) )
       phydro (ij,1   ,l) = phydro (ij,2     ,l) + 0.5_DP * ( rho(ij,2   ,l) + rho(ij,1     ,l) ) &
                                                 * GRAV   * ( z  (ij,2   ,l) - z  (ij,1     ,l) )
       phydro (ij,kdim,l) = phydro (ij,kdim-1,l) - 0.5_DP * ( rho(ij,kdim,l) + rho(ij,kdim-1,l) ) &
                                                 * GRAV   * ( z  (ij,kdim,l) - z  (ij,kdim-1,l) )
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       diff(ij,l) = diff(ij,l) / ( zh(ij,kdim,l) - zh(ij,2,l) )
       diff(ij,l) = max( diff(ij,l), EPS - phydroh(ij,kdim,l) )
    enddo
    enddo
    !$omp end do

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       phydroh(ij,k,l) = phydroh(ij,k,l) + diff(ij,l)
       phydro (ij,k,l) = phydro (ij,k,l) + diff(ij,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    !$acc wait(1)

    return
  end subroutine THRMDYN_phydro_ijkl_DP

  !-----------------------------------------------------------------------------
  !> calculate hydrostatic pressure
  subroutine THRMDYN_phydro_ijkl_SP( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       z,     &
       zh,    &
       phydro )
    use mod_cnst, only: &
       CNST_EPS,  &
       CNST_GRAV
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: ldim
    real(SP), intent(in)  :: tem   (ijdim,kdim,ldim)       ! temperature [K]
    real(SP), intent(in)  :: pre   (ijdim,kdim,ldim)       ! pressure    [Pa]
    real(SP), intent(in)  :: q     (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(SP), intent(in)  :: z     (ijdim,kdim,ldim)       ! altitude (cell center) [m]
    real(SP), intent(in)  :: zh    (ijdim,kdim,ldim)       ! altitude (cell face)   [m]
    real(SP), intent(out) :: phydro(ijdim,kdim,ldim)       ! hydrostatic pressure [Pa]

    real(SP) :: rho    (ijdim,kdim,ldim) ! density     [kg/m3]
    real(SP) :: ein    (ijdim,kdim,ldim) ! internal energy [J]
    real(SP) :: phydroh(ijdim,kdim,ldim) ! hydrostatic pressure (cell face) [Pa]
    real(SP) :: diff   (ijdim,ldim)
    real(SP) :: EPS, GRAV

    integer  :: ij, k, l
    !---------------------------------------------------------------------------

    call THRMDYN_rhoein( ijdim,         & ! [IN]
                         kdim,          & ! [IN]
                         ldim,          & ! [IN]
                         tem (:,:,:),   & ! [IN]
                         pre (:,:,:),   & ! [IN]
                         q   (:,:,:,:), & ! [IN]
                         rho (:,:,:),   & ! [OUT]
                         ein (:,:,:)    ) ! [OUT]

    EPS  = CNST_EPS
    GRAV = CNST_GRAV

    !$omp parallel default(none),private(ij,k,l), &
    !$omp shared(ijdim,kdim,ldim,phydroh,phydro,rho,pre,z,zh,diff,EPS,GRAV)

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       phydroh(ij,kdim  ,l) = pre    (ij,kdim-1,l) - rho(ij,kdim-1,l) * GRAV * ( zh(ij,kdim,l) - z (ij,kdim-1,l) )
       phydroh(ij,kdim-1,l) = phydroh(ij,kdim  ,l) + rho(ij,kdim-1,l) * GRAV * ( zh(ij,kdim,l) - zh(ij,kdim-1,l) )
       phydro (ij,kdim-1,l) = pre    (ij,kdim-1,l)

       diff(ij,l) = 0.0_SP
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do k  = kdim-2, 2, -1
    do ij = 1, ijdim
       phydroh(ij,k,l) = phydroh(ij,k+1,l) + rho(ij,k,l) * GRAV * ( zh(ij,k+1,l) - zh(ij,k,l) )
       phydro (ij,k,l) = phydro (ij,k+1,l) + 0.5_SP * ( rho(ij,k+1,l) + rho(ij,k,l) ) * GRAV * ( z(ij,k+1,l) - z(ij,k,l) )

       diff(ij,l) = diff(ij,l) + ( pre(ij,k,l) - phydro(ij,k,l) ) * ( zh(ij,k+1,l) - zh(ij,k,l) )
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       phydroh(ij,1   ,l) = phydroh(ij,2     ,l) + rho(ij,1,l) * GRAV * ( zh(ij,2,l) - zh(ij,1,l) )
       phydro (ij,1   ,l) = phydro (ij,2     ,l) + 0.5_SP * ( rho(ij,2   ,l) + rho(ij,1     ,l) ) &
                                                 * GRAV   * ( z  (ij,2   ,l) - z  (ij,1     ,l) )
       phydro (ij,kdim,l) = phydro (ij,kdim-1,l) - 0.5_SP * ( rho(ij,kdim,l) + rho(ij,kdim-1,l) ) &
                                                 * GRAV   * ( z  (ij,kdim,l) - z  (ij,kdim-1,l) )
    enddo
    enddo
    !$omp end do

    !$omp do
    do l  = 1, ldim
    do ij = 1, ijdim
       diff(ij,l) = diff(ij,l) / ( zh(ij,kdim,l) - zh(ij,2,l) )
       diff(ij,l) = max( diff(ij,l), EPS - phydroh(ij,kdim,l) )
    enddo
    enddo
    !$omp end do

    !$omp do collapse(2)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       phydroh(ij,k,l) = phydroh(ij,k,l) + diff(ij,l)
       phydro (ij,k,l) = phydro (ij,k,l) + diff(ij,l)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    !$acc wait(1)

    return
  end subroutine THRMDYN_phydro_ijkl_SP

end module mod_thrmdyn
