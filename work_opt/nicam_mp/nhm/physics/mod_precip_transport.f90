!-------------------------------------------------------------------------------
!
!+  precipitation transport module
!
!-------------------------------------------------------------------------------
module mod_precip_transport
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module contains subroutines w.r.t. precipitation transport
  !
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.03                 Imported from Satoh's nonhydro model.
  !                04-12-23   change dimension of precip (ijdim,2)
  !                           add output var: precip_energy
  !                05-10-07   M.Satoh: add precip_transport_2water
  !                06-09-25   Bug fix for drop of momentum with rain. (H.Tomita)
  !                07-01-26   H.Tomita :
  !                           Support the all EIN_TYPEs including 'SIMPLE2'.
  !                           Add  sub[precip_transport_1water].
  !                07-11-26   S.Iga :
  !                           Bug fix in precip_transport_1water.
  !                           (Terms including fprecr_rhoe are deleted)
  !                08-03-06   H.Tomita: Add precip_transport_2water_d for WSM3
  !                08-04-12   T.Mitsui: Add precip_transport_4wq_4wn_inc
  !                                     qr,qi,qs,qg + nr,ni,ns,ng + soluble tracers
  !                08-05-24   T.Mitsui: Add precip_transport_5wq_5wn_inc
  !                                     qc,qr,qi,qs,qg + nc,nr,ni,ns,ng + soluble tracers
  !                08-05-30   T.Mitsui: bug-fix in precip_transport_5wq_5wn_inc
  !                09-04-14   T.Mitsui: del incloud scavenging from 5wq_5wn_inc, 4wq_4wn_inc
  !                09-07-10   H.Tomita: Integrate the precip_transport routines
  !                                     to 'precip_transport_nwater'.
  !                09-07-14   C.Kodama: Bug fix
  !                09-07-28   H.Tomita: Bug fix for energy correction
  !                                     and suffix of vad1d_getflux
  !                09-08-03   S.Iga     Bug fix ( initial value definition of drhokin_h_prec)
  !                10-08-03   T.Mitsui: application for passive tracer(assumption is followings)
  !                                     if NQW_STR <= nq <= NQW_END then nq are active tracers.
  !                                     else nq are passive tracers.
  !                                     (I've forgotten to commit this module)
  !                12-02-01   T.Seiki : add output of aerosol wash-out by optional variables
  !                12-05-18   H.Yashiro : Speedup Flux-Semi-Lagrangian & Simple Upstream-Euler
  !                12-06-06   H.Yashiro : bugfix
  !                12-09-24   FJ : Tuned precip_transport_new
  !                14-12-22   C.Kodama : bug fix
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
  public :: precip_transport
  public :: precip_transport_DP

  interface precip_transport
     module procedure precip_transport_SP
     module procedure precip_transport_DP
  end interface precip_transport

  public :: precip_transport_new
  public :: precip_transport_euler

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine precip_transport_SP( &
       ijdim,              &
       rhog,               &
       rhogvx,             &
       rhogvy,             &
       rhogvz,             &
       rhogw,              &
       rhoge,              &
       rhogq,              &
       rho,                &
       tem,                &
       pre,                &
       vx,                 &
       vy,                 &
       vz,                 &
       w,                  &
       q,                  &
       qd,                 &
       z,                  &
       zh,                 &
       WPREC,              &
       precipitating_flag, &
       precip,             &
       precip_rhoe,        &
       precip_lh_heat,     &
       precip_rhophi,      &
       precip_rhokin,      &
       precip_flux,        &
       gsgam2,             &
       gsgam2h,            &
       rgs,                &
       rgsh,               &
       ix,                 &
       iy,                 &
       iz,                 &
       jx,                 &
       jy,                 &
       jz,                 &
       dt,                 &
       precip_trc          )
    use mod_adm, only: &
       ADM_PRC_ME,       &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR, NQW_END,  &
       I_QC,              &
       I_QR,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW,               &
       LHF,               &
       PRCIP_TRN_ECORRECT
    use mod_vadv1d, only: &
       vadv1d_getflux
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh, &
       GRD_afact, &
       GRD_bfact, &
       GRD_cfact, &
       GRD_dfact
    use mod_cnst, only: &
       CNST_GRAV
    use mod_thrmdyn, only: &
       THRMDYN_qd,        &
       THRMDYN_tempre
    use mod_cnvvar, only: &
       CNVVAR_rhogkin_in
    implicit none

    integer,  intent(in)    :: ijdim
    real(SP), intent(inout) :: rhog              (ijdim,kdim)
    real(SP), intent(inout) :: rhogvx            (ijdim,kdim)
    real(SP), intent(inout) :: rhogvy            (ijdim,kdim)
    real(SP), intent(inout) :: rhogvz            (ijdim,kdim)
    real(SP), intent(inout) :: rhogw             (ijdim,kdim)
    real(SP), intent(inout) :: rhoge             (ijdim,kdim)
    real(SP), intent(inout) :: rhogq             (ijdim,kdim,nqmax)
    real(SP), intent(inout) :: rho               (ijdim,kdim)
    real(SP), intent(inout) :: tem               (ijdim,kdim)
    real(SP), intent(inout) :: pre               (ijdim,kdim)
    real(SP), intent(in)    :: vx                (ijdim,kdim)
    real(SP), intent(in)    :: vy                (ijdim,kdim)
    real(SP), intent(in)    :: vz                (ijdim,kdim)
    real(SP), intent(in)    :: w                 (ijdim,kdim)
    real(SP), intent(inout) :: q                 (ijdim,kdim,nqmax)
    real(SP), intent(out)   :: qd                (ijdim,kdim)
    real(SP), intent(in)    :: WPREC             (ijdim,kdim,nqmax)
    logical,  intent(in)    :: precipitating_flag(nqmax)
    real(SP), intent(in)    :: z                 (ijdim,kdim) ! height
    real(SP), intent(in)    :: zh                (ijdim,kdim)
    real(SP), intent(out)   :: precip            (ijdim,2)    ! rain or snow
    real(SP), intent(out)   :: precip_rhoe       (ijdim)
    real(SP), intent(out)   :: precip_lh_heat    (ijdim)
    real(SP), intent(out)   :: precip_rhophi     (ijdim)
    real(SP), intent(out)   :: precip_rhokin     (ijdim)
    real(SP), intent(out)   :: precip_flux       (ijdim,kdim,2)
    real(SP), intent(in)    :: gsgam2            (ijdim,kdim)
    real(SP), intent(in)    :: gsgam2h           (ijdim,kdim)
    real(SP), intent(in)    :: rgs               (ijdim,kdim)
    real(SP), intent(in)    :: rgsh              (ijdim,kdim)
    real(SP), intent(in)    :: ix                (ijdim)
    real(SP), intent(in)    :: iy                (ijdim)
    real(SP), intent(in)    :: iz                (ijdim)
    real(SP), intent(in)    :: jx                (ijdim)
    real(SP), intent(in)    :: jy                (ijdim)
    real(SP), intent(in)    :: jz                (ijdim)
    real(DP), intent(in)    :: dt
    ![Add] 2012/02/01 T.Seiki
    real(SP), intent(inout), optional :: precip_trc(ijdim,nqmax)

    real(SP) :: fprec(ijdim,kdim,1:nqmax)
    real(SP) :: drho_prec(ijdim,kdim,1:nqmax)
    real(SP) :: rhoeq_prec(ijdim,kdim)
    real(SP) :: fprec_rhoe(ijdim,kdim,1:nqmax)
    real(SP) :: drhoe_prec(ijdim,kdim)
    real(SP) :: rhophiq_prec(ijdim,kdim)
    real(SP) :: rhokin_h_prec(ijdim,kdim)
    real(SP) :: rhokin_v_prec(ijdim,kdim)
    real(SP) :: fprec_rhophi(ijdim,kdim,1:nqmax)
    real(SP) :: fprec_rhokin_h(ijdim,kdim,1:nqmax)
    real(SP) :: fprec_rhokin_v(ijdim,kdim,1:nqmax)
    real(SP) :: drhophi_prec(ijdim,kdim)
    real(SP) :: drhokin_h_prec(ijdim,kdim)
    real(SP) :: drhokin_v_prec(ijdim,kdim)
    real(SP) :: rhouq(ijdim,kdim)
    real(SP) :: rhovq(ijdim,kdim)
    real(SP) :: rhowq(ijdim,kdim)
    real(SP) :: fprec_rhou(ijdim,kdim)
    real(SP) :: fprec_rhov(ijdim,kdim)
    real(SP) :: fprec_rhow(ijdim,kdim)
    real(SP) :: drhogu_prec(ijdim,kdim)
    real(SP) :: drhogv_prec(ijdim,kdim)
    real(SP) :: drhogw_prec(ijdim,kdim)
    real(SP) :: rhogkin(ijdim,kdim)
    real(SP) :: rhogkin_h(ijdim,kdim)
    real(SP) :: rhogkin_v(ijdim,kdim)
    real(SP) :: rhog_h(ijdim,kdim)
    real(SP) :: kin_h0(ijdim,kdim),kin_h(ijdim,kdim)
    real(SP) :: kin_v0(ijdim,kdim),kin_v(ijdim,kdim)
    real(SP) :: vx_t(ijdim,kdim),vy_t(ijdim,kdim),vz_t(ijdim,kdim),w_t(ijdim,kdim)
    real(SP) :: wq(ijdim,kdim), ein(ijdim,kdim)
    real(SP) :: tmp(ijdim,kdim)
    real(SP) :: tmp2(ijdim,kdim)
    real(SP) :: tmp_h         (ijdim,kdim)
    real(SP) :: tmp_v         (ijdim,kdim)
    real(SP) :: C2Wfact       (ijdim,kdim,2)
    real(SP) :: W2Cfact       (ijdim,kdim,2)
    real(SP) :: dz            (ijdim,kdim)
    real(SP) :: dzh           (ijdim,kdim)

    integer :: k, nq, ij
    !---------------------------------------------------------------------------

    call PROF_rapstart('____Precip_Transport',2)

!OCL XFILL
    do k  = kmin, kmax
    do ij = 1, ijdim
       dz (ij,k) = zh(ij,k+1) - zh(ij,k)
       dzh(ij,k) = z (ij,k)   - z (ij,k-1)
    enddo
    enddo
!OCL XFILL
    do ij = 1, ijdim
       dzh(ij,kmin) = z(ij,kmin) - zh(ij,kmin)
    enddo

    do k = kmin, kmax+1
       C2Wfact(:,k,1) = GRD_afact(k) / GSGAM2(:,k  ) * GSGAM2H(:,k)
       C2Wfact(:,k,2) = GRD_bfact(k) / GSGAM2(:,k-1) * GSGAM2H(:,k)
    enddo
    C2Wfact(:,kmin-1,1) = 0.0_RP
    C2Wfact(:,kmin-1,2) = 0.0_RP

    do k = kmin-1, kmax
       W2Cfact(:,k,1) = GRD_cfact(k) * GSGAM2(:,k) / GSGAM2H(:,k+1)
       W2Cfact(:,k,2) = GRD_dfact(k) * GSGAM2(:,k) / GSGAM2H(:,k  )
    enddo
    W2Cfact(:,kmax+1,1) = 0.0_RP
    W2Cfact(:,kmax+1,2) = 0.0_RP

    call CNVVAR_rhogkin_in( ijdim,            & ! [IN]
                            kdim,             & ! [IN]
                            rhog     (:,:),   & ! [IN]
                            rhogvx   (:,:),   & ! [IN]
                            rhogvy   (:,:),   & ! [IN]
                            rhogvz   (:,:),   & ! [IN]
                            rhogw    (:,:),   & ! [IN]
                            C2Wfact  (:,:,:), & ! [IN]
                            W2Cfact  (:,:,:), & ! [IN]
                            rhogkin  (:,:),   & ! [OUT]
                            rhogkin_h(:,:),   & ! [OUT]
                            rhogkin_v(:,:)    ) ! [OUT]

    !--- by rain
    ! [mod] 10/08/03 T.Mitsui
!!$ do nq = NQW_STR,NQW_END
    do nq = 1, nqmax
       if (.NOT.precipitating_flag(nq)) then
          fprec(:,:,nq) = 0.0D0
          fprec_rhoe(:,:,nq) = 0.0D0
          fprec_rhophi(:,:,nq) = 0.0D0
          fprec_rhokin_h(:,:,nq) = 0.0D0
          fprec_rhokin_v(:,:,nq) = 0.0D0
          cycle
       endif
       !----- mass
       fprec(:,:,nq) = 0.d0
       call vadv1d_getflux(             &
            ijdim,                  & !--- in
            kdim,                   &
            kmin,                   &
            kmax,                   &
            dz(:,:),                &
            zh(:,:),                &
            rhogq(:,:,nq)*rgs(:,:), & !--- in
            WPREC(:,:,nq),          & !--- in
            fprec(:,:,nq),          & !--- out
            dt )                      !--- in
       ! [Add] 10/08/03 T.Mitsui, hereafter nq should be active tracers
       if (nq >= NQW_STR .AND. nq <= NQW_END ) then
          !--- internal energy
          rhoeq_prec(:,:) = rhog(:,:) &
               * q(:,:,nq) * CVW(nq) * tem(:,:) * rgs(:,:)
          call vadv1d_getflux(        &
               ijdim,             & !--- in
               kdim,              &
               kmin,              &
               kmax,              &
               dz(:,:),           &
               zh(:,:),           &
               rhoeq_prec(:,:),   & !--- in
               WPREC(:,:,nq),     & !--- in
               fprec_rhoe(:,:,nq),& !--- out
               dt )
          !--- potential energy
          rhophiq_prec(:,:) = rhog(:,:) &
               * q(:,:,nq) * CNST_GRAV * z(:,:) * rgs(:,:)
          call vadv1d_getflux(           &
               ijdim,                & !--- in
               kdim,                 &
               kmin,                 &
               kmax,                 &
               dz(:,:),              &
               zh(:,:),              &
               rhophiq_prec(:,:),    & !--- in
               WPREC(:,:,nq),        & !--- in
               fprec_rhophi(:,:,nq), & !--- out
               dt )
          !--- horizontal kinetic energy
          rhokin_h_prec(:,:) = q(:,:,nq) * rhogkin_h(:,:) * rgs(:,:)
          call vadv1d_getflux(             &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               zh(:,:),                &
               rhokin_h_prec(:,:),     & !--- in
               WPREC(:,:,nq),          & !--- in
               fprec_rhokin_h(:,:,nq), & !--- out
               dt )
          !--- vertical kinetic energy
          rhokin_v_prec(:,:) = 0.0D0
          do k = kmin+1,kmax
             rhokin_v_prec(:,k) = 0.5_RP* (q(:,k,nq)+q(:,k-1,nq)) * rhogkin_v(:,k) * rgsh(:,k)
          enddo
          wq(:,:) = 0.0D0
          do k = kmin+1,kmax-1
             wq(:,k) = 0.5_RP*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
          enddo
          !<--- half -> full
          fprec_rhokin_v(:,:,nq) = 0.0D0
          call vadv1d_getflux(                      &
               ijdim,                           & !--- in
               kdim-1,                          &
               kmin,                            &
               kmax-1,                          &
               dzh(:,kmin:kmax+1),              &
               z(:,kmin-1:kmax),                &
               rhokin_v_prec(:,kmin:kmax+1),    & !--- in
               wq(:,kmin:kmax+1),               & !--- in
               fprec_rhokin_v(:,kmin-1:kmax,nq),& !--- out
               dt )
       endif
    enddo
    !--- tendency of density and energy
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq=NQW_STR,NQW_END
    do nq=1, nqmax
       drho_prec(:,kmin-1,nq) = 0.d0
       drho_prec(:,kmax+1,nq) = 0.d0
    enddo
    drhoe_prec(:,kmin-1) = 0.d0
    drhoe_prec(:,kmax+1) = 0.d0
    drhophi_prec(:,kmin-1) = 0.d0
    drhophi_prec(:,kmax+1) = 0.d0
    ! [Mod] 10/08/03 T.Mitsui
!!$    do k = kmin, kmax
!!$       do nq=NQW_STR,NQW_END
!!$          drho_prec(:,k,nq) &
!!$               = - ( fprec(:,k+1,nq) - fprec(:,k,nq) ) * GRD_rdgz(k)
!!$       enddo
!!$    enddo
    do nq=1, nqmax
       do k = kmin, kmax
          drho_prec(:,k,nq) &
               = - ( fprec(:,k+1,nq) - fprec(:,k,nq) ) * GRD_rdgz(k)
       enddo
    enddo
    do k= kmin, kmax
       drhoe_prec(:,k) = 0.0D0
       drhophi_prec(:,k) = 0.0D0
       drhokin_h_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhoe_prec(:,k) = drhoe_prec(:,k)&
               - ( fprec_rhoe(:,k+1,nq) - fprec_rhoe(:,k,nq) )*GRD_rdgz(k)
          drhophi_prec(:,k) = drhophi_prec(:,k)&
               - ( fprec_rhophi(:,k+1,nq) - fprec_rhophi(:,k,nq) )*GRD_rdgz(k)
          drhokin_h_prec(:,k) = drhokin_h_prec(:,k)&
               - ( fprec_rhokin_h(:,k+1,nq) - fprec_rhokin_h(:,k,nq) )*GRD_rdgz(k)
       enddo
    enddo
    drhokin_h_prec(:,kmax+1) = 0.0D0 ! Iga 090803
    drhokin_h_prec(:,kmin-1) = 0.0D0 ! Iga 090803
    drhokin_v_prec(:,kmin)=0.0D0
    drhokin_v_prec(:,kmin-1)=0.0D0
    drhokin_v_prec(:,kmax+1)=0.0D0
    do k = kmin+1, kmax
       drhokin_v_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhokin_v_prec(:,k) = drhokin_v_prec(:,k)&
               - ( fprec_rhokin_v(:,k,nq) - fprec_rhokin_v(:,k-1,nq) )*GRD_rdgzh(k)
       enddo
    enddo
    !--- momentum transport: new values due to rain
    !<------ use old values of rho and qr for consistency
    do nq=NQW_STR,NQW_END
       if (.NOT.precipitating_flag(nq)) then
          cycle
       endif
       do k = 1,kdim
          rhouq(:,k)                        &
               = ( rhog(:,k)*vx(:,k)*ix(:)  &
                 + rhog(:,k)*vy(:,k)*iy(:)  &
!                 + rhog(:,k)*vz(:,k)*iy(:) )& [bugfix] H.Yashiro 20120518
                 + rhog(:,k)*vz(:,k)*iz(:) )&
                 *q(:,k,nq)*rgs(:,k)
          rhovq(:,k)                 &
               = ( rhog(:,k)*vx(:,k)*jx(:)  &
                 + rhog(:,k)*vy(:,k)*jy(:)  &
!                 + rhog(:,k)*vz(:,k)*jy(:) )& [bugfix] H.Yashiro 20120518
                 + rhog(:,k)*vz(:,k)*jz(:) )&
                 *q(:,k,nq)*rgs(:,k)
       enddo
       call vadv1d_getflux(       &
            ijdim,            & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dz(:,:),          &
            zh(:,:),          &
            rhouq(:,:),       & !--- in
            WPREC(:,:,nq),    & !--- in
            fprec_rhou(:,:),  & !--- out
            dt                & !--- in
            )
       call vadv1d_getflux(  &
            ijdim,               & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dz(:,:),          &
            zh(:,:),          &
            rhovq(:,:),       & !--- in
            WPREC(:,:,nq),    & !--- in
            fprec_rhov(:,:),  & !--- out
            dt                   & !--- in
            )
       !--- update rhogvx, rhogvy, rhogvz
       drhogu_prec(:,kmin-1) = 0.d0
       drhogu_prec(:,kmax+1) = 0.d0
       drhogv_prec(:,kmin-1) = 0.d0
       drhogv_prec(:,kmax+1) = 0.d0
       do k = kmin, kmax
          drhogu_prec(:,k) = - ( fprec_rhou(:,k+1) - fprec_rhou(:,k) ) * GRD_rdgz(k)
          drhogv_prec(:,k) = - ( fprec_rhov(:,k+1) - fprec_rhov(:,k) ) * GRD_rdgz(k)
       enddo
       do k = kmin, kmax
          rhogvx(:,k) = rhogvx(:,k)    &
               + drhogu_prec(:,k)*ix(:)&
               + drhogv_prec(:,k)*jx(:)
          rhogvy(:,k) = rhogvy(:,k)    &
               + drhogu_prec(:,k)*iy(:)&
               + drhogv_prec(:,k)*jy(:)
          rhogvz(:,k) = rhogvz(:,k)    &
               + drhogu_prec(:,k)*iz(:)&
               + drhogv_prec(:,k)*jz(:)
       enddo
       !=====================
       rhowq(:,:) = 0.0D0
       do k = kmin+1,kmax
          rhowq(:,k) &
               = rhogw(:,k)* 0.5_RP* (q(:,k,nq)+q(:,k-1,nq)) *rgsh(:,k)
       enddo
       wq(:,:) = 0.0D0
       do k = kmin+1,kmax-1
          wq(:,k) = 0.5_RP*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
       enddo
       call vadv1d_getflux(            &
            ijdim,                     &
            kdim-1,                    &
            kmin,                      &
            kmax-1,                    &
            dzh(:,kmin:kmax+1),        &
            z(:,kmin-1:kmax),          &
            rhowq(:,kmin:kmax+1),      &
            wq(:,kmin:kmax+1),         &
            fprec_rhow(:,kmin-1:kmax), &
            dt )
       !--- update rhogw
       drhogw_prec(:,:) = 0.0D0
       do k = kmin+1, kmax
          drhogw_prec(:,k) = - ( fprec_rhow(:,k) - fprec_rhow(:,k-1) ) &
               * GRD_rdgzh(k)
       enddo
       rhogw(:,:) = rhogw(:,:) + drhogw_prec(:,:)
    enddo
    ! new values due to rain:
    ! Change in internal energy comes from precipitation of rain and dissipation
    ! of kinetic energy due to drag force.
    ! See Ooyama(2001) (3.13)
    !--- update rhogqr
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq=NQW_STR,NQW_END
    do nq=1, nqmax
       rhogq(:,:,nq) = rhogq(:,:,nq) + drho_prec(:,:,nq)
    enddo
    ! [comment] 10/08/03 T.Mitsui, number never change air mass here.
    !--- update rhog, rho
    do nq=NQW_STR,NQW_END
       rhog(:,:) = rhog(:,:) + drho_prec(:,:,nq)
    enddo
    rho(:,:) = rhog(:,:)/gsgam2(:,:)
    !--- update rhoge, rhogkin
    rhoge(:,:) = rhoge(:,:) + drhoe_prec(:,:) + drhophi_prec(:,:)
    do nq=NQW_STR,NQW_END
       rhoge(:,:) = rhoge(:,:) &
            - CNST_GRAV * z(:,:) * ( drho_prec(:,:,nq) )
    enddo

    if (adm_prc_me==1) then
       do ij=1,ijdim
!          write(*,*) 'a',ij,rhogkin_h(ij,kmin-1),rhogkin_h(ij,kmax+1)
!          write(*,*) 'c',ij,drhokin_h_prec(ij,kmin-1),drhokin_h_prec(ij,kmax+1)
       enddo
    endif
!    do k=1,kdim
!       do ij=1,ijdim

!          if (.NOT.((-1d20<drhokin_v_prec(ij,k)).AND.(drhokin_v_prec(ij,k)<1d20))) then
!             write(*,*) 'a',ij,k,drhokin_v_prec(ij,k),rhogkin_v(ij,k)
!          endif
!          if (.NOT.((-1d20<rhogkin_v(ij,k)).AND.(rhogkin_v(ij,k)<1d20))) then
!             write(*,*) 'b',ij,k,rhogkin_v(ij,k)
!          endif
!       enddo
!    enddo



    rhogkin_h(:,:) = rhogkin_h(:,:)+  drhokin_h_prec(:,:)
    rhogkin_v(:,:) = rhogkin_v(:,:)+  drhokin_v_prec(:,:)
!    rhogkin_h(:,kmin:kmax) = rhogkin_h(:,kmin:kmax)+  drhokin_h_prec(:,kmin:kmax)
!    rhogkin_v(:,kmin:kmax) = rhogkin_v(:,kmin:kmax)+  drhokin_v_prec(:,kmin:kmax)
!    rhogkin_v(:,:) = rhogkin_v(:,:)+1
    !--- update q
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq = NQW_STR,NQW_END
    do nq=1, nqmax
       q(:,:,nq) = rhogq(:,:,nq) / rhog(:,:)
    enddo
    !--- update qd
    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    if (trim(PRCIP_TRN_ECORRECT)=='KIN2EIN') then
!       tmp2(:,kmin-1) = 0.0D0
!       tmp2(:,kmax+1) = 0.0D0
       tmp2(:,kmin-1) = 1.0_RP  ! [mod] C.Kodama 2014.12.22
       tmp2(:,kmax+1) = 1.0_RP  !
       do k = kmin, kmax
          tmp2(:,k) = rhogkin_h(:,k)&
               + (GRD_cfact(k  )*rhogkin_v(:,k+1) &
               +  GRD_dfact(k  )*rhogkin_v(:,k  ) )
       enddo
       call CNVVAR_rhogkin_in( ijdim,          & ! [IN]
                               kdim,           & ! [IN]
                               rhog   (:,:),   & ! [IN]
                               rhogvx (:,:),   & ! [IN]
                               rhogvy (:,:),   & ! [IN]
                               rhogvz (:,:),   & ! [IN]
                               rhogw  (:,:),   & ! [IN]
                               C2Wfact(:,:,:), & ! [IN]
                               W2Cfact(:,:,:), & ! [IN]
                               tmp    (:,:),   & ! [OUT]
                               tmp_h  (:,:),   & ! [OUT]
                               tmp_v  (:,:)    ) ! [OUT]
       rhoge(:,:) = rhoge(:,:) + ( tmp2(:,:) - tmp(:,:) )
    elseif(trim(PRCIP_TRN_ECORRECT)=='KIN2KIN') then
       !--- fix kinetic energy
       ! Modify C.Kodama 09.07.14 =>
       kin_h0(:,kmin:kmax) = rhogkin_h(:,kmin:kmax)/rhog(:,kmin:kmax)
       vx_t(:,kmin:kmax) = rhogvx(:,kmin:kmax)/rhog(:,kmin:kmax)
       vy_t(:,kmin:kmax) = rhogvy(:,kmin:kmax)/rhog(:,kmin:kmax)
       vz_t(:,kmin:kmax) = rhogvz(:,kmin:kmax)/rhog(:,kmin:kmax)
       kin_h(:,kmin:kmax) = ((vx_t(:,kmin:kmax)**2)+(vy_t(:,kmin:kmax)**2)+(vz_t(:,kmin:kmax)**2))*0.5_RP
       where(kin_h(:,kmin:kmax)>1.0D-20)
          vx_t(:,kmin:kmax) = vx_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
          vy_t(:,kmin:kmax) = vy_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
          vz_t(:,kmin:kmax) = vz_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
       end where
       rhogvx(:,kmin:kmax) = rhog(:,kmin:kmax) * vx_t(:,kmin:kmax)
       rhogvy(:,kmin:kmax) = rhog(:,kmin:kmax) * vy_t(:,kmin:kmax)
       rhogvz(:,kmin:kmax) = rhog(:,kmin:kmax) * vz_t(:,kmin:kmax)
       rhog_h(:,kmin:kmax+1) = 1.0D0
       do k =kmin,kmax+1
          rhog_h(:,k) = &
               ( GRD_afact(k) * rhog(:,k  )/gsgam2(:,k  ) &
               + GRD_bfact(k) * rhog(:,k-1)/gsgam2(:,k-1) &
               ) * gsgam2h(:,k)
       enddo
       kin_v0(:,kmin:kmax+1) = rhogkin_v(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       w_t(:,kmin:kmax+1) = rhogw(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       kin_v(:,kmin:kmax+1) = (w_t(:,kmin:kmax+1)**2) * 0.5_RP
       where(kin_v(:,kmin:kmax+1)>1.0D-20)
          w_t(:,kmin:kmax+1) = w_t(:,kmin:kmax+1) *  sqrt(abs(kin_v0(:,kmin:kmax+1)/kin_v(:,kmin:kmax+1)))
       end where
       rhogw(:,kmin:kmax+1) = rhog_h(:,kmin:kmax+1) * w_t(:,kmin:kmax+1)
       !<= Modify C.Kodama 09.07.14
    else
       write(*,*) 'Error in precip_transport'
    endif
    !--- update temerature & pressure
    ein(:,:) =  rhoge(:,:) / rhog(:,:)
    call THRMDYN_tempre( ijdim, kdim, ein(:,:), rho(:,:), q(:,:,:), tem(:,:), pre(:,:) )

    precip     (:,:)   = 0.0D0
    precip_flux(:,:,:) = 0.0D0
    do nq=NQW_STR,NQW_END
       if ( (nq==I_QC).OR.(nq==I_QR) ) then
          precip(:,1) = precip(:,1) -fprec(:,kmin,nq)/dt
          precip_flux(:,:,1) = precip_flux(:,:,1) -fprec(:,:,nq)/dt
       elseif( (nq==I_QI).OR.(nq==I_QS).OR.(nq==I_QG) ) then
          precip(:,2) = precip(:,2) -fprec(:,kmin,nq)/dt
          precip_flux(:,:,2) = precip_flux(:,:,2) -fprec(:,:,nq)/dt
       endif
    enddo
    precip_rhoe(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhoe(:) = precip_rhoe(:) - fprec_rhoe(:,kmin,nq)/dt
    enddo
    precip_lh_heat(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       if ( (nq==I_QI).OR.(nq==I_QS).OR.(nq==I_QG) ) then
          precip_lh_heat(:) = precip_lh_heat(:) + fprec(:,kmin,nq)*LHF/dt
       endif
    enddo
    precip_rhophi(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhophi(:) = precip_rhophi(:) - fprec_rhophi(:,kmin,nq)/dt
    enddo
    precip_rhokin(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhokin(:) = precip_rhokin(:) &
            - fprec_rhokin_h(:,kmin,nq)/dt &
            - fprec_rhokin_v(:,kmin,nq)/dt
    enddo

    ! [Add] 2012/02/01 T.Seiki
    if ( present(precip_trc) ) then
       do nq=1, nqmax
          precip_trc(:,nq) = precip_trc(:,nq) -fprec(:,kmin,nq)/dt
       enddo
    endif

    call PROF_rapend  ('____Precip_Transport',2)

    return
  end subroutine precip_transport_SP

  !-----------------------------------------------------------------------------
  subroutine precip_transport_DP( &
       ijdim,              &
       rhog,               &
       rhogvx,             &
       rhogvy,             &
       rhogvz,             &
       rhogw,              &
       rhoge,              &
       rhogq,              &
       rho,                &
       tem,                &
       pre,                &
       vx,                 &
       vy,                 &
       vz,                 &
       w,                  &
       q,                  &
       qd,                 &
       z,                  &
       zh,                 &
       WPREC,              &
       precipitating_flag, &
       precip,             &
       precip_rhoe,        &
       precip_lh_heat,     &
       precip_rhophi,      &
       precip_rhokin,      &
       precip_flux,        &
       gsgam2,             &
       gsgam2h,            &
       rgs,                &
       rgsh,               &
       ix,                 &
       iy,                 &
       iz,                 &
       jx,                 &
       jy,                 &
       jz,                 &
       dt,                 &
       precip_trc          )
    use mod_adm, only: &
       ADM_PRC_ME,       &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR, NQW_END,  &
       I_QC,              &
       I_QR,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW,               &
       LHF,               &
       PRCIP_TRN_ECORRECT
    use mod_vadv1d, only: &
       vadv1d_getflux_0s_DP, &
       vadv1d_getflux_1s_DP, &
       vadv1d_getflux_2s_DP, &
       vadv1d_getflux
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh, &
       GRD_afact, &
       GRD_bfact, &
       GRD_cfact, &
       GRD_dfact
    use mod_cnst, only: &
       CNST_GRAV
    use mod_thrmdyn, only: &
       THRMDYN_qd,        &
       THRMDYN_tempre
    use mod_cnvvar, only: &
       CNVVAR_rhogkin_in
    implicit none

    integer,  intent(in)    :: ijdim
    real(DP), intent(inout) :: rhog              (ijdim,kdim)
    real(DP), intent(inout) :: rhogvx            (ijdim,kdim)
    real(DP), intent(inout) :: rhogvy            (ijdim,kdim)
    real(DP), intent(inout) :: rhogvz            (ijdim,kdim)
    real(DP), intent(inout) :: rhogw             (ijdim,kdim)
    real(DP), intent(inout) :: rhoge             (ijdim,kdim)
    real(DP), intent(inout) :: rhogq             (ijdim,kdim,nqmax)
    real(DP), intent(inout) :: rho               (ijdim,kdim)
    real(DP), intent(inout) :: tem               (ijdim,kdim)
    real(DP), intent(inout) :: pre               (ijdim,kdim)
    real(DP), intent(in)    :: vx                (ijdim,kdim)
    real(DP), intent(in)    :: vy                (ijdim,kdim)
    real(DP), intent(in)    :: vz                (ijdim,kdim)
    real(DP), intent(in)    :: w                 (ijdim,kdim)
    real(DP), intent(inout) :: q                 (ijdim,kdim,nqmax)
    real(DP), intent(out)   :: qd                (ijdim,kdim)
    real(DP), intent(in)    :: WPREC             (ijdim,kdim,nqmax)
    logical,  intent(in)    :: precipitating_flag(nqmax)
    real(DP), intent(in)    :: z                 (ijdim,kdim) ! height
    real(DP), intent(in)    :: zh                (ijdim,kdim)
    real(DP), intent(out)   :: precip            (ijdim,2)    ! rain or snow
    real(DP), intent(out)   :: precip_rhoe       (ijdim)
    real(DP), intent(out)   :: precip_lh_heat    (ijdim)
    real(DP), intent(out)   :: precip_rhophi     (ijdim)
    real(DP), intent(out)   :: precip_rhokin     (ijdim)
    real(DP), intent(out)   :: precip_flux       (ijdim,kdim,2)
    real(DP), intent(in)    :: gsgam2            (ijdim,kdim)
    real(DP), intent(in)    :: gsgam2h           (ijdim,kdim)
    real(DP), intent(in)    :: rgs               (ijdim,kdim)
    real(DP), intent(in)    :: rgsh              (ijdim,kdim)
    real(DP), intent(in)    :: ix                (ijdim)
    real(DP), intent(in)    :: iy                (ijdim)
    real(DP), intent(in)    :: iz                (ijdim)
    real(DP), intent(in)    :: jx                (ijdim)
    real(DP), intent(in)    :: jy                (ijdim)
    real(DP), intent(in)    :: jz                (ijdim)
    real(DP), intent(in)    :: dt
    ![Add] 2012/02/01 T.Seiki
    real(DP), intent(inout), optional :: precip_trc(ijdim,nqmax)

    real(DP) :: fprec(ijdim,kdim,1:nqmax)
    real(DP) :: drho_prec(ijdim,kdim,1:nqmax)
    real(DP) :: rhoeq_prec(ijdim,kdim)
    real(DP) :: fprec_rhoe(ijdim,kdim,1:nqmax)
    real(DP) :: drhoe_prec(ijdim,kdim)
    real(DP) :: rhophiq_prec(ijdim,kdim)
    real(DP) :: rhokin_h_prec(ijdim,kdim)
    real(DP) :: rhokin_v_prec(ijdim,kdim)
    real(DP) :: fprec_rhophi(ijdim,kdim,1:nqmax)
    real(DP) :: fprec_rhokin_h(ijdim,kdim,1:nqmax)
    real(DP) :: fprec_rhokin_v(ijdim,kdim,1:nqmax)
    real(DP) :: drhophi_prec(ijdim,kdim)
    real(DP) :: drhokin_h_prec(ijdim,kdim)
    real(DP) :: drhokin_v_prec(ijdim,kdim)
    real(DP) :: rhouq(ijdim,kdim)
    real(DP) :: rhovq(ijdim,kdim)
    real(DP) :: rhowq(ijdim,kdim)
    real(DP) :: fprec_rhou(ijdim,kdim)
    real(DP) :: fprec_rhov(ijdim,kdim)
    real(DP) :: fprec_rhow(ijdim,kdim)
    real(DP) :: drhogu_prec(ijdim,kdim)
    real(DP) :: drhogv_prec(ijdim,kdim)
    real(DP) :: drhogw_prec(ijdim,kdim)
    real(DP) :: rhogkin(ijdim,kdim)
    real(DP) :: rhogkin_h(ijdim,kdim)
    real(DP) :: rhogkin_v(ijdim,kdim)
    real(DP) :: rhog_h(ijdim,kdim)
    real(DP) :: kin_h0(ijdim,kdim),kin_h(ijdim,kdim)
    real(DP) :: kin_v0(ijdim,kdim),kin_v(ijdim,kdim)
    real(DP) :: vx_t(ijdim,kdim),vy_t(ijdim,kdim),vz_t(ijdim,kdim),w_t(ijdim,kdim)
    real(DP) :: wq(ijdim,kdim), ein(ijdim,kdim)
    real(DP) :: tmp(ijdim,kdim)
    real(DP) :: tmp2(ijdim,kdim)
    real(DP) :: tmp_h         (ijdim,kdim)
    real(DP) :: tmp_v         (ijdim,kdim)
    real(DP) :: C2Wfact       (ijdim,kdim,2)
    real(DP) :: W2Cfact       (ijdim,kdim,2)
    real(DP) :: dz            (ijdim,kdim)
    real(DP) :: dzh           (ijdim,kdim)

    integer  :: kcell     (ijdim,kdim,max(NQW_STR,1):MAX(NQW_END,1))
    real(DP) :: zdis      (ijdim,kdim,max(NQW_STR,1):MAX(NQW_END,1))
    integer  :: kcell_max (kdim,max(NQW_STR,1):MAX(NQW_END,1))
    integer  :: kcell_min (kdim,max(NQW_STR,1):MAX(NQW_END,1))

    integer :: k, nq, ij
    logical, save :: ofirst = .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('____Precip_Transport',2)

!OCL XFILL
    do k  = kmin, kmax
    do ij = 1, ijdim
       dz (ij,k) = zh(ij,k+1) - zh(ij,k)
       dzh(ij,k) = z (ij,k)   - z (ij,k-1)
    enddo
    enddo
!OCL XFILL
    do ij = 1, ijdim
       dzh(ij,kmin) = z(ij,kmin) - zh(ij,kmin)
    enddo

    do k = kmin, kmax+1
       C2Wfact(:,k,1) = GRD_afact(k) / GSGAM2(:,k  ) * GSGAM2H(:,k)
       C2Wfact(:,k,2) = GRD_bfact(k) / GSGAM2(:,k-1) * GSGAM2H(:,k)
    enddo
    C2Wfact(:,kmin-1,1) = 0.0_RP
    C2Wfact(:,kmin-1,2) = 0.0_RP

    do k = kmin-1, kmax
       W2Cfact(:,k,1) = GRD_cfact(k) * GSGAM2(:,k) / GSGAM2H(:,k+1)
       W2Cfact(:,k,2) = GRD_dfact(k) * GSGAM2(:,k) / GSGAM2H(:,k  )
    enddo
    W2Cfact(:,kmax+1,1) = 0.0_RP
    W2Cfact(:,kmax+1,2) = 0.0_RP

    if(ofirst)then
      call vadv1d_getflux_0s_DP( &
            ijdim,            & !--- IN : number of horizontal grids
            kdim,             & !--- IN
            kmin,             & !--- IN
            kmax,             & !--- IN
            dz,               & !--- IN
            zh,               & !--- IN
            dt                & !--- IN : time interval
            )
        ofirst = .false.
    endif

    call CNVVAR_rhogkin_in( ijdim,            & ! [IN]
                            kdim,             & ! [IN]
                            rhog     (:,:),   & ! [IN]
                            rhogvx   (:,:),   & ! [IN]
                            rhogvy   (:,:),   & ! [IN]
                            rhogvz   (:,:),   & ! [IN]
                            rhogw    (:,:),   & ! [IN]
                            C2Wfact  (:,:,:), & ! [IN]
                            W2Cfact  (:,:,:), & ! [IN]
                            rhogkin  (:,:),   & ! [OUT]
                            rhogkin_h(:,:),   & ! [OUT]
                            rhogkin_v(:,:)    ) ! [OUT]

    !--- by rain
    ! [mod] 10/08/03 T.Mitsui
!!$ do nq = NQW_STR,NQW_END
    do nq = 1, nqmax
       if (.NOT.precipitating_flag(nq)) then
          fprec(:,:,nq) = 0.0D0
          fprec_rhoe(:,:,nq) = 0.0D0
          fprec_rhophi(:,:,nq) = 0.0D0
          fprec_rhokin_h(:,:,nq) = 0.0D0
          fprec_rhokin_v(:,:,nq) = 0.0D0
          cycle
       endif
       !----- mass
       fprec(:,:,nq) = 0.d0
       if (nq >= NQW_STR .AND. nq <= NQW_END ) then
          call vadv1d_getflux_1s_DP( &
                ijdim,            & !--- IN : number of horizontal grids
                kdim,             & !--- IN
                kmin,             & !--- IN
                kmax,             & !--- IN
                dz,               & !--- IN
                zh,               & !--- IN
                WPREC(:,:,nq),    & !--- in
                dt                & !--- IN : time interval
                , kcell     (:,:,nq) & !--- OUT
                , zdis      (:,:,nq) & !--- OUT
                , kcell_max (:,nq)   & !--- OUT
                , kcell_min (:,nq)   & !--- OUT
                )

          call vadv1d_getflux_2s_DP(   &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               rhogq(:,:,nq)*rgs(:,:), & !--- in
               fprec(:,:,nq)           & !--- out
               , kcell     (:,:,nq)    & !--- in
               , zdis      (:,:,nq)    & !--- in
               , kcell_max (:,nq)      & !--- in
               , kcell_min (:,nq)      & !--- in
               )
          !--- internal energy
          rhoeq_prec(:,:) = rhog(:,:) &
               * q(:,:,nq) * CVW(nq) * tem(:,:) * rgs(:,:)
          call vadv1d_getflux_2s_DP(   &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               rhoeq_prec(:,:),        & !--- in
               fprec_rhoe(:,:,nq)      & !--- out
               , kcell     (:,:,nq)    & !--- in
               , zdis      (:,:,nq)    & !--- in
               , kcell_max (:,nq)      & !--- in
               , kcell_min (:,nq)      & !--- in
               )
          !--- potential energy
          rhophiq_prec(:,:) = rhog(:,:) &
               * q(:,:,nq) * CNST_GRAV * z(:,:) * rgs(:,:)
          call vadv1d_getflux_2s_DP(             &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               rhophiq_prec(:,:),      & !--- in
               fprec_rhophi(:,:,nq)    & !--- out
               , kcell     (:,:,nq)    & !--- in
               , zdis      (:,:,nq)    & !--- in
               , kcell_max (:,nq)      & !--- in
               , kcell_min (:,nq)      & !--- in
               )
          !--- horizontal kinetic energy
          rhokin_h_prec(:,:) = q(:,:,nq) * rhogkin_h(:,:) * rgs(:,:)
          call vadv1d_getflux_2s_DP(             &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               rhokin_h_prec(:,:),     & !--- in
               fprec_rhokin_h(:,:,nq)  & !--- out
               , kcell     (:,:,nq)    & !--- in
               , zdis      (:,:,nq)    & !--- in
               , kcell_max (:,nq)      & !--- in
               , kcell_min (:,nq)      & !--- in
               )
          !--- vertical kinetic energy
          rhokin_v_prec(:,:) = 0.0D0
          do k = kmin+1,kmax
             rhokin_v_prec(:,k) = 0.5_RP* (q(:,k,nq)+q(:,k-1,nq)) * rhogkin_v(:,k) * rgsh(:,k)
          enddo
          wq(:,:) = 0.0D0
          do k = kmin+1,kmax-1
             wq(:,k) = 0.5_RP*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
          enddo
          !<--- half -> full
          fprec_rhokin_v(:,:,nq) = 0.0D0
          call vadv1d_getflux(                      &
               ijdim,                           & !--- in
               kdim-1,                          &
               kmin,                            &
               kmax-1,                          &
               dzh(:,kmin:kmax+1),              &
               z(:,kmin-1:kmax),                &
               rhokin_v_prec(:,kmin:kmax+1),    & !--- in
               wq(:,kmin:kmax+1),               & !--- in
               fprec_rhokin_v(:,kmin-1:kmax,nq),& !--- out
               dt )
       else
          call vadv1d_getflux(             &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dz(:,:),                &
               zh(:,:),                &
               rhogq(:,:,nq)*rgs(:,:), & !--- in
               WPREC(:,:,nq),          & !--- in
               fprec(:,:,nq),          & !--- out
               dt )                      !--- in
       endif

    enddo
    !--- tendency of density and energy
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq=NQW_STR,NQW_END
    do nq=1, nqmax
       drho_prec(:,kmin-1,nq) = 0.d0
       drho_prec(:,kmax+1,nq) = 0.d0
    enddo
    drhoe_prec(:,kmin-1) = 0.d0
    drhoe_prec(:,kmax+1) = 0.d0
    drhophi_prec(:,kmin-1) = 0.d0
    drhophi_prec(:,kmax+1) = 0.d0
    ! [Mod] 10/08/03 T.Mitsui
!!$    do k = kmin, kmax
!!$       do nq=NQW_STR,NQW_END
!!$          drho_prec(:,k,nq) &
!!$               = - ( fprec(:,k+1,nq) - fprec(:,k,nq) ) * GRD_rdgz(k)
!!$       enddo
!!$    enddo
    do nq=1, nqmax
       do k = kmin, kmax
          drho_prec(:,k,nq) &
               = - ( fprec(:,k+1,nq) - fprec(:,k,nq) ) * GRD_rdgz(k)
       enddo
    enddo
    do k= kmin, kmax
       drhoe_prec(:,k) = 0.0D0
       drhophi_prec(:,k) = 0.0D0
       drhokin_h_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhoe_prec(:,k) = drhoe_prec(:,k)&
               - ( fprec_rhoe(:,k+1,nq) - fprec_rhoe(:,k,nq) )*GRD_rdgz(k)
          drhophi_prec(:,k) = drhophi_prec(:,k)&
               - ( fprec_rhophi(:,k+1,nq) - fprec_rhophi(:,k,nq) )*GRD_rdgz(k)
          drhokin_h_prec(:,k) = drhokin_h_prec(:,k)&
               - ( fprec_rhokin_h(:,k+1,nq) - fprec_rhokin_h(:,k,nq) )*GRD_rdgz(k)
       enddo
    enddo
    drhokin_h_prec(:,kmax+1) = 0.0D0 ! Iga 090803
    drhokin_h_prec(:,kmin-1) = 0.0D0 ! Iga 090803
    drhokin_v_prec(:,kmin)=0.0D0
    drhokin_v_prec(:,kmin-1)=0.0D0
    drhokin_v_prec(:,kmax+1)=0.0D0
    do k = kmin+1, kmax
       drhokin_v_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhokin_v_prec(:,k) = drhokin_v_prec(:,k)&
               - ( fprec_rhokin_v(:,k,nq) - fprec_rhokin_v(:,k-1,nq) )*GRD_rdgzh(k)
       enddo
    enddo
    !--- momentum transport: new values due to rain
    !<------ use old values of rho and qr for consistency
    do nq=NQW_STR,NQW_END
       if (.NOT.precipitating_flag(nq)) then
          cycle
       endif
       do k = 1,kdim
          rhouq(:,k)                        &
               = ( rhog(:,k)*vx(:,k)*ix(:)  &
                 + rhog(:,k)*vy(:,k)*iy(:)  &
!                 + rhog(:,k)*vz(:,k)*iy(:) )& [bugfix] H.Yashiro 20120518
                 + rhog(:,k)*vz(:,k)*iz(:) )&
                 *q(:,k,nq)*rgs(:,k)
          rhovq(:,k)                 &
               = ( rhog(:,k)*vx(:,k)*jx(:)  &
                 + rhog(:,k)*vy(:,k)*jy(:)  &
!                 + rhog(:,k)*vz(:,k)*jy(:) )& [bugfix] H.Yashiro 20120518
                 + rhog(:,k)*vz(:,k)*jz(:) )&
                 *q(:,k,nq)*rgs(:,k)
       enddo
       call vadv1d_getflux_2s_DP(             &
            ijdim,                  & !--- in
            kdim,                   &
            kmin,                   &
            kmax,                   &
            dz(:,:),                &
            rhouq(:,:),             & !--- in
            fprec_rhou(:,:)         & !--- out
            , kcell     (:,:,nq) & !--- in
            , zdis      (:,:,nq) & !--- in
            , kcell_max (:,nq)   & !--- in
            , kcell_min (:,nq)   & !--- in
            )
       call vadv1d_getflux_2s_DP(             &
            ijdim,                  & !--- in
            kdim,                   &
            kmin,                   &
            kmax,                   &
            dz(:,:),                &
            rhovq(:,:),             & !--- in
            fprec_rhov(:,:)         & !--- out
            , kcell     (:,:,nq) & !--- in
            , zdis      (:,:,nq) & !--- in
            , kcell_max (:,nq)   & !--- in
            , kcell_min (:,nq)   & !--- in
            )
       !--- update rhogvx, rhogvy, rhogvz
       drhogu_prec(:,kmin-1) = 0.d0
       drhogu_prec(:,kmax+1) = 0.d0
       drhogv_prec(:,kmin-1) = 0.d0
       drhogv_prec(:,kmax+1) = 0.d0
       do k = kmin, kmax
          drhogu_prec(:,k) = - ( fprec_rhou(:,k+1) - fprec_rhou(:,k) ) * GRD_rdgz(k)
          drhogv_prec(:,k) = - ( fprec_rhov(:,k+1) - fprec_rhov(:,k) ) * GRD_rdgz(k)
       enddo
       do k = kmin, kmax
          rhogvx(:,k) = rhogvx(:,k)    &
               + drhogu_prec(:,k)*ix(:)&
               + drhogv_prec(:,k)*jx(:)
          rhogvy(:,k) = rhogvy(:,k)    &
               + drhogu_prec(:,k)*iy(:)&
               + drhogv_prec(:,k)*jy(:)
          rhogvz(:,k) = rhogvz(:,k)    &
               + drhogu_prec(:,k)*iz(:)&
               + drhogv_prec(:,k)*jz(:)
       enddo
       !=====================
       rhowq(:,:) = 0.0D0
       do k = kmin+1,kmax
          rhowq(:,k) &
               = rhogw(:,k)* 0.5_RP* (q(:,k,nq)+q(:,k-1,nq)) *rgsh(:,k)
       enddo
       wq(:,:) = 0.0D0
       do k = kmin+1,kmax-1
          wq(:,k) = 0.5_RP*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
       enddo
       call vadv1d_getflux(            &
            ijdim,                     &
            kdim-1,                    &
            kmin,                      &
            kmax-1,                    &
            dzh(:,kmin:kmax+1),        &
            z(:,kmin-1:kmax),          &
            rhowq(:,kmin:kmax+1),      &
            wq(:,kmin:kmax+1),         &
            fprec_rhow(:,kmin-1:kmax), &
            dt )
       !--- update rhogw
       drhogw_prec(:,:) = 0.0D0
       do k = kmin+1, kmax
          drhogw_prec(:,k) = - ( fprec_rhow(:,k) - fprec_rhow(:,k-1) ) &
               * GRD_rdgzh(k)
       enddo
       rhogw(:,:) = rhogw(:,:) + drhogw_prec(:,:)
    enddo
    ! new values due to rain:
    ! Change in internal energy comes from precipitation of rain and dissipation
    ! of kinetic energy due to drag force.
    ! See Ooyama(2001) (3.13)
    !--- update rhogqr
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq=NQW_STR,NQW_END
    do nq=1, nqmax
       rhogq(:,:,nq) = rhogq(:,:,nq) + drho_prec(:,:,nq)
    enddo
    ! [comment] 10/08/03 T.Mitsui, number never change air mass here.
    !--- update rhog, rho
    do nq=NQW_STR,NQW_END
       rhog(:,:) = rhog(:,:) + drho_prec(:,:,nq)
    enddo
    rho(:,:) = rhog(:,:)/gsgam2(:,:)
    !--- update rhoge, rhogkin
    rhoge(:,:) = rhoge(:,:) + drhoe_prec(:,:) + drhophi_prec(:,:)
    do nq=NQW_STR,NQW_END
       rhoge(:,:) = rhoge(:,:) &
            - CNST_GRAV * z(:,:) * ( drho_prec(:,:,nq) )
    enddo

    if (adm_prc_me==1) then
       do ij=1,ijdim
!          write(*,*) 'a',ij,rhogkin_h(ij,kmin-1),rhogkin_h(ij,kmax+1)
!          write(*,*) 'c',ij,drhokin_h_prec(ij,kmin-1),drhokin_h_prec(ij,kmax+1)
       enddo
    endif
!    do k=1,kdim
!       do ij=1,ijdim

!          if (.NOT.((-1d20<drhokin_v_prec(ij,k)).AND.(drhokin_v_prec(ij,k)<1d20))) then
!             write(*,*) 'a',ij,k,drhokin_v_prec(ij,k),rhogkin_v(ij,k)
!          endif
!          if (.NOT.((-1d20<rhogkin_v(ij,k)).AND.(rhogkin_v(ij,k)<1d20))) then
!             write(*,*) 'b',ij,k,rhogkin_v(ij,k)
!          endif
!       enddo
!    enddo



    rhogkin_h(:,:) = rhogkin_h(:,:)+  drhokin_h_prec(:,:)
    rhogkin_v(:,:) = rhogkin_v(:,:)+  drhokin_v_prec(:,:)
!    rhogkin_h(:,kmin:kmax) = rhogkin_h(:,kmin:kmax)+  drhokin_h_prec(:,kmin:kmax)
!    rhogkin_v(:,kmin:kmax) = rhogkin_v(:,kmin:kmax)+  drhokin_v_prec(:,kmin:kmax)
!    rhogkin_v(:,:) = rhogkin_v(:,:)+1
    !--- update q
    ! [Mod] 10/08/03 T.Mitsui
!!$ do nq = NQW_STR,NQW_END
    do nq=1, nqmax
       q(:,:,nq) = rhogq(:,:,nq) / rhog(:,:)
    enddo
    !--- update qd
    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    if (trim(PRCIP_TRN_ECORRECT)=='KIN2EIN') then
!       tmp2(:,kmin-1) = 0.0D0
!       tmp2(:,kmax+1) = 0.0D0
       tmp2(:,kmin-1) = 1.0_RP  ! [mod] C.Kodama 2014.12.22
       tmp2(:,kmax+1) = 1.0_RP  !
       do k = kmin, kmax
          tmp2(:,k) = rhogkin_h(:,k)&
               + (GRD_cfact(k  )*rhogkin_v(:,k+1) &
               +  GRD_dfact(k  )*rhogkin_v(:,k  ) )
       enddo
       call CNVVAR_rhogkin_in( ijdim,          & ! [IN]
                               kdim,           & ! [IN]
                               rhog   (:,:),   & ! [IN]
                               rhogvx (:,:),   & ! [IN]
                               rhogvy (:,:),   & ! [IN]
                               rhogvz (:,:),   & ! [IN]
                               rhogw  (:,:),   & ! [IN]
                               C2Wfact(:,:,:), & ! [IN]
                               W2Cfact(:,:,:), & ! [IN]
                               tmp    (:,:),   & ! [OUT]
                               tmp_h  (:,:),   & ! [OUT]
                               tmp_v  (:,:)    ) ! [OUT]
       rhoge(:,:) = rhoge(:,:) + ( tmp2(:,:) - tmp(:,:) )
    elseif(trim(PRCIP_TRN_ECORRECT)=='KIN2KIN') then
       !--- fix kinetic energy
       ! Modify C.Kodama 09.07.14 =>
       kin_h0(:,kmin:kmax) = rhogkin_h(:,kmin:kmax)/rhog(:,kmin:kmax)
       vx_t(:,kmin:kmax) = rhogvx(:,kmin:kmax)/rhog(:,kmin:kmax)
       vy_t(:,kmin:kmax) = rhogvy(:,kmin:kmax)/rhog(:,kmin:kmax)
       vz_t(:,kmin:kmax) = rhogvz(:,kmin:kmax)/rhog(:,kmin:kmax)
       kin_h(:,kmin:kmax) = ((vx_t(:,kmin:kmax)**2)+(vy_t(:,kmin:kmax)**2)+(vz_t(:,kmin:kmax)**2))*0.5_RP
       where(kin_h(:,kmin:kmax)>1.0D-20)
          vx_t(:,kmin:kmax) = vx_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
          vy_t(:,kmin:kmax) = vy_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
          vz_t(:,kmin:kmax) = vz_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
       end where
       rhogvx(:,kmin:kmax) = rhog(:,kmin:kmax) * vx_t(:,kmin:kmax)
       rhogvy(:,kmin:kmax) = rhog(:,kmin:kmax) * vy_t(:,kmin:kmax)
       rhogvz(:,kmin:kmax) = rhog(:,kmin:kmax) * vz_t(:,kmin:kmax)
       rhog_h(:,kmin:kmax+1) = 1.0D0
       do k =kmin,kmax+1
          rhog_h(:,k) = &
               ( GRD_afact(k) * rhog(:,k  )/gsgam2(:,k  ) &
               + GRD_bfact(k) * rhog(:,k-1)/gsgam2(:,k-1) &
               ) * gsgam2h(:,k)
       enddo
       kin_v0(:,kmin:kmax+1) = rhogkin_v(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       w_t(:,kmin:kmax+1) = rhogw(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       kin_v(:,kmin:kmax+1) = (w_t(:,kmin:kmax+1)**2) * 0.5_RP
       where(kin_v(:,kmin:kmax+1)>1.0D-20)
          w_t(:,kmin:kmax+1) = w_t(:,kmin:kmax+1) *  sqrt(abs(kin_v0(:,kmin:kmax+1)/kin_v(:,kmin:kmax+1)))
       end where
       rhogw(:,kmin:kmax+1) = rhog_h(:,kmin:kmax+1) * w_t(:,kmin:kmax+1)
       !<= Modify C.Kodama 09.07.14
    else
       write(*,*) 'Error in precip_transport'
    endif
    !--- update temerature & pressure
    ein(:,:) =  rhoge(:,:) / rhog(:,:)
    call THRMDYN_tempre( ijdim, kdim, ein(:,:), rho(:,:), q(:,:,:), tem(:,:), pre(:,:) )

    precip     (:,:)   = 0.0D0
    precip_flux(:,:,:) = 0.0D0
    do nq=NQW_STR,NQW_END
       if ( (nq==I_QC).OR.(nq==I_QR) ) then
          precip(:,1) = precip(:,1) -fprec(:,kmin,nq)/dt
          precip_flux(:,:,1) = precip_flux(:,:,1) -fprec(:,:,nq)/dt
       elseif( (nq==I_QI).OR.(nq==I_QS).OR.(nq==I_QG) ) then
          precip(:,2) = precip(:,2) -fprec(:,kmin,nq)/dt
          precip_flux(:,:,2) = precip_flux(:,:,2) -fprec(:,:,nq)/dt
       endif
    enddo
    precip_rhoe(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhoe(:) = precip_rhoe(:) - fprec_rhoe(:,kmin,nq)/dt
    enddo
    precip_lh_heat(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       if ( (nq==I_QI).OR.(nq==I_QS).OR.(nq==I_QG) ) then
          precip_lh_heat(:) = precip_lh_heat(:) + fprec(:,kmin,nq)*LHF/dt
       endif
    enddo
    precip_rhophi(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhophi(:) = precip_rhophi(:) - fprec_rhophi(:,kmin,nq)/dt
    enddo
    precip_rhokin(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhokin(:) = precip_rhokin(:) &
            - fprec_rhokin_h(:,kmin,nq)/dt &
            - fprec_rhokin_v(:,kmin,nq)/dt
    enddo

    ! [Add] 2012/02/01 T.Seiki
    if ( present(precip_trc) ) then
       do nq=1, nqmax
          precip_trc(:,nq) = precip_trc(:,nq) -fprec(:,kmin,nq)/dt
       enddo
    endif

    call PROF_rapend  ('____Precip_Transport',2)

    return
  end subroutine precip_transport_DP

  !-----------------------------------------------------------------------------
  subroutine precip_transport_new( &
       ijdim,              &
       rhog,               &
       rhogvx,             &
       rhogvy,             &
       rhogvz,             &
       rhogw,              &
       rhoge,              &
       rhogq,              &
       rho,                &
       tem,                &
       pre,                &
       vx,                 &
       vy,                 &
       vz,                 &
       w,                  &
       q,                  &
       qd,                 &
       z,                  &
       zh,                 &
       WPREC,              &
       precipitating_flag, &
       precip,             &
       precip_ISO1,        & ! [Add] M.Tanoue 20200908
       precip_ISO2,        & ! [Add] M.Tanoue 20200908
       precip_rhoe,        &
       precip_lh_heat,     &
       precip_rhophi,      &
       precip_rhokin,      &
       precip_flux,        &
       gsgam2,             &
       gsgam2h,            &
       rgs,                &
       rgsh,               &
       ix,                 &
       iy,                 &
       iz,                 &
       jx,                 &
       jy,                 &
       jz,                 &
       dt,                 &
       precip_trc          )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_EPS,  &
       CNST_GRAV
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh, &
       GRD_afact, &
       GRD_bfact, &
       GRD_cfact, &
       GRD_dfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QC,              &
       I_QR,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       opt_water_isotope, & ! [Add] M.Tanoue 20200908
       ISO1_QC,           & ! [Add] M.Tanoue 20200908
       ISO2_QC,           & ! [Add] M.Tanoue 20200908
       ISO1_QR,           & ! [Add] M.Tanoue 20200908
       ISO2_QR,           & ! [Add] M.Tanoue 20200908
       ISO1_QS,           & ! [Add] M.Tanoue 20200908
       ISO2_QS,           & ! [Add] M.Tanoue 20200908
       ISO1_QI,           & ! [Add] M.Tanoue 20200908
       ISO2_QI,           & ! [Add] M.Tanoue 20200908
       ISO1_QG,           & ! [Add] M.Tanoue 20200908
       ISO2_QG,           & ! [Add] M.Tanoue 20200908
       CVW,               &
       LHF,               &
       PRCIP_TRN_ECORRECT
    use mod_cnvvar, only: &
       CNVVAR_rhogkin_in, &
       cnvvar_wait_on => acc_wait_on
    use mod_thrmdyn, only: &
       THRMDYN_qd,    &
       THRMDYN_tempre,&
       thrmdyn_wait_on => acc_wait_on
    use mod_vadv1d, only: &
       vadv1d_setup,      &
       vadv1d_prep,       &
       vadv1d_getflux_new
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(inout) :: rhog              (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz            (ijdim,kdim)
    real(RP), intent(inout) :: rhogw             (ijdim,kdim)
    real(RP), intent(inout) :: rhoge             (ijdim,kdim)
    real(RP), intent(inout) :: rhogq             (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: rho               (ijdim,kdim)
    real(RP), intent(inout) :: tem               (ijdim,kdim)
    real(RP), intent(inout) :: pre               (ijdim,kdim)
    real(RP), intent(in)    :: vx                (ijdim,kdim)
    real(RP), intent(in)    :: vy                (ijdim,kdim)
    real(RP), intent(in)    :: vz                (ijdim,kdim)
    real(RP), intent(in)    :: w                 (ijdim,kdim)
    real(RP), intent(inout) :: q                 (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd                (ijdim,kdim)
    real(RP), intent(in)    :: z                 (ijdim,kdim)
    real(RP), intent(in)    :: zh                (ijdim,kdim)
    real(RP), intent(in)    :: WPREC             (ijdim,kdim,nqmax)
    logical,  intent(in)    :: precipitating_flag(nqmax)
    real(RP), intent(out)   :: precip            (ijdim,2)
    real(RP), intent(out)   :: precip_ISO1       (ijdim,2) ! [Add] M.Tanoue 20200908
    real(RP), intent(out)   :: precip_ISO2       (ijdim,2) ! [Add] M.Tanoue 20200908
    real(RP), intent(out)   :: precip_rhoe       (ijdim)
    real(RP), intent(out)   :: precip_lh_heat    (ijdim)
    real(RP), intent(out)   :: precip_rhophi     (ijdim)
    real(RP), intent(out)   :: precip_rhokin     (ijdim)
    real(RP), intent(out)   :: precip_flux       (ijdim,kdim,2)
    real(RP), intent(in)    :: gsgam2            (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h           (ijdim,kdim)
    real(RP), intent(in)    :: rgs               (ijdim,kdim)
    real(RP), intent(in)    :: rgsh              (ijdim,kdim)
    real(RP), intent(in)    :: ix                (ijdim)
    real(RP), intent(in)    :: iy                (ijdim)
    real(RP), intent(in)    :: iz                (ijdim)
    real(RP), intent(in)    :: jx                (ijdim)
    real(RP), intent(in)    :: jy                (ijdim)
    real(RP), intent(in)    :: jz                (ijdim)
    real(DP), intent(in)    :: dt
    real(RP), intent(inout), optional :: precip_trc(ijdim,nqmax) ! [Add] 2012/02/01 T.Seiki

    real(RP) :: rhogkin       (ijdim,kdim)
    real(RP) :: rhogkin_h     (ijdim,kdim)
    real(RP) :: rhogkin_v     (ijdim,kdim)

    real(RP) :: rhoq          (ijdim,kdim)
    real(RP) :: rhoqe         (ijdim,kdim)
    real(RP) :: rhoqphi       (ijdim,kdim)
    real(RP) :: rhoqkin_h     (ijdim,kdim)
    real(RP) :: rhoqkin_v     (ijdim,kdim)
    real(RP) :: rhoqu         (ijdim,kdim)
    real(RP) :: rhoqv         (ijdim,kdim)
    real(RP) :: rhoqw         (ijdim,kdim)

    real(RP) :: flux_rhoq     (ijdim,kdim)
    real(RP) :: flux_rhoqe    (ijdim,kdim)
    real(RP) :: flux_rhoqphi  (ijdim,kdim)
    real(RP) :: flux_rhoqkin_h(ijdim,kdim)
    real(RP) :: flux_rhoqkin_v(ijdim,kdim)
    real(RP) :: flux_rhoqu    (ijdim,kdim)
    real(RP) :: flux_rhoqv    (ijdim,kdim)
    real(RP) :: flux_rhoqw    (ijdim,kdim)

    real(RP) :: drhogq        (ijdim,kdim,nqmax)
    real(RP) :: drhoge        (ijdim,kdim)
    real(RP) :: drhogphi      (ijdim,kdim)
    real(RP) :: drhogkin_h    (ijdim,kdim)
    real(RP) :: drhogkin_v    (ijdim,kdim)
    real(RP) :: drhogu        (ijdim,kdim)
    real(RP) :: drhogv        (ijdim,kdim)
    real(RP) :: drhogw        (ijdim,kdim)

    real(RP) :: kin_h0        (ijdim,kdim)
    real(RP) :: kin_h         (ijdim,kdim)
    real(RP) :: vx_t          (ijdim,kdim)
    real(RP) :: vy_t          (ijdim,kdim)
    real(RP) :: vz_t          (ijdim,kdim)
    real(RP) :: kin_v0        (ijdim,kdim)
    real(RP) :: kin_v         (ijdim,kdim)
    real(RP) :: w_t           (ijdim,kdim)

    real(RP) :: zdis0         (ijdim,kdim)
    integer  :: kcell         (ijdim,kdim)
    integer  :: kcell_lim     (kdim)
    integer  :: kcell_max     (kdim)
    real(RP) :: zdis0h        (ijdim,kdim)
    integer  :: kcellh        (ijdim,kdim)
    integer  :: kcellh_lim    (kdim)
    integer  :: kcellh_max    (kdim)

    real(RP) :: WPRECH        (ijdim,kdim)
    real(RP) :: qh
    real(RP) :: rhogh         (ijdim,kdim)
    real(RP) :: ein           (ijdim,kdim)
    real(RP) :: tmp           (ijdim,kdim)
    real(RP) :: tmp2          (ijdim,kdim)
    real(RP) :: tmp_h         (ijdim,kdim)
    real(RP) :: tmp_v         (ijdim,kdim)
    real(RP) :: C2Wfact       (ijdim,kdim,2)
    real(RP) :: W2Cfact       (ijdim,kdim,2)

    real(RP) :: dz            (ijdim,kdim)
    real(RP) :: dzh           (ijdim,kdim)
    real(RP) :: dz2           (ijdim,kdim)
    real(RP) :: dzh2          (ijdim,kdim)
    real(RP) :: z_shift       (ijdim,kdim)

    real(RP) :: GRAV, EPS
    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data &
    !$acc pcreate(rhogkin,rhogkin_h,rhogkin_v) &
    !$acc pcreate(rhoq,rhoqe,rhoqphi,rhoqkin_h,rhoqkin_v,rhoqu,rhoqv,rhoqw) &
    !$acc pcreate(flux_rhoq,flux_rhoqe,flux_rhoqphi,flux_rhoqkin_h,flux_rhoqkin_v,flux_rhoqu,flux_rhoqv,flux_rhoqw) &
    !$acc pcreate(drhogq,drhoge,drhogphi,drhogkin_h,drhogkin_v,drhogu,drhogv,drhogw) &
    !$acc pcreate(kin_h0,kin_h,vx_t,vy_t,vz_t,kin_v0,kin_v,w_t) &
    !$acc pcreate(zdis0,kcell,kcell_max,zdis0h,kcellh,kcellh_max) &
    !$acc pcreate(WPRECH,qh,rhogh,ein,tmp,tmp2,tmp_h,tmp_v,C2Wfact,W2Cfact,dz,dzh,dz2,dzh2,z_shift) &
    !$acc pcopyin(vx,vy,vz,w,z,zh,WPREC,precipitating_flag) &
    !$acc pcopyin(gsgam2,gsgam2h,rgs,rgsh,ix,iy,iz,jx,jy,jz) &
    !$acc pcopyin(GRD_rdgz,GRD_rdgzh,GRD_afact,GRD_bfact,GRD_cfact,GRD_dfact,CVW) &
    !$acc pcreate(kcell_lim,kcellh_lim)

    call PROF_rapstart('____Precip_Transport',2)

    cnvvar_wait_on = .false.
    thrmdyn_wait_on = .false.

    EPS  = CNST_EPS
    GRAV = CNST_GRAV

#ifdef MORETIMER
call PROF_rapstart('____Precip_Transport_01')
#endif

!OCL XFILL
    !$acc kernels pcopy(dz,dzh) pcopyin(z,zh) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,dz,dzh,zh,z)
    do k  = kmin, kmax
    do ij = 1, ijdim
       dz (ij,k) = zh(ij,k+1) - zh(ij,k)
       dzh(ij,k) = z (ij,k)   - z (ij,k-1)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(dzh) pcopyin(zh) async(1)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,dzh,zh,z)
    do ij = 1, ijdim
       dzh(ij,kmin) = z(ij,kmin) - zh(ij,kmin)
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(dz2,dzh2) pcopyin(dz,dzh) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,dz2,dzh2,dz,dzh)
    do k  = kmin+1, kmax
    do ij = 1, ijdim
       dz2 (ij,k) = dz (ij,k) + dz (ij,k-1)
       dzh2(ij,k) = dzh(ij,k) + dzh(ij,k-1)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(z_shift) pcopyin(z) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,z_shift,z)
    do k  = kmin, kmax+1
    do ij = 1, ijdim
       z_shift(ij,k) = z(ij,k-1)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(C2Wfact) pcopyin(GRD_afact,GRD_bfact,gsgam2,gsgam2h) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,C2Wfact,GRD_afact,GRD_bfact,gsgam2,gsgam2h)
    do k  = kmin, kmax+1
    do ij = 1, ijdim
       C2Wfact(ij,k,1) = GRD_afact(k) / gsgam2(ij,k  ) * gsgam2h(ij,k)
       C2Wfact(ij,k,2) = GRD_bfact(k) / gsgam2(ij,k-1) * gsgam2h(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(C2Wfact) async(1)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,C2Wfact)
    do ij = 1, ijdim
       C2Wfact(ij,kmin-1,1) = 0.0_RP
       C2Wfact(ij,kmin-1,2) = 0.0_RP
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(W2Cfact) pcopyin(GRD_cfact,GRD_dfact,gsgam2,gsgam2h) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,W2Cfact,GRD_cfact,GRD_dfact,gsgam2,gsgam2h)
    do k  = kmin-1, kmax
    do ij = 1, ijdim
       W2Cfact(ij,k,1) = GRD_cfact(k) * gsgam2(ij,k) / gsgam2h(ij,k+1)
       W2Cfact(ij,k,2) = GRD_dfact(k) * gsgam2(ij,k) / gsgam2h(ij,k  )
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(W2Cfact) async(1)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmax,W2Cfact)
    do ij = 1, ijdim
       W2Cfact(ij,kmax+1,1) = 0.0_RP
       W2Cfact(ij,kmax+1,2) = 0.0_RP
    enddo
    !$omp end parallel do
    !$acc end kernels

    call vadv1d_setup( ijdim, kdim, kmin, kmax, & ! [IN]
                       dz       (:,:),          & ! [IN]
                       zh       (:,:),          & ! [IN]
                       kcell_lim(:),            & ! [OUT]
                       dt                       ) ! [IN]

    call vadv1d_setup( ijdim, kdim, kmin+1, kmax, & ! [IN]
                       dzh       (:,:),           & ! [IN]
                       z_shift   (:,:),           & ! [IN]
                       kcellh_lim(:),             & ! [OUT]
                       dt                         ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_01')
call PROF_rapstart('____Precip_Transport_02z')
#endif

    call CNVVAR_rhogkin_in( ijdim,            & ! [IN]
                            kdim,             & ! [IN]
                            rhog     (:,:),   & ! [IN]
                            rhogvx   (:,:),   & ! [IN]
                            rhogvy   (:,:),   & ! [IN]
                            rhogvz   (:,:),   & ! [IN]
                            rhogw    (:,:),   & ! [IN]
                            C2Wfact  (:,:,:), & ! [IN]
                            W2Cfact  (:,:,:), & ! [IN]
                            rhogkin  (:,:),   & ! [OUT]
                            rhogkin_h(:,:),   & ! [OUT]
                            rhogkin_v(:,:)    ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_02z')
call PROF_rapstart('____Precip_Transport_03')
#endif

!OCL XFILL
    !$acc kernels async(1) pcopy(drhogq)
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,nqmax,drhogq)
!OCL NOPREFETCH
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       drhogq(ij,k,nq) = 0.0_RP
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(drhoge,drhogphi,drhogkin_h,drhogkin_v,drhogu,drhogv,drhogw,precip_flux)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,drhoge,drhogphi,drhogkin_h,drhogkin_v,drhogu,drhogv,drhogw,precip_flux)
!OCL NOPREFETCH
    do k  = 1, kdim
    do ij = 1, ijdim
       drhoge    (ij,k) = 0.0_RP
       drhogphi  (ij,k) = 0.0_RP
       drhogkin_h(ij,k) = 0.0_RP
       drhogkin_v(ij,k) = 0.0_RP
       drhogu    (ij,k) = 0.0_RP
       drhogv    (ij,k) = 0.0_RP
       drhogw    (ij,k) = 0.0_RP

       precip_flux(ij,k,1) = 0.0_RP
       precip_flux(ij,k,2) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(precip,precip_rhoe,precip_lh_heat,precip_rhophi,precip_rhokin)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,precip,precip_rhoe,precip_lh_heat,precip_rhophi,precip_rhokin)
!OCL NOPREFETCH
    do ij = 1, ijdim
       precip        (ij,1) = 0.0_RP ! rain
       precip        (ij,2) = 0.0_RP ! snow
       precip_rhoe   (ij)   = 0.0_RP
       precip_lh_heat(ij)   = 0.0_RP
       precip_rhophi (ij)   = 0.0_RP
       precip_rhokin (ij)   = 0.0_RP
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(precip_ISO1,precip_ISO2)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,precip_ISO1,precip_ISO2)
    do ij = 1, ijdim
       precip_ISO1(ij,1) = 0.0_RP ! rain (for H218O)
       precip_ISO2(ij,1) = 0.0_RP ! rain (for HDO)
       precip_ISO1(ij,2) = 0.0_RP ! snow (for H218O)
       precip_ISO2(ij,2) = 0.0_RP ! snow (for HDO)
    enddo
    !$omp end parallel do
    !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_03')
#endif

    do nq = 1, nqmax

       if ( .NOT. precipitating_flag(nq) ) cycle

#ifdef MORETIMER
call PROF_rapstart('____Precip_Transport_04')
#endif

!OCL XFILL
       !$acc kernels async(1) pcopy(WPRECH) pcopyin(GRD_afact,GRD_bfact,WPREC)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kmin,kmax,WPRECH,WPREC,GRD_afact,GRD_bfact)
       do k  = kmin+1, kmax
       do ij = 1, ijdim
          WPRECH(ij,k) = ( GRD_afact(k) * WPREC(ij,k  ,nq) &
                         + GRD_bfact(k) * WPREC(ij,k-1,nq) )
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

!OCL XFILL
       !$acc kernels pcopy(WPRECH) pcopyin(WPREC) async(1)
       !$omp parallel do default(none),private(ij), &
       !$omp shared(nq,ijdim,kmin,kmax,WPRECH,WPREC)
       do ij = 1, ijdim
          WPRECH(ij,kmin-1) = 0.0_RP
          WPRECH(ij,kmin  ) = WPREC(ij,kmin,nq)
          WPRECH(ij,kmax+1) = WPREC(ij,kmax,nq)
       enddo
       !$omp end parallel do
       !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_04')
call PROF_rapstart('____Precip_Transport_05z')
#endif

       call vadv1d_prep( ijdim, kdim, kmin, kmax, & ! [IN]
                         dz       (:,:),          & ! [IN]
                         dz2      (:,:),          & ! [IN]
                         zh       (:,:),          & ! [IN]
                         WPRECH   (:,:),          & ! [IN]
                         kcell_lim(:),            & ! [IN]
                         zdis0    (:,:),          & ! [OUT]
                         kcell    (:,:),          & ! [OUT]
                         kcell_max(:),            & ! [OUT]
                         dt                       ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_05z')
call PROF_rapstart('____Precip_Transport_06')
#endif

       !----- mass
!OCL XFILL
       !$acc kernels pcopy(rhoq) pcopyin(rhogq,rgs) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kdim,rhoq,rhogq,rgs)
       do k  = 1, kdim
       do ij = 1, ijdim
          rhoq(ij,k) = rhogq(ij,k,nq) * rgs(ij,k)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_06')
call PROF_rapstart('____Precip_Transport_07z')
#endif

       call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                dz       (:,:),          & ! [IN]
                                rhoq     (:,:),          & ! [IN]
                                zdis0    (:,:),          & ! [IN]
                                kcell    (:,:),          & ! [IN]
                                kcell_max(:),            & ! [IN]
                                flux_rhoq(:,:)           ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_07z')
call PROF_rapstart('____Precip_Transport_08')
#endif

!OCL XFILL
       !$acc kernels pcopy(drhogq,flux_rhoq) pcopyin(GRD_rdgz) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kmin,kmax,drhogq,flux_rhoq,GRD_rdgz)
!OCL NOPREFETCH
       do k  = kmin, kmax
       do ij = 1, ijdim
          drhogq(ij,k,nq) = -( flux_rhoq(ij,k+1)-flux_rhoq(ij,k) ) * GRD_rdgz(k) ! G^1/2 / dz = 1 / dgsi
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_08')
#endif

       ! Only for mass tracer
       if (       nq >= NQW_STR &
            .AND. nq <= NQW_END ) then

          !--- internal energy
          !--- potential energy
          !--- horizontal kinetic energy
          !--- momentum u
          !--- momentum v

#ifdef MORETIMER
call PROF_rapstart('____Precip_Transport_09')
#endif

!OCL XFILL
          !$acc kernels pcopy(rhogkin_h,rhoq,rhoqe,rhoqphi,rhoqkin_h,rhoqu,rhoqv) async(1) &
          !$acc pcopyin(rhog,tem,vx,vy,vz,ix,iy,iz,jx,jy,jz)
          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhoqe,rhoqphi,rhoqkin_h,rhoqu,rhoqv,           &
          !$omp        rhoq,tem,vx,vy,vz,z,ix,iy,iz,jx,jy,jz,rhogkin_h,rhog,CVW,GRAV)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhoqe    (ij,k) = rhoq(ij,k) * CVW(nq) * tem(ij,k)
             rhoqphi  (ij,k) = rhoq(ij,k) * GRAV * z(ij,k)
             rhoqkin_h(ij,k) = rhoq(ij,k) * rhogkin_h(ij,k) / rhog(ij,k)
             rhoqu    (ij,k) = rhoq(ij,k) * ( vx(ij,k) * ix(ij) &
                                            + vy(ij,k) * iy(ij) &
                                            + vz(ij,k) * iz(ij) )
             rhoqv    (ij,k) = rhoq(ij,k) * ( vx(ij,k) * jx(ij) &
                                            + vy(ij,k) * jy(ij) &
                                            + vz(ij,k) * jz(ij) )
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_09')
call PROF_rapstart('____Precip_Transport_10z')
#endif

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                   dz       (:,:),          & ! [IN]
                                   rhoqe    (:,:),          & ! [IN]
                                   zdis0    (:,:),          & ! [IN]
                                   kcell    (:,:),          & ! [IN]
                                   kcell_max(:),            & ! [IN]
                                   flux_rhoqe(:,:)          ) ! [OUT]

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                   dz       (:,:),          & ! [IN]
                                   rhoqphi  (:,:),          & ! [IN]
                                   zdis0    (:,:),          & ! [IN]
                                   kcell    (:,:),          & ! [IN]
                                   kcell_max(:),            & ! [IN]
                                   flux_rhoqphi(:,:)        ) ! [OUT]

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                   dz       (:,:),          & ! [IN]
                                   rhoqkin_h(:,:),         & ! [IN]
                                   zdis0    (:,:),          & ! [IN]
                                   kcell    (:,:),          & ! [IN]
                                   kcell_max(:),            & ! [IN]
                                   flux_rhoqkin_h(:,:)      ) ! [OUT]

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                   dz       (:,:),          & ! [IN]
                                   rhoqu    (:,:),          & ! [IN]
                                   zdis0    (:,:),          & ! [IN]
                                   kcell    (:,:),          & ! [IN]
                                   kcell_max(:),            & ! [IN]
                                   flux_rhoqu(:,:)          ) ! [OUT]

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & ! [IN]
                                   dz       (:,:),          & ! [IN]
                                   rhoqv    (:,:),          & ! [IN]
                                   zdis0    (:,:),          & ! [IN]
                                   kcell    (:,:),          & ! [IN]
                                   kcell_max(:),            & ! [IN]
                                   flux_rhoqv(:,:)          ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_10z')
call PROF_rapstart('____Precip_Transport_11')
#endif

          !$acc kernels pcopy(drhoge,drhogphi,drhogkin_h,drhogu,drhogv) async(1) &
          !$acc pcopyin(flux_rhoqe,flux_rhoqphi,flux_rhoqkin_h,flux_rhoqu,flux_rhoqv,GRD_rdgz)
          !$omp parallel default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kmin,kmax,drhoge,drhogphi,drhogkin_h,drhogu,drhogv,         &
          !$omp        flux_rhoqe,flux_rhoqphi,flux_rhoqkin_h,flux_rhoqu,flux_rhoqv,GRD_rdgz)
          do k = kmin, kmax
             !$omp do
             do ij = 1, ijdim
                drhoge    (ij,k) = drhoge    (ij,k) -( flux_rhoqe    (ij,k+1)-flux_rhoqe    (ij,k) ) * GRD_rdgz(k)
#ifndef _OPENACC
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
#endif
                drhogphi  (ij,k) = drhogphi  (ij,k) -( flux_rhoqphi  (ij,k+1)-flux_rhoqphi  (ij,k) ) * GRD_rdgz(k)
#ifndef _OPENACC
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
#endif
                drhogkin_h(ij,k) = drhogkin_h(ij,k) -( flux_rhoqkin_h(ij,k+1)-flux_rhoqkin_h(ij,k) ) * GRD_rdgz(k)
#ifndef _OPENACC
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
#endif
                drhogu    (ij,k) = drhogu    (ij,k) -( flux_rhoqu    (ij,k+1)-flux_rhoqu    (ij,k) ) * GRD_rdgz(k)
#ifndef _OPENACC
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
#endif
                drhogv    (ij,k) = drhogv    (ij,k) -( flux_rhoqv    (ij,k+1)-flux_rhoqv    (ij,k) ) * GRD_rdgz(k)
             enddo
             !$omp end do
          enddo
          !$omp end parallel
          !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_11')
call PROF_rapstart('____Precip_Transport_12z')
#endif

          !--- vertical kinetic energy
          !--- moment w
          !--- half level
          call vadv1d_prep( ijdim, kdim, kmin+1, kmax, & ! [IN]
                            dzh         (:,:),         & ! [IN]
                            dzh2        (:,:),         & ! [IN]
                            z_shift     (:,:),         & ! [IN]
                            WPREC       (:,:,nq),      & ! [IN]
                            kcellh_lim  (:),           & ! [IN]
                            zdis0h      (:,:),         & ! [OUT]
                            kcellh      (:,:),         & ! [OUT]
                            kcellh_max  (:),           & ! [OUT]
                            dt                         ) ! [IN]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_12z')
call PROF_rapstart('____Precip_Transport_13')
#endif

!OCL XFILL
          !$acc kernels pcopy(rhoqkin_v,rhoqw) pcopyin(rhogw,rgsh) async(1)
          !$omp parallel do default(none),private(ij,k,qh), &
          !$omp shared(nq,ijdim,kmin,kmax,q,rhoqkin_v,rhoqw,rhogkin_v,rhogw,rgsh,GRD_afact,GRD_bfact)
          do k  = kmin+1, kmax
          do ij = 1, ijdim
             qh = ( GRD_afact(k) * q(ij,k  ,nq) &
                  + GRD_bfact(k) * q(ij,k-1,nq) )

             rhoqkin_v(ij,k) = qh * rhogkin_v(ij,k) * rgsh(ij,k)
             rhoqw    (ij,k) = qh * rhogw    (ij,k) * rgsh(ij,k)
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_13')
call PROF_rapstart('____Precip_Transport_14z')
#endif

          call vadv1d_getflux_new( ijdim, kdim, kmin+1, kmax, & ! [IN]
                                   dzh       (:,:),           & ! [IN]
                                   rhoqkin_v (:,:),           & ! [IN]
                                   zdis0h    (:,:),           & ! [IN]
                                   kcellh    (:,:),           & ! [IN]
                                   kcellh_max(:),             & ! [IN]
                                   flux_rhoqkin_v(:,:)        ) ! [OUT]

          call vadv1d_getflux_new( ijdim, kdim, kmin+1, kmax, & ! [IN]
                                   dzh       (:,:),           & ! [IN]
                                   rhoqw     (:,:),           & ! [IN]
                                   zdis0h    (:,:),           & ! [IN]
                                   kcellh    (:,:),           & ! [IN]
                                   kcellh_max(:),             & ! [IN]
                                   flux_rhoqw(:,:)            ) ! [OUT]

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_14z')
call PROF_rapstart('____Precip_Transport_15')
#endif

          !$acc kernels pcopy(drhogkin_v,drhogw) pcopyin(flux_rhoqkin_v,flux_rhoqw,GRD_rdgzh) async(1)
          !$omp parallel default(none),private(ij,k), &
          !$omp shared(ijdim,kmin,kmax,drhogkin_v,drhogw,flux_rhoqkin_v,flux_rhoqw,GRD_rdgzh)
          do k  = kmin+1, kmax
             !$omp do
             do ij = 1, ijdim
                drhogkin_v(ij,k) = drhogkin_v(ij,k) -( flux_rhoqkin_v(ij,k)-flux_rhoqkin_v(ij,k-1) ) * GRD_rdgzh(k)
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
                drhogw    (ij,k) = drhogw    (ij,k) -( flux_rhoqw    (ij,k)-flux_rhoqw    (ij,k-1) ) * GRD_rdgzh(k)
             enddo
             !$omp end do
          enddo
          !$omp end parallel
          !$acc end kernels

          ! precipitation on the ground
          if ( nq == I_QC .OR. nq == I_QR ) then

             !$omp parallel default(none),private(ij,k), &
             !$omp shared(ijdim,kmin,kmax,precip_flux,precip,flux_rhoq,dt)

             !$acc kernels pcopy(precip) pcopyin(flux_rhoq) async(1)
             !$omp do
             do ij = 1, ijdim
                precip(ij,1) = precip(ij,1) - flux_rhoq(ij,kmin) / dt
             enddo
             !$omp end do nowait
             !$acc end kernels

             !$acc kernels pcopy(precip_flux) pcopyin(flux_rhoq) async(1)
             !$omp do
             do k  = kmin, kmax
             do ij = 1, ijdim
                precip_flux(ij,k,1) = precip_flux(ij,k,1) - flux_rhoq(ij,k) / dt
             enddo
             enddo
             !$omp end do nowait
             !$acc end kernels

             !$omp end parallel

          elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then

             !$omp parallel default(none),private(ij,k), &
             !$omp shared(ijdim,kmin,kmax,precip_flux,precip,precip_lh_heat,flux_rhoq,dt,LHF)

             !$acc kernels pcopy(precip,precip_lh_heat) pcopyin(flux_rhoq) async(1)
             !$omp do
             do ij = 1, ijdim
                precip(ij,2) = precip(ij,2) - flux_rhoq(ij,kmin) / dt

                precip_lh_heat(ij) = precip_lh_heat(ij) + flux_rhoq(ij,kmin) * LHF / dt
             enddo
             !$omp end do nowait
             !$acc end kernels

             !$acc kernels pcopy(precip_flux) pcopyin(flux_rhoq) async(1)
             !$omp do
             do k  = kmin, kmax
             do ij = 1, ijdim
                precip_flux(ij,k,2) = precip_flux(ij,k,2) - flux_rhoq(ij,k) / dt
             enddo
             enddo
             !$omp end do nowait
             !$acc end kernels

             !$omp end parallel

          endif

          if ( present(precip_trc) ) then ! [Add] 2012/02/01 T.Seiki
             !$acc kernels pcopy(precip_trc) pcopyin(flux_rhoq) async(1)
             !$omp parallel do default(none),private(ij), &
             !$omp shared(nq,ijdim,kmin,precip_trc,flux_rhoq,dt)
             do ij = 1, ijdim
                precip_trc(ij,nq) = precip_trc(ij,nq) - flux_rhoq(ij,kmin) / dt
             enddo
             !$omp end parallel do
             !$acc end kernels
          endif

          !$acc kernels pcopy(precip_rhoe,precip_rhophi,precip_rhokin) async(1) &
          !$acc pcopyin(flux_rhoqe,flux_rhoqphi,flux_rhoqkin_h,flux_rhoqkin_v)
          !$omp parallel do default(none),private(ij), &
          !$omp shared(ijdim,kmin,precip_rhoe,precip_rhophi,precip_rhokin,     &
          !$omp        flux_rhoqe,flux_rhoqphi,flux_rhoqkin_h,flux_rhoqkin_v,dt)
          do ij = 1, ijdim
             precip_rhoe  (ij) = precip_rhoe  (ij) - flux_rhoqe    (ij,kmin) / dt
             precip_rhophi(ij) = precip_rhophi(ij) - flux_rhoqphi  (ij,kmin) / dt
             precip_rhokin(ij) = precip_rhokin(ij) - flux_rhoqkin_h(ij,kmin) / dt &
                                                   - flux_rhoqkin_v(ij,kmin) / dt
          enddo
          !$omp end parallel do
          !$acc end kernels

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_15')
#endif

       endif

       ! precipitation isotopes on the ground [Add] M.Tanoue 20200907
       if ( opt_water_isotope ) then
          if    ( nq == ISO1_QC .OR. nq == ISO1_QR ) then
             do ij = 1, ijdim
                precip_ISO1(ij,1) = precip_ISO1(ij,1) - flux_rhoq(ij,kmin) / dt
             enddo
          elseif( nq == ISO2_QC .OR. nq == ISO2_QR ) then
             do ij = 1, ijdim
                precip_ISO2(ij,1) = precip_ISO2(ij,1) - flux_rhoq(ij,kmin) / dt
             enddo
          elseif( nq == ISO1_QI .OR. nq == ISO1_QS .OR. nq == ISO1_QG ) then
             do ij = 1, ijdim
                precip_ISO1(ij,2) = precip_ISO1(ij,2) - flux_rhoq(ij,kmin) / dt
             enddo
          elseif( nq == ISO2_QI .OR. nq == ISO2_QS .OR. nq == ISO2_QG ) then
             do ij = 1, ijdim
                precip_ISO2(ij,2) = precip_ISO2(ij,2) - flux_rhoq(ij,kmin) / dt
             enddo
          endif
       endif
       !!$acc wait(1)
    enddo ! tracer LOOP

#ifdef MORETIMER
call PROF_rapstart('____Precip_Transport_16')
#endif

    ! Change in internal energy comes from precipitation and dissipation of kinetic energy due to drag force.
    ! See Ooyama(2001) (3.13)

    !$acc kernels pcopy(rhoge,rhogkin_h,rhogkin_v,rhogvx,rhogvy,rhogvz,rhogw) async(1) &
    !$acc pcopyin(drhoge,drhogphi,drhogkin_h,drhogkin_v,drhogu,drhogv,drhogw,ix,iy,iz,jx,jy,jz)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhoge,rhogkin_h,rhogkin_v,rhogvx,rhogvy,rhogvz,rhogw,            &
    !$omp        drhoge,drhogphi,drhogkin_h,drhogkin_v,drhogu,drhogv,drhogw,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       rhoge    (ij,k) = rhoge    (ij,k) + drhoge    (ij,k) + drhogphi(ij,k)
       rhogkin_h(ij,k) = rhogkin_h(ij,k) + drhogkin_h(ij,k)
       rhogkin_v(ij,k) = rhogkin_v(ij,k) + drhogkin_v(ij,k)
       rhogvx   (ij,k) = rhogvx   (ij,k) + drhogu    (ij,k) * ix(ij) + drhogv(ij,k) * jx(ij)
       rhogvy   (ij,k) = rhogvy   (ij,k) + drhogu    (ij,k) * iy(ij) + drhogv(ij,k) * jy(ij)
       rhogvz   (ij,k) = rhogvz   (ij,k) + drhogu    (ij,k) * iz(ij) + drhogv(ij,k) * jz(ij)
       rhogw    (ij,k) = rhogw    (ij,k) + drhogw    (ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    do nq = 1, nqmax
       if ( nq >= NQW_STR .AND. nq <= NQW_END ) then
          !$acc kernels pcopy(rhog,rhogq,rhoge) pcopyin(drhogq,z) async(1)
          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhog,rhogq,rhoge,drhogq,z,GRAV)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhogq(ij,k,nq)
             rhog (ij,k)    = rhog (ij,k)    + drhogq(ij,k,nq)
             rhoge(ij,k)    = rhoge(ij,k)    - drhogq(ij,k,nq) * GRAV * z(ij,k)
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels
       else
          !$acc kernels pcopy(rhogq) pcopyin(drhogq) async(1)
          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhogq,drhogq)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhogq(ij,k,nq)
          enddo
          enddo
          !$omp end parallel do
          !$acc end kernels
       endif
       !!$acc wait(1)
    enddo

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_16')
call PROF_rapstart('____Precip_Transport_17')
#endif

    if ( PRCIP_TRN_ECORRECT == 'KIN2EIN' ) then

!OCL XFILL
       !$acc kernels async(1) pcopy(tmp2) pcopyin(rhogkin_h,rhogkin_v,GRD_cfact,GRD_dfact)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,tmp2,rhogkin_h,rhogkin_v,GRD_cfact,GRD_dfact)
       do k  = kmin, kmax
       do ij = 1, ijdim
          tmp2(ij,k) = rhogkin_h(ij,k) + ( GRD_cfact(k) * rhogkin_v(ij,k+1) &
                                         + GRD_dfact(k) * rhogkin_v(ij,k  ) )
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       call CNVVAR_rhogkin_in( ijdim,          & ! [IN]
                               kdim,           & ! [IN]
                               rhog   (:,:),   & ! [IN]
                               rhogvx (:,:),   & ! [IN]
                               rhogvy (:,:),   & ! [IN]
                               rhogvz (:,:),   & ! [IN]
                               rhogw  (:,:),   & ! [IN]
                               C2Wfact(:,:,:), & ! [IN]
                               W2Cfact(:,:,:), & ! [IN]
                               tmp    (:,:),   & ! [OUT]
                               tmp_h  (:,:),   & ! [OUT]
                               tmp_v  (:,:)    ) ! [OUT]

       !$acc kernels pcopy(rhoge) pcopyin(tmp2,tmp) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhoge,tmp2,tmp)
       do k  = kmin, kmax
       do ij = 1, ijdim
          rhoge(ij,k) = rhoge(ij,k) + ( tmp2(ij,k) - tmp(ij,k) )
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

    elseif( PRCIP_TRN_ECORRECT == 'KIN2KIN' ) then

       !$acc kernels async(1) pcopy(kin_h0,kin_h,vx_t,vy_t,vz_t) pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogkin_h)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,kin_h0,kin_h,vx_t,vy_t,vz_t,rhog,rhogvx,rhogvy,rhogvz,rhogkin_h)
       do k  = kmin, kmax
       do ij = 1, ijdim
          kin_h0(ij,k) = rhogkin_h(ij,k) / rhog(ij,k)
          vx_t  (ij,k) = rhogvx   (ij,k) / rhog(ij,k)
          vy_t  (ij,k) = rhogvy   (ij,k) / rhog(ij,k)
          vz_t  (ij,k) = rhogvz   (ij,k) / rhog(ij,k)

          kin_h (ij,k) = 0.5_RP * ( vx_t(ij,k)*vx_t(ij,k) &
                                  + vy_t(ij,k)*vy_t(ij,k) &
                                  + vz_t(ij,k)*vz_t(ij,k) )
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       !$acc kernels async(1) pcopy(vx_t,vy_t,vz_t) pcopyin(kin_h0,kin_h)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,kin_h0,kin_h,vx_t,vy_t,vz_t,EPS)
       do k  = kmin, kmax
       do ij = 1, ijdim
          if ( kin_h(ij,k) > EPS ) then
             vx_t(ij,k) = vx_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
             vy_t(ij,k) = vy_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
             vz_t(ij,k) = vz_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
          endif
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

!OCL XFILL
       !$acc kernels async(1) pcopy(rhogvx,rhogvy,rhogvz) pcopyin(vx_t,vy_t,vz_t,rhog)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogvx,rhogvy,rhogvz,vx_t,vy_t,vz_t,rhog)
       do k  = kmin, kmax
       do ij = 1, ijdim
          rhogvx(ij,k) = vx_t(ij,k) * rhog(ij,k)
          rhogvy(ij,k) = vy_t(ij,k) * rhog(ij,k)
          rhogvz(ij,k) = vz_t(ij,k) * rhog(ij,k)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

!OCL XFILL
       !$acc kernels pcopy(rhogh) pcopyin(rhog,C2Wfact) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogh,rhog,C2Wfact)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          rhogh(ij,k) = ( C2Wfact(ij,k,1) * rhog(ij,k  ) &
                        + C2Wfact(ij,k,2) * rhog(ij,k-1) )
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       !$acc kernels pcopy(kin_v0,kin_v,w_t) present(rhogh) pcopyin(rhogw,rhogkin_v) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,kin_v0,kin_v,w_t,rhogh,rhogw,rhogkin_v)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          kin_v0(ij,k) = rhogkin_v(ij,k) / rhogh(ij,k)
          w_t   (ij,k) = rhogw    (ij,k) / rhogh(ij,k)
          kin_v (ij,k) = 0.5_RP * w_t(ij,k)**2
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

       !$acc kernels pcopy(w_t) pcopyin(kin_v0,kin_v) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,w_t,kin_v0,kin_v,EPS)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          if ( kin_v(ij,k) > EPS ) then
             w_t(ij,k) = w_t(ij,k) * sqrt( abs( kin_v0(ij,k) / kin_v(ij,k) ) )
          endif
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

!OCL XFILL
       !$acc kernels pcopy(rhogw) pcopyin(w_t,rhogh) async(1)
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogw,w_t,rhogh)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          rhogw(ij,k) = w_t(ij,k) * rhogh(ij,k)
       enddo
       enddo
       !$omp end parallel do
       !$acc end kernels

    else
       write(*,*) 'Error in PRCIP_TRN_ECORRECT: ', trim(PRCIP_TRN_ECORRECT)
    endif

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_17')
call PROF_rapstart('____Precip_Transport_18')
#endif

!OCL XFILL
    !$acc kernels pcopy(q) pcopyin(rhog,rhogq) async(1)
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,nqmax,q,rhog,rhogq)
!OCL NOPREFETCH
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(rho,ein) pcopyin(rhog,rhoge,gsgam2) async(1)
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,ein,rhog,rhoge,gsgam2)
!OCL NOPREFETCH
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog (ij,k) / gsgam2(ij,k)
       ein(ij,k) = rhoge(ij,k) / rhog(ij,k)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels
    call THRMDYN_qd    ( ijdim, kdim, q(:,:,:), qd(:,:) )
    call THRMDYN_tempre( ijdim, kdim, ein(:,:), rho(:,:), q(:,:,:), tem(:,:), pre(:,:) )

    cnvvar_wait_on = .true.
    thrmdyn_wait_on = .true.

#ifdef MORETIMER
call PROF_rapend  ('____Precip_Transport_18')
#endif

    call PROF_rapend  ('____Precip_Transport',2)

    !$acc end data
    return
  end subroutine precip_transport_new

  !-----------------------------------------------------------------------------
  subroutine precip_transport_euler( &
       ijdim,              &
       rhog,               &
       rhogvx,             &
       rhogvy,             &
       rhogvz,             &
       rhogw,              &
       rhoge,              &
       rhogq,              &
       rho,                &
       tem,                &
       pre,                &
       vx,                 &
       vy,                 &
       vz,                 &
       w,                  &
       q,                  &
       qd,                 &
       z,                  &
       zh,                 &
       WPREC,              &
       precipitating_flag, &
       precip,             &
       precip_rhoe,        &
       precip_lh_heat,     &
       precip_rhophi,      &
       precip_rhokin,      &
       precip_flux,        &
       gsgam2,             &
       gsgam2h,            &
       rgs,                &
       rgsh,               &
       ix,                 &
       iy,                 &
       iz,                 &
       jx,                 &
       jy,                 &
       jz,                 &
       dt,                 &
       precip_trc          )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_cnst, only: &
       CNST_EPS,  &
       CNST_GRAV
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh, &
       GRD_afact, &
       GRD_bfact, &
       GRD_cfact, &
       GRD_dfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX, &
       NQW_STR,           &
       NQW_END,           &
       I_QC,              &
       I_QR,              &
       I_QI,              &
       I_QS,              &
       I_QG,              &
       CVW,               &
       LHF,               &
       PRCIP_TRN_ECORRECT
    use mod_cnvvar, only: &
       CNVVAR_rhogkin_in
    use mod_thrmdyn, only: &
       THRMDYN_qd,    &
       THRMDYN_tempre
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(inout) :: rhog              (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz            (ijdim,kdim)
    real(RP), intent(inout) :: rhogw             (ijdim,kdim)
    real(RP), intent(inout) :: rhoge             (ijdim,kdim)
    real(RP), intent(inout) :: rhogq             (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: rho               (ijdim,kdim)
    real(RP), intent(inout) :: tem               (ijdim,kdim)
    real(RP), intent(inout) :: pre               (ijdim,kdim)
    real(RP), intent(in)    :: vx                (ijdim,kdim)
    real(RP), intent(in)    :: vy                (ijdim,kdim)
    real(RP), intent(in)    :: vz                (ijdim,kdim)
    real(RP), intent(in)    :: w                 (ijdim,kdim)
    real(RP), intent(inout) :: q                 (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd                (ijdim,kdim)
    real(RP), intent(in)    :: z                 (ijdim,kdim)
    real(RP), intent(in)    :: zh                (ijdim,kdim)
    real(RP), intent(in)    :: WPREC             (ijdim,kdim,nqmax)
    logical,  intent(in)    :: precipitating_flag(nqmax)
    real(RP), intent(out)   :: precip            (ijdim,2)
    real(RP), intent(out)   :: precip_rhoe       (ijdim)
    real(RP), intent(out)   :: precip_lh_heat    (ijdim)
    real(RP), intent(out)   :: precip_rhophi     (ijdim)
    real(RP), intent(out)   :: precip_rhokin     (ijdim)
    real(RP), intent(out)   :: precip_flux       (ijdim,kdim,2)
    real(RP), intent(in)    :: gsgam2            (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h           (ijdim,kdim)
    real(RP), intent(in)    :: rgs               (ijdim,kdim)
    real(RP), intent(in)    :: rgsh              (ijdim,kdim)
    real(RP), intent(in)    :: ix                (ijdim)
    real(RP), intent(in)    :: iy                (ijdim)
    real(RP), intent(in)    :: iz                (ijdim)
    real(RP), intent(in)    :: jx                (ijdim)
    real(RP), intent(in)    :: jy                (ijdim)
    real(RP), intent(in)    :: jz                (ijdim)
    real(DP), intent(in)    :: dt
    real(RP), intent(inout), optional :: precip_trc(ijdim,nqmax) ! [Add] 2012/02/01 T.Seiki

    real(RP) :: rhogkin       (ijdim,kdim)
    real(RP) :: rhogkin_h     (ijdim,kdim)
    real(RP) :: rhogkin_v     (ijdim,kdim)

    real(RP) :: flux_rhoq     (ijdim,kdim)
    real(RP) :: flux_rhoqe    (ijdim,kdim)
    real(RP) :: flux_rhoqphi  (ijdim,kdim)
    real(RP) :: flux_rhoqkin_h(ijdim,kdim)
    real(RP) :: flux_rhoqkin_v(ijdim,kdim)
    real(RP) :: flux_rhoqu    (ijdim,kdim)
    real(RP) :: flux_rhoqv    (ijdim,kdim)
    real(RP) :: flux_rhoqw    (ijdim,kdim)

    real(RP) :: drhogq        (ijdim,kdim,nqmax)
    real(RP) :: drhoge        (ijdim,kdim)
    real(RP) :: drhogphi      (ijdim,kdim)
    real(RP) :: drhogkin_h    (ijdim,kdim)
    real(RP) :: drhogkin_v    (ijdim,kdim)
    real(RP) :: drhogu        (ijdim,kdim)
    real(RP) :: drhogv        (ijdim,kdim)
    real(RP) :: drhogw        (ijdim,kdim)

    real(RP) :: kin_h0        (ijdim,kdim)
    real(RP) :: kin_h         (ijdim,kdim)
    real(RP) :: vx_t          (ijdim,kdim)
    real(RP) :: vy_t          (ijdim,kdim)
    real(RP) :: vz_t          (ijdim,kdim)
    real(RP) :: kin_v0        (ijdim,kdim)
    real(RP) :: kin_v         (ijdim,kdim)
    real(RP) :: w_t           (ijdim,kdim)

    real(RP) :: qh
    real(RP) :: rhogh         (ijdim,kdim)
    real(RP) :: ein           (ijdim,kdim)
    real(RP) :: tmp           (ijdim,kdim)
    real(RP) :: tmp2          (ijdim,kdim)
    real(RP) :: tmp_h         (ijdim,kdim)
    real(RP) :: tmp_v         (ijdim,kdim)
    real(RP) :: C2Wfact       (ijdim,kdim,2)
    real(RP) :: W2Cfact       (ijdim,kdim,2)

    real(RP) :: dz            (ijdim,kdim)
    real(RP) :: dzh           (ijdim,kdim)
    real(RP) :: z_shift       (ijdim,kdim)

    real(RP) :: dz_min, wprec_max
    real(RP) :: dt_precip    ! dt for small step
    real(RP) :: factor
    integer  :: nstep_precip ! #  of small step
    integer  :: step

    real(RP) :: GRAV, EPS
    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    call PROF_rapstart('____Precip_Transport',2)

    EPS  = CNST_EPS
    GRAV = CNST_GRAV

!OCL XFILL
    do k  = kmin, kmax
    do ij = 1, ijdim
       dz (ij,k) = zh(ij,k+1) - zh(ij,k)
       dzh(ij,k) = z (ij,k)   - z (ij,k-1)
    enddo
    enddo
!OCL XFILL
    do ij = 1, ijdim
       dzh(ij,kmin) = z(ij,kmin) - zh(ij,kmin)
    enddo
!OCL XFILL
    do k  = kmin, kmax+1
    do ij = 1, ijdim
       z_shift(ij,k) = z(ij,k-1)
    enddo
    enddo

    dz_min = maxval(dz(:,:))
    do k  = kmin, kmax
    do ij = 1, ijdim
       dz_min = min( dz_min, dz(ij,k), dzh(ij,k) )
    enddo
    enddo

    wprec_max = 0.0_RP
    do nq = 1, nqmax
       if ( .NOT. precipitating_flag(nq) ) cycle

       do k  = kmin, kmax
       do ij = 1, ijdim
          wprec_max = max( wprec_max, -WPREC(ij,k,nq) ) ! w-precip is always negative
       enddo
       enddo
    enddo

    nstep_precip = int( dt * wprec_max / dz_min / 0.5_RP ) + 1 ! CFL=0.5
    factor       = 1.0_DP / real(nstep_precip,kind=DP)
    dt_precip    = dt * factor

    if ( nstep_precip > 1 ) then
       write(IO_FID_LOG,'(1x,4(A,F6.1))') '*** dt(MP) => dt(PREC) :', dt,        &
                                                               ' =>', dt_precip, &
                                                        ', DZ min =', dz_min,    &
                                                        ', Vt max =', wprec_max
    endif

    do k = kmin, kmax+1
       C2Wfact(:,k,1) = GRD_afact(k) / GSGAM2(:,k  ) * GSGAM2H(:,k)
       C2Wfact(:,k,2) = GRD_bfact(k) / GSGAM2(:,k-1) * GSGAM2H(:,k)
    enddo
    C2Wfact(:,kmin-1,1) = 0.0_RP
    C2Wfact(:,kmin-1,2) = 0.0_RP

    do k = kmin-1, kmax
       W2Cfact(:,k,1) = GRD_cfact(k) * GSGAM2(:,k) / GSGAM2H(:,k+1)
       W2Cfact(:,k,2) = GRD_dfact(k) * GSGAM2(:,k) / GSGAM2H(:,k  )
    enddo
    W2Cfact(:,kmax+1,1) = 0.0_RP
    W2Cfact(:,kmax+1,2) = 0.0_RP

!OCL XFILL
    precip        (:,:)   = 0.0_RP ! rain/snow
!OCL XFILL
    precip_rhoe   (:)     = 0.0_RP
!OCL XFILL
    precip_lh_heat(:)     = 0.0_RP
!OCL XFILL
    precip_rhophi (:)     = 0.0_RP
!OCL XFILL
    precip_rhokin (:)     = 0.0_RP
!OCL XFILL
    precip_flux   (:,:,:) = 0.0_RP ! rain/snow

    do step = 1, nstep_precip

       call CNVVAR_rhogkin_in( ijdim,            & ! [IN]
                               kdim,             & ! [IN]
                               rhog     (:,:),   & ! [IN]
                               rhogvx   (:,:),   & ! [IN]
                               rhogvy   (:,:),   & ! [IN]
                               rhogvz   (:,:),   & ! [IN]
                               rhogw    (:,:),   & ! [IN]
                               C2Wfact  (:,:,:), & ! [IN]
                               W2Cfact  (:,:,:), & ! [IN]
                               rhogkin  (:,:),   & ! [OUT]
                               rhogkin_h(:,:),   & ! [OUT]
                               rhogkin_v(:,:)    ) ! [OUT]

!OCL XFILL
       drhoge    (:,:) = 0.0_RP
!OCL XFILL
       drhogphi  (:,:) = 0.0_RP
!OCL XFILL
       drhogkin_h(:,:) = 0.0_RP
!OCL XFILL
       drhogkin_v(:,:) = 0.0_RP
!OCL XFILL
       drhogu    (:,:) = 0.0_RP
!OCL XFILL
       drhogv    (:,:) = 0.0_RP
!OCL XFILL
       drhogw    (:,:) = 0.0_RP

       do nq = 1, nqmax

          if ( .NOT. precipitating_flag(nq) ) then
!OCL XFILL
             drhogq(:,:,nq) = 0.0_RP
             cycle
          endif
          !----- mass
!OCL XFILL
          do k  = kmin+1, kmax
          do ij = 1, ijdim
             flux_rhoq(ij,k) = ( GRD_afact(k) * WPREC(ij,k  ,nq) &
                               + GRD_bfact(k) * WPREC(ij,k-1,nq) ) * rhogq(ij,k,nq) * rgs(ij,k)
          enddo
          enddo
!OCL XFILL
          do ij = 1, ijdim
             flux_rhoq(ij,kmin  ) = WPREC(ij,kmin,nq) * rhogq(ij,kmin,nq) * rgs(ij,kmin)
             flux_rhoq(ij,kmin-1) = 0.0_RP
             flux_rhoq(ij,kmax+1) = 0.0_RP
          enddo

!OCL XFILL
          do k = kmin, kmax
          do ij = 1,   ijdim
             drhogq(ij,k,nq) = - dt_precip * ( flux_rhoq(ij,k+1)-flux_rhoq(ij,k) ) * GRD_rdgz(k) ! G^1/2 / dz = 1 / dgsi
          enddo
          enddo

          ! Only for mass tracer
          if (       nq >= NQW_STR &
               .AND. nq <= NQW_END ) then

             !--- internal energy
             !--- potential energy
             !--- horizontal kinetic energy
             !--- momentum u
             !--- momentum v
!OCL XFILL
             do k  = kmin, kmax
             do ij = 1, ijdim
                flux_rhoqe    (ij,k) = flux_rhoq(ij,k) * CVW(nq) * tem(ij,k+1)
                flux_rhoqphi  (ij,k) = flux_rhoq(ij,k) * GRAV * z(ij,k+1)
                flux_rhoqkin_h(ij,k) = flux_rhoq(ij,k) * rhogkin_h(ij,k+1) / rhog(ij,k)
                flux_rhoqu    (ij,k) = flux_rhoq(ij,k) * ( vx(ij,k+1) * ix(ij) &
                                                         + vy(ij,k+1) * iy(ij) &
                                                         + vz(ij,k+1) * iy(ij) )
                flux_rhoqv    (ij,k) = flux_rhoq(ij,k) * ( vx(ij,k+1) * jx(ij) &
                                                         + vy(ij,k+1) * jy(ij) &
                                                         + vz(ij,k+1) * jy(ij) )
             enddo
             enddo
!OCL XFILL
             do ij = 1, ijdim
                flux_rhoqe    (ij,kmax+1) = 0.0_RP
                flux_rhoqphi  (ij,kmax+1) = 0.0_RP
                flux_rhoqkin_h(ij,kmax+1) = 0.0_RP
                flux_rhoqu    (ij,kmax+1) = 0.0_RP
                flux_rhoqv    (ij,kmax+1) = 0.0_RP
             enddo

             do k  = kmin, kmax
             do ij = 1, ijdim
                drhoge    (ij,k) = drhoge    (ij,k) - dt_precip * ( flux_rhoqe    (ij,k+1)-flux_rhoqe    (ij,k) ) * GRD_rdgz(k)
                drhogphi  (ij,k) = drhogphi  (ij,k) - dt_precip * ( flux_rhoqphi  (ij,k+1)-flux_rhoqphi  (ij,k) ) * GRD_rdgz(k)
                drhogkin_h(ij,k) = drhogkin_h(ij,k) - dt_precip * ( flux_rhoqkin_h(ij,k+1)-flux_rhoqkin_h(ij,k) ) * GRD_rdgz(k)
                drhogu    (ij,k) = drhogu    (ij,k) - dt_precip * ( flux_rhoqu    (ij,k+1)-flux_rhoqu    (ij,k) ) * GRD_rdgz(k)
                drhogv    (ij,k) = drhogv    (ij,k) - dt_precip * ( flux_rhoqv    (ij,k+1)-flux_rhoqv    (ij,k) ) * GRD_rdgz(k)
             enddo
             enddo

             !--- vertical kinetic energy
             !--- moment w
             !--- half level
!OCL XFILL
             do k  = kmin, kmax
             do ij = 1, ijdim
                qh = ( GRD_afact(k) * q(ij,k  ,nq) &
                     + GRD_bfact(k) * q(ij,k-1,nq) )

                flux_rhoqkin_v(ij,k) = WPREC(ij,k,nq) * qh * rhogkin_v(ij,k+1) * rgsh(ij,k)
                flux_rhoqw    (ij,k) = WPREC(ij,k,nq) * qh * rhogw    (ij,k+1) * rgsh(ij,k)
             enddo
             enddo

             do k  = kmin+1, kmax
             do ij = 1, ijdim
                drhogkin_v(ij,k) = drhogkin_v(ij,k) - dt_precip * ( flux_rhoqkin_v(ij,k)-flux_rhoqkin_v(ij,k-1) ) * GRD_rdgzh(k)
                drhogw    (ij,k) = drhogw    (ij,k) - dt_precip * ( flux_rhoqw    (ij,k)-flux_rhoqw    (ij,k-1) ) * GRD_rdgzh(k)
             enddo
             enddo

             ! ground precipitation
             if ( nq == I_QC .OR.  nq == I_QR ) then

                do ij = 1, ijdim
                   precip(ij,1) = precip(ij,1) - flux_rhoq(ij,kmin) * factor
                enddo

                do k  = kmin, kmax
                do ij = 1, ijdim
                   precip_flux(ij,k,1) = precip_flux(ij,k,1) - flux_rhoq(ij,k) * factor
                enddo
                enddo

             elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then

                do ij = 1, ijdim
                   precip(ij,2) = precip(ij,2) - flux_rhoq(ij,kmin) * factor
                   precip_lh_heat(ij) = precip_lh_heat(ij) + flux_rhoq(ij,kmin) * LHF * factor
                enddo

                do k  = kmin, kmax
                do ij = 1, ijdim
                   precip_flux(ij,k,2) = precip_flux(ij,k,2) - flux_rhoq(ij,k) * factor
                enddo
                enddo

             endif

             if ( present(precip_trc) ) then ! [Add] 2012/02/01 T.Seiki
                do ij = 1, ijdim
                   precip_trc(ij,nq) = precip_trc(ij,nq) - flux_rhoq(ij,kmin) * factor
                enddo
             endif

             do ij = 1, ijdim
                precip_rhoe  (ij) = precip_rhoe  (ij) - flux_rhoqe    (ij,kmin) * factor
                precip_rhophi(ij) = precip_rhophi(ij) - flux_rhoqphi  (ij,kmin) * factor
                precip_rhokin(ij) = precip_rhokin(ij) - flux_rhoqkin_h(ij,kmin) * factor &
                                                      - flux_rhoqkin_v(ij,kmin) * factor
             enddo

       endif

       enddo ! tracer LOOP

       ! Change in internal energy comes from precipitation and dissipation of kinetic energy due to drag force.
       ! See Ooyama(2001) (3.13)
       do k  = kmin, kmax
       do ij = 1, ijdim
          rhoge    (ij,k) = rhoge    (ij,k) + drhoge    (ij,k) + drhogphi(ij,k)
          rhogkin_h(ij,k) = rhogkin_h(ij,k) + drhogkin_h(ij,k)
          rhogkin_v(ij,k) = rhogkin_v(ij,k) + drhogkin_v(ij,k)
          rhogvx   (ij,k) = rhogvx   (ij,k) + drhogu    (ij,k) * ix(ij) + drhogv(ij,k) * jx(ij)
          rhogvy   (ij,k) = rhogvy   (ij,k) + drhogu    (ij,k) * iy(ij) + drhogv(ij,k) * jy(ij)
          rhogvz   (ij,k) = rhogvz   (ij,k) + drhogu    (ij,k) * iz(ij) + drhogv(ij,k) * jz(ij)
          rhogw    (ij,k) = rhogw    (ij,k) + drhogw    (ij,k)
       enddo
       enddo

       do nq = 1, nqmax
          if ( nq >= NQW_STR .AND. nq <= NQW_END ) then
             do k  = kmin, kmax
             do ij = 1, ijdim
                rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhogq(ij,k,nq)
                rhog (ij,k)    = rhog (ij,k)    + drhogq(ij,k,nq)
                rhoge(ij,k)    = rhoge(ij,k)    - drhogq(ij,k,nq) * GRAV * z(ij,k)
             enddo
             enddo
          else
             do k  = kmin, kmax
             do ij = 1, ijdim
                rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhogq(ij,k,nq)
             enddo
             enddo
          endif
       enddo

       if ( PRCIP_TRN_ECORRECT == 'KIN2EIN' ) then

!OCL XFILL
          do k  = kmin, kmax
          do ij = 1, ijdim
             tmp2(ij,k) = rhogkin_h(ij,k) + ( GRD_cfact(k) * rhogkin_v(ij,k+1) &
                                            + GRD_dfact(k) * rhogkin_v(ij,k  ) )
          enddo
          enddo

          call CNVVAR_rhogkin_in( ijdim,          & ! [IN]
                                  kdim,           & ! [IN]
                                  rhog   (:,:),   & ! [IN]
                                  rhogvx (:,:),   & ! [IN]
                                  rhogvy (:,:),   & ! [IN]
                                  rhogvz (:,:),   & ! [IN]
                                  rhogw  (:,:),   & ! [IN]
                                  C2Wfact(:,:,:), & ! [IN]
                                  W2Cfact(:,:,:), & ! [IN]
                                  tmp    (:,:),   & ! [OUT]
                                  tmp_h  (:,:),   & ! [OUT]
                                  tmp_v  (:,:)    ) ! [OUT]

          do k  = kmin, kmax
          do ij = 1, ijdim
             rhoge(ij,k) = rhoge(ij,k) + ( tmp2(ij,k) - tmp(ij,k) )
          enddo
          enddo

       elseif( PRCIP_TRN_ECORRECT == 'KIN2KIN' ) then

          do k  = kmin, kmax
          do ij = 1, ijdim
             kin_h0(ij,k) = rhogkin_h(ij,k) / rhog(ij,k)
             vx_t  (ij,k) = rhogvx   (ij,k) / rhog(ij,k)
             vy_t  (ij,k) = rhogvy   (ij,k) / rhog(ij,k)
             vz_t  (ij,k) = rhogvz   (ij,k) / rhog(ij,k)

             kin_h (ij,k) = 0.5_RP * ( vx_t(ij,k)*vx_t(ij,k) &
                                    + vy_t(ij,k)*vy_t(ij,k) &
                                    + vz_t(ij,k)*vz_t(ij,k) )
          enddo
          enddo

          do k  = kmin, kmax
          do ij = 1, ijdim
             if ( kin_h(ij,k) > EPS ) then
                vx_t(ij,k) = vx_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
                vy_t(ij,k) = vy_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
                vz_t(ij,k) = vz_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
             endif
          enddo
          enddo
!OCL XFILL
          do k  = kmin, kmax
          do ij = 1, ijdim
             rhogvx(ij,k) = vx_t(ij,k) * rhog(ij,k)
             rhogvy(ij,k) = vy_t(ij,k) * rhog(ij,k)
             rhogvz(ij,k) = vz_t(ij,k) * rhog(ij,k)
          enddo
          enddo
!OCL XFILL
          do k  = kmin, kmax+1
          do ij = 1, ijdim
             rhogh(ij,k) = ( C2Wfact(ij,k,1) * rhog(ij,k  ) &
                           + C2Wfact(ij,k,2) * rhog(ij,k-1) )
          enddo
          enddo
!OCL XFILL
          do k  = kmin, kmax+1
          do ij = 1, ijdim
             kin_v0(ij,k) = rhogkin_v(ij,k) / rhogh(ij,k)
             w_t   (ij,k) = rhogw    (ij,k) / rhogh(ij,k)
             kin_v (ij,k) = 0.5_RP * w_t(ij,k)**2
          enddo
          enddo

          do k  = kmin, kmax+1
          do ij = 1, ijdim
             if ( kin_v(ij,k) > EPS ) then
                w_t(ij,k) = w_t(ij,k) * sqrt( abs( kin_v0(ij,k) / kin_v(ij,k) ) )
             endif
          enddo
          enddo
!OCL XFILL
          do k  = kmin, kmax+1
          do ij = 1, ijdim
             rhogw(ij,k) = w_t(ij,k) * rhogh(ij,k)
          enddo
          enddo

       else
          write(*,*) 'Error in PRCIP_TRN_ECORRECT: ', trim(PRCIP_TRN_ECORRECT)
       endif

!OCL XFILL
       do nq = 1, nqmax
       do k  = 1, kdim
       do ij = 1, ijdim
          q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
       enddo
       enddo
       enddo
!OCL XFILL
       do k  = 1, kdim
       do ij = 1, ijdim
          rho(ij,k) = rhog (ij,k) / gsgam2(ij,k)
          ein(ij,k) = rhoge(ij,k) / rhog(ij,k)
       enddo
       enddo

       call THRMDYN_qd    ( ijdim, kdim, q(:,:,:), qd(:,:) )
       call THRMDYN_tempre( ijdim, kdim, ein(:,:), rho(:,:), q(:,:,:), tem(:,:), pre(:,:) )

    enddo ! small step LOOP

    call PROF_rapend  ('____Precip_Transport',2)

    return
  end subroutine precip_transport_euler

end module mod_precip_transport
