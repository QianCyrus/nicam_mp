!-------------------------------------------------------------------------------
!> Module variable conversion
!!
!! @par Description
!!         Conversion tools for prognostic variables
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_cnvvar
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

  use mod_runconf, only: &
     PRG_vmax0,  &
     I_RHOG,     &
     I_RHOGVX,   &
     I_RHOGVY,   &
     I_RHOGVZ,   &
     I_RHOGW,    &
     I_RHOGE,    &
     I_RHOGQstr, &
     I_RHOGQend, &
     DIAG_vmax0, &
     I_pre,      &
     I_tem,      &
     I_vx,       &
     I_vy,       &
     I_vz,       &
     I_w,        &
     I_qstr,     &
     I_qend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
!MINIAPP  public :: CNVVAR_prg2diag
!MINIAPP  public :: CNVVAR_diag2prg
  public :: CNVVAR_prg2diag_in
  public :: CNVVAR_diag2prg_in
!MINIAPP  public :: CNVVAR_rhogkin
  public :: CNVVAR_rhogkin_in
!MINIAPP  public :: CNVVAR_uv2vh
!MINIAPP  public :: CNVVAR_vh2uv
  public :: CNVVAR_vh2uv_in
!MINIAPP  public :: CNVVAR_uv2vh_2D
!MINIAPP  public :: CNVVAR_vh2uv_2D

!MINIAPP  interface CNVVAR_diag2prg
!MINIAPP     module procedure CNVVAR_diag2prg_SP
!MINIAPP     module procedure CNVVAR_diag2prg_DP
!MINIAPP  end interface CNVVAR_diag2prg

  interface CNVVAR_prg2diag_in
     module procedure CNVVAR_prg2diag_in_SP
     module procedure CNVVAR_prg2diag_in_DP
  end interface CNVVAR_prg2diag_in

  interface CNVVAR_diag2prg_in
     module procedure CNVVAR_diag2prg_in_SP
     module procedure CNVVAR_diag2prg_in_DP
  end interface CNVVAR_diag2prg_in

!MINIAPP  interface CNVVAR_rhogkin
!MINIAPP     module procedure CNVVAR_rhogkin_SP
!MINIAPP     module procedure CNVVAR_rhogkin_DP
!MINIAPP  end interface CNVVAR_rhogkin

  interface CNVVAR_rhogkin_in
     module procedure CNVVAR_rhogkin_in_SP
     module procedure CNVVAR_rhogkin_in_DP
  end interface CNVVAR_rhogkin_in

!MINIAPP  interface CNVVAR_uv2vh
!MINIAPP     module procedure CNVVAR_uv2vh_SP
!MINIAPP     module procedure CNVVAR_uv2vh_DP
!MINIAPP  end interface CNVVAR_uv2vh

!MINIAPP  interface CNVVAR_vh2uv
!MINIAPP     module procedure CNVVAR_vh2uv_SP
!MINIAPP     module procedure CNVVAR_vh2uv_DP
!MINIAPP  end interface CNVVAR_vh2uv

  interface CNVVAR_vh2uv_in
     module procedure CNVVAR_vh2uv_in_SP
     module procedure CNVVAR_vh2uv_in_DP
  end interface CNVVAR_vh2uv_in

!MINIAPP  interface CNVVAR_uv2vh_2D
!MINIAPP     module procedure CNVVAR_uv2vh_2D_SP
!MINIAPP     module procedure CNVVAR_uv2vh_2D_DP
!MINIAPP  end interface CNVVAR_uv2vh_2D

!MINIAPP  interface CNVVAR_vh2uv_2D
!MINIAPP     module procedure CNVVAR_vh2uv_2D_SP
!MINIAPP     module procedure CNVVAR_vh2uv_2D_DP
!MINIAPP  end interface CNVVAR_vh2uv_2D

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_prg2diag(&
!MINIAPP       prg,  prg_pl, &
!MINIAPP       diag, diag_pl )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_RGSGAM2,    &
!MINIAPP       VMTR_RGSGAM2_pl, &
!MINIAPP       VMTR_C2Wfact,    &
!MINIAPP       VMTR_C2Wfact_pl
!MINIAPP    use mod_runconf, only: &
!MINIAPP       PRG_vmax,  &
!MINIAPP       DIAG_vmax, &
!MINIAPP       TRC_vmax
!MINIAPP    use mod_thrmdyn, only: &
!MINIAPP       THRMDYN_tempre
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(RP), intent(in)  :: prg    (ADM_gall   ,ADM_kall,ADM_lall   ,PRG_vmax )
!MINIAPP    real(RP), intent(in)  :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
!MINIAPP    real(RP), intent(out) :: diag   (ADM_gall   ,ADM_kall,ADM_lall   ,DIAG_vmax)
!MINIAPP    real(RP), intent(out) :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)
!MINIAPP
!MINIAPP    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(RP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(RP) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(RP) :: rhog_h   (ADM_gall,   ADM_kall)
!MINIAPP    real(RP) :: rhog_h_pl(ADM_gall_pl,ADM_kall)
!MINIAPP
!MINIAPP    integer  :: g, k, l, iq
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(rho,ein,rhog_h) pcopyin(prg,VMTR_RGSGAM2,VMTR_C2Wfact) if(acc_is_present(diag))
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(diag,rho,ein) pcopyin(prg,VMTR_RGSGAM2) if(acc_is_present(diag))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       rho (g,k,l)      = prg(g,k,l,I_RHOG  ) * VMTR_RGSGAM2(g,k,l)
!MINIAPP       diag(g,k,l,I_vx) = prg(g,k,l,I_RHOGVX) / prg(g,k,l,I_RHOG)
!MINIAPP       diag(g,k,l,I_vy) = prg(g,k,l,I_RHOGVY) / prg(g,k,l,I_RHOG)
!MINIAPP       diag(g,k,l,I_vz) = prg(g,k,l,I_RHOGVZ) / prg(g,k,l,I_RHOG)
!MINIAPP       ein (g,k,l)      = prg(g,k,l,I_RHOGE ) / prg(g,k,l,I_RHOG)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(diag) pcopyin(prg) if(acc_is_present(diag))
!MINIAPP    do iq = 1, TRC_vmax
!MINIAPP    do l  = 1, ADM_lall
!MINIAPP    do k  = 1, ADM_kall
!MINIAPP    do g  = 1, ADM_gall
!MINIAPP       diag(g,k,l,DIAG_vmax0+iq) = prg(g,k,l,PRG_vmax0+iq) / prg(g,k,l,I_RHOG)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    call THRMDYN_tempre( ADM_gall,                  & ! [IN]
!MINIAPP                         ADM_kall,                  & ! [IN]
!MINIAPP                         ADM_lall,                  & ! [IN]
!MINIAPP                         ein (:,:,:),               & ! [IN]
!MINIAPP                         rho (:,:,:),               & ! [IN]
!MINIAPP                         diag(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                         diag(:,:,:,I_tem),         & ! [OUT]
!MINIAPP                         diag(:,:,:,I_pre)          ) ! [OUT]
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(diag,rhog_h) pcopyin(prg,VMTR_C2Wfact) if(acc_is_present(diag))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       !------ interpolation of rhog_h
!MINIAPP       do k = 2, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,k) = ( VMTR_C2Wfact(g,k,1,l) * prg(g,k  ,l,I_RHOG) &
!MINIAPP                        + VMTR_C2Wfact(g,k,2,l) * prg(g,k-1,l,I_RHOG) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,1) = rhog_h(g,2)
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          diag(g,k,l,I_w) = prg(g,k,l,I_RHOGW) / rhog_h(g,k)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc end data
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          rho_pl (g,k,l)      = prg_pl(g,k,l,I_RHOG  ) * VMTR_RGSGAM2_pl(g,k,l)
!MINIAPP          diag_pl(g,k,l,I_vx) = prg_pl(g,k,l,I_RHOGVX) / prg_pl(g,k,l,I_RHOG)
!MINIAPP          diag_pl(g,k,l,I_vy) = prg_pl(g,k,l,I_RHOGVY) / prg_pl(g,k,l,I_RHOG)
!MINIAPP          diag_pl(g,k,l,I_vz) = prg_pl(g,k,l,I_RHOGVZ) / prg_pl(g,k,l,I_RHOG)
!MINIAPP          ein_pl (g,k,l)      = prg_pl(g,k,l,I_RHOGE ) / prg_pl(g,k,l,I_RHOG)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do iq = 1, TRC_vmax
!MINIAPP       do l  = 1, ADM_lall_pl
!MINIAPP       do k  = 1, ADM_kall
!MINIAPP       do g  = 1, ADM_gall_pl
!MINIAPP          diag_pl(g,k,l,DIAG_vmax0+iq) = prg_pl(g,k,l,PRG_vmax0+iq) / prg_pl(g,k,l,I_RHOG)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       call THRMDYN_tempre( ADM_gall_pl,                  & ! [IN]
!MINIAPP                            ADM_kall,                     & ! [IN]
!MINIAPP                            ADM_lall_pl,                  & ! [IN]
!MINIAPP                            ein_pl (:,:,:),               & ! [IN]
!MINIAPP                            rho_pl (:,:,:),               & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_tem),         & ! [OUT]
!MINIAPP                            diag_pl(:,:,:,I_pre)          ) ! [OUT]
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          !------ interpolation of rhog_h
!MINIAPP          do k = 2, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,k) = ( VMTR_C2Wfact_pl(g,k,1,l) * prg_pl(g,k  ,l,I_RHOG) &
!MINIAPP                              + VMTR_C2Wfact_pl(g,k,2,l) * prg_pl(g,k-1,l,I_RHOG) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,1) = rhog_h_pl(g,2)
!MINIAPP          enddo
!MINIAPP
!MINIAPP          do k = 1, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             diag_pl(g,k,l,I_w) = prg_pl(g,k,l,I_RHOGW) / rhog_h_pl(g,k)
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP    endif
!MINIAPP    !$acc wait(1)
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_prg2diag

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_diag2prg_SP( &
!MINIAPP       prg,  prg_pl, &
!MINIAPP       diag, diag_pl )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_GSGAM2,     &
!MINIAPP       VMTR_GSGAM2_pl,  &
!MINIAPP       VMTR_C2Wfact,    &
!MINIAPP       VMTR_C2Wfact_pl
!MINIAPP    use mod_runconf, only: &
!MINIAPP       PRG_vmax,  &
!MINIAPP       DIAG_vmax, &
!MINIAPP       TRC_vmax
!MINIAPP    use mod_thrmdyn, only: &
!MINIAPP       THRMDYN_rhoein
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(out) :: prg    (ADM_gall   ,ADM_kall,ADM_lall   ,PRG_vmax )
!MINIAPP    real(SP), intent(out) :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
!MINIAPP    real(SP), intent(in)  :: diag   (ADM_gall   ,ADM_kall,ADM_lall   ,DIAG_vmax)
!MINIAPP    real(SP), intent(in)  :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)
!MINIAPP
!MINIAPP    real(SP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP) :: rhog_h   (ADM_gall,   ADM_kall)
!MINIAPP    real(SP) :: rhog_h_pl(ADM_gall_pl,ADM_kall)
!MINIAPP
!MINIAPP    integer  :: g, k, l, iq
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(rho,ein,rhog_h) pcopyin(diag,VMTR_GSGAM2,VMTR_C2Wfact) if(acc_is_present(prg))
!MINIAPP
!MINIAPP    call THRMDYN_rhoein( ADM_gall,                  & ! [IN]
!MINIAPP                         ADM_kall,                  & ! [IN]
!MINIAPP                         ADM_lall,                  & ! [IN]
!MINIAPP                         diag(:,:,:,I_tem),         & ! [IN]
!MINIAPP                         diag(:,:,:,I_pre),         & ! [IN]
!MINIAPP                         diag(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                         rho (:,:,:),               & ! [OUT]
!MINIAPP                         ein (:,:,:)                ) ! [OUT]
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg) pcopyin(rho,ein,diag,VMTR_GSGAM2) if(acc_is_present(prg))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       prg(g,k,l,I_RHOG  ) = rho(g,k,l) * VMTR_GSGAM2(g,k,l)
!MINIAPP       prg(g,k,l,I_RHOGVX) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vx)
!MINIAPP       prg(g,k,l,I_RHOGVY) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vy)
!MINIAPP       prg(g,k,l,I_RHOGVZ) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vz)
!MINIAPP       prg(g,k,l,I_RHOGE ) = prg(g,k,l,I_RHOG) * ein (g,k,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg) pcopyin(diag) if(acc_is_present(prg))
!MINIAPP    do iq = 1, TRC_vmax
!MINIAPP    do l  = 1, ADM_lall
!MINIAPP    do k  = 1, ADM_kall
!MINIAPP    do g  = 1, ADM_gall
!MINIAPP       prg(g,k,l,PRG_vmax0+iq) = prg(g,k,l,I_RHOG) * diag(g,k,l,DIAG_vmax0+iq)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg,rhog_h) pcopyin(diag,VMTR_C2Wfact) if(acc_is_present(prg))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       !------ interpolation of rhog_h
!MINIAPP       do k = 2, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,k) = ( VMTR_C2Wfact(g,k,1,l) * prg(g,k  ,l,I_RHOG) &
!MINIAPP                        + VMTR_C2Wfact(g,k,2,l) * prg(g,k-1,l,I_RHOG) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,1) = rhog_h(g,2)
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          prg(g,k,l,I_RHOGW) = rhog_h(g,k) * diag(g,k,l,I_w)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc end data
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP
!MINIAPP       call THRMDYN_rhoein( ADM_gall_pl,                  & ! [IN]
!MINIAPP                            ADM_kall,                     & ! [IN]
!MINIAPP                            ADM_lall_pl,                  & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_tem),         & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_pre),         & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                            rho_pl (:,:,:),               & ! [OUT]
!MINIAPP                            ein_pl (:,:,:)                ) ! [OUT]
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          prg_pl(g,k,l,I_RHOG  ) = rho_pl(g,k,l) * VMTR_GSGAM2_pl(g,k,l)
!MINIAPP          prg_pl(g,k,l,I_RHOGVX) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vx)
!MINIAPP          prg_pl(g,k,l,I_RHOGVY) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vy)
!MINIAPP          prg_pl(g,k,l,I_RHOGVZ) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vz)
!MINIAPP          prg_pl(g,k,l,I_RHOGE ) = prg_pl(g,k,l,I_RHOG) * ein_pl (g,k,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do iq = 1, TRC_vmax
!MINIAPP       do l  = 1, ADM_lall_pl
!MINIAPP       do k  = 1, ADM_kall
!MINIAPP       do g  = 1, ADM_gall_pl
!MINIAPP          prg_pl(g,k,l,PRG_vmax0+iq) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,DIAG_vmax0+iq)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          !------ interpolation of rhog_h
!MINIAPP          do k = 2, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,k) = ( VMTR_C2Wfact_pl(g,k,1,l) * prg_pl(g,k  ,l,I_RHOG) &
!MINIAPP                              + VMTR_C2Wfact_pl(g,k,2,l) * prg_pl(g,k-1,l,I_RHOG) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,1) = rhog_h_pl(g,2)
!MINIAPP          enddo
!MINIAPP
!MINIAPP          do k = 1, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             prg_pl(g,k,l,I_RHOGW) = rhog_h_pl(g,k) * diag_pl(g,k,l,I_w)
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_diag2prg_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_diag2prg_DP( &
!MINIAPP       prg,  prg_pl, &
!MINIAPP       diag, diag_pl )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_GSGAM2,     &
!MINIAPP       VMTR_GSGAM2_pl,  &
!MINIAPP       VMTR_C2Wfact,    &
!MINIAPP       VMTR_C2Wfact_pl
!MINIAPP    use mod_runconf, only: &
!MINIAPP       PRG_vmax,  &
!MINIAPP       DIAG_vmax, &
!MINIAPP       TRC_vmax
!MINIAPP    use mod_thrmdyn, only: &
!MINIAPP       THRMDYN_rhoein
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(out) :: prg    (ADM_gall   ,ADM_kall,ADM_lall   ,PRG_vmax )
!MINIAPP    real(DP), intent(out) :: prg_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,PRG_vmax )
!MINIAPP    real(DP), intent(in)  :: diag   (ADM_gall   ,ADM_kall,ADM_lall   ,DIAG_vmax)
!MINIAPP    real(DP), intent(in)  :: diag_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,DIAG_vmax)
!MINIAPP
!MINIAPP    real(DP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP) :: ein_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP) :: rhog_h   (ADM_gall,   ADM_kall)
!MINIAPP    real(DP) :: rhog_h_pl(ADM_gall_pl,ADM_kall)
!MINIAPP
!MINIAPP    integer  :: g, k, l, iq
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(rho,ein,rhog_h) pcopyin(diag,VMTR_GSGAM2,VMTR_C2Wfact) if(acc_is_present(prg))
!MINIAPP
!MINIAPP    call THRMDYN_rhoein( ADM_gall,                  & ! [IN]
!MINIAPP                         ADM_kall,                  & ! [IN]
!MINIAPP                         ADM_lall,                  & ! [IN]
!MINIAPP                         diag(:,:,:,I_tem),         & ! [IN]
!MINIAPP                         diag(:,:,:,I_pre),         & ! [IN]
!MINIAPP                         diag(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                         rho (:,:,:),               & ! [OUT]
!MINIAPP                         ein (:,:,:)                ) ! [OUT]
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg) pcopyin(rho,ein,diag,VMTR_GSGAM2) if(acc_is_present(prg))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       prg(g,k,l,I_RHOG  ) = rho(g,k,l) * VMTR_GSGAM2(g,k,l)
!MINIAPP       prg(g,k,l,I_RHOGVX) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vx)
!MINIAPP       prg(g,k,l,I_RHOGVY) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vy)
!MINIAPP       prg(g,k,l,I_RHOGVZ) = prg(g,k,l,I_RHOG) * diag(g,k,l,I_vz)
!MINIAPP       prg(g,k,l,I_RHOGE ) = prg(g,k,l,I_RHOG) * ein (g,k,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg) pcopyin(diag) if(acc_is_present(prg))
!MINIAPP    do iq = 1, TRC_vmax
!MINIAPP    do l  = 1, ADM_lall
!MINIAPP    do k  = 1, ADM_kall
!MINIAPP    do g  = 1, ADM_gall
!MINIAPP       prg(g,k,l,PRG_vmax0+iq) = prg(g,k,l,I_RHOG) * diag(g,k,l,DIAG_vmax0+iq)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(prg,rhog_h) pcopyin(diag,VMTR_C2Wfact) if(acc_is_present(prg))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       !------ interpolation of rhog_h
!MINIAPP       do k = 2, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,k) = ( VMTR_C2Wfact(g,k,1,l) * prg(g,k  ,l,I_RHOG) &
!MINIAPP                        + VMTR_C2Wfact(g,k,2,l) * prg(g,k-1,l,I_RHOG) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          rhog_h(g,1) = rhog_h(g,2)
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall
!MINIAPP          prg(g,k,l,I_RHOGW) = rhog_h(g,k) * diag(g,k,l,I_w)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !$acc end data
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP
!MINIAPP       call THRMDYN_rhoein( ADM_gall_pl,                  & ! [IN]
!MINIAPP                            ADM_kall,                     & ! [IN]
!MINIAPP                            ADM_lall_pl,                  & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_tem),         & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_pre),         & ! [IN]
!MINIAPP                            diag_pl(:,:,:,I_qstr:I_qend), & ! [IN]
!MINIAPP                            rho_pl (:,:,:),               & ! [OUT]
!MINIAPP                            ein_pl (:,:,:)                ) ! [OUT]
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          prg_pl(g,k,l,I_RHOG  ) = rho_pl(g,k,l) * VMTR_GSGAM2_pl(g,k,l)
!MINIAPP          prg_pl(g,k,l,I_RHOGVX) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vx)
!MINIAPP          prg_pl(g,k,l,I_RHOGVY) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vy)
!MINIAPP          prg_pl(g,k,l,I_RHOGVZ) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,I_vz)
!MINIAPP          prg_pl(g,k,l,I_RHOGE ) = prg_pl(g,k,l,I_RHOG) * ein_pl (g,k,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do iq = 1, TRC_vmax
!MINIAPP       do l  = 1, ADM_lall_pl
!MINIAPP       do k  = 1, ADM_kall
!MINIAPP       do g  = 1, ADM_gall_pl
!MINIAPP          prg_pl(g,k,l,PRG_vmax0+iq) = prg_pl(g,k,l,I_RHOG) * diag_pl(g,k,l,DIAG_vmax0+iq)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          !------ interpolation of rhog_h
!MINIAPP          do k = 2, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,k) = ( VMTR_C2Wfact_pl(g,k,1,l) * prg_pl(g,k  ,l,I_RHOG) &
!MINIAPP                              + VMTR_C2Wfact_pl(g,k,2,l) * prg_pl(g,k-1,l,I_RHOG) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhog_h_pl(g,1) = rhog_h_pl(g,2)
!MINIAPP          enddo
!MINIAPP
!MINIAPP          do k = 1, ADM_kall
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             prg_pl(g,k,l,I_RHOGW) = rhog_h_pl(g,k) * diag_pl(g,k,l,I_w)
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP       enddo
!MINIAPP
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_diag2prg_DP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_prg2diag_in_SP(&
       ijdim,       &
       rho,         &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       phi,         &
       gsqgam2,     &
       gsqgam2h,    &
       set_boundary )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_tempre
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(out) :: rho     (ijdim,kdim)
    real(SP), intent(out) :: pre     (ijdim,kdim)
    real(SP), intent(out) :: tem     (ijdim,kdim)
    real(SP), intent(out) :: vx      (ijdim,kdim)
    real(SP), intent(out) :: vy      (ijdim,kdim)
    real(SP), intent(out) :: vz      (ijdim,kdim)
    real(SP), intent(out) :: w       (ijdim,kdim)
    real(SP), intent(out) :: q       (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: rhog    (ijdim,kdim)
    real(SP), intent(in)  :: rhogvx  (ijdim,kdim)
    real(SP), intent(in)  :: rhogvy  (ijdim,kdim)
    real(SP), intent(in)  :: rhogvz  (ijdim,kdim)
    real(SP), intent(in)  :: rhogw   (ijdim,kdim)
    real(SP), intent(in)  :: rhoge   (ijdim,kdim)
    real(SP), intent(in)  :: rhogq   (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: phi     (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2h(ijdim,kdim)
    logical,  intent(in)  :: set_boundary

    real(SP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data pcreate(ein) pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq) &
    !$acc pcopyin(phi,gsqgam2,gsqgam2h,GRD_afact,GRD_bfact) if(acc_is_present(rho))

!OCL XFILL
    !$acc kernels async(1) pcopy(rho,vx,vy,vz,ein) pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhoge,gsqgam2) if(acc_is_present(rho))
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog  (ij,k) / gsqgam2(ij,k)
       vx (ij,k) = rhogvx(ij,k) / rhog(ij,k)
       vy (ij,k) = rhogvy(ij,k) / rhog(ij,k)
       vz (ij,k) = rhogvz(ij,k) / rhog(ij,k)
       ein(ij,k) = rhoge (ij,k) / rhog(ij,k)
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(q) pcopyin(rhog,rhogq) if(acc_is_present(q))
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$acc end kernels

    call THRMDYN_tempre( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         ein(:,:),   & ! [IN]
                         rho(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         tem(:,:),   & ! [OUT]
                         pre(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(w) pcopyin(rhogw,rho,GRD_afact,GRD_bfact,gsqgam2h) if(acc_is_present(w))
    do k  = kmin-1, kmax+1
       if ( k == kmin-1 ) then
          do ij = 1, ijdim
             w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k+1) * rho(ij,k+1) &
                                       + GRD_bfact(k+1) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
          enddo
       else
          do ij = 1, ijdim
             w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k) * rho(ij,k  ) &
                                       + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
          enddo
       endif
    enddo
    !$acc end kernels

    if ( set_boundary ) then
       call BNDCND_thermo( ijdim, & ! [IN]
                           tem,   & ! [INOUT]
                           rho,   & ! [INOUT]
                           pre,   & ! [INOUT]
                           phi    ) ! [IN]

!OCL XFILL
       !$acc kernels async(1) pcopy(vx,vy,vz) if(acc_is_present(vx))
       do ij = 1, ijdim
          vx(ij,kmax+1) = vx(ij,kmax)
          vy(ij,kmax+1) = vy(ij,kmax)
          vz(ij,kmax+1) = vz(ij,kmax)
          vx(ij,kmin-1) = vx(ij,kmin)
          vy(ij,kmin-1) = vy(ij,kmin)
          vz(ij,kmin-1) = vz(ij,kmin)
       enddo
       !$acc end kernels

!OCL XFILL
       !$acc kernels async(1) pcopy(q) if(acc_is_present(q))
       do nq = 1, nqmax
       do ij = 1, ijdim
          q(ij,kmax+1,nq) = 0.0_SP
          q(ij,kmin-1,nq) = 0.0_SP
       enddo
       enddo
       !$acc end kernels
    endif

    !$acc wait(1)
    !$acc end data
    return
  end subroutine CNVVAR_prg2diag_in_SP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_prg2diag_in_DP(&
       ijdim,       &
       rho,         &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       phi,         &
       gsqgam2,     &
       gsqgam2h,    &
       set_boundary )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_tempre
    use mod_bndcnd, only: &
       BNDCND_thermo
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(out) :: rho     (ijdim,kdim)
    real(DP), intent(out) :: pre     (ijdim,kdim)
    real(DP), intent(out) :: tem     (ijdim,kdim)
    real(DP), intent(out) :: vx      (ijdim,kdim)
    real(DP), intent(out) :: vy      (ijdim,kdim)
    real(DP), intent(out) :: vz      (ijdim,kdim)
    real(DP), intent(out) :: w       (ijdim,kdim)
    real(DP), intent(out) :: q       (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: rhog    (ijdim,kdim)
    real(DP), intent(in)  :: rhogvx  (ijdim,kdim)
    real(DP), intent(in)  :: rhogvy  (ijdim,kdim)
    real(DP), intent(in)  :: rhogvz  (ijdim,kdim)
    real(DP), intent(in)  :: rhogw   (ijdim,kdim)
    real(DP), intent(in)  :: rhoge   (ijdim,kdim)
    real(DP), intent(in)  :: rhogq   (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: phi     (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2h(ijdim,kdim)
    logical,  intent(in)  :: set_boundary

    real(DP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data pcreate(ein) pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq) &
    !$acc pcopyin(phi,gsqgam2,gsqgam2h,GRD_afact,GRD_bfact) if(acc_is_present(rho))

!OCL XFILL
    !$acc kernels async(1) pcopy(rho,vx,vy,vz,ein) pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhoge,gsqgam2) if(acc_is_present(rho))
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog  (ij,k) / gsqgam2(ij,k)
       vx (ij,k) = rhogvx(ij,k) / rhog(ij,k)
       vy (ij,k) = rhogvy(ij,k) / rhog(ij,k)
       vz (ij,k) = rhogvz(ij,k) / rhog(ij,k)
       ein(ij,k) = rhoge (ij,k) / rhog(ij,k)
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(q) pcopyin(rhog,rhogq) if(acc_is_present(q))
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$acc end kernels

    call THRMDYN_tempre( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         ein(:,:),   & ! [IN]
                         rho(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         tem(:,:),   & ! [OUT]
                         pre(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(w) pcopyin(rhogw,rho,GRD_afact,GRD_bfact,gsqgam2h) if(acc_is_present(w))
    do k  = kmin-1, kmax+1
       if ( k == kmin-1 ) then
          do ij = 1, ijdim
             w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k+1) * rho(ij,k+1) &
                                       + GRD_bfact(k+1) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
          enddo
       else
          do ij = 1, ijdim
             w(ij,k) = rhogw(ij,k) / ( ( GRD_afact(k) * rho(ij,k  ) &
                                       + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
          enddo
       endif
    enddo
    !$acc end kernels

    if ( set_boundary ) then
       call BNDCND_thermo( ijdim, & ! [IN]
                           tem,   & ! [INOUT]
                           rho,   & ! [INOUT]
                           pre,   & ! [INOUT]
                           phi    ) ! [IN]

!OCL XFILL
       !$acc kernels async(1) pcopy(vx,vy,vz) if(acc_is_present(vx))
       do ij = 1, ijdim
          vx(ij,kmax+1) = vx(ij,kmax)
          vy(ij,kmax+1) = vy(ij,kmax)
          vz(ij,kmax+1) = vz(ij,kmax)
          vx(ij,kmin-1) = vx(ij,kmin)
          vy(ij,kmin-1) = vy(ij,kmin)
          vz(ij,kmin-1) = vz(ij,kmin)
       enddo
       !$acc end kernels

!OCL XFILL
       !$acc kernels async(1) pcopy(q) if(acc_is_present(q))
       do nq = 1, nqmax
       do ij = 1, ijdim
          q(ij,kmax+1,nq) = 0.0_DP
          q(ij,kmin-1,nq) = 0.0_DP
       enddo
       enddo
       !$acc end kernels
    endif

    !$acc wait(1)
    !$acc end data
    return
  end subroutine CNVVAR_prg2diag_in_DP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_diag2prg_in_SP(&
       ijdim,       &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       gsqgam2,     &
       gsqgam2h     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_rhoein
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(in)  :: pre     (ijdim,kdim)
    real(SP), intent(in)  :: tem     (ijdim,kdim)
    real(SP), intent(in)  :: vx      (ijdim,kdim)
    real(SP), intent(in)  :: vy      (ijdim,kdim)
    real(SP), intent(in)  :: vz      (ijdim,kdim)
    real(SP), intent(in)  :: w       (ijdim,kdim)
    real(SP), intent(in)  :: q       (ijdim,kdim,nqmax)
    real(SP), intent(out) :: rhog    (ijdim,kdim)
    real(SP), intent(out) :: rhogvx  (ijdim,kdim)
    real(SP), intent(out) :: rhogvy  (ijdim,kdim)
    real(SP), intent(out) :: rhogvz  (ijdim,kdim)
    real(SP), intent(out) :: rhogw   (ijdim,kdim)
    real(SP), intent(out) :: rhoge   (ijdim,kdim)
    real(SP), intent(out) :: rhogq   (ijdim,kdim,nqmax)
    real(SP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(SP), intent(in)  :: gsqgam2h(ijdim,kdim)

    real(SP) :: rho(ijdim,kdim)
    real(SP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data pcreate(rho,ein) pcopyin(pre,tem,vx,vy,vz,w,q) &
    !$acc pcopyin(gsqgam2,gsqgam2h,GRD_afact,GRD_bfact) if(acc_is_present(rhog))

    call THRMDYN_rhoein( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         tem(:,:),   & ! [IN]
                         pre(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         rho(:,:),   & ! [OUT]
                         ein(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(rhog,rhogvx,rhogvy,rhogvz,rhoge) pcopyin(rho,vx,vy,vz,ein,gsqgam2) if(acc_is_present(rhog))
    do k  = 1, kdim
    do ij = 1, ijdim
       rhog  (ij,k) = rho(ij,k) * gsqgam2(ij,k)
       rhogvx(ij,k) = vx (ij,k) * rhog(ij,k)
       rhogvy(ij,k) = vy (ij,k) * rhog(ij,k)
       rhogvz(ij,k) = vz (ij,k) * rhog(ij,k)
       rhoge (ij,k) = ein(ij,k) * rhog(ij,k)
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogq) pcopyin(q,rhog) if(acc_is_present(rhogq))
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,nq) = q(ij,k,nq) * rhog(ij,k)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1) pcopy(rhogw) pcopyin(w,rho,GRD_afact,GRD_bfact,gsqgam2h) if(acc_is_present(rhogw))
    do k  = kmin-1, kmax+1
       if ( k == kmin-1 ) then
          do ij = 1, ijdim
             rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(k+1) * rho(ij,k+1) &
                                       + GRD_bfact(k+1) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
          enddo
       else
          do ij = 1, ijdim
             rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(k) * rho(ij,k  ) &
                                       + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
          enddo
       endif
    enddo
    !$acc end kernels

    !$acc wait(1)
    !$acc end data
    return
  end subroutine CNVVAR_diag2prg_in_SP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_diag2prg_in_DP(&
       ijdim,       &
       pre,         &
       tem,         &
       vx,          &
       vy,          &
       vz,          &
       w,           &
       q,           &
       rhog,        &
       rhogvx,      &
       rhogvy,      &
       rhogvz,      &
       rhogw,       &
       rhoge,       &
       rhogq,       &
       gsqgam2,     &
       gsqgam2h     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use mod_grd, only: &
       GRD_afact, &
       GRD_bfact
    use mod_runconf, only: &
       nqmax => TRC_VMAX
    use mod_thrmdyn, only: &
       THRMDYN_rhoein
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(in)  :: pre     (ijdim,kdim)
    real(DP), intent(in)  :: tem     (ijdim,kdim)
    real(DP), intent(in)  :: vx      (ijdim,kdim)
    real(DP), intent(in)  :: vy      (ijdim,kdim)
    real(DP), intent(in)  :: vz      (ijdim,kdim)
    real(DP), intent(in)  :: w       (ijdim,kdim)
    real(DP), intent(in)  :: q       (ijdim,kdim,nqmax)
    real(DP), intent(out) :: rhog    (ijdim,kdim)
    real(DP), intent(out) :: rhogvx  (ijdim,kdim)
    real(DP), intent(out) :: rhogvy  (ijdim,kdim)
    real(DP), intent(out) :: rhogvz  (ijdim,kdim)
    real(DP), intent(out) :: rhogw   (ijdim,kdim)
    real(DP), intent(out) :: rhoge   (ijdim,kdim)
    real(DP), intent(out) :: rhogq   (ijdim,kdim,nqmax)
    real(DP), intent(in)  :: gsqgam2 (ijdim,kdim)
    real(DP), intent(in)  :: gsqgam2h(ijdim,kdim)

    real(DP) :: rho(ijdim,kdim)
    real(DP) :: ein(ijdim,kdim)

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------
    !$acc data pcreate(rho,ein) pcopyin(pre,tem,vx,vy,vz,w,q) &
    !$acc pcopyin(gsqgam2,gsqgam2h,GRD_afact,GRD_bfact) if(acc_is_present(rhog))

    call THRMDYN_rhoein( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         tem(:,:),   & ! [IN]
                         pre(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         rho(:,:),   & ! [OUT]
                         ein(:,:)    ) ! [OUT]

    !$acc kernels async(1) pcopy(rhog,rhogvx,rhogvy,rhogvz,rhoge) pcopyin(rho,vx,vy,vz,ein,gsqgam2) if(acc_is_present(rhog))
    do k  = 1, kdim
    do ij = 1, ijdim
       rhog  (ij,k) = rho(ij,k) * gsqgam2(ij,k)
       rhogvx(ij,k) = vx (ij,k) * rhog(ij,k)
       rhogvy(ij,k) = vy (ij,k) * rhog(ij,k)
       rhogvz(ij,k) = vz (ij,k) * rhog(ij,k)
       rhoge (ij,k) = ein(ij,k) * rhog(ij,k)
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogq) pcopyin(q,rhog) if(acc_is_present(rhogq))
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,nq) = q(ij,k,nq) * rhog(ij,k)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels async(1) pcopy(rhogw) pcopyin(w,rho,GRD_afact,GRD_bfact,gsqgam2h) if(acc_is_present(rhogw))
    k = kmin-1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(kmin) * rho(ij,k+1) &
                                 + GRD_bfact(kmin) * rho(ij,k  ) ) * gsqgam2h(ij,k+1) )
    enddo

    do k  = kmin, kmax+1
    do ij = 1, ijdim
       rhogw(ij,k) = w(ij,k) * ( ( GRD_afact(k) * rho(ij,k  ) &
                                 + GRD_bfact(k) * rho(ij,k-1) ) * gsqgam2h(ij,k) )
    enddo
    enddo
    !$acc end kernels

    !$acc wait(1)
    !$acc end data
    return
  end subroutine CNVVAR_diag2prg_in_DP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_rhogkin_SP( &
!MINIAPP       rhog,    rhog_pl,   &
!MINIAPP       rhogvx,  rhogvx_pl, &
!MINIAPP       rhogvy,  rhogvy_pl, &
!MINIAPP       rhogvz,  rhogvz_pl, &
!MINIAPP       rhogw,   rhogw_pl,  &
!MINIAPP       rhogkin, rhogkin_pl,&
!MINIAPP       rhogkin_h, rhogkin_v ) ! working arrays for OenACC
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_kall,    &
!MINIAPP       ADM_kmax,    &
!MINIAPP       ADM_kmin
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_C2Wfact,    &
!MINIAPP       VMTR_C2Wfact_pl, &
!MINIAPP       VMTR_W2Cfact,    &
!MINIAPP       VMTR_W2Cfact_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(in)  :: rhog      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 )
!MINIAPP    real(SP), intent(in)  :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: rhogvx    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vx
!MINIAPP    real(SP), intent(in)  :: rhogvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: rhogvy    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vy
!MINIAPP    real(SP), intent(in)  :: rhogvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: rhogvz    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vz
!MINIAPP    real(SP), intent(in)  :: rhogvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: rhogw     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X w
!MINIAPP    real(SP), intent(in)  :: rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: rhogkin   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X kin
!MINIAPP    real(SP), intent(out) :: rhogkin_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    real(SP) :: rhogkin_h   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
!MINIAPP    real(SP) :: rhogkin_h_pl(ADM_gall_pl,ADM_kall)
!MINIAPP    real(SP) :: rhogkin_v   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (vertical)
!MINIAPP    real(SP) :: rhogkin_v_pl(ADM_gall_pl,ADM_kall)
!MINIAPP
!MINIAPP    integer  :: gall, kmin, kmax, lall
!MINIAPP
!MINIAPP    integer  :: g, k, l
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    call PROF_rapstart('CNV_rhogkin',3)
!MINIAPP
!MINIAPP    !!$acc data pcreate(rhogkin_h,rhogkin_v) &
!MINIAPP    !!$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_W2Cfact) if(acc_is_present(rhogkin))
!MINIAPP
!MINIAPP    gall = ADM_gall
!MINIAPP    kmin = ADM_kmin
!MINIAPP    kmax = ADM_kmax
!MINIAPP    lall = ADM_lall
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(rhogkin) present(rhogkin_h,rhogkin_v) &
!MINIAPP    !$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_W2Cfact) if(acc_is_present(rhogkin))
!MINIAPP    !$omp parallel default(none),private(g,k,l), &
!MINIAPP    !$omp shared(gall,kmin,kmax,lall,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
!MINIAPP    !$omp        rhogkin_h,rhogkin_v,VMTR_C2Wfact,VMTR_W2Cfact)
!MINIAPP    do l = 1, lall
!MINIAPP       !--- horizontal kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_h(g,k) = 0.5_SP * ( rhogvx(g,k,l) * rhogvx(g,k,l) &
!MINIAPP                                    + rhogvy(g,k,l) * rhogvy(g,k,l) &
!MINIAPP                                    + rhogvz(g,k,l) * rhogvz(g,k,l) ) / rhog(g,k,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP       !--- vertical kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin+1, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_v(g,k) = 0.5_SP * ( rhogw(g,k,l) * rhogw(g,k,l) ) &
!MINIAPP                         / ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
!MINIAPP                           + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_v(g,kmin  ) = 0.0_SP
!MINIAPP          rhogkin_v(g,kmax+1) = 0.0_SP
!MINIAPP       enddo
!MINIAPP       !$omp end do
!MINIAPP
!MINIAPP       !--- total kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin(g,k,l) = rhogkin_h(g,k)                             & ! horizontal
!MINIAPP                         + ( VMTR_W2Cfact(g,k,1,l) * rhogkin_v(g,k+1) & ! vertical
!MINIAPP                           + VMTR_W2Cfact(g,k,2,l) * rhogkin_v(g,k  ) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin(g,kmin-1,l) = 0.0_SP
!MINIAPP          rhogkin(g,kmax+1,l) = 0.0_SP
!MINIAPP       enddo
!MINIAPP       !$omp end do
!MINIAPP    enddo
!MINIAPP    !$omp end parallel
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !!$acc end data
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          !--- horizontal kinetic energy
!MINIAPP          do k = ADM_kmin, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_h_pl(g,k) = 0.5_SP * ( rhogvx_pl(g,k,l) * rhogvx_pl(g,k,l) &
!MINIAPP                                          + rhogvy_pl(g,k,l) * rhogvy_pl(g,k,l) &
!MINIAPP                                          + rhogvz_pl(g,k,l) * rhogvz_pl(g,k,l) ) / rhog_pl(g,k,l)
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP
!MINIAPP          !--- vertical kinetic energy
!MINIAPP          do k = ADM_kmin+1, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_v_pl(g,k) = 0.5_SP * ( rhogw_pl(g,k,l) * rhogw_pl(g,k,l) ) &
!MINIAPP                               / ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
!MINIAPP                                 + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          rhogkin_v_pl(:,ADM_kmin  ) = 0.0_SP
!MINIAPP          rhogkin_v_pl(:,ADM_kmax+1) = 0.0_SP
!MINIAPP
!MINIAPP          !--- total kinetic energy
!MINIAPP          do k = ADM_kmin, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_pl(g,k,l) = rhogkin_h_pl(g,k)                                & ! horizontal
!MINIAPP                               + ( VMTR_W2Cfact_pl(g,k,1,l) * rhogkin_v_pl(g,k+1) & ! vertical
!MINIAPP                                 + VMTR_W2Cfact_pl(g,k,2,l) * rhogkin_v_pl(g,k  ) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          rhogkin_pl(:,ADM_kmin-1,l) = 0.0_SP
!MINIAPP          rhogkin_pl(:,ADM_kmax+1,l) = 0.0_SP
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP    if ( acc_wait_on ) then
!MINIAPP      !$acc wait(1)
!MINIAPP    endif
!MINIAPP    call PROF_rapend  ('CNV_rhogkin',3)
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_rhogkin_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_rhogkin_DP( &
!MINIAPP       rhog,    rhog_pl,   &
!MINIAPP       rhogvx,  rhogvx_pl, &
!MINIAPP       rhogvy,  rhogvy_pl, &
!MINIAPP       rhogvz,  rhogvz_pl, &
!MINIAPP       rhogw,   rhogw_pl,  &
!MINIAPP       rhogkin, rhogkin_pl, &
!MINIAPP       rhogkin_h, rhogkin_v ) ! working arrays for OenACC
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_kall,    &
!MINIAPP       ADM_kmax,    &
!MINIAPP       ADM_kmin
!MINIAPP    use mod_vmtr, only: &
!MINIAPP       VMTR_C2Wfact,    &
!MINIAPP       VMTR_C2Wfact_pl, &
!MINIAPP       VMTR_W2Cfact,    &
!MINIAPP       VMTR_W2Cfact_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(in)  :: rhog      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 )
!MINIAPP    real(DP), intent(in)  :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: rhogvx    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vx
!MINIAPP    real(DP), intent(in)  :: rhogvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: rhogvy    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vy
!MINIAPP    real(DP), intent(in)  :: rhogvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: rhogvz    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X vz
!MINIAPP    real(DP), intent(in)  :: rhogvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: rhogw     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X w
!MINIAPP    real(DP), intent(in)  :: rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: rhogkin   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho X ( G^1/2 X gamma2 ) X kin
!MINIAPP    real(DP), intent(out) :: rhogkin_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    real(DP) :: rhogkin_h   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
!MINIAPP    real(DP) :: rhogkin_h_pl(ADM_gall_pl,ADM_kall)
!MINIAPP    real(DP) :: rhogkin_v   (ADM_gall,   ADM_kall) ! rho X ( G^1/2 X gamma2 ) X kin (vertical)
!MINIAPP    real(DP) :: rhogkin_v_pl(ADM_gall_pl,ADM_kall)
!MINIAPP
!MINIAPP    integer  :: gall, kmin, kmax, lall
!MINIAPP
!MINIAPP    integer  :: g, k, l
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    call PROF_rapstart('CNV_rhogkin',3)
!MINIAPP
!MINIAPP    !!$acc data pcreate(rhogkin_h,rhogkin_v) &
!MINIAPP    !!$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_W2Cfact) if(acc_is_present(rhogkin))
!MINIAPP
!MINIAPP    gall = ADM_gall
!MINIAPP    kmin = ADM_kmin
!MINIAPP    kmax = ADM_kmax
!MINIAPP    lall = ADM_lall
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(rhogkin,rhogkin_h,rhogkin_v) &
!MINIAPP    !$acc pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,VMTR_C2Wfact,VMTR_W2Cfact) if(acc_is_present(rhogkin))
!MINIAPP    !$omp parallel default(none),private(g,k,l), &
!MINIAPP    !$omp shared(gall,kmin,kmax,lall,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
!MINIAPP    !$omp        rhogkin_h,rhogkin_v,VMTR_C2Wfact,VMTR_W2Cfact)
!MINIAPP    do l = 1, lall
!MINIAPP       !--- horizontal kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_h(g,k) = 0.5_DP * ( rhogvx(g,k,l) * rhogvx(g,k,l) &
!MINIAPP                                    + rhogvy(g,k,l) * rhogvy(g,k,l) &
!MINIAPP                                    + rhogvz(g,k,l) * rhogvz(g,k,l) ) / rhog(g,k,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP       !--- vertical kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin+1, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_v(g,k) = 0.5_DP * ( rhogw(g,k,l) * rhogw(g,k,l) ) &
!MINIAPP                         / ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
!MINIAPP                           + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin_v(g,kmin  ) = 0.0_DP
!MINIAPP          rhogkin_v(g,kmax+1) = 0.0_DP
!MINIAPP       enddo
!MINIAPP       !$omp end do
!MINIAPP
!MINIAPP       !--- total kinetic energy
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do k = kmin, kmax
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin(g,k,l) = rhogkin_h(g,k)                             & ! horizontal
!MINIAPP                         + ( VMTR_W2Cfact(g,k,1,l) * rhogkin_v(g,k+1) & ! vertical
!MINIAPP                           + VMTR_W2Cfact(g,k,2,l) * rhogkin_v(g,k  ) )
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       !$omp end do nowait
!MINIAPP
!MINIAPP!OCL XFILL
!MINIAPP       !$omp do
!MINIAPP       do g = 1, gall
!MINIAPP          rhogkin(g,kmin-1,l) = 0.0_DP
!MINIAPP          rhogkin(g,kmax+1,l) = 0.0_DP
!MINIAPP       enddo
!MINIAPP       !$omp end do
!MINIAPP    enddo
!MINIAPP    !$omp end parallel
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    !!$acc end data
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          !--- horizontal kinetic energy
!MINIAPP          do k = ADM_kmin, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_h_pl(g,k) = 0.5_DP * ( rhogvx_pl(g,k,l) * rhogvx_pl(g,k,l) &
!MINIAPP                                          + rhogvy_pl(g,k,l) * rhogvy_pl(g,k,l) &
!MINIAPP                                          + rhogvz_pl(g,k,l) * rhogvz_pl(g,k,l) ) / rhog_pl(g,k,l)
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP
!MINIAPP          !--- vertical kinetic energy
!MINIAPP          do k = ADM_kmin+1, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_v_pl(g,k) = 0.5_DP * ( rhogw_pl(g,k,l) * rhogw_pl(g,k,l) ) &
!MINIAPP                               / ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
!MINIAPP                                 + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          rhogkin_v_pl(:,ADM_kmin  ) = 0.0_DP
!MINIAPP          rhogkin_v_pl(:,ADM_kmax+1) = 0.0_DP
!MINIAPP
!MINIAPP          !--- total kinetic energy
!MINIAPP          do k = ADM_kmin, ADM_kmax
!MINIAPP          do g = 1, ADM_gall_pl
!MINIAPP             rhogkin_pl(g,k,l) = rhogkin_h_pl(g,k)                                & ! horizontal
!MINIAPP                               + ( VMTR_W2Cfact_pl(g,k,1,l) * rhogkin_v_pl(g,k+1) & ! vertical
!MINIAPP                                 + VMTR_W2Cfact_pl(g,k,2,l) * rhogkin_v_pl(g,k  ) )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          rhogkin_pl(:,ADM_kmin-1,l) = 0.0_DP
!MINIAPP          rhogkin_pl(:,ADM_kmax+1,l) = 0.0_DP
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP    if ( acc_wait_on ) then
!MINIAPP       !$acc wait(1)
!MINIAPP    endif
!MINIAPP    call PROF_rapend  ('CNV_rhogkin',3)
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_rhogkin_DP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_rhogkin_in_SP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact,   &
       W2Cfact,   &
       rhogkin,   &
       rhogkin_h, &
       rhogkin_v  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(SP), intent(in)  :: rhog     (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 )
    real(SP), intent(in)  :: rhogvx   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vx
    real(SP), intent(in)  :: rhogvy   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vy
    real(SP), intent(in)  :: rhogvz   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vz
    real(SP), intent(in)  :: rhogw    (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(in)  :: C2Wfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(in)  :: W2Cfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(SP), intent(out) :: rhogkin  (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin
    real(SP), intent(out) :: rhogkin_h(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(SP), intent(out) :: rhogkin_v(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (vertical)

    integer  :: gall, kmin, kmax

    integer  :: g, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,C2Wfact,W2Cfact) if(acc_is_present(rhogkin))

    call PROF_rapstart('CNV_rhogkin',3)

    gall = ijdim
    kmin = 2
    kmax = kdim-1

    !$omp parallel default(none),private(g,k),                           &
    !$omp shared(gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,C2Wfact,W2Cfact)

    !--- horizontal kinetic energy
!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_h) pcopyin(rhog,rhogvx,rhogvy,rhogvz) if(acc_is_present(rhogkin_h))
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin_h(g,k) = 0.5_SP * ( rhogvx(g,k) * rhogvx(g,k) &
                                 + rhogvy(g,k) * rhogvy(g,k) &
                                 + rhogvz(g,k) * rhogvz(g,k) ) / rhog(g,k)
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_h) if(acc_is_present(rhogkin_h))
    !$omp do
    do g = 1, gall
       rhogkin_h(g,kmin-1) = 0.0_SP
       rhogkin_h(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do nowait
    !$acc end kernels

    !--- vertical kinetic energy
    !$acc kernels async(1) pcopy(rhogkin_v) pcopyin(rhog,rhogw,C2Wfact) if(acc_is_present(rhogkin_v))
    !$omp do
    do k = kmin+1, kmax
    do g = 1, gall
       rhogkin_v(g,k) = 0.5_SP * ( rhogw(g,k) * rhogw(g,k) ) &
                      / ( C2Wfact(g,k,1) * rhog(g,k  ) &
                        + C2Wfact(g,k,2) * rhog(g,k-1) )
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_v) if(acc_is_present(rhogkin_v))
    !$omp do
    do g = 1, gall
       rhogkin_v(g,kmin-1) = 0.0_SP
       rhogkin_v(g,kmin  ) = 0.0_SP
       rhogkin_v(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do
    !$acc end kernels

    !--- total kinetic energy
!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin) pcopyin(rhogkin_h,rhogkin_v,W2Cfact) if(acc_is_present(rhogkin))
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin(g,k) = rhogkin_h(g,k)                      & ! horizontal
                    + ( W2Cfact(g,k,1) * rhogkin_v(g,k+1) & ! vertical
                      + W2Cfact(g,k,2) * rhogkin_v(g,k  ) )
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin) if(acc_is_present(rhogkin))
    !$omp do
    do g = 1, gall
       rhogkin(g,kmin-1) = 0.0_SP
       rhogkin(g,kmax+1) = 0.0_SP
    enddo
    !$omp end do
    !$acc end kernels

    !$omp end parallel

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif
    call PROF_rapend  ('CNV_rhogkin',3)

    !!$acc end data
    return
  end subroutine CNVVAR_rhogkin_in_SP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_rhogkin_in_DP( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact,   &
       W2Cfact,   &
       rhogkin,   &
       rhogkin_h, &
       rhogkin_v  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(DP), intent(in)  :: rhog     (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 )
    real(DP), intent(in)  :: rhogvx   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vx
    real(DP), intent(in)  :: rhogvy   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vy
    real(DP), intent(in)  :: rhogvz   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vz
    real(DP), intent(in)  :: rhogw    (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(in)  :: C2Wfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(in)  :: W2Cfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(DP), intent(out) :: rhogkin  (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin
    real(DP), intent(out) :: rhogkin_h(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(DP), intent(out) :: rhogkin_v(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (vertical)

    integer  :: gall, kmin, kmax

    integer  :: g, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(rhog,rhogvx,rhogvy,rhogvz,rhogw,C2Wfact,W2Cfact) if(acc_is_present(rhogkin))

    call PROF_rapstart('CNV_rhogkin',3)

    gall = ijdim
    kmin = 2
    kmax = kdim-1

    !$omp parallel default(none),private(g,k),                           &
    !$omp shared(gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,C2Wfact,W2Cfact)

    !--- horizontal kinetic energy
!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_h) pcopyin(rhog,rhogvx,rhogvy,rhogvz) if(acc_is_present(rhogkin_h))
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin_h(g,k) = 0.5_DP * ( rhogvx(g,k) * rhogvx(g,k) &
                                 + rhogvy(g,k) * rhogvy(g,k) &
                                 + rhogvz(g,k) * rhogvz(g,k) ) / rhog(g,k)
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_h) if(acc_is_present(rhogkin_h))
    !$omp do
    do g = 1, gall
       rhogkin_h(g,kmin-1) = 0.0_DP
       rhogkin_h(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do nowait
    !$acc end kernels

    !--- vertical kinetic energy
    !$acc kernels async(1) pcopy(rhogkin_v) pcopyin(rhog,rhogw,C2Wfact) if(acc_is_present(rhogkin_v))
    !$omp do
    do k = kmin+1, kmax
    do g = 1, gall
       rhogkin_v(g,k) = 0.5_DP * ( rhogw(g,k) * rhogw(g,k) ) &
                      / ( C2Wfact(g,k,1) * rhog(g,k  ) &
                        + C2Wfact(g,k,2) * rhog(g,k-1) )
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin_v) if(acc_is_present(rhogkin_v))
    !$omp do
    do g = 1, gall
       rhogkin_v(g,kmin-1) = 0.0_DP
       rhogkin_v(g,kmin  ) = 0.0_DP
       rhogkin_v(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do
    !$acc end kernels

    !--- total kinetic energy
!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin) pcopyin(rhogkin_h,rhogkin_v,W2Cfact) if(acc_is_present(rhogkin))
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin(g,k) = rhogkin_h(g,k)                      & ! horizontal
                    + ( W2Cfact(g,k,1) * rhogkin_v(g,k+1) & ! vertical
                      + W2Cfact(g,k,2) * rhogkin_v(g,k  ) )
    enddo
    enddo
    !$omp end do nowait
    !$acc end kernels

!OCL XFILL
    !$acc kernels async(1) pcopy(rhogkin) if(acc_is_present(rhogkin))
    !$omp do
    do g = 1, gall
       rhogkin(g,kmin-1) = 0.0_DP
       rhogkin(g,kmax+1) = 0.0_DP
    enddo
    !$omp end do
    !$acc end kernels

    !$omp end parallel

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif
    call PROF_rapend  ('CNV_rhogkin',3)

    !!$acc end data
    return
  end subroutine CNVVAR_rhogkin_in_DP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_uv2vh_SP( &
!MINIAPP       ucos, ucos_pl, &
!MINIAPP       vcos, vcos_pl, &
!MINIAPP       vx,   vx_pl,   &
!MINIAPP       vy,   vy_pl,   &
!MINIAPP       vz,   vz_pl,   &
!MINIAPP       withcos        )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(in)  :: ucos   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vcos   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, k, l, k0
!MINIAPP
!MINIAPP    real(SP) :: sw
!MINIAPP    real(SP) :: u
!MINIAPP    real(SP) :: v
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(ucos,vcos,GRD_LAT,GMTR_p) if(acc_is_present(vx))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_SP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_SP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(vx))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(ucos,vcos,coslat,GMTR_p) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       sw = 0.5_SP + sign(0.5_SP,real(-abs(coslat(g,l)),kind=SP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP       u = ucos(g,k,l) * ( 1.0_SP - sw ) / ( coslat(g,l) - sw )
!MINIAPP       v = vcos(g,k,l) * ( 1.0_SP - sw ) / ( coslat(g,l) - sw )
!MINIAPP
!MINIAPP       vx(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JX)
!MINIAPP       vy(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JY)
!MINIAPP       vz(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IZ) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JZ)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          sw = 0.5_SP + sign(0.5_SP,real(-abs(coslat_pl(g,l)),kind=SP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP          u = ucos_pl(g,k,l) * ( 1.0_SP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP          v = vcos_pl(g,k,l) * ( 1.0_SP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP
!MINIAPP          vx_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JX)
!MINIAPP          vy_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JY)
!MINIAPP          vz_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IZ) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JZ)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_uv2vh_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_uv2vh_DP( &
!MINIAPP       ucos, ucos_pl, &
!MINIAPP       vcos, vcos_pl, &
!MINIAPP       vx,   vx_pl,   &
!MINIAPP       vy,   vy_pl,   &
!MINIAPP       vz,   vz_pl,   &
!MINIAPP       withcos        )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(in)  :: ucos   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vcos   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vx     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vy     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vz     (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, k, l, k0
!MINIAPP
!MINIAPP    real(DP) :: sw
!MINIAPP    real(DP) :: u
!MINIAPP    real(DP) :: v
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(ucos,vcos,GRD_LAT,GMTR_p) if(acc_is_present(vx))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_DP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_DP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(vx))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(ucos,vcos,coslat,GMTR_p) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       sw = 0.5_DP + sign(0.5_DP,real(-abs(coslat(g,l)),kind=DP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP       u = ucos(g,k,l) * ( 1.0_DP - sw ) / ( coslat(g,l) - sw )
!MINIAPP       v = vcos(g,k,l) * ( 1.0_DP - sw ) / ( coslat(g,l) - sw )
!MINIAPP
!MINIAPP       vx(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JX)
!MINIAPP       vy(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JY)
!MINIAPP       vz(g,k,l) = u * GMTR_p(g,k0,l,GMTR_p_IZ) &
!MINIAPP                 + v * GMTR_p(g,k0,l,GMTR_p_JZ)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          sw = 0.5_DP + sign(0.5_DP,real(-abs(coslat_pl(g,l)),kind=DP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP          u = ucos_pl(g,k,l) * ( 1.0_DP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP          v = vcos_pl(g,k,l) * ( 1.0_DP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP
!MINIAPP          vx_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JX)
!MINIAPP          vy_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JY)
!MINIAPP          vz_pl(g,k,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IZ) &
!MINIAPP                       + v * GMTR_p_pl(g,k0,l,GMTR_p_JZ)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_uv2vh_DP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_vh2uv_SP( &
!MINIAPP       u,  u_pl,  &
!MINIAPP       v,  v_pl,  &
!MINIAPP       vx, vx_pl, &
!MINIAPP       vy, vy_pl, &
!MINIAPP       vz, vz_pl, &
!MINIAPP       withcos    )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(out) :: u    (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: v    (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, k, l, k0
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(vx,vy,vz,GRD_LAT,GMTR_p) if(acc_is_present(u))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_SP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_SP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(u))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,coslat,GMTR_p) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       u(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
!MINIAPP       v(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
!MINIAPP                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
!MINIAPP                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          u_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
!MINIAPP          v_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
!MINIAPP                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
!MINIAPP                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_vh2uv_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_vh2uv_DP( &
!MINIAPP       u,  u_pl,  &
!MINIAPP       v,  v_pl,  &
!MINIAPP       vx, vx_pl, &
!MINIAPP       vy, vy_pl, &
!MINIAPP       vz, vz_pl, &
!MINIAPP       withcos    )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl, &
!MINIAPP       ADM_kall
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(out) :: u    (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: v    (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, k, l, k0
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(vx,vy,vz,GRD_LAT,GMTR_p) if(acc_is_present(u))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_DP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_DP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(u))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,coslat,GMTR_p) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do k = 1, ADM_kall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       u(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
!MINIAPP       v(g,k,l) = ( vx(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
!MINIAPP                  + vy(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
!MINIAPP                  + vz(g,k,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do k = 1, ADM_kall
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          u_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
!MINIAPP          v_pl(g,k,l) = ( vx_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
!MINIAPP                        + vy_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
!MINIAPP                        + vz_pl(g,k,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_vh2uv_DP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_vh2uv_in_SP( &
       ijdim, &
       u,     &
       v,     &
       vx,    &
       vy,    &
       vz,    &
       ix,    &
       iy,    &
       iz,    &
       jx,    &
       jy,    &
       jz     )
    use mod_adm, only: &
       kdim => ADM_kall
    implicit none

    integer,  intent(in)  :: ijdim
    real(SP), intent(out) :: u (ijdim,kdim)
    real(SP), intent(out) :: v (ijdim,kdim)
    real(SP), intent(in)  :: vx(ijdim,kdim)
    real(SP), intent(in)  :: vy(ijdim,kdim)
    real(SP), intent(in)  :: vz(ijdim,kdim)
    real(SP), intent(in)  :: ix(ijdim)
    real(SP), intent(in)  :: iy(ijdim)
    real(SP), intent(in)  :: iz(ijdim)
    real(SP), intent(in)  :: jx(ijdim)
    real(SP), intent(in)  :: jy(ijdim)
    real(SP), intent(in)  :: jz(ijdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(vx,vy,vz,ix,iy,iz) if(acc_is_present(u))

!OCL XFILL
    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,ix,iy,iz) if(acc_is_present(u))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,u,v,vx,vy,vz,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       u(ij,k) = vx(ij,k) * ix(ij) &
               + vy(ij,k) * iy(ij) &
               + vz(ij,k) * iz(ij)
       v(ij,k) = vx(ij,k) * jx(ij) &
               + vy(ij,k) * jy(ij) &
               + vz(ij,k) * jz(ij)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif
    !!$acc end data
    return
  end subroutine CNVVAR_vh2uv_in_SP

  !-----------------------------------------------------------------------------
  subroutine CNVVAR_vh2uv_in_DP( &
       ijdim, &
       u,     &
       v,     &
       vx,    &
       vy,    &
       vz,    &
       ix,    &
       iy,    &
       iz,    &
       jx,    &
       jy,    &
       jz     )
    use mod_adm, only: &
       kdim => ADM_kall
    implicit none

    integer,  intent(in)  :: ijdim
    real(DP), intent(out) :: u (ijdim,kdim)
    real(DP), intent(out) :: v (ijdim,kdim)
    real(DP), intent(in)  :: vx(ijdim,kdim)
    real(DP), intent(in)  :: vy(ijdim,kdim)
    real(DP), intent(in)  :: vz(ijdim,kdim)
    real(DP), intent(in)  :: ix(ijdim)
    real(DP), intent(in)  :: iy(ijdim)
    real(DP), intent(in)  :: iz(ijdim)
    real(DP), intent(in)  :: jx(ijdim)
    real(DP), intent(in)  :: jy(ijdim)
    real(DP), intent(in)  :: jz(ijdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------
    !!$acc data pcopyin(vx,vy,vz,ix,iy,iz) if(acc_is_present(u))

!OCL XFILL
    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,ix,iy,iz) if(acc_is_present(u))
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,u,v,vx,vy,vz,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       u(ij,k) = vx(ij,k) * ix(ij) &
               + vy(ij,k) * iy(ij) &
               + vz(ij,k) * iz(ij)
       v(ij,k) = vx(ij,k) * jx(ij) &
               + vy(ij,k) * jy(ij) &
               + vz(ij,k) * jz(ij)
    enddo
    enddo
    !$omp end parallel do
    !$acc end kernels

    if ( acc_wait_on ) then
       !$acc wait(1)
    endif
    !!$acc end data
    return
  end subroutine CNVVAR_vh2uv_in_DP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_uv2vh_2D_SP( &
!MINIAPP       ucos, ucos_pl, &
!MINIAPP       vcos, vcos_pl, &
!MINIAPP       vx,   vx_pl,   &
!MINIAPP       vy,   vy_pl,   &
!MINIAPP       vz,   vz_pl,   &
!MINIAPP       withcos        )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(in)  :: ucos   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vcos   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vx     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vy     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: vz     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, l, k0
!MINIAPP
!MINIAPP    real(SP) :: sw
!MINIAPP    real(SP) :: u
!MINIAPP    real(SP) :: v
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(ucos,vcos,GRD_LAT,GMTR_p) if(acc_is_present(vx))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_SP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_SP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(vx))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(ucos,vcos,coslat,GMTR_p) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       sw = 0.5_SP + sign(0.5_SP,real(-abs(coslat(g,l)),kind=SP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP       u = ucos(g,k0,l) * ( 1.0_SP - sw ) / ( coslat(g,l) - sw )
!MINIAPP       v = vcos(g,k0,l) * ( 1.0_SP - sw ) / ( coslat(g,l) - sw )
!MINIAPP
!MINIAPP       vx(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JX)
!MINIAPP       vy(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JY)
!MINIAPP       vz(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IZ) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JZ)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          sw = 0.5_SP + sign(0.5_SP,real(-abs(coslat_pl(g,l)),kind=SP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP          u = ucos_pl(g,k0,l) * ( 1.0_SP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP          v = vcos_pl(g,k0,l) * ( 1.0_SP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP
!MINIAPP          vx_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JX)
!MINIAPP          vy_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JY)
!MINIAPP          vz_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IZ) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JZ)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_uv2vh_2D_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_uv2vh_2D_DP( &
!MINIAPP       ucos, ucos_pl, &
!MINIAPP       vcos, vcos_pl, &
!MINIAPP       vx,   vx_pl,   &
!MINIAPP       vy,   vy_pl,   &
!MINIAPP       vz,   vz_pl,   &
!MINIAPP       withcos        )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(in)  :: ucos   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vcos   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vx     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vy     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: vz     (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, l, k0
!MINIAPP
!MINIAPP    real(DP) :: sw
!MINIAPP    real(DP) :: u
!MINIAPP    real(DP) :: v
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(ucos,vcos,GRD_LAT,GMTR_p) if(acc_is_present(vx))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_DP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_DP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(vx))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(vx,vy,vz) pcopyin(ucos,vcos,coslat,GMTR_p) if(acc_is_present(vx))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       sw = 0.5_DP + sign(0.5_DP,real(-abs(coslat(g,l)),kind=DP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP       u = ucos(g,k0,l) * ( 1.0_DP - sw ) / ( coslat(g,l) - sw )
!MINIAPP       v = vcos(g,k0,l) * ( 1.0_DP - sw ) / ( coslat(g,l) - sw )
!MINIAPP
!MINIAPP       vx(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JX)
!MINIAPP       vy(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JY)
!MINIAPP       vz(g,k0,l) = u * GMTR_p(g,k0,l,GMTR_p_IZ) &
!MINIAPP                  + v * GMTR_p(g,k0,l,GMTR_p_JZ)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          sw = 0.5_DP + sign(0.5_DP,real(-abs(coslat_pl(g,l)),kind=DP)) ! if (coslat == 0), u=v=0
!MINIAPP
!MINIAPP          u = ucos_pl(g,k0,l) * ( 1.0_DP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP          v = vcos_pl(g,k0,l) * ( 1.0_DP - sw ) / ( coslat_pl(g,l) - sw )
!MINIAPP
!MINIAPP          vx_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JX)
!MINIAPP          vy_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JY)
!MINIAPP          vz_pl(g,k0,l) = u * GMTR_p_pl(g,k0,l,GMTR_p_IZ) &
!MINIAPP                        + v * GMTR_p_pl(g,k0,l,GMTR_p_JZ)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_uv2vh_2D_DP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_vh2uv_2D_SP( &
!MINIAPP       u,  u_pl,  &
!MINIAPP       v,  v_pl,  &
!MINIAPP       vx, vx_pl, &
!MINIAPP       vy, vy_pl, &
!MINIAPP       vz, vz_pl, &
!MINIAPP       withcos    )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(SP), intent(out) :: u    (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: u_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(out) :: v    (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(out) :: v_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vx   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vy   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(SP), intent(in)  :: vz   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(SP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(SP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(SP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, l, k0
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(vx,vy,vz,GRD_LAT,GMTR_p) if(acc_is_present(u))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_SP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_SP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(u))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,coslat,GMTR_p) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       u(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
!MINIAPP       v(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
!MINIAPP                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
!MINIAPP                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          u_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
!MINIAPP          v_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
!MINIAPP                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
!MINIAPP                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_vh2uv_2D_SP

  !-----------------------------------------------------------------------------
  !--- MINIAPP: whole-grid variant removed: it needs VMTR_* / GMTR_*, i.e. mod_vmtr + mod_gmtr + mod_oprt + mod_comm.  The microphysics chain only calls the '_in' (region, gall_in) variants, which need nothing but GRD_afact/GRD_bfact. ---
!MINIAPP  subroutine CNVVAR_vh2uv_2D_DP( &
!MINIAPP       u,  u_pl,  &
!MINIAPP       v,  v_pl,  &
!MINIAPP       vx, vx_pl, &
!MINIAPP       vy, vy_pl, &
!MINIAPP       vz, vz_pl, &
!MINIAPP       withcos    )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_KNONE,   &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_lall,    &
!MINIAPP       ADM_lall_pl, &
!MINIAPP       ADM_gall,    &
!MINIAPP       ADM_gall_pl
!MINIAPP    use mod_grd, only: &
!MINIAPP       GRD_LAT,   &
!MINIAPP       GRD_LAT_pl
!MINIAPP    use mod_gmtr, only: &
!MINIAPP       GMTR_p_IX, &
!MINIAPP       GMTR_p_IY, &
!MINIAPP       GMTR_p_IZ, &
!MINIAPP       GMTR_p_JX, &
!MINIAPP       GMTR_p_JY, &
!MINIAPP       GMTR_p_JZ, &
!MINIAPP       GMTR_p,    &
!MINIAPP       GMTR_p_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(DP), intent(out) :: u    (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: u_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(out) :: v    (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(out) :: v_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vx   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vy   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP    real(DP), intent(in)  :: vz   (ADM_gall   ,ADM_KNONE,ADM_lall   )
!MINIAPP    real(DP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
!MINIAPP
!MINIAPP    logical,  intent(in), optional :: withcos
!MINIAPP
!MINIAPP    real(DP) :: coslat   (ADM_gall   ,ADM_lall   )
!MINIAPP    real(DP) :: coslat_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    integer  :: g, l, k0
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP    !$acc data pcreate(coslat) pcopyin(vx,vy,vz,GRD_LAT,GMTR_p) if(acc_is_present(u))
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(coslat) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       coslat(g,l) = 1.0_DP
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP    coslat_pl(:,:) = 1.0_DP
!MINIAPP
!MINIAPP    if ( present(withcos) ) then
!MINIAPP       if ( withcos ) then
!MINIAPP          !$acc kernels async(1) pcopy(coslat) pcopyin(GRD_LAT) if(acc_is_present(u))
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do g = 1, ADM_gall
!MINIAPP             coslat(g,l) = cos(GRD_LAT(g,l))
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          !$acc end kernels
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             coslat_pl(:,:) = cos(GRD_LAT_pl(:,:))
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc kernels async(1) pcopy(u,v) pcopyin(vx,vy,vz,coslat,GMTR_p) if(acc_is_present(u))
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       u(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IX) &
!MINIAPP                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IY) &
!MINIAPP                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_IZ) ) * coslat(g,l)
!MINIAPP       v(g,k0,l) = ( vx(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JX) &
!MINIAPP                   + vy(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JY) &
!MINIAPP                   + vz(g,k0,l) * GMTR_p(g,k0,l,GMTR_p_JZ) ) * coslat(g,l)
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP    !$acc end kernels
!MINIAPP
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          u_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IX) &
!MINIAPP                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IY) &
!MINIAPP                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_IZ) ) * coslat_pl(g,l)
!MINIAPP          v_pl(g,k0,l) = ( vx_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JX) &
!MINIAPP                         + vy_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JY) &
!MINIAPP                         + vz_pl(g,k0,l) * GMTR_p_pl(g,k0,l,GMTR_p_JZ) ) * coslat_pl(g,l)
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !$acc wait(1)
!MINIAPP    !$acc end data
!MINIAPP    return
!MINIAPP  end subroutine CNVVAR_vh2uv_2D_DP

end module mod_cnvvar
