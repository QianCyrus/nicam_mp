!-------------------------------------------------------------------------------
!> module PRECISION
!!
!! @par Description
!!         Precision definition
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_precision
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: SP = kind(0.E0) ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = kind(0.D0) ! Double Precision: kind(0.D0)

  integer, public, parameter :: SP_PREC = precision(0.E0)
  integer, public, parameter :: DP_PREC = precision(0.D0)

  ! Default setting

#ifdef SINGLE
  integer, public, parameter :: RP      = SP
  integer, public, parameter :: RP_PREC = SP_PREC
#else
  integer, public, parameter :: RP      = DP
  integer, public, parameter :: RP_PREC = DP_PREC
#endif

  ! Individual setting

#ifdef SINGLE_DYN
  integer, public, parameter :: RP_DYN      = SP
  integer, public, parameter :: RP_DYN_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_DYN      = SP
  integer, public, parameter :: RP_DYN_PREC = SP_PREC
#else
  integer, public, parameter :: RP_DYN      = DP
  integer, public, parameter :: RP_DYN_PREC = DP_PREC
#endif

#ifdef SINGLE_PHY
  integer, public, parameter :: RP_PHY      = SP
  integer, public, parameter :: RP_PHY_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_PHY      = SP
  integer, public, parameter :: RP_PHY_PREC = SP_PREC
#else
  integer, public, parameter :: RP_PHY      = DP
  integer, public, parameter :: RP_PHY_PREC = DP_PREC
#endif

#ifdef SINGLE_CP
  integer, public, parameter :: RP_CP      = SP
  integer, public, parameter :: RP_CP_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_CP      = SP
  integer, public, parameter :: RP_CP_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_CP      = SP
  integer, public, parameter :: RP_CP_PREC = SP_PREC
#else
  integer, public, parameter :: RP_CP      = DP
  integer, public, parameter :: RP_CP_PREC = DP_PREC
#endif

#ifdef SINGLE_MP
  integer, public, parameter :: RP_MP      = SP
  integer, public, parameter :: RP_MP_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_MP      = SP
  integer, public, parameter :: RP_MP_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_MP      = SP
  integer, public, parameter :: RP_MP_PREC = SP_PREC
#else
  integer, public, parameter :: RP_MP      = DP
  integer, public, parameter :: RP_MP_PREC = DP_PREC
#endif

#ifdef SINGLE_RD
  integer, public, parameter :: RP_RD      = SP
  integer, public, parameter :: RP_RD_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_RD      = SP
  integer, public, parameter :: RP_RD_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_RD      = SP
  integer, public, parameter :: RP_RD_PREC = SP_PREC
#else
  integer, public, parameter :: RP_RD      = DP
  integer, public, parameter :: RP_RD_PREC = DP_PREC
#endif

#ifdef SINGLE_TB
  integer, public, parameter :: RP_TB      = SP
  integer, public, parameter :: RP_TB_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_TB      = SP
  integer, public, parameter :: RP_TB_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_TB      = SP
  integer, public, parameter :: RP_TB_PREC = SP_PREC
#else
  integer, public, parameter :: RP_TB      = DP
  integer, public, parameter :: RP_TB_PREC = DP_PREC
#endif

#ifdef SINGLE_SF
  integer, public, parameter :: RP_SF      = SP
  integer, public, parameter :: RP_SF_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_SF      = SP
  integer, public, parameter :: RP_SF_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_SF      = SP
  integer, public, parameter :: RP_SF_PREC = SP_PREC
#else
  integer, public, parameter :: RP_SF      = DP
  integer, public, parameter :: RP_SF_PREC = DP_PREC
#endif

#ifdef SINGLE_AE
  integer, public, parameter :: RP_AE      = SP
  integer, public, parameter :: RP_AE_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_AE      = SP
  integer, public, parameter :: RP_AE_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_AE      = SP
  integer, public, parameter :: RP_AE_PREC = SP_PREC
#else
  integer, public, parameter :: RP_AE      = DP
  integer, public, parameter :: RP_AE_PREC = DP_PREC
#endif

#ifdef SINGLE_CH
  integer, public, parameter :: RP_CH      = SP
  integer, public, parameter :: RP_CH_PREC = SP_PREC
#elif SINGLE_PHY
  integer, public, parameter :: RP_CH      = SP
  integer, public, parameter :: RP_CH_PREC = SP_PREC
#elif SINGLE
  integer, public, parameter :: RP_CH      = SP
  integer, public, parameter :: RP_CH_PREC = SP_PREC
#else
  integer, public, parameter :: RP_CH      = DP
  integer, public, parameter :: RP_CH_PREC = DP_PREC
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module mod_precision
