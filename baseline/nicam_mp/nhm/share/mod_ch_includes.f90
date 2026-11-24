!-------------------------------------------------------------------------------
!> Module chem includes
!!
!! @par Description
!!          This module contains the chemistry tracer variables.
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_ch_includes
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
  public :: CHEM_tracerindex

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef DISABLE_CHEM
  integer, public, parameter :: CH_nreac_photo = 1
#else
  integer, public, parameter :: CH_nreac_photo = 41
#endif
  integer, public, parameter :: CH_TRC_VMAX    = 36

  character(len=H_SHORT), public :: CH_TRC_name(CH_TRC_VMAX) ! short name  of tracers
  character(len=H_MID),   public :: CH_TRC_desc(CH_TRC_VMAX) ! description of tracers

  character(len=H_SHORT), public :: CH_SLS_name(CH_TRC_VMAX) ! short name  of short-lived species
  character(len=H_MID),   public :: CH_SLS_desc(CH_TRC_VMAX) ! description of short-lived species

  real(DP),               public :: CH_emisfactorA(CH_TRC_VMAX)

  data CH_TRC_name / 'OX',       &
                     'NOX',      &
                     'N2O5',     &
                     'HNO3',     &
                     'HNO4',     &
                     'H2O2',     &
                     'CO',       &
                     'C2H6',     &
                     'C3H8',     &
                     'C2H4',     &
                     'C3H6',     &
                     'ONMV',     &
                     'C5H8',     &
                     'C10H16',   &
                     'CH3COCH3', &
                     'CH2O',     &
                     'CH3CHO',   &
                     'CH3OH',    &
                     'NALD',     &
                     'MGLY',     &
                     'HACET',    &
                     'MACR',     &
                     'PAN',      &
                     'MPAN',     &
                     'ISON',     &
                     'CH3OOH',   &
                     'C2H5OOH',  &
                     'C3H7OOH',  &
                     'ISOOH',    &
                     'HOROOH',   &
                     'CH3COOOH', &
                     'MACROOH',  &
                     'OXS',      &
                     'SO2',      &
                     'DMS',      &
                     'SO4'       /

  data CH_TRC_desc / 'Ozone family(Ox)',           &
                     'NOx',                        &
                     'Nitrogen pentoxide',         &
                     'Nitric acid',                &
                     'Peroxy nitric acid ',        &
                     'Hydrogen peroxide',          &
                     'Carbon monoxide',            &
                     'Ethane',                     &
                     'Propane',                    &
                     'Ethene',                     &
                     'Propene',                    &
                     'Other VOCs',                 &
                     'Isoprene',                   &
                     'Monoterpenes',               &
                     'Acetone',                    &
                     'Formaldehyde',               &
                     'Acetaldehyde',               &
                     'Methanol',                   &
                     'Nitrooxy acetaldehyde',      &
                     'Methyl glyoxal',             &
                     'Hydroxy acetone',            &
                     'Methacrolein',               &
                     'Peroxyacetyl nitrate',       &
                     'Peroxyacyl nitrates',        &
                     'Hydroxyalkyl nitrates',      &
                     'Methyl hydroperoxide',       &
                     'Ethyl  hydroperoxide',       &
                     'Propyl hydroperoxide',       &
                     'ISOP-origin hydroperoxides', &
                     'Alkane-origin hydroxy hydroper.', &
                     'Peracetic acid',             &
                     'MACR-origin hydroperoxides', &
                     'O3 from stratosphere',       &
                     'Sulfur dioxide ',            &
                     'Dimethyl sulfide',           &
                     'Sulfate SO4(--)'             /

  data CH_emisfactorA / 1.D0,       & ! OX
                        3.285714D0, & ! NOX
                        1.D0,       & ! N2O5
                        1.D0,       & ! HNO3
                        1.D0,       & ! NO2HO2
                        1.D0,       & ! HOOH
                        1.D0,       & ! CO
                        1.250000D0, & ! C2H6
                        1.222222D0, & ! C3H8
                        1.166667D0, & ! C2H4
                        1.166667D0, & ! C3H6
                        1.D0,       & ! ONMV
                        1.D0,       & ! C5H8
                        1.D0,       & ! C10H16
                        1.611111D0, & ! CH3COCH3
                        1.D0,       & ! HCHO
                        1.D0,       & ! CH3CHO
                        1.D0,       & ! CH3OH
                        1.D0,       & ! NALD
                        1.D0,       & ! CH3COCHO
                        1.D0,       & ! CH3COCH2OH
                        1.D0,       & ! MACR
                        1.D0,       & ! PAN
                        1.D0,       & ! MPAN
                        1.D0,       & ! ISOPN
                        1.D0,       & ! CH3OOH
                        1.D0,       & ! C2H5OOH
                        1.D0,       & ! C3H7OOH
                        1.D0,       & ! HOROOH
                        1.D0,       & ! ISOPOOH
                        1.D0,       & ! CH3COOOH
                        1.D0,       & ! MACROOH
                        1.D0,       & ! OXS
                        2.D0,       & ! SO2
                        1.D0,       & ! DMS
                        1.D0        / ! SO4

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
  !> Return tracer ID
  subroutine CHEM_tracerindex( &
       varname, &
       varid    )
    implicit none

    character(len=*), intent(in)  :: varname
    integer,          intent(out) :: varid

    integer :: n
    !---------------------------------------------------------------------------

    !--- search existing item
    varid = -1
    do n = 1, CH_TRC_VMAX
       if ( varname == CH_TRC_name(n) ) then ! match existing item
          varid = n
          return
       endif
    enddo

    return
  end subroutine CHEM_tracerindex

end module mod_ch_includes
