!-------------------------------------------------------------------------------
!> Module GRID
!!
!! @par Description
!!         Management of the icosahedral grid system
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_grd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_stdio

  use mod_adm, only: &
     ADM_nxyz,    &
     ADM_TI,      &
     ADM_TJ,      &
     ADM_AI,      &
     ADM_AIJ,     &
     ADM_AJ,      &
     ADM_KNONE,   &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: GRD_setup
  !--- MINIAPP: horizontal-grid entry points removed (see below). ---
!MINIAPP  public :: GRD_input_hgrid
!MINIAPP  public :: GRD_output_hgrid
  public :: GRD_input_vgrid
  public :: GRD_output_vgrid
!MINIAPP  public :: GRD_makelatlon
!MINIAPP  public :: GRD_setlatlon
!MINIAPP  public :: GRD_scaling

  interface GRD_input_vgrid
    module procedure GRD_input_vgrid_SP
    module procedure GRD_input_vgrid_DP
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Indentifiers for the directions in the Cartesian coordinate
  integer,  public, parameter :: GRD_XDIR = 1
  integer,  public, parameter :: GRD_YDIR = 2
  integer,  public, parameter :: GRD_ZDIR = 3

  ! Indentifiers for the directions in the spherical coordinate
  integer,  public, parameter :: I_LAT = 1
  integer,  public, parameter :: I_LON = 2

  !====== Horizontal Grid ======
  !
  ! Grid points ( X: CELL CENTER )
  !          .___.
  !         /     \
  !        .   p   .
  !         \ ___ /
  !          '   '
  !
  ! Grid points ( Xt: CELL VERTEX )
  !          p___p
  !         /     \
  !        p       p
  !         \ ___ /
  !          p   p
  !
  ! Grid points ( Xr: CELL ARC )
  !          ._p_.
  !         p     p
  !        .       .
  !         p _ _ p
  !          ' p '
  !
  real(RP), public              :: GRD_rscale ! scaling factor for the radius of the sphere

#ifdef _FIXEDINDEX_
  real(RP), public              :: GRD_x     (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz)
  real(RP), public              :: GRD_x_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
  real(RP), public              :: GRD_xt    (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz)
  real(RP), public              :: GRD_xt_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
  real(RP), public              :: GRD_xr    (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,ADM_nxyz)
  real(RP), public              :: GRD_xr_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)

  real(RP), public              :: GRD_s     (ADM_gall   ,ADM_KNONE,ADM_lall   ,              2)
  real(RP), public              :: GRD_s_pl  (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              2)
  real(RP), public              :: GRD_st    (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,2)
  real(RP), public              :: GRD_st_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              2)

  real(RP), public              :: GRD_LAT   (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GRD_LAT_pl(ADM_gall_pl,ADM_lall_pl)
  real(RP), public              :: GRD_LON   (ADM_gall   ,ADM_lall   )
  real(RP), public              :: GRD_LON_pl(ADM_gall_pl,ADM_lall_pl)
#else
  real(RP), public, allocatable :: GRD_x     (:,:,:,:)
  real(RP), public, allocatable :: GRD_x_pl  (:,:,:,:)
  real(RP), public, allocatable :: GRD_xt    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_xt_pl (:,:,:,:)
  real(RP), public, allocatable :: GRD_xr    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_xr_pl (:,:,:,:)

  real(RP), public, allocatable :: GRD_s     (:,:,:,:)
  real(RP), public, allocatable :: GRD_s_pl  (:,:,:,:)
  real(RP), public, allocatable :: GRD_st    (:,:,:,:,:)
  real(RP), public, allocatable :: GRD_st_pl (:,:,:,:)

  real(RP), public, allocatable :: GRD_LAT   (:,:)
  real(RP), public, allocatable :: GRD_LAT_pl(:,:)
  real(RP), public, allocatable :: GRD_LON   (:,:)
  real(RP), public, allocatable :: GRD_LON_pl(:,:)
#endif



  !====== Topography ======
  integer,  public, parameter   :: GRD_ZSFC = 1
  integer,  public, parameter   :: GRD_ZSD  = 2

#ifdef _FIXEDINDEX_
  real(RP), public              :: GRD_zs   (ADM_gall   ,ADM_KNONE,ADM_lall   ,GRD_ZSFC:GRD_ZSD)
  real(RP), public              :: GRD_zs_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,GRD_ZSFC:GRD_ZSD)
#else
  real(RP), public, allocatable :: GRD_zs   (:,:,:,:)
  real(RP), public, allocatable :: GRD_zs_pl(:,:,:,:)
#endif



  !====== Vertical Grid ======
  !
  ! z coordinate ( actual height )
  !
  !         GRD_vz(1:ADM_gall,  &
  !                1:ADM_kall,  &
  !                1:ADM_lall,  &
  !                GRD_Z:GRD_ZH )
  !
  !         GRD_vz_pl(1:ADM_gall_pl, &
  !                   1:ADM_kall,    &
  !                   1:ADM_lall_pl, &
  !                   GRD_Z:GRD_ZH   )
  !

  integer,  public, parameter   :: GRD_Z  = 1
  integer,  public, parameter   :: GRD_ZH = 2

  real(RP), public              :: GRD_htop ! model top height [m]

#ifdef _FIXEDINDEX_
  real(RP), public              :: GRD_gz   (ADM_kall)
  real(RP), public              :: GRD_gzh  (ADM_kall)
  real(RP), public              :: GRD_dgz  (ADM_kall)
  real(RP), public              :: GRD_dgzh (ADM_kall)
  real(RP), public              :: GRD_rdgz (ADM_kall)
  real(RP), public              :: GRD_rdgzh(ADM_kall)

  real(RP), public              :: GRD_afact(ADM_kall)
  real(RP), public              :: GRD_bfact(ADM_kall)
  real(RP), public              :: GRD_cfact(ADM_kall)
  real(RP), public              :: GRD_dfact(ADM_kall)

  real(RP), public              :: GRD_vz   (ADM_gall   ,ADM_kall,ADM_lall   ,GRD_Z:GRD_ZH)
  real(RP), public              :: GRD_vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_Z:GRD_ZH)
#else
  real(RP), public, allocatable :: GRD_gz   (:) ! gsi (z-star) coordinate
  real(RP), public, allocatable :: GRD_gzh  (:) ! gsi (z-star) coordinate at the half point
  real(RP), public, allocatable :: GRD_dgz  (:) ! d(gsi)
  real(RP), public, allocatable :: GRD_dgzh (:) ! d(gsi) at the half point
  real(RP), public, allocatable :: GRD_rdgz (:)
  real(RP), public, allocatable :: GRD_rdgzh(:)

  real(RP), public, allocatable :: GRD_afact(:) ! From the cell center value to the cell wall value
  real(RP), public, allocatable :: GRD_bfact(:) !  A(k-1/2) = ( afact(k) A(k) + bfact(k) * A(k-1) )
  real(RP), public, allocatable :: GRD_cfact(:) ! From the cell wall value to the cell center value
  real(RP), public, allocatable :: GRD_dfact(:) !  A(k) = ( cfact(k) A(k+1/2) + dfact(k) * A(k-1/2) )

  real(RP), public, allocatable :: GRD_vz   (:,:,:,:)
  real(RP), public, allocatable :: GRD_vz_pl(:,:,:,:)
#endif

  integer,  public, parameter   :: GRD_grid_type_on_sphere = 1
  integer,  public, parameter   :: GRD_grid_type_on_plane  = 2
  integer,  public              :: GRD_grid_type           = GRD_grid_type_on_sphere
  logical,  public              :: GRD_channel_1dmd        = .false.

  !%%%%% Obsolete variables for backward compatibility %%%%%
  real(RP), public, allocatable :: GRD_afac(:) !=> 2 x GRD_afact(:)
  real(RP), public, allocatable :: GRD_bfac(:) !=> 2 x GRD_bfact(:)
  real(RP), public, allocatable :: GRD_cfac(:) !=> 2 x GRD_cfact(:)
  real(RP), public, allocatable :: GRD_dfac(:) !=> 2 x GRD_dfact(:)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  real(RP), public              :: GRD_s1_sleve            = 15000.0_RP ! e-folding height for SLEVE vertical coordinate [m]
  real(RP), public              :: GRD_s2_sleve            = 2500.0_RP  ! e-folding height for SLEVE vertical coordinate [m]
  real(RP), public              :: GRD_ns_sleve            = 1.3_RP     ! SLEVE vertical coordinate parameter
  real(RP), public              :: GRD_hflat_sleve         = -999.0_RP  ! [m]
  character(len=H_SHORT), public:: vgrid_scheme            = 'LINEAR'   ! vertical coordinate scheme
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
!MINIAPP  private :: GRD_gen_plgrid
!MINIAPP  private :: GRD_makearc

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: hgrid_io_mode       = 'LEGACY'
  character(len=H_SHORT), private :: topo_io_mode        = 'LEGACY'
  logical,                private :: hgrid_direct_access = .false.
  logical,                private :: topo_direct_access  = .false.
  character(len=H_FILE),  private :: hgrid_fname         = ''         ! horizontal grid file
  character(len=H_FILE),  private :: topo_fname          = ''         ! topography file
  character(len=H_FILE),  private :: toposd_fname        = ''         ! topography file

  character(len=H_FILE),  private :: vgrid_fname         = ''         ! vertical grid file
  real(RP),               private :: h_efold             = 10000.0_RP ! e-folding height for hybrid vertical coordinate [m]
  real(RP),               private :: hflat               =  -999.0_RP ! [m]
  logical,                private :: output_vgrid        = .false.    ! output verical grid file?

  logical,                private :: hgrid_comm_flg      = .true.     ! communicate GRD_x?          [add] T.Ohno 110722
  real(RP),               private :: triangle_size       = 0.0_RP     ! length of sides of triangle [add] T.Ohno 110722

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  !--- MINIAPP: ------------------------------------------------------------
  !--- The original GRD_setup reads the horizontal grid (FIO), the topography,
  !--- builds GRD_vz by the terrain-following transformation and fills the halo
  !--- with COMM_data_transfer.  The microphysics chain reads NONE of that: the
  !--- terrain-following heights z / zh arrive as arguments (from the dump), and
  !--- the only GRD_* symbols reachable are the 1-D vertical arrays
  !--- GRD_dgz/dgzh/rdgz/rdgzh/afact/bfact/cfact/dfact + GRD_gz/gzh, which are
  !--- pure arithmetic on the vgrid file.  The vertical block below is copied
  !--- verbatim from the original; everything else is dropped.
  !--- The original body is preserved, commented out, after this routine.
  !-------------------------------------------------------------------------
  subroutine GRD_setup
    use mod_adm, only: &
       ADM_proc_stop,    &
       ADM_KNONE,        &
       ADM_lall,         &
       ADM_kall,         &
       ADM_kmin,         &
       ADM_kmax,         &
       ADM_IsMaster
    implicit none

    logical :: channel_1dmd = .false.

    namelist / GRDPARAM / &
       GRD_grid_type,       &
       channel_1dmd,        &
       hgrid_io_mode,       &
       topo_io_mode,        &
       hgrid_direct_access, &
       topo_direct_access,  &
       hgrid_fname,         &
       topo_fname,          &
       toposd_fname,        &
       vgrid_fname,         &
       vgrid_scheme,        &
       h_efold,             &
       hflat,               &
       output_vgrid,        &
       hgrid_comm_flg,      &
       triangle_size,       &
       GRD_s1_sleve,        &
       GRD_s2_sleve,        &
       GRD_ns_sleve,        &
       GRD_hflat_sleve

    integer  :: ierr
    integer  :: k
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[grd]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=GRDPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** GRDPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist GRDPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist GRDPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=GRDPARAM)

    GRD_channel_1dmd = channel_1dmd

    !---< vertical coordinate >--- (verbatim from the original GRD_setup)
    if ( ADM_kall /= ADM_KNONE ) then
       allocate( GRD_gz   (ADM_kall) )
       allocate( GRD_gzh  (ADM_kall) )
       allocate( GRD_dgz  (ADM_kall) )
       allocate( GRD_dgzh (ADM_kall) )
       allocate( GRD_rdgz (ADM_kall) )
       allocate( GRD_rdgzh(ADM_kall) )

       allocate( GRD_afact(ADM_kall) )
       allocate( GRD_bfact(ADM_kall) )
       allocate( GRD_cfact(ADM_kall) )
       allocate( GRD_dfact(ADM_kall) )

       if ( trim(vgrid_fname) /= 'NONE' ) then
          call GRD_input_vgrid( vgrid_fname, & ! [IN]
                                ADM_kall,    & ! [IN]
                                GRD_htop,    & ! [IN]
                                GRD_gz (:),  & ! [OUT]
                                GRD_gzh(:)   ) ! [OUT]
       else
          ! dummy
          do k = 1, ADM_kall
             GRD_gz (k) = dble(k)
             GRD_gzh(k) = dble(k) - 0.5_RP
          end do
       end if

       ! calculation of grid intervals ( cell center )
       do k = ADM_kmin-1, ADM_kmax
          GRD_dgz(k) = GRD_gzh(k+1) - GRD_gzh(k)
       enddo
       GRD_dgz(ADM_kmax+1) = GRD_dgz(ADM_kmax)

       ! calculation of grid intervals ( cell wall )
       do k = ADM_kmin, ADM_kmax+1
          GRD_dgzh(k) = GRD_gz(k) - GRD_gz(k-1)
       enddo
       GRD_dgzh(ADM_kmin-1) = GRD_dgzh(ADM_kmin)

       do k = 1, ADM_kall
          GRD_rdgz (k) = 1.0_RP / grd_dgz (k)
          GRD_rdgzh(k) = 1.0_RP / grd_dgzh(k)
       enddo

       ! vertical interpolation factor
       do k = ADM_kmin, ADM_kmax+1
          GRD_afact(k) = ( GRD_gzh(k) - GRD_gz(k-1) ) &
                       / ( GRD_gz (k) - GRD_gz(k-1) )
       enddo
       GRD_afact(ADM_kmin-1) = 1.0_RP

       GRD_bfact(:) = 1.0_RP - GRD_afact(:)

       do k = ADM_kmin, ADM_kmax
          GRD_cfact(k) = ( GRD_gz (k  ) - GRD_gzh(k) ) &
                       / ( GRD_gzh(k+1) - GRD_gzh(k) )
       enddo
       GRD_cfact(ADM_kmin-1) = 1.0_RP
       GRD_cfact(ADM_kmax+1) = 0.0_RP

       GRD_dfact(:) = 1.0_RP - GRD_cfact(:)

       !--- MINIAPP: the terrain-following GRD_vz(:,:,:,GRD_Z/GRD_ZH) transform
       !--- and its COMM_data_transfer halo fill are dropped; the mini-app takes
       !--- z / zh straight from the dump (mpF_z.bin / mpF_zh.bin).
    else
       allocate( GRD_gz (ADM_KNONE) )
       allocate( GRD_gzh(ADM_KNONE) )
       GRD_gz (:) = 1.0_RP
       GRD_gzh(:) = 1.0_RP
    endif

    !--- output information about grid.
    if ( ADM_kall /= ADM_KNONE ) then
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,'(1x,A)')             '|======      Vertical Coordinate [m]      ======|'
       write(IO_FID_LOG,'(1x,A)')             '|                                               |'
       write(IO_FID_LOG,'(1x,A)')             '|          -GRID CENTER-       -GRID INTERFACE- |'
       write(IO_FID_LOG,'(1x,A)')             '|  k        gz     d(gz)      gzh    d(gzh)   k |'
       write(IO_FID_LOG,'(1x,A)')             '|                                               |'
       k = ADM_kmax + 1
       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | TOA'
       k = ADM_kmax
       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmax'
       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
       do k = ADM_kmax-1, ADM_kmin+1, -1
       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        |'
       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
       enddo
       k = ADM_kmin
       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmin'
       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | ground'
       k = ADM_kmin-1
       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
       write(IO_FID_LOG,'(1x,A)')             '|===============================================|'
       write(IO_FID_LOG,*) '--- Vertical layer scheme = ', trim(vgrid_scheme)
       if ( vgrid_scheme == 'HYBRID' ) then
          write(IO_FID_LOG,*) '--- e-folding height = ', h_efold
       endif

       if ( output_vgrid ) then
          if ( ADM_IsMaster ) then
             call GRD_output_vgrid('./vgrid_used.dat')
          endif
       endif
    else
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- vartical layer = 1'
    endif

    !%%%%% Obsolete variables for backward compatibility %%%%%
    if ( ADM_kall /= ADM_KNONE ) then
       allocate( GRD_afac(ADM_kall) )
       allocate( GRD_bfac(ADM_kall) )
       allocate( GRD_cfac(ADM_kall) )
       allocate( GRD_dfac(ADM_kall) )
       GRD_afac(:) = 2.0_RP * GRD_afact(:)
       GRD_bfac(:) = 2.0_RP * GRD_bfact(:)
       GRD_cfac(:) = 2.0_RP * GRD_cfact(:)
       GRD_dfac(:) = 2.0_RP * GRD_dfact(:)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return
  end subroutine GRD_setup

  !--- MINIAPP: original GRD_setup, preserved verbatim but commented out. ----
!MINIAPP  subroutine GRD_setup
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_IsMaster, &
!MINIAPP       ADM_proc_stop
!MINIAPP    use mod_cnst, only: &
!MINIAPP       UNDEF  => CNST_UNDEF,  &
!MINIAPP       RADIUS => CNST_RADIUS
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_flg_scm,      &
!MINIAPP       ADM_HGRID_SYSTEM, & ! [add] 20210322 Junya D-NICAM
!MINIAPP       ADM_have_pl,      &
!MINIAPP       ADM_lall,         &
!MINIAPP       ADM_gall,         &
!MINIAPP       ADM_gall_1d,      &
!MINIAPP       ADM_gmin,         &
!MINIAPP       ADM_gmax,         &
!MINIAPP       ADM_gslf_pl,      &
!MINIAPP       ADM_kmin,         &
!MINIAPP       ADM_kmax
!MINIAPP    use mod_comm, only: &
!MINIAPP       COMM_data_transfer
!MINIAPP    use mod_fio, only: &
!MINIAPP       FIO_FORMAT,    &
!MINIAPP       FIO_set_layer
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    logical :: channel_1dmd = .false.
!MINIAPP
!MINIAPP    namelist / GRDPARAM / &
!MINIAPP       GRD_grid_type,       &
!MINIAPP       channel_1dmd,        &
!MINIAPP       hgrid_io_mode,       &
!MINIAPP       topo_io_mode,        &
!MINIAPP       hgrid_direct_access, &
!MINIAPP       topo_direct_access,  &
!MINIAPP       hgrid_fname,         &
!MINIAPP       topo_fname,          &
!MINIAPP       toposd_fname,        &
!MINIAPP       vgrid_fname,         &
!MINIAPP       vgrid_scheme,        &
!MINIAPP       h_efold,             &
!MINIAPP       hflat,               &
!MINIAPP       output_vgrid,        &
!MINIAPP       hgrid_comm_flg,      &
!MINIAPP       triangle_size,       &
!MINIAPP       GRD_s1_sleve,        &
!MINIAPP       GRD_s2_sleve,        &
!MINIAPP       GRD_ns_sleve,        &
!MINIAPP       GRD_hflat_sleve
!MINIAPP
!MINIAPP    real(DP), allocatable :: zaxis(:) ! = GRD_gz
!MINIAPP
!MINIAPP    real(RP) :: htop
!MINIAPP    integer  :: nstart, nend
!MINIAPP    integer  :: kflat
!MINIAPP
!MINIAPP    integer  :: ierr
!MINIAPP    integer  :: n, k, l, k0
!MINIAPP
!MINIAPP    integer  :: suf, i, j
!MINIAPP    suf(i,j) = ADM_gall_1d * (j-1) + i
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    !--- read parameters
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,*) '+++ Module[grd]/Category[common share]'
!MINIAPP    rewind(IO_FID_CONF)
!MINIAPP    read(IO_FID_CONF,nml=GRDPARAM,iostat=ierr)
!MINIAPP    if ( ierr < 0 ) then
!MINIAPP       write(IO_FID_LOG,*) '*** GRDPARAM is not specified. use default.'
!MINIAPP    elseif( ierr > 0 ) then
!MINIAPP       write(*         ,*) 'xxx Not appropriate names in namelist GRDPARAM. STOP.'
!MINIAPP       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist GRDPARAM. STOP.'
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=GRDPARAM)
!MINIAPP
!MINIAPP    GRD_channel_1dmd = channel_1dmd
!MINIAPP
!MINIAPP
!MINIAPP
!MINIAPP    !---< horizontal grid >---
!MINIAPP#ifndef _FIXEDINDEX_
!MINIAPP    allocate( GRD_x     (ADM_gall   ,k0,ADM_lall   ,              ADM_nxyz) )
!MINIAPP    allocate( GRD_x_pl  (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )
!MINIAPP    allocate( GRD_xt    (ADM_gall   ,k0,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz) )
!MINIAPP    allocate( GRD_xt_pl (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )
!MINIAPP    allocate( GRD_xr    (ADM_gall   ,k0,ADM_lall   ,ADM_AI:ADM_AJ,ADM_nxyz) )
!MINIAPP    allocate( GRD_xr_pl (ADM_gall_pl,k0,ADM_lall_pl,              ADM_nxyz) )
!MINIAPP
!MINIAPP    allocate( GRD_s     (ADM_gall   ,k0,ADM_lall   ,              2) )
!MINIAPP    allocate( GRD_s_pl  (ADM_gall_pl,k0,ADM_lall_pl,              2) )
!MINIAPP    allocate( GRD_st    (ADM_gall   ,k0,ADM_lall   ,ADM_TI:ADM_TJ,2) )
!MINIAPP    allocate( GRD_st_pl (ADM_gall_pl,k0,ADM_lall_pl,              2) )
!MINIAPP
!MINIAPP    allocate( GRD_LAT   (ADM_gall   ,ADM_lall   ) )
!MINIAPP    allocate( GRD_LAT_pl(ADM_gall_pl,ADM_lall_pl) )
!MINIAPP    allocate( GRD_LON   (ADM_gall   ,ADM_lall   ) )
!MINIAPP    allocate( GRD_LON_pl(ADM_gall_pl,ADM_lall_pl) )
!MINIAPP#endif
!MINIAPP    GRD_x     (:,:,:,:)   = UNDEF
!MINIAPP    GRD_x_pl  (:,:,:,:)   = UNDEF
!MINIAPP    GRD_xt    (:,:,:,:,:) = UNDEF
!MINIAPP    GRD_xt_pl (:,:,:,:)   = UNDEF
!MINIAPP    GRD_xr    (:,:,:,:,:) = UNDEF
!MINIAPP    GRD_xr_pl (:,:,:,:)   = UNDEF
!MINIAPP
!MINIAPP    GRD_s     (:,:,:,:)   = UNDEF
!MINIAPP    GRD_s_pl  (:,:,:,:)   = UNDEF
!MINIAPP    GRD_st    (:,:,:,:,:) = UNDEF
!MINIAPP    GRD_st_pl (:,:,:,:)   = UNDEF
!MINIAPP
!MINIAPP    GRD_LAT   (:,:)       = UNDEF
!MINIAPP    GRD_LAT_pl(:,:)       = UNDEF
!MINIAPP    GRD_LON   (:,:)       = UNDEF
!MINIAPP    GRD_LON_pl(:,:)       = UNDEF
!MINIAPP
!MINIAPP    if ( .NOT. ADM_flg_scm ) then
!MINIAPP       call GRD_input_hgrid( hgrid_fname,                        & ![IN]
!MINIAPP                             input_vertex  = .true.,             & ![IN]
!MINIAPP                             io_mode       = hgrid_io_mode,      & ![IN]
!MINIAPP                             direct_access = hgrid_direct_access ) ![IN]
!MINIAPP
!MINIAPP       ! fail safe for floating overflow
!MINIAPP       if ( GRD_grid_type /= GRD_grid_type_on_plane ) then
!MINIAPP       GRD_x (:,:,:,:)   = max( GRD_x (:,:,:,:)  , -9.99_RP )
!MINIAPP       GRD_xt(:,:,:,:,:) = max( GRD_xt(:,:,:,:,:), -9.99_RP )
!MINIAPP
!MINIAPP       if ( ADM_have_pl ) then
!MINIAPP          GRD_x_pl (:,:,:,:) = max( GRD_x_pl (:,:,:,:), -9.99_RP )
!MINIAPP          GRD_xt_pl(:,:,:,:) = max( GRD_xt_pl(:,:,:,:), -9.99_RP )
!MINIAPP       endif
!MINIAPP       endif
!MINIAPP    else
!MINIAPP       GRD_x (:,:,:,:)   = 0.0_RP
!MINIAPP       GRD_xt(:,:,:,:,:) = 0.0_RP
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! data transfer for GRD_x (note: do not communicate GRD_xt)
!MINIAPP    if( hgrid_comm_flg ) call COMM_data_transfer( GRD_x, GRD_x_pl )
!MINIAPP
!MINIAPP    ! scaling
!MINIAPP    if ( GRD_grid_type == GRD_grid_type_on_plane ) then
!MINIAPP       call GRD_scaling(triangle_size)
!MINIAPP    else
!MINIAPP       call GRD_scaling(RADIUS)
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! calc latitude/longitude of each grid point
!MINIAPP    if ( .NOT. ADM_flg_scm ) then
!MINIAPP       call GRD_makelatlon
!MINIAPP    else
!MINIAPP       GRD_s(:,:,:,:) = 0.0_RP
!MINIAPP    endif
!MINIAPP
!MINIAPP    call GRD_setlatlon
!MINIAPP
!MINIAPP    GRD_LAT   (:,:) = GRD_s   (:,k0,:,I_LAT)
!MINIAPP    GRD_LON   (:,:) = GRD_s   (:,k0,:,I_LON)
!MINIAPP    GRD_LAT_pl(:,:) = GRD_s_pl(:,k0,:,I_LAT)
!MINIAPP    GRD_LON_pl(:,:) = GRD_s_pl(:,k0,:,I_LON)
!MINIAPP
!MINIAPP    ! calc position of cell arc
!MINIAPP    if ( .NOT. ADM_flg_scm ) then
!MINIAPP       call GRD_makearc
!MINIAPP    else
!MINIAPP       GRD_xr(:,:,:,:,:) = 0.0_RP
!MINIAPP    endif
!MINIAPP
!MINIAPP
!MINIAPP
!MINIAPP    !---< surface height >---
!MINIAPP#ifndef _FIXEDINDEX_
!MINIAPP    allocate( GRD_zs   (ADM_gall,   k0,ADM_lall,   GRD_ZSFC:GRD_ZSD) )
!MINIAPP    allocate( GRD_zs_pl(ADM_gall_pl,k0,ADM_lall_pl,GRD_ZSFC:GRD_ZSD) )
!MINIAPP#endif
!MINIAPP    GRD_zs   (:,:,:,:) = 0.0_RP
!MINIAPP    GRD_zs_pl(:,:,:,:) = 0.0_RP
!MINIAPP
!MINIAPP    call GRD_input_topograph( topo_fname,   & ![IN]
!MINIAPP                              toposd_fname, & ![IN]
!MINIAPP                              topo_io_mode  ) ![IN]
!MINIAPP
!MINIAPP
!MINIAPP
!MINIAPP    !---< vertical coordinate >---
!MINIAPP    if ( ADM_kall /= ADM_KNONE ) then
!MINIAPP#ifndef _FIXEDINDEX_
!MINIAPP       allocate( GRD_gz   (ADM_kall) )
!MINIAPP       allocate( GRD_gzh  (ADM_kall) )
!MINIAPP       allocate( GRD_dgz  (ADM_kall) )
!MINIAPP       allocate( GRD_dgzh (ADM_kall) )
!MINIAPP       allocate( GRD_rdgz (ADM_kall) )
!MINIAPP       allocate( GRD_rdgzh(ADM_kall) )
!MINIAPP
!MINIAPP       allocate( GRD_afact(ADM_kall) )
!MINIAPP       allocate( GRD_bfact(ADM_kall) )
!MINIAPP       allocate( GRD_cfact(ADM_kall) )
!MINIAPP       allocate( GRD_dfact(ADM_kall) )
!MINIAPP
!MINIAPP       allocate( GRD_vz   (ADM_gall   ,ADM_kall,ADM_lall   ,GRD_Z:GRD_ZH) )
!MINIAPP       allocate( GRD_vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_Z:GRD_ZH) )
!MINIAPP#endif
!MINIAPP       if ( trim(vgrid_fname) /= 'NONE' ) then
!MINIAPP          call GRD_input_vgrid( vgrid_fname, & ! [IN]
!MINIAPP                                ADM_kall,    & ! [IN]
!MINIAPP                                GRD_htop,    & ! [IN]
!MINIAPP                                GRD_gz (:),  & ! [OUT]
!MINIAPP                                GRD_gzh(:)   ) ! [OUT]
!MINIAPP       else
!MINIAPP          ! dummy
!MINIAPP          do k = 1, ADM_kall
!MINIAPP             GRD_gz (k) = dble(k)
!MINIAPP             GRD_gzh(k) = dble(k) - 0.5_RP
!MINIAPP          end do
!MINIAPP       end if
!MINIAPP
!MINIAPP       ! calculation of grid intervals ( cell center )
!MINIAPP       do k = ADM_kmin-1, ADM_kmax
!MINIAPP          GRD_dgz(k) = GRD_gzh(k+1) - GRD_gzh(k)
!MINIAPP       enddo
!MINIAPP       GRD_dgz(ADM_kmax+1) = GRD_dgz(ADM_kmax)
!MINIAPP
!MINIAPP       ! calculation of grid intervals ( cell wall )
!MINIAPP       do k = ADM_kmin, ADM_kmax+1
!MINIAPP          GRD_dgzh(k) = GRD_gz(k) - GRD_gz(k-1)
!MINIAPP       enddo
!MINIAPP       GRD_dgzh(ADM_kmin-1) = GRD_dgzh(ADM_kmin)
!MINIAPP
!MINIAPP       do k = 1, ADM_kall
!MINIAPP          GRD_rdgz (k) = 1.0_RP / grd_dgz (k)
!MINIAPP          GRD_rdgzh(k) = 1.0_RP / grd_dgzh(k)
!MINIAPP       enddo
!MINIAPP
!MINIAPP       ! vertical interpolation factor
!MINIAPP       do k = ADM_kmin, ADM_kmax+1
!MINIAPP          GRD_afact(k) = ( GRD_gzh(k) - GRD_gz(k-1) ) &
!MINIAPP                       / ( GRD_gz (k) - GRD_gz(k-1) )
!MINIAPP       enddo
!MINIAPP       GRD_afact(ADM_kmin-1) = 1.0_RP
!MINIAPP
!MINIAPP       GRD_bfact(:) = 1.0_RP - GRD_afact(:)
!MINIAPP
!MINIAPP       do k = ADM_kmin, ADM_kmax
!MINIAPP          GRD_cfact(k) = ( GRD_gz (k  ) - GRD_gzh(k) ) &
!MINIAPP                       / ( GRD_gzh(k+1) - GRD_gzh(k) )
!MINIAPP       enddo
!MINIAPP       GRD_cfact(ADM_kmin-1) = 1.0_RP
!MINIAPP       GRD_cfact(ADM_kmax+1) = 0.0_RP
!MINIAPP
!MINIAPP       GRD_dfact(:) = 1.0_RP - GRD_cfact(:)
!MINIAPP
!MINIAPP       ! setup z-coordinate
!MINIAPP       if ( ADM_flg_scm ) then
!MINIAPP          nstart = 1
!MINIAPP          nend   = 1
!MINIAPP       elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! [add] 20210322 Junya D-NICAM
!MINIAPP          nstart = suf(ADM_gmin-1,ADM_gmin-1)
!MINIAPP          nend   = suf(ADM_gmax+1,ADM_gmax+1)
!MINIAPP       else ! normal
!MINIAPP          nstart = suf(ADM_gmin,ADM_gmin)
!MINIAPP          nend   = suf(ADM_gmax,ADM_gmax)
!MINIAPP       endif
!MINIAPP
!MINIAPP       select case(vgrid_scheme)
!MINIAPP       case('LINEAR')
!MINIAPP
!MINIAPP          ! linear transfromation : (Gal-Chen & Sommerville(1975)
!MINIAPP          !    gz = H(z-zs)/(H-zs) -> z = (H-zs)/H * gz + zs
!MINIAPP
!MINIAPP          kflat = -1
!MINIAPP          if ( hflat > 0.0_RP ) then ! default = -999.0
!MINIAPP             do k = ADM_kmin+1, ADM_kmax+1
!MINIAPP                if ( hflat < GRD_gzh(k) ) then
!MINIAPP                   kflat = k
!MINIAPP                   exit
!MINIAPP                endif
!MINIAPP             enddo
!MINIAPP          endif
!MINIAPP
!MINIAPP          if ( kflat == -1 ) then
!MINIAPP             kflat = ADM_kmax + 1
!MINIAPP             htop  = GRD_htop
!MINIAPP          else
!MINIAPP             htop  = GRD_gzh(kflat) - GRD_gzh(ADM_kmin)
!MINIAPP          endif
!MINIAPP
!MINIAPP          do l = 1, ADM_lall
!MINIAPP             do k = ADM_kmin-1, kflat
!MINIAPP             do n = nstart, nend
!MINIAPP                GRD_vz(n,k,l,GRD_Z ) = GRD_zs(n,k0,l,GRD_ZSFC) &
!MINIAPP                                     + ( htop - GRD_zs(n,k0,l,GRD_ZSFC) ) / htop * GRD_gz(k)
!MINIAPP                GRD_vz(n,k,l,GRD_ZH) = GRD_zs(n,k0,l,GRD_ZSFC) &
!MINIAPP                                     + ( htop - GRD_zs(n,k0,l,GRD_ZSFC) ) / htop * GRD_gzh(k)
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP
!MINIAPP             if ( kflat < ADM_kmax+1 ) then
!MINIAPP                do k = kflat+1, ADM_kmax+1
!MINIAPP                do n = nstart, nend
!MINIAPP                   GRD_vz(n,k,l,GRD_Z ) = GRD_gz (k)
!MINIAPP                   GRD_vz(n,k,l,GRD_ZH) = GRD_gzh(k)
!MINIAPP                enddo
!MINIAPP                enddo
!MINIAPP             endif
!MINIAPP          enddo
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             n = ADM_gslf_pl
!MINIAPP             do l = 1, ADM_lall_pl
!MINIAPP                do k = ADM_kmin-1, kflat
!MINIAPP                   GRD_vz_pl(n,k,l,GRD_Z)  = GRD_zs_pl(n,k0,l,GRD_ZSFC) &
!MINIAPP                                           + ( htop - GRD_zs_pl(n,k0,l,GRD_ZSFC) ) / htop * GRD_gz(k)
!MINIAPP                   GRD_vz_pl(n,k,l,GRD_ZH) = GRD_zs_pl(n,k0,l,GRD_ZSFC) &
!MINIAPP                                           + ( htop - GRD_zs_pl(n,k0,l,GRD_ZSFC) ) / htop * GRD_gzh(k)
!MINIAPP                enddo
!MINIAPP
!MINIAPP                if ( kflat < ADM_kmax+1 ) then
!MINIAPP                   do k = kflat+1, ADM_kmax+1
!MINIAPP                      GRD_vz_pl(n,k,l,GRD_Z ) = GRD_gz (k)
!MINIAPP                      GRD_vz_pl(n,k,l,GRD_ZH) = GRD_gzh(k)
!MINIAPP                   enddo
!MINIAPP                endif
!MINIAPP             enddo
!MINIAPP          endif
!MINIAPP
!MINIAPP       case('HYBRID')
!MINIAPP
!MINIAPP          ! Hybrid transformation : like as Simmons & Buridge(1981)
!MINIAPP
!MINIAPP          do l = 1, ADM_lall
!MINIAPP          do k = ADM_kmin-1, ADM_kmax+1
!MINIAPP          do n = nstart, nend
!MINIAPP             GRD_vz(n,k,l,GRD_Z)  = GRD_gz(k)                               &
!MINIAPP                                  + GRD_zs(n,k0,l,GRD_ZSFC)                 &
!MINIAPP                                  * sinh( (GRD_htop-GRD_gz (k)) / h_efold ) &
!MINIAPP                                  / sinh(  GRD_htop             / h_efold )
!MINIAPP             GRD_vz(n,k,l,GRD_ZH) = GRD_gzh(k)                              &
!MINIAPP                                  + GRD_zs(n,k0,l,GRD_ZSFC)                 &
!MINIAPP                                  * sinh( (GRD_htop-GRD_gzh(k)) / h_efold ) &
!MINIAPP                                  / sinh(  GRD_htop             / h_efold )
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP          enddo
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             n = ADM_gslf_pl
!MINIAPP             do l = 1, ADM_lall_pl
!MINIAPP             do k = ADM_kmin-1, ADM_kmax+1
!MINIAPP                GRD_vz_pl(n,k,l,GRD_Z)  = GRD_gz(k)                               &
!MINIAPP                                        + GRD_zs_pl(n,k0,l,GRD_ZSFC)              &
!MINIAPP                                        * sinh( (GRD_htop-GRD_gz (k)) / h_efold ) &
!MINIAPP                                        / sinh(  GRD_htop             / h_efold )
!MINIAPP                GRD_vz_pl(n,k,l,GRD_ZH) = GRD_gzh(k)                              &
!MINIAPP                                        + GRD_zs_pl(n,k0,l,GRD_ZSFC)              &
!MINIAPP                                        * sinh( (GRD_htop-GRD_gzh(k)) / h_efold ) &
!MINIAPP                                        / sinh(  GRD_htop             / h_efold )
!MINIAPP             enddo
!MINIAPP             enddo
!MINIAPP          endif
!MINIAPP
!MINIAPP       endselect
!MINIAPP
!MINIAPP       ! fill HALO
!MINIAPP       call COMM_data_transfer( GRD_vz, GRD_vz_pl )
!MINIAPP
!MINIAPP       allocate( zaxis(ADM_kall) )
!MINIAPP       zaxis(:) = real(GRD_gz(:),kind=DP)
!MINIAPP
!MINIAPP       call FIO_set_layer( "xir", "model layer (restart)", "", "m", ADM_kall, zaxis(:) )
!MINIAPP
!MINIAPP       deallocate( zaxis )
!MINIAPP    else
!MINIAPP       allocate( GRD_gz (ADM_KNONE) )
!MINIAPP       allocate( GRD_gzh(ADM_KNONE) )
!MINIAPP       GRD_gz (:) = 1.0_RP
!MINIAPP       GRD_gzh(:) = 1.0_RP
!MINIAPP
!MINIAPP       allocate( zaxis(1) )
!MINIAPP       zaxis(1) = 0.0_DP
!MINIAPP
!MINIAPP       call FIO_set_layer( "zsfc1", "model layer (restart)", "", "NIL", 1, zaxis(:) )
!MINIAPP
!MINIAPP       deallocate( zaxis )
!MINIAPP    endif
!MINIAPP
!MINIAPP    !--- output information about grid.
!MINIAPP    if ( ADM_kall /= ADM_KNONE ) then
!MINIAPP       write(IO_FID_LOG,*)
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|======      Vertical Coordinate [m]      ======|'
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|                                               |'
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|          -GRID CENTER-       -GRID INTERFACE- |'
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|  k        gz     d(gz)      gzh    d(gzh)   k |'
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|                                               |'
!MINIAPP       k = ADM_kmax + 1
!MINIAPP       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
!MINIAPP       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | TOA'
!MINIAPP       k = ADM_kmax
!MINIAPP       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmax'
!MINIAPP       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
!MINIAPP       do k = ADM_kmax-1, ADM_kmin+1, -1
!MINIAPP       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        |'
!MINIAPP       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' |'
!MINIAPP       enddo
!MINIAPP       k = ADM_kmin
!MINIAPP       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | kmin'
!MINIAPP       write(IO_FID_LOG,'(1x,A,2F10.1,I4,A)') '|                      ',GRD_gzh(k),GRD_dgzh(k),k,' | ground'
!MINIAPP       k = ADM_kmin-1
!MINIAPP       write(IO_FID_LOG,'(1x,A,I3,2F10.1,A)') '|',k,GRD_gz(k),GRD_dgz(k), '                        | dummy'
!MINIAPP       write(IO_FID_LOG,'(1x,A)')             '|===============================================|'
!MINIAPP       write(IO_FID_LOG,*) '--- Vertical layer scheme = ', trim(vgrid_scheme)
!MINIAPP       if ( vgrid_scheme == 'HYBRID' ) then
!MINIAPP          write(IO_FID_LOG,*) '--- e-folding height = ', h_efold
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( output_vgrid ) then
!MINIAPP          if ( ADM_IsMaster ) then
!MINIAPP             call GRD_output_vgrid('./vgrid_used.dat')
!MINIAPP          endif
!MINIAPP       endif
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,*)
!MINIAPP       write(IO_FID_LOG,*) '--- vartical layer = 1'
!MINIAPP    endif
!MINIAPP
!MINIAPP    !%%%%% Obsolete variables for backward compatibility %%%%%
!MINIAPP    if ( ADM_kall /= ADM_KNONE ) then
!MINIAPP       allocate( GRD_afac(ADM_kall) )
!MINIAPP       allocate( GRD_bfac(ADM_kall) )
!MINIAPP       allocate( GRD_cfac(ADM_kall) )
!MINIAPP       allocate( GRD_dfac(ADM_kall) )
!MINIAPP       GRD_afac(:) = 2.0_RP * GRD_afact(:)
!MINIAPP       GRD_bfac(:) = 2.0_RP * GRD_bfact(:)
!MINIAPP       GRD_cfac(:) = 2.0_RP * GRD_cfact(:)
!MINIAPP       GRD_dfac(:) = 2.0_RP * GRD_dfact(:)
!MINIAPP    endif
!MINIAPP    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_setup
  !--- MINIAPP: end of original GRD_setup -----------------------------------

  !-----------------------------------------------------------------------------
  !> Input horizontal grid
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_input_hgrid( &
!MINIAPP       basename,     &
!MINIAPP       input_vertex, &
!MINIAPP       io_mode,      &
!MINIAPP       direct_access )
!MINIAPP    use mod_adm, only: &
!MINIAPP       GLOBAL_prefix_dir,    &
!MINIAPP       GLOBAL_extension_ens, &
!MINIAPP       ADM_make_idstr,       &
!MINIAPP       ADM_proc_stop,        &
!MINIAPP       ADM_TI,               &
!MINIAPP       ADM_TJ,               &
!MINIAPP       ADM_have_pl,          &
!MINIAPP       ADM_lall,             &
!MINIAPP       ADM_lall_pl,          &
!MINIAPP       ADM_gall,             &
!MINIAPP       ADM_gall_pl,          &
!MINIAPP       ADM_gall_1d,          &
!MINIAPP       ADM_KNONE
!MINIAPP    use mod_fio_common, only: &
!MINIAPP       FIO_HSHORT, &
!MINIAPP       FIO_REAL8
!MINIAPP    use mod_fio, only: &
!MINIAPP       FIO_FORMAT,       &
!MINIAPP       FIO_input,        &
!MINIAPP       FIO_input_legacy, &
!MINIAPP       FIO_close_legacy
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    character(len=*), intent(in) :: basename                ! input basename
!MINIAPP    logical,          intent(in) :: input_vertex            ! flag of B-grid input
!MINIAPP    character(len=*), intent(in), optional :: io_mode       ! io_mode
!MINIAPP    logical,          intent(in), optional :: direct_access ! direct access? (only for LEGACY)
!MINIAPP
!MINIAPP    character(len=H_SHORT) :: io_mode_
!MINIAPP    logical                :: direct_access_
!MINIAPP    character(len=H_FILE)  :: fname
!MINIAPP
!MINIAPP    real(DP) :: x_DP    (ADM_gall   ,              ADM_nxyz)
!MINIAPP    real(DP) :: x_DP_pl (ADM_gall_pl,ADM_lall_pl,  ADM_nxyz)
!MINIAPP    real(DP) :: xt_DP   (ADM_gall   ,ADM_TI:ADM_TJ,ADM_nxyz)
!MINIAPP    real(DP) :: xt_DP_pl(ADM_gall_pl,ADM_lall_pl,  ADM_nxyz)
!MINIAPP
!MINIAPP    character(len=FIO_HSHORT) :: layername
!MINIAPP
!MINIAPP    integer  :: l
!MINIAPP    integer  :: fid, ierr
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,*) '*** [GRD_input_hgrid] Read horizaontal grid. file = ', trim(basename)
!MINIAPP
!MINIAPP    if ( present(io_mode) ) then
!MINIAPP       io_mode_ = io_mode
!MINIAPP    else
!MINIAPP       io_mode_ = hgrid_io_mode
!MINIAPP    endif
!MINIAPP
!MINIAPP    if ( present(direct_access) ) then
!MINIAPP       direct_access_ = direct_access
!MINIAPP    else
!MINIAPP       direct_access_ = hgrid_direct_access
!MINIAPP    endif
!MINIAPP
!MINIAPP    if ( io_mode_ == 'ADVANCED' ) then
!MINIAPP
!MINIAPP       if ( FIO_FORMAT == 'NETCDF' ) then
!MINIAPP          layername = 'zsfc1'
!MINIAPP       else
!MINIAPP          layername = 'ZSSFC1'
!MINIAPP       endif
!MINIAPP
!MINIAPP       call FIO_input( GRD_x(:,:,:,GRD_XDIR),basename,'grd_x_x',layername,1,1,1 )
!MINIAPP       call FIO_input( GRD_x(:,:,:,GRD_YDIR),basename,'grd_x_y',layername,1,1,1 )
!MINIAPP       call FIO_input( GRD_x(:,:,:,GRD_ZDIR),basename,'grd_x_z',layername,1,1,1 )
!MINIAPP       if ( input_vertex ) then
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),basename,'grd_xt_ix',layername,1,1,1 )
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),basename,'grd_xt_jx',layername,1,1,1 )
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),basename,'grd_xt_iy',layername,1,1,1 )
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),basename,'grd_xt_jy',layername,1,1,1 )
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),basename,'grd_xt_iz',layername,1,1,1 )
!MINIAPP          call FIO_input( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),basename,'grd_xt_jz',layername,1,1,1 )
!MINIAPP       endif
!MINIAPP
!MINIAPP       call GRD_gen_plgrid
!MINIAPP
!MINIAPP    elseif( io_mode_ == 'LEGACY' ) then
!MINIAPP
!MINIAPP       if ( direct_access_ ) then
!MINIAPP
!MINIAPP          call FIO_input_legacy( GRD_x(:,:,:,GRD_XDIR),basename,'grd_x_x',1,1,irec=1,dtype=FIO_REAL8 )
!MINIAPP          call FIO_input_legacy( GRD_x(:,:,:,GRD_YDIR),basename,'grd_x_y',1,1,irec=2,dtype=FIO_REAL8 )
!MINIAPP          call FIO_input_legacy( GRD_x(:,:,:,GRD_ZDIR),basename,'grd_x_z',1,1,irec=3,dtype=FIO_REAL8 )
!MINIAPP          if (input_vertex) then
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),basename,'grd_xt_ix',1,1,irec=4,dtype=FIO_REAL8 )
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),basename,'grd_xt_iy',1,1,irec=5,dtype=FIO_REAL8 )
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),basename,'grd_xt_iz',1,1,irec=6,dtype=FIO_REAL8 )
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),basename,'grd_xt_jx',1,1,irec=7,dtype=FIO_REAL8 )
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),basename,'grd_xt_jy',1,1,irec=8,dtype=FIO_REAL8 )
!MINIAPP             call FIO_input_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),basename,'grd_xt_jz',1,1,irec=9,dtype=FIO_REAL8 )
!MINIAPP          endif
!MINIAPP          call FIO_close_legacy(basename)
!MINIAPP
!MINIAPP          call GRD_gen_plgrid
!MINIAPP
!MINIAPP       else
!MINIAPP
!MINIAPP          do l = 1, ADM_lall
!MINIAPP             call ADM_make_idstr(fname,basename,l)
!MINIAPP             fid = IO_get_available_fid()
!MINIAPP
!MINIAPP             open( unit   = fid,           &
!MINIAPP                   file   = trim(fname),   &
!MINIAPP                   status = 'old',         &
!MINIAPP                   form   = 'unformatted', &
!MINIAPP                   access = 'sequential',  &
!MINIAPP                   iostat = ierr           )
!MINIAPP
!MINIAPP                if ( ierr /= 0 ) then
!MINIAPP                   write(IO_FID_LOG,*) 'xxx [grid] No input grid file. ', trim(fname)
!MINIAPP                   call ADM_proc_stop
!MINIAPP                endif
!MINIAPP
!MINIAPP                read(fid) ADM_gall_1d
!MINIAPP                read(fid) x_DP(:,GRD_XDIR)
!MINIAPP                read(fid) x_DP(:,GRD_YDIR)
!MINIAPP                read(fid) x_DP(:,GRD_ZDIR)
!MINIAPP                GRD_x(:,ADM_KNONE,l,GRD_XDIR) = x_DP(:,GRD_XDIR)
!MINIAPP                GRD_x(:,ADM_KNONE,l,GRD_YDIR) = x_DP(:,GRD_YDIR)
!MINIAPP                GRD_x(:,ADM_KNONE,l,GRD_ZDIR) = x_DP(:,GRD_ZDIR)
!MINIAPP
!MINIAPP                if (input_vertex) then
!MINIAPP                   read(fid) xt_DP(:,:,GRD_XDIR)
!MINIAPP                   read(fid) xt_DP(:,:,GRD_YDIR)
!MINIAPP                   read(fid) xt_DP(:,:,GRD_ZDIR)
!MINIAPP                   GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_XDIR) = xt_DP(:,:,GRD_XDIR)
!MINIAPP                   GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_YDIR) = xt_DP(:,:,GRD_YDIR)
!MINIAPP                   GRD_xt(:,ADM_KNONE,l,ADM_TI:ADM_TJ,GRD_ZDIR) = xt_DP(:,:,GRD_ZDIR)
!MINIAPP                endif
!MINIAPP             close(fid)
!MINIAPP          enddo
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             fname = trim(GLOBAL_prefix_dir)//trim(basename)//trim(GLOBAL_extension_ens)//'.pl'
!MINIAPP             fid = IO_get_available_fid()
!MINIAPP
!MINIAPP             open( unit   = fid,           &
!MINIAPP                   file   = trim(fname),   &
!MINIAPP                   status = 'old',         &
!MINIAPP                   form   = 'unformatted', &
!MINIAPP                   access = 'sequential'   )
!MINIAPP
!MINIAPP                read(fid) x_DP_pl(:,:,:)
!MINIAPP                GRD_x_pl(:,ADM_KNONE,:,:) = x_DP_pl(:,:,:)
!MINIAPP
!MINIAPP                if (input_vertex) then
!MINIAPP                   read(fid) xt_DP_pl(:,:,:)
!MINIAPP                   GRD_xt_pl(:,ADM_KNONE,:,:) = xt_DP_pl(:,:,:)
!MINIAPP                endif
!MINIAPP             close(fid)
!MINIAPP          endif
!MINIAPP
!MINIAPP       endif
!MINIAPP
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,*) 'xxx [grd/GRD_input_hgrid] Invalid io_mode!', trim(io_mode_)
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_input_hgrid

  !-----------------------------------------------------------------------------
  !> Output horizontal grid
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_output_hgrid( &
!MINIAPP       basename,      &
!MINIAPP       output_vertex, &
!MINIAPP       io_mode,       &
!MINIAPP       direct_access  )
!MINIAPP    use mod_adm, only: &
!MINIAPP       GLOBAL_prefix_dir,    &
!MINIAPP       GLOBAL_extension_ens, &
!MINIAPP       ADM_make_idstr,       &
!MINIAPP       ADM_proc_stop,        &
!MINIAPP       ADM_TI,               &
!MINIAPP       ADM_TJ,               &
!MINIAPP       ADM_have_pl,          &
!MINIAPP       ADM_lall,             &
!MINIAPP       ADM_gall_1d
!MINIAPP    use mod_fio_common, only: &
!MINIAPP       FIO_REAL8
!MINIAPP    use mod_fio, only: &
!MINIAPP       FIO_FORMAT,        &
!MINIAPP       FIO_set_layer,     &
!MINIAPP       FIO_output,        &
!MINIAPP       FIO_output_legacy, &
!MINIAPP       FIO_close_legacy
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    character(len=*), intent(in) :: basename                ! output basename
!MINIAPP    logical,          intent(in) :: output_vertex           ! output flag of B-grid
!MINIAPP    character(len=*), intent(in), optional :: io_mode       ! io_mode
!MINIAPP    logical,          intent(in), optional :: direct_access ! direct access? (only for LEGACY)
!MINIAPP
!MINIAPP    character(len=H_SHORT) :: io_mode_
!MINIAPP    logical                :: direct_access_
!MINIAPP    character(len=H_FILE)  :: fname
!MINIAPP    character(len=H_MID)   :: desc = 'HORIZONTAL GRID FILE'
!MINIAPP    character(len=H_SHORT) :: layername
!MINIAPP    real(DP)               :: zaxis(1) = 0.0_DP
!MINIAPP
!MINIAPP    integer  :: fid
!MINIAPP    integer  :: l
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    if ( present(io_mode) ) then
!MINIAPP       io_mode_ = io_mode
!MINIAPP    else
!MINIAPP       io_mode_ = hgrid_io_mode
!MINIAPP    endif
!MINIAPP
!MINIAPP    if ( present(direct_access) ) then
!MINIAPP       direct_access_ = direct_access
!MINIAPP    else
!MINIAPP       direct_access_ = hgrid_direct_access
!MINIAPP    endif
!MINIAPP
!MINIAPP    if ( io_mode_ == 'ADVANCED' ) then
!MINIAPP
!MINIAPP       if ( FIO_FORMAT == 'NETCDF' ) then
!MINIAPP          layername = "zsfc1"
!MINIAPP          call FIO_set_layer( "zsfc1", "zaxis for surface variable", "", "NIL", 1, zaxis(:) )
!MINIAPP       else
!MINIAPP          layername = "ZSSFC1"
!MINIAPP       endif
!MINIAPP
!MINIAPP       call FIO_output( GRD_x(:,:,:,GRD_XDIR),                              & ! [IN]
!MINIAPP                        basename, desc, '',                                 & ! [IN]
!MINIAPP                       'grd_x_x', 'GRD_x (X_DIR)', '',                      & ! [IN]
!MINIAPP                       'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP       call FIO_output( GRD_x(:,:,:,GRD_YDIR),                              & ! [IN]
!MINIAPP                        basename, desc, '',                                 & ! [IN]
!MINIAPP                       'grd_x_y', 'GRD_x (Y_DIR)', '',                      & ! [IN]
!MINIAPP                       'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP       call FIO_output( GRD_x(:,:,:,GRD_ZDIR),                              & ! [IN]
!MINIAPP                        basename, desc, '',                                 & ! [IN]
!MINIAPP                       'grd_x_z', 'GRD_x (Z_DIR)', '',                      & ! [IN]
!MINIAPP                       'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP
!MINIAPP       if ( output_vertex ) then
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_ix', 'GRD_xt (TI,X_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_jx', 'GRD_xt (TJ,X_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_iy', 'GRD_xt (TI,Y_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_jy', 'GRD_xt (TJ,Y_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_iz', 'GRD_xt (TI,Z_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP          call FIO_output( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),                      & ! [IN]
!MINIAPP                           basename, desc, '',                                 & ! [IN]
!MINIAPP                          'grd_xt_jz', 'GRD_xt (TJ,Z_DIR)', '',                & ! [IN]
!MINIAPP                          'NIL', FIO_REAL8, layername, 1, 1, 1, 0.0_DP, 0.0_DP ) ! [IN]
!MINIAPP       endif
!MINIAPP
!MINIAPP    elseif( io_mode_ == 'LEGACY' ) then
!MINIAPP
!MINIAPP       if ( direct_access_ ) then
!MINIAPP
!MINIAPP          call FIO_output_legacy( GRD_x(:,:,:,GRD_XDIR),basename,'grd_x_x',1,1,irec=1,dtype=FIO_REAL8 )
!MINIAPP          call FIO_output_legacy( GRD_x(:,:,:,GRD_YDIR),basename,'grd_x_y',1,1,irec=2,dtype=FIO_REAL8 )
!MINIAPP          call FIO_output_legacy( GRD_x(:,:,:,GRD_ZDIR),basename,'grd_x_z',1,1,irec=3,dtype=FIO_REAL8 )
!MINIAPP
!MINIAPP          if (output_vertex) then
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TI,GRD_XDIR),basename,'grd_xt_ix',1,1,irec=4,dtype=FIO_REAL8 )
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TI,GRD_YDIR),basename,'grd_xt_iy',1,1,irec=5,dtype=FIO_REAL8 )
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TI,GRD_ZDIR),basename,'grd_xt_iz',1,1,irec=6,dtype=FIO_REAL8 )
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_XDIR),basename,'grd_xt_jx',1,1,irec=7,dtype=FIO_REAL8 )
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_YDIR),basename,'grd_xt_jy',1,1,irec=8,dtype=FIO_REAL8 )
!MINIAPP             call FIO_output_legacy( GRD_xt(:,:,:,ADM_TJ,GRD_ZDIR),basename,'grd_xt_jz',1,1,irec=9,dtype=FIO_REAL8 )
!MINIAPP          endif
!MINIAPP          call FIO_close_legacy(basename)
!MINIAPP
!MINIAPP       else ! sequential, to be abolished
!MINIAPP
!MINIAPP          do l = 1, ADM_lall
!MINIAPP             call ADM_make_idstr(fname,basename,l)
!MINIAPP             fid = IO_get_available_fid()
!MINIAPP
!MINIAPP             open( unit   = fid,           &
!MINIAPP                   file   = trim(fname),   &
!MINIAPP                   form   = 'unformatted', &
!MINIAPP                   access = 'sequential'   )
!MINIAPP                write(fid) ADM_gall_1d
!MINIAPP                write(fid) real(GRD_x(:,:,l,GRD_XDIR),kind=DP)
!MINIAPP                write(fid) real(GRD_x(:,:,l,GRD_YDIR),kind=DP)
!MINIAPP                write(fid) real(GRD_x(:,:,l,GRD_ZDIR),kind=DP)
!MINIAPP                if (output_vertex) then
!MINIAPP                   write(fid) real(GRD_xt(:,:,l,:,GRD_XDIR),kind=DP)
!MINIAPP                   write(fid) real(GRD_xt(:,:,l,:,GRD_YDIR),kind=DP)
!MINIAPP                   write(fid) real(GRD_xt(:,:,l,:,GRD_ZDIR),kind=DP)
!MINIAPP                endif
!MINIAPP             close(fid)
!MINIAPP
!MINIAPP          enddo
!MINIAPP
!MINIAPP          if ( ADM_have_pl ) then
!MINIAPP             fname = trim(GLOBAL_prefix_dir)//trim(basename)//trim(GLOBAL_extension_ens)//'.pl'
!MINIAPP             fid = IO_get_available_fid()
!MINIAPP
!MINIAPP             open( unit   = fid,           &
!MINIAPP                   file   = trim(fname),   &
!MINIAPP                   form   = 'unformatted', &
!MINIAPP                   access = 'sequential'   )
!MINIAPP                write(fid) real(GRD_x_pl(:,:,:,:),kind=DP)
!MINIAPP                if (output_vertex) then
!MINIAPP                   write(fid) real(GRD_xt_pl(:,:,:,:),kind=DP)
!MINIAPP                endif
!MINIAPP             close(fid)
!MINIAPP          endif
!MINIAPP
!MINIAPP       endif
!MINIAPP
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,*) 'xxx [grd/GRD_output_hgrid] Invalid io_mode!', trim(io_mode_)
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_output_hgrid

  !-----------------------------------------------------------------------------
  !> Input vertical grid
  subroutine GRD_input_vgrid_SP( &
       fname, &
       kall,  &
       htop,  &
       gz,    &
       gzh    )
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    character(len=*), intent(in)  :: fname ! vertical grid file name
    integer,          intent(in)  :: kall
    real(SP),         intent(out) :: htop
    real(SP),         intent(out) :: gz (kall)
    real(SP),         intent(out) :: gzh(kall)

    integer  :: num_of_layer
    real(DP) :: gz_DP (kall)
    real(DP) :: gzh_DP(kall)

    integer  :: fid, ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** [GRD_input_vgrid_SP] Read vertical grid. file = ', trim(fname)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'old',         &
          form   = 'unformatted', &
          iostat = ierr           )

       if ( ierr /= 0 ) then
          write(IO_FID_LOG,*) 'xxx [grd/GRD_input_vgrid] No vertical grid file.'
          call ADM_proc_stop
       endif

       read(fid) num_of_layer

       if ( num_of_layer+2 /= kall ) then
          write(IO_FID_LOG,*) 'xxx [grd/GRD_input_vgrid] inconsistency in number of vertical layers.'
          call ADM_proc_stop
       endif

       read(fid) gz_DP (:)
       read(fid) gzh_DP(:)

    close(fid)

    gz (:) = real(gz_DP ,kind=SP)
    gzh(:) = real(gzh_DP,kind=SP)

    ! hight top
    htop = gzh(kall)

    return
  end subroutine GRD_input_vgrid_SP

  !-----------------------------------------------------------------------------
  !> Input vertical grid
  subroutine GRD_input_vgrid_DP( &
       fname, &
       kall,  &
       htop,  &
       gz,    &
       gzh    )
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    character(len=*), intent(in)  :: fname ! vertical grid file name
    integer,          intent(in)  :: kall
    real(DP),         intent(out) :: htop
    real(DP),         intent(out) :: gz (kall)
    real(DP),         intent(out) :: gzh(kall)

    integer  :: num_of_layer
    real(DP) :: gz_DP (kall)
    real(DP) :: gzh_DP(kall)

    integer  :: fid, ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** [GRD_input_vgrid_DP] Read vertical grid. file = ', trim(fname)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'old',         &
          form   = 'unformatted', &
          iostat = ierr           )

       if ( ierr /= 0 ) then
          write(IO_FID_LOG,*) 'xxx [grd/GRD_input_vgrid] No vertical grid file.'
          call ADM_proc_stop
       endif

       read(fid) num_of_layer

       if ( num_of_layer+2 /= kall ) then
          write(IO_FID_LOG,*) 'xxx [grd/GRD_input_vgrid] inconsistency in number of vertical layers.'
          call ADM_proc_stop
       endif

       read(fid) gz_DP (:)
       read(fid) gzh_DP(:)

    close(fid)

    gz (:) = real(gz_DP ,kind=DP)
    gzh(:) = real(gzh_DP,kind=DP)

    ! hight top
    htop = gzh(kall)

    return
  end subroutine GRD_input_vgrid_DP

  !-----------------------------------------------------------------------------
  !> Output vertical grid
  subroutine GRD_output_vgrid( &
       basename )
    use mod_adm, only: &
       GLOBAL_prefix_dir,    &
       GLOBAL_extension_ens, &
       ADM_vlayer
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_FILE) :: fname

    integer  :: fid, ierr
    !---------------------------------------------------------------------------

    fname = trim(GLOBAL_prefix_dir)//trim(basename)//trim(GLOBAL_extension_ens)

    write(IO_FID_LOG,*) '*** Write vertical grid file: ', trim(fname)

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'new',         &
          form   = 'unformatted', &
          iostat = ierr           )

       write(fid) ADM_vlayer
       write(fid) real(GRD_gz ,kind=DP)
       write(fid) real(GRD_gzh,kind=DP)

    close(fid)

    return
  end subroutine GRD_output_vgrid

  !-----------------------------------------------------------------------------
  !> Input topography data
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_input_topograph( &
!MINIAPP       topo_basename,   &
!MINIAPP       toposd_basename, &
!MINIAPP       io_mode          )
!MINIAPP    use mod_adm, only: &
!MINIAPP       GLOBAL_prefix_dir,    &
!MINIAPP       GLOBAL_extension_ens, &
!MINIAPP       ADM_make_idstr,       &
!MINIAPP       ADM_proc_stop,        &
!MINIAPP       ADM_have_pl,          &
!MINIAPP       ADM_lall,             &
!MINIAPP       ADM_lall_pl,          &
!MINIAPP       ADM_gall,             &
!MINIAPP       ADM_gall_pl,          &
!MINIAPP       ADM_KNONE
!MINIAPP    use mod_fio_common, only: &
!MINIAPP       FIO_HSHORT, &
!MINIAPP       FIO_REAL8
!MINIAPP    use mod_fio, only: &
!MINIAPP       FIO_FORMAT,       &
!MINIAPP       FIO_input,        &
!MINIAPP       FIO_input_legacy, &
!MINIAPP       FIO_close_legacy
!MINIAPP    use mod_comm, only: &
!MINIAPP       COMM_data_transfer, &
!MINIAPP       COMM_var
!MINIAPP!TENTATIVE    use mod_ideal_topo, only: &
!MINIAPP!TENTATIVE       IDEAL_topo
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    character(len=*), intent(in) :: topo_basename   ! input basename (topography)
!MINIAPP    character(len=*), intent(in) :: toposd_basename ! input basename (std.dev. of topography)
!MINIAPP    character(len=*), intent(in) :: io_mode         ! io_mode
!MINIAPP
!MINIAPP    character(len=H_FILE) :: fname
!MINIAPP
!MINIAPP    real(DP) :: val_DP   (ADM_gall)
!MINIAPP    real(DP) :: val_DP_pl(ADM_gall_pl,ADM_lall_pl)
!MINIAPP
!MINIAPP    character(len=FIO_HSHORT) :: layername
!MINIAPP
!MINIAPP    integer  :: fid, ierr
!MINIAPP    integer  :: l
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,*) '*** [GRD_input_topograph] Read topography data. file = ', trim(topo_basename)
!MINIAPP
!MINIAPP    if ( io_mode == 'ADVANCED' ) then
!MINIAPP
!MINIAPP       if ( FIO_FORMAT == 'NETCDF' ) then
!MINIAPP          layername = 'zsfc1'
!MINIAPP       else
!MINIAPP          layername = 'ZSSFC1'
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( topo_basename /= 'NONE' ) then
!MINIAPP          call FIO_input(GRD_zs(:,:,:,GRD_ZSFC),topo_basename,  'topo'       ,layername,1,1,1)
!MINIAPP       endif
!MINIAPP       if ( toposd_basename /= 'NONE' ) then
!MINIAPP          call FIO_input(GRD_zs(:,:,:,GRD_ZSD ),toposd_basename,'topo_stddev',layername,1,1,1)
!MINIAPP       endif
!MINIAPP
!MINIAPP    elseif( topo_io_mode == 'LEGACY' ) then
!MINIAPP
!MINIAPP       if (topo_direct_access) then
!MINIAPP
!MINIAPP          if ( topo_basename /= 'NONE' .AND. topo_basename /= '' ) then
!MINIAPP             call FIO_input_legacy( GRD_zs(:,:,:,GRD_ZSFC),topo_basename,  'topo'       ,1,1,irec=1,dtype=FIO_REAL8 )
!MINIAPP             call FIO_close_legacy(topo_basename)
!MINIAPP          endif
!MINIAPP          if ( toposd_basename /= 'NONE' .AND. toposd_basename /= '' ) then
!MINIAPP             call FIO_input_legacy( GRD_zs(:,:,:,GRD_ZSD ),toposd_basename,'topo_stddev',1,1,irec=1,dtype=FIO_REAL8 )
!MINIAPP             call FIO_close_legacy(toposd_basename)
!MINIAPP          endif
!MINIAPP
!MINIAPP       else !--- sequential access
!MINIAPP
!MINIAPP          if ( topo_basename /= 'NONE' ) then
!MINIAPP             do l = 1, ADM_lall
!MINIAPP                call ADM_make_idstr(fname,trim(topo_basename),l)
!MINIAPP                fid = IO_get_available_fid()
!MINIAPP
!MINIAPP                open( unit   = fid,           &
!MINIAPP                      file   = trim(fname),   &
!MINIAPP                      status = 'old',         &
!MINIAPP                      form   = 'unformatted', &
!MINIAPP                      access = 'sequential',  &
!MINIAPP                      iostat = ierr           )
!MINIAPP
!MINIAPP                   if ( ierr /= 0 ) then
!MINIAPP                      write(IO_FID_LOG,*) 'xxx [grid] No topography file.'
!MINIAPP                      return
!MINIAPP                   endif
!MINIAPP
!MINIAPP                   read(fid) val_DP(:)
!MINIAPP                   GRD_zs(:,ADM_KNONE,l,GRD_ZSFC) = val_DP(:)
!MINIAPP
!MINIAPP                close(fid)
!MINIAPP             enddo
!MINIAPP
!MINIAPP             if ( ADM_have_pl ) then
!MINIAPP                fname = trim(GLOBAL_prefix_dir)//trim(topo_basename)//trim(GLOBAL_extension_ens)//'.pl'
!MINIAPP                fid = IO_get_available_fid()
!MINIAPP
!MINIAPP                open( unit   = fid,           &
!MINIAPP                      file   = trim(fname),   &
!MINIAPP                      status = 'old',         &
!MINIAPP                      form   = 'unformatted', &
!MINIAPP                      access = 'sequential'   )
!MINIAPP
!MINIAPP                    read(fid) val_DP_pl(:,:)
!MINIAPP                    GRD_zs_pl(:,ADM_KNONE,:,GRD_ZSFC) = val_DP_pl(:,:)
!MINIAPP                close(fid)
!MINIAPP             endif
!MINIAPP          endif
!MINIAPP
!MINIAPP          if ( toposd_basename /= 'NONE' ) then
!MINIAPP             do l = 1, ADM_lall
!MINIAPP                call ADM_make_idstr(fname,trim(toposd_basename),l)
!MINIAPP                fid = IO_get_available_fid()
!MINIAPP
!MINIAPP                open( unit   = fid,           &
!MINIAPP                      file   = trim(fname),   &
!MINIAPP                      status = 'old',         &
!MINIAPP                      form   = 'unformatted', &
!MINIAPP                      access = 'sequential',  &
!MINIAPP                      iostat = ierr           )
!MINIAPP
!MINIAPP                   if ( ierr /= 0 ) then
!MINIAPP                      write(IO_FID_LOG,*) 'xxx [grid] No topography (stdev) file.'
!MINIAPP                      return
!MINIAPP                   endif
!MINIAPP
!MINIAPP                   read(fid) val_DP(:)
!MINIAPP                   GRD_zs(:,ADM_KNONE,l,GRD_ZSD) = val_DP(:)
!MINIAPP
!MINIAPP                close(fid)
!MINIAPP             enddo
!MINIAPP
!MINIAPP             if ( ADM_have_pl ) then
!MINIAPP                fname = trim(GLOBAL_prefix_dir)//trim(toposd_basename)//trim(GLOBAL_extension_ens)//'.pl'
!MINIAPP                fid = IO_get_available_fid()
!MINIAPP
!MINIAPP                open( unit   = fid,           &
!MINIAPP                      file   = trim(fname),   &
!MINIAPP                      status = 'old',         &
!MINIAPP                      form   = 'unformatted', &
!MINIAPP                      access = 'sequential'   )
!MINIAPP
!MINIAPP                    read(fid) val_DP_pl(:,:)
!MINIAPP                    GRD_zs_pl(:,ADM_KNONE,:,GRD_ZSD) = val_DP_pl(:,:)
!MINIAPP
!MINIAPP                close(fid)
!MINIAPP             endif
!MINIAPP          endif
!MINIAPP
!MINIAPP          call COMM_data_transfer(GRD_zs,GRD_zs_pl)
!MINIAPP
!MINIAPP       endif !--- direct/sequencial
!MINIAPP
!MINIAPP!TENTATIVE    elseif( io_mode == 'IDEAL' ) then
!MINIAPP!TENTATIVE
!MINIAPP!TENTATIVE       write(IO_FID_LOG,*) '*** make ideal topography'
!MINIAPP!TENTATIVE
!MINIAPP!TENTATIVE       call IDEAL_topo( GRD_s (:,:,:,I_LAT), & ! [IN]
!MINIAPP!TENTATIVE                        GRD_s (:,:,:,I_LON), & ! [IN]
!MINIAPP!TENTATIVE                        GRD_zs(:,:,:,GRD_ZSFC) ) ! [OUT]
!MINIAPP!TENTATIVE
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,*) 'xxx [grd/GRD_input_topograph] Invalid io_mode!', trim(io_mode)
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif ! io_mode
!MINIAPP
!MINIAPP    call COMM_var( GRD_zs, GRD_zs_pl, ADM_KNONE, 2 )
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_input_topograph

  !-----------------------------------------------------------------------------
  !> Communicate grid data for pole region: This routine is NOT same as COMM_var
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_gen_plgrid
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_COMM_WORLD,     &
!MINIAPP       ADM_proc_stop,      &
!MINIAPP       I_N,                &
!MINIAPP       I_S,                &
!MINIAPP       I_NPL,              &
!MINIAPP       I_SPL,              &
!MINIAPP       I_RGNID,            &
!MINIAPP       I_prc,              &
!MINIAPP       RGNMNG_vert_tab_pl, &
!MINIAPP       RGNMNG_lp2r,        &
!MINIAPP       RGNMNG_r2lp,        &
!MINIAPP       RGNMNG_r2p_pl,      &
!MINIAPP       ADM_vlink,          &
!MINIAPP       ADM_prc_me,         &
!MINIAPP       ADM_gall_1d,        &
!MINIAPP       ADM_gmax,           &
!MINIAPP       ADM_gmin,           &
!MINIAPP       ADM_gslf_pl
!MINIAPP    use mod_comm, only: &
!MINIAPP       COMM_var
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer  :: prctab   (ADM_vlink)
!MINIAPP    integer  :: rgntab   (ADM_vlink)
!MINIAPP    integer  :: sreq     (ADM_vlink)
!MINIAPP    integer  :: rreq     (ADM_vlink)
!MINIAPP    logical  :: send_flag(ADM_vlink)
!MINIAPP
!MINIAPP    real(RP) :: vsend_pl (ADM_nxyz,ADM_vlink)
!MINIAPP    real(RP) :: vrecv_pl (ADM_nxyz,ADM_vlink)
!MINIAPP
!MINIAPP    integer  :: datatype
!MINIAPP
!MINIAPP    integer  :: istat(MPI_STATUS_SIZE)
!MINIAPP    integer  :: n, nn, l, ierr
!MINIAPP
!MINIAPP    integer  :: suf, i, j
!MINIAPP    suf(i,j) = ADM_gall_1d * (j-1) + i
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    if( GRD_grid_type == GRD_grid_type_on_plane ) return
!MINIAPP
!MINIAPP    if ( RP == DP ) then
!MINIAPP       datatype = MPI_DOUBLE_PRECISION
!MINIAPP    elseif( RP == SP ) then
!MINIAPP       datatype = MPI_REAL
!MINIAPP    else
!MINIAPP       write(*,*) 'xxx [grd/GRD_gen_plgrid] unsupportd precision! RP=', RP
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP
!MINIAPP    !--- send information of grid around north pole from regular region
!MINIAPP
!MINIAPP    ! find region which has the north pole
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       nn = n + 1
!MINIAPP       if( nn > ADM_vlink ) nn = 1
!MINIAPP
!MINIAPP       rgntab(n) = RGNMNG_vert_tab_pl(I_RGNID,I_NPL,nn)
!MINIAPP       prctab(n) = RGNMNG_r2lp(I_prc,rgntab(n))
!MINIAPP    enddo
!MINIAPP
!MINIAPP    call MPI_Barrier(ADM_COMM_WORLD,ierr)
!MINIAPP
!MINIAPP    send_flag(:) = .false.
!MINIAPP
!MINIAPP    ! send grid position from regular region
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       do l = 1, ADM_lall
!MINIAPP          if ( RGNMNG_lp2r(l,ADM_prc_me) == rgntab(n) ) then
!MINIAPP             vsend_pl(:,n) = GRD_xt(suf(ADM_gmin,ADM_gmax),ADM_KNONE,l,ADM_TJ,:) ! [mod] H.Yashiro 20120525
!MINIAPP
!MINIAPP             call MPI_ISEND( vsend_pl(:,n),          &
!MINIAPP                             3,                      &
!MINIAPP                             datatype,               &
!MINIAPP                             RGNMNG_r2p_pl(I_NPL)-1, &
!MINIAPP                             rgntab(n),              &
!MINIAPP                             ADM_COMM_WORLD,         &
!MINIAPP                             sreq(n),                &
!MINIAPP                             ierr                    )
!MINIAPP
!MINIAPP             send_flag(n) = .true.
!MINIAPP          endif
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_prc_me == RGNMNG_r2p_pl(I_NPL) ) then
!MINIAPP       do n = 1, ADM_vlink
!MINIAPP          call MPI_IRECV( vrecv_pl(:,n),  &
!MINIAPP                          3,              &
!MINIAPP                          datatype,       &
!MINIAPP                          prctab(n)-1,    &
!MINIAPP                          rgntab(n),      &
!MINIAPP                          ADM_COMM_WORLD, &
!MINIAPP                          rreq(n),        &
!MINIAPP                          ierr            )
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! wait and store
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       if ( send_flag(n) ) then
!MINIAPP          call MPI_WAIT(sreq(n),istat,ierr)
!MINIAPP       endif
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_prc_me == RGNMNG_r2p_pl(I_NPL) ) then
!MINIAPP       do n = 1, ADM_vlink
!MINIAPP          call MPI_WAIT(rreq(n),istat,ierr)
!MINIAPP          GRD_xt_pl(n+1,ADM_KNONE,I_NPL,:) = vrecv_pl(:,n) ! [mod] H.Yashiro 20120525
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    !--- send information of grid around south pole from regular region
!MINIAPP
!MINIAPP    ! find region which has the south pole
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       rgntab(n) = RGNMNG_vert_tab_pl(I_RGNID,I_SPL,n)
!MINIAPP       prctab(n) = RGNMNG_r2lp(I_prc,rgntab(n))
!MINIAPP    enddo
!MINIAPP
!MINIAPP    call MPI_Barrier(ADM_COMM_WORLD,ierr)
!MINIAPP
!MINIAPP    send_flag(:) = .false.
!MINIAPP
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       do l =1, ADM_lall
!MINIAPP          if (RGNMNG_lp2r(l,ADM_prc_me) == rgntab(n) ) then
!MINIAPP             vsend_pl(:,n) = GRD_xt(suf(ADM_gmax,ADM_gmin),ADM_KNONE,l,ADM_TI,:) ! [mod] H.Yashiro 20120525
!MINIAPP
!MINIAPP             call MPI_ISEND( vsend_pl(:,n),          &
!MINIAPP                             3,                      &
!MINIAPP                             datatype,               &
!MINIAPP                             RGNMNG_r2p_pl(I_SPL)-1, &
!MINIAPP                             rgntab(n),              &
!MINIAPP                             ADM_COMM_WORLD,         &
!MINIAPP                             sreq(n),                &
!MINIAPP                             ierr                    )
!MINIAPP
!MINIAPP             send_flag(n) = .true.
!MINIAPP          endif
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_prc_me == RGNMNG_r2p_pl(I_SPL) ) then
!MINIAPP       do n = 1, ADM_vlink
!MINIAPP          call MPI_IRECV( vrecv_pl(:,n),  &
!MINIAPP                          3,              &
!MINIAPP                          datatype,       &
!MINIAPP                          prctab(n)-1,    &
!MINIAPP                          rgntab(n),      &
!MINIAPP                          ADM_COMM_WORLD, &
!MINIAPP                          rreq(n),        &
!MINIAPP                          ierr            )
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! wait and store
!MINIAPP    do n = 1, ADM_vlink
!MINIAPP       if ( send_flag(n) ) then
!MINIAPP          call MPI_WAIT(sreq(n),istat,ierr)
!MINIAPP       endif
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_prc_me == RGNMNG_r2p_pl(I_SPL) ) then
!MINIAPP       do n = 1, ADM_vlink
!MINIAPP          call MPI_WAIT(rreq(n),istat,ierr)
!MINIAPP          GRD_xt_pl(n+1,ADM_KNONE,I_SPL,:) = vrecv_pl(:,n) ! [mod] H.Yashiro 20120525
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    ! grid point communication
!MINIAPP    call COMM_var( GRD_x, GRD_x_pl, ADM_KNONE, 3 )
!MINIAPP
!MINIAPP    if (      ADM_prc_me == RGNMNG_r2p_pl(I_NPL) &
!MINIAPP         .OR. ADM_prc_me == RGNMNG_r2p_pl(I_SPL) ) then
!MINIAPP       GRD_xt_pl(ADM_gslf_pl,:,:,:) = GRD_x_pl(ADM_gslf_pl,:,:,:)
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_gen_plgrid

  !-----------------------------------------------------------------------------
  !> scaling grid position on the plane/sphere
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_scaling( fact )
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(RP), intent(in) :: fact ! scaling factor
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    GRD_x (:,:,:,:)   = GRD_x (:,:,:,:)   * fact
!MINIAPP    GRD_xt(:,:,:,:,:) = GRD_xt(:,:,:,:,:) * fact
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       GRD_x_pl (:,:,:,:)   = GRD_x_pl (:,:,:,:) * fact
!MINIAPP       GRD_xt_pl(:,:,:,:)   = GRD_xt_pl(:,:,:,:) * fact
!MINIAPP    endif
!MINIAPP
!MINIAPP    if ( GRD_grid_type == GRD_grid_type_on_plane ) then
!MINIAPP       ! do nothing
!MINIAPP    else
!MINIAPP       ! setting the sphere radius
!MINIAPP       GRD_rscale = fact
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_scaling

  !-----------------------------------------------------------------------------
  !> calculate longitude and latitude
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_makelatlon
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl
!MINIAPP    use mod_vector, only: &
!MINIAPP       VECTR_xyz2latlon
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer  :: g, k0, l
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    do l = 1, ADM_lall
!MINIAPP    do g = 1, ADM_gall
!MINIAPP       call VECTR_xyz2latlon( GRD_x(g,k0,l,GRD_XDIR), & ! [IN]
!MINIAPP                              GRD_x(g,k0,l,GRD_YDIR), & ! [IN]
!MINIAPP                              GRD_x(g,k0,l,GRD_ZDIR), & ! [IN]
!MINIAPP                              GRD_s(g,k0,l,I_LAT),    & ! [OUT]
!MINIAPP                              GRD_s(g,k0,l,I_LON)     ) ! [OUT]
!MINIAPP
!MINIAPP       call VECTR_xyz2latlon( GRD_xt(g,k0,l,ADM_TI,GRD_XDIR), & ! [IN]
!MINIAPP                              GRD_xt(g,k0,l,ADM_TI,GRD_YDIR), & ! [IN]
!MINIAPP                              GRD_xt(g,k0,l,ADM_TI,GRD_ZDIR), & ! [IN]
!MINIAPP                              GRD_st(g,k0,l,ADM_TI,I_LAT),    & ! [OUT]
!MINIAPP                              GRD_st(g,k0,l,ADM_TI,I_LON)     ) ! [OUT]
!MINIAPP
!MINIAPP       call VECTR_xyz2latlon( GRD_xt(g,k0,l,ADM_TJ,GRD_XDIR), & ! [IN]
!MINIAPP                              GRD_xt(g,k0,l,ADM_TJ,GRD_YDIR), & ! [IN]
!MINIAPP                              GRD_xt(g,k0,l,ADM_TJ,GRD_ZDIR), & ! [IN]
!MINIAPP                              GRD_st(g,k0,l,ADM_TJ,I_LAT),    & ! [OUT]
!MINIAPP                              GRD_st(g,k0,l,ADM_TJ,I_LON)     ) ! [OUT]
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP       do g = 1, ADM_gall_pl
!MINIAPP          call VECTR_xyz2latlon( GRD_x_pl(g,k0,l,GRD_XDIR), & ! [IN]
!MINIAPP                                 GRD_x_pl(g,k0,l,GRD_YDIR), & ! [IN]
!MINIAPP                                 GRD_x_pl(g,k0,l,GRD_ZDIR), & ! [IN]
!MINIAPP                                 GRD_s_pl(g,k0,l,I_LAT),    & ! [OUT]
!MINIAPP                                 GRD_s_pl(g,k0,l,I_LON)     ) ! [OUT]
!MINIAPP
!MINIAPP          call VECTR_xyz2latlon( GRD_xt_pl(g,k0,l,GRD_XDIR), & ! [IN]
!MINIAPP                                 GRD_xt_pl(g,k0,l,GRD_YDIR), & ! [IN]
!MINIAPP                                 GRD_xt_pl(g,k0,l,GRD_ZDIR), & ! [IN]
!MINIAPP                                 GRD_st_pl(g,k0,l,I_LAT),    & ! [OUT]
!MINIAPP                                 GRD_st_pl(g,k0,l,I_LON)     ) ! [OUT]
!MINIAPP       enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_makelatlon

  !-----------------------------------------------------------------------------
  !> set longitude and latitude for IO
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_setlatlon
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl,     &
!MINIAPP       ADM_have_sgp,    &
!MINIAPP       ADM_lall,        &
!MINIAPP       ADM_lall_pl,     &
!MINIAPP       ADM_gall,        &
!MINIAPP       ADM_gall_pl,     &
!MINIAPP       ADM_gall_1d,     &
!MINIAPP       ADM_gmin,        &
!MINIAPP       ADM_gmax,        &
!MINIAPP       ADM_gslf_pl,     &
!MINIAPP       ADM_gmin_pl,     &
!MINIAPP       ADM_gmax_pl,     &
!MINIAPP       ADM_KNONE,       &
!MINIAPP       RGNMNG_l2r,      &
!MINIAPP       RGNMNG_edge_tab, &
!MINIAPP       I_DIR,           &
!MINIAPP       I_SW,            &
!MINIAPP       I_NW,            &
!MINIAPP       I_NE,            &
!MINIAPP       I_SE,            &
!MINIAPP       I_NPL,           &
!MINIAPP       I_SPL
!MINIAPP    use mod_fio, only: &
!MINIAPP       FIO_set_latlon
!MINIAPP    use mod_comm, only: &
!MINIAPP       COMM_data_transfer
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    real(RP) :: latlon_with_around   (ADM_gall   ,7,ADM_lall   ,2) ! lat-lon point, center and neighbors
!MINIAPP    real(RP) :: latlon_with_around_pl(ADM_gall_pl,7,ADM_lall_pl,2) ! lat-lon point, center and neighbors
!MINIAPP
!MINIAPP    integer  :: ij
!MINIAPP    integer  :: ip1j, ijp1, ip1jp1
!MINIAPP    integer  :: im1j, ijm1, im1jm1
!MINIAPP
!MINIAPP    real(RP) :: tmp(7)
!MINIAPP    integer  :: rgnid
!MINIAPP
!MINIAPP    integer  :: i, j, k0, l, d, n, v, vv
!MINIAPP
!MINIAPP    integer  :: suf
!MINIAPP    suf(i,j) = ADM_gall_1d * (j-1) + i
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    do d = I_LAT, I_LON
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       do j = ADM_gmin, ADM_gmax
!MINIAPP       do i = ADM_gmin, ADM_gmax
!MINIAPP          ij     = suf(i  ,j  )
!MINIAPP          ip1j   = suf(i+1,j  )
!MINIAPP          ip1jp1 = suf(i+1,j+1)
!MINIAPP          ijp1   = suf(i  ,j+1)
!MINIAPP          im1j   = suf(i-1,j  )
!MINIAPP          im1jm1 = suf(i-1,j-1)
!MINIAPP          ijm1   = suf(i  ,j-1)
!MINIAPP
!MINIAPP          ! 1 center and 6 neighbors
!MINIAPP          latlon_with_around(ij,1,l,d) = GRD_s(ij    ,k0,l,d)
!MINIAPP          latlon_with_around(ij,2,l,d) = GRD_s(ijm1  ,k0,l,d)
!MINIAPP          latlon_with_around(ij,3,l,d) = GRD_s(ip1j  ,k0,l,d)
!MINIAPP          latlon_with_around(ij,4,l,d) = GRD_s(ip1jp1,k0,l,d)
!MINIAPP          latlon_with_around(ij,5,l,d) = GRD_s(ijp1  ,k0,l,d)
!MINIAPP          latlon_with_around(ij,6,l,d) = GRD_s(im1j  ,k0,l,d)
!MINIAPP          latlon_with_around(ij,7,l,d) = GRD_s(im1jm1,k0,l,d)
!MINIAPP       enddo ! i loop
!MINIAPP       enddo ! j loop
!MINIAPP
!MINIAPP       if ( ADM_have_sgp(l) ) then ! pentagon
!MINIAPP          latlon_with_around(suf(ADM_gmin,ADM_gmin),7,l,d) = latlon_with_around(suf(ADM_gmin,ADM_gmin),2,l,d)
!MINIAPP       endif
!MINIAPP    enddo ! l loop
!MINIAPP    enddo ! d loop
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       n = ADM_gslf_pl
!MINIAPP
!MINIAPP       do d = I_LAT, I_LON
!MINIAPP          latlon_with_around_pl(n,1,I_NPL,d) = GRD_s_pl(n,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,2,I_NPL,d) = GRD_s_pl(5,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,3,I_NPL,d) = GRD_s_pl(4,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,4,I_NPL,d) = GRD_s_pl(3,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,5,I_NPL,d) = GRD_s_pl(2,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,6,I_NPL,d) = GRD_s_pl(6,k0,I_NPL,d)
!MINIAPP          latlon_with_around_pl(n,7,I_NPL,d) = GRD_s_pl(6,k0,I_NPL,d)
!MINIAPP
!MINIAPP          latlon_with_around_pl(n,1,I_SPL,d) = GRD_s_pl(n,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,2,I_SPL,d) = GRD_s_pl(5,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,3,I_SPL,d) = GRD_s_pl(4,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,4,I_SPL,d) = GRD_s_pl(3,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,5,I_SPL,d) = GRD_s_pl(2,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,6,I_SPL,d) = GRD_s_pl(2,k0,I_SPL,d)
!MINIAPP          latlon_with_around_pl(n,7,I_SPL,d) = GRD_s_pl(6,k0,I_SPL,d)
!MINIAPP       enddo ! d loop
!MINIAPP    endif
!MINIAPP
!MINIAPP    call COMM_data_transfer( latlon_with_around, latlon_with_around_pl )
!MINIAPP
!MINIAPP    do d = I_LAT, I_LON
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       rgnid = RGNMNG_l2r(l)
!MINIAPP
!MINIAPP       ! SW halo (Southern Hemisphere, Edge of diamond)
!MINIAPP       if ( RGNMNG_edge_tab(I_DIR,I_SW,rgnid) == I_SE ) then
!MINIAPP          do j = ADM_gmin-1, ADM_gmin-1
!MINIAPP          do i = ADM_gmin-1, ADM_gmax+1
!MINIAPP             ij = suf(i,j)
!MINIAPP
!MINIAPP             tmp(:) = latlon_with_around(ij,:,l,d)
!MINIAPP
!MINIAPP             do v = 2, 7
!MINIAPP                vv = v + 1
!MINIAPP                if( vv > 7 ) vv = 2
!MINIAPP
!MINIAPP                latlon_with_around(ij,vv,l,d) = tmp(v)
!MINIAPP             enddo
!MINIAPP          enddo ! i loop
!MINIAPP          enddo ! j loop
!MINIAPP       endif
!MINIAPP
!MINIAPP       ! NW halo (Northern Hemisphere, Edge of diamond)
!MINIAPP       if ( RGNMNG_edge_tab(I_DIR,I_NW,rgnid) == I_NE ) then
!MINIAPP          do j = ADM_gmin  , ADM_gmax
!MINIAPP          do i = ADM_gmin-1, ADM_gmin-1
!MINIAPP             ij = suf(i,j)
!MINIAPP
!MINIAPP             tmp(:) = latlon_with_around(ij,:,l,d)
!MINIAPP
!MINIAPP             do v = 2, 7
!MINIAPP                vv = v - 1
!MINIAPP                if( vv < 2 ) vv = 7
!MINIAPP
!MINIAPP                latlon_with_around(ij,vv,l,d) = tmp(v)
!MINIAPP             enddo
!MINIAPP          enddo ! i loop
!MINIAPP          enddo ! j loop
!MINIAPP       endif
!MINIAPP
!MINIAPP       ! NE halo (Northern Hemisphere, Edge of diamond)
!MINIAPP       if ( RGNMNG_edge_tab(I_DIR,I_NE,rgnid) == I_NW ) then
!MINIAPP          do j = ADM_gmax+1, ADM_gmax+1
!MINIAPP          do i = ADM_gmin-1, ADM_gmax+1
!MINIAPP             ij = suf(i,j)
!MINIAPP
!MINIAPP             tmp(:) = latlon_with_around(ij,:,l,d)
!MINIAPP
!MINIAPP             do v = 2, 7
!MINIAPP                vv = v + 1
!MINIAPP                if( vv > 7 ) vv = 2
!MINIAPP
!MINIAPP                latlon_with_around(ij,vv,l,d) = tmp(v)
!MINIAPP             enddo
!MINIAPP          enddo ! i loop
!MINIAPP          enddo ! j loop
!MINIAPP       endif
!MINIAPP
!MINIAPP       ! SE halo (Southern Hemisphere, Edge of diamond)
!MINIAPP       if ( RGNMNG_edge_tab(I_DIR,I_SE,rgnid) == I_SW ) then
!MINIAPP          do j = ADM_gmin  , ADM_gmax+1
!MINIAPP          do i = ADM_gmax+1, ADM_gmax+1
!MINIAPP             ij = suf(i,j)
!MINIAPP
!MINIAPP             tmp(:) = latlon_with_around(ij,:,l,d)
!MINIAPP
!MINIAPP             do v = 2, 7
!MINIAPP                vv = v - 1
!MINIAPP                if( vv < 2 ) vv = 7
!MINIAPP
!MINIAPP                latlon_with_around(ij,vv,l,d) = tmp(v)
!MINIAPP             enddo
!MINIAPP          enddo ! i loop
!MINIAPP          enddo ! j loop
!MINIAPP       endif
!MINIAPP    enddo
!MINIAPP    enddo
!MINIAPP
!MINIAPP    call FIO_set_latlon( latlon_with_around(:,:,:,I_LAT), & ![IN]
!MINIAPP                         latlon_with_around(:,:,:,I_LON)  ) ![IN]
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_setlatlon

  !-----------------------------------------------------------------------------
  !> calculate location of the mid-point of cell arc
  !--- MINIAPP: horizontal-grid / topography / lat-lon routines removed: they need mod_comm, mod_fio and mod_vector, and the column-local MP chain never reads GRD_x/GRD_s/GRD_zs. ---
!MINIAPP  subroutine GRD_makearc
!MINIAPP    use mod_adm, only: &
!MINIAPP       ADM_have_pl, &
!MINIAPP       ADM_gall_1d, &
!MINIAPP       ADM_gmin,    &
!MINIAPP       ADM_gmax,    &
!MINIAPP       ADM_gslf_pl, &
!MINIAPP       ADM_gmin_pl, &
!MINIAPP       ADM_gmax_pl
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    integer  :: ij
!MINIAPP    integer  :: im1j, ijm1
!MINIAPP
!MINIAPP    integer  :: nstart,nend
!MINIAPP    integer  :: n, l, v, k0
!MINIAPP
!MINIAPP    integer  :: suf, i, j
!MINIAPP    suf(i,j) = ADM_gall_1d * (j-1) + i
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    k0 = ADM_KNONE
!MINIAPP
!MINIAPP    do l = 1, ADM_lall
!MINIAPP       nstart = suf(ADM_gmin-1,ADM_gmin  )
!MINIAPP       nend   = suf(ADM_gmax  ,ADM_gmax  )
!MINIAPP
!MINIAPP       do n = nstart, nend
!MINIAPP          ij     = n
!MINIAPP          ijm1   = n     - ADM_gall_1d
!MINIAPP
!MINIAPP          GRD_xr(n,k0,l,ADM_AI ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_XDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_XDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AI ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_YDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_YDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AI ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ijm1,k0,l,ADM_TJ,GRD_ZDIR) + GRD_xt(ij,k0,l,ADM_TI,GRD_ZDIR) )
!MINIAPP       enddo
!MINIAPP
!MINIAPP       nstart = suf(ADM_gmin-1,ADM_gmin-1)
!MINIAPP       nend   = suf(ADM_gmax  ,ADM_gmax  )
!MINIAPP
!MINIAPP       do n = nstart, nend
!MINIAPP          ij     = n
!MINIAPP
!MINIAPP          GRD_xr(n,k0,l,ADM_AIJ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_XDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_XDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AIJ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_YDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_YDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AIJ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TI,GRD_ZDIR) + GRD_xt(ij,k0,l,ADM_TJ,GRD_ZDIR) )
!MINIAPP       enddo
!MINIAPP
!MINIAPP       nstart = suf(ADM_gmin  ,ADM_gmin-1)
!MINIAPP       nend   = suf(ADM_gmax  ,ADM_gmax  )
!MINIAPP
!MINIAPP       do n = nstart, nend
!MINIAPP          ij     = n
!MINIAPP          im1j   = n - 1
!MINIAPP
!MINIAPP          GRD_xr(n,k0,l,ADM_AJ ,GRD_XDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_XDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_XDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AJ ,GRD_YDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_YDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_YDIR) )
!MINIAPP          GRD_xr(n,k0,l,ADM_AJ ,GRD_ZDIR) = 0.5_RP * ( GRD_xt(ij,k0,l,ADM_TJ,GRD_ZDIR) + GRD_xt(im1j,k0,l,ADM_TI,GRD_ZDIR) )
!MINIAPP       enddo
!MINIAPP    enddo
!MINIAPP
!MINIAPP    if ( ADM_have_pl ) then
!MINIAPP       n = ADM_gslf_pl
!MINIAPP
!MINIAPP       do l = 1, ADM_lall_pl
!MINIAPP          do v = ADM_gmin_pl, ADM_gmax_pl
!MINIAPP             ij   = v
!MINIAPP             ijm1 = v - 1
!MINIAPP             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl
!MINIAPP
!MINIAPP             GRD_xr_pl(v,k0,l,GRD_XDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_XDIR)+GRD_xt_pl(ij,k0,l,GRD_XDIR))
!MINIAPP             GRD_xr_pl(v,k0,l,GRD_YDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_YDIR)+GRD_xt_pl(ij,k0,l,GRD_YDIR))
!MINIAPP             GRD_xr_pl(v,k0,l,GRD_ZDIR) = 0.5_RP * (GRD_xt_pl(ijm1,k0,l,GRD_ZDIR)+GRD_xt_pl(ij,k0,l,GRD_ZDIR))
!MINIAPP          enddo
!MINIAPP       enddo
!MINIAPP    endif
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine GRD_makearc

end module mod_grd
