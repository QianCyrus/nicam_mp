!-------------------------------------------------------------------------------
!> Module ADMINISTRATION
!!
!! @par Description
!!         Management of process and region
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_adm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
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
  public :: ADM_openparam
  public :: ADM_proc_init
  public :: ADM_proc_stop
  public :: ADM_proc_finish

  public :: ADM_setup
  public :: ADM_mk_suffix

  public :: ADM_MPIbarrier
  public :: ADM_MPItime
  public :: ADM_MPItimestat

  public :: ADM_make_idstr
  public :: ADM_set_prefix_extension

  public :: ADM_proc_submodel
  public :: ADM_setup_single

  public :: RGNMNG_generate_ico
  public :: RGNMNG_vertex_walkaround

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !====== Basic definition & information ======

  integer,  public, parameter :: I_l   = 1   !< local region
  integer,  public, parameter :: I_prc = 2   !< process

  integer,  public, parameter :: I_RGNID = 1 !< region id
  integer,  public, parameter :: I_DIR   = 2 !< direction

  ! Identifiers of directions of region edges
  integer,  public, parameter :: I_SW = 1    !< south west
  integer,  public, parameter :: I_NW = 2    !< north west
  integer,  public, parameter :: I_NE = 3    !< north east
  integer,  public, parameter :: I_SE = 4    !< south east

  ! Identifiers of directions of region vertices
  integer,  public, parameter :: I_W = 1     !< west
  integer,  public, parameter :: I_N = 2     !< north
  integer,  public, parameter :: I_E = 3     !< east
  integer,  public, parameter :: I_S = 4     !< south

  ! Identifier of poles (north pole or south pole)
  integer,  public, parameter :: I_NPL = 1   !< north pole
  integer,  public, parameter :: I_SPL = 2   !< south pole

  ! Identifier of triangle element (i-axis-side or j-axis side)
  integer,  public, parameter :: ADM_TI = 1
  integer,  public, parameter :: ADM_TJ = 2

  ! Identifier of arc element (i-axis-side, ij-axis side, or j-axis side)
  integer,  public, parameter :: ADM_AI  = 1
  integer,  public, parameter :: ADM_AIJ = 2
  integer,  public, parameter :: ADM_AJ  = 3

  ! Identifier of 1 variable
  integer,  public, parameter :: ADM_KNONE = 1

  !#############################################################################
  ! Global Parameters
  !#############################################################################
  integer,  public            :: GLOBAL_glevel                ! grid     division level
  integer,  public            :: GLOBAL_rlevel                ! region   division level
  integer,  public            :: GLOBAL_elevel                ! ensemble division level
  integer,  public            :: GLOBAL_COMM_WORLD            ! communication world per member
  integer,  public            :: GLOBAL_prc_all               ! number of MPI process
  integer,  public            :: GLOBAL_prc_me                ! my process ID
  integer,  public, parameter :: GLOBAL_prc_master  = 1       ! master process ID

  logical,  public            :: GLOBAL_ensemble_on = .false. ! ensemble run is on?

  character(len=H_LONG), public              :: GLOBAL_prefix_dir    = '' !< prefix,    directory name
  character(len=H_LONG), public              :: GLOBAL_extension_ens = '' !< extension, ensemble name
  character(len=H_LONG), public, allocatable :: GLOBAL_extension_rgn(:)   !< extension, region name

#ifdef _FIXEDINDEX_
  include 'inc_index.h'
#else
  !#############################################################################
  ! Basic Index Parameters
  !#############################################################################

  ! main parameter
  integer,  public            :: ADM_glevel           ! grid   division level
  integer,  public            :: ADM_rlevel           ! region division level
  integer,  public            :: ADM_vlayer           ! number of vertical layer
  integer,  public            :: ADM_DMD              ! number of diamond
  integer,  public            :: ADM_prc_all          ! number of MPI process

  ! region
  integer,  public            :: ADM_rgn_nmax         ! number of regular region
  integer,  public            :: ADM_lall             ! number of regular region per process
  integer,  public, parameter :: ADM_rgn_nmax_pl =  2 ! number of pole    region
  integer,  public, parameter :: ADM_lall_pl     =  2 ! number of pole    region per process

  ! horizontal grid
  integer,  public            :: ADM_gall             ! number of horizontal grid per regular region
  integer,  public            :: ADM_gall_in          ! number of horizontal grid (inner part)
  integer,  public            :: ADM_gall_1d          ! number of horizontal grid (1D)
  integer,  public            :: ADM_gmin             ! start index of 1D horizontal grid
  integer,  public            :: ADM_gmax             ! end   index of 1D horizontal grid

  integer,  public            :: ADM_vlink       = -1 ! maximum number of vertex linkage, ICO:5, PSP:6, LCP, MLCP:k
  integer,  public            :: ADM_gall_pl          ! number of horizontal grid for pole region
  integer,  public, parameter :: ADM_gslf_pl     =  1 ! index for pole point
  integer,  public, parameter :: ADM_gmin_pl     =  2 ! start index of grid around the pole point
  integer,  public            :: ADM_gmax_pl          ! end   index of grid around the pole point

  ! vertical grid
  integer,  public            :: ADM_kall             ! number of vertical grid
  integer,  public            :: ADM_kmin             ! start index of vertical grid
  integer,  public            :: ADM_kmax             ! end   index of vertical grid

  ! List vectors
  integer,  public            :: ADM_IooJoo_nmax
  integer,  public            :: ADM_IooJmo_nmax
  integer,  public            :: ADM_IooJop_nmax
  integer,  public            :: ADM_IooJmp_nmax
  integer,  public            :: ADM_ImoJoo_nmax
  integer,  public            :: ADM_ImoJmo_nmax
  integer,  public            :: ADM_ImoJop_nmax
  integer,  public            :: ADM_ImoJmp_nmax
  integer,  public            :: ADM_IopJoo_nmax
  integer,  public            :: ADM_IopJmo_nmax
  integer,  public            :: ADM_IopJop_nmax
  integer,  public            :: ADM_IopJmp_nmax
  integer,  public            :: ADM_ImpJoo_nmax
  integer,  public            :: ADM_ImpJmo_nmax
  integer,  public            :: ADM_ImpJop_nmax
  integer,  public            :: ADM_ImpJmp_nmax
#endif

  !
  !====== Information for processes ======
  !
  integer,  public            :: ADM_COMM_WORLD          ! local communicator (each member)
  integer,  public            :: ADM_prc_me              ! my process ID
  logical,  public            :: ADM_IsMaster  = .false. ! master process in local communicator?
  logical,  public            :: ADM_mpi_alive = .false. ! MPI is alive?

  integer,  public, parameter :: ADM_prc_master = 1      ! master process ID
  integer,  public            :: ADM_prc_pl              ! process ID which manages the pole regions
  logical,  public            :: ADM_have_pl             ! this ID manages pole region?

  !
  !====== Information for devices ======
  !
  integer,  public            :: ADM_device_all          ! number of accelerator device (local node)
  integer,  public            :: ADM_device_me           ! my device ID (local node)

  !
  !====== Information for processes-region relationship ======
  !
  integer,  public, parameter   :: RGNMNG_llim = 2560          !< maximum number of region per process.

  integer,  public, allocatable :: RGNMNG_edge_tab   (:,:,:)   !< region link information (for 4 edges)

  integer,  public, allocatable :: RGNMNG_vert_num   (:,:)     !< number of region around the vertex (4 vertexes)
  integer,  public, allocatable :: RGNMNG_vert_tab   (:,:,:,:) !< region link information (for 4 vertexes)
  integer,  public, allocatable :: RGNMNG_vert_tab_pl(:,:,:)   !< region link information (for 4 vertexes)

  integer,  public, allocatable :: RGNMNG_lnum(:)              !< number of region in each process
  integer,  public, allocatable :: RGNMNG_lp2r(:,:)            !< l,prc       => rgn
  integer,  public, allocatable :: RGNMNG_r2lp(:,:)            !< rgn         => l,prc
  integer,  public, allocatable :: RGNMNG_l2r (:)              !< l,prc_me    => rgn

  integer,  public              :: RGNMNG_r2p_pl(I_NPL:I_SPL)  !< process ID which have the pole regions
  integer,  public              :: RGNMNG_rgn4pl(I_NPL:I_SPL)  !< region, having pole data in the halo

  !
  !====== Information for regions ======
  !
  integer,  public              :: ADM_l_me        ! Present Local region number
  logical,  public, allocatable :: ADM_have_sgp(:) ! region have singlar point?

  !
  !====== Information for grids ======
  !
  integer,  public, parameter   :: ADM_nxyz = 3 ! dimension of the spacial vector

  ! Identifiers of grid points around a grid point
  integer,  public, parameter   :: ADM_GIJ_nmax = 7
  integer,  public, parameter   :: ADM_GIoJo    = 1
  integer,  public, parameter   :: ADM_GIpJo    = 2
  integer,  public, parameter   :: ADM_GIpJp    = 3
  integer,  public, parameter   :: ADM_GIoJp    = 4
  integer,  public, parameter   :: ADM_GImJo    = 5
  integer,  public, parameter   :: ADM_GImJm    = 6
  integer,  public, parameter   :: ADM_GIoJm    = 7

  integer,  public, allocatable :: ADM_IooJoo(:,:)
  integer,  public, allocatable :: ADM_IooJmo(:,:)
  integer,  public, allocatable :: ADM_IooJop(:,:)
  integer,  public, allocatable :: ADM_IooJmp(:,:)
  integer,  public, allocatable :: ADM_ImoJoo(:,:)
  integer,  public, allocatable :: ADM_ImoJmo(:,:)
  integer,  public, allocatable :: ADM_ImoJop(:,:)
  integer,  public, allocatable :: ADM_ImoJmp(:,:)
  integer,  public, allocatable :: ADM_IopJoo(:,:)
  integer,  public, allocatable :: ADM_IopJmo(:,:)
  integer,  public, allocatable :: ADM_IopJop(:,:)
  integer,  public, allocatable :: ADM_IopJmp(:,:)
  integer,  public, allocatable :: ADM_ImpJoo(:,:)
  integer,  public, allocatable :: ADM_ImpJmo(:,:)
  integer,  public, allocatable :: ADM_ImpJop(:,:)
  integer,  public, allocatable :: ADM_ImpJmp(:,:)

  character(len=H_SHORT), public :: ADM_HGRID_SYSTEM = 'ICO' ! [XTMS] Horizontal Grid type
                                                     ! 'ICO'      icosahedral
                                                     ! 'ICO-XTMS' icosahedral but XTMS is used in oprt
                                                     ! 'LCP'      Lambert-cornial (including PSP)
                                                     ! 'MLCP'     Mercator+Lambert-cornial
                                                     ! 'MLCP-OLD' OLD vergion (only for s=1)

  integer,                public :: ADM_XTMS_MLCP_S  = 1 ! [XTMS] Number of segment for MLCP

  logical,                public :: ADM_flg_scm      = .false.
  integer,                public :: ADM_idx_aid4stat = 0

  !
  !====== Information for ensemble ======
  !
  integer,  public              :: ENS_COMM_WORLD     ! communication world for ensemble
  integer,  public              :: ENS_prc_all        ! number of MPI process in the ensemble
  integer,  public              :: ENS_prc_me         ! my process ID in the ensemble
  integer,  public, parameter   :: ENS_prc_master = 1 ! master process ID in the emsemble
  logical,  public              :: ENS_have_pl        ! this ID manages pole grid?

  integer,  public              :: ENS_elevel         ! ensemble division level
  integer,  public              :: ENS_div_nmax       ! number of divided part
  integer,  public              :: ENS_gall           ! number of horizontal grid
  integer,  public, parameter   :: ENS_gall_pl = 2    ! number of horizontal grid, pole point
  integer,  public              :: ENS_eall           ! number of ensemble
  integer,  public              :: ENS_eall_total     ! number of ensemble (total)

  !
  !====== Information for processes-division relationship ======
  !
  integer,  public              :: ENS_ediv           ! number of divided part per process
  integer,  public              :: ENS_gdiv           ! number of horizontal grid per divided part
  integer,  public              :: ENS_gdiv_1d        ! number of horizontal grid(1D) per divided part

  integer,  public, allocatable :: ENS_d2e(:)         ! division ID -> ensemble member ID
  integer,  public, allocatable :: ENS_d2l(:)         ! division ID -> local l
  integer,  public, allocatable :: ENS_imin(:)        ! start i-index of 1D horizontal grid for devide
  integer,  public, allocatable :: ENS_imax(:)        ! end   i-index of 1D horizontal grid for divide
  integer,  public, allocatable :: ENS_jmin(:)        ! start j-index of 1D horizontal grid for devide
  integer,  public, allocatable :: ENS_jmax(:)        ! end   j-index of 1D horizontal grid for divide

  !%%%%% Obsolete variables for backward compatibility %%%%%
  integer,  public, parameter   :: ADM_NSYS     = H_SHORT
  integer,  public, parameter   :: ADM_MAXFNAME = H_FILE
  integer,  public              :: ADM_LOG_FID
  integer,  public              :: ADM_CTL_FID
  integer,  public, allocatable :: ADM_prc_tab(:,:) ! => RGNMNG_lp2r
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: ADM_member_setup
  private :: ADM_ensemble_setup

  private :: RGNMNG_setup
  private :: RGNMNG_input
  private :: RGNMNG_output
  private :: output_info

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_FILE), private            :: RGNMNG_in_fname    = ''  !< input  file name for region management file
  character(len=H_FILE), private            :: RGNMNG_out_fname   = ''  !< output file name for region management file

  character(len=2),      private, parameter :: RGNMNG_edgename(4) = (/'SW','NW','NE','SE'/) !< name table
  character(len=2),      private, parameter :: RGNMNG_vertname(4) = (/'W ','N ','E ','S '/) !< name table

  character(len=H_LONG), private            :: prefix_dir    = ''
  character(len=H_LONG), private            :: extension_ens = 'mem'

  logical, private :: debug = .false. !< debug option

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Open config file
  subroutine ADM_openparam( &
       param_fname )
    implicit none

    character(len=*), intent(inout) :: param_fname ! namelist file name

    character(len=H_FILE) :: fname !< filename

    integer  :: ierr
    !---------------------------------------------------------------------------

    ! try to get config file from command argument
    if ( COMMAND_ARGUMENT_COUNT() == 1 ) then
       call get_command_argument(1,fname)
       param_fname = trim(fname)
    endif

    ! open config file
    IO_FID_CONF = IO_get_available_fid()
    open( unit   = IO_FID_CONF,       &
          file   = trim(param_fname), &
          form   = 'formatted',       &
          status = 'old',             &
          iostat = ierr               )

    if ( ierr /= 0 ) then
       write(*,*) 'xxx Cannot open parameter control file!'
       write(*,*) 'xxx filename:', trim(param_fname)
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    return
  end subroutine ADM_openparam

  !-----------------------------------------------------------------------------
  !> Start MPI (global scope)
  subroutine ADM_proc_init
    !--- MINIAPP: mod_coupler (Jcup) is not part of the microphysics chain and
    !--- drags in the whole Coupler tree; the CPL_jcup_active branch is removed
    !--- and only the plain-MPI path is kept. ---
!MINIAPP    use mod_coupler, only: &
!MINIAPP       CPL_jcup_active,    &
!MINIAPP       CPL_mpi_init,       &
!MINIAPP       CPL_prc_all,        &
!MINIAPP       CPL_my_rank,        &
!MINIAPP       CPL_my_comm
    implicit none

    integer  :: prc_all, my_rank, my_comm
    integer  :: ierr
    !---------------------------------------------------------------------------

!MINIAPP    if ( CPL_jcup_active ) then
!MINIAPP       call CPL_mpi_init
!MINIAPP       my_comm = CPL_my_comm
!MINIAPP       prc_all = CPL_prc_all
!MINIAPP       my_rank = CPL_my_rank
!MINIAPP    else
       call MPI_Init(ierr)
       my_comm = MPI_COMM_WORLD
       call MPI_Comm_size(MPI_COMM_WORLD,prc_all,ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
!MINIAPP   endif

    GLOBAL_COMM_WORLD = my_comm
    GLOBAL_prc_all    = prc_all
    GLOBAL_prc_me     = my_rank + 1

    ADM_mpi_alive = .true.

#ifdef _OPENACC
    ADM_device_all = acc_get_num_devices(acc_device_nvidia)
    ADM_device_me  = mod(my_rank,ADM_device_all)
    call acc_set_device_num(ADM_device_me,acc_device_nvidia)
#endif

    return
  end subroutine ADM_proc_init

  !-----------------------------------------------------------------------------
  !> get global MPI parameters (act as subcomponent)
  subroutine ADM_proc_submodel( &
       my_comm, &
       prc_all, &
       my_rank  )
    implicit none

    integer, intent(in) :: my_comm
    integer, intent(in) :: prc_all
    integer, intent(in) :: my_rank
    !---------------------------------------------------------------------------

    GLOBAL_COMM_WORLD = my_comm
    GLOBAL_prc_all    = prc_all
    GLOBAL_prc_me     = my_rank + 1

    ADM_mpi_alive     = .true.

    return
  end subroutine ADM_proc_submodel

  !-----------------------------------------------------------------------------
  !> Abort MPI process
  subroutine ADM_proc_stop
    implicit none

    integer  :: ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '+++ Abort MPI, rank(GLOBAL,ADM,ENS)=', GLOBAL_prc_me-1, ADM_prc_me-1, ENS_prc_me-1
    write(*,*)          '+++ Abort MPI, rank(GLOBAL,ADM,ENS)=', GLOBAL_prc_me-1, ADM_prc_me-1, ENS_prc_me-1

    call IO_close_allfile

    ! Abort MPI
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)

    stop
  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------------
  !> Finish MPI process
  subroutine ADM_proc_finish
    !--- MINIAPP: mod_coupler (Jcup) removed, see ADM_proc_init. ---
!MINIAPP    use mod_coupler, only: &
!MINIAPP       CPL_jcup_active,  &
!MINIAPP       CPL_mpi_finalize
    implicit none

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( ADM_mpi_alive ) then

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++ finalize MPI'
       call MPI_Barrier(ADM_COMM_WORLD,ierr)

!MINIAPP       if ( CPL_jcup_active ) then
!MINIAPP          call CPL_mpi_finalize
!MINIAPP       else
          call MPI_Finalize(ierr)
!MINIAPP       endif

       write(IO_FID_LOG,*) '*** MPI is peacefully finalized'
    else
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '+++ stop serial process.'
    endif

    call IO_close_allfile

    return
  end subroutine ADM_proc_finish

  !-----------------------------------------------------------------------------
  !--- MINIAPP: ------------------------------------------------------------
  !--- The original ADM_setup builds the whole process/region management from
  !--- &GLOBALPARAM + &IOPARAM + &ADMPARAM and the region-management file
  !--- (rl<RL>-prc<NP>.info), because NICAM needs the region topology for the
  !--- halo exchange.  The microphysics chain is column-local: it never calls
  !--- COMM, never touches a neighbour cell and never uses the _pl poles, so
  !--- the topology is dead weight and the region count of the mini-app is a
  !--- property of the DUMP (10 for 'full', 1 for 'region0'), not of the mesh.
  !--- The body below is therefore a MINIAPP init that opens the log file and
  !--- derives the very same index set from &ADMPARAM (glevel/rlevel/vlayer)
  !--- plus a mini-app-only 'lall' knob.  The index arithmetic is copied
  !--- verbatim from ADM_member_setup (kept below, now unreferenced).
  !--- The original body is preserved, commented out with the !MINIAPP tag,
  !--- immediately after this routine.
  !-------------------------------------------------------------------------
  subroutine ADM_setup( &
       param_fname, &
       msg_base     )
    implicit none

    character(len=*), intent(in)           :: param_fname ! namelist file name
    character(len=*), intent(in), optional :: msg_base    ! [add] H.Yashiro 20110701

    integer               :: glevel      = -1
    integer               :: rlevel      = -1
    integer               :: vlayer      =  1
    integer               :: lall        = -1 !--- MINIAPP: number of column blocks in the dump
    character(len=H_FILE) :: rgnmngfname = ''

    namelist / ADMPARAM / &
       glevel,           &
       rlevel,           &
       vlayer,           &
       lall,             &
       rgnmngfname,      &
       ADM_HGRID_SYSTEM, &
       ADM_vlink,        &
       ADM_XTMS_MLCP_S,  &
       debug

    character(len=H_LONG)  :: log_basename          = 'msg'
    character(len=H_SHORT) :: namelist_out_type     = 'LOG'
    character(len=H_FILE)  :: namelist_out_filename = 'namelist.cnf'

    namelist / IOPARAM / &
       log_basename,         &
       namelist_out_type,    &
       namelist_out_filename

    character(len=H_FILE) :: LOG_fname
    character(len=H_FILE) :: NML_fname
    character(len=H_LONG) :: msg

    integer  :: nmax, dmd
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( IO_FID_CONF <= 0 ) then ! config file is not opened yet
       IO_FID_CONF = IO_get_available_fid()
       open( unit   = IO_FID_CONF,       &
             file   = trim(param_fname), &
             form   = 'formatted',       &
             status = 'old',             &
             iostat = ierr               )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx Cannot open parameter control file!'
          write(*,*) 'xxx filename:', trim(param_fname)
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif
    endif

    !--- MINIAPP: no ensemble, no member split: one communicator, one world.
    GLOBAL_ensemble_on   = .false.
    GLOBAL_glevel        = -1
    GLOBAL_rlevel        = -1
    GLOBAL_elevel        = 0
    ADM_COMM_WORLD       = GLOBAL_COMM_WORLD
    ADM_prc_all          = GLOBAL_prc_all
    ADM_prc_me           = GLOBAL_prc_me
    ENS_COMM_WORLD       = -1
    ENS_prc_all          = 0
    ENS_prc_me           = -1
    ENS_eall_total       = 0
    GLOBAL_prefix_dir    = ''
    GLOBAL_extension_ens = ''

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=IOPARAM,iostat=ierr)
    if ( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist IOPARAM. STOP.'
       call ADM_proc_stop
    endif
    IO_LOG_BASENAME = trim(log_basename)

    msg = trim(IO_LOG_BASENAME)
    if ( present(msg_base) ) msg = trim(msg_base)

    ! open log file
    IO_FID_LOG = IO_get_available_fid()
    write(LOG_fname,'(A,A,A,I6.6)') trim(GLOBAL_prefix_dir), trim(msg), '.pe', GLOBAL_prc_me-1
    open( unit   = IO_FID_LOG,      &
          file   = trim(LOG_fname), &
          form   = 'formatted',     &
          status = 'replace',       &
          iostat = ierr             )

    write(IO_FID_LOG,*) '############################################################'
    write(IO_FID_LOG,*) '#                                                          #'
    write(IO_FID_LOG,*) '#   NICAM : NSW6 cloud-microphysics MINI-APP               #'
    write(IO_FID_LOG,*) '#                                                          #'
    write(IO_FID_LOG,*) '############################################################'
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)')        '====== Process management info. (global) ======'
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Total number of process            : ', GLOBAL_prc_all
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Communicator ID                    : ', GLOBAL_COMM_WORLD
    write(IO_FID_LOG,'(1x,A,I10)')    '--- My Process rank (+1)               : ', GLOBAL_prc_me
#ifdef _OPENACC
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Number of device per node          : ', ADM_device_all
    write(IO_FID_LOG,'(1x,A,I10)')    '--- My device rank                     : ', ADM_device_me
#endif
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)')        '====== Basic I/O info. ======'
    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Name of log    file                : ', trim(LOG_fname)
    write(IO_FID_LOG,'(1x,A,I10)')    '--- FID  of log    file                : ', IO_FID_LOG
    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Name of config file                : ', trim(param_fname)
    write(IO_FID_LOG,'(1x,A,I10)')    '--- FID  of config file                : ', IO_FID_CONF

    if    ( namelist_out_type == 'NONE' ) then
       IO_FID_NML = -1
       write(IO_FID_LOG,'(1x,A,A)')   '--- Name of namelist file              : surpressed'
    elseif( namelist_out_type == 'LOG'  ) then
       IO_FID_NML = IO_FID_LOG
       write(IO_FID_LOG,'(1x,A,A)'  ) '--- Name of namelist file              : merged to logfile'
    elseif( namelist_out_type == 'FILE' ) then
       if ( ADM_prc_me == ADM_prc_master ) then
          IO_FID_NML = IO_get_available_fid()
          NML_fname = trim(GLOBAL_prefix_dir)//trim(namelist_out_filename)
          open( unit   = IO_FID_NML,      &
                file   = trim(NML_fname), &
                form   = 'formatted',     &
                status = 'replace',       &
                iostat = ierr             )
       else
          IO_FID_NML = -1
       endif
       write(IO_FID_LOG,'(1x,A,A)'  ) '--- Name of namelist file              : ', trim(namelist_out_filename)
       write(IO_FID_LOG,'(1x,A,I10)') '--- FID  of namelist file              : ', IO_FID_NML
    else
       write(IO_FID_LOG,*) 'xxx Not appropriate names is set to namelist_out_type = ', trim(namelist_out_type)
       write(*         ,*) 'xxx Not appropriate names is set to namelist_out_type = ', trim(namelist_out_type)
       call ADM_proc_stop
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)')        '====== Precision info. ======'
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Single precision kind ~= byte      : ', SP
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Double precision kind ~= byte      : ', DP
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (global default)    : ', RP
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (physics)           : ', RP_PHY
    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (-microphysics)     : ', RP_MP

    !--- read parameters (mirrors ADM_member_setup)
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[adm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*         ,*) 'xxx Not found namelist ADMPARAM! STOP.'
       write(IO_FID_LOG,*) 'xxx Not found namelist ADMPARAM! STOP.'
       call ADM_proc_stop
    elseif ( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=ADMPARAM)

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*         ,*) 'xxx glevel is not appropriate :', glevel
       write(IO_FID_LOG,*) 'xxx glevel is not appropriate :', glevel
       call ADM_proc_stop
    endif
    if ( rlevel < 0 ) then
       write(*         ,*) 'xxx rlevel is not appropriate :', rlevel
       write(IO_FID_LOG,*) 'xxx rlevel is not appropriate :', rlevel
       call ADM_proc_stop
    endif

    !--- MINIAPP: only the icosahedral index set is reachable from the MP chain.
    ADM_HGRID_SYSTEM = 'ICO'
    ADM_vlink   = 5
    dmd         = 10
    ADM_prc_pl  = 1
    ADM_gall_pl = ADM_vlink + 1
    ADM_gmax_pl = ADM_vlink + 1

    ADM_glevel   = glevel
    ADM_rlevel   = rlevel
    ADM_vlayer   = vlayer
    ADM_DMD      = dmd

    ADM_rgn_nmax = 2**ADM_rlevel * 2**ADM_rlevel * ADM_DMD
    !--- MINIAPP: lall comes from the dump, not from rgn_nmax/prc_all.
    if ( lall < 1 ) lall = ADM_rgn_nmax / ADM_prc_all
    ADM_lall     = lall

    nmax         = 2**(ADM_glevel-ADM_rlevel)
    ADM_gall_1d  = 1 + nmax + 1
    ADM_gmin     = 1 + 1
    ADM_gmax     = 1 + nmax

    ADM_gall     = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in  = (   nmax+1 ) * (   nmax+1 )

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif

    if ( ADM_prc_me == ADM_prc_master ) then
       ADM_IsMaster = .true.
    else
       ADM_IsMaster = .false.
    endif

    !--- MINIAPP: no region management, no pole region, no singular point.
    ADM_have_pl   = .false.
    ADM_l_me      = 0
    RGNMNG_r2p_pl(:) = -1
    RGNMNG_rgn4pl(:) = -1
    allocate( ADM_have_sgp(ADM_lall) )
    ADM_have_sgp(:) = .false.
    allocate( RGNMNG_lnum(ADM_prc_all)          ); RGNMNG_lnum(:)   = ADM_lall
    allocate( RGNMNG_lp2r(ADM_lall,ADM_prc_all) )
    allocate( RGNMNG_l2r (ADM_lall)             )
    do dmd = 1, ADM_lall
       RGNMNG_lp2r(dmd,:) = dmd
       RGNMNG_l2r (dmd)   = dmd
    enddo

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)')     '====== Grid dimension info. (MINIAPP) ======'
    write(IO_FID_LOG,'(1x,A,I10)') '--- grid   division level (grid)       : ', ADM_glevel
    write(IO_FID_LOG,'(1x,A,I10)') '--- region division level (rgn )       : ', ADM_rlevel
    write(IO_FID_LOG,'(1x,A,I10)') '--- number of vertical layer           : ', ADM_vlayer
    write(IO_FID_LOG,'(1x,A,I10)') '--- number of local region             : ', ADM_lall
    write(IO_FID_LOG,'(1x,A,I10)') '--- horizontal grid (1D)               : ', ADM_gall_1d
    write(IO_FID_LOG,'(1x,A,I10)') '--- horizontal grid (whole)            : ', ADM_gall
    write(IO_FID_LOG,'(1x,A,I10)') '--- horizontal grid (inner, used by MP): ', ADM_gall_in
    write(IO_FID_LOG,'(1x,A,I10)') '--- vertical   grid (whole)            : ', ADM_kall
    write(IO_FID_LOG,'(1x,A,2I10)')'--- vertical   grid (kmin,kmax)        : ', ADM_kmin, ADM_kmax

    !%%%%% Obsolete variables for backward compatibility %%%%%
    ADM_CTL_FID      = IO_FID_CONF
    ADM_LOG_FID      = IO_FID_LOG
    allocate( ADM_prc_tab(ADM_lall,ADM_prc_all) )
    ADM_prc_tab(:,:) = RGNMNG_lp2r(:,:)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    return
  end subroutine ADM_setup

  !--- MINIAPP: original ADM_setup, preserved verbatim but commented out. ----
!MINIAPP  subroutine ADM_setup( &
!MINIAPP       param_fname, &
!MINIAPP       msg_base     )
!MINIAPP    use mod_coupler, only: &
!MINIAPP       CPL_jcup_active,  &
!MINIAPP       CPL_modelname_me, &
!MINIAPP       CPL_exists_coco,  &
!MINIAPP       CPL_exists_io,    &
!MINIAPP       CPL_exists_pair,  &
!MINIAPP       CPL_exists_python
!MINIAPP    implicit none
!MINIAPP
!MINIAPP    character(len=*), intent(in)           :: param_fname ! namelist file name
!MINIAPP    character(len=*), intent(in), optional :: msg_base    ! [add] H.Yashiro 20110701
!MINIAPP
!MINIAPP    integer :: glevel         = -1
!MINIAPP    integer :: rlevel         = -1
!MINIAPP    integer :: elevel         = -1
!MINIAPP    integer :: prc_member     = -1
!MINIAPP    integer :: prc_ensemble   = -1
!MINIAPP    integer :: ensemble_total = -1
!MINIAPP    integer :: nmem_offset    = 0
!MINIAPP
!MINIAPP    namelist / GLOBALPARAM / &
!MINIAPP       glevel,         &
!MINIAPP       rlevel,         &
!MINIAPP       elevel,         &
!MINIAPP       prc_member,     &
!MINIAPP       prc_ensemble,   &
!MINIAPP       prefix_dir,     &
!MINIAPP       extension_ens,  &
!MINIAPP       ensemble_total, &
!MINIAPP       nmem_offset
!MINIAPP
!MINIAPP    character(len=H_LONG)  :: log_basename          = 'msg'
!MINIAPP    character(len=H_SHORT) :: namelist_out_type     = 'LOG'
!MINIAPP    character(len=H_FILE)  :: namelist_out_filename = 'namelist.cnf'
!MINIAPP
!MINIAPP    namelist / IOPARAM / &
!MINIAPP       log_basename,         &
!MINIAPP       namelist_out_type,    &
!MINIAPP       namelist_out_filename
!MINIAPP
!MINIAPP    character(len=H_FILE) :: CTL_fname
!MINIAPP    character(len=H_FILE) :: LOG_fname
!MINIAPP    character(len=H_FILE) :: NML_fname
!MINIAPP    character(len=H_LONG) :: msg ! [add] H.Yashiro 20110701
!MINIAPP
!MINIAPP    integer  :: member_me, ensemble_me
!MINIAPP    integer  :: my_rank, prc_all
!MINIAPP    integer  :: ierr
!MINIAPP    !---------------------------------------------------------------------------
!MINIAPP
!MINIAPP    if ( IO_FID_CONF <= 0 ) then ! config file is not opened yet
!MINIAPP       IO_FID_CONF = IO_get_available_fid()
!MINIAPP       open( unit   = IO_FID_CONF,       &
!MINIAPP             file   = trim(param_fname), &
!MINIAPP             form   = 'formatted',       &
!MINIAPP             status = 'old',             &
!MINIAPP             iostat = ierr               )
!MINIAPP
!MINIAPP       if ( ierr /= 0 ) then
!MINIAPP          write(*,*) 'xxx Cannot open parameter control file!'
!MINIAPP          write(*,*) 'xxx filename:', trim(param_fname)
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP    endif
!MINIAPP
!MINIAPP    rewind(IO_FID_CONF)
!MINIAPP    read(IO_FID_CONF,nml=GLOBALPARAM,iostat=ierr)
!MINIAPP    if ( ierr == 0 ) then
!MINIAPP       GLOBAL_ensemble_on = .true.
!MINIAPP    elseif ( ierr < 0 ) then
!MINIAPP       GLOBAL_ensemble_on = .false.
!MINIAPP    elseif ( ierr > 0 ) then
!MINIAPP       write(*,*) 'xxx Not appropriate names in namelist GLOBALPARAM. STOP.'
!MINIAPP       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP    endif
!MINIAPP
!MINIAPP    GLOBAL_glevel = glevel
!MINIAPP    GLOBAL_rlevel = rlevel
!MINIAPP    GLOBAL_elevel = elevel
!MINIAPP
!MINIAPP    if ( GLOBAL_ensemble_on ) then
!MINIAPP       ! check process number
!MINIAPP       if ( GLOBAL_prc_all /= prc_member *  prc_ensemble ) then
!MINIAPP          write(*,*) 'xxx total process must be prc_member *  prc_ensemble :', &
!MINIAPP                     GLOBAL_prc_all, prc_member, prc_ensemble
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       ! check glevel & rlevel & elevel
!MINIAPP       if ( glevel < 1 ) then
!MINIAPP          write(*,*) 'xxx glevel is not appropriate :', glevel
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( rlevel < 0 ) then
!MINIAPP          write(*,*) 'xxx rlevel is not appropriate :', rlevel
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( elevel < 0 ) then
!MINIAPP          write(*,*) 'xxx elevel is not appropriate :', elevel
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( glevel < rlevel + elevel ) then
!MINIAPP          write(*,*) 'xxx glevel must be rlevel + elevel :', glevel, rlevel, elevel
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       ensemble_me = ( GLOBAL_prc_me - 1 ) / prc_member + 1
!MINIAPP       member_me   = mod(GLOBAL_prc_me,prc_member)
!MINIAPP
!MINIAPP       !--- split to each members
!MINIAPP       call MPI_Comm_split(GLOBAL_COMM_WORLD, ensemble_me, GLOBAL_prc_me, ADM_COMM_WORLD, ierr)
!MINIAPP       call MPI_Comm_size (ADM_COMM_WORLD,prc_all,ierr)
!MINIAPP       call MPI_Comm_rank (ADM_COMM_WORLD,my_rank,ierr)
!MINIAPP       ADM_prc_all = prc_all
!MINIAPP       ADM_prc_me  = my_rank + 1
!MINIAPP
!MINIAPP       !--- split to local ensemble groups
!MINIAPP       call MPI_Comm_split(GLOBAL_COMM_WORLD, prc_ensemble+member_me, GLOBAL_prc_me, ENS_COMM_WORLD, ierr)
!MINIAPP       call MPI_Comm_size (ENS_COMM_WORLD,prc_all,ierr)
!MINIAPP       call MPI_Comm_rank (ENS_COMM_WORLD,my_rank,ierr)
!MINIAPP       ENS_prc_all = prc_all
!MINIAPP       ENS_prc_me  = my_rank + 1
!MINIAPP
!MINIAPP       if ( prefix_dir == '' ) then
!MINIAPP          GLOBAL_prefix_dir = ''
!MINIAPP       else
!MINIAPP          write(GLOBAL_prefix_dir,'(A,I6.6,A)') trim(prefix_dir), ENS_prc_me-1+nmem_offset, '/'
!MINIAPP       endif
!MINIAPP
!MINIAPP       if ( extension_ens == '' ) then
!MINIAPP          GLOBAL_extension_ens = ''
!MINIAPP       else
!MINIAPP          write(GLOBAL_extension_ens,'(A,A,I6.6)') '.', trim(extension_ens), ENS_prc_me-1+nmem_offset
!MINIAPP       endif
!MINIAPP
!MINIAPP       ! close global config file
!MINIAPP       close(IO_FID_CONF)
!MINIAPP       ! open config file for each member
!MINIAPP       CTL_fname = trim(GLOBAL_prefix_dir)//trim(param_fname)//trim(GLOBAL_extension_ens)
!MINIAPP       open( unit   = IO_FID_CONF,     &
!MINIAPP             file   = trim(CTL_fname), &
!MINIAPP             form   = 'formatted',     &
!MINIAPP             status = 'old',           &
!MINIAPP             iostat = ierr             )
!MINIAPP
!MINIAPP       if ( ierr /= 0 ) then
!MINIAPP          write(*,*) 'xxx Cannot open parameter control file!'
!MINIAPP          write(*,*) 'xxx filename:', trim(CTL_fname)
!MINIAPP          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
!MINIAPP       endif
!MINIAPP
!MINIAPP       ENS_eall_total = ensemble_total
!MINIAPP
!MINIAPP    else
!MINIAPP
!MINIAPP       GLOBAL_elevel  = 0
!MINIAPP
!MINIAPP       ADM_COMM_WORLD = GLOBAL_COMM_WORLD
!MINIAPP       ADM_prc_all    = GLOBAL_prc_all
!MINIAPP       ADM_prc_me     = GLOBAL_prc_me
!MINIAPP
!MINIAPP       ENS_COMM_WORLD = -1
!MINIAPP       ENS_prc_all    = 0
!MINIAPP       ENS_prc_me     = -1
!MINIAPP
!MINIAPP       ENS_eall_total = 0
!MINIAPP
!MINIAPP       GLOBAL_prefix_dir    = ''
!MINIAPP       GLOBAL_extension_ens = ''
!MINIAPP
!MINIAPP       ! use current config file
!MINIAPP       CTL_fname = param_fname
!MINIAPP
!MINIAPP    endif
!MINIAPP
!MINIAPP    rewind(IO_FID_CONF)
!MINIAPP    read(IO_FID_CONF,nml=IOPARAM,iostat=ierr)
!MINIAPP    if ( ierr > 0 ) then
!MINIAPP       write(*,*) 'xxx Not appropriate names in namelist IOPARAM. STOP.'
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP    IO_LOG_BASENAME = trim(log_basename)
!MINIAPP
!MINIAPP    msg = trim(IO_LOG_BASENAME)
!MINIAPP    if ( present(msg_base) ) msg = trim(msg_base)
!MINIAPP
!MINIAPP    ! open log file
!MINIAPP    IO_FID_LOG = IO_get_available_fid()
!MINIAPP    write(LOG_fname,'(A,A,A,I6.6)') trim(GLOBAL_prefix_dir), trim(msg), '.pe', GLOBAL_prc_me-1
!MINIAPP    open( unit   = IO_FID_LOG,     &
!MINIAPP          file   = trim(LOG_fname), &
!MINIAPP          form   = 'formatted',     &
!MINIAPP          status = 'replace',       &
!MINIAPP          iostat = ierr             )
!MINIAPP
!MINIAPP    write(IO_FID_LOG,*) '############################################################'
!MINIAPP    write(IO_FID_LOG,*) '#                                                          #'
!MINIAPP    write(IO_FID_LOG,*) '#   NICAM : Nonhydrostatic ICosahedal Atmospheric Model    #'
!MINIAPP    write(IO_FID_LOG,*) '#                                                          #'
!MINIAPP    write(IO_FID_LOG,*) '############################################################'
!MINIAPP#ifdef DEBUG
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,'(1x,A)')        '--- DEBUG mode'
!MINIAPP#endif
!MINIAPP#ifdef QUICKDEBUG
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,'(1x,A)')        '--- QUICKDEBUG mode'
!MINIAPP#endif
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,'(1x,A)')        '====== Process management info. (global) ======'
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Total number of process            : ', GLOBAL_prc_all
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Communicator ID                    : ', GLOBAL_COMM_WORLD
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- My Process rank (+1)               : ', GLOBAL_prc_me
!MINIAPP#ifdef _OPENACC
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Number of device per node          : ', ADM_device_all
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- My device rank                     : ', ADM_device_me
!MINIAPP#endif
!MINIAPP    if ( CPL_jcup_active ) then
!MINIAPP       write(IO_FID_LOG,'(1x,A)')     '--- Use coupling library Jcup?         :        YES'
!MINIAPP       write(IO_FID_LOG,'(1x,A,A10)') '--- My modelname                       : ', trim(CPL_modelname_me)
!MINIAPP       write(IO_FID_LOG,'(1x,A,L10)') '--- COCO       exists?                 : ', CPL_exists_coco
!MINIAPP       write(IO_FID_LOG,'(1x,A,L10)') '--- CoupledIO  exists?                 : ', CPL_exists_io
!MINIAPP       write(IO_FID_LOG,'(1x,A,L10)') '--- NICAM-pair exists?                 : ', CPL_exists_pair
!MINIAPP       write(IO_FID_LOG,'(1x,A,L10)') '--- CoupledPY  exists?                 : ', CPL_exists_python
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,'(1x,A)')     '--- Use coupling library Jcup?         :         NO'
!MINIAPP    endif
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,'(1x,A)')        '====== Basic I/O info. ======'
!MINIAPP    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Prefix (dir)                       : ', trim(GLOBAL_prefix_dir)
!MINIAPP    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Extension for ensemble             : ', trim(GLOBAL_extension_ens)
!MINIAPP    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Name of log    file                : ', trim(LOG_fname)
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- FID  of log    file                : ', IO_FID_LOG
!MINIAPP    write(IO_FID_LOG,'(1x,A,A)'  )    '--- Name of config file                : ', trim(CTL_fname)
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- FID  of config file                : ', IO_FID_CONF
!MINIAPP
!MINIAPP    if    ( namelist_out_type == "NONE" ) then
!MINIAPP       IO_FID_NML = -1
!MINIAPP       write(IO_FID_LOG,'(1x,A,A)')   '--- Name of namelist file              : surpressed'
!MINIAPP    elseif( namelist_out_type == "LOG"  ) then
!MINIAPP       IO_FID_NML = IO_FID_LOG
!MINIAPP       write(IO_FID_LOG,'(1x,A,A)'  ) '--- Name of namelist file              : merged to logfile'
!MINIAPP    elseif( namelist_out_type == "FILE" ) then
!MINIAPP       if ( ADM_prc_me == ADM_prc_master ) then
!MINIAPP          IO_FID_NML = IO_get_available_fid()
!MINIAPP          NML_fname = trim(GLOBAL_prefix_dir)//trim(namelist_out_filename)
!MINIAPP          open( unit   = IO_FID_NML,      &
!MINIAPP                file   = trim(NML_fname), &
!MINIAPP                form   = 'formatted',     &
!MINIAPP                status = 'replace',       &
!MINIAPP                iostat = ierr             )
!MINIAPP       else
!MINIAPP          IO_FID_NML = -1
!MINIAPP       endif
!MINIAPP
!MINIAPP       write(IO_FID_LOG,'(1x,A,A)'  ) '--- Name of namelist file              : ', trim(namelist_out_filename)
!MINIAPP       write(IO_FID_LOG,'(1x,A,I10)') '--- FID  of namelist file              : ', IO_FID_NML
!MINIAPP    else
!MINIAPP       write(IO_FID_LOG,*) 'xxx Not appropriate names is set to namelist_out_type = ', trim(namelist_out_type)
!MINIAPP       write(IO_FID_LOG,*) 'xxx Please set NONE, LOG, or FILE'
!MINIAPP       write(*         ,*) 'xxx Not appropriate names is set to namelist_out_type = ', trim(namelist_out_type)
!MINIAPP       write(*         ,*) 'xxx Please set NONE, LOG, or FILE'
!MINIAPP       call ADM_proc_stop
!MINIAPP    endif
!MINIAPP
!MINIAPP    write(IO_FID_LOG,*)
!MINIAPP    write(IO_FID_LOG,'(1x,A)')        '====== Precision info. ======'
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Single precision kind ~= byte      : ', SP
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Double precision kind ~= byte      : ', DP
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Single precision digit             : ', SP_PREC
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Double precision digit             : ', DP_PREC
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (global default)    : ', RP
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (dynamics)          : ', RP_DYN
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (physics)           : ', RP_PHY
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (-microphysics)     : ', RP_MP
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (-radiation)        : ', RP_RD
!MINIAPP    write(IO_FID_LOG,'(1x,A,I10)')    '--- Precision kind (-turbulence)       : ', RP_TB
!MINIAPP
!MINIAPP    call ADM_member_setup
!MINIAPP
!MINIAPP    call ADM_ensemble_setup
!MINIAPP
!MINIAPP    !%%%%% Obsolete variables for backward compatibility %%%%%
!MINIAPP    ADM_CTL_FID      = IO_FID_CONF
!MINIAPP    ADM_LOG_FID      = IO_FID_LOG
!MINIAPP    allocate( ADM_prc_tab(ADM_lall,ADM_prc_all) )
!MINIAPP    ADM_prc_tab(:,:) = RGNMNG_lp2r(:,:)
!MINIAPP    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!MINIAPP
!MINIAPP    return
!MINIAPP  end subroutine ADM_setup
  !--- MINIAPP: end of original ADM_setup -----------------------------------

  !-----------------------------------------------------------------------------
  subroutine ADM_setup_single( &
       msg_base )
    implicit none

    character(len=*), intent(in) :: msg_base

    character(len=H_FILE) :: log_basename

    namelist / IOPARAM / &
       log_basename

    character(len=H_FILE) :: LOG_fname

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( IO_FID_CONF <= 0 ) then ! config file is not opened yet
       write(*,*) 'xxx ADM_openparam is not called yet! STOP.'
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif

    ADM_COMM_WORLD = GLOBAL_COMM_WORLD
    ADM_prc_all    = GLOBAL_prc_all ! will be replaced after
    ADM_prc_me     = GLOBAL_prc_me  ! will be replaced after

    GLOBAL_elevel        = 0
    ENS_COMM_WORLD       = -1
    ENS_prc_all          = 0
    ENS_prc_me           = -1
    ENS_eall_total       = 0
    GLOBAL_prefix_dir    = ''
    GLOBAL_extension_ens = ''

    log_basename = trim(msg_base)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=IOPARAM,iostat=ierr)
    if ( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist IOPARAM. STOP.'
       call ADM_proc_stop
    endif

    ! open log file
    IO_FID_LOG = IO_get_available_fid()
    write(LOG_fname,'(A,A,I6.6)') trim(log_basename), '.pe', GLOBAL_prc_me-1
    open( unit   = IO_FID_LOG,     &
          file   = trim(LOG_fname), &
          form   = 'formatted',     &
          status = 'replace',       &
          iostat = ierr             )

    write(IO_FID_LOG,*) '############################################################'
    write(IO_FID_LOG,*) '#                                                          #'
    write(IO_FID_LOG,*) '#   Simple MPI execution for NICAM tools                   #'
    write(IO_FID_LOG,*) '#                                                          #'
    write(IO_FID_LOG,*) '############################################################'
#ifdef DEBUG
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)' )    '--- DEBUG mode'
#endif
#ifdef QUICKDEBUG
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)' )    '--- QUICKDEBUG mode'
#endif
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Process management info. (global) ======'
    write(IO_FID_LOG,'(1x,A,I11)') '--- Total number of process            : ', GLOBAL_prc_all
    write(IO_FID_LOG,'(1x,A,I11)') '--- Communicator ID                    : ', GLOBAL_COMM_WORLD
    write(IO_FID_LOG,'(1x,A,I11)') '--- My Process rank (+1)               : ', GLOBAL_prc_me
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Basic I/O info. ======'
    write(IO_FID_LOG,'(1x,A,A)'  ) '--- Name of log    file                : ', trim(LOG_fname)
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- FID  of log    file                : ', IO_FID_LOG
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- FID  of config file                : ', IO_FID_CONF

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Precision info. ======'
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- Single precision kind ~= byte      : ', SP
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- Double precision kind ~= byte      : ', DP
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- Single precision digit             : ', SP_PREC
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- Double precision digit             : ', DP_PREC
    write(IO_FID_LOG,'(1x,A,I7)' ) '--- Precision kind (global default)    : ', RP

    return
  end subroutine ADM_setup_single

  !-----------------------------------------------------------------------------
  subroutine ADM_member_setup
    implicit none

    integer               :: glevel      = -1
    integer               :: rlevel      = -1
    integer               :: vlayer      =  1
    character(len=H_FILE) :: rgnmngfname = ''

    namelist / ADMPARAM / &
       glevel,           &
       rlevel,           &
       vlayer,           &
       rgnmngfname,      &
       ADM_HGRID_SYSTEM, &
       ADM_vlink,        &
       ADM_XTMS_MLCP_S,  &
       debug

    integer  :: nmax, dmd
    integer  :: l, rgnid
    integer  :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[adm]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ADMPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(*         ,*) 'xxx Not found namelist! STOP.'
       write(IO_FID_LOG,*) 'xxx Not found namelist! STOP.'
       call ADM_proc_stop
    elseif ( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist ADMPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=ADMPARAM)

    if ( ADM_HGRID_SYSTEM == 'SCM' ) then
       write(IO_FID_LOG,*) '*** glevel & rlevel are modified for SCM.'
       glevel = 1
       rlevel = 0
    endif

    ! Error if glevel & rlevel are not defined
    if ( glevel < 1 ) then
       write(*         ,*) 'xxx glevel is not appropriate :', glevel
       write(IO_FID_LOG,*) 'xxx glevel is not appropriate :', glevel
       call ADM_proc_stop
    endif
    if ( rlevel < 0 ) then
       write(*         ,*) 'xxx rlevel is not appropriate :', rlevel
       write(IO_FID_LOG,*) 'xxx rlevel is not appropriate :', rlevel
       call ADM_proc_stop
    endif

    ! Error if glevel & rlevel are different from global parameter
    if ( GLOBAL_ensemble_on ) then
       if ( glevel /= GLOBAL_glevel ) then
          write(*         ,*) 'xxx glevel is not equal (global,local) :', GLOBAL_glevel, glevel
          write(IO_FID_LOG,*) 'xxx glevel is not equal (global,local) :', GLOBAL_glevel, glevel
          call ADM_proc_stop
       endif
       if ( rlevel /= GLOBAL_rlevel ) then
          write(*         ,*) 'xxx rlevel is not equal (global,local) :', GLOBAL_rlevel, rlevel
          write(IO_FID_LOG,*) 'xxx rlevel is not equal (global,local) :', GLOBAL_rlevel, rlevel
          call ADM_proc_stop
       endif
    endif

    if ( ADM_HGRID_SYSTEM == 'ICO' ) then
       ADM_vlink  = 5
       dmd        = 10
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'LCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 4 * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'MLCP-OLD' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = 2 * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'MLCP' ) then
       if( ADM_vlink == -1 ) ADM_vlink = 6
       dmd        = (1+ADM_XTMS_MLCP_S)  * ADM_vlink
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'PERIODIC-1DMD' ) then ! T.Ohno 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'PARALLELOGRAM' ) then ! T.Ohno 250418
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' ) then ! M.Hara 110721
       ADM_vlink  = 5
       dmd        = 1
       ADM_prc_pl = -999
    elseif( ADM_HGRID_SYSTEM == 'ICO-XTMS' ) then
       ADM_vlink  = 5
       dmd        = 10
       ADM_prc_pl = 1
    elseif( ADM_HGRID_SYSTEM == 'SCM' ) then
       call setup_for_scm

       return
    else
       write(*         ,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       write(IO_FID_LOG,*) 'xxx Name of ADM_HGRID_SYSTEM is wrong. STOP.'
       call ADM_proc_stop
    endif

    ADM_gall_pl = ADM_vlink + 1
    ADM_gmax_pl = ADM_vlink + 1

    ADM_glevel   = glevel
    ADM_rlevel   = rlevel
    ADM_vlayer   = vlayer
    ADM_DMD      = dmd

    ADM_rgn_nmax = 2**ADM_rlevel * 2**ADM_rlevel * ADM_DMD
    ADM_lall     = ADM_rgn_nmax / ADM_prc_all

    nmax         = 2**(ADM_glevel-ADM_rlevel)
    ADM_gall_1d  = 1 + nmax + 1
    ADM_gmin     = 1 + 1
    ADM_gmax     = 1 + nmax

    ADM_gall     = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in  = (   nmax+1 ) * (   nmax+1 )

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif

    if ( ADM_prc_me == ADM_prc_master ) then
       ADM_IsMaster = .true.
    else
       ADM_IsMaster = .false.
    endif



    call RGNMNG_setup( rgnmngfname )



    allocate( GLOBAL_extension_rgn(ADM_lall) )
    allocate( ADM_have_sgp        (ADM_lall) )
    GLOBAL_extension_rgn(:) = ''
    ADM_have_sgp        (:) = .false.

    do l = 1, ADM_lall
       rgnid = RGNMNG_lp2r(l,ADM_prc_me)
       if (ADM_rlevel > 6 ) then
          write(GLOBAL_extension_rgn(l),'(A,I6.6)') '.rgn', rgnid-1
       else
          write(GLOBAL_extension_rgn(l),'(A,I5.5)') '.rgn', rgnid-1
       end if

       if ( RGNMNG_vert_num(I_W,rgnid) == 3 ) then
          ADM_have_sgp(l) = .true.
       endif

       if (       ADM_HGRID_SYSTEM == '1DMD-ON-SPHERE' &
            .AND. rgnid == 1                           ) then ! [mod] 20210322 Junya D-NICAM
          ADM_have_sgp(l) = .true.
       endif
    enddo

    if ( ADM_prc_me == ADM_prc_pl ) then
       ADM_have_pl = .true.
    else
       ADM_have_pl = .false.
    endif

    ADM_l_me = 0

    !--- make suffix for list-vector loop.
    call ADM_mk_suffix

    call output_info

    return
    !---------------------------------------------------------------------------
  contains
    !---------------------------------------------------------------------------
    subroutine setup_for_scm
      implicit none
      !-------------------------------------------------------------------------

      !Tentative (T.Ohno)
      ADM_flg_scm      = .true.
      ADM_idx_aid4stat = 1
      ADM_prc_pl       = -999
      dmd              = 1

      ADM_gall_pl  = 1

      ADM_glevel   = glevel
      ADM_rlevel   = rlevel
      ADM_vlayer   = vlayer
      ADM_DMD      = dmd

      ADM_rgn_nmax = 1
      ADM_lall     = 1

      ADM_gall_1d  = 1
      ADM_gmin     = 1
      ADM_gmax     = 0

      ADM_gall     = 1
      ADM_gall_in  = 1

      if ( ADM_vlayer == 1 ) then
         ADM_kall = 1
         ADM_kmin = 1
         ADM_kmax = 1
      else
         ADM_kall = 1 + ADM_vlayer + 1
         ADM_kmin = 1 + 1
         ADM_kmax = 1 + ADM_vlayer
      endif

      ADM_IsMaster = .true.

      allocate( RGNMNG_edge_tab(I_RGNID:I_DIR,I_SW:I_SE,ADM_rgn_nmax) )
      RGNMNG_edge_tab(:,:,:) = -1

      allocate( RGNMNG_lp2r(ADM_lall,ADM_prc_all) )
      allocate( RGNMNG_l2r (ADM_lall)             )
      RGNMNG_lp2r(1,1) = 1
      RGNMNG_l2r (1)   = 1

      allocate( GLOBAL_extension_rgn(ADM_lall) )
      GLOBAL_extension_rgn(:) = ''

      l = 1
      rgnid = RGNMNG_lp2r(l,ADM_prc_me)
      write(GLOBAL_extension_rgn(l),'(A,I5.5)') '.rgn', rgnid-1

      ADM_have_pl = .false.
      ADM_l_me    = 0

      ADM_IopJop_nmax = 1
      allocate( ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax) )
      ADM_IopJop = 1

      allocate( ADM_have_sgp(ADM_lall) )
      ADM_have_sgp(1) = .false.

      return
    end subroutine setup_for_scm

  end subroutine ADM_member_setup

  !-----------------------------------------------------------------------------
  subroutine ADM_mk_suffix
    implicit none

    integer  :: isize ! i-direction grid size without halo
    integer  :: jsize ! i-direction grid size without halo
    integer  :: n

    integer  :: i, j, suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    isize = ADM_gmax - ADM_gmin + 1
    jsize = ADM_gmax - ADM_gmin + 1

    ADM_IooJoo_nmax = ( isize     ) * ( jsize     )
    ADM_IooJmo_nmax = ( isize     ) * ( jsize + 1 )
    ADM_IooJop_nmax = ( isize     ) * ( jsize + 1 )
    ADM_IooJmp_nmax = ( isize     ) * ( jsize + 2 )
    ADM_ImoJoo_nmax = ( isize + 1 ) * ( jsize     )
    ADM_ImoJmo_nmax = ( isize + 1 ) * ( jsize + 1 )
    ADM_ImoJop_nmax = ( isize + 1 ) * ( jsize + 1 )
    ADM_ImoJmp_nmax = ( isize + 1 ) * ( jsize + 2 )
    ADM_IopJoo_nmax = ( isize + 1 ) * ( jsize     )
    ADM_IopJmo_nmax = ( isize + 1 ) * ( jsize + 1 )
    ADM_IopJop_nmax = ( isize + 1 ) * ( jsize + 1 )
    ADM_IopJmp_nmax = ( isize + 1 ) * ( jsize + 2 )
    ADM_ImpJoo_nmax = ( isize + 2 ) * ( jsize     )
    ADM_ImpJmo_nmax = ( isize + 2 ) * ( jsize + 1 )
    ADM_ImpJop_nmax = ( isize + 2 ) * ( jsize + 1 )
    ADM_ImpJmp_nmax = ( isize + 2 ) * ( jsize + 2 )

    allocate( ADM_IooJoo(ADM_IooJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJmo(ADM_IooJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJop(ADM_IooJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IooJmp(ADM_IooJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJoo(ADM_ImoJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJmo(ADM_ImoJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJop(ADM_ImoJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImoJmp(ADM_ImoJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJoo(ADM_IopJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJmo(ADM_IopJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJop(ADM_IopJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_IopJmp(ADM_IopJmp_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJoo(ADM_ImpJoo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJmo(ADM_ImpJmo_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJop(ADM_ImpJop_nmax,ADM_GIJ_nmax) )
    allocate( ADM_ImpJmp(ADM_ImpJmp_nmax,ADM_GIJ_nmax) )

    !--- ADM_IooJoo
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ADM_IooJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin  , ADM_gmax
       ADM_IooJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJop
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax
       ADM_IooJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IooJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin  , ADM_gmax
       ADM_IooJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IooJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IooJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IooJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IooJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IooJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IooJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJoo
    n = 1
    do j = ADM_gmin  , ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJop
    n = 1
    do j = ADM_gmin  , ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImoJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax
       ADM_ImoJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImoJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImoJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImoJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImoJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImoJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImoJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJoo
    n = 1
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin  , ADM_gmax+1
       ADM_IopJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJop
    n = 1
    do j = ADM_gmin, ADM_gmax+1
    do i = ADM_gmin, ADM_gmax+1
       ADM_IopJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_IopJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin  , ADM_gmax+1
       ADM_IopJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_IopJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_IopJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_IopJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_IopJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_IopJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_IopJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJoo
    n = 1
    do j = ADM_gmin  , ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJoo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJoo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJoo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJoo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJoo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJoo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJoo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmo
    n = 1
    do j = ADM_gmin-1, ADM_gmax
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmo(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmo(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmo(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmo(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmo(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmo(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmo(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJop
    n = 1
    do j = ADM_gmin  , ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJop(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJop(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJop(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJop(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJop(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJop(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJop(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    !--- ADM_ImpJmp
    n = 1
    do j = ADM_gmin-1, ADM_gmax+1
    do i = ADM_gmin-1, ADM_gmax+1
       ADM_ImpJmp(n,ADM_GIoJo) = suf(i  ,j  )
       ADM_ImpJmp(n,ADM_GIpJo) = suf(i+1,j  )
       ADM_ImpJmp(n,ADM_GIpJp) = suf(i+1,j+1)
       ADM_ImpJmp(n,ADM_GIoJp) = suf(i  ,j+1)
       ADM_ImpJmp(n,ADM_GImJo) = suf(i-1,j  )
       ADM_ImpJmp(n,ADM_GImJm) = suf(i-1,j-1)
       ADM_ImpJmp(n,ADM_GIoJm) = suf(i  ,j-1)
       n = n + 1
    enddo
    enddo

    return
  end subroutine ADM_mk_suffix

  !-----------------------------------------------------------------------------
  subroutine RGNMNG_setup( &
       rgnmngfname )
    implicit none

    character(len=*), intent(in) :: rgnmngfname

    namelist / RGNMNGPARAM / &
       RGNMNG_in_fname,  &
       RGNMNG_out_fname, &
       debug

    integer  :: l, p
    integer  :: ierr
    !---------------------------------------------------------------------------

    RGNMNG_in_fname = trim(rgnmngfname)

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[rgnmng]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=RGNMNGPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** RGNMNGPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist RGNMNGPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=RGNMNGPARAM)

    if ( ADM_lall > RGNMNG_llim ) then
       write(*,*) 'xxx limit exceed! local region: ', ADM_lall, RGNMNG_llim
       call ADM_proc_stop
    endif

    ! Global information (Each process has all the information)
    allocate( RGNMNG_edge_tab(I_RGNID:I_DIR,I_SW:I_SE,ADM_rgn_nmax) )
    allocate( RGNMNG_lnum    (ADM_prc_all) )
    allocate( RGNMNG_lp2r    (ADM_lall,ADM_prc_all) )

    if ( RGNMNG_in_fname /= '' ) then
       call RGNMNG_input( RGNMNG_in_fname,        & ! [IN]
                          ADM_rgn_nmax,           & ! [IN]
                          ADM_prc_all,            & ! [IN]
                          ADM_lall,               & ! [IN]
                          RGNMNG_edge_tab(:,:,:), & ! [OUT]
                          RGNMNG_lnum    (:),     & ! [OUT]
                          RGNMNG_lp2r    (:,:)    ) ! [OUT]
    else
       write(IO_FID_LOG,*) '*** input file is not specified.'

       if ( ADM_HGRID_SYSTEM == 'ICO' ) then
          call RGNMNG_generate_ico( ADM_rlevel,             & ! [IN]
                                    ADM_rgn_nmax,           & ! [IN]
                                    ADM_prc_all,            & ! [IN]
                                    ADM_lall,               & ! [IN]
                                    RGNMNG_edge_tab(:,:,:), & ! [OUT]
                                    RGNMNG_lnum    (:),     & ! [OUT]
                                    RGNMNG_lp2r    (:,:)    ) ! [OUT]
       else
          write(*,*) 'xxx auto generation of mnginfo is not supported for ADM_HGRID_SYSTEM = ', ADM_HGRID_SYSTEM
          call ADM_proc_stop
       endif
    endif

    if ( RGNMNG_out_fname /= '' ) then
       call RGNMNG_output( RGNMNG_out_fname,       & ! [IN]
                           ADM_rgn_nmax,           & ! [IN]
                           ADM_prc_all,            & ! [IN]
                           ADM_lall,               & ! [IN]
                           RGNMNG_edge_tab(:,:,:), & ! [IN]
                           RGNMNG_lnum    (:),     & ! [IN]
                           RGNMNG_lp2r    (:,:)    ) ! [IN]
    endif

    !--- additional table (reversal,local)
    allocate( RGNMNG_r2lp(2,ADM_rgn_nmax) )
    allocate( RGNMNG_l2r (ADM_lall) )

    do p = 1, ADM_prc_all
    do l = 1, ADM_lall
       RGNMNG_r2lp(I_l,  RGNMNG_lp2r(l,p)) = l
       RGNMNG_r2lp(I_prc,RGNMNG_lp2r(l,p)) = p
    enddo
    enddo

    do l = 1, ADM_lall
       RGNMNG_l2r(l) = RGNMNG_lp2r(l,ADM_prc_me)
    enddo

    !--- region connection chains around the diamond vertexes
    allocate( RGNMNG_vert_num   ( I_W:I_S,ADM_rgn_nmax ) )
    allocate( RGNMNG_vert_tab   ( I_RGNID:I_DIR,   &
                                  I_W:I_S,         &
                                  ADM_rgn_nmax,    &
                                  ADM_vlink        ) )
    allocate( RGNMNG_vert_tab_pl( I_RGNID:I_DIR,   &
                                  ADM_rgn_nmax_pl, &
                                  ADM_vlink        ) )

    call RGNMNG_vertex_walkaround( ADM_rgn_nmax,             & ! [IN]
                                   ADM_rgn_nmax_pl,          & ! [IN]
                                   ADM_vlink,                & ! [IN]
                                   RGNMNG_edge_tab(:,:,:),   & ! [IN]
                                   RGNMNG_vert_num(:,:),     & ! [OUT]
                                   RGNMNG_vert_tab(:,:,:,:), & ! [OUT]
                                   RGNMNG_vert_tab_pl(:,:,:) ) ! [OUT]

    !--- tables for pole
    RGNMNG_rgn4pl(I_NPL) = RGNMNG_vert_tab_pl(I_RGNID,I_NPL,ADM_vlink)
    RGNMNG_rgn4pl(I_SPL) = RGNMNG_vert_tab_pl(I_RGNID,I_SPL,ADM_vlink)

    RGNMNG_r2p_pl(I_NPL) = ADM_prc_pl
    RGNMNG_r2p_pl(I_SPL) = ADM_prc_pl

    return
  end subroutine RGNMNG_setup

  !-----------------------------------------------------------------------------
  !> Input mnginfo file
  subroutine RGNMNG_input( &
       in_fname, &
       rall,     &
       pall,     &
       lall,     &
       edge_tab, &
       lnum,     &
       lp2r      )
    implicit none

    character(len=*), intent(in)  :: in_fname           !< input file
    integer,          intent(in)  :: rall               !< total number of region
    integer,          intent(in)  :: pall               !< total number of process
    integer,          intent(in)  :: lall               !< limit number of local region
    integer,          intent(out) :: edge_tab(2,4,rall) !< region link information (for 4 edges)
    integer,          intent(out) :: lnum(pall)         !< number of local region
    integer,          intent(out) :: lp2r(lall,pall)    !< l,prc => region

    integer  :: num_of_rgn !< number of region

    namelist / rgn_info / &
       num_of_rgn

    integer  :: rgnid      !< region ID
    integer  :: sw(2) = -1 !< south-west region info
    integer  :: nw(2) = -1 !< north-west region info
    integer  :: ne(2) = -1 !< north-east region info
    integer  :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
       rgnid, &
       sw,    &
       nw,    &
       ne,    &
       se

    integer  :: num_of_proc !< number of processes

    namelist / proc_info / &
       num_of_proc

    integer  :: peid            !< process ID
    integer  :: num_of_mng      !< number of regions be managed
    integer  :: mng_rgnid(RGNMNG_llim) !< managed region ID

    namelist / rgn_mng_info / &
       peid,       &
       num_of_mng, &
       mng_rgnid

    integer  :: fid, ierr
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** input mnginfo. filename = ',trim(in_fname)

    fid = IO_get_available_fid()
    open( unit   = fid,            &
          file   = trim(in_fname), &
          form   = 'formatted',    &
          status = 'old',          &
          iostat = ierr            )

       ! ERROR if filename are not defined
       if ( ierr /= 0 ) then
          write(IO_FID_LOG,*) 'xxx mnginfo file is not found! STOP. ', trim(in_fname)
          call ADM_proc_stop
       endif

       read(fid,nml=rgn_info)

       if ( num_of_rgn /= rall ) then
          write(IO_FID_LOG,*) 'xxx No match for region number! STOP.'
          write(IO_FID_LOG,*) 'xxx rall= ', rall, ' num_of_rgn=', num_of_rgn
          call ADM_proc_stop
       endif

       do r = 1, rall
          read(fid,nml=rgn_link_info)

          edge_tab(I_RGNID:I_DIR,I_SW,rgnid) = sw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NW,rgnid) = nw(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_NE,rgnid) = ne(I_RGNID:I_DIR)
          edge_tab(I_RGNID:I_DIR,I_SE,rgnid) = se(I_RGNID:I_DIR)
       enddo

       read(fid,nml=proc_info)
       if ( num_of_proc /= pall ) then
          write(IO_FID_LOG,*) ' xxx No match for  process number! STOP.'
          write(IO_FID_LOG,*) ' xxx prc_all= ', pall, ' num_of_proc=', num_of_proc
          call ADM_proc_stop
       endif

       do p = 1, pall
          read(fid,nml=rgn_mng_info)

          lnum(peid) = num_of_mng

          if ( lnum(peid) /= ADM_lall ) then
             write(*,*) 'xxx local region number is not match: ', peid, lnum(peid), ADM_lall
             call ADM_proc_stop
          endif

          lp2r(:,peid) = -1 ! initialize
          do l = 1, lnum(peid)
             lp2r(l,peid) = mng_rgnid(l)
          enddo
       enddo

    close(fid)

    return
  end subroutine RGNMNG_input

  !-----------------------------------------------------------------------------
  !> Output mnginfo file
  subroutine RGNMNG_output( &
       out_fname, &
       rall,      &
       pall,      &
       lall,      &
       edge_tab,  &
       lnum,      &
       lp2r       )
    implicit none

    character(len=*), intent(in) :: out_fname          !< output file
    integer,          intent(in) :: rall               !< total number of region
    integer,          intent(in) :: pall               !< total number of process
    integer,          intent(in) :: lall               !< limit number of local region
    integer,          intent(in) :: edge_tab(2,4,rall) !< region link information (for 4 edges)
    integer,          intent(in) :: lnum(pall)         !< number of local region
    integer,          intent(in) :: lp2r(lall,pall)    !< l,prc => region

    integer  :: num_of_rgn !< number of region

    namelist / rgn_info / &
       num_of_rgn

    integer  :: rgnid      !< region ID
    integer  :: sw(2) = -1 !< south-west region info
    integer  :: nw(2) = -1 !< north-west region info
    integer  :: ne(2) = -1 !< north-east region info
    integer  :: se(2) = -1 !< south-east region info

    namelist / rgn_link_info / &
       rgnid, &
       sw,    &
       nw,    &
       ne,    &
       se

    integer  :: num_of_proc !< number of processes

    namelist / proc_info / &
       num_of_proc

    integer  :: peid            !< process ID
    integer  :: num_of_mng      !< number of regions be managed
    integer  :: mng_rgnid(RGNMNG_llim) !< managed region ID

    namelist / rgn_mng_info / &
       peid,       &
       num_of_mng, &
       mng_rgnid

    integer  :: fid
    integer  :: r, p, l
    !---------------------------------------------------------------------------

    if ( ADM_IsMaster ) then
       write(IO_FID_LOG,*) '*** output mnginfo: ',trim(out_fname)

       fid = IO_get_available_fid()
       open( unit   = fid,             &
             file   = trim(out_fname), &
             form   = 'formatted'      )

          num_of_rgn = rall
          write(fid,nml=rgn_info)

          do r = 1, rall
             rgnid = r
             sw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SW,rgnid)
             nw(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NW,rgnid)
             ne(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_NE,rgnid)
             se(I_RGNID:I_DIR) = edge_tab(I_RGNID:I_DIR,I_SE,rgnid)

             write(fid,nml=rgn_link_info)
          enddo

          num_of_proc = pall
          write(fid,nml=proc_info)

          do p = 1, pall
             peid = p
             num_of_mng = lnum(peid)

             mng_rgnid(:) = -1
             do l = 1, lnum(peid)
                mng_rgnid(l) = lp2r(l,peid)
             enddo

             write(fid,nml=rgn_mng_info)
          enddo

       close(fid)
    else
       write(IO_FID_LOG,*) '*** output mnginfo: only in the master process'
    endif

    return
  end subroutine RGNMNG_output

  !-----------------------------------------------------------------------------
  !> Generate region management info
  subroutine RGNMNG_generate_ico( &
       rlevel,   &
       rall,     &
       pall,     &
       lall,     &
       edge_tab, &
       lnum,     &
       lp2r      )
    implicit none

    integer, intent(in)  :: rlevel
    integer, intent(in)  :: rall
    integer, intent(in)  :: pall
    integer, intent(in)  :: lall
    integer, intent(out) :: edge_tab(2,4,rall)
    integer, intent(out) :: lnum(pall)
    integer, intent(out) :: lp2r(lall,pall)

    integer, parameter :: nmax_dmd = 10
    integer  :: dmd_data(4,nmax_dmd)

    integer  :: rall_1d, rall_1dmd

    integer  :: d_nb, i_nb, j_nb, rgnid_nb, direction
    integer  :: d, i, j, rgnid
    integer  :: l, p
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** generate region management table'
    write(IO_FID_LOG,*) '*** Topology: default icosahedral'

    dmd_data(:, 1) = (/  6, 5, 2,10 /)
    dmd_data(:, 2) = (/ 10, 1, 3, 9 /)
    dmd_data(:, 3) = (/  9, 2, 4, 8 /)
    dmd_data(:, 4) = (/  8, 3, 5, 7 /)
    dmd_data(:, 5) = (/  7, 4, 1, 6 /)
    dmd_data(:, 6) = (/  7, 5, 1,10 /)
    dmd_data(:, 7) = (/  8, 4, 5, 6 /)
    dmd_data(:, 8) = (/  9, 3, 4, 7 /)
    dmd_data(:, 9) = (/ 10, 2, 3, 8 /)
    dmd_data(:,10) = (/  6, 1, 2, 9 /)

    rall_1d   = 2**rlevel
    rall_1dmd = rall_1d*rall_1d

    !--- make region link table
    do d = 1, nmax_dmd
    do i = 1, rall_1d
    do j = 1, rall_1d
       rgnid = (d-1)*rall_1dmd + (j-1)*rall_1d + i

       !--- I_SW
       if ( j == 1 ) then
          if ( d <= 5 ) then
             i_nb = i
             j_nb = rall_1d
             d_nb = dmd_data(I_SW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_SW,d)
             direction = I_SE
          endif
       else
          i_nb = i
          j_nb = j-1
          d_nb = d
          direction = I_NE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SW,rgnid) = direction

       !--- I_NW
       if ( i == 1 ) then
          if ( d <= 5 ) then
             i_nb = rall_1d+1-j
             j_nb = rall_1d
             d_nb = dmd_data(I_NW,d)
             direction = I_NE
          else
             i_nb = rall_1d
             j_nb = j
             d_nb = dmd_data(I_NW,d)
             direction = I_SE
          endif
       else
          i_nb = i-1
          j_nb = j
          d_nb = d
          direction = I_SE
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NW,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NW,rgnid) = direction

       !--- I_NE
       if ( j == rall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = rall_1d+1-i
             d_nb = dmd_data(I_NE,d)
             direction = I_NW
          else
             i_nb = i
             j_nb = 1
             d_nb = dmd_data(I_NE,d)
             direction = I_SW
          endif
       else
          i_nb = i
          j_nb = j+1
          d_nb = d
          direction = I_SW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_NE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_NE,rgnid) = direction

       !--- I_SE
       if ( i == rall_1d ) then
          if ( d <= 5 ) then
             i_nb = 1
             j_nb = j
             d_nb = dmd_data(I_SE,d)
             direction = I_NW
          else
             i_nb = rall_1d+1-j
             j_nb = 1
             d_nb = dmd_data(I_SE,d)
             direction = I_SW
          endif
       else
          i_nb = i+1
          j_nb = j
          d_nb = d
          direction = I_NW
       endif
       rgnid_nb = (d_nb-1)*rall_1dmd + (j_nb-1)*rall_1d + i_nb

       edge_tab(I_RGNID,I_SE,rgnid) = rgnid_nb
       edge_tab(I_DIR,  I_SE,rgnid) = direction

    enddo
    enddo
    enddo

    if ( mod(rall,pall) /= 0 ) then
       write(*,*) 'invalid number of process!', rall, pall
       call ADM_proc_stop
    else
       if ( rall/pall /= lall ) then
          write(*,*) 'xxx local region number is not match: ', rall/pall, lall
          call ADM_proc_stop
       endif

       do p = 1, pall
          lnum(p) = rall / pall
       enddo
    endif

    !--- make region-pe relationship
    lp2r(:,:) = -1
    rgnid = 0

    do p = 1, pall
    do l = 1, lnum(p)
       rgnid = rgnid + 1

       lp2r(l,p) = rgnid
    enddo
    enddo

    return
  end subroutine RGNMNG_generate_ico

  !-----------------------------------------------------------------------------
  !> Search region with walking around the diamond vertexes
  subroutine RGNMNG_vertex_walkaround( &
       rall,       &
       rall_pl,    &
       vlink_nmax, &
       edge_tab,   &
       vert_num,   &
       vert_tab,   &
       vert_tab_pl )
    implicit none

    integer, intent(in)  :: rall
    integer, intent(in)  :: rall_pl
    integer, intent(in)  :: vlink_nmax
    integer, intent(in)  :: edge_tab(2,4,rall)

    integer, intent(out) :: vert_num   (4,rall)
    integer, intent(out) :: vert_tab   (2,4,rall,vlink_nmax)
    integer, intent(out) :: vert_tab_pl(2,rall_pl,vlink_nmax)

    integer  :: rgnid, dir, dir0
    integer  :: rgnid_next, dir_next

    integer  :: r, d, v
    !---------------------------------------------------------------------------

    vert_num(:,:)     = -1
    vert_tab(:,:,:,:) = -1

    do r = 1, rall
    do d = I_W, I_S

       rgnid = r
       select case(d)
       case(I_W)
          dir = I_SW
       case(I_N)
          dir = I_NW
       case(I_E)
          dir = I_NE
       case(I_S)
          dir = I_SE
       endselect
       dir0 = dir

       v = 0
       do
          rgnid_next = edge_tab(I_RGNID,dir,rgnid)
          dir_next   = edge_tab(I_DIR,  dir,rgnid) - 1

          if( dir_next == 0 ) dir_next = 4
          v = v + 1
          vert_tab(I_RGNID,d,r,v) = rgnid
          vert_tab(I_DIR,  d,r,v) = dir

          rgnid = rgnid_next
          dir   = dir_next

          if( rgnid == r .and. dir == dir0 ) exit
       enddo
       vert_num(d,r) = v

    enddo
    enddo

    vert_tab_pl(:,:,:) = -1

    do r = rall, 1, -1
       if ( vert_num(I_N,r) == vlink_nmax ) then
          do v = 1, vlink_nmax
             vert_tab_pl(I_RGNID,I_NPL,v) = vert_tab(I_RGNID,I_N,r,v)
             vert_tab_pl(I_DIR,  I_NPL,v) = vert_tab(I_DIR,  I_N,r,v)
          enddo
          exit
       endif
    enddo

    do r = rall, 1, -1
       if ( vert_num(I_S,r) == vlink_nmax ) then
          do v = 1, vlink_nmax
             vert_tab_pl(I_RGNID,I_SPL,v) = vert_tab(I_RGNID,I_S,r,v)
             vert_tab_pl(I_DIR,  I_SPL,v) = vert_tab(I_DIR,  I_S,r,v)
          enddo
          exit
       endif
    enddo

    return
  end subroutine RGNMNG_vertex_walkaround

  !-----------------------------------------------------------------------------
  subroutine output_info
    implicit none

    integer          :: rgnid, rgnid_next
    character(len=2) :: dstr, dstr_next

    integer  :: l, d, v
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Process management info. ======'
    write(IO_FID_LOG,'(1x,A,I10)') '--- Total number of process            : ', ADM_prc_all
    write(IO_FID_LOG,'(1x,A,I10)') '--- Communicator ID                    : ', ADM_COMM_WORLD
    write(IO_FID_LOG,'(1x,A,I10)') '--- My Process number = (my rank + 1)  : ', ADM_prc_me
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Region/Grid topology info. ======'
    write(IO_FID_LOG,'(1x,A,A10)') '--- Grid sysytem                       : ', trim(ADM_HGRID_SYSTEM)
    write(IO_FID_LOG,'(1x,A,I10)') '--- #  of diamond                      : ', ADM_DMD
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Region management info. ======'
    write(IO_FID_LOG,'(1x,A,I10)') '--- Region level (RL)                  : ', ADM_rlevel
    write(IO_FID_LOG,'(1x,A,I10,3(A,I4),A)') &
                                   '--- Total number of region             : ', ADM_rgn_nmax, &
                                        ' (', 2**ADM_rlevel, ' x', 2**ADM_rlevel, ' x', ADM_DMD, ' )'
    write(IO_FID_LOG,'(1x,A,I10)') '--- #  of region per process           : ', ADM_lall
    write(IO_FID_LOG,'(1x,A)'    ) '--- ID of region in my process         : '
    write(IO_FID_LOG,*)           RGNMNG_lp2r(1:ADM_lall,ADM_prc_me)

    write(IO_FID_LOG,'(1x,A,I10)') '--- Region ID, contains north pole     : ', RGNMNG_rgn4pl(I_NPL)
    write(IO_FID_LOG,'(1x,A,I10)') '--- Region ID, contains south pole     : ', RGNMNG_rgn4pl(I_SPL)
    write(IO_FID_LOG,'(1x,A,I10)') '--- Process rank, managing north pole  : ', RGNMNG_r2p_pl(I_NPL)
    write(IO_FID_LOG,'(1x,A,I10)') '--- Process rank, managing south pole  : ', RGNMNG_r2p_pl(I_SPL)
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== Grid management info. ======'
    write(IO_FID_LOG,'(1x,A,I10)') '--- Grid level (GL)                    : ', ADM_glevel
    write(IO_FID_LOG,'(1x,A,I10,2(A,I4),A,I7,A)') &
                                   '--- Total number of grid (horizontal)  : ',  &
                                        4**(ADM_glevel-ADM_rlevel)*ADM_rgn_nmax, &
                                        ' (', 2**(ADM_glevel-ADM_rlevel),        &
                                        ' x', 2**(ADM_glevel-ADM_rlevel),        &
                                        ' x', ADM_rgn_nmax, ' )'
    write(IO_FID_LOG,'(1x,A,I10)') '--- Number of vertical layer           : ', ADM_kmax-ADM_kmin+1

    if ( debug ) then
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '====== region management information ======'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,'(1x,A,I7)') '--- # of region in this node : ', ADM_lall

       write(IO_FID_LOG,*) '--- (l,prc_me) => (rgn)'
       do l = 1, ADM_lall
          rgnid = RGNMNG_l2r(l)
          write(IO_FID_LOG,'(1x,A,I4,A,I6,A,I6,A)') '--- (',l,',',ADM_prc_me,') => (',rgnid,') '
       enddo

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- Link information'
       do l = 1, ADM_lall
          rgnid = RGNMNG_l2r(l)

          write(IO_FID_LOG,*)
          write(IO_FID_LOG,*) '--- edge link: (rgn,direction)'
          do d = I_SW, I_SE
             rgnid_next = RGNMNG_edge_tab(I_RGNID,d,rgnid)
             dstr       = RGNMNG_edgename(d)
             dstr_next  = RGNMNG_edgename(RGNMNG_edge_tab(I_DIR,d,rgnid))
             write(IO_FID_LOG,'(5x,A,I6,A,A,A,I6,A,A,A)') '(',rgnid,',',dstr,') -> (', rgnid_next,',', dstr_next,')'
          enddo

          write(IO_FID_LOG,*) '--- vertex link: (rgn)'
          do d = I_W, I_S
             dstr = RGNMNG_vertname(d)
             write(IO_FID_LOG,'(5x,A,I6,A,A,A)',advance='no') '(',rgnid,',',dstr,')'
             do v = 2, RGNMNG_vert_num(d,rgnid)
                dstr = RGNMNG_vertname(RGNMNG_vert_tab(I_DIR,d,rgnid,v))
                write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',RGNMNG_vert_tab(I_RGNID,d,rgnid,v),',',dstr,')'
             enddo
             write(IO_FID_LOG,*)
          enddo
       enddo

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- Pole information (in the global scope)'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- region, having north pole data : ', RGNMNG_rgn4pl(I_NPL)
       write(IO_FID_LOG,*) '--- vertex link: (north pole)'
       do v = 1, ADM_vlink
          rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_NPL,v)
          dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_NPL,v))
          write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       enddo
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- process, managing north pole : ', RGNMNG_r2p_pl(I_NPL)

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- region, having south pole data : ', RGNMNG_rgn4pl(I_SPL)
       write(IO_FID_LOG,*) '--- vertex link: (south pole)'
       do v = 1, ADM_vlink
          rgnid = RGNMNG_vert_tab_pl(I_RGNID,I_SPL,v)
          dstr  = RGNMNG_vertname(RGNMNG_vert_tab_pl(I_DIR,I_SPL,v))
          write(IO_FID_LOG,'(A,I6,A,A,A)',advance='no') ' -> (',rgnid,',',dstr,')'
       enddo
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '--- process, managing south pole : ', RGNMNG_r2p_pl(I_SPL)
    endif

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  subroutine ADM_ensemble_setup
    implicit none

    integer  :: d_1d, di, dj, gdiv_1d
    integer  :: gss, ges, gse, gee
    integer  :: d, l

    integer  :: i,j,suf
    suf(i,j) = ADM_gall_1d * (j-1) + i
    !---------------------------------------------------------------------------

    if ( GLOBAL_ensemble_on ) then

       do l = 1, ADM_lall
          if (      ADM_prc_me == RGNMNG_r2p_pl(I_NPL) &
               .OR. ADM_prc_me == RGNMNG_r2p_pl(I_SPL) ) then
             ENS_have_pl = .true.
          else
             ENS_have_pl = .false.
          endif
       enddo

       ! management of spatial division
       ENS_elevel   = GLOBAL_elevel
       ENS_div_nmax = 2**ENS_elevel * 2**ENS_elevel * ADM_lall

       if ( mod(ENS_div_nmax,ENS_prc_all) /= 0 ) then
          write(*         ,*) 'xxx ENS_div_nmax must be divided by ENS_prc_all :', ENS_div_nmax, ENS_prc_all
          write(IO_FID_LOG,*) 'xxx ENS_div_nmax must be divided by ENS_prc_all :', ENS_div_nmax, ENS_prc_all
          call ADM_proc_stop
       endif

       ENS_eall = ENS_prc_all
       ENS_ediv = ENS_div_nmax / ENS_eall

       allocate( ENS_d2e(ENS_div_nmax) )
       allocate( ENS_d2l(ENS_div_nmax) )

       do d = 1, ENS_div_nmax
          ENS_d2e(d) = (d-1) / ENS_ediv + 1
          ENS_d2l(d) = (d-1) / (4**ENS_elevel) + 1
       enddo

       ! management of grid for spatial division
       allocate( ENS_imin(ENS_div_nmax) )
       allocate( ENS_imax(ENS_div_nmax) )
       allocate( ENS_jmin(ENS_div_nmax) )
       allocate( ENS_jmax(ENS_div_nmax) )

       d_1d    = 2**ENS_elevel
       gdiv_1d = 2**(ADM_glevel-ADM_rlevel-ENS_elevel)

       ! gdiv have halo
       ENS_gdiv_1d = gdiv_1d+1
       ENS_gdiv    = (gdiv_1d+1) * (gdiv_1d+1)

       do d = 1, ENS_div_nmax
          di = mod((d-1)     ,d_1d) + 1
          dj = mod((d-1)/d_1d,d_1d) + 1
          ENS_imin(d) = 1 + (di-1)*gdiv_1d + 1
          ENS_imax(d) = 1 + (di-1)*gdiv_1d + gdiv_1d + 1
          ENS_jmin(d) = 1 + (dj-1)*gdiv_1d + 1
          ENS_jmax(d) = 1 + (dj-1)*gdiv_1d + gdiv_1d + 1
       enddo

       ENS_gall = ENS_ediv * ENS_gdiv

       !--- report
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,'(1x,A)'    ) '====== Process management info. (local ensemble) ======'
       write(IO_FID_LOG,'(1x,A,I11)') '--- Total number of process         : ', ENS_prc_all
       write(IO_FID_LOG,'(1x,A,I11)') '--- Communicator ID                 : ', ENS_COMM_WORLD
       write(IO_FID_LOG,'(1x,A,I11)') '--- My Process rank (+1)            : ', ENS_prc_me
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,'(1x,A)'   ) '====== Ensemble management info. ======'
       write(IO_FID_LOG,'(1x,A,I7)') '--- Ensemble level (EL)             : ', ENS_elevel
       write(IO_FID_LOG,'(1x,A,I7)') '--- Ensemble member                 : ', ENS_eall
       write(IO_FID_LOG,'(1x,A,I7)') '--- Number of horizontal grid       : ', ENS_gall
       write(IO_FID_LOG,'(1x,A,L2)') '--- This process have pole value?   : ', ENS_have_pl
!       write(IO_FID_LOG,'(1x,A,I7)') '--- Number of horizontal grid, pole : ', ENS_gall_pl
       write(IO_FID_LOG,'(1x,A,I7,2(A,I4),A)') &
                                      '--- Total number of divided part    : ', ENS_div_nmax, &
                                      ' (', d_1d, &
                                      ' x', d_1d, ' )'
       write(IO_FID_LOG,'(1x,A,I7)') '--- # of division in this node      : ', ENS_ediv
       write(IO_FID_LOG,'(1x,A,I7,2(A,I4),A)') &
                                      '--- horizontal grid per division    : ', ENS_gdiv, &
                                      ' (', ENS_gdiv_1d, &
                                      ' x', ENS_gdiv_1d, ' )'

       if ( debug ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,*) '====== division management information ======'
          write(IO_FID_LOG,*)

          write(IO_FID_LOG,*) '--- (div) => (e)'
          do d = 1, ENS_div_nmax
             write(IO_FID_LOG,'(1x,2(A,I6),A)') '--- (',d,') => (',ENS_d2e(d),')'
             gss = suf(ENS_imin(d),ENS_jmin(d))
             ges = suf(ENS_imax(d),ENS_jmin(d))
             gse = suf(ENS_imin(d),ENS_jmax(d))
             gee = suf(ENS_imax(d),ENS_jmax(d))
             write(IO_FID_LOG,'(1x,A,I4)') '    region position :', ENS_d2l(d)
             write(IO_FID_LOG,'(1x,A)')    '    grid position :'
             write(IO_FID_LOG,'(1x,A,I6,2(A,I4),A,A,I6,2(A,I4),A)') &
                                            '    ',gss,'(',ENS_imin(d),',',ENS_jmin(d),')', &
                                             ' - ',ges,'(',ENS_imax(d),',',ENS_jmin(d),')'
             write(IO_FID_LOG,'(1x,A)')    '        |           |'
             write(IO_FID_LOG,'(1x,A,I6,2(A,I4),A,A,I6,2(A,I4),A)') &
                                            '    ',gse,'(',ENS_imin(d),',',ENS_jmax(d),')', &
                                             ' - ',gee,'(',ENS_imax(d),',',ENS_jmax(d),')'
          enddo
       endif

    endif

    return
  end subroutine ADM_ensemble_setup

  !-----------------------------------------------------------------------------
  !> Barrier MPI
  subroutine ADM_MPIbarrier
    implicit none

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( ADM_mpi_alive ) then
       call MPI_barrier(GLOBAL_COMM_WORLD,ierr)
    endif

  end subroutine ADM_MPIbarrier

  !-----------------------------------------------------------------------------
  !> Get MPI time
  !> @return time
  function ADM_MPItime() result(time)
    implicit none

    real(DP) :: time
    !---------------------------------------------------------------------------

    if ( ADM_mpi_alive ) then
       time = real(MPI_WTIME(),kind=DP)
    else
       call CPU_TIME(time)
    endif

  end function ADM_MPItime

  !-----------------------------------------------------------------------------
  !> Calc global statistics for timer
  subroutine ADM_MPItimestat( &
      avgvar, &
      maxvar, &
      minvar, &
      maxidx, &
      minidx, &
      var     )
    implicit none

    real(DP), intent(out) :: avgvar(:) !< average
    real(DP), intent(out) :: maxvar(:) !< maximum
    real(DP), intent(out) :: minvar(:) !< minimum
    integer,  intent(out) :: maxidx(:) !< index of maximum
    integer,  intent(out) :: minidx(:) !< index of minimum
    real(DP), intent(in)  :: var(:)    !< values for statistics

    real(DP), allocatable :: statval(:,:)
    integer               :: vsize

    real(DP) :: totalvar
    integer  :: ierr
    integer  :: v, p
    !---------------------------------------------------------------------------

    call ADM_MPIbarrier

    vsize = size(var(:))

    allocate( statval(vsize,0:ADM_prc_all-1) )
    statval(:,:) = 0.0_DP

    if ( ADM_mpi_alive ) then

!        do p = 0, ADM_prc_all-1
!           call MPI_Bcast( statval(1,p),         &
!                           vsize,                &
!                           MPI_DOUBLE_PRECISION, &
!                           p,                    &
!                           ADM_COMM_WORLD,       &
!                           ierr                  )
!        enddo

       call MPI_Allgather( var(:),               &
                           vsize,                &
                           MPI_DOUBLE_PRECISION, &
                           statval(:,:),         &
                           vsize,                &
                           MPI_DOUBLE_PRECISION, &
                           ADM_COMM_WORLD,       &
                           ierr                  )

       do v = 1, vsize
          totalvar = 0.0_DP
          do p = 0, ADM_prc_all-1
             totalvar = totalvar + statval(v,p)
          enddo
          avgvar(v) = totalvar / real(ADM_prc_all,kind=DP)

          maxvar(v)   = maxval(statval(v,:))
          minvar(v)   = minval(statval(v,:))
          maxidx(v:v) = maxloc(statval(v,:))
          minidx(v:v) = minloc(statval(v,:))
       enddo
    else
       do v = 1, vsize
          avgvar(v)   = statval(v,ADM_prc_me-1)
          maxvar(v)   = statval(v,ADM_prc_me-1)
          minvar(v)   = statval(v,ADM_prc_me-1)
          maxidx(v:v) = ADM_prc_me-1
          minidx(v:v) = ADM_prc_me-1
       enddo
    endif

    call ADM_MPIbarrier

    deallocate( statval )

    return
  end subroutine ADM_MPItimestat

  !-----------------------------------------------------------------------------
  subroutine ADM_make_idstr( &
       fname,    &
       besename, &
       l         )
    implicit none

    character(len=*), intent(out) :: fname    ! combined string (file name)
    character(len=*), intent(in)  :: besename ! basename
    integer,          intent(in)  :: l        ! local region number
    !---------------------------------------------------------------------------

    fname =  trim(GLOBAL_prefix_dir)      &
          // trim(besename)               &
          // trim(GLOBAL_extension_ens)   &
          // trim(GLOBAL_extension_rgn(l) )

    return
  end subroutine ADM_make_idstr

  !-----------------------------------------------------------------------------
  subroutine ADM_set_prefix_extension( &
       nmem_offset )
    implicit none

    integer, intent(in)  :: nmem_offset ! offset for data
    !---------------------------------------------------------------------------

    if ( GLOBAL_prefix_dir /= '' .AND. prefix_dir /= '' ) then
       write(GLOBAL_prefix_dir ,'(A,I6.6,A)') trim(prefix_dir), ENS_prc_me-1+nmem_offset, '/'
    endif

    if ( GLOBAL_extension_ens /= '' .AND. extension_ens /= '' ) then
       write(GLOBAL_extension_ens,'(A,A,I6.6)') '.', trim(extension_ens), ENS_prc_me-1+nmem_offset
    endif

    return
  end subroutine ADM_set_prefix_extension

end module mod_adm
