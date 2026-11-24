!-------------------------------------------------------------------------------
!> module PROFILER
!!
!! @par Description
!!         Time counter & FLOP counter(PAPI) toolbox + simple array checker
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_prof
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: PROF_setup
  public :: PROF_setprefx
  public :: PROF_rapstart
  public :: PROF_rapend
  public :: PROF_rapreport

#ifdef _PAPI_
  public :: PROF_PAPI_rapstart
  public :: PROF_PAPI_rapstop
  public :: PROF_PAPI_rapreport
#endif

  public :: PROF_valcheck

  interface PROF_valcheck
     module procedure PROF_valcheck_SP_1D
     module procedure PROF_valcheck_SP_2D
     module procedure PROF_valcheck_SP_3D
     module procedure PROF_valcheck_SP_4D
     module procedure PROF_valcheck_SP_5D
     module procedure PROF_valcheck_SP_6D
     module procedure PROF_valcheck_DP_1D
     module procedure PROF_valcheck_DP_2D
     module procedure PROF_valcheck_DP_3D
     module procedure PROF_valcheck_DP_4D
     module procedure PROF_valcheck_DP_5D
     module procedure PROF_valcheck_DP_6D
     module procedure PROF_valcheck_INT_1D
     module procedure PROF_valcheck_INT_2D
     module procedure PROF_valcheck_INT_3D
     module procedure PROF_valcheck_INT_4D
     module procedure PROF_valcheck_INT_5D
     module procedure PROF_valcheck_INT_6D
  end interface PROF_valcheck

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: get_rapid
  private :: get_grpid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                  private, parameter :: PROF_rapnlimit = 500
  character(len=H_SHORT),   private            :: PROF_prefix    = ''
  integer,                  private            :: PROF_rapnmax   = 0
  character(len=H_SHORT*2), private            :: PROF_rapname (PROF_rapnlimit)
  integer,                  private            :: PROF_grpnmax   = 0
  character(len=H_SHORT),   private            :: PROF_grpname (PROF_rapnlimit)
  integer,                  private            :: PROF_grpid   (PROF_rapnlimit)
  real(DP),                 private            :: PROF_raptstr (PROF_rapnlimit)
  real(DP),                 private            :: PROF_rapttot (PROF_rapnlimit)
  integer,                  private            :: PROF_rapnstr (PROF_rapnlimit)
  integer,                  private            :: PROF_rapnend (PROF_rapnlimit)
  integer,                  private            :: PROF_raplevel(PROF_rapnlimit)

  integer,                  private, parameter :: PROF_default_rap_level = 2
  integer,                  private            :: PROF_rap_level         = 2
  logical,                  private            :: PROF_mpi_barrier       = .false.

#ifdef _PAPI_
  integer(DP),              private            :: PROF_PAPI_flops     = 0   !> total floating point operations since the first call
  real(SP),                 private            :: PROF_PAPI_real_time = 0.0 !> total realtime since the first PROF_PAPI_flops() call
  real(SP),                 private            :: PROF_PAPI_proc_time = 0.0 !> total process time since the first PROF_PAPI_flops() call
  real(SP),                 private            :: PROF_PAPI_mflops    = 0.0 !> Mflop/s achieved since the previous call
  integer,                  private            :: PROF_PAPI_check
#endif

  !--- MINIAPP: opt-in GPU fence at every rapstart/rapend. The rap timers are
  !    host walls around async(1) device work, so per-section attribution is
  !    approximate; PROF_rap_sync=.true. inserts an `acc wait` at both edges of
  !    every section -> honest per-section GPU time (totals slightly perturbed
  !    because overlap is serialized). Default off = original behavior.
  logical,                  public             :: PROF_rap_sync = .false.

  character(len=7),         private            :: PROF_header
  character(len=16),        private            :: PROF_item
  real(DP),                 private            :: PROF_max
  real(DP),                 private            :: PROF_min
  real(DP),                 private            :: PROF_sum

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine PROF_setup
    use mod_adm, only: &
       ADM_proc_stop
    implicit none

    namelist / PROFPARAM / &
       PROF_rap_level,  &
       PROF_rap_sync,   & !--- MINIAPP
       PROF_mpi_barrier

    integer  :: ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[prof]/Category[common share]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PROFPARAM,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
       write(IO_FID_LOG,*) '*** PROFPARAM is not specified. use default.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*         ,*) 'xxx Not appropriate names in namelist PROFPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist PROFPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=PROFPARAM)

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*)           '====== Time profiler info. ======'
    write(IO_FID_LOG,'(1x,A,I1)') '--- Rap output level              = ', PROF_rap_level
    write(IO_FID_LOG,'(1x,A,L1)') '--- Add MPI_barrier in every rap? = ', PROF_mpi_barrier

    PROF_prefix = ''

    return
  end subroutine PROF_setup

  !-----------------------------------------------------------------------------
  subroutine PROF_setprefx( &
       prefxname )
    implicit none

    character(len=*), intent(in) :: prefxname !< prefix
    !---------------------------------------------------------------------------

    if ( prefxname == '' ) then !--- no prefix
       PROF_prefix = ''
    else
       PROF_prefix = trim(prefxname)//'_'
    endif

    return
  end subroutine PROF_setprefx

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine PROF_rapstart( rapname_base, level )
    use mod_adm, only: &
       ADM_MPIbarrier, &
       ADM_MPItime
#ifdef _OPENACC
       use nvtx
#endif
    implicit none

    character(len=*), intent(in) :: rapname_base    !< name of item

    integer,          intent(in), optional :: level !< level of item (default is 2)

    character(len=H_SHORT*2) :: rapname             !< name of item with prefix
    integer                  :: id, level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       level_ = level
    else
       level_ = PROF_default_rap_level
    endif

    if( level_ > PROF_rap_level ) return

    rapname = trim(PROF_prefix)//trim(rapname_base)

#ifdef _OPENACC
    call nvtxStartRange(rapname)
    if ( PROF_rap_sync ) then !--- MINIAPP: honest section start
       !$acc wait
    endif
#endif

    id = get_rapid( rapname, level_ )

    if(PROF_mpi_barrier) call ADM_MPIbarrier

    PROF_raptstr(id) = ADM_MPItime()
    PROF_rapnstr(id) = PROF_rapnstr(id) + 1

#if defined DEBUG || QUICKDEBUG
    !write(IO_FID_LOG,*) '<DEBUG> [PROF] ', rapname, PROF_rapnstr(id)
#endif

#if defined DEBUG || QUICKDEBUG
    call flush(IO_FID_LOG)
#endif

#ifdef _FAPP_
    call FAPP_START( trim(rapname), id, level_ )
#endif

    return
  end subroutine PROF_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine PROF_rapend( rapname_base, level )
    use mod_adm, only: &
       ADM_MPIbarrier, &
       ADM_MPItime
#ifdef _OPENACC
       use nvtx
#endif
    implicit none

    character(len=*), intent(in) :: rapname_base    !< name of item

    integer,          intent(in), optional :: level !< level of item

    character(len=H_SHORT*2) :: rapname             !< name of item with prefix
    integer                  :: id, level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       if( level > PROF_rap_level ) return
    endif

#ifdef _OPENACC
    call nvtxEndRange()
    if ( PROF_rap_sync ) then !--- MINIAPP: honest section end
       !$acc wait
    endif
#endif

    rapname = trim(PROF_prefix)//trim(rapname_base)

    id = get_rapid( rapname, level_ )

    if( level_ > PROF_rap_level ) return

    if(PROF_mpi_barrier) call ADM_MPIbarrier

    PROF_rapttot(id) = PROF_rapttot(id) + ( ADM_MPItime()-PROF_raptstr(id) )
    PROF_rapnend(id) = PROF_rapnend(id) + 1

#ifdef _FAPP_
    call FAPP_STOP( trim(rapname), id, level_ )
#endif

    return
  end subroutine PROF_rapend

  !-----------------------------------------------------------------------------
  !> Report raptime
  subroutine PROF_rapreport
    use mod_adm, only: &
       ADM_MPItimestat, &
       ADM_IsMaster
    implicit none

    real(DP) :: avgvar(PROF_rapnlimit)
    real(DP) :: maxvar(PROF_rapnlimit)
    real(DP) :: minvar(PROF_rapnlimit)
    integer  :: maxidx(PROF_rapnlimit)
    integer  :: minidx(PROF_rapnlimit)

    integer  :: id, gid
    integer  :: fid
    !---------------------------------------------------------------------------

    do id = 1, PROF_rapnmax
       if ( PROF_rapnstr(id) /= PROF_rapnend(id) ) then
           write(*,*) '*** Mismatch Report',id,PROF_rapname(id),PROF_rapnstr(id),PROF_rapnend(id)
       endif
    enddo

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Computational Time Report'
    write(IO_FID_LOG,*) '*** Rap level is ', PROF_rap_level

    if ( .false. ) then ! report for each node

       do gid = 1, PROF_rapnmax
       do id  = 1, PROF_rapnmax
          if (       PROF_raplevel(id) <= PROF_rap_level &
               .AND. PROF_grpid(id)    == gid            ) then
             write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I9)') &
                  '*** ID=',id,' : ',PROF_rapname(id),' T=',PROF_rapttot(id),' N=',PROF_rapnstr(id)
          endif
       enddo
       enddo

    else

       call ADM_MPItimestat( avgvar      (1:PROF_rapnmax), & ! [OUT]
                             maxvar      (1:PROF_rapnmax), & ! [OUT]
                             minvar      (1:PROF_rapnmax), & ! [OUT]
                             maxidx      (1:PROF_rapnmax), & ! [OUT]
                             minidx      (1:PROF_rapnmax), & ! [OUT]
                             PROF_rapttot(1:PROF_rapnmax)  ) ! [IN]

       fid = -1
       if ( .false. ) then ! report to STDOUT
          if ( ADM_IsMaster ) then
             write(*,*) '*** Computational Time Report'
             fid = 6 ! master node
          endif
       else
          fid = IO_FID_LOG
       endif

       do gid = 1, PROF_rapnmax
       do id  = 1, PROF_rapnmax
          if (       PROF_raplevel(id) <= PROF_rap_level &
               .AND. PROF_grpid(id)    == gid            &
               .AND. fid > 0                             ) then
             write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I9)') &
                  '*** ID=',id,' : ',PROF_rapname(id), &
                  ' T(avg)=',avgvar(id), &
                  ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                  ', T(min)=',minvar(id),'[',minidx(id),']', &
                  ' N=',PROF_rapnstr(id)
          endif
       enddo
       enddo

    endif

    return
  end subroutine PROF_rapreport

#ifdef _PAPI_
  !-----------------------------------------------------------------------------
  !> Start flop counter
  subroutine PROF_PAPI_rapstart
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PROF_PAPI_real_time, PROF_PAPI_proc_time, PROF_PAPI_flops, PROF_PAPI_mflops, PROF_PAPI_check )

    return
  end subroutine PROF_PAPI_rapstart

  !-----------------------------------------------------------------------------
  !> Stop flop counter
  subroutine PROF_PAPI_rapstop
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PROF_PAPI_real_time, PROF_PAPI_proc_time, PROF_PAPI_flops, PROF_PAPI_mflops, PROF_PAPI_check )

    return
  end subroutine PROF_PAPI_rapstop

  !-----------------------------------------------------------------------------
  !> Report flop
  subroutine PROF_PAPI_rapreport
    use mod_adm, only: &
       ADM_MPItimestat, &
       ADM_prc_all,     &
       ADM_IsMaster
    implicit none

    real(DP) :: avgvar(3)
    real(DP) :: maxvar(3)
    real(DP) :: minvar(3)
    integer  :: maxidx(3)
    integer  :: minidx(3)

    real(DP) :: PROF_PAPI_gflop
    real(DP) :: statistics(3)
    !---------------------------------------------------------------------------

    PROF_PAPI_gflop = real(PROF_PAPI_flops,kind=DP) / 1024.0_DP**3

    if ( IO_LOG_ALLNODE ) then ! report for each node

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '*** PAPI Report [Local PE information]'
       write(IO_FID_LOG,'(1x,A,F15.3)') '*** Real time          [sec] : ', PROF_PAPI_real_time
       write(IO_FID_LOG,'(1x,A,F15.3)') '*** CPU  time          [sec] : ', PROF_PAPI_proc_time
       write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP             [GFLOP] : ', PROF_PAPI_gflop
       write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOPS by PAPI   [GFLOPS] : ', PROF_PAPI_mflops/1024.0_DP
       write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP / CPU Time [GFLOPS] : ', PROF_PAPI_gflop/PROF_PAPI_proc_time

    else
       statistics(1) = real(PROF_PAPI_real_time,kind=DP)
       statistics(2) = real(PROF_PAPI_proc_time,kind=DP)
       statistics(3) = PROF_PAPI_gflop

       call ADM_MPItimestat( avgvar    (1:3), & ! [OUT]
                             maxvar    (1:3), & ! [OUT]
                             minvar    (1:3), & ! [OUT]
                             maxidx    (1:3), & ! [OUT]
                             minidx    (1:3), & ! [OUT]
                             statistics(1:3)  ) ! [IN]

       write(IO_FID_LOG,*)
       write(IO_FID_LOG,*) '*** PAPI Report'
       write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
       write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
       write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
       write(IO_FID_LOG,*)
       write(IO_FID_LOG,'(1x,A,F15.3,A,I6,A)') &
                  '*** TOTAL FLOP    [GFLOP] : ', avgvar(3)*ADM_prc_all, '(',ADM_prc_all,' PEs)'
       write(IO_FID_LOG,'(1x,A,F15.3)') &
                  '*** FLOPS        [GFLOPS] : ', avgvar(3)*ADM_prc_all/maxvar(2)
       write(IO_FID_LOG,'(1x,A,F15.3)') &
                  '*** FLOPS per PE [GFLOPS] : ', avgvar(3)/maxvar(2)
       write(IO_FID_LOG,*)

       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( ADM_IsMaster ) then ! master node
             write(*,*)
             write(*,*) '*** PAPI Report'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
             write(*,*)
             write(*,'(1x,A,F15.3,A,I6,A)') &
                  '*** TOTAL FLOP    [GFLOP] : ', avgvar(3)*ADM_prc_all, '(',ADM_prc_all,' PEs)'
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS        [GFLOPS] : ', avgvar(3)*ADM_prc_all/maxvar(2)
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS per PE [GFLOPS] : ', avgvar(3)/maxvar(2)
          endif
       endif
    endif

    return
  end subroutine PROF_PAPI_rapreport
#endif

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function get_rapid( rapname, level ) result(id)
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item
    integer,          intent(inout) :: level   !< level of item
    integer                         :: id

    character(len=H_SHORT*2) :: trapname
    character(len=H_SHORT)   :: trapname2
    !---------------------------------------------------------------------------

    trapname  = trim(rapname)
    trapname2 = trim(rapname)

    do id = 1, PROF_rapnmax
       if ( trapname == PROF_rapname(id) ) then
          level = PROF_raplevel(id)
          return
       endif
    enddo

    PROF_rapnmax     = PROF_rapnmax + 1
    id               = PROF_rapnmax
    PROF_rapname(id) = trapname

    PROF_rapnstr(id) = 0
    PROF_rapnend(id) = 0
    PROF_rapttot(id) = 0.0_DP

    PROF_grpid   (id) = get_grpid(trapname2)
    PROF_raplevel(id) = level

    return
  end function get_rapid

  !-----------------------------------------------------------------------------
  !> Get group ID
  function get_grpid( rapname ) result(gid)
    implicit none

    character(len=*), intent(in) :: rapname !< name of item
    integer                      :: gid

    character(len=H_SHORT) :: grpname
    integer                :: idx
    !---------------------------------------------------------------------------

    idx = index(rapname,' ')
    if ( idx > 1 ) then
       grpname = rapname(1:idx-1)
    else
       grpname = rapname
    endif

    do gid = 1, PROF_grpnmax
       if( grpname == PROF_grpname(gid) ) return
    enddo

    PROF_grpnmax      = PROF_grpnmax + 1
    gid               = PROF_grpnmax
    PROF_grpname(gid) = grpname

    return
  end function get_grpid

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_1D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, i1
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1)
          sp_min = var(1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1))
             sp_min = min(sp_min,  var(i1))
             sp_sum =     sp_sum + var(i1)
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_2D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, n2, i1, i2
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1,1)
          sp_min = var(1,1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i2 = 1, n2
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1,i2))
             sp_min = min(sp_min,  var(i1,i2))
             sp_sum =     sp_sum + var(i1,i2)
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_3D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, n2, n3, i1, i2, i3
    !---------------------------------------------------------------------------
    stats_acc = .false.
    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1,1,1)
          sp_min = var(1,1,1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1,i2,i3))
             sp_min = min(sp_min,  var(i1,i2,i3))
             sp_sum =     sp_sum + var(i1,i2,i3)
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_4D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, n2, n3, n4, i1, i2, i3, i4
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1,1,1,1)
          sp_min = var(1,1,1,1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1,i2,i3,i4))
             sp_min = min(sp_min,  var(i1,i2,i3,i4))
             sp_sum =     sp_sum + var(i1,i2,i3,i4)
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_5D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, n2, n3, n4, n5, i1, i2, i3, i4, i5
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1,1,1,1,1)
          sp_min = var(1,1,1,1,1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1,i2,i3,i4,i5))
             sp_min = min(sp_min,  var(i1,i2,i3,i4,i5))
             sp_sum =     sp_sum + var(i1,i2,i3,i4,i5)
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_6D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    real(SP) :: sp_sum, sp_max, sp_min
    integer  :: n1, n2, n3, n4, n5, n6, i1, i2, i3, i4, i5, i6
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          n6 = size(var,6)
          sp_sum = 0.0_SP

          !$acc kernels
          sp_max = var(1,1,1,1,1,1)
          sp_min = var(1,1,1,1,1,1)
          !$acc loop reduction(max:sp_max) reduction(min:sp_min) reduction(+:sp_sum)
          do i6 = 1, n6
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             sp_max = max(sp_max,  var(i1,i2,i3,i4,i5,i6))
             sp_min = min(sp_min,  var(i1,i2,i3,i4,i5,i6))
             sp_sum =     sp_sum + var(i1,i2,i3,i4,i5,i6)
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(sp_max,kind=DP)
          PROF_min = real(sp_min,kind=DP)
          PROF_sum = real(sp_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_SP_6D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_1D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, i1
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1)
          PROF_min = var(1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1))
             PROF_min = min(PROF_min,  var(i1))
             PROF_sum =     PROF_sum + var(i1)
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_2D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, n2, i1, i2
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1,1)
          PROF_min = var(1,1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i2 = 1, n2
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1,i2))
             PROF_min = min(PROF_min,  var(i1,i2))
             PROF_sum =     PROF_sum + var(i1,i2)
          enddo
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_3D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, n2, n3, i1, i2, i3
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1,1,1)
          PROF_min = var(1,1,1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1,i2,i3))
             PROF_min = min(PROF_min,  var(i1,i2,i3))
             PROF_sum =     PROF_sum + var(i1,i2,i3)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_4D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, n2, n3, n4, i1, i2, i3, i4
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1,1,1,1)
          PROF_min = var(1,1,1,1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1,i2,i3,i4))
             PROF_min = min(PROF_min,  var(i1,i2,i3,i4))
             PROF_sum =     PROF_sum + var(i1,i2,i3,i4)
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_5D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, n2, n3, n4, n5, i1, i2, i3, i4, i5
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1,1,1,1,1)
          PROF_min = var(1,1,1,1,1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1,i2,i3,i4,i5))
             PROF_min = min(PROF_min,  var(i1,i2,i3,i4,i5))
             PROF_sum =     PROF_sum + var(i1,i2,i3,i4,i5)
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_6D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: n1, n2, n3, n4, n5, n6, i1, i2, i3, i4, i5, i6
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          n6 = size(var,6)
          PROF_sum = 0.0_DP

          !$acc kernels
          PROF_max = var(1,1,1,1,1,1)
          PROF_min = var(1,1,1,1,1,1)
          !$acc loop reduction(max:PROF_max) reduction(min:PROF_min) reduction(+:PROF_sum)
          do i6 = 1, n6
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 2, n3
          do i2 = 1, n2
          do i1 = 1, n1
             PROF_max = max(PROF_max,  var(i1,i2,i3,i4,i5,i6))
             PROF_min = min(PROF_min,  var(i1,i2,i3,i4,i5,i6))
             PROF_sum =     PROF_sum + var(i1,i2,i3,i4,i5,i6)
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_DP_6D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_1D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, i1
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          int_sum = 0

          !$acc kernels
          int_max = var(1)
          int_min = var(1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i1 = 1, n1
             int_max = max(int_max,  var(i1))
             int_min = min(int_min,  var(i1))
             int_sum =     int_sum + var(i1)
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_2D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, n2, i1, i2
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          int_sum = 0

          !$acc kernels
          int_max = var(1,1)
          int_min = var(1,1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i2 = 1, n2
          do i1 = 1, n1
             int_max = max(int_max,  var(i1,i2))
             int_min = min(int_min,  var(i1,i2))
             int_sum =     int_sum + var(i1,i2)
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_3D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, n2, n3, i1, i2, i3
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          int_sum = 0

          !$acc kernels
          int_max = var(1,1,1)
          int_min = var(1,1,1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             int_max = max(int_max,  var(i1,i2,i3))
             int_min = min(int_min,  var(i1,i2,i3))
             int_sum =     int_sum + var(i1,i2,i3)
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_4D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, n2, n3, n4, i1, i2, i3, i4
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          int_sum = 0

          !$acc kernels
          int_max = var(1,1,1,1)
          int_min = var(1,1,1,1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             int_max = max(int_max,  var(i1,i2,i3,i4))
             int_min = min(int_min,  var(i1,i2,i3,i4))
             int_sum =     int_sum + var(i1,i2,i3,i4)
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_5D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, n2, n3, n4, n5, i1, i2, i3, i4, i5
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          int_sum = 0

          !$acc kernels
          int_max = var(1,1,1,1,1)
          int_min = var(1,1,1,1,1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             int_max = max(int_max,  var(i1,i2,i3,i4,i5))
             int_min = min(int_min,  var(i1,i2,i3,i4,i5))
             int_sum =     int_sum + var(i1,i2,i3,i4,i5)
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_INT_6D( &
       header,  &
       varname, &
       var,     &
       acc      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: var(:,:,:,:,:,:)
    character(len=*), intent(in), optional :: acc

    logical  :: stats_acc = .false.
    integer  :: int_sum, int_max, int_min
    integer  :: n1, n2, n3, n4, n5, n6, i1, i2, i3, i4, i5, i6
    !---------------------------------------------------------------------------

    if ( present(acc) ) then
       if ( acc == 'D2H' ) then
          stats_acc = .true.

          n1 = size(var,1)
          n2 = size(var,2)
          n3 = size(var,3)
          n4 = size(var,4)
          n5 = size(var,5)
          n6 = size(var,6)
          int_sum = 0

          !$acc kernels
          int_max = var(1,1,1,1,1,1)
          int_min = var(1,1,1,1,1,1)
          !$acc loop reduction(max:int_max) reduction(min:int_min) reduction(+:int_sum)
          do i6 = 1, n6
          do i5 = 1, n5
          do i4 = 1, n4
          do i3 = 1, n3
          do i2 = 1, n2
          do i1 = 1, n1
             int_max = max(int_max,  var(i1,i2,i3,i4,i5,i6))
             int_min = min(int_min,  var(i1,i2,i3,i4,i5,i6))
             int_sum =     int_sum + var(i1,i2,i3,i4,i5,i6)
          enddo
          enddo
          enddo
          enddo
          enddo
          enddo
          !$acc end kernels
          PROF_max = real(int_max,kind=DP)
          PROF_min = real(int_min,kind=DP)
          PROF_sum = real(int_sum,kind=DP)
       endif
    endif

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    if ( .NOT. stats_acc ) then
       PROF_max = real(maxval(var),kind=DP)
       PROF_min = real(minval(var),kind=DP)
       PROF_sum = real(sum   (var),kind=DP)
    endif

    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,', min=',PROF_min,', sum=',PROF_sum
    call flush(IO_FID_LOG)

    return
  end subroutine PROF_valcheck_INT_6D

end module mod_prof
