!-------------------------------------------------------------------------------
!> module STANDARD I/O
!!
!! @par Description
!!         Standard, common I/O
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_stdio
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
  public :: IO_get_available_fid
  public :: IO_close_allfile

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,               public, parameter :: H_SHORT         =   16    ! character length for short var.
  integer,               public, parameter :: H_MID           =   64    ! character length for middle var.
  integer,               public, parameter :: H_LONG          =  256    ! character length for long var.
  integer,               public, parameter :: H_FILE          = 1024    ! character length for long var.

  integer,               public            :: IO_FID_CONF     = -1      ! Config file ID
  integer,               public            :: IO_FID_LOG      = -1      ! Log file ID
  integer,               public            :: IO_FID_NML      = -1      ! NAMELIST file ID
  character(len=H_LONG), public            :: IO_LOG_BASENAME = 'msg'   ! basename of logfile

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: IO_MINFID =  10 !< minimum available fid
  integer, private, parameter :: IO_MAXFID = 900 !< maximum available fid

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> search & get available file ID
  !> @return fid
  function IO_get_available_fid() result(fid)
    implicit none

    integer :: fid      ! file ID
    logical :: i_opened ! file ID is already used or not?
    !---------------------------------------------------------------------------

    do fid = IO_MINFID, IO_MAXFID
       inquire(fid,OPENED=i_opened)
       if( .NOT. i_opened ) return
    enddo

    if ( fid >= IO_MAXFID ) then ! over limit
       write(IO_FID_LOG,*) 'xxx Used I/O unit number reached to the limit! STOP'
       write(*,*)          'xxx Used I/O unit number reached to the limit! STOP'
       stop 1
    endif

  end function IO_get_available_fid

  !-----------------------------------------------------------------------------
  !> search & close all file ID
  subroutine IO_close_allfile
    implicit none

    integer :: fid      ! file ID
    logical :: i_opened ! file ID is already used or not?
    !---------------------------------------------------------------------------

    do fid = IO_MINFID, IO_MAXFID
       inquire(fid,OPENED=i_opened)
       if( i_opened ) close(fid)
    enddo

    return
  end subroutine IO_close_allfile

end module mod_stdio
