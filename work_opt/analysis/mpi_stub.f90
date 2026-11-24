! Minimal 1-rank MPI module stub for local macOS gfortran verification (no MPI needed)
module mpi
  implicit none
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: MPI_SUCCESS = 0
  integer, parameter :: MPI_STATUS_SIZE = 5
  integer, parameter :: MPI_DOUBLE_PRECISION = 1
  integer, parameter :: MPI_REAL = 2
  integer, parameter :: MPI_INTEGER = 3
  integer, parameter :: MPI_CHARACTER = 4
  integer, parameter :: MPI_LOGICAL = 5
  integer, parameter :: MPI_SUM = 1
  integer, parameter :: MPI_MAX = 2
  integer, parameter :: MPI_MIN = 3
contains
  subroutine MPI_Init(ierr)
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine MPI_Init
  subroutine MPI_Initialized(flag, ierr)
    logical, intent(out) :: flag
    integer, intent(out) :: ierr
    flag = .true.; ierr = 0
  end subroutine MPI_Initialized
  subroutine MPI_Comm_size(comm, sz, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: sz, ierr
    sz = 1; ierr = 0
  end subroutine MPI_Comm_size
  subroutine MPI_Comm_rank(comm, rk, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: rk, ierr
    rk = 0; ierr = 0
  end subroutine MPI_Comm_rank
  subroutine MPI_Barrier(comm, ierr)
    integer, intent(in) :: comm
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine MPI_Barrier
  subroutine MPI_Finalize(ierr)
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine MPI_Finalize
  subroutine MPI_Abort(comm, errcode, ierr)
    integer, intent(in) :: comm, errcode
    integer, intent(out) :: ierr
    ierr = 0
    stop 1
  end subroutine MPI_Abort
  double precision function MPI_Wtime()
    integer(8) :: c, r
    call system_clock(c, r)
    MPI_Wtime = dble(c) / dble(r)
  end function MPI_Wtime
end module mpi

subroutine mpi_allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm, ierr)
  implicit none
  real(8), intent(in) :: sendbuf(*)
  real(8), intent(out) :: recvbuf(*)
  integer, intent(in) :: sendcnt, sendtype, recvcnt, recvtype, comm
  integer, intent(out) :: ierr
  recvbuf(1:sendcnt) = sendbuf(1:sendcnt)
  ierr = 0
end subroutine mpi_allgather
