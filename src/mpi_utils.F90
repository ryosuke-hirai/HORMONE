module mpi_utils
#ifdef MPI
  use mpi
#endif
  implicit none

  integer :: myrank, nprocs

  contains

  subroutine init_mpi()
#ifdef MPI
    integer :: ierr
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
#endif
  end subroutine init_mpi

  subroutine finalize_mpi()
#ifdef MPI
    integer :: ierr
    call MPI_FINALIZE(ierr)
#endif
  end subroutine finalize_mpi

  subroutine allreduce_mpi(op_str, value)
    character(len=*), intent(in) :: op_str
    real(8), intent(inout) :: value

#ifdef MPI
    real(8) :: value_reduced
    integer :: op, ierr

    select case (op_str)
    case ('sum')
      op = MPI_SUM
    case ('max')
      op = MPI_MAX
    case ('min')
      op = MPI_MIN
    case default
      error stop "Invalid operation specified"
    end select

    call MPI_ALLREDUCE(value, value_reduced, 1, MPI_DOUBLE_PRECISION, op, MPI_COMM_WORLD, ierr)
    value = value_reduced
#endif

  end subroutine allreduce_mpi

end module mpi_utils
