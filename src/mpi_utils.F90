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

end module mpi_utils
