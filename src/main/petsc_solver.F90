module petsc_solver_mod
#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
  use petsc
  use mpi_utils, only: myrank, stop_mpi
  implicit none
  ! private

  public :: init_petsc, finalise_petsc

  PetscErrorCode :: ierr
  Mat :: A_petsc
  Vec :: x_petsc, b_petsc
  KSP :: ksp
  PC :: pc

  contains

  subroutine init_petsc
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) stop 'PETSc initialization failed'
  end subroutine init_petsc

  subroutine finalise_petsc
    call PetscFinalize(ierr)
  end subroutine finalise_petsc

#endif
end module petsc_solver_mod
