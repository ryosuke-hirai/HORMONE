module petsc_solver_mod
#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
  use petsc
  use mpi_utils, only: myrank, stop_mpi
#endif
  implicit none
  private

  public :: init_petsc, finalise_petsc, solve_system_petsc
  public :: A_petsc

#ifdef USE_PETSC
  PetscErrorCode :: ierr
  Mat :: A_petsc
  Vec :: x_petsc, b_petsc
  KSP :: ksp
  PC :: pc
#else
  ! Dummy declarations for non-PETSc builds
  integer :: ierr
  integer :: A_petsc, x_petsc, b_petsc, ksp, pc
#endif

  contains

  subroutine init_petsc
#ifdef USE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) stop 'PETSc initialization failed'
#endif
  end subroutine init_petsc

  subroutine finalise_petsc
#ifdef USE_PETSC
    call PetscFinalize(ierr)
#endif
  end subroutine finalise_petsc

  subroutine solve_system_petsc(cgsrc, x)
    !-------------------------------------------------------------------
    ! Solves the linear system A*x = b using PETSc.
    !
    ! Input:
    !   cgsrc : real(8) array for the right-hand side vector b.
    !
    ! In/Out:
    !   x     : real(8) array where the solution is placed.
    !
    ! PETSc objects are created locally and then destroyed.
    !-------------------------------------------------------------------
    real(8), intent(in) :: cgsrc(:)
    real(8), intent(inout) :: x(:)
#ifdef USE_PETSC
    integer :: ierr, i, lmax
    PetscInt :: N
    real(8), pointer :: x_array(:)

    ! System size
    lmax = size(cgsrc)
    N = lmax

    ! Create vector for right-hand side.
    call VecCreate(PETSC_COMM_WORLD, b_petsc, ierr)
    call VecSetSizes(b_petsc, PETSC_DECIDE, N, ierr)
    call VecSetFromOptions(b_petsc, ierr)

    ! Create vector for solution.
    call VecCreate(PETSC_COMM_WORLD, x_petsc, ierr)
    call VecSetSizes(x_petsc, PETSC_DECIDE, N, ierr)
    call VecSetFromOptions(x_petsc, ierr)

    ! Insert values from cgsrc into b_petsc.
    do i = 0, N-1
      call VecSetValue(b_petsc, i, cgsrc(i+1), INSERT_VALUES, ierr)
    end do
    call VecAssemblyBegin(b_petsc, ierr)
    call VecAssemblyEnd(b_petsc, ierr)

    ! Create KSP solver context.
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, A_petsc, A_petsc, ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
    call KSPSetFromOptions(ksp, ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetup(pc,ierr)
    call KSPSetup(ksp, ierr)

    ! Solve the linear system.
    call KSPSolve(ksp, b_petsc, x_petsc, ierr)

    ! Copy solution back to x.
    call VecGetArrayF90(x_petsc, x_array, ierr)

    do i = 1, lmax
      x(i) = x_array(i)
    end do
    call VecRestoreArrayF90(x_petsc, x_array, ierr)

    ! Clean up PETSc objects.
    call KSPDestroy(ksp, ierr)
    call VecDestroy(x_petsc, ierr)
    call VecDestroy(b_petsc, ierr)

#endif
  end subroutine solve_system_petsc

end module petsc_solver_mod
