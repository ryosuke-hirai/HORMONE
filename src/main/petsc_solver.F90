module petsc_solver_mod
#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
  use petsc
  use mpi_utils, only: myrank, stop_mpi
#endif
  implicit none
  private

  public :: init_petsc, finalise_petsc, solve_system_petsc

#ifdef USE_PETSC
  PetscErrorCode :: ierr
  Vec :: x_petsc, b_petsc
  KSP :: ksp
  PC :: pc
#else
  ! Dummy declarations for non-PETSc builds
  integer :: ierr
  integer :: x_petsc, b_petsc, ksp, pc
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

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SETUP_PETSC_A
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: PETSc branch for assembling the Laplacian matrix.
!          For each grid point, we compute the stencil using compute_coeffs,
!          then insert the values directly into the PETSc matrix (A_petsc).
!          Preconditioning is handled by PETSc, so we do not compute it here.

subroutine write_A_petsc(system)
  integer, intent(in) :: system
#ifdef USE_PETSC
  use grid,only:is=>gis,ie=>gie,js=>gjs,je=>gje,ks=>gks,ke=>gke
  use petsc_solver_mod,only:A_petsc_grv, A_petsc_rad
  integer :: in, jn, kn, dim, i, j, k, l, ncoeff, m
  real(8) :: coeffs(5)
  integer :: row, col
  integer :: ierr
  PetscInt    :: n
  PetscScalar :: val
  type(Mat) :: A_petsc
  integer :: lmax

  ! Select the appropriate PETSc matrix based on the system.
  if (system == igrv) then
    A_petsc = A_petsc_grv
  else if (system == irad) then
    A_petsc = A_petsc_rad
  else
    print *, 'Error in setup_petsc_A: unsupported system'
    stop
  end if

  ! Compute grid dimensions and total number of points.
  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1
  lmax = in * jn * kn

  if (system == igrv) then
    lmax_grv = lmax
  else if (system == irad) then
    lmax_rad = lmax
  end if

  ! Determine problem dimension based on input indices.
  if(ie > is .and. je > js .and. ke > ks) then
    dim = 3
    ncoeff = 5
  else if(ie > is .and. (je > js .or. ke > ks)) then
    dim = 2
    ncoeff = 3
  else if(ie > is .and. (je == js .and. ke == ks)) then
    dim = 1
    ncoeff = 2
  else
    print *, 'Error in setup_petsc_A: unsupported dimension'
    stop
  end if

  ! Set up PETSc matrix
  call MatCreate(PETSC_COMM_WORLD, A_petsc, ierr)
  n = lmax
  call MatSetSizes(A_petsc, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  call MatSetFromOptions(A_petsc, ierr)
  call MatSetUp(A_petsc, ierr)

  ! Loop over all grid points.
  do l = 1, lmax
    call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
    call compute_coeffs(system, dim, i, j, k, ie, in, jn, kn, coeffs)

    ! As the dimension increases, additional diagonals are added, but the
    ! offsets of existing diagonals are not changed.
    ! In 1D, the offsets are at 0, 1.
    ! In 2D, the offsets are at 0, 1, in
    ! In 3D, the offsets are at 0, 1, in, in*jn, in*jn*(kn-1)
    ! These are the same for gravity and radiation
    do m = 1, ncoeff
      select case(m)
      case(1)
        row = l - 1
        col = l - 1
      case(2)
        row = l - 1
        col = l - 1 + 1
      case(3)
        row = l - 1
        col = l - 1 + in
      case(4)
        row = l - 1
        col = l - 1 + in * jn
      case(5)
        row = l - 1
        col = l - 1 + in * jn * (kn - 1)
      end select

      if(row >= 0 .and. row < lmax .and. col >= 0 .and. col < lmax) then
        val = coeffs(m)
        ! Symmetric matrix
                       call MatSetValue(A_petsc, row, col, val, INSERT_VALUES, ierr)
        if(row /= col) call MatSetValue(A_petsc, col, row, val, INSERT_VALUES, ierr)
      end if
    enddo
  end do

  ! Finalize PETSc matrix assembly.
  call MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY, ierr)

#endif
end subroutine write_A_petsc

  subroutine solve_system_petsc(A_petsc, cgsrc, x)
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
#ifdef USE_PETSC
    Mat, intent(in) :: A_petsc
#else
    integer, intent(in) :: A_petsc
#endif
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
