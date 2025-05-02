module petsc_solver_mod
#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
  use petsc
  use mpi_utils, only: myrank, stop_mpi
#endif
  implicit none
  private

  public :: init_petsc, finalize_petsc
#ifdef USE_PETSC
  public :: solve_system_petsc, write_A_petsc, setup_petsc
#endif

  contains

  subroutine init_petsc
#ifdef USE_PETSC
    PetscErrorCode :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) stop 'PETSc initialization failed'
#endif
  end subroutine init_petsc

  subroutine finalize_petsc
#ifdef USE_PETSC
    PetscErrorCode :: ierr
    call PetscFinalize(ierr)
#endif
  end subroutine finalize_petsc


#ifdef USE_PETSC
subroutine setup_petsc(A_petsc, b_petsc, x_petsc, ksp, lmax_petsc)
  use grid, only: is, ie, js, je, ks, ke
  integer :: in, jn, kn
  Mat :: A_petsc
  Vec :: b_petsc, x_petsc
  KSP :: ksp
  PetscInt :: lmax_petsc
  PetscErrorCode :: ierr

  ! Compute grid dimensions and total number of points.
  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1
  lmax_petsc = in * jn * kn

  ! Set up matrix A
  call MatCreate(PETSC_COMM_WORLD, A_petsc, ierr)
  call MatSetSizes(A_petsc, PETSC_DECIDE, PETSC_DECIDE, lmax_petsc, lmax_petsc, ierr)
  call MatSetFromOptions(A_petsc, ierr)
  call MatSetUp(A_petsc, ierr)
  ! Create vector for right-hand side.
  call VecCreate(PETSC_COMM_WORLD, b_petsc, ierr)
  call VecSetSizes(b_petsc, PETSC_DECIDE, lmax_petsc, ierr)
  call VecSetFromOptions(b_petsc, ierr)
  ! Create vector for solution.
  call VecCreate(PETSC_COMM_WORLD, x_petsc, ierr)
  call VecSetSizes(x_petsc, PETSC_DECIDE, lmax_petsc, ierr)
  call VecSetFromOptions(x_petsc, ierr)
  ! Create the KSP solver context.
  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
  call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
  call KSPSetFromOptions(ksp, ierr)

end subroutine setup_petsc
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE WRITE_A_PETSC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: PETSc branch for assembling the Laplacian matrix.
!          For each grid point, we compute the stencil using compute_coeffs,
!          then insert the values directly into the PETSc matrix (A_petsc).
!          Preconditioning is handled by PETSc, so we do not compute it here.

#ifdef USE_PETSC
subroutine write_A_petsc(A_petsc, ksp, lmax_petsc, system)
  use grid,only:is=>gis,ie=>gie,js=>gjs,je=>gje,ks=>gks,ke=>gke
  use matrix_coeffs_mod,only:compute_coeffs
  use miccg_mod,only:ijk_from_l ! TODO: move this elsewhere
  Mat, intent(in) :: A_petsc
  KSP, intent(in) :: ksp
  PetscInt, intent(in) :: lmax_petsc
  integer, intent(in) :: system
  integer :: in, jn, kn, dim, i, j, k, l, ncoeff, m
  real(8) :: coeffs(5)
  integer :: row, col
  PC             :: pc
  PetscScalar    :: val
  PetscErrorCode :: ierr

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

  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1

  ! Loop over all grid points.
  do l = 1, lmax_petsc
    call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
    call compute_coeffs(system, dim, i, j, k, coeffs)

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

      if(row >= 0 .and. row < lmax_petsc .and. col >= 0 .and. col < lmax_petsc) then
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

  ! Set operators for the KSP solver.
  call KSPSetOperators(ksp, A_petsc, A_petsc, ierr)

  ! Set preconditioner.
  call KSPGetPC(ksp,pc,ierr)
  call PCSetup(pc,ierr)
  call KSPSetup(ksp, ierr)

end subroutine write_A_petsc
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SOLVE_SYSTEM_PETSC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solves the linear system A*x = b using PETSc.

#ifdef USE_PETSC
subroutine solve_system_petsc(x_petsc, b_petsc, ksp, lmax_petsc, x, b)
  Vec, intent(in) :: x_petsc
  Vec, intent(in) :: b_petsc
  KSP, intent(in) :: ksp
  PetscInt, intent(in) :: lmax_petsc
  real(8), intent(inout) :: x(:)
  real(8), intent(in)    :: b(:)
  integer :: ierr, i
  real(8), pointer :: x_array(:)

  ! Insert values from b into b_petsc.
  do i = 0, lmax_petsc-1
    call VecSetValue(b_petsc, i, b(i+1), INSERT_VALUES, ierr)
  end do
  call VecAssemblyBegin(b_petsc, ierr)
  call VecAssemblyEnd(b_petsc, ierr)

  ! Solve the linear system.
  call KSPSolve(ksp, b_petsc, x_petsc, ierr)

  ! Copy solution back to x.
  call VecGetArrayF90(x_petsc, x_array, ierr)
  do i = 1, lmax_petsc
    x(i) = x_array(i)
  end do
  call VecRestoreArrayF90(x_petsc, x_array, ierr)

end subroutine solve_system_petsc
#endif

end module petsc_solver_mod
