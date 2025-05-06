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
subroutine setup_petsc(petsc_matrix)
  use grid, only: is, ie, js, je, ks, ke
  use matrix_vars, only: petsc_set
  type(petsc_set), intent(out) :: petsc_matrix
  PC  :: pc
  PetscErrorCode :: ierr

  ! Compute grid dimensions and total number of points.
  petsc_matrix%is = is; petsc_matrix%ie = ie
  petsc_matrix%js = js; petsc_matrix%je = je
  petsc_matrix%ks = ks; petsc_matrix%ke = ke

  petsc_matrix%in = ie - is + 1
  petsc_matrix%jn = je - js + 1
  petsc_matrix%kn = ke - ks + 1

  petsc_matrix%lmax = petsc_matrix%in * petsc_matrix%jn * petsc_matrix%kn

  ! Set up matrix A
  call MatCreate(PETSC_COMM_WORLD, petsc_matrix%A, ierr)
  call MatSetSizes(petsc_matrix%A, PETSC_DECIDE, PETSC_DECIDE, petsc_matrix%lmax, petsc_matrix%lmax, ierr)
  call MatSetFromOptions(petsc_matrix%A, ierr)
  call MatSetUp(petsc_matrix%A, ierr)
  ! Create vector for right-hand side.
  call VecCreate(PETSC_COMM_WORLD, petsc_matrix%b, ierr)
  call VecSetSizes(petsc_matrix%b, PETSC_DECIDE, petsc_matrix%lmax, ierr)
  call VecSetFromOptions(petsc_matrix%b, ierr)
  ! Create vector for solution.
  call VecCreate(PETSC_COMM_WORLD, petsc_matrix%x, ierr)
  call VecSetSizes(petsc_matrix%x, PETSC_DECIDE, petsc_matrix%lmax, ierr)
  call VecSetFromOptions(petsc_matrix%x, ierr)
  ! Create the KSP solver context.
  call KSPCreate(PETSC_COMM_WORLD, petsc_matrix%ksp, ierr)
  call KSPSetInitialGuessNonzero(petsc_matrix%ksp, PETSC_TRUE, ierr)
  ! Equivalent to options: -ksp_type cg -pc_type bjacobi -ksp_rtol 1.e-16
  call KSPSetType(petsc_matrix%ksp, KSPCG, ierr)
  call KSPGetPC(petsc_matrix%ksp, pc, ierr)
  call PCSetType(pc, PCBJACOBI, ierr)
  call KSPSetTolerances(petsc_matrix%ksp, 1.d-16, -1.d0, -1.d0, PETSC_DECIDE, ierr)
  ! Allow command line options/overrides
  call KSPSetFromOptions(petsc_matrix%ksp, ierr)

end subroutine setup_petsc
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE WRITE_A_PETSC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: PETSc branch for assembling the Laplacian matrix.
!          For each grid point, we compute the stencil using compute_coeffs,
!          then insert the values directly into the PETSc matrix.
!          Preconditioning is handled by PETSc, so we do not compute it here.

#ifdef USE_PETSC
subroutine write_A_petsc(petsc_matrix, system)
  use utils, only: get_dim
  use grid, only:is=>gis,ie=>gie,js=>gjs,je=>gje,ks=>gks,ke=>gke
  use matrix_coeffs_mod, only:compute_coeffs, get_matrix_offsets
  use matrix_utils, only:ijk_from_l
  use matrix_vars, only: petsc_set
  type(petsc_set), intent(inout) :: petsc_matrix
  integer, intent(in) :: system
  integer :: in, jn, kn, dim, i, j, k, l, ncoeff, m
  real(8) :: coeffs(5)
  integer :: offsets(5)
  integer :: row, col
  PC             :: pc
  PetscScalar    :: val
  PetscErrorCode :: ierr

  ! Determine problem dimension based on input indices.
  call get_dim(is, ie, js, je, ks, ke, dim)

  select case (dim)
  case(1)
    ncoeff = 2
  case(2)
    ncoeff = 3
  case(3)
    ncoeff = 5
  end select

  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1

  call get_matrix_offsets(dim, offsets)

  ! Loop over all grid points.
  do l = 1, petsc_matrix%lmax
    call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
    call compute_coeffs(system, dim, i, j, k, coeffs)

    do m = 1, ncoeff
      row = (l - 1)
      col = row + offsets(m)

      if(row >= 0 .and. row < petsc_matrix%lmax .and. col >= 0 .and. col < petsc_matrix%lmax) then
        val = coeffs(m)
        ! Symmetric matrix
                       call MatSetValue(petsc_matrix%A, row, col, val, INSERT_VALUES, ierr)
        if(row /= col) call MatSetValue(petsc_matrix%A, col, row, val, INSERT_VALUES, ierr)
      end if
    enddo
  end do

  ! Finalize PETSc matrix assembly.
  call MatAssemblyBegin(petsc_matrix%A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(petsc_matrix%A, MAT_FINAL_ASSEMBLY, ierr)

  ! Set operators for the KSP solver.
  call KSPSetOperators(petsc_matrix%ksp, petsc_matrix%A, petsc_matrix%A, ierr)

  ! Set preconditioner.
  call KSPGetPC(petsc_matrix%ksp, pc, ierr)
  call PCSetup(pc, ierr)
  call KSPSetup(petsc_matrix%ksp, ierr)

end subroutine write_A_petsc
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SOLVE_SYSTEM_PETSC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solves the linear system A*x = b using PETSc.

#ifdef USE_PETSC
subroutine solve_system_petsc(petsc_matrix, x, b)
  use matrix_vars, only: petsc_set
  type(petsc_set), intent(in) :: petsc_matrix
  real(8), intent(inout) :: x(:)
  real(8), intent(in)    :: b(:)
  integer :: ierr, i
  real(8), pointer :: x_array(:)

  do i = 0, petsc_matrix%lmax-1
    call VecSetValue(petsc_matrix%b, i, b(i+1), INSERT_VALUES, ierr)
  end do
  call VecAssemblyBegin(petsc_matrix%b, ierr)
  call VecAssemblyEnd(petsc_matrix%b, ierr)

  ! Solve the linear system.
  call KSPSolve(petsc_matrix%ksp, petsc_matrix%b, petsc_matrix%x, ierr)

  ! Copy solution back to x.
  call VecGetArrayF90(petsc_matrix%x, x_array, ierr)
  do i = 1, petsc_matrix%lmax
    x(i) = x_array(i)
  end do
  call VecRestoreArrayF90(petsc_matrix%x, x_array, ierr)

end subroutine solve_system_petsc
#endif

end module petsc_solver_mod
