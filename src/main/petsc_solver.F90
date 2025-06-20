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
subroutine setup_petsc(is, ie, js, je, ks, ke, is_global, ie_global, js_global, je_global, ks_global, ke_global, pm)
  use settings, only: cgerr
  use utils, only: get_dim
  use matrix_utils, only: l_from_ijk, ijk_from_l
  use matrix_vars, only: petsc_set
  integer, intent(in) :: is, ie, js, je, ks, ke
  integer, intent(in) :: is_global, ie_global, js_global, je_global, ks_global, ke_global
  type(petsc_set), intent(out) :: pm
  integer :: i, j, k, l, ll
  PC  :: pc
  PetscErrorCode :: ierr
  PetscInt :: ls, le

  ! Compute grid dimensions and total number of points.
  pm%is = is; pm%ie = ie
  pm%js = js; pm%je = je
  pm%ks = ks; pm%ke = ke

  ! in, jn, kn are the number of grid points in each direction locally
  pm%in = ie - is + 1
  pm%jn = je - js + 1
  pm%kn = ke - ks + 1

  pm%lmax = pm%in * pm%jn * pm%kn

  pm%is_global = is_global; pm%ie_global = ie_global
  pm%js_global = js_global; pm%je_global = je_global
  pm%ks_global = ks_global; pm%ke_global = ke_global

  pm%in_global = ie_global - is_global + 1
  pm%jn_global = je_global - js_global + 1
  pm%kn_global = ke_global - ks_global + 1

  pm%lmax_global = pm%in_global * pm%jn_global * pm%kn_global

  call get_dim(pm%is_global, pm%ie_global, pm%js_global, pm%je_global, pm%ks_global, pm%ke_global, pm%dim)

  ! Set up matrix A
  call MatCreate(PETSC_COMM_WORLD, pm%A, ierr)
  ! Divide the matrix between processes
  call MatSetSizes(pm%A, pm%lmax, pm%lmax, pm%lmax_global, pm%lmax_global, ierr)
  call MatSetFromOptions(pm%A, ierr)
  call MatSetUp(pm%A, ierr)

  ! Create vector for right-hand side.
  call VecCreate(PETSC_COMM_WORLD, pm%b, ierr)
  call VecSetSizes(pm%b, pm%lmax, pm%lmax_global, ierr)
  call VecSetFromOptions(pm%b, ierr)

  ! Create vector for solution.
  call VecCreate(PETSC_COMM_WORLD, pm%x, ierr)
  call VecSetSizes(pm%x, pm%lmax, pm%lmax_global, ierr)
  call VecSetFromOptions(pm%x, ierr)

  ! Get the ownership ranges for matrix
  call MatGetOwnershipRange(pm%A, ls, le, ierr)

  ! print*, "Rank ", myrank, ":  Matrix ownership range: ", ls, " to ", le, " (global size: ", pm%lmax_global, ")"

  ! Create the KSP solver context.
  call KSPCreate(PETSC_COMM_WORLD, pm%ksp, ierr)

  ! Re-use the solution from the previous iteration as the initial guess.
  call KSPSetInitialGuessNonzero(pm%ksp, PETSC_TRUE, ierr)

  ! Equivalent to options: -ksp_type cg -pc_type bjacobi -ksp_rtol [cgerr]
  call KSPSetType(pm%ksp, KSPCG, ierr)
  call KSPGetPC(pm%ksp, pc, ierr)
  call PCSetType(pc, PCBJACOBI, ierr)
  call KSPSetTolerances(pm%ksp, cgerr, -1.d0, -1.d0, PETSC_DECIDE, ierr)

  ! Allow command line options/overrides
  call KSPSetFromOptions(pm%ksp, ierr)

  ! Define the local rows for this task.
  allocate(pm%my_rows(pm%lmax))

  ll = 0
  do k = ks, ke
    do j = js, je
      do i = is, ie
        ll = ll + 1
        l = l_from_ijk(i, j, k, pm%is_global, pm%js_global, pm%ks_global, pm%in_global, pm%jn_global)
        pm%my_rows(ll) = l
      end do
    end do
  end do

  if (ll /= pm%lmax) then
    print *, "Error: Number of local rows does not match expected count."
    call stop_mpi(1)
  end if

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
subroutine write_A_petsc(system, pm)
  use settings,only:crdnt,solve_j,solve_k
  use utils, only: get_dim
  use matrix_coeffs, only:compute_coeffs, get_matrix_offsets
  use matrix_utils, only:l_from_ijk, get_raddim, ijk_from_l
  use matrix_vars, only: petsc_set, irad, igrv
  integer, intent(in) :: system
  type(petsc_set), intent(inout) :: pm
  integer :: dim, i, j, k, l, ll, ncoeff, m
  real(8) :: coeffs(5)
  integer :: offsets(5)
  integer :: row, col
  PC             :: pc
  PetscScalar    :: val
  PetscErrorCode :: ierr

  ncoeff = pm%dim + 1
  if(crdnt==1.and.solve_j)ncoeff=ncoeff+1
  if(crdnt==2.and.solve_k)ncoeff=ncoeff+1

  call get_matrix_offsets(pm%dim, offsets)

  ! The radiation coefficients depend on the direction of the problem in 1D and 2D
  if (system == irad) then
    dim = get_raddim(pm%in, pm%jn, pm%kn)
  elseif (system == igrv) then
    dim = pm%dim
  endif

  ! Loop over all grid points.
  ! l is the row number in the matrix
  do ll = 1, pm%lmax
    l = pm%my_rows(ll)
    call ijk_from_l(l, pm%is_global, pm%js_global, pm%ks_global, pm%in_global, pm%jn_global, i, j, k)
    call compute_coeffs(system, dim, i, j, k, coeffs)

    do m = 1, ncoeff
      row = (l - 1)
      col = row + offsets(m)

      if(row >= 0 .and. row < pm%lmax_global .and. col >= 0 .and. col < pm%lmax_global) then
        val = coeffs(m)
        ! Symmetric matrix
                       call MatSetValue(pm%A, row, col, val, INSERT_VALUES, ierr)
        if(row /= col) call MatSetValue(pm%A, col, row, val, INSERT_VALUES, ierr) ! Not guaranteed to be in the current task's ownership range
      end if
    enddo
  end do

  ! Finalize PETSc matrix assembly.
  call MatAssemblyBegin(pm%A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(pm%A, MAT_FINAL_ASSEMBLY, ierr)

  ! Set operators for the KSP solver.
  call KSPSetOperators(pm%ksp, pm%A, pm%A, ierr)

  ! Set preconditioner.
  call KSPGetPC(pm%ksp, pc, ierr)
  call PCSetup(pc, ierr)
  call KSPSetup(pm%ksp, ierr)

end subroutine write_A_petsc
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SOLVE_SYSTEM_PETSC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solves the linear system A*x = b using PETSc.

#ifdef USE_PETSC
subroutine solve_system_petsc(pm, b, x)
  use matrix_vars, only: petsc_set
  type(petsc_set), intent(in) :: pm
  real(8), intent(in)  :: b(pm%lmax)
  real(8), intent(out) :: x(pm%lmax)
  integer :: ierr, l, ll
  real(8), pointer :: x_array(:)
  PetscInt, allocatable :: idx_array(:)
  IS :: is_from
  Vec :: x_gather
  VecScatter :: scatter_ctx

  ! Set the right-hand side vector b.
  do ll = 1, pm%lmax
    l = pm%my_rows(ll)
    call VecSetValue(pm%b, l-1, b(ll), INSERT_VALUES, ierr)
  end do
  call VecAssemblyBegin(pm%b, ierr)
  call VecAssemblyEnd(pm%b, ierr)

  ! Solve the linear system.
  call KSPSolve(pm%ksp, pm%b, pm%x, ierr)

  ! Create index set with only our owned indices
  allocate(idx_array(pm%lmax))
  idx_array(:) = pm%my_rows(:) - 1  ! Convert to zero-based indexing
  call ISCreateGeneral(PETSC_COMM_SELF, pm%lmax, idx_array, &
                       PETSC_COPY_VALUES, is_from, ierr)

  ! Create small vector to receive only our owned values
  call VecCreate(PETSC_COMM_SELF, x_gather, ierr)
  call VecSetSizes(x_gather, PETSC_DECIDE, pm%lmax, ierr)
  call VecSetFromOptions(x_gather, ierr)

  ! Create scatter context to extract only our owned values
  call VecScatterCreate(pm%x, is_from, x_gather, PETSC_NULL_IS, scatter_ctx, ierr)

  ! Perform scatter to extract only our owned values
  call VecScatterBegin(scatter_ctx, pm%x, x_gather, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(scatter_ctx, pm%x, x_gather, INSERT_VALUES, SCATTER_FORWARD, ierr)

  ! Copy solution back to Fortran array x.
  ! Leave values in pm%x to re-use in the next iteration.
  call VecGetArrayF90(x_gather, x_array, ierr)
  x(1:pm%lmax) = x_array(1:pm%lmax)
  call VecRestoreArrayF90(x_gather, x_array, ierr)

end subroutine solve_system_petsc
#endif

end module petsc_solver_mod
