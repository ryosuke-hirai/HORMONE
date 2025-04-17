module matrix_solver
  use miccg_mod, only: cg_set, ijk_from_l, get_preconditioner
  use matrix_coeffs_mod, only: compute_coeffs_gravity
#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
  use petsc, only: MatSetValue, MatAssemblyBegin, MatAssemblyEnd, &
  INSERT_VALUES, MAT_FINAL_ASSEMBLY, PETSC_DECIDE, PETSC_COMM_WORLD
#endif
  implicit none
  public :: setup_grvA

  integer, public :: lmax

  contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SETUP_GRVA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Top-level assembly routine for the matrix A
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.
!          - The MICCG branch uses the compute_coeffs to fill the cg object.
!          - The PETSc branch directly fills the PETSc matrix A_petsc.
!            The preconditioner is handled by PETSc.
!
! TODO: This needs to be generalised to handle the radiation matrix as well

subroutine setup_grvA

#ifdef USE_PETSC
  call setup_petsc_A
#else
  call setup_cg
#endif

end subroutine setup_grvA

subroutine solve_system(cgsrc, x)
  use petsc_solver_mod
  real(8), allocatable, intent(in) :: cgsrc(:)
  real(8), allocatable, intent(inout) :: x(:)

#ifdef USE_PETSC
  ! PETSc solver
  call solve_system_petsc(cgsrc, x)
#else
  ! Original MICCG solver
  call miccg(cg, cgsrc, x)
#endif

end subroutine solve_system

#ifdef USE_PETSC
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SETUP_PETSC_A
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: PETSc branch for assembling the Laplacian matrix.
!          For each grid point, we compute the stencil using compute_coeffs,
!          then insert the values directly into the PETSc matrix (A_petsc).
!          Preconditioning is handled by PETSc, so we do not compute it here.

subroutine setup_petsc_A
  use grid,only:is=>gis,ie=>gie,js=>gjs,je=>gje,ks=>gks,ke=>gke
  use petsc_solver_mod,only:A_petsc
  integer :: in, jn, kn, dim, i, j, k, l, ncoeff, m
  real(8) :: coeffs(5)
  integer :: row, col
  integer :: ierr
  PetscInt    :: n
  PetscScalar :: val

  ! Compute grid dimensions and total number of points.
  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1
  lmax = in * jn * kn

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
    call compute_coeffs_gravity(dim, i, j, k, ie, kn, coeffs)

    ! As the dimension increases, additional diagonals are added, but the
    ! offsets of existing diagonals are not changed.
    ! In 1D, the offsets are at 0, 1.
    ! In 2D, the offsets are at 0, 1, in
    ! In 3D, the offsets are at 0, 1, in, in*jn, in*jn*(kn-1)
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

end subroutine setup_petsc_A
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SETUP_CG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: This is the assembly routine for the MICCG branch.
!          The coefficients are computed using compute_coeffs and
!          the preconditioner is computed as well.
!          The resulting data is stored in the cg object.

subroutine setup_cg
  use grid,only:is=>gis,ie=>gie,js=>gjs,je=>gje,ks=>gks,ke=>gke
  use miccg_mod, only: cg=>cg_grv
  integer :: in, jn, kn, dim, i, j, k, l
  real(8), allocatable :: coeffs(:)

  ! Set grid parameters.
  cg%is = is;  cg%ie = ie
  cg%js = js;  cg%je = je
  cg%ks = ks;  cg%ke = ke

  in = ie - is + 1
  jn = je - js + 1
  kn = ke - ks + 1
  cg%in = in;  cg%jn = jn;  cg%kn = kn
  lmax = in * jn * kn
  cg%lmax = lmax

  ! Determine problem dimensionality.
  if(ie > is .and. je > js .and. ke > ks) then
    dim = 3
  else if(ie > is .and. (je > js .or. ke > ks)) then
    dim = 2
  else if(ie > is .and. (je == js .and. ke == ks)) then
    dim = 1
  else
    print *, 'Error in setup_grvcg_new: Unsupported dimension'
    stop
  end if

  ! Set number of diagonals and coefficient count.
  select case(dim)
  case(1)
    cg%Adiags = 2
  case(2)
    cg%Adiags = 3
  case(3)
    cg%Adiags = 5
  end select

  ! Allocate the equation matrix and set the offsets.
  allocate(cg%ia(1:cg%Adiags))
  allocate(cg%A(1:cg%Adiags, 1:lmax))
  select case(dim)
  case(1)
    cg%ia(1) = 0
    cg%ia(2) = 1
  case(2)
    cg%ia(1) = 0
    cg%ia(2) = 1
    cg%ia(3) = in
  case(3)
    cg%ia(1) = 0
    cg%ia(2) = 1
    cg%ia(3) = in
    cg%ia(4) = in * jn
    cg%ia(5) = in * jn * (kn - 1)
  end select

  ! Loop over all grid points. For each point, compute the stencil
  ! coefficients from the physics via compute_coeffs.
  allocate(coeffs(cg%Adiags))
  do l = 1, lmax
    call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
    call compute_coeffs_gravity(dim, i, j, k, ie, kn, coeffs)
    ! Save the computed coefficients into the equation matrix.
    cg%A(:, l) = coeffs
  end do
  deallocate(coeffs)

  ! Set up the preconditioner.
  select case(dim)
  case(1)
    cg%cdiags = 2
  case(2)
    cg%cdiags = 4
  case(3)
    cg%cdiags = 5
  end select

  allocate(cg%ic(1:cg%cdiags))
  allocate(cg%c(1:cg%cdiags, 1:lmax))

  select case(dim)
  case(1)
    cg%ic(1) = 0
    cg%ic(2) = 1
    cg%alpha = 0.99d0
  case(2)
    cg%ic(1) = 0
    cg%ic(2) = 1
    cg%ic(3) = in - 1
    cg%ic(4) = in
    cg%alpha = 0.99d0
  case(3)
    cg%ic(1) = 0
    cg%ic(2) = 1
    cg%ic(3) = in
    cg%ic(4) = in * jn
    cg%ic(5) = in * jn * (kn - 1)
    cg%alpha = 0.999d0
  end select

  call get_preconditioner(cg)

end subroutine setup_cg

end module matrix_solver