module matrix_utils
  use miccg_mod, only: cg_set, ijk_from_l, get_preconditioner
  use gravity_miccg_mod, only: compute_coeffs

  implicit none
  public :: setup_grvA

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
!
! TODO: not fully implemented yet

subroutine setup_petsc_A

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
  integer :: in, jn, kn, lmax, dim, i, j, k, l
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
    call compute_coeffs(dim, i, j, k, ie, kn, coeffs)
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

end module matrix_utils