module matrix_utils
  use miccg_mod, only: cg_set, ijk_from_l, get_preconditioner  ! setup_grvcg already exists
  use gravity_miccg_mod, only: compute_coeffs, compute_precond
#ifdef USE_PETSC
  use petsc, only: MatSetValue, MatAssemblyBegin, MatAssemblyEnd, &
                     MAT_INSERT_VALUES, MAT_FINAL_ASSEMBLY
#endif
  implicit none
  public :: setup_grvA
contains

subroutine setup_grvA(is, ie, js, je, ks, ke)
    !-----------------------------------------------------------------------
    ! Top-level assembly routine.
    ! It simply calls a branch:
    !   - For MICCG: setup_grvcg (which builds the CDS-based cg object)
    !   - For PETSc: setup_petsc_A (which writes directly into PETSc arrays)
    !-----------------------------------------------------------------------
    integer, intent(in) :: is, ie, js, je, ks, ke

#ifdef USE_PETSC
        call setup_petsc_A(is, ie, js, je, ks, ke)
#else
        call setup_cg(is, ie, js, je, ks, ke)
#endif
end subroutine setup_grvA

#ifdef USE_PETSC
  subroutine setup_petsc_A(is, ie, js, je, ks, ke)
    !-----------------------------------------------------------------------
    ! PETSc branch for assembling the Laplacian matrix.
    ! For each grid point, we compute the stencil using compute_coeffs,
    ! then insert the values directly into the PETSc matrix (A_petsc).
    ! Preconditioning is handled by PETSc, so we do not compute it here.
    !-----------------------------------------------------------------------
    integer, intent(in) :: is, ie, js, je, ks, ke
    integer :: in, jn, kn, lmax, dim, i, j, k, l, ncoeff
    real(8), allocatable :: coeffs(:)
    integer :: ierr

    ! Compute grid dimensions and total number of points.
    in = ie - is + 1
    jn = je - js + 1
    kn = ke - ks + 1
    lmax = in * jn * kn

    ! Determine problem dimension based on input indices.
    if(ie > is .and. je > js .and. ke > ks) then
      dim = 3
    else if(ie > is .and. (je > js .or. ke > ks)) then
      dim = 2
    else if(ie > is .and. (je == js .and. ke == ks)) then
      dim = 1
    else
      print *, 'Error in setup_petsc_A: unsupported dimension'
      stop
    end if

    ! Loop over all grid points.
    do l = 1, lmax
      call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
      select case(dim)
      case(1)
         ncoeff = 2
      case(2)
         ncoeff = 3
      case(3)
         ncoeff = 5
      end select
      allocate(coeffs(ncoeff))
      call compute_coeffs(dim, i, j, k, ie, kn, coeffs)
      ! Insert the computed coefficients into the PETSc matrix.
      ! Here we assume each grid point corresponds to a row (l-1).
      ! Insert the main (diagonal) coefficient.
      call MatSetValue(A_petsc, l-1, l-1, coeffs(1), MAT_INSERT_VALUES, ierr)
      if(ncoeff >= 2) then
         ! For example, assume offset of 1 for the first off-diagonal.
         if(l > 1) then
           call MatSetValue(A_petsc, l-1, l-2, coeffs(2), MAT_INSERT_VALUES, ierr)
         end if
      end if
      if(ncoeff >= 3) then
         ! Example using offset equal to in for the second off-diagonal.
         if(l + in <= lmax) then
           call MatSetValue(A_petsc, l-1, (l-1) + in, coeffs(3), MAT_INSERT_VALUES, ierr)
         end if
      end if
      if(ncoeff >= 4) then
         ! For 3D: further off-diagonals (offset = in*jn).
         if(l + in*jn <= lmax) then
           call MatSetValue(A_petsc, l-1, (l-1) + in*jn, coeffs(4), MAT_INSERT_VALUES, ierr)
         end if
      end if
      if(ncoeff >= 5) then
         ! And the last off-diagonal (offset = in*jn*(kn-1))
         if(l + in*jn*(kn-1) <= lmax) then
           call MatSetValue(A_petsc, l-1, (l-1) + in*jn*(kn-1), coeffs(5), MAT_INSERT_VALUES, ierr)
         end if
      end if
      deallocate(coeffs)
    end do

    ! Finalize PETSc matrix assembly.
    call MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A_petsc, MAT_FINAL_ASSEMBLY, ierr)
  end subroutine setup_petsc_A
#endif

subroutine setup_cg(is, ie, js, je, ks, ke)
    !-------------------------------------------------------------
    ! New MICCG assembly routine using compute_coeffs and
    ! compute_precond. This routine fills the CDS-based cg object.
    !
    ! Arguments:
    !   is,ie,js,je,ks,ke : grid indices in the three directions
    !-------------------------------------------------------------
    use miccg_mod, only: cg=>cg_grv
    integer, intent(in) :: is, ie, js, je, ks, ke
    integer :: in, jn, kn, lmax, dim, i, j, k, l, ncoeff
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
      ncoeff = 2
    case(2)
      cg%Adiags = 3
      ncoeff = 3
    case(3)
      cg%Adiags = 5
      ncoeff = 5
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
    do l = 1, lmax
       call ijk_from_l(l, is, js, ks, in, jn, i, j, k)
       allocate(coeffs(ncoeff))
       call compute_coeffs(dim, i, j, k, ie, kn, coeffs)
       ! Save the computed coefficients into the equation matrix.
       cg%A(:, l) = coeffs
       deallocate(coeffs)
    end do

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
    call compute_precond(dim, in, jn, kn, cg%cdiags, cg%ic, cg%alpha)
    call get_preconditioner(cg)

  end subroutine setup_cg

end module matrix_utils