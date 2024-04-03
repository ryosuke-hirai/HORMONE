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

  subroutine domain_decomp()
    use grid
#ifdef MPI
    integer :: nx, ny, nz
    integer :: ierr
    integer :: dims(3), coords(3)
    logical :: periods(3)
    integer :: cart_comm, mycoords(3)
#endif

    is_global = is
    ie_global = ie
    js_global = js
    je_global = je
    ks_global = ks
    ke_global = ke

    gis_global = gis
    gie_global = gie
    gjs_global = gjs
    gje_global = gje
    gks_global = gks
    gke_global = gke

#ifdef MPI
    nx = ie - is + 1
    ny = je - js + 1
    nz = ke - ks + 1

    dims = [nprocs, 1, 1] ! Trivial decomposition, for now... TODO
    periods = [.false., .false., .false.] ! TODO

    call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, .false., cart_comm, ierr)
    call MPI_CART_COORDS(cart_comm, myrank, 3, mycoords, ierr)

    ! Compute local domain, including offset if starting index is not 1
    is = mycoords(1) * (nx / dims(1)) + 1 - (is_global - 1)
    ie = (mycoords(1) + 1) * (nx / dims(1)) - (is_global - 1)
    js = mycoords(2) * (ny / dims(2)) + 1 - (js_global - 1)
    je = (mycoords(2) + 1) * (ny / dims(2)) - (js_global - 1)
    ks = mycoords(3) * (nz / dims(3)) + 1 - (ks_global - 1)
    ke = (mycoords(3) + 1) * (nz / dims(3)) - (ks_global - 1)

    ! Special treatment of final proc, in case nx/ny/nz is not divisible by nprocs
    if (mycoords(1) == dims(1) - 1) ie = nx - (is_global - 1)
    if (mycoords(2) == dims(2) - 1) je = ny - (js_global - 1)
    if (mycoords(3) == dims(3) - 1) ke = nz - (ks_global - 1)

#endif

  end subroutine domain_decomp

end module mpi_utils
