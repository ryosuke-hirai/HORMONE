
module mpi_domain
#ifdef MPI
   use mpi
#endif
   use mpi_utils, only: myrank, nprocs

   implicit none

   integer :: cart_comm

   contains

   subroutine domain_decomp()
      use grid
#ifdef MPI
      integer :: nx, ny, nz
      integer :: ierr
      integer :: dims(3), coords(3)
      logical :: periods(3)
      integer :: mycoords(3)
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

      ! Trivial decomposition for now
      ! Slice along the dimension with most cells to ensure that all 1D problems work
      ! TODO: decompose along all 3 dimensions
      if (nx >= ny .and. nx >= nz) then
         dims = [nprocs, 1, 1]
      else if (ny >= nx .and. ny >= nz) then
         dims = [1, nprocs, 1]
      else
         dims = [1, 1, nprocs]
      endif

      periods = [.true., .true., .true.] ! Always set to periodic and allow boundary conditions to override

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

      if (myrank == 0) then
         print*, 'Global domain: ', is_global, ie_global, js_global, je_global, ks_global, ke_global
      endif

      print*, 'Rank ', myrank, ' has domain ', is, ie, js, je, ks, ke

#endif

   end subroutine domain_decomp

   subroutine exchange_mpi
      use settings
      use grid
      use physval
      integer :: i, left_rank, right_rank, ierr
      ! Exchange data between MPI domains
      ! Scalar quantities: d, p, phi, spc
      ! Vector quantities: v1, v2, v3, b1, b2, b3

      ! Always perform periodic exchange, and allow the boundary condition
      ! routine to override if necessary

      call exchange_scalar(d)
      call exchange_scalar(p)
      call exchange_scalar(phi)
      do i = 1, spn
         call exchange_scalar(spc(i,:,:,:)) ! TODO: This is inefficient, but works for now
      enddo
      call exchange_scalar(v1)
      call exchange_scalar(v2)
      call exchange_scalar(v3)
      call exchange_scalar(b1)
      call exchange_scalar(b2)
      call exchange_scalar(b3)

      call exchange_scalar(e)
      call exchange_scalar(eint)

   end subroutine exchange_mpi

   subroutine exchange_scalar(val)
      use grid
      real(8), intent(in) :: val(is-2:ie+2,js-2:je+2,ks-2:ke+2)

#ifdef MPI
      integer :: i, j, k, d, ghost
      integer :: left_rank, right_rank, ierr

      do d = 1, 3

         if (d == 1) then
            if (is_global == ie_global) cycle
            if (is==is_global .and. ie==ie_global) cycle
            i = 1; j = 0; k = 0
         else if (d == 2) then
            if (js_global == je_global) cycle
            if (js==js_global .and. je==je_global) cycle
            i = 0; j = 1; k = 0
         else if (d == 3) then
            if (ks_global == ke_global) cycle
            if (ks==ks_global .and. ke==ke_global) cycle
            i = 0; j = 0; k = 1
         endif

         call MPI_CART_SHIFT(cart_comm, d-1, 1, left_rank, right_rank, ierr)

         do ghost = 1, 2
            ! As long as the first-depth ghost gets sent first, the second-depth ghost will
            ! always propagate the correct value, even if it originates from two neighbors deep
            ! because the nearest MPI neighbour will have the correct value after the first update
            ! i.e. there is no need to communicate directly with the second-nearest MPI neighbour, even
            ! if the required value is two MPI neighbours away

            ! Send to the left, receive from the right
            call MPI_SENDRECV(val(is + (ghost-1)*i, js + (ghost-1)*j, ks + (ghost-1)*k), 1, MPI_REAL8, left_rank,  0, &
                              val(ie + (ghost  )*i, js + (ghost  )*j, ks + (ghost  )*k), 1, MPI_REAL8, right_rank, 0, &
                              cart_comm, MPI_STATUS_IGNORE, ierr)

            ! Send to the right, receive from the left
            call MPI_SENDRECV(val(ie - (ghost-1)*i, js - (ghost-1)*j, ks - (ghost-1)*k), 1, MPI_REAL8, right_rank, 0, &
                              val(is - (ghost  )*i, js - (ghost  )*j, ks - (ghost  )*k), 1, MPI_REAL8, left_rank,  0, &
                              cart_comm, MPI_STATUS_IGNORE, ierr)
         enddo
      enddo
#endif
   end subroutine exchange_scalar

end module mpi_domain
