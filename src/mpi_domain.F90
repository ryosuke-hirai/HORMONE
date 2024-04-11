
module mpi_domain
#ifdef MPI
   use mpi
#endif
   use mpi_utils, only: myrank, nprocs

   implicit none

#ifdef MPI
   integer :: cart_comm

   ! Indices for the real and ghost zones involved in the exchange
   ! (2nd index j is for exchange in the j-direction)
   integer :: l_ghost(3,3), r_ghost(3,3), l_real(3,3), r_real(3,3)

   ! MPI subarray datatype for exchange
   integer :: subarray(3)

   ! Left and right neighbours in each direction
   integer :: left_rank(3), right_rank(3)
#endif

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

      call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, .false., cart_comm, ierr)
      call MPI_Cart_coords(cart_comm, myrank, 3, mycoords, ierr)

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

      call setup_mpi_exchange

#endif

   end subroutine domain_decomp

   subroutine setup_mpi_exchange
#ifdef MPI
      use grid

      ! Indices and sizes for MPI subarray datatypes used for exchange in each direction
      integer :: sizes(3), subsizes(3), starts(3)
      integer :: d
      integer :: ierr

      ! Set up the subarrays for MPI exchange, which only needs to be done once at the start

      ! The exchange populates ghost cells with the real cells of neighbours,
      ! and send real cells to the ghost cells of neighbours
      ! E.g. In the x-direction:
      ! Send val(is:is+1,:,:) to the left rank's val(ie+1:ie+2,:,:)
      ! Recv val(ie+1:ie+2,:,:) from the right rank's val(is:is+1,:,:)
      ! Send val(ie-1:ie,:,:) to the right rank's val(is-2:is-1,:,:)
      ! Recv val(is-2:is-1,:,:) from the left rank's val(ie-1:ie,:,:)

      ! Rather than slicing the array manually into many MPI_SENDRECV calls,
      ! create a subarray datatype to send a subset of the array

      ! The send and recv buffers in MPI_Sendrecv must be disjoint, so the subarrays have to be created
      ! relative to the array indices passed in, not the start of the array

      ! Get the left and right neighbours in each direction
      do d = 1, 3
         call MPI_Cart_shift(cart_comm, d-1, 1, left_rank(d), right_rank(d), ierr)
      enddo

      ! Size of the array on this task
      sizes = [ie - is + 5, je - js + 5, ke - ks + 5]

      ! --- x-1 direction ---
      subsizes = [2, je-js+1, ke-ks+1] ! Size of the ghost cells to send
      starts   = [0, 2, 2] ! Offset relative to the address passed to MPI_Sendrecv
      call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL8, subarray(1), ierr)
      call MPI_Type_commit(subarray(1), ierr)

      ! Starting indices of the real and ghost zones involved in the exchange
      l_real (:,1) = [is,   js-2, ks-2]
      r_real (:,1) = [ie-1, js-2, ks-2]
      l_ghost(:,1) = [is-2, js-2, ks-2]
      r_ghost(:,1) = [ie+1, js-2, ks-2]

      ! --- x-2 direction ---
      subsizes = [ie-is+1, 2, ke-ks+1] ! Size of the ghost cells to send
      starts   = [2, 0, 2] ! Offset relative to the address passed to MPI_Sendrecv
      call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL8, subarray(2), ierr)
      call MPI_Type_commit(subarray(2), ierr)

      ! Starting indices of the real and ghost zones involved in the exchange
      l_real (:,2) = [is-2, js,   ks-2]
      r_real (:,2) = [is-2, je-1, ks-2]
      l_ghost(:,2) = [is-2, js-2, ks-2]
      r_ghost(:,2) = [is-2, je+1, ks-2]

      ! --- x-3 direction ---
      subsizes = [ie-is+1, je-js+1, 2] ! Size of the ghost cells to send
      starts   = [2, 2, 0] ! Offset relative to the address passed to MPI_Sendrecv
      call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL8, subarray(3), ierr)
      call MPI_Type_commit(subarray(3), ierr)

      ! Starting indices of the real and ghost zones involved in the exchange
      l_real (:,3) = [is-2, js-2, ks  ]
      r_real (:,3) = [is-2, js-2, ke-1]
      l_ghost(:,3) = [is-2, js-2, ks-2]
      r_ghost(:,3) = [is-2, js-2, ke+1]

#endif
   end subroutine setup_mpi_exchange

   subroutine exchange_mpi
      use settings
      use grid
      use physval
      integer :: i
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
      real(8), intent(inout) :: val(is-2:ie+2,js-2:je+2,ks-2:ke+2)

#ifdef MPI
      integer :: d
      integer :: ierr

      do d = 1, 3
         if (d == 1) then
            if (is_global == ie_global) cycle
            if (is==is_global .and. ie==ie_global) cycle
         else if (d == 2) then
            if (js_global == je_global) cycle
            if (js==js_global .and. je==je_global) cycle
         else if (d == 3) then
            if (ks_global == ke_global) cycle
            if (ks==ks_global .and. ke==ke_global) cycle
         endif

         ! Send left real cells to left neighbour's right ghost cells
         call MPI_Sendrecv(val(l_real (1,d), l_real (2,d), l_real (3,d)), 1, subarray(d), left_rank (d), 0, &
                           val(r_ghost(1,d), r_ghost(2,d), r_ghost(3,d)), 1, subarray(d), right_rank(d), 0, &
                           cart_comm, MPI_STATUS_IGNORE, ierr)

         ! Send right real cells to right neighbour's left ghost cells
         call MPI_Sendrecv(val(r_real (1,d), r_real (2,d), r_real (3,d)), 1, subarray(d), right_rank(d), 0, &
                           val(l_ghost(1,d), l_ghost(2,d), l_ghost(3,d)), 1, subarray(d), left_rank (d), 0, &
                           cart_comm, MPI_STATUS_IGNORE, ierr)
      enddo
#endif
   end subroutine exchange_scalar

end module mpi_domain
