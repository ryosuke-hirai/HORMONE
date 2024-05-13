
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

   subroutine domain_decomp
      use grid
#ifdef MPI
      use settings, only:gravswitch

      integer :: nx, ny, nz
      integer :: ierr
      integer, allocatable :: factors(:)
      integer :: num_factors, axis
      real(8) :: n_tmp(3)
      integer, dimension(3) :: dims
      logical :: periods(3)
      integer :: mycoords(3)
      integer :: i

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

      ! Prime factors of the number of MPI tasks
      call prime_factors(nprocs, factors, num_factors)

      ! For each factor, split the domain along the dimension with the most cells
      dims = [1, 1, 1]
      n_tmp = [real(nx), real(ny), real(nz)] ! Real, because cells may not divide exactly
      do i = 1, num_factors
         axis = maxloc(n_tmp, 1) ! Index of the largest value in n_tmp
         dims(axis) = dims(axis) * factors(i)
         n_tmp(axis) = n_tmp(axis) / real(factors(i))
      enddo

      ! Always set to periodic and allow boundary conditions to override
      periods = [.true., .true., .true.]

      call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, .false., cart_comm, ierr)
      call MPI_Cart_coords(cart_comm, myrank, 3, mycoords, ierr)

      ! Compute local domain, including offset if starting index is not 1
      is = mycoords(1) * (nx / dims(1)) + 1 + (is_global - 1)
      ie = (mycoords(1) + 1) * (nx / dims(1)) + (is_global - 1)
      js = mycoords(2) * (ny / dims(2)) + 1 + (js_global - 1)
      je = (mycoords(2) + 1) * (ny / dims(2)) + (js_global - 1)
      ks = mycoords(3) * (nz / dims(3)) + 1 + (ks_global - 1)
      ke = (mycoords(3) + 1) * (nz / dims(3)) + (ks_global - 1)

      ! Special treatment of final proc, in case nx/ny/nz is not divisible by nprocs
      if (mycoords(1) == dims(1) - 1) ie = nx + (is_global - 1)
      if (mycoords(2) == dims(2) - 1) je = ny + (js_global - 1)
      if (mycoords(3) == dims(3) - 1) ke = nz + (ks_global - 1)

      ! Copy the hydro grid decomposition to the gravity grid
      gis = is
      gie = ie
      gjs = js
      gje = je
      gks = ks
      gke = ke

      ! Check for gravity grid extension and adjust if necessary
      if (gravswitch >= 1) then
         if (mycoords(1) == 0) then
            if (gis_global /= is_global) gis = gis_global
         endif
         if (mycoords(1) == dims(1) - 1) then
            if (gie_global /= ie_global) gie = gie_global
         endif
         if (mycoords(2) == 0) then
            if (gjs_global /= js_global) gjs = gjs_global
         endif
         if (mycoords(2) == dims(2) - 1) then
            if (gje_global /= je_global) gje = gje_global
         endif
         if (mycoords(3) == 0) then
            if (gks_global /= ks_global) gks = gks_global
         endif
         if (mycoords(3) == dims(3) - 1) then
            if (gke_global /= ke_global) gke = gke_global
         endif
      endif

      call setup_mpi_exchange
      call setup_mpi_io
      call print_domain_decomposition

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
      use profiler_mod
      use mpi_utils

      integer :: i

      ! Timing: measure the time spent waiting for other tasks to catch up
      call start_clock(wtwai)
      call barrier_mpi
      call stop_clock(wtwai)

      call start_clock(wtmpi)

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

      call stop_clock(wtmpi)

   end subroutine exchange_mpi

   subroutine exchange_scalar(val)
      use grid

      real(8), intent(inout) :: val(is-2:ie+2,js-2:je+2,ks-2:ke+2)

#ifdef MPI
      integer :: d
      integer :: ierr
      integer :: i, n

      do d = 1, 3
         n = 1
         if (d == 1) then
            if (is_global == ie_global) cycle
            if (is==is_global .and. ie==ie_global) cycle
            if (is == ie) n = 2
         else if (d == 2) then
            if (js_global == je_global) cycle
            if (js==js_global .and. je==je_global) cycle
            if (js == je) n = 2
         else if (d == 3) then
            if (ks_global == ke_global) cycle
            if (ks==ks_global .and. ke==ke_global) cycle
            if (ks == ke) n = 2
         endif

         ! If the domain is only one cell wide in this direction, exchange twice
         ! so that the ghost cells (two deep) are propagated correctly
         do i = 1, n

         ! Send left real cells to left neighbour's right ghost cells
         call MPI_Sendrecv(val(l_real (1,d), l_real (2,d), l_real (3,d)), 1, subarray(d), left_rank (d), 0, &
                           val(r_ghost(1,d), r_ghost(2,d), r_ghost(3,d)), 1, subarray(d), right_rank(d), 0, &
                           cart_comm, MPI_STATUS_IGNORE, ierr)

         ! Send right real cells to right neighbour's left ghost cells
         call MPI_Sendrecv(val(r_real (1,d), r_real (2,d), r_real (3,d)), 1, subarray(d), right_rank(d), 0, &
                           val(l_ghost(1,d), l_ghost(2,d), l_ghost(3,d)), 1, subarray(d), left_rank (d), 0, &
                           cart_comm, MPI_STATUS_IGNORE, ierr)
         enddo
      enddo
#endif
   end subroutine exchange_scalar

   subroutine setup_mpi_io
#ifdef MPI
      use mpi
      use mpi_utils, only: mpi_subarray_default, mpi_subarray_spc, mpi_subarray_gravity, mpi_subarray_extgrtv
      use grid
      use settings, only: spn, compswitch, include_sinks, gravswitch, include_extgrv
      integer, dimension(3) :: sizes3, subsizes3, starts3
      integer, dimension(4) :: sizes4, subsizes4, starts4
      integer :: ierr, istart, iend, jstart, jend, kstart, kend

      ! Set up the subarray which selects only the real cells for I/O
      sizes3 = [ie_global-is_global+1,je_global-js_global+1,ke_global-ks_global+1]
      subsizes3 = [ie-is+1,je-js+1,ke-ks+1]
      starts3 = [is-is_global,js-js_global,ks-ks_global]
      call mpi_type_create_subarray(3, sizes3, subsizes3, starts3, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_subarray_default, ierr)
      call mpi_type_commit(mpi_subarray_default, ierr)

      if (compswitch>=2) then
         ! Set up the subarray for spc
         sizes4 = [spn,ie_global-is_global+1,je_global-js_global+1,ke_global-ks_global+1]
         subsizes4 = [spn,ie-is+1,je-js+1,ke-ks+1]
         starts4 = [0,is-is_global,js-js_global,ks-ks_global]
         call mpi_type_create_subarray(4, sizes4, subsizes4, starts4, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_subarray_spc, ierr)
         call mpi_type_commit(mpi_subarray_spc, ierr)
      endif

      if (include_sinks) then
         call create_sink_type_mpi
      endif

      if (gravswitch>=1) then
         sizes3 = [gie_global-gis_global+1,gje_global-gjs_global+1,gke_global-gks_global+1]
         subsizes3 = [gie-gis+1,gje-gjs+1,gke-gks+1]
         starts3 = [gis-gis_global,gjs-gjs_global,gks-gks_global]
         call mpi_type_create_subarray(3, sizes3, subsizes3, starts3, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_subarray_gravity, ierr)
         call mpi_type_commit(mpi_subarray_gravity, ierr)
      end if

      if (include_extgrv) then
         istart = gis; iend = gie
         jstart = gjs; jend = gje
         kstart = gks; kend = gke
         if (gis==gis_global) istart = gis-2
         if (gie==gie_global) iend = gie+2
         if (gjs==gjs_global) jstart = gjs-2
         if (gje==gje_global) jend = gje+2
         if (gks==gks_global) kstart = gks-2
         if (gke==gke_global) kend = gke+2
         sizes3 = [(gie_global+2)-(gis_global-2)+1, (gje_global+2)-(gjs_global-2)+1, (gke_global+2)-(gks_global-2)+1]
         subsizes3 = [iend-istart+1, jend-jstart+1, kend-kstart+1]
         starts3 = [istart-(gis_global-2), jstart-(gjs_global-2), kstart-(gks_global-2)]
         call mpi_type_create_subarray(3, sizes3, subsizes3, starts3, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_subarray_extgrtv, ierr)
         call mpi_type_commit(mpi_subarray_extgrtv, ierr)
      endif

#endif
   end subroutine setup_mpi_io

   subroutine create_sink_type_mpi
#ifdef MPI
      use mpi_utils, only: mpi_type_sink_prop
      use sink_mod, only: sink_prop
      type(sink_prop) :: sink
      integer, parameter :: nattr = 12
      integer, dimension(nattr) :: types, blocklengths
      integer(kind=MPI_ADDRESS_KIND), dimension(nattr) :: offsets
      integer :: ierr

      ! Get the addresses of each component of sink_prop
      call MPI_Get_address(sink%i, offsets(1), ierr)
      call MPI_Get_address(sink%j, offsets(2), ierr)
      call MPI_Get_address(sink%k, offsets(3), ierr)
      call MPI_Get_address(sink%mass, offsets(4), ierr)
      call MPI_Get_address(sink%softfac, offsets(5), ierr)
      call MPI_Get_address(sink%lsoft, offsets(6), ierr)
      call MPI_Get_address(sink%locres, offsets(7), ierr)
      call MPI_Get_address(sink%dt, offsets(8), ierr)
      call MPI_Get_address(sink%x, offsets(9), ierr)
      call MPI_Get_address(sink%v, offsets(10), ierr)
      call MPI_Get_address(sink%a, offsets(11), ierr)
      call MPI_Get_address(sink%xpol, offsets(12), ierr)

      ! Compute offsets as relative to the start of sink
      offsets = offsets - offsets(1)

      ! Set the blocklengths (number of elements in each block)
      blocklengths = (/1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3/)

      ! Set the types (type of each block)
      types = (/MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, \
                MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, \
                MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION/)

      ! Create the custom datatype
      call MPI_Type_create_struct(nattr, blocklengths, offsets, types, mpi_type_sink_prop, ierr)
      call MPI_Type_commit(mpi_type_sink_prop, ierr)
#endif
   end subroutine create_sink_type_mpi

   subroutine prime_factors(n, factors, num_factors)
      ! Compute the prime factors of a number
      implicit none
      integer, intent(in) :: n
      integer, dimension(:), allocatable, intent(out) :: factors
      integer, intent(out) :: num_factors
      integer :: i, count, temp

      ! Count the number of factors
      temp = n
      count = 0
      do i = 2, n
         do while (mod(temp, i) == 0)
            temp = temp / i
            count = count + 1
         enddo
        if (temp == 1) exit
      end do

      ! Repeat, saving the factors
      allocate(factors(count))
      temp = n
      count = 0
      do i = 2, n
         do while (mod(temp, i) == 0)
            temp = temp / i
            count = count + 1
            factors(count) = i
         enddo
        if (temp == 1) exit
      end do

      num_factors = count
   end subroutine prime_factors

   subroutine calculate_decomposition_efficiency(eff)
      use settings
      use grid
      real(8), intent(out) :: eff
      ! Calculate the ratio: (real cells) / (total cells including ghost)
      eff = 0.d0
      if (solve_i) then
         eff = eff + real(2 * (je - js + 1) * (ke - ks + 1))
      endif
      if (solve_j) then
         eff = eff + real(2 * (ie - is + 1) * (ke - ks + 1))
      endif
      if (solve_k) then
         eff = eff + real(2 * (ie - is + 1) * (je - js + 1))
      endif
      eff = 1.d0 - (eff / real((ie - is + 5) * (je - js + 5) * (ke - ks + 5)))
   end subroutine calculate_decomposition_efficiency

   subroutine print_domain_decomposition
#ifdef MPI
      use grid
      real(8) :: eff
      integer :: i, ierr

      call calculate_decomposition_efficiency(eff)

      if (myrank == 0) then
         write(*,'(7(A,I0),A)') 'Global domain (', is_global, ':', ie_global, ', ', js_global, ':', je_global, ', ', ks_global, ':', ke_global, ') split between ', nprocs, ' MPI ranks:'
      endif
      call MPI_Barrier(cart_comm, ierr)

      do i = 1, nprocs
         if (myrank == i-1) then
            write(*,'(7(A,I0),A,F6.2,A)') '  Rank ', myrank, ' has domain (', is, ':', ie, ', ', js, ':', je, ', ', ks, ':', ke, '), volume efficiency=', eff*100.d0, '%'
         endif
         call MPI_Barrier(cart_comm, ierr)
      enddo

      if (myrank == 0) then
         write(*, *)
      endif
      call MPI_Barrier(cart_comm, ierr)
#endif
   end subroutine print_domain_decomposition

end module mpi_domain
