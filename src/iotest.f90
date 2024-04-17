module iotest_mod
  implicit none

  contains
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                           SUBROUTINE IOTEST
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: To test the read and write routines

  subroutine iotest
    use mpi_utils
    use grid
    use physval
    use output_mod
    use profiler_mod
    use readbin_mod

    integer:: i,j,k,numerr
    real(8):: err

    if (myrank == 0) print*, 'Running I/O test'

    ! MHD unused
    ! TODO: test these quantities
    b1 = 0.d0
    b2 = 0.d0
    b3 = 0.d0

    ! Initialize the arrays
    do i = is, ie
      do j = js, je
        do k = ks, ke
          d(i,j,k)  = 1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v1(i,j,k) = 2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v2(i,j,k) = 3.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v3(i,j,k) = 4.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          e(i,j,k)  = 1.d3 + 1.d-2*i + 1.d-4*j + 1.d-6*k
        end do
      end do
    end do

    ! Override the profiler time to prevent a divide by zero error during output
    wtime(wtlop) = 1.d0


    ! Write the arrays to file
    call output

    ! Reset the arrays in memory
    d = 0.d0
    v1 = 0.d0
    v2 = 0.d0
    v3 = 0.d0
    e = 0.d0

    ! Read the arrays from file
    call readbin('data/bin00000000000s.dat')

    ! Check the arrays
    numerr = 0
    do i = is, ie
      do j = js, je
        do k = ks, ke
          err = 0.d0
          err = err + abs(d(i,j,k)  - (1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(v1(i,j,k) - (2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(v2(i,j,k) - (3.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(v3(i,j,k) - (4.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(e(i,j,k)  - (1.d3 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          if (err > 0.d0) then
            print*, 'Error in array values at i=',i,'j=',j,'k=',k,'err=',err
            numerr = numerr + 1
          endif
        end do
      end do
    end do

    call allreduce_mpi('sum', numerr)

    if (myrank == 0) then
      if (numerr > 0) then
        print*, 'Number of errors = ',numerr
        print*, 'Test failed'
        error stop 1
      else
        print*, 'Test passed'
      endif
    endif

    call finalize_mpi

    stop

    return
  end subroutine iotest

end module iotest_mod
