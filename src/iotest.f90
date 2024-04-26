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
    use settings, only: spn
    use sink_mod, only: sink, nsink, sink_prop

    integer:: a,i,j,k,numerr
    real(8):: err
    character(len=4) :: speci

    if (myrank == 0) print*, 'Running I/O test'

    ! Initialize the arrays
    do i = is, ie
      do j = js, je
        do k = ks, ke
          d(i,j,k)   = 1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v1(i,j,k)  = 2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v2(i,j,k)  = 3.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v3(i,j,k)  = 4.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          b1(i,j,k)  = 5.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          b2(i,j,k)  = 6.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          b3(i,j,k)  = 7.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          phi(i,j,k) = 8.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          e(i,j,k)   = 1.d3 + 1.d-2*i + 1.d-4*j + 1.d-6*k
        end do
      end do
    end do

    do a = 1, spn
      do i = is, ie
        do j = js, je
          do k = ks, ke
            spc(a,i,j,k) = 8.d0 + 1.d0*a + 1.d-2*i + 1.d-4*j + 1.d-6*k
          end do
        end do
      end do
      write(speci, '(I0)') a
      species(a) = 'spc'//trim(speci)
    end do

    do i=1,nsink
      sink(i)%i = 1*i
      sink(i)%j = 10*i
      sink(i)%k = 100*i
      sink(i)%mass = 1.d0*i
      sink(i)%softfac = 10.d0*i
      sink(i)%lsoft = 100.d0*i
      sink(i)%locres = 1000.d0*i
      sink(i)%dt = 10000.d0*i
      sink(i)%x = (/1.d0*i,10.d0*i,100.d0*i/)
      sink(i)%v = (/2.d0*i,20.d0*i,200.d0*i/)
      sink(i)%a = (/3.d0*i,30.d0*i,300.d0*i/)
      sink(i)%xpol = (/4.d0*i,40.d0*i,400.d0*i/)
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
    b1 = 0.d0
    b2 = 0.d0
    b3 = 0.d0
    phi = 0.d0
    e = 0.d0

    spc = 0.d0
    species = ''

    do i = 1, size(sink)
      sink(i)%i = 0
      sink(i)%j = 0
      sink(i)%k = 0
      sink(i)%mass = 0.d0
      sink(i)%softfac = 0.d0
      sink(i)%lsoft = 0.d0
      sink(i)%locres = 0.d0
      sink(i)%dt = 0.d0
      sink(i)%x = 0.d0
      sink(i)%v = 0.d0
      sink(i)%a = 0.d0
      sink(i)%xpol = 0.d0
    end do

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
          err = err + abs(b1(i,j,k) - (5.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(b2(i,j,k) - (6.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(b3(i,j,k) - (7.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(phi(i,j,k)- (8.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          err = err + abs(e(i,j,k)  - (1.d3 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          if (err > 0.d0) then
            print*, 'Error in array values at i=',i,'j=',j,'k=',k,'err=',err
            numerr = numerr + 1
          endif
        end do
      end do
    end do

    do a = 1,spn
      do i = is, ie
        do j = js, je
          do k = ks, ke
            err = 0.d0
            err = err + abs(spc(a,i,j,k) - (8.d0 + 1.d0*a + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            if (err > 0.d0) then
              print*, 'Error in spc values at a=',a,'i=',i,'j=',j,'k=',k,'err=',err
              numerr = numerr + 1
            endif
          end do
        end do
      end do
      write(speci, '(I0)') a
      if (trim(species(a)) /= 'spc'//trim(speci)) then
        print*, 'Error in species names at a=',a,'species(a)=',trim(species(a))
        numerr = numerr + 1
      end if
    end do

    do i = 1, nsink
      err = 0.d0
      err = err + abs(sink(i)%i - 1*i)
      err = err + abs(sink(i)%j - 10*i)
      err = err + abs(sink(i)%k - 100*i)
      err = err + abs(sink(i)%mass - 1.d0*i)
      err = err + abs(sink(i)%softfac - 10.d0*i)
      err = err + abs(sink(i)%lsoft - 100.d0*i)
      err = err + abs(sink(i)%locres - 1000.d0*i)
      err = err + abs(sink(i)%dt - 10000.d0*i)
      err = err + sum(abs(sink(i)%x - (/1.d0*i, 10.d0*i, 100.d0*i/)))
      err = err + sum(abs(sink(i)%v - (/2.d0*i, 20.d0*i, 200.d0*i/)))
      err = err + sum(abs(sink(i)%a - (/3.d0*i, 30.d0*i, 300.d0*i/)))
      err = err + sum(abs(sink(i)%xpol - (/4.d0*i, 40.d0*i, 400.d0*i/)))
      if (err > 0.d0) then
        print*, 'Error in sink array values at i=', i, 'err=', err
        numerr = numerr + 1
      end if
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
