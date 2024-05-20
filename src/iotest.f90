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
    use settings, only: spn, include_extgrv
    use sink_mod, only: sink, nsink, sink_prop
    use gravmod, only: grvphi, grvpsi, dt_old
    use utils, only: isequal
    use settings, only: include_sinks, compswitch, mag_on, gravswitch

    integer:: a,i,j,k,numerr=0
    real(8):: err
    character(len=4) :: speci
    character(len=30) :: filename

    if (myrank == 0) print*, 'Running I/O test'

    tn = 123
    time = 456.d0

    call iotest_grid(numerr)

    ! Initialize the arrays
    do i = is, ie
      do j = js, je
        do k = ks, ke
          d(i,j,k)   = 1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v1(i,j,k)  = 2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v2(i,j,k)  = 3.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          v3(i,j,k)  = 4.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          e(i,j,k)   = 1.d3 + 1.d-2*i + 1.d-4*j + 1.d-6*k
        end do
      end do
    end do

    if (mag_on) then
      do i = is, ie
        do j = js, je
          do k = ks, ke
            b1(i,j,k)  = 5.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
            b2(i,j,k)  = 6.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
            b3(i,j,k)  = 7.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
            phi(i,j,k) = 8.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          end do
        end do
      end do
    endif

    if (compswitch>=2) then
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
    endif

    if (include_sinks) then
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
      call open_sinkfile
    endif

    if (gravswitch>=2) then
      do i = gis,gie
        do j = gjs,gje
          do k = gks,gke
            grvphi(i,j,k) = 1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          end do
        end do
      end do
    endif

    if (gravswitch==3) then
      do i = gis,gie
        do j = gjs,gje
          do k = gks,gke
            grvpsi(i,j,k) = 2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
          end do
        end do
      end do
      dt_old = 789.d0
    endif

    ! Override the profiler time to prevent a divide by zero error during output
    wtime(wtlop) = 1.d0

    ! Write the arrays to file
    call output

    ! Get filename, based off time
    call set_file_name('bin', tn, time, filename)

    ! Reset the arrays in memory
    tn = 0
    time = 0.d0
    d = 0.d0
    v1 = 0.d0
    v2 = 0.d0
    v3 = 0.d0
    b1 = 0.d0
    b2 = 0.d0
    b3 = 0.d0
    phi = 0.d0
    e = 0.d0
    grvphi = 0.d0
    grvpsi = 0.d0
    dt_old = 0.d0

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
    if (myrank == 0) print*, 'Reading: ', filename
    call readbin(trim(filename))

    if (tn /= 123) then
      print*, 'Error in tn value, tn=',tn,'should be:',123
      numerr = numerr + 1
    endif

    if ( .not. isequal(time, 456.d0)) then
      print*, 'Error in time value, time=',time,'should be:',456.d0
      numerr = numerr + 1
    endif

    ! Check the arrays
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

    if (mag_on) then
      do i = is, ie
        do j = js, je
          do k = ks, ke
            err = 0.d0
            err = err + abs(b1(i,j,k) - (5.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            err = err + abs(b2(i,j,k) - (6.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            err = err + abs(b3(i,j,k) - (7.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            err = err + abs(phi(i,j,k)- (8.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            if (err > 0.d0) then
              print*, 'Error in array values at i=',i,'j=',j,'k=',k,'err=',err
              numerr = numerr + 1
            endif
          end do
        end do
      end do
    endif

    if (compswitch>=2) then
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
    endif

    if (include_sinks) then
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
    endif

    if (gravswitch>=2) then
      do i = gis,gie
        do j = gjs,gje
          do k = gks,gke
            err = 0.d0
            err = err + abs(grvphi(i,j,k) - (1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            if (err > 0.d0) then
              print*, 'Error in grvphi values at i=',i,'j=',j,'k=',k,'err=',err
              numerr = numerr + 1
            endif
          end do
        end do
      end do
    endif

    if (gravswitch==3) then
      do i = gis,gie
        do j = gjs,gje
          do k = gks,gke
            err = 0.d0
            err = err + abs(grvpsi(i,j,k) - (2.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
            if (err > 0.d0) then
              print*, 'Error in grvpsi values at i=',i,'j=',j,'k=',k,'err=',err
              numerr = numerr + 1
            endif
          end do
        end do
      end do
      if (.not. isequal(dt_old, 789.d0)) then
        print*, 'Error in dt_old value, dt_old=', dt_old, 'should be:', 789.d0
        numerr = numerr + 1
      endif
    endif

    if (include_extgrv) call iotest_extgrv(numerr)

    call allreduce_mpi('sum', numerr)

    if (myrank == 0) then
      if (numerr > 0) then
        print*, 'Number of errors = ',numerr
        print*, 'Test failed'
        error stop 1
      else
        call cleanup_files(tn=123, time=456.d0)
        print*, 'Test passed'
      endif
    endif

    call finalize_mpi

    stop

    return
  end subroutine iotest

  subroutine iotest_extgrv(numerr)
    use grid
    use gravmod, only: extgrv,mc
    use output_mod, only: write_extgrv
    use readbin_mod, only: read_extgrv
    use mpi_utils, only: myrank
    integer, intent(inout) :: numerr
    integer :: i,j,k
    real(8) :: err
    integer :: istart, iend, jstart, jend, kstart, kend

    if (is==is_global) then
      mc(is-1) = 1.d0
    endif

    istart = gis; iend = gie
    jstart = gjs; jend = gje
    kstart = gks; kend = gke
    if (gis==gis_global) istart = gis-2
    if (gie==gie_global) iend = gie+2
    if (gjs==gjs_global) jstart = gjs-2
    if (gje==gje_global) jend = gje+2
    if (gks==gks_global) kstart = gks-2
    if (gke==gke_global) kend = gke+2

    do i = istart, iend
      do j = jstart, jend
        do k = kstart, kend
          extgrv(i,j,k) = 1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k
        enddo
      enddo
    enddo

    call write_extgrv

    mc = -999.d0
    extgrv = -999.d0

    if (myrank == 0) print*, 'Reading: extgrv.bin'
    call read_extgrv('data/extgrv.bin')

    if (is==is_global) then
      err = abs(mc(is-1) - 1.d0)
      if ( err > 0.d0 ) then
        print*, 'Error in mc values at i=',is-1,'err=',err
        numerr = numerr + 1
      endif
    endif

    do i = istart, iend
      do j = jstart, jend
        do k = kstart, kend
          err = abs(extgrv(i,j,k) - (1.d0 + 1.d-2*i + 1.d-4*j + 1.d-6*k))
          if ( err > 0.d0 ) then
            print*, 'Error in extgrv values at i=',i,'j=',j,'k=',k,'err=',err
            numerr = numerr + 1
          endif
        enddo
      enddo
    enddo

  end subroutine iotest_extgrv

  subroutine iotest_grid(numerr)
    use grid
    use output_mod, only: write_grid
    use readbin_mod, only: readgrid
    use utils, only: isequal
    integer, intent(inout) :: numerr
    integer :: i
    real(8),allocatable,dimension(:):: x1_old, xi1_old, dx1_old, dxi1_old
    real(8),allocatable,dimension(:):: x2_old, xi2_old, dx2_old, dxi2_old
    real(8),allocatable,dimension(:):: x3_old, xi3_old, dx3_old, dxi3_old

    allocate(x1_old,mold=x1)
    allocate(x2_old,mold=x2)
    allocate(x3_old,mold=x3)
    allocate(xi1_old,mold=xi1)
    allocate(xi2_old,mold=xi2)
    allocate(xi3_old,mold=xi3)
    allocate(dx1_old,mold=dx1)
    allocate(dx2_old,mold=dx2)
    allocate(dx3_old,mold=dx3)
    allocate(dxi1_old,mold=dxi1)
    allocate(dxi2_old,mold=dxi2)
    allocate(dxi3_old,mold=dxi3)

    ! backup arrays
    x1_old(:) = x1(:)
    x2_old(:) = x2(:)
    x3_old(:) = x3(:)
    xi1_old(:) = xi1(:)
    xi2_old(:) = xi2(:)
    xi3_old(:) = xi3(:)
    dx1_old(:) = dx1(:)
    dx2_old(:) = dx2(:)
    dx3_old(:) = dx3(:)
    dxi1_old(:) = dxi1(:)
    dxi2_old(:) = dxi2(:)
    dxi3_old(:) = dxi3(:)

    call write_grid

    ! reset values
    x1(:) = -999.d0
    x2(:) = -999.d0
    x3(:) = -999.d0
    xi1(:) = -999.d0
    xi2(:) = -999.d0
    xi3(:) = -999.d0
    dx1(:) = -999.d0
    dx2(:) = -999.d0
    dx3(:) = -999.d0
    dxi1(:) = -999.d0
    dxi2(:) = -999.d0
    dxi3(:) = -999.d0

    call readgrid('data/gridfile.bin')

    ! Check the arrays
    do i = gis_global-2, gie_global+2
      if ( .not. isequal(x1(i), x1_old(i))) then
        print*, 'Error in x1 values at i=',i,'x1(i)=',x1(i),'should be:',x1_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(xi1(i), xi1_old(i))) then
        print*, 'Error in xi1 values at i=',i,'xi1(i)=',xi1(i),'should be:',xi1_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dx1(i), dx1_old(i))) then
        print*, 'Error in dx1 values at i=',i,'dx1(i)=',dx1(i),'should be:',dx1_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dxi1(i), dxi1_old(i))) then
        print*, 'Error in dxi1 values at i=',i,'dxi1(i)=',dxi1(i),'should be:',dxi1_old(i)
        numerr = numerr + 1
      endif
    end do

    do i = gjs_global-2, gje_global+2
      if ( .not. isequal(x2(i), x2_old(i))) then
        print*, 'Error in x2 values at i=',i,'x2(i)=',x2(i),'should be:',x2_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(xi2(i), xi2_old(i))) then
        print*, 'Error in xi2 values at i=',i,'xi2(i)=',xi2(i),'should be:',xi2_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dx2(i), dx2_old(i))) then
        print*, 'Error in dx2 values at i=',i,'dx2(i)=',dx2(i),'should be:',dx2_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dxi2(i), dxi2_old(i))) then
        print*, 'Error in dxi2 values at i=',i,'dxi2(i)=',dxi2(i),'should be:',dxi2_old(i)
        numerr = numerr + 1
      endif
    end do

    do i = gks_global-2, gke_global+2
      if ( .not. isequal(x3(i), x3_old(i))) then
        print*, 'Error in x3 values at i=',i,'x3(i)=',x3(i),'should be:',x3_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(xi3(i), xi3_old(i))) then
        print*, 'Error in xi3 values at i=',i,'xi3(i)=',xi3(i),'should be:',xi3_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dx3(i), dx3_old(i))) then
        print*, 'Error in dx3 values at i=',i,'dx3(i)=',dx3(i),'should be:',dx3_old(i)
        numerr = numerr + 1
      endif
      if ( .not. isequal(dxi3(i), dxi3_old(i))) then
        print*, 'Error in dxi3 values at i=',i,'dxi3(i)=',dxi3(i),'should be:',dxi3_old(i)
        numerr = numerr + 1
      endif
    end do

  end subroutine iotest_grid

  subroutine cleanup_files(tn, time)
    use output_mod, only: set_file_name
    use mpi_utils, only: myrank
    integer, intent(in) :: tn
    real(8), intent(in) :: time
    integer :: fh, stat
    character(len=60) :: filename

    if (myrank/=0) return

    print*, 'Cleaning up files...'

    open(newunit=fh, iostat=stat, file='data/extgrv.bin', status='old')
    if (stat == 0) close(fh, status='delete')

    open(newunit=fh, iostat=stat, file='data/gridfile.bin', status='old')
    if (stat == 0) close(fh, status='delete')

    open(newunit=fh, iostat=stat, file='data/gridfile.dat', status='old')
    if (stat == 0) close(fh, status='delete')

    open(newunit=fh, iostat=stat, file='data/sinks.dat', status='old')
    if (stat == 0) close(fh, status='delete')

    call set_file_name('bin', tn, time, filename)
    open(newunit=fh, iostat=stat, file=trim(filename), status='old')
    if (stat == 0) close(fh, status='delete')

    call set_file_name('plt', tn, time, filename)
    open(newunit=fh, iostat=stat, file=trim(filename), status='old')
    if (stat == 0) close(fh, status='delete')

    open(newunit=fh, iostat=stat, file='walltime.dat', status='old')
    if (stat == 0) close(fh, status='delete')

  end subroutine cleanup_files

end module iotest_mod
