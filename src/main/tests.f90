module tests_mod
 implicit none

 public:: check_testlist,test
 private:: testfilename,print_errors,ratio,&
           max_norm_error,L1_norm_error,L2_norm_error

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE CHECK_TESTLIST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To check whether the test data is available

 subroutine check_testlist(simtype)

  character(len=*),intent(in)::simtype
  logical:: test_available
  character(40):: filename

!-----------------------------------------------------------------------------

  filename = testfilename(simtype)
  inquire(file=filename,exist=test_available)

  if(.not.test_available)then
   print*,'Test data is not available for the following simtype'
   print'(3a)',' simtype = "',trim(simtype),'"'
   stop
  end if

  return
 end subroutine check_testlist

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE TEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compare results against pre-computed data

 subroutine test(passed)

  use utils,only:isequal
  use mpi_utils,only:myrank
  use settings,only:simtype,mag_on,test_tol,Mach_tol,compswitch,spn,radswitch
  use grid,only:is,ie,js,je,ks,ke,dim
  use physval
  use gravmod
  use readbin_mod,only:readbin
  use timestep_mod,only:timestep

  logical, intent(out) :: passed
  integer:: n
  integer,parameter:: nn = 11
  character(40):: testfile
  real(8),allocatable:: error(:)
  real(8),allocatable,dimension(:,:,:,:):: val,valorg,scale
  real(8),allocatable,dimension(:,:,:):: vmag,bmag
  character(len=10),allocatable:: label(:)

!-----------------------------------------------------------------------------

  allocate(val(is:ie,js:je,ks:ke,nn+spn),label(1:nn+spn),error(1:nn+spn))
  allocate(valorg,scale,mold=val)
  allocate(vmag(is:ie,js:je,ks:ke),bmag(is:ie,js:je,ks:ke))
  error = 0d0

  label( 1) = 'density'
  label( 2) = 'energy'
  label( 3) = 'v1'
  label( 4) = 'v2'
  label( 5) = 'v3'
  label( 6) = 'b1'
  label( 7) = 'b2'
  label( 8) = 'b3'
  label( 9) = 'divB'
  label(10) = 'gravity'
  label(11) = 'erad'
  if(compswitch>=2)then ! Chemical elements
   do n = 1, spn
    label(nn+n) = trim(species(n))
   end do
  end if
  if(.not.mag_on)       label(6:9) = 'aaa'   ! No magnetic field
  if(mag_on.and.dim==1) label(9) = 'aaa'     ! 1D MHD
  if(gravswitch==0)     label(10) = 'aaa'    ! No gravity
  if(radswitch==0)       label(11) = 'aaa'    ! No radiation

! First record simulated variables
  val(:,:,:,1) = d (is:ie,js:je,ks:ke)
  val(:,:,:,2) = e (is:ie,js:je,ks:ke)
  val(:,:,:,3) = v1(is:ie,js:je,ks:ke)
  val(:,:,:,4) = v2(is:ie,js:je,ks:ke)
  val(:,:,:,5) = v3(is:ie,js:je,ks:ke)
  val(:,:,:,6) = b1(is:ie,js:je,ks:ke)
  val(:,:,:,7) = b2(is:ie,js:je,ks:ke)
  val(:,:,:,8) = b3(is:ie,js:je,ks:ke)
  if(mag_on.and.dim>=2) val(:,:,:, 9) = phi   (is:ie,js:je,ks:ke)
  if(gravswitch>0)      val(:,:,:,10) = grvphi(is:ie,js:je,ks:ke)
  if(radswitch>0)       val(:,:,:,11) = erad  (is:ie,js:je,ks:ke)
  if(compswitch>=2)then
   do n = 1, spn
    val(:,:,:,10+n) = spc(n,is:ie,js:je,ks:ke)
   end do
  end if

! Open test data file
  testfile = testfilename(simtype)

  call readbin(testfile)

! Next record variables from pre-computed file
  valorg(:,:,:,1) = d (is:ie,js:je,ks:ke)
  valorg(:,:,:,2) = e (is:ie,js:je,ks:ke)
  valorg(:,:,:,3) = v1(is:ie,js:je,ks:ke)
  valorg(:,:,:,4) = v2(is:ie,js:je,ks:ke)
  valorg(:,:,:,5) = v3(is:ie,js:je,ks:ke)
  valorg(:,:,:,6) = b1(is:ie,js:je,ks:ke)
  valorg(:,:,:,7) = b2(is:ie,js:je,ks:ke)
  valorg(:,:,:,8) = b3(is:ie,js:je,ks:ke)
  if(mag_on.and.dim>=2) valorg(:,:,:, 9) = phi   (is:ie,js:je,ks:ke)
  if(gravswitch>0)      valorg(:,:,:,10) = grvphi(is:ie,js:je,ks:ke)
  if(radswitch>0)       valorg(:,:,:,11) = erad  (is:ie,js:je,ks:ke)
  if(compswitch>=2)then
   do n = 1, spn
    valorg(:,:,:,10+n) = spc(n,is:ie,js:je,ks:ke)
   end do
  end if

  ! Calculate the magnitude of velocity and magnetic field
  vmag(is:ie,js:je,ks:ke) = sqrt(v1(is:ie,js:je,ks:ke)**2+v2(is:ie,js:je,ks:ke)**2+v3(is:ie,js:je,ks:ke)**2)
  bmag(is:ie,js:je,ks:ke) = sqrt(b1(is:ie,js:je,ks:ke)**2+b2(is:ie,js:je,ks:ke)**2+b3(is:ie,js:je,ks:ke)**2)
  ! Calculate the sound speed
  call timestep

  ! Magnitude of this quantity used to normalise the error
  do n = 1, nn+spn
    ! If the variable is a component of velocity or magnetic field, use the magnitude
    if (n>=3.and.n<=5) then
      ! Down-weigh velocity errors when Mach<<1
      scale(is:ie,js:je,ks:ke,n) = max(vmag(is:ie,js:je,ks:ke), &
                                       cs  (is:ie,js:je,ks:ke) * Mach_tol )
    else if (n>=6.and.n<=8) then
      scale(is:ie,js:je,ks:ke,n) = bmag(is:ie,js:je,ks:ke)
    else
      scale(is:ie,js:je,ks:ke,n) = abs(valorg(is:ie,js:je,ks:ke,n))
    end if

    ! If the value is smaller than epsilon, set the scale to epsilon
    scale(is:ie,js:je,ks:ke,n) = max(scale(is:ie,js:je,ks:ke,n),epsilon(valorg(is:ie,js:je,ks:ke,n)))

  end do

! Calculate max norm errors
  call print_errors('max',max_norm_error,label,val,valorg,scale,error)

! Calculate L1 norm errors
  call print_errors('L1',L1_norm_error,label,val,valorg,scale,error)

! Calculate L2 norm errors
  call print_errors('L2',L2_norm_error,label,val,valorg,scale,error)

! Check if maximum L2 norm error is within acceptable bounds
  if (myrank==0) print*, 'Test tolerance (L2) =', test_tol
  if(maxval(error)<test_tol)then
   if (myrank==0) then
      print*,trim(simtype),' test: passed'
      if(maxval(error)<=0d0) print*,trim(simtype),'     : Identical!'
   endif
   passed = .true.
  else
   if (myrank==0) print*,trim(simtype),' test: failed'
   passed = .false.
  end if

  return
 end subroutine test

 function testfilename(simtype) result(file)
  use settings,only:flux_limiter
  character(len=*),intent(in)::simtype
  character(40)::file

  if(flux_limiter=='flat')then
   ! If using flat reconstruction
   file = '../tests/'//trim(simtype)//'_flat.bin'
  else
   file = '../tests/'//trim(simtype)//'.bin'
  end if

 end function testfilename

 function max_norm_error(var,var0,scale) result(norm)
  use mpi_utils,only:allreduce_mpi
  real(8),dimension(:,:,:),intent(in):: var,var0,scale
  real(8):: norm, pos

  pos = maxval(var)*minval(var)

  if(pos>0d0)then ! for strictly positive quantities
   norm = maxval(abs(var/var0-1d0))
  else ! for quantities that can contain zeroes
   norm = maxval(abs(var-var0)/max(abs(var0),scale))
  end if

  ! Reduce error across MPI tasks
  call allreduce_mpi('max', norm)

 end function max_norm_error

 function L1_norm_error(var,var0,scale) result(norm)
  use mpi_utils,only:allreduce_mpi
  use grid,only:is,ie,js,je,ks,ke,dvol
  real(8),dimension(:,:,:),intent(in):: var,var0,scale
  real(8):: norm, pos, errsum, wsum
  real(8),allocatable:: w(:,:,:)

  pos = maxval(var)*minval(var)

! Define weight function = volume
  allocate(w,mold=var)
  w(:,:,:) = dvol(is:ie,js:je,ks:ke)

  if(pos>0d0)then ! for strictly positive quantities
   errsum = sum(abs(var/var0-1d0)*w)
   wsum = sum(w)
  else ! for quantities that can contain zeroes
   errsum = sum(abs(var-var0)/max(abs(var0),scale)*w)
   wsum = sum(w)
  end if

  ! Reduce error across MPI tasks
  call allreduce_mpi('sum', errsum)
  call allreduce_mpi('sum', wsum)
  norm = sqrt(errsum)/sqrt(wsum)

 end function L1_norm_error

 function L2_norm_error(var,var0,scale) result(norm)
  use mpi_utils,only:allreduce_mpi
  use grid,only:is,ie,js,je,ks,ke,dvol
  real(8),dimension(:,:,:),intent(in):: var,var0,scale
  real(8):: norm, pos, jump, base, denom, floor, errsum, wsum
  real(8),allocatable:: relerr(:,:,:),w(:,:,:)
  integer:: i,j,k,il,jl,kl,iu,ju,ku,disco_range

  pos = maxval(var)*minval(var)

  ! var and var0, when passed through as regular array arguments,
  ! have the same shape as the original arrays, but indices starting from 1,
  ! not is, js, ks. When relerr and w are allocated with mold=var, they will
  ! also have indices starting from 1.

  allocate(relerr,w,mold=var)
  w(:,:,:) = dvol(is:ie,js:je,ks:ke)

! Weigh down the error if there is a discontinuity within this number of cells
  disco_range = 3
  do k = lbound(var,3), ubound(var,3)
   kl=min(disco_range,k-lbound(var,3));ku=min(disco_range,ubound(var,3)-k)
   do j = lbound(var,2), ubound(var,2)
    jl=min(disco_range,j-lbound(var,2));ju=min(disco_range,ubound(var,2)-j)
    do i = lbound(var,1), ubound(var,1)
     il=min(disco_range,i-lbound(var,1));iu=min(disco_range,ubound(var,1)-i)
     base = var0(i,j,k)
     denom = maxval(abs(var0(i-il:i+iu,j-jl:j+ju,k-kl:k+ku)))
     floor = 0d0
     if(pos<=0d0)then ! Use floor value if variable is allowed to be zero
      base = max(abs(base),scale(i,j,k))
      denom = max(denom,scale(i,j,k))
      floor = scale(i,j,k)
     end if
     jump = maxval(ratio(var0(i-il:i+iu,j-jl:j+ju,k-kl:k+ku),base,floor))

     relerr(i,j,k) = (var(i,j,k)-var0(i,j,k))/denom*exp(1d0-jump)
    end do
   end do
  end do

  ! Reduce error across MPI tasks
  errsum = sum(relerr**2*w)
  wsum = sum(w)
  call allreduce_mpi('sum', errsum)
  call allreduce_mpi('sum', wsum)
  norm = sqrt(errsum)/sqrt(wsum)

 end function L2_norm_error

 pure elemental function ratio(x,y,floor)
  real(8),intent(in)::x,y,floor
  real(8):: ratio
  ratio = max(abs(x)/max(abs(y),floor),abs(y)/max(abs(x),floor),1d0)
 end function ratio

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE PRINT_ERRORS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To print errors defined by a given norm

 subroutine print_errors(name,f,label,val,valorg,scale,error)
  use mpi_utils,only:myrank

  interface
   real(8) function f(var,var0,scale_n)
    real(8),dimension(:,:,:),intent(in):: var,var0,scale_n
    real(8):: norm
   end function f
  end interface
  character(len=*),intent(in):: name
  character(len=10),intent(in):: label(:)
  real(8),allocatable,intent(inout),dimension(:,:,:,:):: val,valorg,scale
  real(8),intent(inout)::error(:)
  integer:: n

!-----------------------------------------------------------------------------

  if (myrank == 0) print*,trim(name),' norm errors:'
  do n = 1, size(error)
   if(trim(label(n))=='aaa')cycle
   error(n) = f(val(:,:,:,n),valorg(:,:,:,n),scale(:,:,:,n))
   if (myrank == 0) print*,'  ',label(n),' =',error(n)
  end do

  return
 end subroutine print_errors

end module tests_mod
