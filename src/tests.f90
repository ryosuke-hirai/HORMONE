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
  character(30):: filename

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

  use settings,only:simtype,mag_on
  use grid,only:is,ie,js,je,ks,ke
  use physval
  use gravmod
  use readbin_mod,only:readbin

  logical, intent(out) :: passed
  integer,parameter:: nn = 9
  character(30):: testfile
  real(8),parameter:: tol=1d-4
  real(8):: error(nn)
  real(8),allocatable,dimension(:,:,:,:):: val,valorg
  character(len=10):: label(nn)

!-----------------------------------------------------------------------------

  allocate(val(is:ie,js:je,ks:ke,nn))
  allocate(valorg,mold=val)
  error = 0d0

  label(1) = 'density'
  label(2) = 'energy'
  label(3) = 'v1'
  label(4) = 'v2'
  label(5) = 'v3'
  label(6) = 'b1'
  label(7) = 'b2'
  label(8) = 'b3'
  label(9) = 'gravity'

  if(.not.mag_on) label(6:8) = 'aaa'

! First record simulated variables
  val(:,:,:,1) = d (is:ie,js:je,ks:ke)
  val(:,:,:,2) = e (is:ie,js:je,ks:ke)
  val(:,:,:,3) = v1(is:ie,js:je,ks:ke)
  val(:,:,:,4) = v2(is:ie,js:je,ks:ke)
  val(:,:,:,5) = v3(is:ie,js:je,ks:ke)
  val(:,:,:,6) = b1(is:ie,js:je,ks:ke)
  val(:,:,:,7) = b2(is:ie,js:je,ks:ke)
  val(:,:,:,8) = b3(is:ie,js:je,ks:ke)
  if(gravswitch>0)then
   val(:,:,:,9) = grvphi(is:ie,js:je,ks:ke)
  else
   label(9) = 'aaa'
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
  if(gravswitch>0)then
   valorg(:,:,:,9) = grvphi(is:ie,js:je,ks:ke)
  end if

! Calculate max norm errors
  call print_errors('max',max_norm_error,label,val,valorg,tol,error)

! Calculate L1 norm errors
  call print_errors('L1',L1_norm_error,label,val,valorg,tol,error)

! Calculate L2 norm errors
  call print_errors('L2',L2_norm_error,label,val,valorg,tol,error)

! Check if maximum L2 norm error is within acceptable bounds
  if(maxval(error)<tol)then
   print*,trim(simtype),' test: passed'
   if(maxval(error)<=0d0)print*,trim(simtype),'     : Identical!'
   passed = .true.
  else
   print*,trim(simtype),' test: failed'
   passed = .false.
  end if

  return
 end subroutine test

 function testfilename(simtype) result(file)
  character(len=*),intent(in)::simtype
  character(30)::file
  file = '../tests/'//trim(simtype)//'.bin'
 end function testfilename

 function max_norm_error(var,var0,tol) result(norm)
  real(8),dimension(:,:,:),intent(in):: var,var0
  real(8),intent(in):: tol
  real(8):: norm, pos

  pos = maxval(var)*minval(var)

  if(pos>0d0)then ! for strictly positive quantities
   norm = maxval(abs(var/var0-1d0))
  else ! for quantities that can contain zeroes
   norm = maxval(abs(var-var0)/max(abs(var0),tol))
  end if

 end function max_norm_error

 function L1_norm_error(var,var0,tol) result(norm)
  use grid,only:is,ie,js,je,ks,ke,dvol
  real(8),dimension(:,:,:),intent(in):: var,var0
  real(8),intent(in):: tol
  real(8):: norm, pos
  real(8),allocatable:: w(:,:,:)

  pos = maxval(var)*minval(var)

! Define weight function = volume
  allocate(w,mold=var)
  w(:,:,:) = dvol(is:ie,js:je,ks:ke)

  if(pos>0d0)then ! for strictly positive quantities
   norm = sum(abs(var/var0-1d0)*w)/sum(w)
  else ! for quantities that can contain zeroes
   norm = sum(abs(var-var0)/max(abs(var0),tol)*w)/sum(w)
  end if

 end function L1_norm_error

 function L2_norm_error(var,var0,tol) result(norm)
  use grid,only:is,ie,js,je,ks,ke,dvol
  real(8),dimension(:,:,:),intent(in):: var,var0
  real(8),intent(in):: tol
  real(8):: norm, pos, jump, base, denom, floor
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
      base = max(abs(base),tol)
      denom = max(denom,tol)
      floor = tol
     end if
     jump = maxval(ratio(var0(i-il:i+iu,j-jl:j+ju,k-kl:k+ku),base,floor))

     relerr(i,j,k) = (var(i,j,k)-var0(i,j,k))/denom*exp(1d0-jump)
    end do
   end do
  end do

  norm = sqrt(sum(relerr**2*w))/sqrt(sum(w))

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

 subroutine print_errors(name,f,label,val,valorg,tol,error)

  interface
   real(8) function f(var,var0,tol)
    real(8),dimension(:,:,:),intent(in):: var,var0
    real(8),intent(in):: tol
    real(8):: norm
   end function f
  end interface
  character(len=*),intent(in):: name
  character(len=10),intent(in):: label(:)
  real(8),allocatable,intent(inout),dimension(:,:,:,:):: val,valorg
  real(8),intent(in):: tol
  real(8),intent(inout)::error(:)
  integer:: n

!-----------------------------------------------------------------------------

  print*,trim(name),' norm errors:'
  do n = 1, size(error)
   if(trim(label(n))=='aaa')cycle
   error(n) = f(val(:,:,:,n),valorg(:,:,:,n),tol)
   print*,'  ',label(n),' =',error(n)
  end do

  return
 end subroutine print_errors

end module tests_mod
