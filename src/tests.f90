module tests_mod
 implicit none

 public:: check_testlist,test
 private:: testfilename

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

 subroutine test

  use settings,only:simtype,mag_on
  use grid,only:is,ie,js,je,ks,ke
  use physval
  use gravmod
  use readbin_mod,only:readbin

  integer,parameter:: nn = 9
  integer:: n
  character(30):: testfile
  real(8),parameter:: tol=1d-10
  real(8):: error(nn)
  real(8),allocatable,dimension(:,:,:,:):: val,valorg
  character(len=10):: label(nn)
  
!-----------------------------------------------------------------------------

  allocate(val(nn,is:ie,js:je,ks:ke))
  allocate(valorg,mold=val)

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

! First record variables
  val(1,:,:,:) = d (is:ie,js:je,ks:ke) 
  val(2,:,:,:) = e (is:ie,js:je,ks:ke) 
  val(3,:,:,:) = v1(is:ie,js:je,ks:ke) 
  val(4,:,:,:) = v2(is:ie,js:je,ks:ke) 
  val(5,:,:,:) = v3(is:ie,js:je,ks:ke) 
  val(6,:,:,:) = b1(is:ie,js:je,ks:ke) 
  val(7,:,:,:) = b2(is:ie,js:je,ks:ke) 
  val(8,:,:,:) = b3(is:ie,js:je,ks:ke)
  if(gravswitch>0)then
   val(9,:,:,:) = grvphi(is:ie,js:je,ks:ke)
  else
   label(9) = 'aaa'
  end if

! Open test data file
  testfile = testfilename(simtype)
  
  call readbin(testfile)

  valorg(1,:,:,:) = d (is:ie,js:je,ks:ke) 
  valorg(2,:,:,:) = e (is:ie,js:je,ks:ke) 
  valorg(3,:,:,:) = v1(is:ie,js:je,ks:ke) 
  valorg(4,:,:,:) = v2(is:ie,js:je,ks:ke) 
  valorg(5,:,:,:) = v3(is:ie,js:je,ks:ke) 
  valorg(6,:,:,:) = b1(is:ie,js:je,ks:ke) 
  valorg(7,:,:,:) = b2(is:ie,js:je,ks:ke) 
  valorg(8,:,:,:) = b3(is:ie,js:je,ks:ke)
  if(gravswitch>0)then
   valorg(9,:,:,:) = grvphi(is:ie,js:je,ks:ke)
  end if

! Calculate L1 norm errors
  print*,'L1 norm errors:'
  do n = 1, nn
   error(n) = L1_norm_error(val(n,:,:,:),valorg(n,:,:,:),tol)
   if(trim(label(n))/='aaa')print*,label(n),' =',error(n)
  end do

! Calculate L2 norm errors
  print*,'L2 norm errors:'
  do n = 1, nn
   error(n) = L2_norm_error(val(n,:,:,:),valorg(n,:,:,:),tol)
   if(trim(label(n))/='aaa')print*,label(n),' =',error(n)
  end do

!  if(max(derr,eerr,v1err,v2err,v3err,gerr)<tol)then
  if(maxval(error)<tol)then
   print*,trim(simtype),' test: passed'
  else
   print*,trim(simtype),' test: failed'
  end if
  
  return
 end subroutine test

 function testfilename(simtype) result(file)
  character(len=*),intent(in)::simtype
  character(30)::file
  file = '../tests/'//trim(simtype)//'.bin'
 end function testfilename

 function L1_norm_error(var,var0,tol) result(norm)
  real(8),dimension(:,:,:),intent(in):: var,var0
  real(8),intent(in):: tol
  real(8):: norm, pos

  pos = maxval(var)*minval(var)

  if(pos>0d0)then ! for strictly positive quantities
   norm = maxval(abs(var/var0-1d0))
  else ! for quantities that can contain zeroes
   norm = maxval(abs(var-var0)/max(abs(var0),tol))
  end if

 end function L1_norm_error

 function L2_norm_error(var,var0,tol) result(norm)
  real(8),dimension(:,:,:),intent(in):: var,var0
  real(8),intent(in):: tol
  real(8):: norm, pos

  pos = maxval(var)*minval(var)

  if(pos>0d0)then ! for strictly non-zero quantities
   norm = norm2(abs(var/var0-1d0))/sqrt(dble(size(var)))
  else ! for quantities that can contain zeroes
   norm = norm2(max(abs(var),tol)/max(abs(var0),tol)-1d0)/sqrt(dble(size(var)))
  end if

 end function L2_norm_error

end module tests_mod
