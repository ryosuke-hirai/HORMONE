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

  use settings,only:simtype
  use grid,only:is,ie,js,je,ks,ke
  use physval
  use readbin_mod,only:readbin

  character(30):: testfile
  real*8,parameter:: tolerance=1d-10
  real*8:: derr,perr,v1err,v2err,v3err
  
!-----------------------------------------------------------------------------

! First record variables
  d0 (is:ie,js:je,ks:ke) = d (is:ie,js:je,ks:ke) 
  p0 (is:ie,js:je,ks:ke) = p (is:ie,js:je,ks:ke) 
  v10(is:ie,js:je,ks:ke) = v1(is:ie,js:je,ks:ke) 
  v20(is:ie,js:je,ks:ke) = v2(is:ie,js:je,ks:ke) 
  v30(is:ie,js:je,ks:ke) = v3(is:ie,js:je,ks:ke) 
  b10(is:ie,js:je,ks:ke) = b1(is:ie,js:je,ks:ke) 
  b20(is:ie,js:je,ks:ke) = b2(is:ie,js:je,ks:ke) 
  b30(is:ie,js:je,ks:ke) = b3(is:ie,js:je,ks:ke) 

! Open test data file
  testfile = testfilename(simtype)
  
  call readbin(testfile)

  derr = maxval(abs(d0(is:ie,js:je,ks:ke)-d(is:ie,js:je,ks:ke)) &
                / d(is:ie,js:je,ks:ke))
  perr = maxval(abs(p0(is:ie,js:je,ks:ke)-p(is:ie,js:je,ks:ke)) &
                / p(is:ie,js:je,ks:ke))
  v1err = maxval(abs(v10(is:ie,js:je,ks:ke)-v1(is:ie,js:je,ks:ke)) &
                 / max(abs(v1(is:ie,js:je,ks:ke)),tolerance))
  v2err = maxval(abs(v20(is:ie,js:je,ks:ke)-v2(is:ie,js:je,ks:ke)) &
                 / max(abs(v1(is:ie,js:je,ks:ke)),tolerance))
  v3err = maxval(abs(v30(is:ie,js:je,ks:ke)-v3(is:ie,js:je,ks:ke)) &
                 / max(abs(v1(is:ie,js:je,ks:ke)),tolerance))
  
  print*,'Maximum error is =',max(derr,perr,v1err,v2err,v3err)
  if(max(derr,perr,v1err,v2err,v3err)<tolerance)then
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

end module tests_mod
