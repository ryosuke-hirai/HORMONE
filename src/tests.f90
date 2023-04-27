module tests_mod
 implicit none

 public:: check_testlist

 integer,parameter:: ntest=8
! List of available test data
 character(30),dimension(ntest),parameter:: &
  testlist = (/ &
  'sedov_default',&
  'briowushock_x',&
  'briowushock_y',&
  'briowushock_z',&
  'sodshock_x',&
  'sodshock_y',&
  'sodshock_z',&
  'orszagtang_xy' &
  /)

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
  integer:: i

!-----------------------------------------------------------------------------

  test_available = .false.
  do i = 1, ntest
   if(trim(simtype)==trim(testlist(i)))then
    test_available = .true.
    exit
   end if
  end do

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
  integer:: ui,istat
  real*8,parameter:: tolerance=1d-10
  real*8:: derr,perr,v1err,v2err,v3err,berr
  
!-----------------------------------------------------------------------------

! First record variables
  d0  = d
  p0  = p
  v10 = v1
  v20 = v2
  v30 = v3
  b10 = b1
  b20 = b2
  b30 = b3

! Open test data file
  testfile = '../tests/'//trim(simtype)//'.bin'
  
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

end module tests_mod
