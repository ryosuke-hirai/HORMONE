!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SHOCKTUBE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for shock tube problems.

subroutine shocktube(whichaxis, &
                     dl,pl,v1l,v2l,v3l,b1l,b2l,b3l,&
                     dr,pr,v1r,v2r,v3r,b1r,b2r,b3r)

  use grid
  use physval

  implicit none

  integer,intent(in):: whichaxis
  real*8,intent(in):: dl,dr,pl,pr,v1l,v1r,v2l,v2r,v3l,v3r,&
                                  b1l,b1r,b2l,b2r,b3l,b3r

!--------------------------------------------------------------------

  select case(whichaxis)
  case(1)! for shock tube x direction
! left state
   d (is:(is+ie)/2  ,js:je,ks:ke) = dl
   p (is:(is+ie)/2  ,js:je,ks:ke) = pl
   v1(is:(is+ie)/2  ,js:je,ks:ke) = v1l
   v2(is:(is+ie)/2  ,js:je,ks:ke) = v2l
   v3(is:(is+ie)/2  ,js:je,ks:ke) = v3l
   b1(is:(is+ie)/2  ,js:je,ks:ke) = b1l
   b2(is:(is+ie)/2  ,js:je,ks:ke) = b2l
   b3(is:(is+ie)/2  ,js:je,ks:ke) = b3l
! right state
   d ((is+ie)/2+1:ie,js:je,ks:ke) = dr
   p ((is+ie)/2+1:ie,js:je,ks:ke) = pr
   v1((is+ie)/2+1:ie,js:je,ks:ke) = v1r
   v2((is+ie)/2+1:ie,js:je,ks:ke) = v2r
   v3((is+ie)/2+1:ie,js:je,ks:ke) = v3r
   b1((is+ie)/2+1:ie,js:je,ks:ke) = b1r
   b2((is+ie)/2+1:ie,js:je,ks:ke) = b2r
   b3((is+ie)/2+1:ie,js:je,ks:ke) = b3r

  case(2)! for shock tube y direction
! left state
   d (is:ie,js:(js+je)/2  ,ks:ke) = dl
   p (is:ie,js:(js+je)/2  ,ks:ke) = pl
   v1(is:ie,js:(js+je)/2  ,ks:ke) = v1l
   v2(is:ie,js:(js+je)/2  ,ks:ke) = v2l
   v3(is:ie,js:(js+je)/2  ,ks:ke) = v3l
   b1(is:ie,js:(js+je)/2  ,ks:ke) = b1l
   b2(is:ie,js:(js+je)/2  ,ks:ke) = b2l
   b3(is:ie,js:(js+je)/2  ,ks:ke) = b3l
! right state
   d (is:ie,(js+je)/2+1:je,ks:ke) = dr
   p (is:ie,(js+je)/2+1:je,ks:ke) = pr
   v1(is:ie,(js+je)/2+1:je,ks:ke) = v1r
   v2(is:ie,(js+je)/2+1:je,ks:ke) = v2r
   v3(is:ie,(js+je)/2+1:je,ks:ke) = v3r
   b1(is:ie,(js+je)/2+1:je,ks:ke) = b1r
   b2(is:ie,(js+je)/2+1:je,ks:ke) = b2r
   b3(is:ie,(js+je)/2+1:je,ks:ke) = b3r

  case(3)! for shock tube z direction
! left state
   d (is:ie,js:je,ks:(ks+ke)/2  ) = dl
   p (is:ie,js:je,ks:(ks+ke)/2  ) = pl
   v1(is:ie,js:je,ks:(ks+ke)/2  ) = v1l
   v2(is:ie,js:je,ks:(ks+ke)/2  ) = v2l
   v3(is:ie,js:je,ks:(ks+ke)/2  ) = v3l
   b1(is:ie,js:je,ks:(ks+ke)/2  ) = b1l
   b2(is:ie,js:je,ks:(ks+ke)/2  ) = b2l
   b3(is:ie,js:je,ks:(ks+ke)/2  ) = b3l
! right state
   d (is:ie,js:je,(ks+ke)/2+1:ke) = dr
   p (is:ie,js:je,(ks+ke)/2+1:ke) = pr
   v1(is:ie,js:je,(ks+ke)/2+1:ke) = v1r
   v2(is:ie,js:je,(ks+ke)/2+1:ke) = v2r
   v3(is:ie,js:je,(ks+ke)/2+1:ke) = v3r
   b1(is:ie,js:je,(ks+ke)/2+1:ke) = b1r
   b2(is:ie,js:je,(ks+ke)/2+1:ke) = b2r
   b3(is:ie,js:je,(ks+ke)/2+1:ke) = b3r

  case default
   print*,'Provide an index between 1 and 3 for shock tube direction'
   print*,'whichaxis = "',whichaxis,'"'
   stop
  end select

return
end subroutine shocktube
