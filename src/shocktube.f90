module shocktube_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SHOCKTUBE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for shock tube problems.

subroutine shocktube

 use settings,only:mag_on,simtype,extrasfile
 use physval,only:gamma
 use input_mod,only:error_extras,error_nml

 integer::dir,ui,strl,istat
 real(8):: dl,dr,pl,pr,v1l,v1r,v2l,v2r,v3l,v3r,b1l,b1r,b2l,b2r,b3l,b3r
!-----------------------------------------------------------------------------

 namelist /tubecon/ mag_on,gamma,&
                    dl,dr,pl,pr,v1l,v1r,v2l,v2r,v3l,v3r,&
                                b1l,b1r,b2l,b2r,b3l,b3r

! Set default values
 open(newunit=ui,file='../para/extras_shocktube',status='old')
 read(ui,NML=tubecon)
 close(ui)

 strl = len(trim(simtype))

! Find the direction of shock tube
 select case(simtype(strl:strl))
 case('x')
  dir = 1
 case('y')
  dir = 2
 case('z')
  dir = 3
 end select

! Select the type of shock tube problem
 select case(simtype(1:strl-2))
 case('sodshock')
  dl = 1d0 ; dr = 0.125d0
  pl = 1d0 ; pr = 0.1d0

 case('briowushock')
  dl  = 1d0    ; dr  = 0.125d0
  pl  = 1d0    ; pr  = 0.1d0
  select case(dir)
  case(1)
   b1l = 0.75d0 ; b1r = 0.75d0
   b2l = 1d0    ; b2r =-1d0
   b3l = 0d0    ; b3r = 0d0
  case(2)
   b1l = 0d0    ; b1r = 0d0
   b2l = 0.75d0 ; b2r = 0.75d0
   b3l = 1d0    ; b3r =-1d0
  case(3)
   b1l = 1d0    ; b1r =-1d0
   b2l = 0d0    ; b2r = 0d0
   b3l = 0.75d0 ; b3r = 0.75d0
  end select

 case('other_shocktube')
  open(newunit=ui,file=extrasfile,status='old',iostat=istat)
  if(istat/=0)call error_extras('shocktube',extrasfile)
  read(ui,NML=tubecon,iostat=istat)
  if(istat/=0)call error_nml('shocktube',extrasfile)
  close(ui)

 end select


 call setup_shocktube(dir,dl,pl,v1l,v2l,v3l,b1l,b2l,b3l,&
                          dr,pr,v1r,v2r,v3r,b1r,b2r,b3r)

 return
end subroutine shocktube


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SETUP_SHOCKTUBE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for shock tube problems.

subroutine setup_shocktube(whichaxis, &
                           dl,pl,v1l,v2l,v3l,b1l,b2l,b3l,&
                           dr,pr,v1r,v2r,v3r,b1r,b2r,b3r)

  use grid
  use physval

  implicit none

  integer,intent(in):: whichaxis
  real(8),intent(in):: dl,dr,pl,pr,v1l,v1r,v2l,v2r,v3l,v3r,&
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
end subroutine setup_shocktube

end module shocktube_mod
