module hllflux_mod

 implicit none
 
 public:: hlldflux,hllflux
 private:: fstar
 
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE HLLDFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte hlldflux
!          See Miyoshi & Kusano 2005 for details

 subroutine hlldflux(fflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr,el,er,&
                     ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

  use physval,only:ch

  implicit none

  real(8),intent(in):: cfl, cfr, v1l, v1r, v2l, v2r, v3l, v3r,&
                      dl, dr, ptl, ptr, el, er, b1l, b1r, b2l, b2r, b3l, b3r
  real(8):: sl, sr, sm, sla, sra
  real(8):: dla, dra, v2la, v2ra, v3la, v3ra, b2la, b2ra, b3la, b3ra
  real(8):: ela, era, pta
  real(8):: v2aa, v3aa, b2aa, b3aa, elaa, eraa
  real(8):: signsm, sortsla, sortslr, sortsra
  real(8),intent(in):: phil, phir
  real(8):: phi, b1 ! for 9-wave method
  real(8),dimension(9),intent(inout):: fflux

!----------------------------------------------------------------------------

! added for 9wave method -------------------------!
  b1  = 0.5d0 * ( (b1l +b1r ) - (phir-phil)/ch )  !
  phi = 0.5d0 * ( (phil+phir) - (b1r -b1l )*ch )  !
! end 9wave method--------------------------------!


! find left and right waves
  sl  = min(v1r-cfr,v1l-cfl,0.d0)
  sr  = max(v1r+cfr,v1l+cfl,0.d0)
! sr & sl set to 0 for supersonic states
!print *,v1l,cfl
  sm  = ( (sr-v1r)*dr*v1r - (sl-v1l)*dl*v1l - ptr + ptl ) &
        / ( (sr-v1r)*dr - (sl-v1l)*dl )
!print *,(sr-v1r)*dr*v1r, - (sl-v1l)*dl*v1l, - ptr + ptl
  dla = dl * (sl-v1l)/(sl-sm)
  dra = dr * (sr-v1r)/(sr-sm)
!print *,sr,sm,v1r,dra,dla
  sla = sm - abs(b1)/sqrt(dla)
  sra = sm + abs(b1)/sqrt(dra)

  pta = ( (sr-v1r)*dr*ptl-(sl-v1l)*dl*ptr+dl*dr*(sr-v1r)*(sl-v1l)*(v1r-v1l) )&
        / ( (sr-v1r)*dr-(sl-v1l)*dl )

  v2la = v2l - b1*b2l*(sm-v1l) &
             / ( dl*(sl-v1l)*(sl-sm) - b1*b1 )
  v2ra = v2r - b1*b2r*(sm-v1r) &
             / ( dr*(sr-v1r)*(sr-sm) - b1*b1 )

  v3la = v3l - b1*b3l*(sm-v1l) &
             / ( dl*(sl-v1l)*(sl-sm) - b1*b1 )
  v3ra = v3r - b1*b3r*(sm-v1r) &
             / ( dr*(sr-v1r)*(sr-sm) - b1*b1 )


  b2la = b2l * ( dl*(sl-v1l)*(sl-v1l) - b1*b1 ) &
             / ( dl*(sl-v1l)*(sl-sm) -b1*b1 )
  b2ra = b2r * ( dr*(sr-v1r)*(sr-v1r) - b1*b1 ) &
             / ( dr*(sr-v1r)*(sr-sm) -b1*b1 )

  b3la = b3l * ( dl*(sl-v1l)*(sl-v1l) - b1*b1 ) &
             / ( dl*(sl-v1l)*(sl-sm) -b1*b1 )
  b3ra = b3r * ( dr*(sr-v1r)*(sr-v1r) - b1*b1 ) &
             / ( dr*(sr-v1r)*(sr-sm) -b1*b1 )


  ela  = ( (sl-v1l)*el - ptl*v1l + pta*sm + b1 * &
         ( (v1l*b1+v2l*b2l+v3l*b3l) - (sm*b1+v2la*b2la+v3la*b3la)) ) &
         / (sl-sm)
  era  = ( (sr-v1r)*er - ptr*v1r + pta*sm + b1 * &
         ( (v1r*b1+v2r*b2r+v3r*b3r) - (sm*b1+v2ra*b2ra+v3ra*b3ra)) ) &
         / (sr-sm)


  v2aa = ( sqrt(dla)*v2la + sqrt(dra)*v2ra + (b2ra-b2la)*sign(1.d0,b1) ) &
         / ( sqrt(dla) + sqrt(dra) )
  v3aa = ( sqrt(dla)*v3la + sqrt(dra)*v3ra + (b3ra-b3la)*sign(1.d0,b1) ) &
         / ( sqrt(dla) + sqrt(dra) )

  b2aa = ( sqrt(dla)*b2ra+sqrt(dra)*b2la + &
           sqrt(dla*dra)*(v2ra-v2la)*sign(1.d0,b1) ) &
         / ( sqrt(dla) + sqrt(dra) )
  b3aa = ( sqrt(dla)*b3ra+sqrt(dra)*b3la + &
           sqrt(dla*dra)*(v3ra-v3la)*sign(1.d0,b1) ) &
         / ( sqrt(dla) + sqrt(dra) )

  eraa = era + sign(1.d0,b1)*sqrt(dra)* &
         ((v2ra*b2ra+v3ra*b3ra)-(v2aa*b2aa+v3aa*b3aa))
  elaa = ela - sign(1.d0,b1)*sqrt(dla)* &
         ((v2la*b2la+v3la*b3la)-(v2aa*b2aa+v3aa*b3aa))

! start sorting

! set frequently used variables --------!
  signsm  = sign(0.5d0,sm)             !
  sortsla = 0.5d0+sign(0.5d0, sla)     !
  sortsra = 0.5d0+sign(0.5d0,-sra)     !
  sortslr = 0.5d0-sign(0.5d0,sla*sra)  !
!---------------------------------------!
!print *,sortsla,sortsra,sortslr
  dla  =   (0.5d0+signsm) * dla  &
         + (0.5d0-signsm) * dra

  v2la =   sortsla * v2la &
         + sortslr * v2aa &
         + sortsra * v2ra 

  v3la =   sortsla * v3la &
         + sortslr * v3aa &
         + sortsra * v3ra 

  b2la =   sortsla * b2la &
         + sortslr * b2aa &
         + sortsra * b2ra 

  b3la =   sortsla * b3la &
         + sortslr * b3aa &
         + sortsra * b3ra 

  ela  =   sortsla * ela &
         + (0.5d0-sign(0.5d0, sla*sm)) * elaa &
         + (0.5d0+sign(0.5d0,-sm*sra)) * eraa &
         + (0.5d0-sign(0.5d0,sra)) * era

  call fstar(fflux,dla,sm,v2la,v3la,b1,b2la,b3la,ela,pta,phi)

  return
 end subroutine hlldflux



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 subroutine fstar(tmpflux,d,v1,v2,v3,b1,b2,b3,e,p,phi)

  use physval,only:ch

  implicit none

  real(8),dimension(9),intent(out)::tmpflux
  real(8),intent(in)::d,v1,v2,v3,b1,b2,b3,e,p,phi

!----------------------------------------------------------------------------

  tmpflux(1) = d*v1
  tmpflux(2) = d*v1*v1 + p - b1*b1 
  tmpflux(3) = d*v1*v2 - b1*b2 
  tmpflux(4) = d*v1*v3 - b1*b3
  tmpflux(5) = phi
  tmpflux(6) = b2*v1 - b1*v2 
  tmpflux(7) = b3*v1 - b1*v3 
  tmpflux(8) = (e+p)*v1 - b1*(v1*b1+v2*b2+v3*b3 )
! for 9wave method
  tmpflux(9) = ch*ch * b1
! end 9wave method
  return
 end subroutine fstar

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE HLLFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte hllflux

 pure subroutine hllflux(fflux,fl,fr,ul,ur,cfl,cfr,vl,vr)

  implicit none

  real(8),intent(in) :: fl, fr, ul, ur, cfl, cfr, vl, vr
  real(8):: sl, sr
  real(8),intent(inout)::fflux

!--------------------------------------------------------------------

! find left and right waves
  sr  = max(vr+cfr,vl+cfl,0.d0)
  sl  = min(vr-cfr,vl-cfl,0.d0)

  fflux = (sr*fl - sl*fr + sr*sl*(ur-ul)) / (sr-sl)

  return
 end subroutine hllflux

end module hllflux_mod
