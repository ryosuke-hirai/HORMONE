!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE SEDOV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set initial conditions for Sedov-Taylor blast wave test

subroutine sedov(damb,Eexp)

 use grid
 use physval
 use pressure_mod,only:eos_p

 implicit none

 real*8,intent(in):: damb,Eexp
 real*8:: ein, pin, pamb, Tin, imuconst
 integer:: i_inj

!-----------------------------------------------------------------------------

 i_inj = 10
 ein = Eexp/sum(dvol(is:i_inj,js:je,ks:ke))
 Tin = 1d3
 imu = 1d0/muconst
 pin = eos_p(damb,ein,Tin,imuconst)
 pamb = 1d-3*Eexp/sum(dvol(is:ie,js:je,ks:ke)) ! Low energy

 d(is:ie,js:je,ks:ke) = damb
 p(is:ie,js:je,ks:ke) = pamb
 p(is:i_inj,js:je,ks:ke) = pin
 v1(is:ie,js:je,ks:ke) = 0d0
 v2(is:ie,js:je,ks:ke) = 0d0
 v3(is:ie,js:je,ks:ke) = 0d0

return
end subroutine sedov
s
