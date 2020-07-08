!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE COOLING_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set parameters required for cooling

subroutine cooling_setup

 use cooling_mod

 implicit none

 integer k

!-----------------------------------------------------------------------------
! All based on Townsend 2009

! Give a cooling efficiency function in piecewise power-law form
! Not suitable for alph=1

 Tint(0) = 1d4
 Tint(1) = 1.73780082874938d4
 Tint(2) = 3.16227766016838d4
 Tint(3) = 1d5
 Tint(4) = 2.51188643150958d7
 Tint(5) = Tref

 alph(0) =  3.94899234386478d0
 alph(1) = -1.32649386004528d0
 alph(2) =  1.56058212317378d0
 alph(3) = -0.672204056342891d0
 alph(4) =  0.277914373766375d0

 lint(0) = 4.44841730085144d-39*Tint(0)**alph(0)
 lint(1) = 7.83525736700239d-17*Tint(1)**alph(1)
 lint(2) = 7.91899595725077d-30*Tint(2)**alph(2)
 lint(3) = 1.39019574514672d-18*Tint(3)**alph(3)
 lint(4) = 1.60990554337853d-25*Tint(4)**alph(4)
 lint(5) = lint(4)*(Tref/Tint(NN-1))**alph(NN-1)

! Calculate Y_k values
 Yint(NN) = 0d0
 do k = NN-1, 0, -1
  Yint(k) = Yint(k+1) - &
            lint(NN)*Tint(k)/((1d0-alph(k))*lint(k)*Tint(NN))&
            *(1d0-(Tint(k)/Tint(k+1))**(alph(k)-1d0))
 end do

return
end subroutine cooling_setup
