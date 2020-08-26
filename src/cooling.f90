module cooling_mod
 implicit none

 private
 integer,parameter:: NN = 5
 real*8:: Yint(0:NN), Tint(0:NN), lint(0:NN), alph(0:NN-1), tcool, lambda, Y
 real*8,parameter:: Tref = 1d8

 public cooling_setup,cooling
 
 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE COOLING_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set parameters required for cooling

subroutine cooling_setup

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

  
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE COOLING
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To implement cooling.

subroutine cooling

 use grid
 use settings,only:eostype,mag_on
 use physval
 use constants
 use pressure_mod

 implicit none

 real*8 YY

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,n,lambda,tcool,Y,YY)
 do k = ks, ke
  do j = js, je
   do i = is, ie

    if(T(i,j,k)<Tint(0).or.T(i,j,k)>=Tref)then
     cycle
    else
     do n = 0, NN-1
      if(T(i,j,k)>=Tint(n).and.T(i,j,k)<Tint(n+1))then
       lambda = lint(n)*(T(i,j,k)/Tint(n))**alph(n)
       tcool = kbol*amu/((gamma-1d0)*d(i,j,k)) &
             * (3d0*spc(1,i,j,k)+0.5d0*spc(2,i,j,k)+1d0) &
              /(spc(1,i,j,k)*(1d0+spc(1,i,j,k))) &
             * T(i,j,k)/lambda
       Y = Yint(n) + lint(NN)*Tint(n)/((1d0-alph(n))*lint(n)*Tint(NN)) &
                    *(1d0-(Tint(n)/T(i,j,k))**(alph(n)-1d0))
       YY = Y + T(i,j,k)*lint(NN)*dt/(Tref*lambda*tcool)
      end if
     end do
     if(YY>Yint(0))then
      T(i,j,k) = max(Tint(0),Yinv(YY,0))
     elseif(YY<Yint(NN))then
      T(i,j,k) = Yinv(YY,NN-1)
     else
      do n = 0, NN-1
       if(YY>=Yint(n).and.YY<Yint(n+1))then
        T(i,j,k) = Yinv(YY,n)
       end if
      end do
     end if
    end if

    eint(i,j,k) = fac_egas*imu(i,j,k)*d(i,j,k)*T(i,j,k)
    e(i,j,k) = eint(i,j,k) &
             + 0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2)
    if(mag_on) e(i,j,k) = e(i,j,k)+0.5d0*(b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2)
    select case (eostype)
    case(0:1) ! without recombination
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k), &
                      spc(1,i,j,k),spc(2,i,j,k))
    end select
    ptot(i,j,k) = p(i,j,k)

   end do
  end do
 end do
!$omp end parallel do

 call conserve

return

end subroutine cooling

real*8 function Yinv(yy,kk)
 implicit none
 integer,intent(in)::kk
 real*8,intent(in)::yy

 Yinv = Tint(kk)*&
      ( 1d0 - (1d0-alph(kk))*lint(kk)/lint(NN)*Tint(NN)/Tint(kk)&
             *(yy-Yint(kk)) ) **(1d0/(1d0-alph(kk)))

end function Yinv

end module cooling_mod
