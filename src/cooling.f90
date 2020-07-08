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
 use cooling_mod
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

contains

 real*8 function Yinv(yy,kk)
  implicit none
  integer,intent(in)::kk
  real*8,intent(in)::yy

  Yinv = Tint(kk)*&
       ( 1d0 - (1d0-alph(kk))*lint(kk)/lint(NN)*Tint(NN)/Tint(kk)&
              *(yy-Yint(kk)) ) **(1d0/(1d0-alph(kk)))

 end function Yinv

end subroutine cooling
