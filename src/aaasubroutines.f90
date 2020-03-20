module pressure_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE EOS_P_CF
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure and sound speed based on chosen EoS.

subroutine eos_p_cf(d,v1,v2,v3,b1,b2,b3,e,Tini,imu,p,cf)
! PURPOSE: To calculate pressure and sound speed from density and internal energy
 use funcs
 use settings,only:eostype,eoserr
 use constants
 use physval,only:gamma
 implicit none
 real*8,intent( in):: d,v1,v2,v3,b1,b2,b3,e,Tini,imu
 real*8,intent(out):: p,cf
 real*8:: eint, gamma_eff, bsq, T, corr

!-----------------------------------------------------------------------------

 bsq = 0.5d0*(b1*b1+b2*b2+b3*b3)
 
 eint = e - 0.5d0*d*(v1*v1+v2*v2+v3*v3) - bsq
! eint = abs(eint)
! eint = max(eint,1d-4*e)

 if(eostype==0)then ! ideal gas
  p = (gamma-1d0)*eint + bsq
  cf = sqrt(gamma*p/d)

 elseif(eostype==1)then ! ideal gas + radiation pressure
!  T = eint/(d*fac_egas) ! initial guess
  T = Tini
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (eint - (arad*pw(3,T)+d*fac_egas*imu)*T) &
        / (  - 4d0*arad*pw(3,T)-d*fac_egas*imu)
   T = T - corr
  end do
  p = ( fac_pgas*imu*d + arad*T*T*T/3d0 )*T + bsq
  gamma_eff = 1d0+p/eint
  cf = sqrt(gamma_eff*p/d)

 else
  stop 'Error in eostype'
 end if

return
end subroutine eos_p_cf


real*8 function eos_p(d,eint,Tini,imu)
! PURPOSE: To calculate pressure from density and internal energy
 use settings,only:eostype,eoserr
 use constants
 use physval,only:gamma
 implicit none
 real*8,intent(in):: d,eint,Tini,imu
 real*8:: T, corr, bsq

 bsq = 0d0

 if(eostype==0)then
  eos_p = (gamma-1d0)*eint + bsq
 elseif(eostype==1)then
!  T = eint/(d*fac_egas) ! initial guess
  T = Tini
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (eint-(arad*T*T*T+d*fac_egas*imu)*T) &
        / ( -4d0*arad*T*T*T-d*fac_egas*imu)
   T = T - corr
  end do
  eos_p = ( fac_pgas*imu*d + arad*T*T*T/3d0 )*T + bsq
 else
  stop 'Error in eostype'
 end if

end function eos_p


real*8 function eos_e(d,p,Tini,imu)
! PURPOSE: To calculate internal energy from density and pressure
 use settings,only:eostype,eoserr
 use constants
 use physval,only:gamma
 implicit none
 real*8,intent(in):: d,p,Tini,imu
 real*8:: T, corr

 if(eostype==0)then
  eos_e = p/(gamma-1d0)
 elseif(eostype==1)then
!  T = p/(d*fac_pgas) ! initial guess
  T = Tini
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (p - (arad*T*T*T/3d0+d*fac_pgas*imu)*T) &
        / (-4d0*arad*T*T*T    -d*fac_pgas*imu)
   T = T - corr
  end do
  eos_e = ( fac_egas*imu*d + arad*T*T*T )*T
 else
  stop 'Error in eostype'
 end if

end function eos_e

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE PRESSURE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure in all cells based on EoS.

subroutine pressure
 use funcs
 use grid
 use settings,only:eostype, eoserr
 use constants
 use physval
use gravmod

 implicit none

 real*8 vsq, bsq, corr

!-----------------------------------------------------------------------------

if(eostype==0)then ! ideal gas
!$omp parallel do private(i,j,k,vsq,bsq)
 do k = ks,ke
  do j = js,je
   do i = is,ie
    vsq      = pw(2, v1(i,j,k)) + pw(2, v2(i,j,k)) + pw(2, v3(i,j,k))
    bsq      = pw(2, b1(i,j,k)) + pw(2, b2(i,j,k)) + pw(2, b3(i,j,k))
    eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq - 0.5d0*bsq
    T(i,j,k) = eint(i,j,k) / (d(i,j,k)*fac_egas*imu(i,j,k))
    p(i,j,k) = (gamma-1.d0) * eint(i,j,k)

    ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq

   end do
  end do
 end do
!$omp end parallel do

elseif(eostype==1)then ! ideal gas + radiation pressure
!$omp parallel do private(i,j,k,vsq,bsq,corr)
 do k = ks,ke
  do j = js,je
   do i = is,ie
    vsq = v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k)
    bsq = b1(i,j,k)*b1(i,j,k) + b2(i,j,k)*b2(i,j,k) + b3(i,j,k)*b3(i,j,k)
    eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq - 0.5d0*bsq
eint(i,j,k) = max(eint(i,j,k),1d-4*e(i,j,k))
!    T(i,j,k) = eint(i,j,k) / (d(i,j,k)*fac_egas) ! initial guess
!    T(i,j,k) = T(i,j,k)
    corr = 1d99
if(d(i,j,k)<=0d0)print *,i,j,eint(i,j,k),d(i,j,k),e(i,j,k),grvphi(i-1:i+1,j,k-1:k+1)
!if(rungen==2)print *,i,k,eint(i,j,k),d(i,j,k),e(i,j,k)
    do while(abs(corr)>eoserr*T(i,j,k))
     corr = (eint(i,j,k) - ( arad*T(i,j,k)*T(i,j,k)*T(i,j,k) &
                           + d(i,j,k)*fac_egas*imu(i,j,k) ) *T(i,j,k) )&
          / ( -4d0*arad*T(i,j,k)*T(i,j,k)*T(i,j,k) &
              -d(i,j,k)*fac_egas*imu(i,j,k) )
     T(i,j,k) = T(i,j,k) - corr
    end do
    p(i,j,k) = ( d(i,j,k)*fac_pgas*imu(i,j,k) & ! gas pressure
               + arad*pw(3,T(i,j,k))/3d0 ) *T(i,j,k) ! radiation pressure

    ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq
   end do
  end do
 end do
!$omp end parallel do
else
 stop 'Error in eostype'
end if

return
end subroutine pressure

end module pressure_mod


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module minmod_mod

contains
subroutine minmod(mm,u,dx)

  implicit none

  real*8,intent(in ):: u(1:3),dx(1:2)
  real*8,intent(out):: mm
  real*8 x,y,z

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  x = 1.4d0*(u(3) - u(2))*dx(2)
  y = 1.4d0*(u(2) - u(1))*dx(1)
  z = (u(3)-u(1))*dx(1)*dx(2)/sum(dx)

  mm = sign(1.d0,x) * max(0.d0,min(abs(x),sign(1.d0,x)*y,sign(1.d0,x)*z))

return
end subroutine minmod

end module minmod_mod
