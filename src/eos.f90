module pressure_mod
 use recombination_mod
 use constants,only:fac_egas,fac_pgas,arad,huge
 use physval,only:gamma
 use settings,only:eostype,eoserr
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE EOS_P_CF
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure and sound speed based on chosen EoS.
!          output p includes magnetic pressure
subroutine eos_p_cf(d,v1,v2,v3,b1,b2,b3,e,Tini,imu,p,cf,X,Y)
! PURPOSE: To calculate pressure and sound speed from density and internal energy
 implicit none
 real*8,intent( in):: d,v1,v2,v3,b1,b2,b3,e,Tini
 real*8,intent( in),optional:: X,Y
 real*8,intent(inout):: imu
 real*8,intent(out):: p,cf
 real*8:: eint, gamma_eff, bsq, T, corr
 real*8:: erec, imurec, derecdT, dimurecdT, Tdot, logd, dt
 real*8,parameter:: W4err = 1d-2
 integer n

!-----------------------------------------------------------------------------

 bsq = 0.5d0*(b1**2+b2**2+b3**2)
 
 eint = e - 0.5d0*d*(v1**2+v2**2+v3**2) - bsq

 select case (eostype)
 case(0) ! ideal gas
  p = (gamma-1d0)*eint + bsq
  cf = sqrt(gamma*p/d)

 case(1) ! ideal gas + radiation pressure
  T = Tini
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (eint - (arad*T**3+d*fac_egas*imu)*T) &
        / (  - 4d0*arad*T**3-d*fac_egas*imu)
   T = T - corr
  end do
  p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T + bsq
  gamma_eff = 1d0+p/eint
  cf = sqrt(gamma_eff*p/d)

 case(2) ! ideal gas + radiation pressure + recombination energy
  T = Tini
  corr = 1d99; Tdot=0d0;logd=log10(d);dt=0.9d0;n=0
  do while(abs(corr)>eoserr*T)
   call get_erec_imurec(logd,T,X,Y,erec,imurec,derecdT,dimurecdT)
   if(d*erec>=eint)then ! avoid negative thermal energy
    T = 0.9d0*T; Tdot = 0d0; cycle
   end if
   corr = (eint-(arad*T**3+d*fac_egas*imurec)*T-d*erec) &
        / ( -4d0*arad*T**3-d*(fac_egas*(imurec+dimurecdT*T)+derecdT) )
   if(abs(corr)>W4err*T)then
    T = T + Tdot*dt
    Tdot = (1d0-2d0*dt)*Tdot - dt*corr
   else
    T = T-corr
    Tdot = 0d0
   end if
   n=n+1
   if(n>50)dt=0.5d0
  end do
 
  imu = imurec
  p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T + bsq
  gamma_eff = 1d0+p/(eint-d*erec)

  cf = sqrt(gamma_eff*p/d)

 case default
  stop 'Error in eostype'
 end select

return
end subroutine eos_p_cf

! ***************************************************************************
real*8 function eos_p(d,eint,Tini,imu,X,Y)
! PURPOSE: To calculate pressure from density and internal energy without B field
 implicit none
 real*8,intent(in):: d,eint,Tini,imu
 real*8,intent(in),optional:: X,Y
 real*8:: T, corr, erec, imurec, derecdT, dimurecdT, Tdot, logd, dt
 real*8,parameter:: W4err = 1d-2
 integer n

 select case (eostype)
 case(0) ! ideal gas
  eos_p = (gamma-1d0)*eint
  
 case(1) ! ideal gas + radiation
  T = Tini
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (eint-(arad*T**3+d*fac_egas*imu)*T) &
        / ( -4d0*arad*T**3-d*fac_egas*imu)
   T = T - corr
  end do
  eos_p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T

 case(2) ! ideal gas + radiation + recombination
  T = Tini
  corr=1d99;Tdot=0d0;logd=log10(d);dt=0.9d0;n=0
  do while(abs(corr)>eoserr*T)
   call get_erec_imurec(logd,T,X,Y,erec,imurec,derecdT,dimurecdT)
   if(d*erec>=eint)then ! avoid negative thermal energy
    T = 0.9d0*T; Tdot=0d0;cycle
   end if
   corr = (eint-(arad*T**3+d*fac_egas*imurec)*T-d*erec) &
    / ( -4d0*arad*T**3-d*(fac_egas*(imurec+dimurecdT*T)+derecdT) )
   if(abs(corr)>W4err*T)then
    T = T + Tdot*dt
    Tdot = (1d0-2d0*dt)*Tdot - dt*corr
   else
    T = T-corr
    Tdot = 0d0
   end if
   n=n+1
   if(n>50)dt=0.5d0
  end do

  eos_p = ( fac_pgas*imurec*d + arad*T**3/3d0 )*T

 case default
  stop 'Error in eostype'
 end select

end function eos_p

! **************************************************************************

real*8 function eos_e(d,p,T,imu,X,Y)
! PURPOSE: To calculate internal energy from density and pressure
 implicit none
 real*8,intent(in):: d,p
 real*8,intent(in),optional:: X,Y
 real*8,intent(inout):: imu,T
 real*8:: corr, imurec, imurecold, dimurecdT, logd
 
 select case (eostype)
 case(0) ! ideal gas
  eos_e = p/(gamma-1d0)
  
 case(1) ! ideal gas + radiation
  corr = 1d99
  do while(abs(corr)>eoserr*T)
   corr = (p - (arad*T**3/3d0+d*fac_pgas*imu)*T) &
        / (-4d0*arad*T**3/3d0-d*fac_pgas*imu)
   T = T - corr
  end do
  eos_e = ( fac_egas*imu*d + arad*T**3 )*T
  
 case(2) ! ideal gas + radiation + recombination
  corr = 1d99; logd = log10(d)
  do while(abs(corr)>eoserr*T)
   call get_imurec(logd,T,X,Y,imurec,dimurecdT)
   corr = (p - (arad*T**3/3d0+d*fac_pgas*imurec)*T) &
        / (-4d0*arad*T**3/3d0-d*fac_pgas*(imurec+dimurecdT*T))
   T = T - corr
  end do

  imu = imurec
  eos_e = ( fac_egas*imurec*d + arad*T**3 )*T + d*get_erec(logd,T,X,Y)

 case default
  stop 'Error in eostype'
 end select

end function eos_e

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE PRESSURE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure in all cells based on EoS.

subroutine pressure
 use grid
 use physval

 implicit none

 real*8 vsq, bsq, corr, erec, imurec, derecdT, dimurecdT, Tdot, logd, dtau
 real*8,parameter:: W4err=1d-2
 integer nnn

!-----------------------------------------------------------------------------

 select case (eostype)
 case(0) ! ideal gas
!$omp parallel do private(i,j,k,vsq,bsq)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     vsq      = v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2 
     bsq      = b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2 
     eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq - 0.5d0*bsq
     T(i,j,k) = eint(i,j,k) / (d(i,j,k)*fac_egas*imu(i,j,k))
     p(i,j,k) = (gamma-1.d0) * eint(i,j,k)
     
     ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq

    end do
   end do
  end do
!$omp end parallel do

 case(1) ! ideal gas + radiation pressure
!$omp parallel do private(i,j,k,vsq,bsq,corr)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     vsq = v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2
     bsq = b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2
     eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq - 0.5d0*bsq
!eint(i,j,k) = max(eint(i,j,k),1d-4*e(i,j,k))
!    T(i,j,k) = eint(i,j,k) / (d(i,j,k)*fac_egas) ! initial guess
     corr = 1d99
     do while(abs(corr)>eoserr*T(i,j,k))
      corr = (eint(i,j,k) - ( arad*T(i,j,k)**3 &
                            + d(i,j,k)*fac_egas*imu(i,j,k) ) *T(i,j,k) )&
           / ( -4d0*arad*T(i,j,k)**3 &
               -d(i,j,k)*fac_egas*imu(i,j,k) )
      T(i,j,k) = T(i,j,k) - corr
     end do
     p(i,j,k) = ( d(i,j,k)*fac_pgas*imu(i,j,k) & ! gas pressure
                + arad*T(i,j,k)**3/3d0 ) *T(i,j,k) ! radiation pressure

     ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq
    end do
   end do
  end do
!$omp end parallel do
 case (2) ! ideal gas + radiation pressure + recombination energy
!$omp parallel do private(i,j,k,vsq,bsq,corr,erec,imurec,derecdT,dimurecdT,Tdot,logd,nnn,dtau)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     vsq = v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2
     bsq = b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2
     eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq - 0.5d0*bsq

     corr = 1d99;Tdot=0d0;logd=log10(d(i,j,k));dtau=0.9d0;nnn=0
! W4 method to give guess for Newton-Raphson
     do while(abs(corr)>eoserr*T(i,j,k))
      call get_erec_imurec(logd,T(i,j,k),spc(1,i,j,k),spc(2,i,j,k),&
                           erec,imurec,derecdT,dimurecdT)
      if(d(i,j,k)*erec>=eint(i,j,k))then ! avoid negative thermal energy
       T(i,j,k) = 0.9d0*T(i,j,k); Tdot = 0d0; cycle
      end if
      corr = (eint(i,j,k) - ( arad*T(i,j,k)**3 &
                            + d(i,j,k)*fac_egas*imurec ) *T(i,j,k) &
                            - d(i,j,k)*erec )&
                            
           / ( - 4d0*arad*T(i,j,k)**3 &
               - d(i,j,k)*(fac_egas*(imurec+dimurecdT*T(i,j,k))+derecdT) )
      if(abs(corr)>W4err*T(i,j,k))then
       T(i,j,k) = T(i,j,k) + Tdot*dtau
       Tdot = (1d0-2d0*dtau)*Tdot-dtau*corr
      else
       T(i,j,k) = T(i,j,k) - corr
       Tdot = 0d0
      end if
      nnn=nnn+1
      if(nnn>50)dtau=0.5d0
     end do

     imu(i,j,k) = imurec
     p(i,j,k) = ( d(i,j,k)*fac_pgas*imu(i,j,k) & ! gas pressure
                + arad*T(i,j,k)**3/3d0 ) *T(i,j,k) ! radiation pressure

     ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq
    end do
   end do
  end do
!$omp end parallel do
  
 case default
  stop 'Error in eostype'
 end select

 return
end subroutine pressure

end module pressure_mod

