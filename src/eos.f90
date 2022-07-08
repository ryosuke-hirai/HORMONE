module pressure_mod
 use ionization_mod
 use constants,only:fac_egas,fac_pgas,arad,huge
 use physval,only:gamma
 use settings,only:eostype,eoserr
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE EOS_P_CF
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure and sound speed based on chosen EoS.
!          output p includes magnetic pressure
subroutine eos_p_cf(d,b1,b2,b3,eint,Tini,imu,p,cf,X,Y,ierr)
! PURPOSE: To calculate pressure and sound speed from density and internal energy
 implicit none
 real*8,intent( in):: d,b1,b2,b3,eint,Tini
 real*8,intent( in),optional:: X,Y
 real*8,intent(inout):: imu
 real*8,intent(out):: p,cf
 integer,intent(out):: ierr
 real*8:: gamma_eff, bsq, T, corr
 real*8:: erec, imurec, derecdT, dimurecdT, Tdot, logd, dt
 real*8,parameter:: W4err = 1d-2
 integer n

!-----------------------------------------------------------------------------

 ierr=0
 bsq = 0.5d0*(b1**2+b2**2+b3**2)
 if(eint<=0d0.or.d<=0d0)then
  ierr=1
  return
 end if
 
 select case (eostype)
 case(0) ! ideal gas
  p = (gamma-1d0)*eint + bsq
  gamma_eff = 1d0+p/eint

  cf = sqrt(gamma_eff*p/d)

 case(1) ! ideal gas + radiation pressure
  T = Tini
  corr = huge
  do n = 1, 50
   corr = (eint - (arad*T**3+d*fac_egas*imu)*T) &
        / (  - 4d0*arad*T**3-d*fac_egas*imu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>50)then
   ierr=2
   print*,'Error in eos_p_cf, did not converge, eostype=',eostype
   print*,'d=',d,'e=',eint,'mu=',1d0/imu
   return
  end if
  p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T + bsq
  gamma_eff = 1d0+p/eint
  cf = sqrt(gamma_eff*p/d)

 case(2) ! ideal gas + radiation pressure + recombination energy
  T = Tini
  corr = huge; Tdot=0d0;logd=log10(d);dt=0.9d0;n=0
  do n = 1, 500
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
   if(n>50)dt=0.5d0
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>500)then
   ierr=2
   print*,'Error in eos_p_cf, did not converge, eostype=',eostype
   print*,'d=',d,'e=',eint,'mu=',1d0/imu
   return
  end if
 
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
real*8 function eos_p(d,eint,T,imu,X,Y)
! PURPOSE: To calculate pressure from density and internal energy without B field
 implicit none
 real*8,intent(in):: d,eint
 real*8,intent(inout):: T,imu
 real*8,intent(in),optional:: X,Y
 real*8:: corr, erec, derecdT, dimurecdT, Tdot, logd, dt
 real*8,parameter:: W4err = 1d-2
 integer n

 select case (eostype)
 case(0) ! ideal gas
  eos_p = (gamma-1d0)*eint
  
 case(1) ! ideal gas + radiation
  corr = huge
  do n = 1, 50
   corr = (eint-(arad*T**3+d*fac_egas*imu)*T) &
        / ( -4d0*arad*T**3-d*fac_egas*imu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>50)then
   print*,'Error in eos_p, eostype=',eostype
   print*,'d=',d,'eint=',eint,'mu=',1d0/imu
   stop
  end if
  eos_p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T

 case(2) ! ideal gas + radiation + recombination
  corr=huge;Tdot=0d0;logd=log10(d);dt=0.9d0
  do n = 1, 500
   call get_erec_imurec(logd,T,X,Y,erec,imu,derecdT,dimurecdT)
   if(d*erec>=eint)then ! avoid negative thermal energy
    T = 0.9d0*T; Tdot=0d0;cycle
   end if
   corr = (eint-(arad*T**3+d*fac_egas*imu)*T-d*erec) &
    / ( -4d0*arad*T**3-d*(fac_egas*(imu+dimurecdT*T)+derecdT) )
   if(abs(corr)>W4err*T)then
    T = T + Tdot*dt
    Tdot = (1d0-2d0*dt)*Tdot - dt*corr
   else
    T = T-corr
    Tdot = 0d0
   end if
   if(abs(corr)<eoserr*T)exit
   if(n>50)dt=0.5d0
  end do
  if(n>500)then
   print*,'Error in eos_p, eostype=',eostype
   print*,'d=',d,'eint=',eint,'mu=',1d0/imu
   stop
  end if

  eos_p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T

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
 integer n
 
 select case (eostype)
 case(0) ! ideal gas
  eos_e = p/(gamma-1d0)
  
 case(1) ! ideal gas + radiation
  corr = huge
  do n = 1, 50
   corr = (p - (arad*T**3/3d0+d*fac_pgas*imu)*T) &
        / (-4d0*arad*T**3/3d0-d*fac_pgas*imu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>50)then
   print*,'Error in eos_e, eostype=',eostype
   print*,'d=',d,'p=',p,'mu=',1d0/imu
   stop
  end if
  eos_e = ( fac_egas*imu*d + arad*T**3 )*T
  
 case(2) ! ideal gas + radiation + recombination
  corr = huge; logd = log10(d)
  do n = 1, 500
   call get_imurec(logd,T,X,Y,imurec,dimurecdT)
   corr = (p - (arad*T**3/3d0+d*fac_pgas*imurec)*T) &
        / (-4d0*arad*T**3/3d0-d*fac_pgas*(imurec+dimurecdT*T))
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>500)then
   print*,'Error in eos_e, eostype=',eostype
   print*,'d=',d,'p=',p,'mu=',1d0/imu
   stop
  end if
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

 real*8 bsq

!-----------------------------------------------------------------------------

 call internalenergy

 select case (eostype)
 case(0:1) ! EoSs that don't require composition
!$omp parallel do private(i,j,k,bsq)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     bsq = b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k)) ! gets T too

     ptot(i,j,k) = p(i,j,k) + 0.5d0*bsq
    end do
   end do
  end do
!$omp end parallel do

 case (2) ! EoSs that require composition
!$omp parallel do private(i,j,k,bsq)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     bsq = b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k),&
                      spc(1,i,j,k),spc(2,i,j,k)) ! gets T too

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

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE INTERNALENERGY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute eint from etot, v and B over the whole grid

subroutine internalenergy

 use grid
 use physval

 implicit none

 integer:: ierr

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    eint(i,j,k) = get_eint(e(i,j,k),d(i,j,k),&
                           v1(i,j,k),v2(i,j,k),v3(i,j,k),&
                           b1(i,j,k),b2(i,j,k),b3(i,j,k),ierr )
    if(ierr==1)then
     print*,'Error in internalenergy, Negative internal energy'
     print'("i=",i4,"j=",i4,"k=",i4)',i,j,k
     print*,'etot=',e(i,j,k),'eint=',eint(i,j,k)
     stop
    end if
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine internalenergy

! **************************************************************************

function get_eint(etot,d,v1,v2,v3,b1,b2,b3,ierr) result(eint)
!PURPOSE: To calculate eint from etot, v and B
 implicit none
 real*8,intent(in):: etot,d,v1,v2,v3,b1,b2,b3
 integer,intent(out),optional::ierr
 real*8:: eint, vsq, bsq

 ierr=0
 vsq = v1**2 + v2**2 + v3**2
 bsq = b1**2 + b2**2 + b3**2

 eint = etot-0.5d0*d*vsq-0.5d0*bsq
 if(present(ierr).and.eint<=0d0)ierr=1
 
end function get_eint

end module pressure_mod

