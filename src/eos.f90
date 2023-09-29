module pressure_mod
 use constants,only:fac_egas,fac_pgas,arad,huge
 use physval,only:gamma
 use settings,only:eostype,eoserr
 implicit none
 
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GETT_FROM_DE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine getT_from_de(d,eint,T,imu,X,Y,erec_out)
 use ionization_mod,only:get_erec_imurec
! PURPOSE: To calculate temperature from density and internal energy
 implicit none
 real(8),intent(in):: d,eint
 real(8),intent(inout):: T,imu
 real(8),intent(in),optional:: X,Y
 real(8),intent(out),optional:: erec_out
 real(8):: corr, erec, derecdT, dimurecdT, Tdot, logd, dt
 real(8),parameter:: W4err = 1d-2
 integer:: n

 if(T<=0d0) T=1d3

 select case (eostype)
 case(0) ! ideal gas
  T = eint/(fac_egas*imu)
  
 case(1) ! ideal gas + radiation
  corr = huge
  do n = 1, 50
   corr = (eint-(arad*T**3+d*fac_egas*imu)*T) &
        / ( -4d0*arad*T**3-d*fac_egas*imu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>50)then
   print*,'Error in getT_from_de, eostype=',eostype
   print*,'d=',d,'eint=',eint,'mu=',1d0/imu
   stop
  end if

 case(2) ! ideal gas + radiation + recombination
  corr=huge;Tdot=0d0;logd=log10(d);dt=0.9d0
  do n = 1, 500
   call get_erec_imurec(logd,T,X,Y,erec,imu,derecdT,dimurecdT)
   if(d*erec>=eint)then ! avoid negative thermal energy
    T = 0.95d0*T; Tdot=0d0;cycle
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
   print*,'Error in getT_from_de, eostype=',eostype
   print*,'d=',d,'eint=',eint,'mu=',1d0/imu
   stop
  end if
  if(present(erec_out)) erec_out = erec
  
 case default
  stop 'Error in eostype'
 end select

end subroutine getT_from_de


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GETT_FROM_DP
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine getT_from_dp(d,p,T,imu,X,Y,erec)
 use ionization_mod,only:get_erec,get_imurec
! PURPOSE: To calculate temperature from density and pressure
 implicit none
 real(8),intent(in):: d,p
 real(8),intent(in),optional:: X,Y
 real(8),intent(inout):: T,imu
 real(8),intent(out),optional:: erec
 real(8):: corr, imurec, dimurecdT, logd
 integer:: n

 if(T<=0d0) T=1d3
 
 select case (eostype)
 case(0) ! ideal gas
  T = p/(fac_pgas*imu)
  
 case(1) ! ideal gas + radiation
  corr = huge
  do n = 1, 500
   corr = (p - (arad*T**3/3d0+d*fac_pgas*imu)*T) &
        / (-4d0*arad*T**3/3d0-d*fac_pgas*imu)
   T = T - corr
   if(abs(corr)<eoserr*T)exit
  end do
  if(n>500)then
   print*,'Error in getT_from_dp, eostype=',eostype
   print*,'d=',d,'p=',p,'mu=',1d0/imu
   stop
  end if
  
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
   print*,'Error in getT_from_dp, eostype=',eostype
   print*,'d=',d,'p=',p,'mu=',1d0/imu
   stop
  end if
  imu = imurec
  if(present(erec)) erec = get_erec(logd,T,X,Y)

 case default
  stop 'Error in eostype'
 end select

end subroutine getT_from_dp

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE EOS_P_CS
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure and thermal sound speed from d and eint

subroutine eos_p_cs(d,eint,T,imu,p,cs,X,Y,ierr)

 implicit none
 real(8),intent( in):: d,eint
 real(8),intent( in),optional:: X,Y
 real(8),intent(inout):: T,imu
 real(8),intent(out):: p,cs
 integer,intent(out):: ierr
 real(8):: gamma_eff,erec,Ttemp

!-----------------------------------------------------------------------------

 ierr=0

 if(eint<=0d0.or.d<=0d0)then
  ierr=1
  return
 end if
 
 select case (eostype)
 case(0) ! ideal gas
  p = eos_p(d,eint,T,imu)
  gamma_eff = 1d0+p/eint

  cs = sqrt(gamma_eff*p/d)

 case(1) ! ideal gas + radiation pressure
  p = eos_p(d,eint,T,imu)
  gamma_eff = 1d0+p/eint
  
  cs = sqrt(gamma_eff*p/d)

 case(2) ! ideal gas + radiation pressure + recombination energy
  Ttemp = T
  call getT_from_de(d,eint,Ttemp,imu,X,Y,erec)
  p = eos_p(d,eint,T,imu,X,Y)
  gamma_eff = 1d0+p/(eint-d*erec)

  cs = sqrt(gamma_eff*p/d)

 case default
  stop 'Error in eostype'
 end select

return
end subroutine eos_p_cs

! ***************************************************************************
function eos_p(d,eint,T,imu,X,Y)
! PURPOSE: To calculate pressure from density and internal energy without B field
 implicit none
 real(8),intent(in):: d,eint
 real(8),intent(inout):: T,imu
 real(8),intent(in),optional:: X,Y
 real(8):: eos_p

 select case (eostype)
 case(0) ! ideal gas
  call getT_from_de(d,eint,T,imu)
  eos_p = (gamma-1d0)*eint
  
 case(1) ! ideal gas + radiation
  call getT_from_de(d,eint,T,imu)
  eos_p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T

 case(2) ! ideal gas + radiation + recombination
  call getT_from_de(d,eint,T,imu,X,Y)
  eos_p = ( fac_pgas*imu*d + arad*T**3/3d0 )*T

 case default
  stop 'Error in eostype'
 end select

end function eos_p

! **************************************************************************

function eos_e(d,p,T,imu,X,Y)
! PURPOSE: To calculate internal energy from density and pressure
 implicit none
 real(8),intent(in):: d,p
 real(8),intent(in),optional:: X,Y
 real(8),intent(inout):: imu,T
 real(8):: eos_e, erec
 
 select case (eostype)
 case(0) ! ideal gas
  call getT_from_dp(d,p,T,imu)
  eos_e = p/(gamma-1d0)
  
 case(1) ! ideal gas + radiation
  call getT_from_dp(d,p,T,imu)
  eos_e = ( fac_egas*imu*d + arad*T**3 )*T
  
 case(2) ! ideal gas + radiation + recombination
  call getT_from_dp(d,p,T,imu,X,Y,erec)
  eos_e = ( fac_egas*imu*d + arad*T**3 )*T + d*erec

 case default
  stop 'Error in eostype'
 end select

end function eos_e

! ***************************************************************************

pure function get_cf(d,cs,b1,b2,b3) result(cf)
! PURPOSE: To calculate full sound speed from thermal sound speed and B field
!          Note that this is the sound speed in the b1 direction
 real(8),intent(in):: d,cs,b1,b2,b3
 real(8):: cf,asq,bbsq,bb1,bb2,bb3
 asq = cs**2

 bb1  = b1**2 / d
 bb2  = b2**2 / d
 bb3  = b3**2 / d
 bbsq = bb1 + bb2 + bb3

 cf = sqrt( 0.5d0*( asq + bbsq + sqrt( (asq+bbsq)**2 - 4d0*asq*bb1 ) ) )

end function get_cf

! ***************************************************************************

function entropy_from_dT(d,T,imu,X,Y) result(entropy)
 use constants
 use ionization_mod,only:get_xion
! PURPOSE: To calculate entropy from density and temperature
 implicit none
 real(8),intent(in):: d,T
 real(8),intent(in),optional:: X,Y
 real(8),intent(inout):: imu
 real(8):: p,entropy,S_ion,S_rad,S_ele,n_x,n_y,n_z,n_e,fac,eta,xion(1:4)

 select case(eostype)
 case(0) ! ideal gas
  p = fac_pgas*T*imu
  entropy = p/d**gamma

 case(1) ! ideal gas + radiation
  S_ion = fac_pgas*imu*log(T**1.5d0/d)
  S_rad = 4d0*arad*T**3/(3d0*d)
  entropy = (S_ion + S_rad) / fac_pgas
  
 case(2) ! ideal gas + radiation + recombination
  call get_xion(log(d),T,X,Y,xion)
  n_x = d*X/amu
  n_y = d*Y/(4d0*amu)
  n_z = d*(1d0-X-Y)/(12d0*amu)
  n_e = n_x*xion(2)+n_y*sum(xion(3:4))+6d0*n_z
  fac = 2d0*pi*amu*kbol/hplanck**2
  S_ion = X*(log((fac*T)**1.5d0/n_x)+2.5d0) &
        + Y*(log((fac*T)**1.5d0/n_y)+2.5d0)*0.25d0
  S_ele = (2.5d0*(X+1d0)/2d0-eta*n_e/d)*amu
  S_rad = 4d0*arad*T**3/(3d0*d*kbol)*amu
  entropy = S_ion + S_ele + S_rad

! TODO LIST:
! 1. Write a routine to calculate the degeneracy parameter eta given d,T,X,Y
! 2. Add ionization entropy
! 3. Is S_ion correct with the current expression for partial degeneracy?
  
  print*, 'Entropy calculation for eostype=2 still needs work'
  stop

 case default
  stop 'Error in eostype'
  
 end select
  
end function entropy_from_dT

! ***************************************************************************

function entropy_from_dp(d,p,T,imu,X,Y) result(entropy)
 use constants
! PURPOSE: To calculate entropy from density and pressure
 implicit none
 real(8),intent(in):: d,p
 real(8),intent(in),optional:: X,Y
 real(8),intent(inout):: T,imu
 real(8):: entropy

 select case(eostype)
 case(0) ! ideal gas
  call getT_from_dp(d,p,T,imu)
  entropy = p/d**gamma

 case(1) ! ideal gas + radiation
  call getT_from_dp(d,p,T,imu)
  entropy = entropy_from_dT(d,T,imu)
  
 case(2) ! ideal gas + radiation + recombination
  call getT_from_dp(d,p,T,imu,X,Y)
  entropy = entropy_from_dT(d,T,imu,X,Y)

 case default
  stop 'Error in eostype'
  
 end select
  
end function entropy_from_dp

! ***************************************************************************

function entropy_from_de(d,e,T,imu,X,Y) result(entropy)
 use constants
! PURPOSE: To calculate entropy from density and internal energy
 implicit none
 real(8),intent(in):: d,e
 real(8),intent(in),optional:: X,Y
 real(8),intent(inout):: T,imu
 real(8):: entropy

 select case(eostype)
 case(0) ! ideal gas
  call getT_from_de(d,e,T,imu)
  entropy = (gamma-1d0)*e/d**gamma

 case(1) ! ideal gas + radiation
  call getT_from_de(d,e,T,imu)
  entropy = entropy_from_dT(d,T,imu)
  
 case(2) ! ideal gas + radiation + recombination
  call getT_from_de(d,e,T,imu,X,Y)
  entropy = entropy_from_dT(d,T,imu,X,Y)

 case default
  stop 'Error in eostype'
  
 end select
  
end function entropy_from_de

! **************************************************************************

function get_d_from_ps(p,S,imu,X,Y) result(d)
! PURPOSE: Calculate density given pressure and entropy​
 implicit none

 real(8),intent(in):: p,S
 real(8),intent(inout):: imu
 real(8),intent(in),optional:: X,Y
 real(8):: d,corr,corr0,dp,S0,Sp,dSdd,T,ddot,dt
 real(8),parameter:: dfac=1d-8, Serr_rel_eoserr=7d0, W4err=1d-2
 integer:: n

!-----------------------------------------------------------------------------

 select case (eostype)
 case(0) ! ideal gas
  d = (p/S)**(1d0/gamma)

 case(1) ! ideal gas + radiation (fully ionized, uniform composition)
  d = 1d-8 ! initial guess
  T = p/(fac_pgas*d*imu)
  corr=1d99
  do n = 1, 5000 ! Newton-Raphson iteration to get density
   dp = d*(1d0+dfac)
   S0 = entropy_from_dp(d,p,T,imu)
   Sp = entropy_from_dp(dp,p,T,imu)
   dSdd = (Sp-S0)/(dp-d) ! Numerical differentiation
   corr = min((S0-S)/dSdd,d*0.9d0)
   d = d-corr
   if(abs((S0-S)/dSdd)<Serr_rel_eoserr*eoserr*d)exit
  end do
  if(n>5000)then
   print*,'Error in get_d_from_ps, eostype=',eostype
   print*,'d=',d,'S=',S,'mu=',1d0/imu
   stop
  end if
  
 case(2) ! ideal gas + radiation + recombination
  d = 1d-8 ! initial guess
  T = 1d3
  corr=huge;ddot=0d0;dt=0.9d0
  do n = 1,500
   dp = d*(1d0+dfac)
   S0 = entropy_from_dp(d,p,T,imu,X,Y)
   Sp = entropy_from_dp(dp,p,T,imu,X,Y)
   dSdd = (Sp-S0)/(dp-d) ! Numerical differentiation
   corr0 = (S0-S)/dSdd
   corr = min(corr0,d*0.9d0)
   if(abs(corr)>W4err*d)then
    d = d + ddot*dt
    ddot = (1d0-2d0*dt)*ddot - dt*corr
   else
    d = d-corr
    ddot = 0d0
   end if
   if(abs(corr)<Serr_rel_eoserr*eoserr*d)exit
   if(n>50)dt=0.5d0
  end do
  if(n>500)then
   print*,'Error in get_d_from_ps, eostype=',eostype
   print*,'d=',d,'S=',S,'mu=',1d0/imu
   stop
  end if

 case default
  stop 'Error in eostype'
 end select

return
end function get_d_from_ps

! **************************************************************************

function get_e_from_ds(d,S,imu,X,Y) result(e)
! PURPOSE: Calculate internal energy given density and entropy​
 implicit none

 real(8),intent(in):: d,S
 real(8),intent(inout):: imu
 real(8),intent(in),optional:: X,Y
 real(8):: e,corr,corr0,ep,S0,Sp,dSde,T,edot,dt
 real(8),parameter:: dfac=1d-8, Serr_rel_eoserr=1d2, W4err=1d-2
 integer:: n

!-----------------------------------------------------------------------------

 select case (eostype)
 case(0) ! ideal gas
  e = (S*d**gamma)/(gamma-1d0)

 case(1) ! ideal gas + radiation (fully ionized, uniform composition)
  e = 1d-8 ! initial guess
  T = e/(fac_egas*d*imu)
  corr=1d99
  do n = 1, 500 ! Newton-Raphson iteration to get internal energy
   ep = e*(1d0+dfac)
   S0 = entropy_from_de(d,e,T,imu)
   Sp = entropy_from_de(d,ep,T,imu)
   dSde = (Sp-S0)/(ep-e) ! Numerical differentiation
   corr = min((S0-S)/dSde,e*0.9d0)
   e = e-corr
   if(abs((S0-S)/dSde)<Serr_rel_eoserr*eoserr*e)exit
  end do
  if(n>500)then
   print*,'Error in get_e_from_ds, eostype=',eostype
   print*,'d=',d,'S=',S,'mu=',1d0/imu
   stop
  end if
  
 case(2) ! ideal gas + radiation + recombination
  e = 1d-8 ! initial guess
  T = 1d3
  corr=huge;edot=0d0;dt=0.9d0
  do n = 1,500
   ep = e*(1d0+dfac)
   S0 = entropy_from_de(d,e,T,imu,X,Y)
   Sp = entropy_from_de(d,ep,T,imu,X,Y)
   dSde = (Sp-S0)/(ep-e) ! Numerical differentiation
   corr0 = (S0-S)/dSde
   corr = min(corr0,e*0.9d0)
   if(abs(corr)>W4err*e)then
    e = e + edot*dt
    edot = (1d0-2d0*dt)*edot - dt*corr
   else
    e = e-corr
    edot = 0d0
   end if
   if(abs(corr)<Serr_rel_eoserr*eoserr*e)exit
   if(n>50)dt=0.5d0
  end do
  if(n>500)then
   print*,'Error in getd_from_ps, eostype=',eostype
   print*,'d=',d,'S=',S,'mu=',1d0/imu
   stop
  end if

 case default
  stop 'Error in eostype'
 end select

return
end function get_e_from_ds

! **************************************************************************

function get_p_from_ds(d,S,imu,X,Y) result(p)
! PURPOSE: Calculate internal energy given density and entropy​
 implicit none

 real(8),intent(in):: d,S
 real(8),intent(inout):: imu
 real(8),intent(in),optional:: X,Y
 real(8):: eint,p,T

!-----------------------------------------------------------------------------
 
 select case(eostype)
 case(0) ! ideal gas
  p = S*d**gamma

 case(1) ! ideal gas + radiation (fully ionized, uniform composition)
  eint = get_e_from_ds(d,S,imu)
  p = eos_p(d,eint,T,imu)

 case(2) ! ideal gas + radiation + recombination
  eint = get_e_from_ds(d,S,imu,X,Y)
  p = eos_p(d,eint,T,imu,X,Y)

 case default
  stop 'Error in eostype'
  
 end select

 return
end function get_p_from_ds

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE PRESSURE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate pressure in all cells based on EoS.

subroutine pressure

 use settings,only:mag_on,eostype
 use grid,only:is,ie,js,je,ks,ke
 use physval,only:p,ptot,b1,b2,b3,d,eint,T,imu,spc

 real(8):: bsq=0d0
 integer:: i,j,k

!-----------------------------------------------------------------------------

 call internalenergy

 select case (eostype)
 case(0:1) ! EoSs that don't require composition
!$omp parallel do private(i,j,k,bsq) collapse(3)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k)) ! gets T too
     ptot(i,j,k) = p(i,j,k)

     if(.not.mag_on)cycle
     bsq = b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2
     ptot(i,j,k) = ptot(i,j,k) + 0.5d0*bsq
    end do
   end do
  end do
!$omp end parallel do

 case (2) ! EoSs that require composition
!$omp parallel do private(i,j,k,bsq) collapse(3)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k),&
                      spc(1,i,j,k),spc(2,i,j,k)) ! gets T too
     ptot(i,j,k) = p(i,j,k)

     if(.not.mag_on)cycle
     bsq = b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2
     ptot(i,j,k) = ptot(i,j,k) + 0.5d0*bsq
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

 use grid,only:is,ie,js,je,ks,ke
 use physval,only:d,e,eint,v1,v2,v3,b1,b2,b3

 implicit none

 logical:: err=.false.
 integer:: i,j,k,ierr,erri,errj,errk

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,ierr) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    eint(i,j,k) = get_eint(e (i,j,k),d (i,j,k),&
                           v1(i,j,k),v2(i,j,k),v3(i,j,k),&
                           b1(i,j,k),b2(i,j,k),b3(i,j,k),ierr )
    if(ierr==1)then
     err = .true.
     erri = i; errj = j; errk = k
    end if
   end do
  end do
 end do
!$omp end parallel do

 if(err)then
  i = erri; j = errj; k = errk
  print*,'Error in internalenergy, Negative internal energy'
  print'("i=",i4,", j=",i4,", k=",i4)',i,j,k
  print*,'etot=',e(i,j,k),'eint=',eint(i,j,k)
  stop
 end if
 
return
end subroutine internalenergy

! **************************************************************************

function get_eint(etot,d,v1,v2,v3,b1,b2,b3,ierr) result(eint)
!PURPOSE: To calculate eint from etot, v and B
 use settings,only:mag_on
 implicit none
 real(8),intent(in):: etot,d,v1,v2,v3,b1,b2,b3
 integer,intent(out),optional::ierr
 real(8):: eint, vsq, bsq

 ierr=0
 vsq = v1**2 + v2**2 + v3**2

 eint = etot-0.5d0*d*vsq
 if(mag_on)then
  bsq = b1**2 + b2**2 + b3**2
  eint = eint - 0.5d0*bsq
 end if

 if(present(ierr).and.eint<=0d0)ierr=1
 
end function get_eint

end module pressure_mod

