module radtests_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE MATRAD_COUPLING
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To do matter-radiation coupling test

subroutine matrad_coupling

 use constants
 use grid
 use physval
 use radiation_mod,only:rad_heat_cool
 use pressure_mod,only:Trad

 real(8):: d_0,erad0,eint0,dt0,eint_equil,Tgas,TTrad,error
 real(8),parameter:: tolerance=1d-12

!-----------------------------------------------------------------------------

 d_0   = 1d-7
 erad0 = 1d12
 dt0   = 1d-12

! Heating test ===============================================================
 eint0 = 1d2 ! start from low gas temperature

 d=d_0;erad=erad0-eint0;eint=eint0;dt=dt0
 time=0d0
 do while(time<=1d-6)
  call rad_heat_cool
  time = time + dt
!  write(60,*)time,eint(1,1,1)/(fac_egas*d(1,1,1)*imu(1,1,1)),Trad(erad(1,1,1))
 end do
 eint_equil = fac_egas*d_0*imu(1,1,1)*Trad(erad(1,1,1))
 Tgas=eint(1,1,1)/(fac_egas*d(1,1,1)*imu(1,1,1))
 TTrad=Trad(erad(1,1,1))

 print*,'================================================='
 print*,'Heating test'
 print*,'================================================='
 print*,'Tgas/K =',Tgas
 print*,'Trad/K =',TTrad
 error = abs(Tgas/TTrad-1d0)
 print*,'Equilibriation level is: abs(e_int/e_equil-1) =',error
 if(error>tolerance)then
  print*,'Test failed!'
 else
  print*,'Passed!'
 end if


! Cooling test ===============================================================
 eint0 = 1d10 ! start from high gas temperature

 d=d_0;erad=erad0-eint0;eint=eint0;dt=dt0
 time=0d0
 do while(time<=1d-6)
  call rad_heat_cool
  time = time + dt
!  write(70,*)time,eint(1,1,1)/(fac_egas*d(1,1,1)*imu(1,1,1)),Trad(erad(1,1,1))
 end do
 stop
 eint_equil = fac_egas*d_0*imu(1,1,1)*Trad(erad(1,1,1))
 Tgas=eint(1,1,1)/(fac_egas*d(1,1,1)*imu(1,1,1))
 TTrad=Trad(erad(1,1,1))

 print*,'================================================='
 print*,'Cooling test'
 print*,'================================================='
 print*,'Tgas/K =',eint(1,1,1)/(fac_egas*d(1,1,1)*imu(1,1,1))
 print*,'Trad/K =',TTrad
 error = abs(Tgas/TTrad-1d0)
 print*,'Equilibriation level is: abs(e_int/e_equil-1) =',error
 if(error>tolerance)then
  print*,'Test failed!'
 else
  print*,'Passed!'
 end if

 stop
 
return
end subroutine matrad_coupling

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE RADIATION_SHOCK
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial conditions for radiation-dominated shock

subroutine radshock

 use constants,only:arad,Rgas
 use grid
 use physval

 real(8):: Teq

!-----------------------------------------------------------------------------

 Teq = 10d0 ! K
 d = 7.78d-10
 v1 = -6d5
 p = Rgas*d*Teq/muconst
 erad = arad*Teq**4

return
end subroutine radshock

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE DIFFUSION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Linear diffusion test

subroutine diffusion1d

 use settings,only:simtype
 use grid
 use physval
 use output_mod
 use radiation_mod,only:radiation

 integer:: strl, imid
 real(8):: Etilde

!-----------------------------------------------------------------------------

 strl = len(trim(simtype))

 Etilde = 1d5
 d = 1d0
 p = 1d0
 eint = p/(gamma-1d0)
 e = eint
 erad = 1d0
 
! Find the direction of shock tube
 select case(simtype(strl:strl))
 case('x')
  imid = (is_global+ie_global)/2
  erad(imid,js,ks) = Etilde/dxi1(imid)
 case('y')
  imid = (js_global+je_global)/2
  erad(is,imid,ks) = Etilde/dxi2(imid)
 case('z')
  imid = (ks_global+ke_global)/2
  erad(is,js,imid) = Etilde/dxi3(imid)
 end select

 return
end subroutine diffusion1d

end module radtests_mod
