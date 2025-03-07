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
 use radiation_mod,only:radiation
 use pressure_mod,only:Trad

 real(8):: d_0,erad0,eint0,dt0,eint_equil,Tgas,TTrad,error
 real(8),parameter:: tolerance=1d-5

!-----------------------------------------------------------------------------

 d_0   = 1d-7
 erad0 = 1d12
 dt0   = 1d-19

! Heating test ===============================================================
 eint0 = 1d6 ! start from low gas temperature

 d=d_0;erad=erad0-eint0;eint=eint0;e=eint0;T=eint/(d*Cv*imu);dt=dt0
 p=(gamma-1d0)*eint
 time=0d0
 do while(time<=1d-6)
  call radiation
  time = time + dt
  dt = dt * 1.01d0
!  write(60,*)time,eint(1,1,1)/(Cv*d(1,1,1)*imu(1,1,1)),Trad(erad(1,1,1))
 end do

 eint_equil = Cv*d_0*imu(1,1,1)*Trad(erad(1,1,1))
 Tgas=eint(1,1,1)/(Cv*d(1,1,1)*imu(1,1,1))
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

 d=d_0;erad=erad0-eint0;eint=eint0;e=eint0;T=eint/(d*Cv*imu);dt=dt0
 time=0d0
 do while(time<=1d-6)
  call radiation
  time = time + dt
  dt = dt * 1.01d0
!  write(70,*)time,eint(1,1,1)/(Cv*d(1,1,1)*imu(1,1,1)),Trad(erad(1,1,1))
 end do

 eint_equil = Cv*d_0*imu(1,1,1)*Trad(erad(1,1,1))
 Tgas=eint(1,1,1)/(Cv*d(1,1,1)*imu(1,1,1))
 TTrad=Trad(erad(1,1,1))

 print*,'================================================='
 print*,'Cooling test'
 print*,'================================================='
 print*,'Tgas/K =',eint(1,1,1)/(Cv*d(1,1,1)*imu(1,1,1))
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

 use settings,only:simtype,extrasfile
 use constants,only:arad,Rgas
 use grid
 use physval
 use opacity_mod

 integer:: ui,strl,istat
 real(8):: d_back,Teq,v_piston,kapparho

!-----------------------------------------------------------------------------

 namelist /rshocon/ d_back,Teq,v_piston,kapparho

! Set default values
 open(newunit=ui,file='../para/extras_radshock',status='old')
 read(ui,NML=rshocon)
 close(ui)

! Override with user-set values if present
 open(newunit=ui,file=extrasfile,status='old',iostat=istat)
 if(istat==0)read(ui,NML=rshocon,iostat=istat)
 close(ui)

 d(is:ie,js:je,ks:ke) = d_back
 p(is:ie,js:je,ks:ke) = Rgas*d(is:ie,js:je,ks:ke)*Teq/muconst
 erad(is:ie,js:je,ks:ke) = arad*Teq**4
 c_kappa_p = kapparho/d_back
 c_kappa_r = c_kappa_p

! Find the direction of shock tube
 strl = len(trim(simtype))
 select case(simtype(strl:strl))
 case('x')
  v1 = -v_piston
 case('y')
  v2 = -v_piston
 case('z')
  v3 = -v_piston
 end select

 return
end subroutine radshock

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE DIFFUSION
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
 d(is:ie,js:je,ks:ke) = 1d0
 p(is:ie,js:je,ks:ke) = 1d0
 eint(is:ie,js:je,ks:ke) = p(is:ie,js:je,ks:ke)/(gamma-1d0)
 e(is:ie,js:je,ks:ke) = eint(is:ie,js:je,ks:ke)
 erad(is:ie,js:je,ks:ke) = 1d0
 
! Find the direction of shock tube
 select case(simtype(strl:strl))
 case('x')
  imid = (is_global+ie_global)/2
  erad(imid,js:je,ks:ke) = Etilde/dxi1(imid)
 case('y')
  imid = (js_global+je_global)/2
  erad(is:ie,imid,ks:ke) = Etilde/dxi2(imid)
 case('z')
  imid = (ks_global+ke_global)/2
  erad(is:ie,js:je,imid) = Etilde/dxi3(imid)
 end select

 return
end subroutine diffusion1d

end module radtests_mod
