module checksetup_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE CHECKSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
subroutine checksetup

! purpose: To check whether the setups are consistent.

  use settings
  use grid
  use physval
  use tests_mod

!-----------------------------------------------------------------------------

! Select which directions to solve
  solve_i=.false.;solve_j=.false.;solve_k=.false.
  if(ie>is)solve_i=.true.
  if(je>js)solve_j=.true.
  if(ke>ks)solve_k=.true.

! Reading dimension of grid
  dim=0
  if(solve_i)dim=dim+1
  if(solve_j)dim=dim+1
  if(solve_k)dim=dim+1
  if(dim==0)then
   print*,"Error from gridset"
   stop
  end if
  if(dim/=2) write_other_vel = .false.
  if(dim/=3) write_other_slice = .false.
  if(crdnt/=2) write_mc = .false.

! Count number of equations
  ufnmax = 0
  call add_equation(icnt,ufnmax) ! Continuity equation
  call add_equation(imo1,ufnmax) ! Momentum equation 1
  call add_equation(imo2,ufnmax) ! Momentum equation 2
  call add_equation(imo3,ufnmax) ! Momentum equation 3
  call add_equation(iene,ufnmax) ! Energy equation
  if(mag_on)then
   call add_equation(img1,ufnmax) ! Magnetic field equation 1
   call add_equation(img2,ufnmax) ! Magnetic field equation 2
   call add_equation(img3,ufnmax) ! Magnetic field equation 3
   call add_equation(i9wv,ufnmax) ! Divergence cleaning (9-wave method)
  end if

! Set uniform mesh if that dimension is not used
  if(.not.solve_i)imesh=0
  if(.not.solve_j)jmesh=0
  if(.not.solve_k)kmesh=0

! To warn the periodic boundary condition
  if(bc1is*bc1os==0.and.bc1is+bc1os/=0)then
   print *,"Error from boundary condition 'bc1is,bc1os'"
   print *,"Both should be periodic"
   stop
  elseif(bc2is*bc2os==0.and.bc2is+bc2os/=0)then
   print *,"Error from boundary condition 'bc2is,bc2os'"
   print *,"Both should be periodic"
   stop
  elseif(bc3is*bc3os==0.and.bc3is+bc3os/=0)then
   print *,"Error from boundary condition 'bc3is,bc3os'"
   print *,"Both should be periodic"
   stop
  elseif(bc1iv*bc1ov==0.and.bc1iv+bc1ov/=0)then
   print *,"Error from boundary condition 'bc1iv,bc1ov'"
   print *,"Both should be periodic"
   stop
  elseif(bc2iv*bc2ov==0.and.bc2iv+bc2ov/=0)then
   print *,"Error from boundary condition 'bc2iv,bc2ov'"
   print *,"Both should be periodic"
   stop
  elseif(bc3iv*bc3ov==0.and.bc3iv+bc3ov/=0)then
   print *,"Error from boundary condition 'bc3iv,bc3ov'"
   print *,"Both should be periodic"
   stop
  end if

! Setting boundary conditions for 1D and 2D simulations
  if(.not.solve_i)then
   bc1is=2 ; bc1os=2 ; bc1iv=2 ; bc1ov=2
  end if

  if(.not.solve_j)then
   bc2is=3 ; bc2os=3 ; bc2iv=3 ; bc2ov=3
  end if

  if(.not.solve_k)then
   bc3is=3 ; bc3os=3 ; bc3iv=3 ; bc3ov=3
  end if

! Set boundary condition for polar coordinates
  if(crdnt==1)then
   bc2is=0 ; bc2os=0 ; bc2iv=0 ; bc2ov=0
  elseif(crdnt==2)then
   bc2is=1 ; bc2os=1 ; bc2iv=1 ; bc2ov=1
   bc3is=0 ; bc3os=0 ; bc3iv=0 ; bc3ov=0
  end if

! Setting boundary condition for equatorial symmetry
  if(eq_sym.and.crdnt==1)then
   bc3is=1
   bc3iv=1
  elseif(eq_sym.and.crdnt==2)then
   bc2os=1
   bc2ov=1
  end if

! Set flag if any boundary is set to Dirichlet boundary condition
  dirichlet_on = .false.
  if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or.&
     bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9)then
   dirichlet_on = .true.
  end if

! Set flag if any boundary is set to flux boundary condition
  fluxbound_on = .false.
  if(bc1is==10.or.bc1os==10.or.bc2is==10.or.bc2os==10.or.bc3is==10.or.bc3os==10.or.&
     bc1iv==10.or.bc1ov==10.or.bc2iv==10.or.bc2ov==10.or.bc3iv==10.or.bc3ov==10)then
   fluxbound_on = .true.
  end if

! check CFL condition
  if(courant>1d0)then
   print *,"Error from courant number, courant = ",courant
   print *,"courant number should be <= 1.0"
   stop
  end if

! check gravity
  if(dim==1.and.gravswitch==2)then
   print *,"Self-gravity calculations are not necessary for 1D calculations"
   stop
  end if
  if(crdnt/=2.and.gravswitch==1)then
   print *,"Point source gravity routine only for spherical coordinates"
   stop
  end if

  if(gravswitch/=0.and.crdnt==1)then
   gis = is ; gjs = js ; gje = je
  elseif(gravswitch/=0.and.crdnt==2)then
   gis = is ; gjs = js ; gje = je ; gks = ks ; gke = ke
  elseif(gravswitch==0)then
   gis = is ; gie = ie ; gjs = js ; gje = je ; gks = ks ; gke = ke
  end if

  gis = min(gis,is) ; gjs = min(gjs,js) ; gks = min(gks,ks)
  gie = max(gie,ie) ; gje = max(gje,je) ; gke = max(gke,ke)

  select case(dt_unit)
  case('yr')
   dt_unit_in_sec = 3600d0*24d0*365.25d0
  case('d')
   dt_unit_in_sec = 3600*24d0
  case('hr')
   dt_unit_in_sec = 3600d0
  case('ks')
   dt_unit_in_sec = 1000d0
  case('min')
   dt_unit_in_sec = 60d0
  case('s')
   dt_unit_in_sec = 1d0
  case('ms')
   dt_unit_in_sec = 1d-3
  case('ns')
   dt_unit_in_sec = 1d-9
  case default
   print*,'Error: Set a valid unit for output interval, dt_unit = ',dt_unit
   stop
  end select
  dt_out = dt_out*dt_unit_in_sec
  t_out = dt_out
  t_end = t_end*dt_unit_in_sec

! Spherical composition only for spherical coordinates
  if(compswitch==1.and.(crdnt/=2.or.max(je-js,ke-ks)>0))then
   print *,"compswitch is not consistent with coordinate", compswitch
   stop
  end if

! Cooling only for ideal gas EoS
  if(include_cooling.and.eostype>0)then
   stop "Cooling is only for optically thin equation of state"
  end if

! Recombination only if composition is provided
  if(compswitch<=1.and.eostype==2)then
   print *,'Provide composition to use recombination module'
  end if

! Provide more than 2 species for composition
  if(compswitch>=1) spn = max(spn,2)

! Don't allocate spc if compswitch == 0
  if(compswitch==0) spn = 0

! Check if test data is available
  if(is_test)call check_testlist(simtype)

! Check if fmr_max is reasonable for the grid resolution
  if(2**(fmr_max-1) > max(je,ke))then
   print*,'fmr_max should be smaller so that 2**(fmr_max-1) >= max(je,ke)'
   print*,'fmr_max=',fmr_max,', max(je,ke)=',max(je,ke)
   stop
  end if

! EoS can only be specific types if radiation is switched on
  if(radswitch>0.and.eostype==1)then
   print*,'Cannot select gas+radiation EoS if radswitch>0'
   print*,'eostype=',eostype
   stop
  end if

! Output temperature if eostype>=1
  if(eostype>=1)write_temp = .true.

  return
 end subroutine checksetup

 subroutine add_equation(i,ufnmax)
  integer,intent(out)::i
  integer,intent(inout)::ufnmax
  ufnmax = ufnmax + 1
  i = ufnmax
 end subroutine add_equation

end module checksetup_mod
