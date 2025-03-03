module setup_mod
 implicit none

 public:: read_startfile,read_default,read_parameters
 private::is_it_test
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE READ_STARTFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read startfile

 subroutine read_startfile

  use settings,only:start,simtype,parafile,is_test

!-----------------------------------------------------------------------------

  namelist /startcon/ start, simtype, parafile

  open(unit=1,file='startfile',status='old')
  read(1,NML=startcon)
  close(1)

  call is_it_test(simtype,is_test)
  if(is_test) parafile = ''

  return
 end subroutine read_startfile

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE IS_IT_TEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To check whether the simulation is a test to be compared against pre-computed data

 subroutine is_it_test(simtype,is_test)

  character(len=*),intent(inout)::simtype
  logical,intent(out):: is_test
  integer:: strl

!-----------------------------------------------------------------------------

  strl = len(trim(simtype))
  if(simtype(strl-4:strl)==' test')then
   is_test = .true.
   simtype = simtype(1:strl-5)
  else
   is_test = .false.
  end if

  return
 end subroutine is_it_test


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE READ_DEFAULT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Find and read default parameter file for a given simtype

subroutine read_default

 use settings
 use grid

 character(len=50):: filename, basefile

!-----------------------------------------------------------------------------
! Read base file to set all parameters to default
 basefile='../para/parameters_base'
 call read_parameters(basefile)

! Read default file relevant for the simtype
 select case(simtype)
 case('eostest')
  filename='../para/parameters_eostest'
 case('orszagtang_xy')
  filename='../para/parameters_orszagtang_xy'
 case('orszagtang_yz')
  filename='../para/parameters_orszagtang_yz'
 case('orszagtang_xz')
  filename='../para/parameters_orszagtang_xz'
 case('sodshock_x')
  filename='../para/parameters_sodshock_x'
 case('sodshock_y')
  filename='../para/parameters_sodshock_y'
 case('sodshock_z')
  filename='../para/parameters_sodshock_z'
 case('briowushock_x')
  filename='../para/parameters_briowushock_x'
 case('briowushock_y')
  filename='../para/parameters_briowushock_y'
 case('briowushock_z')
  filename='../para/parameters_briowushock_z'
 case('other_shocktube_x')
  filename='../para/parameters_othershock_x'
 case('other_shocktube_y')
  filename='../para/parameters_othershock_y'
 case('other_shocktube_z')
  filename='../para/parameters_othershock_z'
 case('sedov_default','sedov_other')
  filename='../para/parameters_sedov_default'
 case('eruption')
  filename='../para/parameters_eruption'
 case('KHinstability')
  filename='../para/parameters_KHinstability'
 case('rad_box')
  filename='../para/parameters_rad_box'
 case('star_sph')
  filename='../para/parameters_star_sph'
 case('rsg_sph')
  filename='../para/parameters_rsg_sph'
 case('polytrope')
  filename='../para/parameters_polytrope'
 case('agndisk')
  filename='../para/parameters_agndisk'
 case('windtunnel')
  filename='../para/parameters_windtunnel'
 case('stellarcollision')
  filename='../para/parameters_stellarcollision'
 case('sn2022jli')
  filename='../para/parameters_sn2022jli'
 case('modify')
  filename='../para/parameters_rsg_sph'!temporary
 case('smearing')
  filename='../para/parameters_smearing'
 case('iotest')
  filename='../para/parameters_iotest'
 case default
  print*,'This simtype does not exist yet, simtype ="',trim(simtype),'"'
  stop
 end select

 call read_parameters(filename)

return
end subroutine read_default

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE READ_PARAMETERS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read namelists in parameters file

subroutine read_parameters(filename)

 use settings
 use grid
 use physval

 integer:: ui,istat
 character(len=*),intent(in)::filename

!-----------------------------------------------------------------------------

 namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                    is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                    x1min, x2min, x3min, &
                    fmr_max, fmr_lvl
 namelist /out_con/ outdir, outstyle, endstyle, &
                    tnlim, t_end, dt_out, tn_out, tn_evo, dt_unit, &
                    sigfig, outres, output_ascii, write_other_vel, &
                    write_shock, write_evo, write_other_slice, write_temp, &
                    write_mc
 namelist /eos_con/ eostype, gamma, eoserr, compswitch, muconst, spn
 namelist /simucon/ crdnt, courant, rktype, mag_on, flux_limiter, alpha9wave,&
                    include_cooling, include_extforce, frame, extrasfile
 namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                    bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
 namelist /gravcon/ gravswitch, grvsrctype, grverr, cgerr, &
                    HGfac, hgcfl, maxtngrv, gbtype, grktype, alphagrv, &
                    grav_init_other, grav_init_relax, include_extgrv,&
                    gis, gie, gjs, gje, gks, gke
 namelist /sinkcon/ include_sinks, nsink, include_accretion
 namelist /rad_con/ radswitch, opacitytype, lambdatype
 namelist /partcon/ include_particles, maxptc
 namelist /testcon/ test_tol, Mach_tol

 if(filename=='')return

 open(newunit=ui,file=filename,status='old',iostat=istat)
 if(istat/=0)then
  print*,'Error: Simulation parameter file cannot be found'
  print'(3a)','Missing file name = "',trim(filename),'"'
  stop
 end if
 read(ui,NML=gridcon,iostat=istat);rewind(ui)
 read(ui,NML=out_con,iostat=istat);rewind(ui)
 read(ui,NML=eos_con,iostat=istat);rewind(ui)
 read(ui,NML=simucon,iostat=istat);rewind(ui)
 read(ui,NML=bouncon,iostat=istat);rewind(ui)
 read(ui,NML=gravcon,iostat=istat);rewind(ui)
 read(ui,NML=sinkcon,iostat=istat);rewind(ui)
 read(ui,NML=rad_con,iostat=istat);rewind(ui)
 read(ui,NML=partcon,iostat=istat);rewind(ui)
 read(ui,NML=testcon,iostat=istat)
 close(ui)

return
end subroutine read_parameters

end module setup_mod
