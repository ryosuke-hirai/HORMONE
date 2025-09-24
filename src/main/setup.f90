module setup_mod
 implicit none

 public:: read_startfile,read_default,read_parameters
 private::is_it_test
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READ_STARTFILE
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
!                            SUBROUTINE IS_IT_TEST
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
!                           SUBROUTINE READ_DEFAULT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Find and read default parameter file for a given simtype

subroutine read_default

 use settings
 use grid

 character(len=50):: filename, basefile
 logical:: parafile_available

!-----------------------------------------------------------------------------
! Read base file to set all parameters to default
 basefile='../para/parameters_base'
 call read_parameters(basefile)

! Read default file relevant for the simtype
 filename = '../para/parameters_'//trim(simtype)
 inquire(file=filename,exist=parafile_available)

 if(.not.parafile_available)then
  print*,'This simtype does not exist yet'
  print'(3a)',' simtype = "',trim(simtype),'"'
  stop
 end if

 call read_parameters(filename)

return
end subroutine read_default

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READ_PARAMETERS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read namelists in parameters file

subroutine read_parameters(filename)

 use settings
 use grid
 use physval
 use opacity_mod

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
 namelist /eos_con/ eostype, gamma, eoserr, compswitch, muconst, spn, &
                    X_uniform, Y_uniform, Z_uniform
 namelist /simucon/ crdnt, courant, rktype, mag_on, flux_limiter, alpha9wave,&
                    include_cooling, include_extforce, frame, extrasfile, &
                    solve_hydro, fixed_dt
 namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                    bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
 namelist /gravcon/ gravswitch, grvsrctype, grverr, cgerr, &
                    HGfac, hgcfl, maxtngrv, gbtype, grktype, alphagrv, &
                    grav_init_other, grav_init_relax, include_extgrv,&
                    gis, gie, gjs, gje, gks, gke
 namelist /sinkcon/ include_sinks, nsink, include_accretion, p_acc
 namelist /rad_con/ radswitch, opacitytype, c_kappa_p, c_kappa_r, c_kappa_f, &
                    lambdatype, rbtype
 namelist /partcon/ include_particles, maxptc
 namelist /mat_con/ matrix_solver
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
 read(ui,NML=mat_con,iostat=istat);rewind(ui)
 read(ui,NML=testcon,iostat=istat)
 close(ui)

return
end subroutine read_parameters

end module setup_mod
