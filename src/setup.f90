module setup_mod
 implicit none
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE READ_STARTFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read startfile

 subroutine read_startfile

  use settings,only:start,simtype,parafile

!-----------------------------------------------------------------------------

  namelist /startcon/ start, simtype, parafile
  
  open(unit=1,file='startfile',status='old')
  read(1,NML=startcon)
  close(1)
  
  return
 end subroutine read_startfile


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE READ_DEFAULT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Find and read default parameter file for a given simtype

subroutine read_default

 use settings
 use grid

 integer:: ui
 character*50:: filename, basefile

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
 case('sodshock_x','briowushock_x')
  filename='../para/parameters_shocktube_x'
 case('sodshock_y','briowushock_y')
  filename='../para/parameters_shocktube_y'
 case('sodshock_z','briowushock_z')
  filename='../para/parameters_shocktube_z'
 case('sedov')
  filename='../para/parameters_sedov'
 case('rsg_sph')
  filename='../para/parameters_rsg_sph'
 case default
  print*,'This simtype does not exist yet, simtype ="',trim(simtype),'"'
  stop
 end select

 call read_parameters(filename)

return
end subroutine read_default

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE READ_PARAMETERS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Read namelists in parameters file

subroutine read_parameters(filename)

 use settings
 use grid
 use physval
 use particle_mod,only:maxptc

 integer:: ui
 character*50,intent(in)::filename

!-----------------------------------------------------------------------------

 namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                    is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                    x1min, x2min, x3min, &
                    sphrn, trnsn16, trnsn8, trnsn4, trnsn2
 namelist /out_con/ outstyle, endstyle, tnlim, t_end, dt_out, tn_out, &
                    dt_unit, sigfig, outres, write_other_vel, write_shock, &
                    write_evo
 namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn
 namelist /simucon/ crdnt,courant, rktype, mag_on, flux_limiter, &
                    include_cooling, include_extforce
 namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                    bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
 namelist /gravcon/ gravswitch, grverr, cgerr, HGfac, hgcfl, gbtype, &
                    grav_init_other, include_extgrv, &
                    gis, gie, gjs, gje, gks, gke
 namelist /partcon/ include_particles, maxptc

 if(filename=='')return
 
 open(newunit=ui,file=filename,status='old')
 read(ui,NML=gridcon)
 read(ui,NML=out_con)
 read(ui,NML=eos_con)
 read(ui,NML=simucon)
 read(ui,NML=bouncon)
 read(ui,NML=gravcon)
 read(ui,NML=partcon)
 close(ui)
  
return
end subroutine read_parameters

end module setup_mod
