!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!                                                                           !
!                                                                           !
!                                                                           !
!                             PROGRAM HORMONE                               !
!                                                                           !
!                                                                           !
!                                                                           !
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!                                                                           !
!    High ORder Magnetohydrodynamic cOde with Numerous Enhancements         !
!                                                                           !
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!

! PURPOSE: Main program to solve ideal MHD problems

!############################################################################

program hormone

  use settings
  use grid
  use physval
  use constants
  use pressure_mod
  use particle_mod
  use cooling_mod
  use dirichlet_mod
  use timestep_mod
  use output_mod
  use gravmod
  use composition_mod
  
  implicit none

! set parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                     sphrn, trnsn1, trnsn2, trnsn3
  namelist /out_con/ outstyle, endstyle, tnlim, t_end, dt_out, tn_out, &
                     dt_unit, sigfig, write_other_vel, write_shock
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn, include_cooling
  namelist /simucon/ crdnt,courant, rktype, start, mag_on, flux_limiter
  namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                     bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
  namelist /gravcon/ gravswitch, grverr, cgerr, HGfac, hgcfl, &
                     include_extgrv, gis, gie, gjs, gje, gks, gke
  namelist /partcon/ include_particles, maxptc


!############################## start program ################################

! Initial setups -------------------------------------------------------------

  time = 0.d0; tn = 0

  ! Reading parameters
  open(unit=1,file='parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=out_con)
  read(1,NML=eos_con)
  read(1,NML=simucon)
  read(1,NML=bouncon)
  read(1,NML=gravcon)
  read(1,NML=partcon)
  close(1)

  call checksetup
  call allocations
  call gridset
  call tools
  call metric
  call gravsetup

  call initialcondition
  if(dirichlet_on)call dirichletbound
  call boundarycondition
  call meanmolweight
  call conserve
  call pressure

! Initial output
  if(include_particles.and.tn==0)call particles_setup
!  if(include_particles.and.time==inifile)call particles_setup
  call boundarycondition
  call timestep
  if(gravswitch==3.and.tn==0)dt_old=dt / (courant*HGfac) * hgcfl

  if(tn==0)then
   call gravity
   call output
  end if

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then ! tnlim=0 to just output initial condition
   do while(time < t_end.and.tn<=tnlim)

    call timestep
    print*,tn,time,dt

    if(tn>0)call gravity
    if(dirichlet_on)call dirichletbound
    call shockfind

    do rungen = 1, rktype
     call boundarycondition
     call numflux
     call source
     call rungekutta
    end do

    if(include_cooling)  call cooling
    if(include_particles)call particles

    time = time + dt ; tn = tn + 1
! Output sequence ---------------------- !
    if(outstyle==2)then                  !
     if(tn/=0.and.mod(tn,tn_out)==0)then !
      call boundarycondition             !
      call output                        !
     end if                              !
    elseif(outstyle==1)then              !
     if(time>=t_out)then                 !
      call boundarycondition             !
      call output                        !
      t_out = t_out + dt_out             !
     end if                              !
    else                                 !
     print *,'outstyle out of range'     !
     stop                                !
    end if                               !
! -------------------------------------- !

! End sequence -------------!
    select case (endstyle)  ! 
    case(1) ! time up       ! 
     if(time>=t_end)exit    ! 
    case(2) ! timestep up   !
     if(tn>=tnlim)exit      ! 
    end select              ! 
! --------------------------!

   end do
  end if
! End integration ------------------------------------------------------------

  call output ! To see final state

!------------------------------- end program ---------------------------------

  print *, 'Calculation complete! tn = ',tn

end program hormone
