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
  use merger_mod
  use funcs

  use gravmod

  implicit none

! set parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                     sphrn, trnsn1, trnsn2, trnsn3
  namelist /timecon/ outstyle, endstyle, tnlim, t_end, dt_out, tn_out
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn
  namelist /simucon/ crdnt,courant, rktype, start
  namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                     bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
  namelist /gravcon/ gravswitch, grverr, cgerr, HGfac, hgcfl, &
                     include_extgrv, gis, gie, gjs, gje, gks, gke
  namelist /partcon/ include_particles, maxptc


!############################## start program ################################

! Initial setups -------------------------------------------------------------

  time = 0.d0; tn = 0
  heatdone = .true.
  inifile = 7d5
  inifile2= 7d5

  ! Reading parameters
  open(unit=1,file='parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=timecon)
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
  call boundarycondition
  call meanmolweight
  call conserve
  call pressure
  
! Initial output
  if(bc3os==9)call dirichletbound ! for ejecta
  if(include_particles.and.tn==0)call particles_setup
  if(include_particles.and.time==inifile)call particles_setup
  call boundarycondition
  call timestep
  if(gravswitch==3.and.tn==0)dt_old=dt / (courant*HGfac) * hgcfl

  call merger_setup

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then
   do while(time < t_end.and.tn<=tnlim)

    call timestep

    if(tn>0)call gravity
    if(bc3os==9)call dirichletbound ! for ejecta

    do rungen = 1, rktype
     call boundarycondition
     call numflux
     call source
     call rungekutta
    end do


    if(include_particles)call particles

    call merger

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

    if(endstyle==1)then
     if(time>=t_end)exit
    elseif(endstyle==2)then
     if(tn>=tnlim)exit
    end if

   end do
  end if
! End integration ------------------------------------------------------------

  call output ! To see final state

!------------------------------- end program ---------------------------------

  print *, 'Calculation complete! tn = ',tn

end program hormone
