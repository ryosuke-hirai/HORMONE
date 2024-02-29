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
  use setup_mod
  use allocation_mod
  use tools_mod
  use checksetup_mod
  use gridset_mod
  use pressure_mod
  use conserve_mod
  use timestep_mod
  use metric_mod
  use initialcondition_mod
  use output_mod
  use gravmod
  use sink_mod
  use composition_mod
  use boundary_mod
  use numflux_mod
  use source_mod
  use rungekutta_mod

  use gravity_mod
  use radiation_mod
  use particle_mod
  use cooling_mod
  use dirichlet_mod
  use shockfind_mod
  use tests_mod

  use profiler_mod
  
  implicit none

!############################## start program ################################

! Start profiling
  call init_profiler
  
! Initial setups -------------------------------------------------------------

  call start_clock(wtini)

  time = 0d0; tn = 0

! Read startfile
  call read_startfile

! Read default parameter file
  call read_default

! Reading parameters
  call read_parameters(parafile)

! Start initial setup
  call checksetup
  call allocations
  call gridset
  call metric
  call tools

  call initialcondition
  if(dirichlet_on)call dirichletbound
  call boundarycondition
  call meanmolweight
  call conserve
  call pressure

! Initial output
  call open_evofile
  call open_sinkfile

  if(include_particles.and.tn==0)call particles_setup
  call timestep

  if(gravswitch==3.and.tn==0)dt_old=dt / (courant*HGfac) * hgcfl

  if(tn==0)then
   call gravity
   if(include_sinks)call sinkfield
   call output
  end if

  call stop_clock(wtini)

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then ! tnlim=0 to just output initial condition

   wtime(wtgri) = wtime(wtgrv)
   call reset_clock(wtlop)
   call start_clock(wtlop)

   main_loop:do

    call timestep

    call gravity
    if(include_sinks)call sink_motion

    print'(a,i8,2(3X,a,1PE13.5e2))','tn =',tn,'time =',time,'dt =',dt

    call start_clock(wthyd)

    if(dirichlet_on) call dirichletbound
    call shockfind

    do rungen = 1, rktype
     call boundarycondition
     call numflux
     call source
     call rungekutta
    end do

    call stop_clock(wthyd)

    if(mag_on)           call phidamp
    if(radswitch>0)      call radiation
    if(include_cooling)  call cooling
    if(include_particles)call particles

    time = time + dt ; tn = tn + 1

! Output sequence ---------------------- !
    select case(outstyle)                !
    case(1) ! output by time             !
     if(time>=t_out)then                 !
      call output                        !
      t_out = t_out + dt_out             !
     end if                              !
    case(2) ! output by timestep         !
     if(tn/=0.and.mod(tn,tn_out)==0)then !
      call output                        !
     end if                              !
    case default                         !
     print *,'outstyle out of range'     !
     stop                                !
    end select                           !
    if(write_evo.and.mod(tn,10)==0)then  !
     call evo_output                     !
     call sink_output                    !
    end if                               !
! -------------------------------------- !

! End sequence ------------------- !
    select case (endstyle)         !
    case(1) ! time up              !
     if(time>=t_end)exit main_loop ! 
    case(2) ! timestep up          !
     if(tn>=tnlim)exit main_loop   !
    end select                     !
! -------------------------------- !

   end do main_loop
   call stop_clock(wtlop)

  end if
! End integration ------------------------------------------------------------

  call scaling_output

  if(tn/=0)call output ! To see final state

  if(is_test)call test

!------------------------------- end program ---------------------------------

  print *, 'Calculation complete! tn = ',tn

end program hormone
