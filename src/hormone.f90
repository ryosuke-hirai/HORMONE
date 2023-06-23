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
  use composition_mod
  use boundary_mod
  use numflux_mod
  use source_mod
  use rungekutta_mod

  use gravity_mod
  use particle_mod
  use cooling_mod
  use dirichlet_mod
  use shockfind_mod
  use tests_mod
  use omp_lib
  
  implicit none

!############################## start program ################################

  wtime = 0d0
  wtime(iini) = -omp_get_wtime()
  
! Initial setups -------------------------------------------------------------

  time = 0.d0; tn = 0

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

  if(write_evo)then
   call open_evofile
   call evo_output
  end if

  wtime(iini) = wtime(iini) + omp_get_wtime()

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then ! tnlim=0 to just output initial condition

   wtime(itot) = wtime(itot) - omp_get_wtime()

   main_loop:do

    call timestep
    print'(a,i8,2(3X,a,1PE13.5e2))','tn =',tn,'time =',time,'dt =',dt

    if(tn>0)call gravity
    if(dirichlet_on)call dirichletbound
    call shockfind

    do rungen = 1, rktype
     call boundarycondition
     call numflux
     call source
     call rungekutta
    end do

    if(mag_on)           call phidamp
    if(include_cooling)  call cooling
    if(include_particles)call particles

    time = time + dt ; tn = tn + 1
! Output sequence ---------------------- !
    select case(outstyle)                !
    case(1) ! output by time             !
     if(time>=t_out)then                 !
      call boundarycondition             !
      call output                        !
      t_out = t_out + dt_out             !
     end if                              !
    case(2) ! output by timestep         !
     if(tn/=0.and.mod(tn,tn_out)==0)then !
      call boundarycondition             !
      call output                        !
     end if                              !
    case default                         !
     print *,'outstyle out of range'     !
     stop                                !
    end select                           !
    if(write_evo.and.&                   !
       mod(tn,10)==0)call evo_output     !
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
  end if
! End integration ------------------------------------------------------------

  if(tn/=0)call output ! To see final state

  if(is_test)call test

!------------------------------- end program ---------------------------------

  print *, 'Calculation complete! tn = ',tn

  call scaling_output

end program hormone
