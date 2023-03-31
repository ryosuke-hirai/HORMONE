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
  use pressure_mod
  use particle_mod
  use cooling_mod
  use dirichlet_mod
  use timestep_mod
  use output_mod
  use gravmod
  use composition_mod
  
  implicit none

!############################## start program ################################

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

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then ! tnlim=0 to just output initial condition
   do

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
