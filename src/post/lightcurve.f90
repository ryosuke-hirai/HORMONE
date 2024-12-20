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

! PURPOSE: Program to post-process light curves

!############################################################################

program lightcurve

  use mpi_utils
  use mpi_domain
  use settings
  use grid
  use physval
  use constants
  use allocation_mod
  use boundary_mod
  use tools_mod
  use checksetup_mod
  use setup_mod
  use readbin_mod
  use output_mod
  use gridset_mod
  use profiler_mod
  use analysis_mod

  implicit none

  integer,parameter:: iangle=0, fangle=90, dangle=15
  integer:: outtn, angle, unitn(0:fangle)
  real(8):: outtime, rad, lum
  character(len=70):: file, outfile


!############################## start program ################################

  in_loop = .false.

  call init_mpi

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
  call checksetup

! Decompose domain onto MPI tasks
  call domain_decomp

! Start initial setup
  call allocations

  call readgrid('data/gridfile.bin')

  do angle = iangle, fangle, dangle
   write(outfile,'("data/",i2.2,"deg_lightcurve.dat")')angle
   if(start>0)then
    open(newunit=unitn(angle),file=outfile,status='old',position='append')
   else
    open(newunit=unitn(angle),file=outfile,status='replace')
   end if
  end do

  outtn = start
  outtime = dble(start)*dt_unit_in_sec
  do
   call set_file_name('bin',outtn,outtime,file)
   print*,'Reading ',file
   call readbin(file)
   call boundarycondition

   do angle = iangle, fangle, dangle
    rad = dble(angle)*pi/180d0
    call get_luminosity2(rad,lum)
    print*,time,angle,lum
    write(unitn(angle),'(i5,2(1PE14.6e2))')tn,time,lum
    flush(unitn(angle))
   end do
   
   outtime = outtime + dt_out
   outtn = outtn + tn_out
  end do

  do angle = iangle, fangle, dangle
   close(unitn(angle))
  end do

  call finalize_mpi

!------------------------------- end program ---------------------------------

 end program lightcurve
