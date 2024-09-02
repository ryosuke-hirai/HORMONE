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

program write_ascii

  use mpi_utils
  use mpi_domain
  use settings
  use grid
  use physval
  use allocation_mod
  use tools_mod
  use checksetup_mod
  use setup_mod
  use readbin_mod
  use output_mod
  use gridset_mod
  use profiler_mod

  implicit none

  integer:: outtn
  real(8):: outtime
  character(len=70):: file


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

  outtn = 0
  outtime = 0d0
  do
   call set_file_name('bin',outtn,outtime,file)
   print*,'Reading ',file
   call readbin(file)
   call write_plt
   outtime = outtime + dt_out
   outtn = outtn + tn_out
  end do

  call finalize_mpi

!------------------------------- end program ---------------------------------

 end program write_ascii
