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
  use utils_analysis

  implicit none

  integer,parameter:: iangle=0, fangle=90, dangle=15
  integer:: outtn, angle, unitn(0:fangle), iviewing_angle, n
  real(8):: outtime, rad, lum, Teff, Reff
  real(8),allocatable:: angles(:,:),lmda(:),spec(:)
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

!!$  do angle = iangle, fangle, dangle
!!$   write(outfile,'("data/",i2.2,"deg_lightcurve.dat")')angle
!!$   if(start>0)then
!!$    open(newunit=unitn(angle),file=outfile,status='old',position='append')
!!$   else
!!$    open(newunit=unitn(angle),file=outfile,status='replace')
!!$   end if
!!$  end do

  allocate(lmda(1:6))
  allocate(spec,mold=lmda)
  allocate(angles(1:2,0:8))

  angles(1:2,0) = [0d0      ,      0d0]! +z direction
  angles(1:2,1) = [0.5d0 *pi,      0d0]! +x direction, Equatorial plane 
  angles(1:2,2) = [0.5d0 *pi,       pi]! -x direction, Equatorial plane 
  angles(1:2,3) = [0.5d0 *pi, 0.5d0*pi]! +y direction, Equatorial plane 
  angles(1:2,4) = [0.5d0 *pi,-0.5d0*pi]! -y direction, Equatorial plane 
  angles(1:2,5) = [0.25d0*pi,      0d0]! +x direction, 45 deg angle 
  angles(1:2,6) = [0.25d0*pi,       pi]! -x direction, 45 deg angle 
  angles(1:2,7) = [0.25d0*pi, 0.5d0*pi]! +y direction, 45 deg angle 
  angles(1:2,8) = [0.25d0*pi,-0.5d0*pi]! -y direction, 45 deg angle

  lmda(1:6) = [4000,5000,6000,7000,8000,9000] ! in Angstrom
  lmda=lmda*1d-8 ! convert Angstrom to cm
  do iviewing_angle = 0, 8
   write(outfile,'("data/lightcurve_",i1.1,".dat")')iviewing_angle
   if(start>0)then
    open(newunit=unitn(iviewing_angle),file=outfile,status='old',&
         position='append')
   else
    open(newunit=unitn(iviewing_angle),file=outfile,status='replace')
    write(unitn(iviewing_angle),'("#",2(2X,a,1PE13.5e2))')&
     'theta=',angles(1,iviewing_angle),&
     'phi='  ,angles(2,iviewing_angle)
    write(unitn(iviewing_angle),'(a10)',advance='no')'tn'
    write(unitn(iviewing_angle),'(a14)',advance='no')'time'
    write(unitn(iviewing_angle),'(a14)',advance='no')'lum'
    do n = 1, size(lmda)
     write(unitn(iviewing_angle),'(9X,a1,i4)',advance='no')'I',int(lmda(n)*1d8)
    end do
    write(unitn(iviewing_angle),'(a14)',advance='no')'Teff'
    write(unitn(iviewing_angle),'(a14)',advance='no')'Reff'
    write(unitn(iviewing_angle),'()')
   end if
  end do


  outtn = start
  outtime = dble(start)*dt_unit_in_sec
  do
   call set_file_name('bin',outtn,outtime,file)
   print*,'Reading ',file
   call readbin(file)
   call boundarycondition

!!$   do angle = iangle, fangle, dangle
!!$    rad = dble(angle)*pi/180d0
!!$    call get_luminosity2(angle,lmda,lum,spec)
!!$    print*,time,angle,lum
!!$    flush(6)
!!$    write(unitn(angle),'(i10,2(1PE14.6e2))')tn,time,lum
!!$    flush(unitn(angle))
!!$   end do

   do iviewing_angle = 0, 8
    call get_luminosity3(angles(:,iviewing_angle),lmda,lum,spec)
    call fit_blackbody(lmda,spec,Teff)
    Reff = get_bbradius(Teff,lum)
    print'(1PE13.5e2,i2,3(1PE13.5e2))',time/dt_unit_in_sec,iviewing_angle,lum,Teff,Reff/rsun
    flush(6)
    write(unitn(iviewing_angle),'(i10,99(1PE14.5e3))')tn,time/3600d0/24d0,lum,spec,Teff,Reff/rsun
    flush(unitn(iviewing_angle))
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
