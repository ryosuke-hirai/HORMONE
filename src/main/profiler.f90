module profiler_mod
 implicit none

 public:: init_profiler,profiler_output1,start_clock,stop_clock,reset_clock
 integer,parameter:: n_wt=38 ! number of profiling categories
 real(8):: wtime(0:n_wt),wtime_max(0:n_wt),wtime_min(0:n_wt),wtime_avg(0:n_wt),imbalance(0:n_wt)
 integer,parameter:: &
  wtini=1 ,& ! initial conditions
  wtgri=2 ,& ! initial gravity
  wtou1=3 ,& ! initial output
  wtlop=4 ,& ! main loop
  wthyd=5 ,& ! hydrodynamics
  wtflx=6 ,& ! numerical flux
  wtrng=7 ,& ! Runge-Kutta
  wtbnd=8 ,& ! boundary conditions
  wtsrc=9 ,& ! source terms
  wtint=10,& ! MUSCL interpolation
  wteos=11,& ! equation of state
  wtsmr=12,& ! smearing
  wtsm2=13,& ! MPI sweep for smearing
  wttim=14,& ! time stepping
  wtgrv=15,& ! gravity
  wtgbn=16,& ! gravbound
  wtelg=17,& ! elliptic self-gravity
  wtmig=18,& ! MICCG solver for gravity
  wtpeg=19,& ! PETSc solver for gravity
  wthyp=20,& ! hyperbolic self-gravity
  wtgsm=21,& ! gravity smearing
  wtgs2=22,& ! MPI sweep gravity smearing
  wtsho=23,& ! shockfind
  wtrad=24,& ! radiation
  wtmir=25,& ! MICCG solver for radiation
  wtmra=26,& ! MICCG A matrix
  wtper=27,& ! PETSc solver for radiation
  wtprv=28,& ! PETSc vector assembly
  wtpra=29,& ! PETSc A matrix
  wtprc=30,& ! PETSc A matrix coefficients
  wtprm=31,& ! PETSc A matrix MPI
  wtopc=32,& ! opacity
  wtrfl=33,& ! radiative flux
  wtsnk=34,& ! sink motion
  wtacc=35,& ! sink accretion
  wtout=36,& ! output
  wtmpi=37,& ! mpi exchange
  wtwai=38,& ! mpi wait
  wttot=0    ! total
 integer,public:: parent(0:n_wt),maxlbl
 character(len=30),public:: routine_name(0:n_wt)
 logical,public:: clock_on(0:n_wt)
 real(8),parameter,private:: thres=1d-10

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE INIT_PROFILER
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To initialize variables for profiler

subroutine init_profiler

!-----------------------------------------------------------------------------

 wtime = 0d0

! Default parent=0

 parent(wtini) = wttot ! initial conditions
 parent(wtgri) = wtini ! initial gravity
 parent(wtou1) = wtini ! initial output
 parent(wtlop) = wttot ! main loop
 parent(wthyd) = wtlop ! hydrodynamics
 parent(wtflx) = wthyd ! numerical flux
 parent(wtrng) = wthyd ! Runge-Kutta
 parent(wtbnd) = wthyd ! boundary conditions
 parent(wtsrc) = wthyd ! source terms
 parent(wtint) = wthyd ! MUSCL interpolation
 parent(wteos) = wthyd ! equation of state
 parent(wtsmr) = wthyd ! smearing
 parent(wtsm2) = wtsmr ! MPI sweep for smearing
 parent(wttim) = wtlop ! time stepping
 parent(wtgrv) = wtlop ! gravity
 parent(wtgbn) = wtgrv ! gravbound
 parent(wtelg) = wtgrv ! elliptic self-gravity
 parent(wtmig) = wtelg ! MICCG solver for gravity
 parent(wtpeg) = wtelg ! PETSc solver for gravity
 parent(wthyp) = wtgrv ! hyperbolic self-gravity
 parent(wtgsm) = wthyp ! gravity smearing
 parent(wtgs2) = wtgsm ! MPI sweep gravity smearing
 parent(wtout) = wtlop ! output
 parent(wtsho) = wtlop ! shockfind
 parent(wtrad) = wtlop ! radiation
 parent(wtmir) = wtrad ! MICCG solver for radiation
 parent(wtmra) = wtrad ! MICCG A matrix
 parent(wtper) = wtrad ! PETSc solver for radiation
 parent(wtprv) = wtper ! PETSc vector assembly
 parent(wtpra) = wtrad ! PETSc A assembly
 parent(wtprc) = wtpra ! PETSc A coefficients
 parent(wtprm) = wtpra ! PETSc A MPI assembly
 parent(wtopc) = wtrad ! opacity
 parent(wtrfl) = wtrad ! radiative flux
 parent(wtsnk) = wtlop ! sink motion
 parent(wtacc) = wtlop ! sink accretion
 parent(wtmpi) = wtlop ! mpi exchange
 parent(wtwai) = wtlop ! mpi wait
 parent(wttot) =-1     ! Total

! Make sure to keep routine name short
 routine_name(wtini) = 'Setup'       ! initial conditions
 routine_name(wtgri) = '1st Gravity' ! initial gravity
 routine_name(wtou1) = '1st Output'  ! initial output
 routine_name(wtlop) = 'Main loop'   ! main loop
 routine_name(wthyd) = 'Hydro'       ! hydrodynamics
 routine_name(wtflx) = 'Numflux'     ! numerical flux
 routine_name(wtrng) = 'RungeKutta'  ! Runge-Kutta
 routine_name(wtbnd) = 'Boundary'    ! boundary conditions
 routine_name(wtsrc) = 'Source'      ! source terms
 routine_name(wtint) = 'Interpolate' ! MUSCL interpolation
 routine_name(wteos) = 'EoS'         ! equation of state
 routine_name(wtsmr) = 'Smearing'    ! smearing
 routine_name(wtsm2) = 'MPI sweep'   ! MPI sweep for smearing
 routine_name(wttim) = 'Timestep'    ! time stepping
 routine_name(wtgrv) = 'Gravity'     ! gravity
 routine_name(wtgbn) = 'Gravbound'   ! gravbound
 routine_name(wtelg) = 'Elliptic'    ! elliptic self-gravity
 routine_name(wtmig) = 'MICCG'       ! MICCG solver for gravity
 routine_name(wtpeg) = 'PETSc'       ! PETSc solver for gravity
 routine_name(wthyp) = 'Hyperbolic'  ! hyperbolic self-gravity
 routine_name(wtgsm) = 'Grav smear'  ! gravity smearing
 routine_name(wtgs2) = 'MPI sweep'   ! MPI sweep for gravity smearing
 routine_name(wtout) = 'Output'      ! output
 routine_name(wtsho) = 'Shockfind'   ! shockfind
 routine_name(wtrad) = 'Radiation'   ! radiation
 routine_name(wtmir) = 'MICCG solve' ! MICCG solver for radiation
 routine_name(wtmra) = 'MICCG A'     ! MICCG A matrix
 routine_name(wtper) = 'PETSc solve' ! PETSc solver for radiation
 routine_name(wtprv) = 'PETSc vec'   ! PETSc vector assembly
 routine_name(wtpra) = 'PETSc A'     ! PETSc A matrix
 routine_name(wtprc) = 'PETSc Coeff' ! PETSc A matrix coefficients
 routine_name(wtprm) = 'PETSc MPI'   ! PETSc A matrix MPI
 routine_name(wtopc) = 'Opacity'     ! opacity
 routine_name(wtrfl) = 'Rad flux'    ! radiative flux
 routine_name(wtsnk) = 'Sink motion' ! sink motion
 routine_name(wtacc) = 'Accretion'   ! sink motion
 routine_name(wtmpi) = 'MPI exchange'! MPI exchange
 routine_name(wtwai) = 'MPI wait'    ! MPI wait
 routine_name(wttot) = 'Total'       ! total

 clock_on(0:n_wt) = .false.

 return
end subroutine init_profiler

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE GET_MAXLBL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute maximum number of characters in the label

subroutine get_maxlbl

 integer:: i

!-----------------------------------------------------------------------------

 maxlbl = 0
 do i = 0, n_wt
  if(wtime_max(i)>thres) &
   maxlbl = max(maxlbl,len_trim(routine_name(i))+get_layer(i)*2)
 end do
 maxlbl = maxlbl + 1

return
end subroutine get_maxlbl

! Get layer number \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function get_layer(i) result(ilayer)
 integer,intent(in)::i
 integer:: j,ilayer
 j = i
 ilayer = 0
 if(i>n_wt)return
 do
  j = parent(j)
  if(j>0)then
   ilayer = ilayer + 1
  else
   exit
  end if
 end do
end function get_layer

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE START_CLOCK
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To start recording wall time for a given category

subroutine start_clock(i)

 use omp_lib

 integer,intent(in)::i
 integer::j

!-----------------------------------------------------------------------------

 if(clock_on(i))then
  print*,"Error in start_clock"
  print*,'category: ',routine_name(i)
  stop
 end if
 wtime(i) = wtime(i) - omp_get_wtime()
 clock_on(i) = .true.

! Start clock for all parent categories too
 j=i
 do
  j = parent(j)
  if(j<0)exit
  if(.not.clock_on(j))then
   wtime(j) = wtime(j) - omp_get_wtime()
   clock_on(j) = .true.
  else
   exit
  end if
 end do

return
end subroutine start_clock

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE STOP_CLOCK
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To stop recording wall time for a given category

subroutine stop_clock(i)

 use omp_lib

 integer,intent(in)::i

!-----------------------------------------------------------------------------

 if(.not.clock_on(i))then
  print*,"Error in stop_clock"
  print*,'category: ',routine_name(i)
  stop
 end if

 wtime(i) = wtime(i) + omp_get_wtime()
 clock_on(i) = .false.

return
end subroutine stop_clock

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                            SUBROUTINE RESET_CLOCK
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To reset the clock for a given category

subroutine reset_clock(i)

 use omp_lib

 integer,intent(in)::i
 integer:: n,j

!-----------------------------------------------------------------------------

 wtime(i) = 0d0
 clock_on(i) = .false.

! Reset clock for all child categories too
 all_loop:do n = 1, n_wt
  j = parent(n)
  if(j/=i)then
   find_child:do
    if(j==i)then
     exit find_child
    elseif(j==wttot)then
     cycle all_loop
    end if
    j = parent(j)
   end do find_child
  end if
  wtime(n) = 0d0
  clock_on(n) = .false.
 end do all_loop

return
end subroutine reset_clock

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE REDUCE_CLOCKS_MPI
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To reduce profiling clocks across MPI tasks
subroutine reduce_clocks_mpi
  use mpi_utils, only:allreduce_mpi,nprocs

  integer::wti

  ! The average walltime for each category gives the most meaningful
  ! diagnostic, because it indicates where time is spent. Simply taking the
  ! max may not be useful, because a task that finishes a step quickly will
  ! record less time in that step, but more time in the subsequent step,
  ! meaning the times may add up to more than the total time.

  do wti = 0, n_wt
    wtime_avg(wti) = wtime(wti)
    wtime_max(wti) = wtime(wti)
    wtime_min(wti) = wtime(wti)

    call allreduce_mpi('sum', wtime_avg(wti))
    call allreduce_mpi('max', wtime_max(wti))
    call allreduce_mpi('min', wtime_min(wti))

    wtime_avg(wti) = wtime_avg(wti) / real(nprocs)

    ! Imbalance is the difference between the maximum and minimum wall time
    ! as a fraction of the average wall time
    if(wtime_avg(wti) > 0) then
      imbalance(wti) = (wtime_max(wti) - wtime_min(wti)) / wtime_avg(wti)
    else
      imbalance(wti) = 0
    end if

  end do

end subroutine reduce_clocks_mpi

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE PROFILER_OUTPUT1
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output one line for profiler output

subroutine profiler_output1(ui,wti)

 integer,intent(in):: ui,wti
 character(len=30):: form1,lbl
 integer::i,j,k,l

!-----------------------------------------------------------------------------

 if(wtime_max(wti)<=thres)return

 lbl = trim(routine_name(wti))
 j = get_layer(wti)
 i = maxlbl
 if(j>0)then

  i = i + 4
  if(last_in_category(wti))then
   lbl = "└─"//trim(lbl)
  else
   lbl = "├─"//trim(lbl)
  end if

  if(j>1)then
   l = wti
   do k = j, 2, -1
    l = parent(l)
    if(last_in_category(l))then
     lbl = "  "//trim(lbl)
    else
     lbl = "│ "//trim(lbl)
     i = i+2
    end if
   end do
  end if

 end if

 write(form1,'(a,i2,a)')'(a',i,',":",4(1X,F10.3,a))'
 write(ui,form1)adjustl(lbl),wtime_avg(wti),'s',&
                wtime_avg(wti)/wtime_avg(wtlop)*1d2,'%',&
                wtime_avg(wti)/wtime_avg(wttot)*1d2,'%',&
                imbalance(wti)*1d2,'%'

return
end subroutine profiler_output1

function last_in_category(wti)
! Function to judge whether a profiler element is the last element in a subcategory.
 integer,intent(in):: wti
 logical:: last_in_category
 integer:: k, next

! Find the next category with the same parent and finite wall time
 next = n_wt+1
 if(wti<n_wt)then
  do k = wti+1, n_wt
   if(parent(k)==parent(wti).and.wtime_max(k)>thres)then
    next = k
    exit
   end if
  end do
 end if

 if(next==n_wt+1)then
  last_in_category = .true.
 else
  last_in_category = .false.
 end if

end function last_in_category

end module profiler_mod
