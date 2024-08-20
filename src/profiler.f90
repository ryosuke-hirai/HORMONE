module profiler_mod
 implicit none

 public:: init_profiler,profiler_output1,start_clock,stop_clock,reset_clock
 integer,parameter:: n_wt=27 ! number of profiling categories
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
  wttim=13,& ! time stepping
  wtgrv=14,& ! gravity
  wtgbn=15,& ! gravbound
  wtpoi=16,& ! Poisson solver
  wthyp=17,& ! hyperbolic self-gravity
  wtgsm=18,& ! gravity smearing
  wtsho=19,& ! shockfind
  wtrad=20,& ! radiation
  wtopc=21,& ! opacity
  wtrfl=22,& ! radiative flux
  wtsnk=23,& ! sink motion
  wtacc=24,& ! sink accretion
  wtout=25,& ! output
  wtmpi=26,& ! mpi exchange
  wtwai=27,& ! mpi wait
  wttot=0    ! total
 integer,public:: parent(0:n_wt),maxlbl
 character(len=30),public:: routine_name(0:n_wt)
 logical,public:: clock_on(0:n_wt)

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                       SUBROUTINE INIT_PROFILER
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To initialize variables for profiler

subroutine init_profiler

 integer::i

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
 parent(wttim) = wtlop ! time stepping
 parent(wtgrv) = wtlop ! gravity
 parent(wtgbn) = wtgrv ! gravbound
 parent(wtpoi) = wtgrv ! Poisson solver
 parent(wthyp) = wtgrv ! hyperbolic self-gravity
 parent(wtgsm) = wtgrv ! gravity smearing
 parent(wtout) = wtlop ! output
 parent(wtsho) = wtlop ! shockfind
 parent(wtrad) = wtlop ! radiation
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
 routine_name(wttim) = 'Timestep'    ! time stepping
 routine_name(wtgrv) = 'Gravity'     ! gravity
 routine_name(wtgbn) = 'Gravbound'   ! gravbound
 routine_name(wtpoi) = 'MICCG'       ! Poisson solver
 routine_name(wthyp) = 'Hyperbolic'  ! hyperbolic self-gravity
 routine_name(wtgsm) = 'Grav smear'  ! gravity smearing
 routine_name(wtout) = 'Output'      ! output
 routine_name(wtsho) = 'Shockfind'   ! shockfind
 routine_name(wtrad) = 'Radiation'   ! radiation
 routine_name(wtopc) = 'Opacity'     ! opacity
 routine_name(wtrfl) = 'Rad flux'    ! radiative flux
 routine_name(wtsnk) = 'Sink motion' ! sink motion
 routine_name(wtacc) = 'Accretion'   ! sink motion
 routine_name(wtmpi) = 'MPI exchange'! MPI exchange
 routine_name(wtwai) = 'MPI wait'    ! MPI wait
 routine_name(wttot) = 'Total'       ! total

 do i = 0, n_wt
  maxlbl = max(maxlbl,len_trim(routine_name(i))+get_layer(i)*2)
 end do
 maxlbl = maxlbl + 1

 clock_on(0:n_wt) = .false.

 return
end subroutine init_profiler

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
!                        SUBROUTINE START_CLOCK
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
!                          SUBROUTINE STOP_CLOCK
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
!                          SUBROUTINE RESET_CLOCK
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
!                          SUBROUTINE REDUCE_CLOCKS_MPI
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
!                      SUBROUTINE PROFILER_OUTPUT1
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output one line for profiler output

subroutine profiler_output1(ui,wti)

 integer,intent(in):: ui,wti
 character(len=30):: form1,lbl
 integer::j,k,next
 real(8),parameter:: thres=1d-10

!-----------------------------------------------------------------------------

 if(wtime_max(wti)<=thres)return

 lbl = trim(routine_name(wti))
 j = get_layer(wti)
 if(j>0)then

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

  if(next/=n_wt+1)then
   lbl = "├─"//trim(lbl)
  else
   lbl = "└─"//trim(lbl)
  end if
  if(j>1)then
   do k = j, 2, -1
    lbl = "│ "//trim(lbl)
   end do
  end if
 end if

 select case(get_layer(wti))
 case(0)
  k = maxlbl
 case(1:)
  k = maxlbl + 2 + get_layer(wti)*2
 end select
 write(form1,'(a,i2,a)')'(a',k,',":",4(1X,F10.3,a))'
 write(ui,form1)adjustl(lbl),wtime_avg(wti),'s',&
                wtime_avg(wti)/wtime_avg(wtlop)*1d2,'%',&
                wtime_avg(wti)/wtime_avg(wttot)*1d2,'%',&
                imbalance(wti)*1d2,'%'

return
end subroutine profiler_output1



end module profiler_mod
