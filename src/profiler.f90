module profiler_mod
 implicit none

 public:: init_profiler,profiler_output1,start_clock,stop_clock,reset_clock
 integer,parameter:: n_wt=19 ! number of profiling categories
 real(8):: wtime(0:n_wt)
 integer,parameter:: &
  wtini=1 ,& ! initial conditions
  wtou1=2 ,& ! initial output
  wtlop=3 ,& ! main loop
  wthyd=4 ,& ! hydrodynamics
  wtflx=5 ,& ! numerical flux
  wtrng=6 ,& ! Runge-Kutta
  wtbnd=7 ,& ! boundary conditions
  wtsrc=8 ,& ! source terms
  wtint=9 ,& ! MUSCL interpolation
  wteos=10,& ! equation of state
  wttim=11,& ! time stepping
  wtgrv=12,& ! gravity
  wtgbn=13,& ! gravbound
  wtsho=14,& ! shockfind
  wtrad=15,& ! radiation
  wtopc=16,& ! opacity
  wtrfl=17,& ! radiative flux
  wtsnk=18,& ! sink particles
  wtout=19,& ! output
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
 parent(wtou1) = wtini ! initial output
 parent(wtlop) = wttot ! main loop
 parent(wthyd) = wtlop ! hydrodynamics
 parent(wtflx) = wthyd ! numerical flux
 parent(wtrng) = wthyd ! Runge-Kutta
 parent(wtbnd) = wthyd ! boundary conditions
 parent(wtsrc) = wthyd ! source terms
 parent(wtint) = wthyd ! MUSCL interpolation
 parent(wteos) = wthyd ! equation of state
 parent(wttim) = wthyd ! time stepping
 parent(wtgrv) = wtlop ! gravity
 parent(wtgbn) = wtgrv ! gravbound
 parent(wtout) = wtlop ! output
 parent(wtsho) = wtlop ! shockfind
 parent(wtrad) = wtlop ! radiation
 parent(wtopc) = wtrad ! opacity
 parent(wtrfl) = wtrad ! radiative flux
 parent(wtsnk) = wtlop ! sink particles
 parent(wttot) =-1     ! Total
 
! Make sure to keep routine name short
 routine_name(wtini) = 'Setup'       ! initial conditions
 routine_name(wtou1) = '1st Output'  ! initial output
 routine_name(wtlop) = 'Main loop'   ! main loop
 routine_name(wthyd) = 'Hydro'       ! hydrodynamics
 routine_name(wtflx) = 'Numflux'     ! numerical flux
 routine_name(wtrng) = 'RungeKutta'  ! Runge-Kutta
 routine_name(wtbnd) = 'Boundary'    ! boundary conditions
 routine_name(wtsrc) = 'Source'      ! source terms
 routine_name(wtint) = 'Interpolate' ! MUSCL interpolation
 routine_name(wteos) = 'EoS'         ! equation of state
 routine_name(wttim) = 'Timestep'    ! time stepping
 routine_name(wtgrv) = 'Gravity'     ! gravity
 routine_name(wtgbn) = 'Gravbound'   ! gravbound
 routine_name(wtout) = 'Output'      ! output
 routine_name(wtsho) = 'Shockfind'   ! shockfind
 routine_name(wtrad) = 'Radiation'   ! radiation
 routine_name(wtopc) = 'Opacity'     ! opacity
 routine_name(wtrfl) = 'Rad flux'    ! radiative flux
 routine_name(wtsnk) = 'Sinks'       ! sink particles
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
!                        SUBROUTINE STOP_CLOCK
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
!                        SUBROUTINE RESET_CLOCK
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To reset the clock for a given category

subroutine reset_clock(i)

 use omp_lib

 integer,intent(in)::i

!-----------------------------------------------------------------------------

 wtime(i) = 0d0
 clock_on(i) = .false.

return
end subroutine reset_clock

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

 if(wtime(wti)<=thres)return

 lbl = trim(routine_name(wti))
 j = get_layer(wti)
 if(j>0)then

! Find the next category with the same parent and finite wall time
  next = n_wt+1
  if(wti<n_wt)then
   do k = wti+1, n_wt
    if(parent(k)==parent(wti).and.wtime(k)>thres)then
     next = k
     exit
    end if
   end do
  end if

  if(next/=n_wt+1)then
   lbl = "├─"//lbl
  else
   lbl = "└─"//lbl
  end if
  if(j>1)then
   do k = j, 2, -1
    lbl = "│ "//lbl
   end do
  end if
 end if

 select case(get_layer(wti))
 case(0)
  k = maxlbl
 case(1:)
  k = maxlbl + 2 + get_layer(wti)*2
 end select
 write(form1,'(a,i2,a)')'(a',k,',":",3(1X,F10.3,a))'
 write(ui,form1)adjustl(lbl),wtime(wti),'s',&
                wtime(wti)/wtime(wtlop)*1d2,'%',&
                wtime(wti)/wtime(wttot)*1d2,'%'

return
end subroutine profiler_output1



end module profiler_mod
