module profiler_mod
 implicit none

 public:: init_profiler,profiler_output1,start_clock,stop_clock,reset_clock
 integer,parameter:: n_wt=22 ! number of profiling categories
 real(8):: wtime(0:n_wt)
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
  wttim=12,& ! time stepping
  wtgrv=13,& ! gravity
  wtgbn=14,& ! gravbound
  wtpoi=15,& ! Poisson solver
  wthyp=16,& ! hyperbolic self-gravity
  wtsho=17,& ! shockfind
  wtrad=18,& ! radiation
  wtopc=19,& ! opacity
  wtrfl=20,& ! radiative flux
  wtsnk=21,& ! sink particles
  wtout=22,& ! output
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
 parent(wttim) = wtlop ! time stepping
 parent(wtgrv) = wtlop ! gravity
 parent(wtgbn) = wtgrv ! gravbound
 parent(wtpoi) = wtgrv ! Poisson solver
 parent(wthyp) = wtgrv ! hyperbolic self-gravity
 parent(wtout) = wtlop ! output
 parent(wtsho) = wtlop ! shockfind
 parent(wtrad) = wtlop ! radiation
 parent(wtopc) = wtrad ! opacity
 parent(wtrfl) = wtrad ! radiative flux
 parent(wtsnk) = wtlop ! sink particles
 parent(wttot) =-1     ! Total
 
! Make sure to keep routine name short
 routine_name(wtini) = 'Setup'       ! initial conditions
 routine_name(wtgri) = '1st_Gravity' ! initial gravity
 routine_name(wtou1) = '1st_Output'  ! initial output
 routine_name(wtlop) = 'Main_loop'   ! main loop
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
 routine_name(wtpoi) = 'MICCG'       ! Poisson solver
 routine_name(wthyp) = 'Hyperbolic'  ! hyperbolic self-gravity
 routine_name(wtout) = 'Output'      ! output
 routine_name(wtsho) = 'Shockfind'   ! shockfind
 routine_name(wtrad) = 'Radiation'   ! radiation
 routine_name(wtopc) = 'Opacity'     ! opacity
 routine_name(wtrfl) = 'Rad_flux'    ! radiative flux
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
