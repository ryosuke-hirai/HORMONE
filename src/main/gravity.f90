module gravity_mod

 implicit none
 private

 public:: gravity,gravsetup

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate gravitational forces.

subroutine gravity

 use settings,only:grav_init_relax,grav_init_other,in_loop
 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use gravity_hyperbolic_mod,only:gravity_hyperbolic,hg_boundary_conditions
 use gravity_elliptic_mod,only:gravity_elliptic
 use timestep_mod,only:timestep
 use profiler_mod

 integer:: wtind

 if(gravswitch==0.or.gravswitch==1)return

 if(in_loop)then
  wtind = wtgrv
 else
  wtind = wtgri
 end if

 call start_clock(wtind)

! Set source term for gravity
 call get_gsrc(gsrc)

 if (tn==0) grvtime = 0.d0

 if(grav_init_other.and.tn==0)then
  call hg_boundary_conditions

 elseif(grav_init_relax .and. tn==0) then
  call gravity_relax
  call hg_boundary_conditions

 elseif(gravswitch==2 .or. (gravswitch==3 .and. tn==0))then
  call gravity_elliptic

 elseif(gravswitch==3 .and. tn/=0)then
  call gravity_hyperbolic
 end if

call stop_clock(wtind)

end subroutine gravity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRAVSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the matrix for Poisson equation

subroutine gravsetup
 use settings,only:grav_init_relax
 use utils,only:isequal
 use grid
 use gravmod
 use gravity_hyperbolic_mod,only:setup_grv_hyperbolic
 use matrix_vars,only:igrv
 use matrix_solver_mod,only:setup_matrix,write_A_grv

! set initial x0
 if(tn==0 .and. isequal(maxval(grvphi(gis:gie,gjs:gje,gks:gke)), 0d0))&
  grvphi = -1d3

 if(tn==0) cgrav_old = cgrav

 if(gravswitch==2.or.(gravswitch==3.and..not.grav_init_relax))then
  call setup_matrix(igrv)
  call write_A_grv ! writes to cg_grv or PETSc arrays
 end if

! For Hyperbolic gravity solver ----------------------------------------------
 if(gravswitch==3)then
  call setup_grv_hyperbolic
 end if

 return

end subroutine gravsetup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_GSRC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set gsrc for Poisson's equation

subroutine get_gsrc(gsrc)

 use grid,only:is,ie,js,je,ks,ke

 real(8),allocatable,intent(inout):: gsrc(:,:,:)
 integer:: i,j,k

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    gsrc(i,j,k) = get_gsrc1(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_gsrc

function get_gsrc1(i,j,k) result(gsrc)
 use settings,only:grvsrctype
 use physval,only:d,spc
 integer,intent(in)::i,j,k
 real(8):: gsrc
 select case(grvsrctype)
 case(0)
  gsrc = d(i,j,k) ! default
 case(1)
  gsrc = d(i,j,k)*spc(1,i,j,k) ! For others
 case default
  print*,'Error in grvsrctype, grvsrctype=',grvsrctype
  stop
 end select
end function get_gsrc1

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GRAVITY_RELAX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Use hyperbolic self-gravity as a relaxation method

subroutine gravity_relax

 use settings
 use grid
 use physval
 use timestep_mod
 use gravmod
 use mpi_utils,only:allreduce_mpi
 use mpi_domain,only:myrank
 use gravity_hyperbolic_mod,only:hyperbolic_gravity_step
 use profiler_mod

 real(8):: err
 integer :: grktype_org
 integer, parameter :: maxiter = 10000 ! TEMPORARY: May not be enough to fully converge, but set low for speed
 real(8), parameter :: itertol = 1d-2
 integer :: i

 real(8), allocatable :: mass(:,:,:)
 real(8) :: mtot

!-----------------------------------------------------------------------------

 if (myrank==0) print*, 'Initialising gravity using hyperbolic solver...'

 call timestep
 cgrav_old = cgrav

 grktype_org = grktype
 grktype = 1

 ! Mass of cells
 allocate(mass(is:ie,js:je,ks:ke))
 mass = d(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke)
 mtot = sum(mass)
 call allreduce_mpi('sum', mtot)

 err = 1.d0

 ! Repeat the damped hyperbolic step until the solution is stationary
 do i = 0, maxiter
  call hyperbolic_gravity_step(cgrav,cgrav_old,dtgrav)

  ! Only calculate error every 1000 iterations because it is expensive
  if (mod(i,1000)==0 .or. i==1) then
   ! Error in Poisson equation, weighted by mass
   err = sum( ((lapphi(is:ie,js:je,ks:ke) - hgsrc(is:ie,js:je,ks:ke)) &
             / hgsrc(is:ie,js:je,ks:ke))**2 * mass )
   call allreduce_mpi('sum', err)
   err = sqrt(err)/sqrt(mtot)
   if (myrank==0) print*, 'iteration=', i, 'error=', err

   ! Converged if below tolerance
   if (err < itertol .and. i>0) exit

  endif
 end do

 if (myrank==0) print*, 'Gravity relaxation converged in ', &
                        i, ' iterations, with error', err

 ! Reset grvpsi
 grvpsi = 0.d0

 ! Reset grktype
 grktype = grktype_org

! Reset clocks
 call reset_clock(wthyp)

return
end subroutine gravity_relax

end module gravity_mod
