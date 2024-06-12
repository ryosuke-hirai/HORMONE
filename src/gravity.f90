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

 use settings,only:grav_init_relax
 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use gravity_hyperbolic_mod,only:gravity_hyperbolic,hg_boundary_conditions
 use gravity_miccg_mod,only:gravity_miccg
 use miccg_mod,only:miccg,l_from_ijk,ijk_from_l
 use timestep_mod,only:timestep
 use profiler_mod

 if(gravswitch==0.or.gravswitch==1)return

 call start_clock(wtgrv)

! Set source term for gravity
 call get_gsrc(gsrc)

 if(grav_init_relax .and. tn==0) then
  call gravity_relax
  call hg_boundary_conditions

 elseif(gravswitch==2 .or. (gravswitch==3 .and. tn==0))then
  call gravity_miccg
 endif

if(gravswitch==3 .and. tn/=0)then
 call gravity_hyperbolic
end if

call stop_clock(wtgrv)

end subroutine gravity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRAVSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the matrix for Poisson equation

subroutine gravsetup

 use utils,only:isequal
 use grid
 use gravmod
 use miccg_mod,only:cg_grv
 use gravity_hyperbolic_mod,only:setup_grv_hyperbolic
 use gravity_miccg_mod,only:setup_grvcg

! set initial x0
 if(tn==0 .and. isequal(maxval(grvphi(gis:gie,gjs:gje,gks:gke)), 0d0))&
  grvphi = -1d3

 if(tn==0) cgrav_old = cgrav

 if(gravswitch==2.or.gravswitch==3)then
  call setup_grvcg(gis,gie,gjs,gje,gks,gke,cg_grv)
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
 use timestep_mod
 use gravmod
 use mpi_utils,only:allreduce_mpi
 use mpi_domain,only:myrank
 use gravity_hyperbolic_mod,only:hyperbolic_gravity_step

 real(8):: err_grvphi
 integer :: grktype_org
 integer, parameter :: maxiter = 1000000
 real(8), parameter :: itertol = 1d-9
 integer :: i

 call timestep
 cgrav_old = cgrav

 grktype_org = grktype
 grktype = 1

 if (myrank==0) print*, 'Initialising gravity using hyperbolic solver...'

 ! Repeat the damped hyperbolic step until the solution is stationary
 do i = 1, maxiter
   call hyperbolic_gravity_step(cgrav,cgrav_old,dtgrav)
   err_grvphi = maxval(abs( (grvphi(is:ie,js:je,ks:ke)-grvphiorg(is:ie,js:je,ks:ke,1))/grvphi(is:ie,js:je,ks:ke) ))

   ! Get max error across all MPI tasks
   call allreduce_mpi('max', err_grvphi)

   ! Print error every 1000 steps
   if (mod(i,1000)==0 .and. myrank==0) then
    print*, 'iteration=', i, 'error=', err_grvphi
   endif

   if (i > 2 .and. err_grvphi < itertol) exit
 enddo

 if (myrank==0) print*, 'Gravity relaxation converged in ', i, ' iterations, with error', err_grvphi

 ! Reset grvpsi
 grvpsi = 0.d0

 ! Reset grktype
 grktype = grktype_org

return
end subroutine gravity_relax

end module gravity_mod
