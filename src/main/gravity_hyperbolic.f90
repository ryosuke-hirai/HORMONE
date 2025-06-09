module gravity_hyperbolic_mod

 implicit none
 private

 public :: gravity_hyperbolic, setup_grv_hyperbolic, hg_boundary_conditions, hyperbolic_gravity_step

contains


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                     SUBROUTINE GRAVITY_HYPERBOLIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Wrapper routine for the hyperbolic self-gravity solver

subroutine gravity_hyperbolic

 use grid
 use gravmod
 use timestep_mod,only:timestep
 use profiler_mod
 use gravbound_mod,only:gravbound

 integer:: l,tngrav

!-----------------------------------------------------------------------------

 call start_clock(wthyp)

 tngrav = ceiling((time+dt-grvtime)/dtgrav)
 dtgrav = (time+dt-grvtime)/dble(tngrav)

 if(gbtype==0)call gravbound

 do l = 1, tngrav
  call hyperbolic_gravity_step(cgrav,cgrav_old,dtgrav)

  grvtime = grvtime + dtgrav
  cgrav_old = cgrav

 end do

 call hg_boundary_conditions

 call stop_clock(wthyp)

end subroutine gravity_hyperbolic


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                   SUBROUTINE HYPERBOLIC_GRAVITY_STEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: One time step for the hyperbolic self-gravity solver

subroutine hyperbolic_gravity_step(cgrav_now,cgrav_old,dtg)

 use settings,only:grktype,alphagrv
 use grid
 use gravmod,only:grvphiorg,grvphi,grvpsi,lapphi,hgsrc,gsrc,hgcfl
 use rungekutta_mod,only:get_runge_coeff
 use smear_mod,only:smear
 use mpi_domain,only:exchange_gravity_mpi

 real(8),intent(in):: dtg,cgrav_now,cgrav_old
 integer:: i,j,k,grungen
 real(8):: faco, facn, fact, damp_factor

!-----------------------------------------------------------------------------

 do grungen = 1, grktype

! Perform MPI neighbour exchange
  call exchange_gravity_mpi
  call get_runge_coeff(grungen,grktype,faco,fact,facn)

!$omp parallel

! First set flux and source term
  call get_lapphi_hgsrc(grvphi,gsrc,lapphi,hgsrc)

!$omp do private (i,j,k) collapse(3)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     if(grungen==1)then
      grvphiorg(i,j,k,1) = grvphi(i,j,k)
      grvphiorg(i,j,k,2) = grvpsi(i,j,k)
     end if
     grvphi(i,j,k) = faco*grvphiorg(i,j,k,1) + facn*grvphi(i,j,k) &
                   + fact*dtg*cgrav_now*grvpsi(i,j,k)
     grvpsi(i,j,k) = faco*grvphiorg(i,j,k,2) + facn*grvpsi(i,j,k) &
                   + fact*dtg*cgrav_now*(lapphi(i,j,k) - hgsrc(i,j,k))
    end do
   end do
  end do
!$omp end do
!$omp end parallel

! Smear gravity in central regions --------------------------------------------

  if(fmr_max>0)call smear('grav')

! -----------------------------------------------------------------------------

 end do

! Damping of grvpsi by separation of variables
 damp_factor = exp(-(1d0-cgrav_old/cgrav_now+alphagrv*hgcfl))
!$omp parallel do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    grvpsi(i,j,k) = grvpsi(i,j,k)*damp_factor
   end do
  end do
 end do
!$omp end parallel do

end subroutine hyperbolic_gravity_step


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                   SUBROUTINE HG_BOUNDARY_CONDITIONS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Boundary conditions for the hyperbolic self-gravity solver

subroutine hg_boundary_conditions

 use settings
 use grid
 use gravmod
 use mpi_domain,only:exchange_gravity_mpi

 integer:: i,j,k

!-----------------------------------------------------------------------------

 call exchange_gravity_mpi

 !$omp parallel
 select case(crdnt)
! Cylindrical coordinates ++++++++++++++++++++++++++++++++++++++++++++++++++++
 case(1)
  !$omp do private(j,k) collapse(2)
  do k = ks, ke
   do j = js, je
    if(is==is_global) grvphi(is-1,j,k) = grvphi(is,j,k)
    if(ie==ie_global .and. gbtype==1) grvphi(gie+1,j,k) = grvphi(gie,j,k) &
    *orgdis(gie,j,k)/orgdis(gie+1,j,k)
   end do
  end do
  !$omp end do
  if(je>js)then
   !$omp do private(i,k) collapse(2)
   do k = ks, ke
    do i = is, ie
     if(js==js_global) grvphi(i,js-1,k) = grvphi(i,je,k)
     if(je==je_global) grvphi(i,je+1,k) = grvphi(i,js,k)
    end do
   end do
   !$omp end do
  end if
  !$omp do private(i,j) collapse(2)
  do j = js, je
   do i = is, ie
    if(eq_sym)grvphi(i,j,ks-1) = grvphi(i,j,ks)
    if(gbtype==1)then
     if(ks==ks_global) grvphi(i,j,ks-1) = grvphi(i,j,ks)*orgdis(i,j,ks)/orgdis(i,j,ks-1)
     if(ke==ke_global) grvphi(i,j,ke+1) = grvphi(i,j,ke)*orgdis(i,j,ke)/orgdis(i,j,ke+1)
    end if
   end do
  end do
  !$omp end do

  ! Spherical coordinates ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 case(2)
  !$omp do private(j,k) collapse(2)
  do k = ks, ke
   do j = js, je
    if(is==is_global) grvphi(is-1,j,k) = grvphi(is,j,k)
    if(ie==ie_global .and. gbtype==1) grvphi(ie+1,j,k) = grvphi(ie,j,k)*x1(ie)/x1(ie+1)
   end do
  end do
  !$omp end do
  !$omp do private(i,k) collapse(2)
  do k = ks, ke
   do i = is-1, ie+1
    if(js==js_global) grvphi(i,js-1,k) = grvphi(i,js,k)
    if(je==je_global) grvphi(i,je+1,k) = grvphi(i,je,k)
   end do
  end do
  !$omp end do
  if(ke>ks)then
   ! Note: in MPI mode, cell exchange already implements periodic BCs
   if(ks==ks_global .and. ke==ke_global) then
    !$omp do private(i,j) collapse(2)
    do j = js-1, je+1
     do i = is-1, ie+1
      grvphi(i,j,ks-1) = grvphi(i,j,ke)
      grvphi(i,j,ke+1) = grvphi(i,j,ks)
     end do
    end do
    !$omp end do
   end if
  end if

 end select
 !$omp end parallel

end subroutine hg_boundary_conditions


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE GET_LAPPHI_HGSRC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute Laplacian phi

subroutine get_lapphi_hgsrc(grvphi,gsrc,lapphi,hgsrc)

 use settings,only:bc2is,bc2os,bc3is,bc3os
 use constants,only:pi,G
 use grid
 use gravmod,only:lap_coeff

 real(8),allocatable,dimension(:,:,:),intent(in):: grvphi,gsrc
 real(8),allocatable,dimension(:,:,:),intent(inout):: lapphi,hgsrc
 integer:: i,j,k,jm,jp,km,kp
 real(8):: lap

!-----------------------------------------------------------------------------

!$omp do private(i,j,k,jm,jp,km,kp,lap) collapse(3)
 do k = gks,gke
  do j = gjs,gje
   do i = gis,gie
    lap = lap_coeff(0,i,j,k)*grvphi(i,j,k)

    lap = lap + lap_coeff(1,i-1,j,k)*grvphi(i-1,j,k)
    lap = lap + lap_coeff(1,i  ,j,k)*grvphi(i+1,j,k)
    if(je>js)then
     jm = j-1 ; jp = j+1
     ! If multiple MPI tasks in this dimension, periodic BCs are already handled
     if (js==js_global.and.je==je_global) then
      if(j==gjs_global.and.bc2is==0)jm=je_global
      if(j==gje_global.and.bc2os==0)jp=js_global
     endif
     lap = lap + lap_coeff(2,i,j-1,k)*grvphi(i,jm,k)
     lap = lap + lap_coeff(2,i,j  ,k)*grvphi(i,jp,k)
    end if
    if(ke>ks)then
     km = k-1 ; kp = k+1
     ! If multiple MPI tasks in this dimension, periodic BCs are already handled
     if (ks==ks_global.and.ke==ke_global) then
      if(k==gks_global.and.bc3is==0)km=ke_global
      if(k==gke_global.and.bc3os==0)kp=ks_global
     endif
     lap = lap + lap_coeff(3,i,j,k-1)*grvphi(i,j,km)
     lap = lap + lap_coeff(3,i,j,k  )*grvphi(i,j,kp)
    end if
    lapphi(i,j,k) = lap/lap_coeff(-1,i,j,k)

    hgsrc(i,j,k) = 4d0*pi*G*gsrc(i,j,k)
   end do
  end do
 end do
!$omp end do

end subroutine get_lapphi_hgsrc


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SETUP_GRV_HYPERBOLIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the coefficients for the hyperbolic self-gravity solver

subroutine setup_grv_hyperbolic

 use settings,only:eq_sym
 use utils,only:isequal
 use grid
 use gravmod

 integer:: i,j,k

!-----------------------------------------------------------------------------

  allocate(lap_coeff(-1:3,gis-1:gie,gjs-1:gje,gks-1:gke))

! for axisymmetric cylindrical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(je==js.and.dim==2.and.crdnt==1)then
! TODO: To be replaced with Cartoon mesh method
   j = js
   do k = gks-1, gke
    do i = gis-1, gie
     lap_coeff(0,i,j,k) = -( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
                            +2d0*x1(i)*dxi1(i)/(dx3(k)*dx3(k+1)) ) &
                           * 0.5d0*(dx3(k)+dx3(k+1))
     lap_coeff(1,i,j,k) = xi1(i)*0.5d0*(dx3(k)+dx3(k+1))/dx1(i+1)
     lap_coeff(2,i,j,k) = 0d0
     lap_coeff(3,i,j,k) = x1(i)*dxi1(i)/dx3(k+1)
     if(gbtype==1)then
      if(i==gie_global)then ! for Robin boundary at bc1o
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                          + lap_coeff(1,i,j,k) *rdis(i,k)/rdis(i+1,k)
       lap_coeff(1,i,j,k) = 0d0
      end if
      if(k==gks_global)then ! for Robin boundary at bc3i
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                          + x1(i)*dxi1(i)/dx3(k) *rdis(i,k)/rdis(i,k-1)
      end if
      if(k==gks_global-1)then
       lap_coeff(3,i,j,k) = 0d0
      end if
      if(k==gke_global)then ! for Robin boundary at bc3o
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                          + lap_coeff(3,i,j,k) *rdis(i,k)/rdis(i,k+1)
       lap_coeff(3,i,j,k) = 0d0
      end if
     end if
     if(eq_sym)then
      if(k==gks_global-1)& ! Neumann boundary at bc3i
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) + lap_coeff(3,i,j,k)
      if(k==gks_global-1)&
       lap_coeff(3,i,j,k) = 0d0
     end if
! This normalization is applied after dot(lap,phi) due to stability reasons
     lap_coeff(-1,i,j,k) = x1(i)*dxi1(i)*0.5d0*(dx3(k)+dx3(k+1))
    end do
   end do

! for axisymmetrical spherical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.ke==ks.and.dim==2)then

   k = ks
!$omp parallel do private(i,j) schedule(static) collapse(2)
   do j = gjs-1, gje
    do i = gis-1, gie
     lap_coeff(0,i,j,k) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                           * sinc(j)*dxi2(j) &
                        + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
     lap_coeff(1,i,j,k) = xi1(i)**2 *sinc(j)*dxi2(j)/dx1(i+1)
     lap_coeff(2,i,j,k) = sini(j)*dxi1(i)/dx2(j+1)
     lap_coeff(3,i,j,k) = 0d0
     if(i==gis_global.and.xi1s>0d0)& ! inner boundary
      lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                         + xi1(i-1)**2*sinc(j)*dxi2(j)/dx1(i)
     if(gbtype==1.and.i==gie_global)then ! for Robin boundary at bc1o
      lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                         + lap_coeff(1,i,j,k) *x1(i)/x1(i+1)
      lap_coeff(1,i,j,k) = 0d0
     end if
     if(j==gje_global)then ! Neumann boundary at bc2o
      lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) + lap_coeff(2,i,j,k)
      lap_coeff(2,i,j,k) = 0d0
     end if
! This normalization is applied after dot(lap,phi) due to stability reasons
     lap_coeff(-1,i,j,k) = x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)
    end do
   end do
!$omp end parallel do


! for 3D spherical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.dim==3)then

!$omp parallel do private(i,j,k) collapse(3)
   do k = gks-1, gke
    do j = gjs-1, gje
     do i = gis-1, gie
      lap_coeff(0,i,j,k) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                  * sinc(j)*dxi2(j)*dxi3(k) &
                  + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i)*dxi3(k) &
                  + ( idx3(k)+idx3(k+1) )*dxi1(i)*dxi2(j) )
      lap_coeff(1,i,j,k) = xi1(i)**2 *sinc(j)*dxi2(j)*dxi3(k)/dx1(i+1)
      lap_coeff(2,i,j,k) = sini(j)*dxi1(i)*dxi3(k)/dx2(j+1)
      lap_coeff(3,i,j,k) = dxi1(i)*dxi2(j)/dx3(k+1)

      if(j==gje_global)then ! for Neumann boundary at bc2o
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) + lap_coeff(2,i,j,k)
       lap_coeff(2,i,j,k) = 0d0
      end if
      if(i==gis_global.and.xi1s>0d0)& ! for inner boundary
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                          + xi1(i-1)**2*sinc(j)*dxi2(j)*dxi3(k)/dx1(i)
      if(gbtype==1.and.i==gie_global)then ! for Robin boundary at bc1o
       lap_coeff(0,i,j,k) = lap_coeff(0,i,j,k) &
                          + lap_coeff(1,i,j,k)*x1(i)/x1(i+1)
       lap_coeff(1,i,j,k) = 0d0
      end if
! This normalization is applied after dot(lap,phi) due to stability reasons
      lap_coeff(-1,i,j,k) = x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)*dxi3(k)

     end do
    end do
   end do
!$omp end parallel do

  else
   print *,'coefficients for Hyperbolic gravity not supported for this coordinate system'
   print*,'crdnt=',crdnt,', dim=',dim
   stop
  end if

  hgsrc = 0d0

end subroutine setup_grv_hyperbolic

end module gravity_hyperbolic_mod
