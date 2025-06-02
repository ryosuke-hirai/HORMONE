module gravity_elliptic_mod

 implicit none
 private

 public :: gravity_elliptic

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVITY_ELLIPTIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gravitational forces using the sparse matrix method.

subroutine gravity_elliptic

 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use miccg_mod,only:miccg
 use timestep_mod,only:timestep
 use profiler_mod
 use matrix_utils,only:l_from_ijk,ijk_from_l
 use matrix_solver,only:solve_system_grv
 integer:: i,j,k,l,lmax
 real(8),allocatable,dimension(:):: x, cgsrc

!-------------------------------------------------------------------------

 call start_clock(wtelg)

 lmax = (gie-gis+1)*(gje-gjs+1)*(gke-gks+1)

 allocate( x(1:lmax), cgsrc(1:lmax) )

! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(gbtype==0)call gravbound

! cylindrical (equatorial+axial symmetry) ####################################
 if(je==js.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l)
  do l = 1, lmax
   call ijk_from_l(l,gis,gjs,gks,gie-gis+1,gje-gjs+1,i,j,k)
   x(l) = grvphi(i,j,k)
   cgsrc(l) = 0d0
   if(i>=is)then;if(i<=ie)then;if(k>=ks)then;if(k<=ke)then
    cgsrc(l) = 4d0*pi*G*gsrc(i,j,k)*x1(i)*dxi1(i)*sum(dx3(k:k+1))*0.5d0
   end if;end if;end if;end if
   if(gbtype==0)then
    if(k==gks) cgsrc(l) = cgsrc(l) - x1 (i)*dxi1(i)*idx3(k  )*phi3i(i,k-1)
    if(i==gie) cgsrc(l) = cgsrc(l) - xi1(i)*dxi3(k)*idx1(i+1)*phi1o(i+1,k)
    if(k==gke) cgsrc(l) = cgsrc(l) - x1 (i)*dxi1(i)*idx3(k+1)*phi3o(i,k+1)
   end if
  end do
!$omp end parallel do

! spherical (axial symmetry) #################################################
 elseif(crdnt==2.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l)
  do l=1,lmax
   call ijk_from_l(l,gis,gjs,gks,gie-gis+1,gje-gjs+1,i,j,k)
   x(l) = grvphi(i,j,k)
   cgsrc(l) = 0d0
   if(i>=is)then;if(i<=ie)then;if(j>=js)then;if(j<=je)then
    cgsrc(l) = 4d0*pi*G*gsrc(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)
   end if;end if;end if;end if
   if(gbtype==0)then
    if(i==gie) cgsrc(l) = cgsrc(l) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)&
                                    *phiio(i+1,j)
    if(i==gis) cgsrc(l) = cgsrc(l) - xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)&
                                    *phiii(i-1,j)
   end if
  end do
!$omp end parallel do

! spherical (3D) #############################################################
 elseif(crdnt==2.and.dim==3)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l)
  do l=1,lmax
   call ijk_from_l(l,gis,gjs,gks,gie-gis+1,gje-gjs+1,i,j,k)
   x(l) = grvphi(i,j,k)
   cgsrc(l) = 0d0
   if(i>=is)then;if(i<=ie)then
    cgsrc(l) = 4d0*pi*G*gsrc(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)*dxi3(k)
   end if;end if

   if(gbtype==0)then
    if(i==gie)cgsrc(l)= cgsrc(l) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)*dxi3(k)&
                                  *phiio(i+1,j)
    if(i==gis)cgsrc(l)= cgsrc(l) - xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)*dxi3(k)&
                                  *phiii(i-1,j)
   end if
  end do
!$omp end parallel do
 end if

! ############################################################################

! Solve Poisson equation with the conjugate gradient method
  call solve_system_grv(cgsrc, x)

!-------------------------------------------------------------------------

! convert x to phi
!$omp parallel do private(i,j,k,l) collapse(3)
 do k = gks, gke
  do j = gjs, gje
   do i = gis, gie
    l = l_from_ijk(i,j,k,gis,gjs,gks,gie-gis+1,gje-gjs+1)
    grvphi(i,j,k) = x(l)
   end do
  end do
 end do
!$omp end parallel do


 if(je==js.and.crdnt==1.and.dim==2)then ! for cylindrical coordinates

  do k = gks,gke
   do j = js,je
    grvphi(is-2,j,k) = grvphi(is+1,j,k)
    grvphi(is-1,j,k) = grvphi(is  ,j,k)
   end do
  end do
  grvphi(gis:gie,js:je,gks-2) = grvphi(gis:gie,js:je,gks)
  grvphi(gis:gie,js:je,gks-1) = grvphi(gis:gie,js:je,gks)

 elseif(ke==ks.and.crdnt==2.and.dim==2)then ! for spherical coordinates (2D)

  grvphi(is-1:gie+1,js-2,ks) = grvphi(is-1:gie+1,js+1,ks)
  grvphi(is-1:gie+1,js-1,ks) = grvphi(is-1:gie+1,js,ks)
  grvphi(is-1:gie+1,je+1,ks) = grvphi(is-1:gie+1,je,ks)
  grvphi(is-1:gie+1,je+2,ks) = grvphi(is-1:gie+1,je-1,ks)

  grvphi(is-1,:,:) = grvphi(is,:,:)
  grvphi(is-2,:,:) = grvphi(is+1,:,:)
  grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                  ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)

 elseif(crdnt==2.and.dim==3)then ! for spherical coordinates (3D)

  grvphi(is-1:gie+1,js-2,ks:ke) = grvphi(is-1:gie+1,js,ks:ke)
  grvphi(is-1:gie+1,js-1,ks:ke) = grvphi(is-1:gie+1,js,ks:ke)
  grvphi(is-1:gie+1,je+1,ks:ke) = grvphi(is-1:gie+1,je,ks:ke)
  grvphi(is-1:gie+1,je+2,ks:ke) = grvphi(is-1:gie+1,je,ks:ke)

  grvphi(is-1,:,:) = grvphi(is,:,:)
  grvphi(is-2,:,:) = grvphi(is+1,:,:)
  grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                  ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)

 end if

 if(gravswitch==3)then
  call timestep
  cgrav_old = cgrav
 end if

 call stop_clock(wtelg)
end subroutine gravity_elliptic

end module gravity_elliptic_mod
