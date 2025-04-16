module gravity_miccg_mod

 implicit none
 private

 public :: gravity_miccg
 public :: compute_coeffs

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVITY_MICCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gravitational forces using the MICCG method

subroutine gravity_miccg

 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use miccg_mod,only:cg=>cg_grv,miccg,l_from_ijk,ijk_from_l
 use timestep_mod,only:timestep
 use profiler_mod
#ifdef USE_PETSC
 use petsc_solver_mod,only:init_petsc,finalise_petsc
#endif
 integer:: i,j,k,l
 real(8),allocatable,dimension(:):: x, cgsrc

!-------------------------------------------------------------------------

 allocate( x(1:cg%lmax), cgsrc(1:cg%lmax) )

! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(gbtype==0)call gravbound

 call start_clock(wtpoi)

! cylindrical (equatorial+axial symmetry) ####################################
 if(je==js.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
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
  do l=1,cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
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
  do l=1,cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
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
#ifdef USE_PETSC
  !--- PETSc solver branch ---
  ! call init_petsc()
  ! call solve_system_petsc(cg, cgsrc, x)
  call finalise_petsc()
  stop "petsc solver not implemented"
#else
  !--- Original MICCG solver branch ---
  call miccg(cg, cgsrc, x)
#endif

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

 call stop_clock(wtpoi)
end subroutine gravity_miccg


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE COMPUTE_COEFFS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Computes the stencil coefficients for a given grid point.

subroutine compute_coeffs(dim, i, j, k, ie, ke, coeffs)
  use settings, only: crdnt, eq_sym, gbtype
  use grid,only:gis_global,gie_global,gje_global,gks_global,gke_global,&
                xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,dx3,dxi3,idx3,&
                sini,sinc,rdis

  integer, intent(in) :: dim, i, j, k, ie, ke
  real(8), intent(out), dimension(:) :: coeffs
  real(8) :: sum_dx3

  select case(dim)
  case(1)  ! 1D spherical coordinates (crdnt == 2 assumed)
      ! Expect coeffs(1:2)
      coeffs(1) = - ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) )
      coeffs(2) = xi1(i)**2/dx1(i+1)
      if ( i == gis_global .and. xi1s > 0d0 ) then
        coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j) * dxi2(j) * idx1(i)
      end if
      if ( gbtype == 1 .and. i == gie_global ) then
        coeffs(1) = coeffs(1) + coeffs(2) * x1(i)/x1(i+1)
      end if
      if ( i == ie ) then ! TODO: check MPI
        coeffs(2) = 0d0
      end if

  case(2)  ! 2D
      if ( crdnt == 1 ) then
        ! 2D cylindrical coordinates
        ! Expect coeffs(1:3)
        ! Use sum_dx3 = dx3(k) + dx3(k+1)
        sum_dx3 = dx3(k) + dx3(k+1)
        coeffs(1) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) + &
                    2d0*x1(i)*dxi1(i)/( dx3(k)*dx3(k+1) ) ) * 0.5d0 * sum_dx3
        coeffs(2) = 0.5d0 * xi1(i) * sum_dx3 / dx1(i+1)
        coeffs(3) = x1(i) * dxi1(i) / dx3(k+1)
        if ( gbtype == 1 ) then
          if ( i == ie ) then ! TODO: check MPI
            coeffs(1) = coeffs(1) + coeffs(2) * rdis(i, k)/rdis(i+1, k)
          end if
          if ( k == ke ) then ! TODO: check MPI
            coeffs(1) = coeffs(1) + coeffs(2) * rdis(i, k)/rdis(i, k+1)
          end if
        end if
        if ( i == gie_global ) then
          coeffs(2) = 0d0
        end if
        if ( k == gke_global ) then
          coeffs(3) = 0d0
        end if
        if ( eq_sym .and. k == gks_global ) then
          coeffs(1) = coeffs(1) + x1(i)* dxi1(i) / dx3(k+1)
        end if
      else if ( crdnt == 2 ) then
        ! 2D spherical coordinates
        ! Expect coeffs(1:3)
        coeffs(1) = - ( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) * sinc(j)* dxi2(j) + &
                      ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j) ) * dxi1(i) )
        coeffs(2) = xi1(i)**2 * sinc(j)* dxi2(j) / dx1(i+1)
        coeffs(3) = sini(j) * dxi1(i) / dx2(j+1)
        if ( i == gis_global .and. xi1s > 0d0 ) then
          coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j)* dxi2(j)* idx1(i)
        end if
        if ( gbtype == 1 .and. i == gie_global ) then
          coeffs(1) = coeffs(1) + coeffs(2)* x1(i)/x1(i+1)
        end if
        if ( i == gie_global ) then
          coeffs(2) = 0d0
        end if
        if ( j == gje_global ) then
          coeffs(1) = coeffs(1) + sini(j) * dxi1(i) / dx2(j+1)
          coeffs(3) = 0d0
        end if
      end if

  case(3)  ! 3D spherical coordinates (crdnt == 2 assumed)
      ! Expect coeffs(1:5)
      coeffs(1) = - ( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) * sinc(j)* dxi2(j)* dxi3(k) + &
                    ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j) ) * dxi1(i)* dxi3(k) + &
                    ( idx3(k) + idx3(k-1) ) * dxi1(i)* dxi2(j) )
      coeffs(2) = xi1(i)**2 * sinc(j)* dxi2(j)* dxi3(k) / dx1(i+1)
      coeffs(3) = sini(j) * dxi1(i)* dxi3(k) / dx2(j+1)
      coeffs(4) = dxi1(i)* dxi2(j) / dx3(k+1)
      coeffs(5) = dxi1(i)* dxi2(j) / dx3(k)
      if ( eq_sym .and. j == gje_global ) then
        coeffs(1) = coeffs(1) + coeffs(3)
      end if
      if ( i == gis_global .and. xi1s > 0d0 ) then
        coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j)* dxi2(j)* idx1(i)
      end if
      if ( gbtype == 1 .and. i == gie_global ) then
        coeffs(1) = coeffs(1) + coeffs(2)* x1(i)/x1(i+1)
      end if
      if ( i == gie_global ) then
        coeffs(2) = 0d0
      end if
      if ( j == gje_global ) then
        coeffs(3) = 0d0
      end if
      if ( k == gke_global ) then
        coeffs(4) = 0d0
      end if
      if ( k /= gks_global ) then
        coeffs(5) = 0d0
      end if

  end select

end subroutine compute_coeffs

end module gravity_miccg_mod
