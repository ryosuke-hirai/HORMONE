module gravity_miccg_mod

 implicit none
 private

 public :: gravity_miccg, setup_grvcg

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
  ! call finalise_petsc()
  ! call setup_gravity_matrix_petsc  ! all matrix entries written to A_petsc
  ! call solve_system_petsc(cg, cgsrc, x)
  ! call finalise_petsc()
  stop "finished petsc solver"
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
!                        SUBROUTINE SETUP_GRVCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the A matrix elements for Laplacian

subroutine setup_grvcg(is,ie,js,je,ks,ke,cg)

 use settings,only:crdnt,eq_sym,gbtype
 use grid,only:gis_global,gie_global,gje_global,gks_global,gke_global,&
               xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,dx3,dxi3,idx3,&
               sini,sinc,rdis
 use miccg_mod,only:cg_set,ijk_from_l,get_preconditioner

 integer,intent(in)::is,ie,js,je,ks,ke
 type(cg_set),intent(out):: cg
 integer:: in,jn,kn,lmax,dim,i,j,k,l

!-----------------------------------------------------------------------------

 cg%is=is;cg%ie=ie; cg%js=js;cg%je=je; cg%ks=ks;cg%ke=ke
 in = ie-is+1; jn = je-js+1; kn = ke-ks+1
 cg%in=in; cg%jn=jn; cg%kn=kn
 lmax = in*jn*kn; cg%lmax=lmax
 if(ie>is.and.je>js.and.ke>ks)then
  dim=3
 elseif(ie>is.and.(je>js.or.ke>ks))then
  dim=2
 elseif(ie>is.and.je==js.and.ke==ks.and.crdnt==2)then
  dim=1
 else
  print*,'Error in setup_grvcg, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(1) ! 1D
! 1D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==2)then
   cg%Adiags = 2
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,i,j,k)
    cg%A(1,l) = -( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) )
    cg%A(2,l) = xi1(i)**2/dx1(i+1)
    if(i==gis_global.and.xi1s>0d0)& ! for inner boundary
     cg%A(1,l)=cg%A(1,l)+xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(gbtype==1.and.i==gie_global)& ! for Robin boundary at bc1o
     cg%A(1,l) = cg%A(1,l) + cg%A(2,l)*x1(i)/x1(i+1)
    if(i==ie)cg%A(2,l) = 0d0

   end do

! Pre-conditioner matrix with MICCG(1,1) method
   cg%cdiags = 2
   allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
   cg%ic(1) = 0
   cg%ic(2) = 1
   cg%alpha = 0.99d0

  end if

 case(2) ! 2D
! 2D cylindrical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==1.and.je==js.and.ke>ks)then

   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,i,j,k)
    cg%A(1,l) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
                      + 2d0*x1(i)*dxi1(i) / (dx3(k)*dx3(k+1)) ) &
                    * 0.5d0*sum(dx3(k:k+1))
    cg%A(2,l) = 0.5d0*xi1(i)*sum(dx3(k:k+1))/dx1(i+1)
    cg%A(3,l) = x1(i)*dxi1(i)/dx3(k+1)
    if(gbtype==1)then
     if(i==ie)cg%A(1,l) = cg%A(1,l) + cg%A(2,l)*rdis(i,k)/rdis(i+1,k)
     if(k==ke)cg%A(1,l) = cg%A(1,l) + cg%A(2,l)*rdis(i,k)/rdis(i,k+1)
    end if
    if(i==gie_global)cg%A(2,l) = 0d0
    if(k==gke_global)cg%A(3,l) = 0d0

    if(eq_sym.and.k==gks_global)& ! for Neumann boundary at bc3i
     cg%A(1,l) = cg%A(1,l) + x1(i)*dxi1(i)/dx3(k+1)

   end do

! 2D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.ie>is.and.je>js.and.ke==ks)then

   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,i,j,k)

    cg%A(1,l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                  * sinc(j)*dxi2(j) &
               + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    cg%A(2,l) = xi1(i)**2 *sinc(j)*dxi2(j)/dx1(i+1)
    cg%A(3,l) = sini(j)*dxi1(i)/dx2(j+1)
    if(i==gis_global.and.xi1s>0d0)& ! inner boundary
     cg%A(1,l)=cg%A(1,l)+xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(gbtype==1.and.i==gie_global)& ! for Robin boundary at bc1o
     cg%A(1,l) = cg%A(1,l) + cg%A(2,l)*x1(i)/x1(i+1)
    if(i==gie_global)cg%A(2,l) = 0d0
    if(j==gje_global)then ! Reflection boundary at axis or equatorial plane
     cg%A(1,l) = cg%A(1,l) + sini(j)*dxi1(i)/dx2(j+1)
     cg%A(3,l) = 0d0
    end if

   end do

  else

   write(*,'(8(a,i0))') 'Error (setup_grvcg): dimension is not supported; dim = ',dim,' with crdnt = ',crdnt,' is = ',is,' ie = ',ie,' js = ',js,' je = ',je,' ks = ',ks,' ke = ',ke
   error stop

  end if

! Pre-conditioner matrix with MICCG(1,2) method
  cg%cdiags = 4
  allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
  cg%ic(1) = 0
  cg%ic(2) = 1
  cg%ic(3) = in-1
  cg%ic(4) = in
  cg%alpha = 0.99d0


 case(3) ! 3D
! 3D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==2)then

   cg%Adiags = 5
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in
   cg%ia(4) = in*jn
   cg%ia(5) = in*jn*(kn-1)

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,i,j,k)
    cg%A(1,l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                  * sinc(j)*dxi2(j)*dxi3(k) &
                  + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i)*dxi3(k) &
                  + ( idx3(k)+idx3(k-1) )*dxi1(i)*dxi2(j) )
    cg%A(2,l) = xi1(i)**2 *sinc(j)*dxi2(j)*dxi3(k)/dx1(i+1)
    cg%A(3,l) = sini(j)*dxi1(i)*dxi3(k)/dx2(j+1)
    cg%A(4,l) = dxi1(i)*dxi2(j)/dx3(k+1)
    cg%A(5,l) = dxi1(i)*dxi2(j)/dx3(k  )
    if(eq_sym.and.j==gje_global)& ! for Neumann boundary at bc2o
     cg%A(1,l) = cg%A(1,l) + cg%A(3,l)
    if(i==gis_global.and.xi1s>0d0)& ! for inner boundary
     cg%A(1,l) = cg%A(1,l) + xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(gbtype==1.and.i==gie_global)& ! for Robin boundary at bc1o
     cg%A(1,l) = cg%A(1,l) + cg%A(2,l)*x1(i)/x1(i+1)
    if(i==gie_global)cg%A(2,l) = 0d0
    if(j==gje_global)cg%A(3,l) = 0d0
    if(k==gke_global)cg%A(4,l) = 0d0
    if(k/=gks_global)cg%A(5,l) = 0d0

   end do

! Pre-conditioner matrix with ICCG(1,1) method
   cg%cdiags = 5
   allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
   cg%ic(1) = 0
   cg%ic(2) = 1
   cg%ic(3) = in
   cg%ic(4) = in*jn
   cg%ic(5) = in*jn*(kn-1)
   cg%alpha = 0.999d0 ! No modification

  else
   write(*,'(A,I0,A,I0)') 'Error (setup_grvcg): dimension is not supported; dim = ',dim,' with crdnt = ',crdnt
   error stop

  end if

 end select

 call get_preconditioner(cg)

end subroutine setup_grvcg

end module gravity_miccg_mod
