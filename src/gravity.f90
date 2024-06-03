module gravity_mod
 implicit none

 public:: gravity,gravsetup
 private:: setup_grvcg

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate gravitational forces.

subroutine gravity

 use settings,only:eq_sym,grav_init_relax,grav_init_other
 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use miccg_mod,only:cg=>cg_grv,miccg,l_from_ijk,ijk_from_l
 use timestep_mod,only:timestep
 use profiler_mod

 integer:: i,j,k,n,l, gin, gjn, gkn, tngrav
 real(8):: phih, h
 real(8),allocatable,dimension(:):: x, cgsrc

!-----------------------------------------------------------------------------

 if(gravswitch==0.or.gravswitch==1)return

 call start_clock(wtgrv)

 gin = gie - gis + 1
 gjn = gje - gjs + 1
 gkn = gke - gks + 1

! Set source term for gravity
 call get_gsrc(gsrc)

 if(grav_init_relax .and. tn==0) then
  call gravity_relax

 elseif(gravswitch==2 .or. (gravswitch==3 .and. tn==0))then
  allocate( x(1:cg%lmax), cgsrc(1:cg%lmax) )

  if(grav_init_other.and.gravswitch==3)return
! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if(gbtype==0)call gravbound

  call start_clock(wtpoi)

! cylindrical (equatorial+axial symmetry) ####################################
  if(je==js.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l)
   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
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
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
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
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
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
  call miccg(cg,cgsrc,x)

!-------------------------------------------------------------------------

! convert x to phi
!$omp parallel do private(i,j,k,l) collapse(3)
  do k = gks, gke
   do j = gjs, gje
    do i = gis, gie
     l = l_from_ijk(i,j,k,gis,gjs,gks,gie-gis+1,gje-gjs+1,gke-gks+1)
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

 endif


if(gravswitch==3 .and. tn/=0)then
! Hyperbolic Self-Gravity $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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

end if
call stop_clock(wtgrv)

end subroutine gravity

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

 real(8),intent(in):: dtg,cgrav_now,cgrav_old
 integer:: i,j,k,n, jb, kb, grungen
 real(8):: faco, facn, fact, vol

!-----------------------------------------------------------------------------

!$omp parallel
 do grungen = 1, grktype
!$omp single
  call get_runge_coeff(grungen,grktype,faco,fact,facn)
!$omp end single

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

! Smear gravity in central regions -----------------------------------
  if(fmr_max>0)then
   do n = 1, fmr_max
    if(fmr_lvl(n)==0)cycle
    jb=min(2**(fmr_max-n+1),je)-1 ; kb=min(2**(fmr_max-n+1),ke)-1
    if(n==1)then
     jb=je-js; kb=ke-ks
    end if
!$omp do private(i,j,k,vol) collapse(3)
    do k = ks, ke, kb+1
     do j = js, je, jb+1
      do i = is+sum(fmr_lvl(0:n-1)), is+sum(fmr_lvl(0:n))-1
       vol = sum(dvol(i,j:j+jb,k:k+kb))
       grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb) &
                                      *dvol(i,j:j+jb,k:k+kb)) / vol
       grvpsi(i,j:j+jb,k:k+kb) = sum(grvpsi(i,j:j+jb,k:k+kb) &
                                      *dvol(i,j:j+jb,k:k+kb)) / vol
      end do
     end do
    end do
!$omp end do
   end do
  end if
! --------------------------------------------------------------------

 end do

! Damping of grvpsi by separation of variables
!$omp do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    grvpsi(i,j,k) = grvpsi(i,j,k)*exp(-(1d0-cgrav_old/cgrav_now+alphagrv*hgcfl))
   end do
  end do
 end do
!$omp end do
!$omp end parallel

return
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
 integer:: i,j,k
 !$omp parallel
 select case(crdnt)
  ! Cylindrical coordinates ++++++++++++++++++++++++++++++++++++++++++++++++++++
 case(1)
  !$omp do private(j,k) collapse(2)
  do k = ks, ke
   do j = js, je
    grvphi(is-1,j,k) = grvphi(is,j,k)
    if(gbtype==1)grvphi(gie+1,j,k) = grvphi(gie,j,k) &
    *orgdis(gie,j,k)/orgdis(gie+1,j,k)
   end do
  end do
  !$omp end do
  if(je>js)then
   !$omp do private(i,k) collapse(2)
   do k = ks, ke
    do i = is, ie
     grvphi(i,js-1,k) = grvphi(i,je,k)
     grvphi(i,je+1,k) = grvphi(i,js,k)
    end do
   end do
   !$omp end do
  end if
  !$omp do private(i,j) collapse(2)
  do j = js, je
   do i = is, ie
    if(eq_sym)grvphi(i,j,ks-1) = grvphi(i,j,ks)
    if(gbtype==1)then
     grvphi(i,j,ks-1) = grvphi(i,j,ks)*orgdis(i,j,ks)/orgdis(i,j,ks-1)
     grvphi(i,j,ke+1) = grvphi(i,j,ke)*orgdis(i,j,ke)/orgdis(i,j,ke+1)
    end if
   end do
  end do
  !$omp end do

  ! Spherical coordinates ++++++++++++++++++++++++++++++++++++++++++++++++++++++
 case(2)
  !$omp do private(j,k) collapse(2)
  do k = ks, ke
   do j = js, je
    grvphi(is-1,j,k) = grvphi(is,j,k)
    if(gbtype==1)grvphi(ie+1,j,k) = grvphi(ie,j,k)*x1(ie)/x1(ie+1)
   end do
  end do
  !$omp end do
  !$omp do private(i,k) collapse(2)
  do k = ks, ke
   do i = is, ie
    grvphi(i,js-1,k) = grvphi(i,js,k)
    grvphi(i,je+1,k) = grvphi(i,je,k)
   end do
  end do
  !$omp end do
  if(ke>ks)then
   !$omp do private(i,j) collapse(2)
   do j = js, je
    do i = is, ie
     grvphi(i,j,ks-1) = grvphi(i,j,ke)
     grvphi(i,j,ke+1) = grvphi(i,j,ks)
    end do
   end do
   !$omp end do
  end if

 end select
 ! ============================================================================

 !$omp end parallel

end subroutine hg_boundary_conditions

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SETUP_GRVCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the A matrix elements for Laplacian

subroutine setup_grvcg(is,ie,js,je,ks,ke,cg)

 use settings,only:crdnt,eq_sym,gbtype
 use grid,only:gis_global,gie_global,&
               gjs_global,gje_global,&
               gks_global,gke_global,&
               xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,idx2,dx3,dxi3,idx3,&
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
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
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
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
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
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)

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
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
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

return
end subroutine setup_grvcg

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRAVSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the matrix for Poisson equation

subroutine gravsetup

 use settings,only:courant
 use constants,only:huge
 use utils,only:isequal
 use grid
 use gravmod
 use miccg_mod,only:cg_grv

 integer:: i,j,k
 real(8):: h

!-----------------------------------------------------------------------------

! set initial x0
 if(tn==0 .and. isequal(maxval(grvphi(gis:gie,gjs:gje,gks:gke)), 0d0))&
  grvphi = -1d3

 if(tn==0) cgrav_old = cgrav

 if(gravswitch==2.or.gravswitch==3)then
  call setup_grvcg(gis,gie,gjs,gje,gks,gke,cg_grv)
 end if

! For Hyperbolic gravity solver ----------------------------------------------
 if(gravswitch==3)then

  allocate(lap_coeff(-1:3,gis-1:gie,gjs-1:gje,gks-1:gke))

! for axisymmetric cylindrical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(je==js.and.dim==2.and.crdnt==1)then
! Cartoon mesh method
   j = js
   do k = gks-1, gke
    do i = gis-1, gie

    end do
   end do
   print*,'Cylindrical coordinates not supported yet for hyperbolic self-gravity'
   stop

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
                  + ( idx3(k)+idx3(k-1) )*dxi1(i)*dxi2(j) )
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
!                          SUBROUTINE GET_LAPPHI
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
     if(j==gjs_global.and.bc2is==0)jm=je_global
     if(j==gje_global.and.bc2os==0)jp=js_global
     lap = lap + lap_coeff(2,i,j-1,k)*grvphi(i,jm,k)
     lap = lap + lap_coeff(2,i,j  ,k)*grvphi(i,jp,k)
    end if
    if(ke>ks)then
     km = k-1 ; kp = k+1
     if(k==gks_global.and.bc3is==0)km=ke_global
     if(k==gke_global.and.bc3os==0)kp=ks_global
     lap = lap + lap_coeff(3,i,j,k-1)*grvphi(i,j,km)
     lap = lap + lap_coeff(3,i,j,k  )*grvphi(i,j,kp)
    end if
    lapphi(i,j,k) = lap/lap_coeff(-1,i,j,k)

    hgsrc(i,j,k) = 4d0*pi*G*gsrc(i,j,k)
   end do
  end do
 end do
!$omp end do

return
end subroutine get_lapphi_hgsrc

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

 real(8):: err_grvphi
 integer :: grktype_org
 integer, parameter :: maxiter = 1000000
 real(8), parameter :: itertol = 1d-9
 integer :: i

 call timestep
 cgrav_old = cgrav

 grktype_org = grktype
 grktype = 1

! Repeat the damped hyperbolic step until the solution is stationary
 do i = 1, maxiter
   call hyperbolic_gravity_step(cgrav,cgrav_old,dtgrav)
   err_grvphi = maxval(abs( (grvphi(is:ie,js:je,ks:ke)-grvphiorg(is:ie,js:je,ks:ke,1))/grvphi(is:ie,js:je,ks:ke) ))
   if (i > 2 .and. err_grvphi < itertol) exit
 enddo

 grktype = grktype_org

return
end subroutine gravity_relax

end module gravity_mod
