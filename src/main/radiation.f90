module radiation_mod

 use miccg_mod,only:cg_set
 use opacity_mod
 use radiation_utils

 implicit none

 public :: radiation,radiation_setup,radiative_force,rad_heat_cool
 private:: setup_radcg,get_diffusion_coeff,get_radA,get_radb
 real(8),allocatable,private:: rsrc(:),radK(:,:,:),gradE(:,:,:,:)

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE RADIATION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute radiation transport with flux limited diffusion

subroutine radiation

 use grid
 use physval
 use miccg_mod,only:cg=>cg_rad,miccg,ijk_from_l
 use profiler_mod

 integer:: l,i,j,k
 real(8),allocatable:: x(:)

!-----------------------------------------------------------------------------

 call start_clock(wtrad)

 call rad_boundary

! Advection and radiative acceleration terms are updated in hydro step

! Update heating/cooling term first
 call rad_heat_cool

! Then update the diffusion term
 call get_diffusion_coeff(cg) ! use erad^n for diffusion coefficients
 call get_radA(cg)
 call get_radb(cg)

 allocate( x(1:cg%lmax) )
 !$omp parallel do private(l,i,j,k)
 do l = 1, cg%lmax
  call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
  x(l) = erad(i,j,k)
 end do
!$omp end parallel do

 call miccg(cg,rsrc,x) ! returns erad^{n+1}

! update erad and u
!$omp parallel do private(l,i,j,k)
 do l = 1, cg%lmax
  call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
  erad(i,j,k)   = x(l)
  u(i,j,k,irad) = x(l)
 end do
!$omp end parallel do

 call stop_clock(wtrad)

return
end subroutine radiation


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GET_DIFFUSION_COEFF
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gradE on each cell

subroutine get_diffusion_coeff(cg)

 use constants,only:clight
 use utils,only:get_grad
 use physval,only:erad,d,T,erad

 type(cg_set),intent(in)::cg
 integer:: i,j,k
 real(8):: RR,ll,kappar

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,RR,ll,kappar) collapse(3)
 do k = cg%ks, cg%ke
  do j = cg%js, cg%je
   do i = cg%is, cg%ie
    gradE(1:3,i,j,k) = get_grad(erad,i,j,k)
    kappar = kappa_r(d(i,j,k),T(i,j,k))
    RR = norm2(gradE(1:3,i,j,k)) / (d(i,j,k)*kappar*erad(i,j,k))
    ll = lambda(RR)
    radK(i,j,k) = clight*ll/(d(i,j,k)*kappar)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_diffusion_coeff

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_RADA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Compute divergence term (A) for radiation diffusion equation

subroutine get_radA(cg)

 use utils,only:har_mean
 use grid,only:dt,dvol
 use miccg_mod,only:ijk_from_l,get_preconditioner

 type(cg_set),intent(inout)::cg
 integer:: i,j,k,l,dim

!-----------------------------------------------------------------------------

 if(cg%in>1.and.cg%jn>1.and.cg%kn>1)then
  dim=3
 elseif(cg%in>1.and.cg%jn>1.and.cg%kn==1)then
  dim=22
 elseif(cg%in>1.and.cg%jn==1.and.cg%kn>1)then
  dim=23
 elseif(cg%in>1.and.cg%jn==1.and.cg%kn==1)then
  dim=11
 elseif(cg%in==1.and.cg%jn>1.and.cg%kn==1)then
  dim=12
 elseif(cg%in==1.and.cg%jn==1.and.cg%kn>1)then
  dim=13
 else
  print*,'Error in get_radA, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(11) ! 1D x-direction

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(i>cg%is)& ! gradient term in x- direction
    cg%A(1,l) = cg%A(1,l) + geo(1,i-1,j,k)*har_mean(radK(i-1:i,j,k))
   if(i<cg%ie)& ! gradient term in x+ direction
    cg%A(1,l) = cg%A(1,l) + geo(1,i  ,j,k)*har_mean(radK(i:i+1,j,k))

   cg%A(2,l) = -geo(1,i,j,k)*har_mean(radK(i:i+1,j,k))
   if(i==cg%ie)cg%A(2,l) = 0d0
  end do
!$omp end parallel do

 case(12) ! 1D y-direction

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(j>cg%js)& ! gradient term in y- direction
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j-1,k)*har_mean(radK(i,j-1:j,k))
   if(j<cg%je)& ! gradient term in y+ direction
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j  ,k)*har_mean(radK(i,j:j+1,k))

   cg%A(2,l) = -geo(2,i,j,k)*har_mean(radK(i,j:j+1,k))
   if(j==cg%je)cg%A(2,l) = 0d0
  end do
!$omp end parallel do

 case(13) ! 1D z-direction

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(k>cg%ks)& ! gradient term in z- direction
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j,k-1)*har_mean(radK(i,j,k-1:k))
   if(k<cg%ke)& ! gradient term in z+ direction
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j,k  )*har_mean(radK(i,j,k:k+1))

   cg%A(2,l) = -geo(3,i,j,k)*har_mean(radK(i,j,k:k+1))
   if(k==cg%ke)cg%A(2,l) = 0d0
  end do
!$omp end parallel do

 case(22) ! 2D xy-plane

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(i>cg%is)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i-1,j,k)*har_mean(radK(i-1:i,j,k))
   if(i<cg%ie)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i  ,j,k)*har_mean(radK(i:i+1,j,k))
   if(j>cg%js)&
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j-1,k)*har_mean(radK(i,j-1:j,k))
   if(j<cg%je)&
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j  ,k)*har_mean(radK(i,j:j+1,k))

   cg%A(2,l) = -geo(1,i,j,k)*har_mean(radK(i:i+1,j,k))
   cg%A(3,l) = -geo(2,i,j,k)*har_mean(radK(i,j:j+1,k))
   if(i==cg%ie)cg%A(2,l) = 0d0
   if(j==cg%je)cg%A(3,l) = 0d0
  end do
!$omp end parallel do

 case(23) ! xz-plane

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(i>cg%is)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i-1,j,k)*har_mean(radK(i-1:i,j,k))
   if(i<cg%ie)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i  ,j,k)*har_mean(radK(i:i+1,j,k))
   if(k>cg%ks)&
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j-1,k)*har_mean(radK(i,j,k-1:k))
   if(k<cg%ke)&
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j  ,k)*har_mean(radK(i,j,k:k+1))

   cg%A(2,l) = -geo(1,i,j,k)*har_mean(radK(i:i+1,j,k))
   cg%A(3,l) = -geo(3,i,j,k)*har_mean(radK(i,j,k:k+1))
   if(i==cg%ie)cg%A(2,l) = 0d0
   if(k==cg%ke)cg%A(3,l) = 0d0
  end do
!$omp end parallel do

 case(3) ! 3D

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)

   cg%A(1,l) = dvol(i,j,k)/dt
   if(i>cg%is)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i-1,j,k)*har_mean(radK(i-1:i,j,k))
   if(i<cg%ie)&
    cg%A(1,l) = cg%A(1,l) + geo(1,i  ,j,k)*har_mean(radK(i:i+1,j,k))
   if(j>cg%js)&
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j-1,k)*har_mean(radK(i,j-1:j,k))
   if(j<cg%je)&
    cg%A(1,l) = cg%A(1,l) + geo(2,i,j  ,k)*har_mean(radK(i,j:j+1,k))
   if(k>cg%ks)&
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j-1,k)*har_mean(radK(i,j,k-1:k))
   if(k<cg%ke)&
    cg%A(1,l) = cg%A(1,l) + geo(3,i,j  ,k)*har_mean(radK(i,j,k:k+1))

   cg%A(2,l) = -geo(1,i,j,k)*har_mean(radK(i:i+1,j,k))
   cg%A(3,l) = -geo(2,i,j,k)*har_mean(radK(i,j:j+1,k))
   cg%A(4,l) = -geo(3,i,j,k  )*har_mean(radK(i,j,k:k+1))
   cg%A(5,l) = -geo(3,i,j,k-1)*har_mean(radK(i,j,k-1:k))
   if(i==cg%ie)cg%A(2,l) = 0d0
   if(j==cg%je)cg%A(3,l) = 0d0
   if(k==cg%ke)cg%A(4,l) = 0d0
   if(k/=cg%ks)cg%A(5,l) = 0d0
  end do
!$omp end parallel do

 end select

 call get_preconditioner(cg)

return
end subroutine get_radA

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE GET_RADB
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get source term (b) for the radiation diffusion equation

subroutine get_radb(cg)

 use grid,only:dvol,dt
 use physval,only:erad
 use miccg_mod,only:cg_set,ijk_from_l

 type(cg_set),intent(inout):: cg
 integer:: i,j,k,l

!-----------------------------------------------------------------------------

!$omp parallel do private(l,i,j,k)
  do l = 1, cg%lmax
   call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
! Note: Dirichlet boundary conditions should be included here.
   rsrc(l) = erad(i,j,k)*dvol(i,j,k)/dt
  end do
!$omp end parallel do

return
end subroutine get_radb

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SETUP_RADCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the necessary parameters for the CG method

subroutine setup_radcg(is,ie,js,je,ks,ke,cg)

 use settings,only:crdnt

 integer,intent(in)::is,ie,js,je,ks,ke
 type(cg_set),intent(out):: cg
 integer:: in,jn,kn,lmax,dim

!-----------------------------------------------------------------------------

 cg%is=is;cg%ie=ie; cg%js=js;cg%je=je; cg%ks=ks;cg%ke=ke
 in = ie-is+1; jn = je-js+1; kn = ke-ks+1
 cg%in=in; cg%jn=jn; cg%kn=kn
 lmax = in*jn*kn; cg%lmax=lmax
 if(ie>is.and.je>js.and.ke>ks)then
  dim=3
 elseif(ie>is.and.(je>js.or.ke>ks))then
  dim=2
 elseif(ie>is.and.je==js.and.ke==ks)then
  dim=1
 elseif(ie==is.and.je>js.and.ke==ks)then
  dim=1
 elseif(ie==is.and.je==js.and.ke>ks)then
  dim=1
 else
  print*,'Error in setup_radcg, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(1) ! 1D
! 1D coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cg%Adiags = 2
  allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
  cg%ia(1) = 0
  cg%ia(2) = 1

! Pre-conditioner matrix with MICCG(1,1) method
  cg%cdiags = 2
  allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
  cg%ic(1) = 0
  cg%ic(2) = 1
  cg%alpha = 0.99d0

 case(2) ! 2D
! 2D cylindrical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==1.and.je==js.and.ke>ks)then
   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

! 2D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.ie>is.and.je>js.and.ke==ks)then

   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

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

! Pre-conditioner matrix with ICCG(1,1) method
   cg%cdiags = 5
   allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
   cg%ic(1) = 0
   cg%ic(2) = 1
   cg%ic(3) = in
   cg%ic(4) = in*jn
   cg%ic(5) = in*jn*(kn-1)
   cg%alpha = 0.999d0 ! No modification

  end if

 end select

return
end subroutine setup_radcg


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE RADIATION_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up radiation parameters

subroutine radiation_setup

 use settings,only:radswitch
 use grid,only:is,ie,js,je,ks,ke
 use miccg_mod,only:cg=>cg_rad

!-----------------------------------------------------------------------------

 if(radswitch==1)then
  call setup_radcg(is,ie,js,je,ks,ke,cg)
  call get_geo
  allocate(radK(is-1:ie+1,js-1:je+1,ks-1:ke+1),gradE(1:3,is:ie,js:je,ks:ke),&
           rsrc(1:cg%lmax) )
 end if

return
end subroutine radiation_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE RADIATIVE_FORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute radiative acceleration and associated source terms

subroutine radiative_force

 use utils,only:get_grad
 use grid
 use physval

 integer:: i,j,k,l,m
 real(8):: RR,ll,ff,vdotfrad,Pedd(1:3,1:3),radwork,kappar
 real(8),dimension(1:3):: frad,gradv1,gradv2,gradv3,nn

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,RR,ll,ff,frad,vdotfrad,gradv1,gradv2,gradv3,&
!$omp nn,l,m,radwork,Pedd) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    kappar = kappa_r(d(i,j,k),T(i,j,k))
    if(kappar<=tiny(kappar))cycle
    RR = norm2(gradE(1:3,i,j,k)) / (d(i,j,k)*kappar*erad(i,j,k))
    ll = lambda(RR)
    ff = ll + (ll*RR)**2

    frad(1:3) = -ll*gradE(1:3,i,j,k)
    vdotfrad = frad(1)*v1(i,j,k) + frad(2)*v2(i,j,k) + frad(3)*v3(i,j,k)

    gradv1 = get_grad(v1,i,j,k)
    gradv2 = get_grad(v2,i,j,k)
    gradv3 = get_grad(v3,i,j,k)

    nn(1:3) = gradE(1:3,i,j,k)/max(norm2(gradE(1:3,i,j,k)),epsilon(erad(i,j,k)))
    do m = 1, 3
     do l = 1, 3
      if(l>m)cycle
      Pedd(l,m) = (3d0*ff-1d0)*nn(l)*nn(m)
      if(l==m)Pedd(l,m) = Pedd(l,m) + (1d0-ff)
      Pedd(l,m) = 0.5d0*Pedd(l,m)
     end do
    end do
    radwork = Pedd(1,1)*gradv1(1) + Pedd(1,2)*gradv1(2) + Pedd(1,3)*gradv1(3) &
            + Pedd(1,2)*gradv2(1) + Pedd(2,2)*gradv2(2) + Pedd(2,3)*gradv2(3) &
            + Pedd(1,3)*gradv3(1) + Pedd(2,3)*gradv3(2) + Pedd(3,3)*gradv3(3)

    src(i,j,k,imo1) = src(i,j,k,imo1) + frad(1)
    src(i,j,k,imo2) = src(i,j,k,imo2) + frad(2)
    src(i,j,k,imo3) = src(i,j,k,imo3) + frad(3)
    src(i,j,k,iene) = src(i,j,k,iene) + vdotfrad
    src(i,j,k,irad) = -radwork

   end do
  end do
 end do
!$omp end parallel do

return
end subroutine radiative_force

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE RAD_HEAT_COOL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute the radiative heating and cooling
!          Follows the method of Turner & Stone 2001

subroutine rad_heat_cool

 use constants,only:clight,sigma,fac_egas
 use grid
 use physval
 use pressure_mod,only:get_etot_from_eint,getT_from_de,Trad

 integer:: i,j,k
 real(8):: a1,a2,c1,c2,kappap,eint1

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,a1,a2,c1,c2,kappap,eint1) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    if(kappap<=0d0)cycle

    a1 = 4d0*kappap*sigma/(fac_egas*imu(i,j,k))**4/d(i,j,k)**3*dt
    a2 = clight*kappap*d(i,j,k)*dt
    c1 = (1d0+a2)/a1
    c2 = -c1*eint(i,j,k)-a2/a1*erad(i,j,k)

    eint1 = max(eint(i,j,k),erad(i,j,k))

    call solve_quartic(c1,c2,eint1)
    erad(i,j,k) = (a1*eint1**4+erad(i,j,k))/(1d0+a2)
    eint(i,j,k) = eint1
    e(i,j,k) = get_etot_from_eint(i,j,k)
    call getT_from_de(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))

    u(i,j,k,iene) = e(i,j,k)
    u(i,j,k,irad) = erad(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine rad_heat_cool

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE RAD_BOUNDARY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for erad

subroutine rad_boundary

!-----------------------------------------------------------------------------

! Not needed for now

return
end subroutine rad_boundary

end module radiation_mod
