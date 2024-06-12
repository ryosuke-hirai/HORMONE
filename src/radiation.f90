module radiation_mod
 implicit none

 public:: radiation,radiation_setup
 private::lambda,lambda_LP81,lambda_M78,lambda_K89,kappa_r,kappa_p,get_radflux,get_source_term,get_geo,setup_radcg
 real(8),allocatable,private:: geo(:,:,:,:)
 real(8),allocatable,public:: urad(:,:,:),rsrc(:,:,:)

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
 use miccg_mod,only:cg=>cg_rad,miccg,ijk_from_l,l_from_ijk
 use profiler_mod

 real(8),allocatable::radK(:,:,:,:)

!-----------------------------------------------------------------------------

 call start_clock(wtrad)

 allocate(radK(1:3,cg%is-1:cg%ie,cg%js-1:cg%je,cg%ks-1:cg%ke))

 call get_radflux(urad,d,T,cg,radK)
 call get_radA(radK,cg)
 call get_source_term(d,T,urad,dt,cg,rsrc)

! call miccg(cg,rsrc,x)

 call stop_clock(wtrad)

return
end subroutine radiation

! Get flux limiter for flux-limited diffusion
function lambda(R)
 use settings,only:lambdatype
 real(8),intent(in)::R
 real(8):: lambda
 select case(lambdatype)
 case(1) ! Levermore & Pomraning 1981
  lambda = lambda_LP81(R)
 case(2) ! Minerbo 1978
  lambda = lambda_M78(R)
 case(3) ! Kley 1989
  lambda = lambda_K89(R)
 case default
  print*,'Error from choice of radiation flux limiter'
  print*,'lambdatype=',lambdatype
  stop
 end select
end function lambda

! Levermore & Pomraning 1981
function lambda_LP81(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 lambda = (1d0/tanh(R)-1d0/R)/R
end function lambda_LP81

! Minerbo 1978
function lambda_M78(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 if(R<=1.5d0)then
  lambda = 2d0/(3d0+sqrt(9d0+12d0*R**2))
 else
  lambda = 1d0/(1d0+R+sqrt(1d0+2d0*R))
 end if
end function lambda_M78

! Kley 1989
function lambda_K89(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 if(R<=2d0)then
  lambda = 2d0/(3d0+sqrt(9d0+10d0*R**2))
 else
  lambda = 10d0/(10d0*R+9d0+sqrt(180d0*R+81d0))
 end if
end function lambda_K89

! Get Rosseland mean opacity
function kappa_r(d,T) result(kappa)
 real(8),intent(in)::d,T
 real(8):: kappa
 kappa = 1d0
end function kappa_r

! Get Planck mean opacity
function kappa_p(d,T) result(kappa)
 real(8),intent(in)::d,T
 real(8):: kappa
 kappa = 1d0
end function kappa_p


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE GET_RADFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Compute radiative fluxes at each interface

subroutine get_radflux(urad,d,T,cg,radK)

 use constants,only:c=>clight
 use grid,only: x1,xi1,x2,xi2,x3,xi3,dx1,dx2,dx3,g22,g33
 use miccg_mod,only:cg_set
 use utils,only:intpol

 real(8),allocatable,dimension(:,:,:),intent(in):: urad,d,T
 type(cg_set),intent(in):: cg
 real(8),allocatable,intent(inout):: radK(:,:,:,:)
 real(8):: gradE,R,rho,TT,xi
 integer:: i,j,k

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,gradE,rho,TT,xi,R) collapse(3)
 do k = cg%ks-1, cg%ke
  do j = cg%js-1, cg%je
   do i = cg%is-1, cg%ie
! Flux in x1 direction
    gradE = (urad(i+1,j,k) - urad(i,j,k))/dx1(i+1)
    rho = intpol(x1(i:i+1),d   (i:i+1,j,k),xi1(i))
    TT  = intpol(x1(i:i+1),T   (i:i+1,j,k),xi1(i))
    xi  = intpol(x1(i:i+1),urad(i:i+1,j,k),xi1(i))
    R = abs(gradE)/(kappa_r(rho,TT)*rho*xi)
    radK(1,i,j,k) = c*lambda(R)/(kappa_r(rho,TT)*rho)

! Flux in x2 direction
    gradE = (urad(i,j+1,k) - urad(i,j,k))/(dx2(j+1)*g22(i))
    rho = intpol(x2(j:j+1),d   (i,j:j+1,k),xi2(j))
    TT  = intpol(x2(j:j+1),T   (i,j:j+1,k),xi2(j))
    xi  = intpol(x2(j:j+1),urad(i,j:j+1,k),xi2(j))
    R = abs(gradE)/(kappa_r(rho,TT)*rho*xi)
    radK(2,i,j,k) = c*lambda(R)/(kappa_r(rho,TT)*rho)

! Flux in x3 direction
    gradE = (urad(i,j,k+1) - urad(i,j,k))/(dx3(k+1)*g33(i,j))
    rho = intpol(x3(k:k+1),d   (i,j,k:k+1),xi3(k))
    TT  = intpol(x3(k:k+1),T   (i,j,k:k+1),xi3(k))
    xi  = intpol(x3(k:k+1),urad(i,j,k:k+1),xi3(k))
    R = abs(gradE)/(kappa_r(rho,TT)*rho*xi)
    radK(3,i,j,k) = c*lambda(R)/(kappa_r(rho,TT)*rho)

   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_radflux

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE GET_RADA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Compute A matrix elements for radiation

subroutine get_radA(radK,cg)

 use constants,only:c=>clight,fac_egas,arad
 use settings,only:crdnt
 use grid,only:dt,dvol
 use physval,only:d,T,imu
 use miccg_mod,only:cg_set,ijk_from_l,get_preconditioner

 real(8),allocatable,intent(in):: radK(:,:,:,:)
 type(cg_set),intent(inout)::cg
 integer:: i,j,k,l,dim
 real(8):: kappap

!-----------------------------------------------------------------------------

 if(cg%in>1.and.cg%jn>1.and.cg%kn>1)then
  dim=3
 elseif(cg%in>1.and.(cg%jn>1.or.cg%kn>1))then
  dim=2
 elseif(cg%in>1.and.cg%jn==1.and.cg%kn==1.and.crdnt==2)then
  dim=1
 else
  print*,'Error in get_radA, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(1) ! 1D

  if(cg%jn==1.and.cg%kn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    cg%A(1,l) = dvol(i,j,k)/dt &
     + geo(1,i-1,j,k)*radK(1,i-1,j,k) + geo(1,i,j,k)*radK(1,i,j,k) &
     + dvol(i,j,k)*d(i,j,k)*fac_egas*imu(i,j,k)*kappap*c &
      /(fac_egas*imu(i,j,k)+4d0*kappap*c*arad*T(i,j,k)**3*dt)
    cg%A(2,l) = -geo(1,i,j,k)*radK(1,i,j,k)
    if(i==cg%ie)cg%A(2,l) = 0d0
   end do

  end if

 case(2) ! 2D

  if(cg%kn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    cg%A(1,l) = dvol(i,j,k)/dt &
     + geo(1,i-1,j,k)*radK(1,i-1,j,k) + geo(1,i,j,k)*radK(1,i,j,k) &
     + geo(2,i,j-1,k)*radK(2,i,j-1,k) + geo(2,i,j,k)*radK(2,i,j,k) &
     + dvol(i,j,k)*d(i,j,k)*fac_egas*imu(i,j,k)*kappap*c &
      /(fac_egas*imu(i,j,k)+4d0*kappap*c*arad*T(i,j,k)**3*dt)
    cg%A(2,l) = -geo(1,i,j,k)*radK(1,i,j,k)
    cg%A(3,l) = -geo(2,i,j,k)*radK(2,i,j,k)
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(j==cg%je)cg%A(3,l) = 0d0
   end do

  elseif(cg%jn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    cg%A(1,l) = dvol(i,j,k)/dt &
     + geo(1,i-1,j,k)*radK(1,i-1,j,k) + geo(1,i,j,k)*radK(1,i,j,k) &
     + geo(3,i,j,k-1)*radK(2,i,j-1,k) + geo(3,i,j,k)*radK(2,i,j,k) &
     + dvol(i,j,k)*d(i,j,k)*fac_egas*imu(i,j,k)*kappap*c &
      /(fac_egas*imu(i,j,k)+4d0*kappap*c*arad*T(i,j,k)**3*dt)
    cg%A(2,l) = -geo(1,i,j,k)*radK(1,i,j,k)
    cg%A(3,l) = -geo(3,i,j,k)*radK(3,i,j,k)
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(k==cg%ke)cg%A(3,l) = 0d0
   end do

  end if

 case(3) ! 3D

  if(crdnt==2)then
   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    cg%A(1,l) = dvol(i,j,k)/dt &
    + geo(1,i-1,j,k)*radK(1,i-1,j,k) + geo(1,i,j,k)*radK(1,i,j,k) &
    + geo(2,i,j-1,k)*radK(2,i,j-1,k) + geo(2,i,j,k)*radK(2,i,j,k) &
    + geo(3,i,j,k-1)*radK(3,i,j,k-1) + geo(3,i,j,k)*radK(3,i,j,k) &
    + dvol(i,j,k)*d(i,j,k)*fac_egas*imu(i,j,k)*kappap*c &
     /(fac_egas*imu(i,j,k)+4d0*kappap*c*arad*T(i,j,k)**3*dt)
    cg%A(2,l) = -geo(1,i,j,k)*radK(1,i,j,k)
    cg%A(3,l) = -geo(2,i,j,k)*radK(2,i,j,k)
    cg%A(4,l) = -geo(3,i,j,k)*radK(3,i,j,k)
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(j==cg%je)cg%A(3,l) = 0d0
    if(k==cg%ke)cg%A(4,l) = 0d0
    if(k/=cg%ks)cg%A(5,l) = 0d0
   end do
  end if

 end select

 call get_preconditioner(cg)

return
end subroutine get_radA

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SETUP_RADCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the necessary parameters for the CG method

subroutine setup_radcg(is,ie,js,je,ks,ke,cg)

 use settings,only:crdnt
 use miccg_mod,only:cg_set

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
 elseif(ie>is.and.je==js.and.ke==ks.and.crdnt==2)then
  dim=1
 else
  print*,'Error in setup_radcg, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(1) ! 1D
! 1D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==1)then
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

  end if

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
!                       SUBROUTINE GET_SOURCE_TERM
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute the source term for radiation

subroutine get_source_term(d,T,urad,dt,cg,rsrc)

 use grid,only:dim,crdnt
 use miccg_mod,only:cg_set,ijk_from_l

 real(8),allocatable,dimension(:,:,:),intent(in):: d,T,urad
 real(8),intent(in):: dt
 type(cg_set),intent(inout):: cg
 real(8),allocatable,intent(inout):: rsrc(:,:,:)
 integer:: i,j,k,l

!-----------------------------------------------------------------------------

 select case(dim)
 case(1) ! 1D

  if(cg%jn==1.and.cg%kn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    rsrc(i,j,k) = 0d0
    if(i==cg%ie)cg%A(2,l) = 0d0
   end do

  end if

 case(2) ! 2D

  if(cg%kn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    rsrc(i,j,k) = 0d0
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(j==cg%je)cg%A(3,l) = 0d0
   end do

  elseif(cg%jn==1)then

   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    rsrc(i,j,k) = 0d0
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(k==cg%ke)cg%A(3,l) = 0d0
   end do

  end if

 case(3) ! 3D

  if(crdnt==2)then
   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,i,j,k)
    rsrc(i,j,k) = 0d0
    if(i==cg%ie)cg%A(2,l) = 0d0
    if(j==cg%je)cg%A(3,l) = 0d0
    if(k==cg%ke)cg%A(4,l) = 0d0
    if(k/=cg%ks)cg%A(5,l) = 0d0
   end do
  end if

 end select


return
end subroutine get_source_term

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GET_GEO
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute geometrical coefficients for diffusion equation

subroutine get_geo

 use grid,only:is,ie,js,je,ks,ke,sa1,sa2,sa3,dx1,dx2,dx3,g22,g33

 integer:: i,j,k

!-----------------------------------------------------------------------------

 allocate(geo(1:3,is:ie,js:je,ks:ke))

!$omp parallel do private(i,j,k) collapse(3)
 do k = ks-1, ke
  do j = js-1, je
   do i = is-1, ie
    geo(1,i,j,k) = sa1(i  ,j,k) / dx1(i+1)
    geo(2,i,j,k) = sa2(i,j  ,k) / (dx2(j+1)*g22(i))
    geo(3,i,j,k) = sa3(i,j,k  ) / (dx3(k+1)*g33(i,j))
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_geo

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE RADIATION_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up radiation parameters

subroutine radiation_setup

 use settings,only:radswitch
 use grid,only:is,ie,js,je,ks,ke
 use miccg_mod,only:cg_rad

!-----------------------------------------------------------------------------

 if(radswitch==1)then
  call setup_radcg(is,ie,js,je,ks,ke,cg_rad)
  call get_geo
 end if

return
end subroutine radiation_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE RADIATIVE_FORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute radiative acceleration.

subroutine radiative_force

 use grid

!-----------------------------------------------------------------------------



return
end subroutine radiative_force

end module radiation_mod
