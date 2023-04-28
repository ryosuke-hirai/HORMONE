module boundary_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE BOUNDARYCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary condition

subroutine boundarycondition

 use settings
 use grid
 use physval
 use dirichlet_mod
 use pressure_mod
 use composition_mod

 real(8),allocatable:: plug(:,:,:)

!-------------------------------------------------------------------------

! Note:
!  0: Periodic b.c.
!  1: Reflective b.c.
!  2: Outgoing b.c.
!  3: Free b.c.
!  4: Linear Extrapolation
!  5: Linear Extrapolation but Outgoing (only for vectors)
!  9: Dirichlet b.c. (boundary values should be given elsewhere!)

 call pressure
 if(compswitch==2)call spcboundary
 call meanmolweight

! x1-direction ***********************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! scalar values
 x1_inner_scalar: select case (bc1is)
 case(0) x1_inner_scalar ! periodic
  d (is-2:is-1,js:je,ks:ke) = d (ie-1:ie,js:je,ks:ke)
  p (is-2:is-1,js:je,ks:ke) = p (ie-1:ie,js:je,ks:ke)
  b1(is-2:is-1,js:je,ks:ke) = b1(ie-1:ie,js:je,ks:ke)
  b2(is-2:is-1,js:je,ks:ke) = b2(ie-1:ie,js:je,ks:ke)
  b3(is-2:is-1,js:je,ks:ke) = b3(ie-1:ie,js:je,ks:ke)
  phi(is-2:is-1,js:je,ks:ke) = phi(ie-1:ie,js:je,ks:ke)
 case(1) x1_inner_scalar ! reflective
  d (is-2:is-1,js:je,ks:ke) = d (is+1:is:-1,js:je,ks:ke)
  p (is-2:is-1,js:je,ks:ke) = p (is+1:is:-1,js:je,ks:ke)
  b1(is-2:is-1,js:je,ks:ke) = b1(is+1:is:-1,js:je,ks:ke)
  b2(is-2:is-1,js:je,ks:ke) = b2(is+1:is:-1,js:je,ks:ke)
  b3(is-2:is-1,js:je,ks:ke) = b3(is+1:is:-1,js:je,ks:ke)
 case(2:3) x1_inner_scalar ! outgoing/free
  d (is-2:is-1,js:je,ks:ke) = spread(d (is,js:je,ks:ke),1,2)
  p (is-2:is-1,js:je,ks:ke) = spread(p (is,js:je,ks:ke),1,2)
  b1(is-2:is-1,js:je,ks:ke) = spread(b1(is,js:je,ks:ke),1,2)
  b2(is-2:is-1,js:je,ks:ke) = spread(b2(is,js:je,ks:ke),1,2)
  b3(is-2:is-1,js:je,ks:ke) = spread(b3(is,js:je,ks:ke),1,2)
 case(9) x1_inner_scalar ! Dirichlet
  d (is-2:is-1,js:je,ks:ke) = d0 (is-2:is-1,js:je,ks:ke)
  p (is-2:is-1,js:je,ks:ke) = p0 (is-2:is-1,js:je,ks:ke)
  b1(is-2:is-1,js:je,ks:ke) = b10(is-2:is-1,js:je,ks:ke)
  b2(is-2:is-1,js:je,ks:ke) = b20(is-2:is-1,js:je,ks:ke)
  b3(is-2:is-1,js:je,ks:ke) = b30(is-2:is-1,js:je,ks:ke)
 case default x1_inner_scalar
  print *, "Error from x1 scalar inner boundary condition" ; stop
 end select x1_inner_scalar

! vector values
 x1_inner_vector: select case (bc1iv)
 case (0) x1_inner_vector ! periodic
  v1(is-2:is-1,js:je,ks:ke) = v1(ie-1:ie,js:je,ks:ke)
  v2(is-2:is-1,js:je,ks:ke) = v2(ie-1:ie,js:je,ks:ke)
  v3(is-2:is-1,js:je,ks:ke) = v3(ie-1:ie,js:je,ks:ke)
 case(1) x1_inner_vector ! reflective
  v1(is-2:is-1,js:je,ks:ke) = -v1(is+1:is:-1,js:je,ks:ke)
  v2(is-2:is-1,js:je,ks:ke) =  v2(is+1:is:-1,js:je,ks:ke)
  v3(is-2:is-1,js:je,ks:ke) =  v3(is+1:is:-1,js:je,ks:ke)
 case(2) x1_inner_vector ! outgoing
  allocate( plug(1:1,js:je,ks:ke) )
  plug(1,:,:) = 0.5d0 - sign(0.5d0,v1(is,js:je,ks:ke))
  v1(is-2:is-1,js:je,ks:ke) = spread(v1(is,js:je,ks:ke)*plug(1,:,:),1,2) 
  v2(is-2:is-1,js:je,ks:ke) = spread(v2(is,js:je,ks:ke)*plug(1,:,:),1,2) 
  v3(is-2:is-1,js:je,ks:ke) = spread(v3(is,js:je,ks:ke)*plug(1,:,:),1,2) 
  deallocate( plug )
 case(3) x1_inner_vector ! free
  v1(is-2:is-1,js:je,ks:ke) = spread(v1(is,js:je,ks:ke),1,2) 
  v2(is-2:is-1,js:je,ks:ke) = spread(v2(is,js:je,ks:ke),1,2) 
  v3(is-2:is-1,js:je,ks:ke) = spread(v3(is,js:je,ks:ke),1,2) 
 case(9) x1_inner_vector ! Dirichlet
  v1(is-2:is-1,js:je,ks:ke) = v10(is-2:is-1,js:je,ks:ke)
  v2(is-2:is-1,js:je,ks:ke) = v20(is-2:is-1,js:je,ks:ke)
  v3(is-2:is-1,js:je,ks:ke) = v30(is-2:is-1,js:je,ks:ke)
 case default x1_inner_vector
  print *, "Error from x1 velocity inner boundary condition" ; stop
 end select x1_inner_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ks, ke
  do j = js, je
   do i = is-2, is-1
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(is,j,k)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! scalar values
 x1_outer_scalar: select case (bc1os)
 case(0) x1_outer_scalar ! periodic
  d (ie+1:ie+2,js:je,ks:ke) = d (is:is+1,js:je,ks:ke)
  p (ie+1:ie+2,js:je,ks:ke) = p (is:is+1,js:je,ks:ke)
  b1(ie+1:ie+2,js:je,ks:ke) = b1(is:is+1,js:je,ks:ke)
  b2(ie+1:ie+2,js:je,ks:ke) = b2(is:is+1,js:je,ks:ke)
  b3(ie+1:ie+2,js:je,ks:ke) = b3(is:is+1,js:je,ks:ke)
  phi(ie+1:ie+2,js:je,ks:ke) = phi(is:is+1,js:je,ks:ke)
 case(1) x1_outer_scalar ! reflective
  d (ie+1:ie+2,js:je,ks:ke) = d (ie:ie-1:-1,js:je,ks:ke)
  p (ie+1:ie+2,js:je,ks:ke) = p (ie:ie-1:-1,js:je,ks:ke)
  b1(ie+1:ie+2,js:je,ks:ke) = b1(ie:ie-1:-1,js:je,ks:ke)
  b2(ie+1:ie+2,js:je,ks:ke) = b2(ie:ie-1:-1,js:je,ks:ke)
  b3(ie+1:ie+2,js:je,ks:ke) = b3(ie:ie-1:-1,js:je,ks:ke)
 case(2:3) x1_outer_scalar ! outgoing/free
  d (ie+1:ie+2,js:je,ks:ke) = spread(d (ie,js:je,ks:ke),1,2)
  p (ie+1:ie+2,js:je,ks:ke) = spread(p (ie,js:je,ks:ke),1,2)
  b1(ie+1:ie+2,js:je,ks:ke) = spread(b1(ie,js:je,ks:ke),1,2) 
  b2(ie+1:ie+2,js:je,ks:ke) = spread(b2(ie,js:je,ks:ke),1,2) 
  b3(ie+1:ie+2,js:je,ks:ke) = spread(b3(ie,js:je,ks:ke),1,2)
 case(4:5) x1_outer_scalar ! linear
  d (ie+1,js:je,ks:ke) = d (ie  ,js:je,ks:ke) &
                       + (d (ie,js:je,ks:ke)-d (ie-1,js:je,ks:ke)) &
                       * dx1(ie)/dx1(ie-1)
  d (ie+2,js:je,ks:ke) = d (ie+1,js:je,ks:ke) &
                       + (d (ie+1,js:je,ks:ke)-d (ie,js:je,ks:ke)) &
                       * dx1(ie+1)/dx1(ie)
  p (ie+1,:,:) = p (ie  ,:,:) + (p (ie,:,:)-p (ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  p (ie+2,:,:) = p (ie+1,:,:) + (p (ie+1,:,:)-p (ie,:,:)) * dx1(ie+1)/dx1(ie)
!$omp parallel do private(i,j,k)
  do k = ks, ke
   do j = js, je
    do i = ie+1, ie+2
     if(d(i,j,k)<0d0)then
      d(i,j,k) = d(ie,j,k)
     end if
     if(p(i,j,k)<0d0)then
      p(i,j,k) = p(ie,j,k)
     end if
    end do
   end do
  end do
!$omp end parallel do
  b1(ie+1,:,:) = b1(ie  ,:,:) + (b1(ie,:,:)-b1(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  b1(ie+2,:,:) = b1(ie+1,:,:) + (b1(ie+1,:,:)-b1(ie,:,:)) * dx1(ie+1)/dx1(ie)
  b2(ie+1,:,:) = b2(ie  ,:,:) + (b2(ie,:,:)-b2(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  b2(ie+2,:,:) = b2(ie+1,:,:) + (b2(ie+1,:,:)-b2(ie,:,:)) * dx1(ie+1)/dx1(ie)
  b3(ie+1,:,:) = b3(ie  ,:,:) + (b3(ie,:,:)-b3(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  b3(ie+2,:,:) = b3(ie+1,:,:) + (b3(ie+1,:,:)-b3(ie,:,:)) * dx1(ie+1)/dx1(ie)
 case(9) x1_outer_scalar ! Dirichlet
  d (ie+1:ie+2,js:je,ks:ke) = d0 (ie+1:ie+2,js:je,ks:ke)
  p (ie+1:ie+2,js:je,ks:ke) = p0 (ie+1:ie+2,js:je,ks:ke)
  b1(ie+1:ie+2,js:je,ks:ke) = b10(ie+1:ie+2,js:je,ks:ke)
  b2(ie+1:ie+2,js:je,ks:ke) = b20(ie+1:ie+2,js:je,ks:ke)
  b3(ie+1:ie+2,js:je,ks:ke) = b30(ie+1:ie+2,js:je,ks:ke)
 case default x1_outer_scalar
  print *, "Error from x1 scalar outer boundary condition" ; stop
 end select x1_outer_scalar

! vector values
 x1_outer_vector: select case (bc1ov)
 case(0) x1_outer_vector ! periodic
  v1(ie+1:ie+2,js:je,ks:ke) = v1(is:is+1,js:je,ks:ke)
  v2(ie+1:ie+2,js:je,ks:ke) = v2(is:is+1,js:je,ks:ke)
  v3(ie+1:ie+2,js:je,ks:ke) = v3(is:is+1,js:je,ks:ke)
 case(1) x1_outer_vector ! reflective
  v1(ie+1:ie+2,js:je,ks:ke) = -v1(ie:ie-1:-1,js:je,ks:ke)
  v2(ie+1:ie+2,js:je,ks:ke) =  v2(ie:ie-1:-1,js:je,ks:ke)
  v3(ie+1:ie+2,js:je,ks:ke) =  v3(ie:ie-1:-1,js:je,ks:ke)
 case(2) x1_outer_vector ! outgoing
  allocate( plug(1:1,js:je,ks:ke) )
  plug(1,:,:) = 0.5d0 + sign(0.5d0,v1(ie,js:je,ks:ke))
  v1(ie+1:ie+2,js:je,ks:ke) = spread(v1(ie,js:je,ks:ke)*plug(1,:,:),1,2) 
  v2(ie+1:ie+2,js:je,ks:ke) = spread(v2(ie,js:je,ks:ke)*plug(1,:,:),1,2) 
  v3(ie+1:ie+2,js:je,ks:ke) = spread(v3(ie,js:je,ks:ke)*plug(1,:,:),1,2) 
  deallocate( plug )
 case(3) x1_outer_vector ! free
  v1(ie+1:ie+2,js:je,ks:ke) = spread(v1(ie,js:je,ks:ke),1,2)
  v2(ie+1:ie+2,js:je,ks:ke) = spread(v2(ie,js:je,ks:ke),1,2)
  v3(ie+1:ie+2,js:je,ks:ke) = spread(v3(ie,js:je,ks:ke),1,2)
 case(4) x1_outer_vector ! linear
  v1(ie+1,:,:) = v1(ie  ,:,:) + (v1(ie,:,:)-v1(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v1(ie+2,:,:) = v1(ie+1,:,:) + (v1(ie+1,:,:)-v1(ie,:,:)) * dx1(ie+1)/dx1(ie)
  v2(ie+1,:,:) = v2(ie  ,:,:) + (v2(ie,:,:)-v2(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v2(ie+2,:,:) = v2(ie+1,:,:) + (v2(ie+1,:,:)-v2(ie,:,:)) * dx1(ie+1)/dx1(ie)
  v3(ie+1,:,:) = v3(ie  ,:,:) + (v3(ie,:,:)-v3(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v3(ie+2,:,:) = v3(ie+1,:,:) + (v3(ie+1,:,:)-v3(ie,:,:)) * dx1(ie+1)/dx1(ie)
 case(5) x1_outer_vector ! linear + outgoing
  v1(ie+1,:,:) = v1(ie  ,:,:) + (v1(ie,:,:)-v1(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v1(ie+2,:,:) = v1(ie+1,:,:) + (v1(ie+1,:,:)-v1(ie,:,:)) * dx1(ie+1)/dx1(ie)
  v1(ie+1,:,:) = v1(ie+1,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
  v1(ie+2,:,:) = v1(ie+2,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
  v2(ie+1,:,:) = v2(ie  ,:,:) + (v2(ie,:,:)-v2(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v2(ie+2,:,:) = v2(ie+1,:,:) + (v2(ie+1,:,:)-v2(ie,:,:)) * dx1(ie+1)/dx1(ie)
  v2(ie+1,:,:) = v2(ie+1,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
  v2(ie+2,:,:) = v2(ie+2,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
  v3(ie+1,:,:) = v3(ie  ,:,:) + (v3(ie,:,:)-v3(ie-1,:,:)) * dx1(ie)/dx1(ie-1)
  v3(ie+2,:,:) = v3(ie+1,:,:) + (v3(ie+1,:,:)-v3(ie,:,:)) * dx1(ie+1)/dx1(ie)
  v3(ie+1,:,:) = v3(ie+1,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
  v3(ie+2,:,:) = v3(ie+2,:,:) * (0.5d0+sign(0.5d0,v1(ie,:,:)))
 case(9) x1_outer_vector ! Dirichlet
  v1(ie+1:ie+2,js:je,ks:ke) = v10(ie+1:ie+2,js:je,ks:ke)
  v2(ie+1:ie+2,js:je,ks:ke) = v20(ie+1:ie+2,js:je,ks:ke)
  v3(ie+1:ie+2,js:je,ks:ke) = v30(ie+1:ie+2,js:je,ks:ke)
 case default x1_outer_vector
  print *, "Error from x1 velocity outer boundary condition" ; stop
 end select x1_outer_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ks, ke
  do j = js, je
   do i = ie+1, ie+2
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(ie,j,k)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do
! x2-direction ***********************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(je>1)then
! scalar values
 x2_inner_scalar: select case (bc2is)
 case(0) x2_inner_scalar ! periodic
  d (is:ie,js-2:js-1,ks:ke) = d (is:ie,je-1:je,ks:ke)
  p (is:ie,js-2:js-1,ks:ke) = p (is:ie,je-1:je,ks:ke)
  b1(is:ie,js-2:js-1,ks:ke) = b1(is:ie,je-1:je,ks:ke)
  b2(is:ie,js-2:js-1,ks:ke) = b2(is:ie,je-1:je,ks:ke)
  b3(is:ie,js-2:js-1,ks:ke) = b3(is:ie,je-1:je,ks:ke)
  phi(is:ie,js-2:js-1,ks:ke) = phi(is:ie,je-1:je,ks:ke)
 case(1) x2_inner_scalar ! reflective
  d (is:ie,js-2:js-1,ks:ke) = d (is:ie,js+1:js:-1,ks:ke)
  p (is:ie,js-2:js-1,ks:ke) = p (is:ie,js+1:js:-1,ks:ke)
  b1(is:ie,js-2:js-1,ks:ke) = b1(is:ie,js+1:js:-1,ks:ke)
  b2(is:ie,js-2:js-1,ks:ke) = b2(is:ie,js+1:js:-1,ks:ke)
  b3(is:ie,js-2:js-1,ks:ke) = b3(is:ie,js+1:js:-1,ks:ke)
 case(2:3) x2_inner_scalar ! outgoing/free
!$omp parallel do private (i,j,k)
  do k = ks, ke ; do j = js-2, js-1 ; do i = is, ie
   d (i,j,k) = d (i,js,k) ; p (i,j,k) = p (i,js,k)
   b1(i,j,k) = b1(i,js,k) ; b2(i,j,k) = b2(i,js,k)
   b3(i,j,k) = b3(i,js,k)
  end do; end do ; end do
!$omp end parallel do
!!$   d (is:ie,js-2:js-1,ks:ke) = d (is:ie,js,ks:ke)
!!$   p (is:ie,js-2:js-1,ks:ke) = spread(p (is:ie,js,ks:ke),2,2)
!!$   b1(is:ie,js-2:js-1,ks:ke) = spread(b1(is:ie,js,ks:ke),2,2) 
!!$   b2(is:ie,js-2:js-1,ks:ke) = spread(b2(is:ie,js,ks:ke),2,2) 
!!$   b3(is:ie,js-2:js-1,ks:ke) = spread(b3(is:ie,js,ks:ke),2,2) 
 case(9) x2_inner_scalar ! Dirichlet
  d (is:ie,js-2:js-1,ks:ke) = d0 (is:ie,js-2:js-1,ks:ke)
  p (is:ie,js-2:js-1,ks:ke) = p0 (is:ie,js-2:js-1,ks:ke)
  b1(is:ie,js-2:js-1,ks:ke) = b10(is:ie,js-2:js-1,ks:ke)
  b2(is:ie,js-2:js-1,ks:ke) = b20(is:ie,js-2:js-1,ks:ke)
  b3(is:ie,js-2:js-1,ks:ke) = b30(is:ie,js-2:js-1,ks:ke)
 case default x2_inner_scalar
  print *, "Error from x2 scalar inner boundary condition" ; stop
 end select x2_inner_scalar

! vector values
 x2_inner_vector: select case (bc2iv)
 case(0) x2_inner_vector ! periodic
  v1(is:ie,js-2:js-1,ks:ke) = v1(is:ie,je-1:je,ks:ke)
  v2(is:ie,js-2:js-1,ks:ke) = v2(is:ie,je-1:je,ks:ke)
  v3(is:ie,js-2:js-1,ks:ke) = v3(is:ie,je-1:je,ks:ke)
 case(1) x2_inner_vector ! reflective
  v1(is:ie,js-2:js-1,ks:ke) =  v1(is:ie,js+1:js:-1,ks:ke)
  v2(is:ie,js-2:js-1,ks:ke) = -v2(is:ie,js+1:js:-1,ks:ke)
  v3(is:ie,js-2:js-1,ks:ke) =  v3(is:ie,js+1:js:-1,ks:ke)
 case(2) x2_inner_vector ! outgoing
  allocate( plug(is:ie,1:1,ks:ke) )
  plug(:,1,:) = 0.5d0 - sign(0.5d0,v2(is:ie,js,ks:ke))
  v1(is:ie,js-2:js-1,ks:ke) = spread(v1(is:ie,js,ks:ke)*plug(:,1,:),2,2) 
  v2(is:ie,js-2:js-1,ks:ke) = spread(v2(is:ie,js,ks:ke)*plug(:,1,:),2,2) 
  v3(is:ie,js-2:js-1,ks:ke) = spread(v3(is:ie,js,ks:ke)*plug(:,1,:),2,2) 
  deallocate( plug )
 case(3) x2_inner_vector ! free
  do k = ks, ke
   do j = js-2, js-1
    do i = is, ie
     v1(i,j,k) = v1(i,js,k) ; v2(i,j,k) = v2(i,js,k) ; v3(i,j,k) = v3(i,js,k)
    end do
   end do
  end do
!!$   v1(is:ie,js-2:js-1,ks:ke) = spread(v1(is:ie,js,ks:ke),2,2) 
!!$   v2(is:ie,js-2:js-1,ks:ke) = spread(v2(is:ie,js,ks:ke),2,2) 
!!$   v3(is:ie,js-2:js-1,ks:ke) = spread(v3(is:ie,js,ks:ke),2,2) 
 case(9) x2_inner_vector ! Dirichlet
  v1(is:ie,js-2:js-1,ks:ke) = v10(is:ie,js-2:js-1,ks:ke)
  v2(is:ie,js-2:js-1,ks:ke) = v20(is:ie,js-2:js-1,ks:ke)
  v3(is:ie,js-2:js-1,ks:ke) = v30(is:ie,js-2:js-1,ks:ke)
 case default x2_inner_vector
  print *, "Error from x2 velocity inner boundary condition" ; stop
 end select x2_inner_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ks, ke
  do j = js-2, js-1
   do i = is, ie
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(i,js,k)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! scalar values
 x2_outer_scalar: select case (bc2os)
 case(0) x2_outer_scalar ! periodic
  d (is:ie,je+1:je+2,ks:ke) = d (is:ie,js:js+1,ks:ke)
  p (is:ie,je+1:je+2,ks:ke) = p (is:ie,js:js+1,ks:ke)
  b1(is:ie,je+1:je+2,ks:ke) = b1(is:ie,js:js+1,ks:ke)
  b2(is:ie,je+1:je+2,ks:ke) = b2(is:ie,js:js+1,ks:ke)
  b3(is:ie,je+1:je+2,ks:ke) = b3(is:ie,js:js+1,ks:ke)
  phi(is:ie,je+1:je+2,ks:ke) = phi(is:ie,js:js+1,ks:ke)
 case(1) x2_outer_scalar ! reflective
  d (is:ie,je+1:je+2,ks:ke) = d (is:ie,je:je-1:-1,ks:ke)
  p (is:ie,je+1:je+2,ks:ke) = p (is:ie,je:je-1:-1,ks:ke)
  b1(is:ie,je+1:je+2,ks:ke) = b1(is:ie,je:je-1:-1,ks:ke)
  b2(is:ie,je+1:je+2,ks:ke) = b2(is:ie,je:je-1:-1,ks:ke)
  b3(is:ie,je+1:je+2,ks:ke) = b3(is:ie,je:je-1:-1,ks:ke)
 case(2:3) x2_outer_scalar ! outgoing/free
  do i = is, ie ; do j = je+1, je+2 ; do k = ks, ke
   d (i,j,k) = d (i,je,k)
   p (i,j,k) = p (i,je,k)
   b1(i,j,k) = b1(i,je,k)
   b2(i,j,k) = b2(i,je,k)
   b3(i,j,k) = b3(i,je,k)
  end do; end do ; end do
!!$   d (is:ie,je+1:je+2,ks:ke) = spread(d (is:ie,je,ks:ke),2,2)
!!$   p (is:ie,je+1:je+2,ks:ke) = spread(p (is:ie,je,ks:ke),2,2)
!!$   b1(is:ie,je+1:je+2,ks:ke) = spread(b1(is:ie,je,ks:ke),2,2) 
!!$   b2(is:ie,je+1:je+2,ks:ke) = spread(b2(is:ie,je,ks:ke),2,2) 
!!$   b3(is:ie,je+1:je+2,ks:ke) = spread(b3(is:ie,je,ks:ke),2,2) 
 case(9) x2_outer_scalar ! Dirichlet
  d (is:ie,je+1:je+2,ks:ke) = d0 (is:ie,je+1:je+2,ks:ke)
  p (is:ie,je+1:je+2,ks:ke) = p0 (is:ie,je+1:je+2,ks:ke)
  b1(is:ie,je+1:je+2,ks:ke) = b10(is:ie,je+1:je+2,ks:ke)
  b2(is:ie,je+1:je+2,ks:ke) = b20(is:ie,je+1:je+2,ks:ke)
  b3(is:ie,je+1:je+2,ks:ke) = b30(is:ie,je+1:je+2,ks:ke)
 case default x2_outer_scalar
  print *, "Error from x2 scalar outer boundary condition" ; stop
 end select x2_outer_scalar

! vector values
 x2_outer_vector: select case (bc2ov)
 case(0) x2_outer_vector ! periodic
  v1(is:ie,je+1:je+2,ks:ke) = v1(is:ie,js:js+1,ks:ke)
  v2(is:ie,je+1:je+2,ks:ke) = v2(is:ie,js:js+1,ks:ke)
  v3(is:ie,je+1:je+2,ks:ke) = v3(is:ie,js:js+1,ks:ke)
 case(1) x2_outer_vector ! reflective
  v1(is:ie,je+1:je+2,ks:ke) =  v1(is:ie,je:je-1:-1,ks:ke)
  v2(is:ie,je+1:je+2,ks:ke) = -v2(is:ie,je:je-1:-1,ks:ke)
  v3(is:ie,je+1:je+2,ks:ke) =  v3(is:ie,je:je-1:-1,ks:ke)
 case(2) x2_outer_vector ! outgoing
  allocate( plug(is:ie,1:1,ks:ke) )
  plug(:,1,:) = 0.5d0 + sign(0.5d0,v2(is:ie,je,ks:ke))
  v1(is:ie,je+1:je+2,ks:ke) = spread(v1(is:ie,je,ks:ke)*plug(:,1,:),2,2) 
  v2(is:ie,je+1:je+2,ks:ke) = spread(v2(is:ie,je,ks:ke)*plug(:,1,:),2,2) 
  v3(is:ie,je+1:je+2,ks:ke) = spread(v3(is:ie,je,ks:ke)*plug(:,1,:),2,2) 
  deallocate( plug )
 case(3) x2_outer_vector ! free
  do i = is, ie
   do j = je+1, je+2
    do k = ks, ke
     v1(i,j,k) = v1(i,je,k) ; v2(i,j,k) = v2(i,je,k) ; v3(i,j,k) = v3(i,je,k)
    end do
   end do
  end do
!!$   v1(is:ie,je+1:je+2,ks:ke) = spread(v1(is:ie,je,ks:ke),2,2) 
!!$   v2(is:ie,je+1:je+2,ks:ke) = spread(v2(is:ie,je,ks:ke),2,2) 
!!$   v3(is:ie,je+1:je+2,ks:ke) = spread(v3(is:ie,je,ks:ke),2,2) 
 case(9) x2_outer_vector ! Dirichlet
  v1(is:ie,je+1:je+2,ks:ke) = v10(is:ie,je+1:je+2,ks:ke)
  v2(is:ie,je+1:je+2,ks:ke) = v20(is:ie,je+1:je+2,ks:ke)
  v3(is:ie,je+1:je+2,ks:ke) = v30(is:ie,je+1:je+2,ks:ke)
 case default x2_outer_vector
  print *, "Error from x2 velocity outer boundary condition" ; stop
 end select x2_outer_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ks, ke
  do j = je+1, je+2
   do i = is, ie
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(i,je,k)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do
end if
! x3-direction ***********************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(ke>1)then
! scalar values
 x3_inner_scalar: select case (bc3is)
 case(0) x3_inner_scalar ! periodic
  d (is:ie,js:je,ks-2:ks-1) = d (is:ie,js:je,ke-1:ke)
  p (is:ie,js:je,ks-2:ks-1) = p (is:ie,js:je,ke-1:ke)
  b1(is:ie,js:je,ks-2:ks-1) = b1(is:ie,js:je,ke-1:ke)
  b2(is:ie,js:je,ks-2:ks-1) = b2(is:ie,js:je,ke-1:ke)
  b3(is:ie,js:je,ks-2:ks-1) = b3(is:ie,js:je,ke-1:ke)
  phi(is:ie,js:je,ks-2:ks-1) = phi(is:ie,js:je,ke-1:ke)
 case(1) x3_inner_scalar ! reflective
  d (is:ie,js:je,ks-2:ks-1) = d (is:ie,js:je,ks+1:ks:-1)
  p (is:ie,js:je,ks-2:ks-1) = p (is:ie,js:je,ks+1:ks:-1)
  b1(is:ie,js:je,ks-2:ks-1) = b1(is:ie,js:je,ks+1:ks:-1)
  b2(is:ie,js:je,ks-2:ks-1) = b2(is:ie,js:je,ks+1:ks:-1)
  b3(is:ie,js:je,ks-2:ks-1) = b3(is:ie,js:je,ks+1:ks:-1)
 case(2:3) x3_inner_scalar ! outgoing/free
  d (is:ie,js:je,ks-2:ks-1) = spread(d (is:ie,js:je,ks),3,2)
  p (is:ie,js:je,ks-2:ks-1) = spread(p (is:ie,js:je,ks),3,2)
  b1(is:ie,js:je,ks-2:ks-1) = spread(b1(is:ie,js:je,ks),3,2) 
  b2(is:ie,js:je,ks-2:ks-1) = spread(b2(is:ie,js:je,ks),3,2) 
  b3(is:ie,js:je,ks-2:ks-1) = spread(b3(is:ie,js:je,ks),3,2)
 case(4:5) x3_inner_scalar ! linear
  d (:,:,ks-1) = d (:,:,ks  ) - (d (:,:,ks+1)-d (:,:,ks)) * dx3(ks)/dx3(ks+1)
  d (:,:,ks-2) = d (:,:,ks-1) - (d (:,:,ks)-d (:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  p (:,:,ks-1) = p (:,:,ks  ) - (p (:,:,ks+1)-p (:,:,ks)) * dx3(ks)/dx3(ks+1)
  p (:,:,ks-2) = p (:,:,ks-1) - (p (:,:,ks)-p (:,:,ks-1)) * dx3(ks-1)/dx3(ks)
!$omp parallel do private (i,j,k)
  do k = ks-2, ks-1
   do j = js, je
    do i = is, ie
     if(d(i,j,k)<0d0)then
      d(i,j,k) = d(i,j,ks)
     end if
     if(p(i,j,k)<0d0)then
      p(i,j,k) = p(i,j,ks)
     end if
    end do
   end do
  end do
!$omp end parallel do
  b1(:,:,ks-1) = b1(:,:,ks  ) - (b1(:,:,ks+1)-b1(:,:,ks)) * dx3(ks)/dx3(ks+1)
  b1(:,:,ks-2) = b1(:,:,ks-1) - (b1(:,:,ks)-b1(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  b2(:,:,ks-1) = b2(:,:,ks  ) - (b2(:,:,ks+1)-b2(:,:,ks)) * dx3(ks)/dx3(ks+1)
  b2(:,:,ks-2) = b2(:,:,ks-1) - (b2(:,:,ks)-b2(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  b3(:,:,ks-1) = b3(:,:,ks  ) - (b3(:,:,ks+1)-b3(:,:,ks)) * dx3(ks)/dx3(ks+1)
  b3(:,:,ks-2) = b3(:,:,ks-1) - (b3(:,:,ks)-b3(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
 case(9) x3_inner_scalar ! Dirichlet
  d (is:ie,js:je,ks-2:ks-1) = d0 (is:ie,js:je,ks-2:ks-1)
  p (is:ie,js:je,ks-2:ks-1) = p0 (is:ie,js:je,ks-2:ks-1)
  b1(is:ie,js:je,ks-2:ks-1) = b10(is:ie,js:je,ks-2:ks-1)
  b2(is:ie,js:je,ks-2:ks-1) = b20(is:ie,js:je,ks-2:ks-1)
  b3(is:ie,js:je,ks-2:ks-1) = b30(is:ie,js:je,ks-2:ks-1)
 case default x3_inner_scalar
  print *, "Error from x3 scalar inner boundary condition" ; stop
 end select x3_inner_scalar

! vector values
 x3_inner_vector: select case (bc3iv)
 case(0) x3_inner_vector ! periodic
  v1(is:ie,js:je,ks-2:ks-1) = v1(is:ie,js:je,ke-1:ke)
  v2(is:ie,js:je,ks-2:ks-1) = v2(is:ie,js:je,ke-1:ke)
  v3(is:ie,js:je,ks-2:ks-1) = v3(is:ie,js:je,ke-1:ke)
 case(1) x3_inner_vector ! reflective
  v1(is:ie,js:je,ks-2:ks-1) =  v1(is:ie,js:je,ks+1:ks:-1)
  v2(is:ie,js:je,ks-2:ks-1) =  v2(is:ie,js:je,ks+1:ks:-1)
  v3(is:ie,js:je,ks-2:ks-1) = -v3(is:ie,js:je,ks+1:ks:-1)
 case(2) x3_inner_vector ! outgoing
  allocate( plug(is:ie,js:je,1:1) )
  plug(:,:,1) = 0.5d0 - sign(0.5d0,v3(is:ie,js:je,ks))
  v1(is:ie,js:je,ks-2:ks-1) = spread(v1(is:ie,js:je,ks)*plug(:,:,1),3,2) 
  v2(is:ie,js:je,ks-2:ks-1) = spread(v2(is:ie,js:je,ks)*plug(:,:,1),3,2) 
  v3(is:ie,js:je,ks-2:ks-1) = spread(v3(is:ie,js:je,ks)*plug(:,:,1),3,2) 
  deallocate( plug )
 case(3) x3_inner_vector ! free
  v1(is:ie,js:je,ks-2:ks-1) = spread(v1(is:ie,js:je,ks),3,2)
  v2(is:ie,js:je,ks-2:ks-1) = spread(v2(is:ie,js:je,ks),3,2)
  v3(is:ie,js:je,ks-2:ks-1) = spread(v3(is:ie,js:je,ks),3,2)
 case(4) x3_inner_vector ! linear
  v1(:,:,ks-1) = v1(:,:,ks  ) - (v1(:,:,ks+1)-v1(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v1(:,:,ks-2) = v1(:,:,ks-1) - (v1(:,:,ks)-v1(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  v2(:,:,ks-1) = v2(:,:,ks  ) - (v2(:,:,ks+1)-v2(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v2(:,:,ks-2) = v2(:,:,ks-1) - (v2(:,:,ks)-v2(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  v3(:,:,ks-1) = v3(:,:,ks  ) - (v3(:,:,ks+1)-v3(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v3(:,:,ks-2) = v3(:,:,ks-1) - (v3(:,:,ks)-v3(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
 case(5) x3_inner_vector ! linear + outgoing
  v1(:,:,ks-1) = v1(:,:,ks  ) - (v1(:,:,ks+1)-v1(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v1(:,:,ks-2) = v1(:,:,ks-1) - (v1(:,:,ks)-v1(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  v1(:,:,ks-1) = v1(:,:,ks-1) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
  v1(:,:,ks-2) = v1(:,:,ks-2) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
  v2(:,:,ks-1) = v2(:,:,ks  ) - (v2(:,:,ks+1)-v2(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v2(:,:,ks-2) = v2(:,:,ks-1) - (v2(:,:,ks)-v2(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  v2(:,:,ks-1) = v2(:,:,ks-1) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
  v2(:,:,ks-2) = v2(:,:,ks-2) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
  v3(:,:,ks-1) = v3(:,:,ks  ) - (v3(:,:,ks+1)-v3(:,:,ks)) * dx3(ks)/dx3(ks+1)
  v3(:,:,ks-2) = v3(:,:,ks-1) - (v3(:,:,ks)-v3(:,:,ks-1)) * dx3(ks-1)/dx3(ks)
  v3(:,:,ks-1) = v3(:,:,ks-1) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
  v3(:,:,ks-2) = v3(:,:,ks-2) * (0.5d0-sign(0.5d0,v3(:,:,ks)))
 case(9) x3_inner_vector ! Dirichlet
  v1(is:ie,js:je,ks-2:ks-1) = v10(is:ie,js:je,ks-2:ks-1)
  v2(is:ie,js:je,ks-2:ks-1) = v20(is:ie,js:je,ks-2:ks-1)
  v3(is:ie,js:je,ks-2:ks-1) = v30(is:ie,js:je,ks-2:ks-1)
 case default x3_inner_vector
  print *, "Error from x3 velocity inner boundary condition" ; stop
 end select x3_inner_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ks-2, ks-1
  do j = js, je
   do i = is, ie
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(i,j,ks)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
     e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do
! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! scalar values
 x3_outer_scalar: select case (bc3os)
 case(0) x3_outer_scalar ! periodic
  d (is:ie,js:je,ke+1:ke+2) = d (is:ie,js:je,ks:ks+1)
  p (is:ie,js:je,ke+1:ke+2) = p (is:ie,js:je,ks:ks+1)
  b1(is:ie,js:je,ke+1:ke+2) = b1(is:ie,js:je,ks:ks+1)
  b2(is:ie,js:je,ke+1:ke+2) = b2(is:ie,js:je,ks:ks+1)
  b3(is:ie,js:je,ke+1:ke+2) = b3(is:ie,js:je,ks:ks+1)
  phi(is:ie,js:je,ke+1:ke+2) = phi(is:ie,js:je,ks:ks+1)
 case(1) x3_outer_scalar ! reflective
  d (is:ie,js:je,ke+1:ke+2) = d (is:ie,js:je,ke:ke-1:-1)
  p (is:ie,js:je,ke+1:ke+2) = p (is:ie,js:je,ke:ke-1:-1)
  b1(is:ie,js:je,ke+1:ke+2) = b1(is:ie,js:je,ke:ke-1:-1)
  b2(is:ie,js:je,ke+1:ke+2) = b2(is:ie,js:je,ke:ke-1:-1)
  b3(is:ie,js:je,ke+1:ke+2) = b3(is:ie,js:je,ke:ke-1:-1)
 case(2:3) x3_outer_scalar ! outgoing/free
  d (is:ie,js:je,ke+1:ke+2) = spread(d (is:ie,js:je,ke),3,2)
  p (is:ie,js:je,ke+1:ke+2) = spread(p (is:ie,js:je,ke),3,2)
  b1(is:ie,js:je,ke+1:ke+2) = spread(b1(is:ie,js:je,ke),3,2) 
  b2(is:ie,js:je,ke+1:ke+2) = spread(b2(is:ie,js:je,ke),3,2) 
  b3(is:ie,js:je,ke+1:ke+2) = spread(b3(is:ie,js:je,ke),3,2)
 case(4) x3_outer_scalar ! linear
  d (:,:,ke+1) = d (:,:,ke  ) + (d (:,:,ke)-d (:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  d (:,:,ke+2) = d (:,:,ke+1) + (d (:,:,ke+1)-d (:,:,ke)) * dx3(ke+1)/dx3(ke)
  p (:,:,ke+1) = p (:,:,ke  ) + (p (:,:,ke)-p (:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  p (:,:,ke+2) = p (:,:,ke+1) + (p (:,:,ke+1)-p (:,:,ke)) * dx3(ke+1)/dx3(ke)
!$omp parallel do private (i,j,k)
  do k = ke+1, ke+2
   do j = js, je
    do i = is, ie
     if(d(i,j,k)<0d0)then
      d(i,j,k) = d(i,j,ke)
     end if
     if(p(i,j,k)<0d0)then
      p(i,j,k) = p(i,j,ke)
     end if
    end do
   end do
  end do
!$omp end parallel do
  b1(:,:,ke+1) = b1(:,:,ke  ) + (b1(:,:,ke)-b1(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  b1(:,:,ke+2) = b1(:,:,ke+1) + (b1(:,:,ke+1)-b1(:,:,ke)) * dx3(ke+1)/dx3(ke)
  b2(:,:,ke+1) = b2(:,:,ke  ) + (b2(:,:,ke)-b2(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  b2(:,:,ke+2) = b2(:,:,ke+1) + (b2(:,:,ke+1)-b2(:,:,ke)) * dx3(ke+1)/dx3(ke)
  b3(:,:,ke+1) = b3(:,:,ke  ) + (b3(:,:,ke)-b3(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  b3(:,:,ke+2) = b3(:,:,ke+1) + (b3(:,:,ke+1)-b3(:,:,ke)) * dx3(ke+1)/dx3(ke)
 case(9) x3_outer_scalar ! Dirichlet
  d (is:ie,js:je,ke+1:ke+2) = d0 (is:ie,js:je,ke+1:ke+2)
  p (is:ie,js:je,ke+1:ke+2) = p0 (is:ie,js:je,ke+1:ke+2)
  b1(is:ie,js:je,ke+1:ke+2) = b10(is:ie,js:je,ke+1:ke+2)
  b2(is:ie,js:je,ke+1:ke+2) = b20(is:ie,js:je,ke+1:ke+2)
  b3(is:ie,js:je,ke+1:ke+2) = b30(is:ie,js:je,ke+1:ke+2)
 case default x3_outer_scalar
  print *, "Error from x3 scalar outer boundary condition" ; stop
 end select x3_outer_scalar

! vector values
 x3_outer_vector: select case (bc3ov)
 case(0) x3_outer_vector ! periodic
  v1(is:ie,js:je,ke+1:ke+2) = v1(is:ie,js:je,ks:ks+1)
  v2(is:ie,js:je,ke+1:ke+2) = v2(is:ie,js:je,ks:ks+1)
  v3(is:ie,js:je,ke+1:ke+2) = v3(is:ie,js:je,ks:ks+1)
 case(1) x3_outer_vector ! reflective
  v1(is:ie,js:je,ke+1:ke+2) =  v1(is:ie,js:je,ke:ke-1:-1)
  v2(is:ie,js:je,ke+1:ke+2) =  v2(is:ie,js:je,ke:ke-1:-1)
  v3(is:ie,js:je,ke+1:ke+2) = -v3(is:ie,js:je,ke:ke-1:-1)
 case(2) x3_outer_vector ! outgoing
  allocate( plug(is:ie,js:je,1:1) )
  plug(:,:,1) = 0.5d0 + sign(0.5d0,v3(is:ie,js:je,ke))
  v1(is:ie,js:je,ke+1:ke+2) = spread(v1(is:ie,js:je,ke)*plug(:,:,1),3,2) 
  v2(is:ie,js:je,ke+1:ke+2) = spread(v2(is:ie,js:je,ke)*plug(:,:,1),3,2) 
  v3(is:ie,js:je,ke+1:ke+2) = spread(v3(is:ie,js:je,ke)*plug(:,:,1),3,2) 
  deallocate( plug )
 case(3) x3_outer_vector ! free
  v1(is:ie,js:je,ke+1:ke+2) = spread(v1(is:ie,js:je,ke),3,2)
  v2(is:ie,js:je,ke+1:ke+2) = spread(v2(is:ie,js:je,ke),3,2)
  v3(is:ie,js:je,ke+1:ke+2) = spread(v3(is:ie,js:je,ke),3,2)
 case(4) x3_outer_vector ! linear
  v1(:,:,ke+1) = v1(:,:,ke  ) + (v1(:,:,ke)-v1(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v1(:,:,ke+2) = v1(:,:,ke+1) + (v1(:,:,ke+1)-v1(:,:,ke)) * dx3(ke+1)/dx3(ke)
  v2(:,:,ke+1) = v2(:,:,ke  ) + (v2(:,:,ke)-v2(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v2(:,:,ke+2) = v2(:,:,ke+1) + (v2(:,:,ke+1)-v2(:,:,ke)) * dx3(ke+1)/dx3(ke)
  v3(:,:,ke+1) = v3(:,:,ke  ) + (v3(:,:,ke)-v3(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v3(:,:,ke+2) = v3(:,:,ke+1) + (v3(:,:,ke+1)-v3(:,:,ke)) * dx3(ke+1)/dx3(ke)
 case(5) x3_outer_vector ! linear + outgoing
  v1(:,:,ke+1) = v1(:,:,ke  ) + (v1(:,:,ke)-v1(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v1(:,:,ke+2) = v1(:,:,ke+1) + (v1(:,:,ke+1)-v1(:,:,ke)) * dx3(ke+1)/dx3(ke)
  v1(:,:,ke+1) = v1(:,:,ke+1) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
  v1(:,:,ke+2) = v1(:,:,ke+2) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
  v2(:,:,ke+1) = v2(:,:,ke  ) + (v2(:,:,ke)-v2(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v2(:,:,ke+2) = v2(:,:,ke+1) + (v2(:,:,ke+1)-v2(:,:,ke)) * dx3(ke+1)/dx3(ke)
  v2(:,:,ke+1) = v2(:,:,ke+1) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
  v2(:,:,ke+2) = v2(:,:,ke+2) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
  v3(:,:,ke+1) = v3(:,:,ke  ) + (v3(:,:,ke)-v3(:,:,ke-1)) * dx3(ke)/dx3(ke-1)
  v3(:,:,ke+2) = v3(:,:,ke+1) + (v3(:,:,ke+1)-v3(:,:,ke)) * dx3(ke+1)/dx3(ke)
  v3(:,:,ke+1) = v3(:,:,ke+1) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
  v3(:,:,ke+2) = v3(:,:,ke+2) * (0.5d0+sign(0.5d0,v3(:,:,ke)))
 case(9) x3_outer_vector ! Dirichlet
  v1(is:ie,js:je,ke+1:ke+2) = v10(is:ie,js:je,ke+1:ke+2)
  v2(is:ie,js:je,ke+1:ke+2) = v20(is:ie,js:je,ke+1:ke+2)
  v3(is:ie,js:je,ke+1:ke+2) = v30(is:ie,js:je,ke+1:ke+2)
 case default x3_outer_vector
  print *, "Error from x3 velocity outer boundary condition" ; stop
 end select x3_outer_vector

! set e and ptot
!$omp parallel do private(i,j,k)
 do k = ke+1, ke+2
  do j = js, je
   do i = is, ie
    ptot(i,j,k) = p(i,j,k) &
                + 0.5d0*( b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
    T(i,j,k) = T(i,j,ke)
    select case (eostype)
    case(0:1) ! without recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2) ! with recombination
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    e   (i,j,k) = eint(i,j,k) &
                + 0.5d0*( d(i,j,k)*&
                         ( v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2 )&
                         + b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2 )
   end do
  end do
 end do
!$omp end parallel do
end if

return
contains

 subroutine spcboundary
  
  implicit none

! Free boundaries for chemical composition
  do k = ks, ke
   do j = js-2, js-1
    do i = is, ie
     spc(1:spn,i,j,k) = spc(1:spn,i,js,k)
    end do
   end do
   do j = js, je
    do i = is-2, is-1
     spc(1:spn,i,j,k) = spc(1:spn,is,j,k)
    end do
    do i = ie+1, ie+2
     spc(1:spn,i,j,k) = spc(1:spn,ie,j,k)
    end do
   end do
   do j = je+1, je+2
    do i = is, ie
     spc(1:spn,i,j,k) = spc(1:spn,i,je,k)
    end do
   end do
  end do
  do j = js, je
   do i = is, ie
    do k = ks-2, ks-1
     spc(1:spn,i,j,k) = spc(1:spn,i,j,ks)
    end do
    do k = ke+1, ke+2
     spc(1:spn,i,j,k) = spc(1:spn,i,j,ke)
    end do
   end do
  end do

! temp
  if(bc1os==9)then
   do k = ks, ke
    do j = js, je
     do i = ie+1, ie+2
      spc(1:8,i,j,k) = spc0(1:8,i,j,k)
     end do
    end do
   end do
  end if
  return
 end subroutine spcboundary
 
end subroutine boundarycondition

end module boundary_mod
