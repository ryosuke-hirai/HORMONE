module boundary_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE BOUNDARYCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary condition

subroutine boundarycondition

 use settings
 use grid
 use physval
 use dirichlet_mod
 use pressure_mod
 use composition_mod
 use profiler_mod

 integer:: i,j,k
 real(8):: plug

!-------------------------------------------------------------------------

! Note:
!  0: Periodic b.c.
!  1: Reflective b.c.
!  2: Outgoing b.c.
!  3: Free b.c.
!  4: Linear Extrapolation
!  5: Linear Extrapolation but Outgoing (only for vectors)
!  9: Dirichlet b.c. (boundary values should be given elsewhere!)
! 10: Flux b.c. (flux values should be given elsewhere!)

 call start_clock(wtbnd)

!$omp parallel

! x1-direction ***********************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! If physical boundary condition (not MPI)
 if(is==is_global)then

  ! scalar values
  x1_inner_scalar: select case (bc1is)
  case(0) x1_inner_scalar ! periodic --------------------------------------
  ! When MPI is used, periodic BCs are already applied by the exchange
  ! and applying them here with is and ie will product wrong results.
  ! In serial, this is still necessary.
  if (is==is_global .and. ie==ie_global) then
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(is-2:is-1,j,k) = d(ie-1:ie,j,k)
      p(is-2:is-1,j,k) = p(ie-1:ie,j,k)
      if(mag_on)phi(is-2:is-1,j,k) = phi(ie-1:ie,j,k)
      if(compswitch>=2)spc(1:spn,is-2:is-1,j,k) = spc(1:spn,ie-1:ie,j,k)
    end do
    end do
  !$omp end do
  endif

  case(1) x1_inner_scalar ! reflective ------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(is-2,j,k) = d(is+1,j,k) ; d(is-1,j,k) = d(is,j,k)
      p(is-2,j,k) = p(is+1,j,k) ; p(is-1,j,k) = p(is,j,k)
      if(mag_on)then
      phi(is-2,j,k) = phi(is+1,j,k) ; phi(is-1,j,k) = phi(is,j,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,is-2,j,k) = spc(1:spn,is+1,j,k)
      spc(1:spn,is-1,j,k) = spc(1:spn,is  ,j,k)
      end if
    end do
    end do
  !$omp end do

  case(2:3) x1_inner_scalar ! outgoing/free -------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(is-2:is-1,j,k) = d(is,j,k)
      p(is-2:is-1,j,k) = p(is,j,k)
      if(mag_on)phi(is-2:is-1,j,k) = phi(is,j,k)
      if(compswitch>=2)then
      spc(1:spn,is-2,j,k) = spc(1:spn,is,j,k)
      spc(1:spn,is-1,j,k) = spc(1:spn,is,j,k)
      end if
    end do
    end do
  !$omp end do

  case(9) x1_inner_scalar ! Dirichlet -------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(is-2:is-1,j,k) = d0(is-2:is-1,j,k)
      p(is-2:is-1,j,k) = p0(is-2:is-1,j,k)
      if(compswitch>=2)then
      spc(1:spn,is-2:is-1,j,k) = spc0(1:spn,is-2:is-1,j,k)
      end if
    end do
    end do
  !$omp end do

  case(10) x1_inner_scalar ! Flux -----------------------------------------

  case default x1_inner_scalar ! Error ------------------------------------
    print *, "Error from x1 scalar inner boundary condition" ; stop
  end select x1_inner_scalar

  ! vector values
  x1_inner_vector: select case (bc1iv)
  case (0) x1_inner_vector ! periodic -------------------------------------
  if (is==is_global .and. ie==ie_global) then
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(is-2,j,k) = v1(ie-1,j,k) ; v1(is-1,j,k) = v1(ie,j,k)
      v2(is-2,j,k) = v2(ie-1,j,k) ; v2(is-1,j,k) = v2(ie,j,k)
      v3(is-2,j,k) = v3(ie-1,j,k) ; v3(is-1,j,k) = v3(ie,j,k)
      if(mag_on)then
      b1(is-2,j,k) = b1(ie-1,j,k) ; b1(is-1,j,k) = b1(ie,j,k)
      b2(is-2,j,k) = b2(ie-1,j,k) ; b2(is-1,j,k) = b2(ie,j,k)
      b3(is-2,j,k) = b3(ie-1,j,k) ; b3(is-1,j,k) = b3(ie,j,k)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x1_inner_vector ! reflective ------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(is-2,j,k) =-v1(is+1,j,k) ; v1(is-1,j,k) =-v1(is,j,k)
      v2(is-2,j,k) = v2(is+1,j,k) ; v2(is-1,j,k) = v2(is,j,k)
      v3(is-2,j,k) = v3(is+1,j,k) ; v3(is-1,j,k) = v3(is,j,k)
      if(mag_on)then
      b1(is-2,j,k) =-b1(is+1,j,k) ; b1(is-1,j,k) =-b1(is,j,k)
      b2(is-2,j,k) = b2(is+1,j,k) ; b2(is-1,j,k) = b2(is,j,k)
      b3(is-2,j,k) = b3(is+1,j,k) ; b3(is-1,j,k) = b3(is,j,k)
      end if
    end do
    end do
  !$omp end do

  case(2) x1_inner_vector ! outgoing --------------------------------------
  !$omp do private(plug,i,j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      plug = 0.5d0-sign(0.5d0,v1(is,j,k))
      do i = is-2, is-1
      v1(i,j,k) = min(0d0,v1(is,j,k))
      v2(i,j,k) = v2(is,j,k)
      v3(i,j,k) = v3(is,j,k)
      if(mag_on)then
        b1(i,j,k) = b1(is,j,k)
        b2(i,j,k) = b2(is,j,k)
        b3(i,j,k) = b3(is,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(3) x1_inner_vector ! free ------------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js, je
      do i = is-2, is-1
      v1(i,j,k) = v1(is,j,k)
      v2(i,j,k) = v2(is,j,k)
      v3(i,j,k) = v3(is,j,k)
      if(mag_on)then
        b1(i,j,k) = b1(is,j,k)
        b2(i,j,k) = b2(is,j,k)
        b3(i,j,k) = b3(is,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(9) x1_inner_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js, je
      do i = is-2, is-1
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x1_inner_vector ! Flux -----------------------------------------

  case default x1_inner_vector ! Error ------------------------------------
    print *, "Error from x1 velocity inner boundary condition" ; stop

  end select x1_inner_vector
end if

  ! Set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! ========================================================================

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! If physical boundary condition (not MPI)
 if(ie==ie_global)then

  ! scalar values
  x1_outer_scalar: select case (bc1os)
  case(0) x1_outer_scalar ! periodic --------------------------------------
  if (is==is_global .and. ie==ie_global) then
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(ie+1:ie+2,j,k) = d(is:is+1,j,k)
      p(ie+1:ie+2,j,k) = p(is:is+1,j,k)
      if(mag_on)phi(ie+1:ie+2,j,k) = phi(is:is+1,j,k)
      if(compswitch>=2)spc(1:spn,ie+1:ie+2,j,k) = spc(1:spn,is:is+1,j,k)
    end do
    end do
  !$omp end do
endif

  case(1) x1_outer_scalar ! reflective ------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(ie+1,j,k) = d(ie,j,k) ; d(ie+2,j,k) = d(ie-1,j,k)
      p(ie+1,j,k) = p(ie,j,k) ; p(ie+2,j,k) = p(ie-1,j,k)
      if(mag_on)then
      phi(ie+1,j,k) = phi(ie,j,k) ; phi(ie+2,j,k) = phi(ie-1,j,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,ie+1,j,k) = spc(1:spn,ie  ,j,k)
      spc(1:spn,ie+2,j,k) = spc(1:spn,ie-1,j,k)
      end if
    end do
    end do
  !$omp end do

  case(2:3) x1_outer_scalar ! outgoing/free -------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(ie+1:ie+2,j,k) = d(ie,j,k)
      p(ie+1:ie+2,j,k) = p(ie,j,k)
      if(mag_on)then
      phi(ie+1:ie+2,j,k) = phi(ie,j,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,ie+1,j,k) = spc(1:spn,ie,j,k)
      spc(1:spn,ie+2,j,k) = spc(1:spn,ie,j,k)
      end if
    end do
    end do
  !$omp end do

  case(4:5) x1_outer_scalar ! linear --------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      d(ie+1,j,k) = d(ie  ,j,k) + (d(ie,j,k)-d(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      d(ie+2,j,k) = d(ie+1,j,k) + (d(ie+1,j,k)-d(ie,j,k))*dx1(ie+1)/dx1(ie)
      p(ie+1,j,k) = p(ie  ,j,k) + (p(ie,j,k)-p(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      p(ie+2,j,k) = p(ie+1,j,k) + (p(ie+1,j,k)-p(ie,j,k))*dx1(ie+1)/dx1(ie)
      if(d(ie+2,j,k)<=0d0)then
      d(ie+1:ie+2,j,k) = d(ie,j,k)
      end if
      if(p(ie+2,j,k)<=0d0)then
      p(ie+1:ie+2,j,k) = p(ie,j,k)
      end if
      if(mag_on)then
      phi(ie+1:ie+2,j,k) = phi(ie,j,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,ie+1,j,k) = spc(1:spn,ie,j,k)
      spc(1:spn,ie+2,j,k) = spc(1:spn,ie,j,k)
      end if
    end do
    end do
  !$omp end do

  case(9) x1_outer_scalar ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js, je
      do i = ie+1, ie+2
      d(i,j,k) = d0(i,j,k)
      p(i,j,k) = p0(i,j,k)
      if(compswitch>=2)then
        spc(1:spn,i,j,k) = spc0(1:spn,i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x1_outer_scalar ! Flux -----------------------------------------

  case default x1_outer_scalar ! Error ------------------------------------
    print *, "Error from x1 scalar outer boundary condition" ; stop
  end select x1_outer_scalar

  ! vector values
  x1_outer_vector: select case (bc1ov)
  case(0) x1_outer_vector ! periodic --------------------------------------
  if (is==is_global .and. ie==ie_global) then
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(ie+1,j,k) = v1(is,j,k) ; v1(ie+2,j,k) = v1(is+1,j,k)
      v2(ie+1,j,k) = v2(is,j,k) ; v2(ie+2,j,k) = v2(is+1,j,k)
      v3(ie+1,j,k) = v3(is,j,k) ; v3(ie+2,j,k) = v3(is+1,j,k)
      if(mag_on)then
      b1(ie+1,j,k) = b1(is,j,k) ; b1(ie+2,j,k) = b1(is+1,j,k)
      b2(ie+1,j,k) = b2(is,j,k) ; b2(ie+2,j,k) = b2(is+1,j,k)
      b3(ie+1,j,k) = b3(is,j,k) ; b3(ie+2,j,k) = b3(is+1,j,k)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x1_outer_vector ! reflective ------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(ie+1,j,k) =-v1(ie,j,k) ; v1(ie+2,j,k) =-v1(ie-1,j,k)
      v2(ie+1,j,k) = v2(ie,j,k) ; v2(ie+2,j,k) = v2(ie-1,j,k)
      v3(ie+1,j,k) = v3(ie,j,k) ; v3(ie+2,j,k) = v3(ie-1,j,k)
      if(mag_on)then
      b1(ie+1,j,k) =-b1(ie,j,k) ; b1(ie+2,j,k) =-b1(ie-1,j,k)
      b2(ie+1,j,k) = b2(ie,j,k) ; b2(ie+2,j,k) = b2(ie-1,j,k)
      b3(ie+1,j,k) = b3(ie,j,k) ; b3(ie+2,j,k) = b3(ie-1,j,k)
      end if
    end do
    end do
  !$omp end do

  case(2) x1_outer_vector ! outgoing --------------------------------------
  !$omp do private(plug,i,j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      plug = 0.5d0+sign(0.5d0,v1(ie,j,k))
      do i = ie+1, ie+2
      v1(i,j,k) = v1(ie,j,k)*plug
      v2(i,j,k) = v2(ie,j,k)
      v3(i,j,k) = v3(ie,j,k)
      if(mag_on)then
        b1(i,j,k) = b1(ie,j,k)
        b2(i,j,k) = b2(ie,j,k)
        b3(i,j,k) = b3(ie,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(3) x1_outer_vector ! free ------------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(ie+1:ie+2,j,k) = v1(ie,j,k)
      v2(ie+1:ie+2,j,k) = v2(ie,j,k)
      v3(ie+1:ie+2,j,k) = v3(ie,j,k)
      if(mag_on)then
      b1(ie+1:ie+2,j,k) = b1(ie,j,k)
      b2(ie+1:ie+2,j,k) = b2(ie,j,k)
      b3(ie+1:ie+2,j,k) = b3(ie,j,k)
      end if
    end do
    end do
  !$omp end do

  case(4) x1_outer_vector ! linear ----------------------------------------
  !$omp do private(j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      v1(ie+1,j,k) = v1(ie  ,j,k) + (v1(ie,j,k)-v1(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v1(ie+2,j,k) = v1(ie+1,j,k) + (v1(ie+1,j,k)-v1(ie,j,k))*dx1(ie+1)/dx1(ie)
      v2(ie+1,j,k) = v2(ie  ,j,k) + (v2(ie,j,k)-v2(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v2(ie+2,j,k) = v2(ie+1,j,k) + (v2(ie+1,j,k)-v2(ie,j,k))*dx1(ie+1)/dx1(ie)
      v3(ie+1,j,k) = v3(ie  ,j,k) + (v3(ie,j,k)-v3(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v3(ie+2,j,k) = v3(ie+1,j,k) + (v3(ie+1,j,k)-v3(ie,j,k))*dx1(ie+1)/dx1(ie)
      if(mag_on)then
      b1(ie+1,j,k) = b1(ie  ,j,k) + (b1(ie,j,k)-b1(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b1(ie+2,j,k) = b1(ie+1,j,k) + (b1(ie+1,j,k)-b1(ie,j,k))*dx1(ie+1)/dx1(ie)
      b2(ie+1,j,k) = b2(ie  ,j,k) + (b2(ie,j,k)-b2(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b2(ie+2,j,k) = b2(ie+1,j,k) + (b2(ie+1,j,k)-b2(ie,j,k))*dx1(ie+1)/dx1(ie)
      b3(ie+1,j,k) = b3(ie  ,j,k) + (b3(ie,j,k)-b3(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b3(ie+2,j,k) = b3(ie+1,j,k) + (b3(ie+1,j,k)-b3(ie,j,k))*dx1(ie+1)/dx1(ie)
      end if
    end do
    end do
  !$omp end do

  case(5) x1_outer_vector ! linear + outgoing -----------------------------
  !$omp do private(plug,j,k) collapse(2)
    do k = ks, ke
    do j = js, je
      plug = 0.5d0+sign(0.5d0,v1(ie,j,k))
      v1(ie+1,j,k) = v1(ie  ,j,k) + (v1(ie,j,k)-v1(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v1(ie+2,j,k) = v1(ie+1,j,k) + (v1(ie+1,j,k)-v1(ie,j,k))*dx1(ie+1)/dx1(ie)
      v1(ie+1:ie+2,j,k) = v1(ie+1:ie+2,j,k) * plug
      v2(ie+1,j,k) = v2(ie  ,j,k) + (v2(ie,j,k)-v2(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v2(ie+2,j,k) = v2(ie+1,j,k) + (v2(ie+1,j,k)-v2(ie,j,k))*dx1(ie+1)/dx1(ie)
      v2(ie+1:ie+2,j,k) = v2(ie+1:ie+2,j,k)
      v3(ie+1,j,k) = v3(ie  ,j,k) + (v3(ie,j,k)-v3(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      v3(ie+2,j,k) = v3(ie+1,j,k) + (v3(ie+1,j,k)-v3(ie,j,k))*dx1(ie+1)/dx1(ie)
      v3(ie+1:ie+2,j,k) = v3(ie+1:ie+2,j,k)
      if(mag_on)then
      b1(ie+1,j,k) = b1(ie  ,j,k) + (b1(ie,j,k)-b1(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b1(ie+2,j,k) = b1(ie+1,j,k) + (b1(ie+1,j,k)-b1(ie,j,k))*dx1(ie+1)/dx1(ie)
      b1(ie+1:ie+2,j,k) = b1(ie+1:ie+2,j,k)
      b2(ie+1,j,k) = b2(ie  ,j,k) + (b2(ie,j,k)-b2(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b2(ie+2,j,k) = b2(ie+1,j,k) + (b2(ie+1,j,k)-b2(ie,j,k))*dx1(ie+1)/dx1(ie)
      b2(ie+1:ie+2,j,k) = b2(ie+1:ie+2,j,k)
      b3(ie+1,j,k) = b3(ie  ,j,k) + (b3(ie,j,k)-b3(ie-1,j,k))*dx1(ie)/dx1(ie-1)
      b3(ie+2,j,k) = b3(ie+1,j,k) + (b3(ie+1,j,k)-b3(ie,j,k))*dx1(ie+1)/dx1(ie)
      b3(ie+1:ie+2,j,k) = b3(ie+1:ie+2,j,k)
      end if
    end do
    end do
  !$omp end do

  case(9) x1_outer_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js, je
      do i = ie+1, ie+2
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x1_outer_vector ! Flux -----------------------------------------

  case default x1_outer_vector ! Error ------------------------------------
    print *, "Error from x1 velocity outer boundary condition" ; stop

  end select x1_outer_vector
end if

  ! set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! =======================================================================

! x2-direction ***********************************************************
! If physical boundary condition (not MPI)
if(je>js)then ! TODO: In MPI, this can be false even when x2 is active
! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(js==js_global)then
  ! scalar values
  x2_inner_scalar: select case (bc2is)
  case(0) x2_inner_scalar ! periodic --------------------------------------
  if (js==js_global .and. je==je_global) then
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      d(i,js-2,k) = d(i,je-1,k); d(i,js-1,k) = d(i,je,k)
      p(i,js-2,k) = p(i,je-1,k); p(i,js-1,k) = p(i,je,k)
      if(mag_on)then
      phi(i,js-2,k) = phi(i,je-1,k); phi(i,js-1,k) = phi(i,je,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,js-2,k) = spc(1:spn,i,je-1,k)
      spc(1:spn,i,js-1,k) = spc(1:spn,i,je  ,k)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x2_inner_scalar ! reflective ------------------------------------
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      d(i,js-2,k) = d(i,js+1,k) ; d(i,js-1,k) = d(i,js,k)
      p(i,js-2,k) = p(i,js+1,k) ; p(i,js-1,k) = p(i,js,k)
      if(mag_on)then
      phi(i,js-2,k) = phi(i,js+1,k) ; phi(i,js-1,k) = phi(i,js,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,js-2,k) = spc(1:spn,i,js+1,k)
      spc(1:spn,i,js-1,k) = spc(1:spn,i,js  ,k)
      end if
    end do
    end do
  !$omp end do

  case(2:3) x2_inner_scalar ! outgoing/free -------------------------------
  !$omp do private (i,j,k) collapse(3)
    do k = ks, ke
    do j = js-2, js-1
      do i = is, ie
      d(i,j,k) = d(i,js,k)
      p(i,j,k) = p(i,js,k)
      if(mag_on)phi(i,j,k) = phi(i,js,k)
      if(compswitch>=2)spc(1:spn,i,j,k) = spc(1:spn,i,js,k)
      end do
    end do
    end do
  !$omp end do

  case(9) x2_inner_scalar ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js-2, js-1
      do i = is, ie
      d(i,j,k) = d0(i,j,k)
      p(i,j,k) = p0(i,j,k)
      if(compswitch>=2)then
        spc(1:spn,i,j,k) = spc0(1:spn,i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x2_inner_scalar ! Flux ----------------------------------------

  case default x2_inner_scalar ! Error -----------------------------------
    print *, "Error from x2 scalar inner boundary condition" ; stop
  end select x2_inner_scalar

  ! vector values
  x2_inner_vector: select case (bc2iv)
  case(0) x2_inner_vector ! periodic -------------------------------------
  if (js==js_global .and. je==je_global) then
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      v1(i,js-2,k) = v1(i,je-1,k) ; v1(i,js-1,k) = v1(i,je,k)
      v2(i,js-2,k) = v2(i,je-1,k) ; v2(i,js-1,k) = v2(i,je,k)
      v3(i,js-2,k) = v3(i,je-1,k) ; v3(i,js-1,k) = v3(i,je,k)
      if(mag_on)then
      b1(i,js-2,k) = b1(i,je-1,k) ; b1(i,js-1,k) = b1(i,je,k)
      b2(i,js-2,k) = b2(i,je-1,k) ; b2(i,js-1,k) = b2(i,je,k)
      b3(i,js-2,k) = b3(i,je-1,k) ; b3(i,js-1,k) = b3(i,je,k)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x2_inner_vector ! reflective -----------------------------------
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      v1(i,js-2,k) = v1(i,js+1,k) ; v1(i,js-1,k) = v1(i,js,k)
      v2(i,js-2,k) =-v2(i,js+1,k) ; v2(i,js-1,k) =-v2(i,js,k)
      v3(i,js-2,k) = v3(i,js+1,k) ; v3(i,js-1,k) = v3(i,js,k)
      if(mag_on)then
      b1(i,js-2,k) = b1(i,js+1,k) ; b1(i,js-1,k) = b1(i,js,k)
      b2(i,js-2,k) =-b2(i,js+1,k) ; b2(i,js-1,k) =-b2(i,js,k)
      b3(i,js-2,k) = b3(i,js+1,k) ; b3(i,js-1,k) = b3(i,js,k)
      end if
    end do
    end do
  !$omp end do

  case(2) x2_inner_vector ! outgoing -------------------------------------
  !$omp do private(plug,i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      plug = 0.5d0-sign(0.5d0,v2(i,js,k))
      v1(i,js-2:js-1,k) = v1(i,js,k)
      v2(i,js-2:js-1,k) = v2(i,js,k)*plug
      v3(i,js-2:js-1,k) = v3(i,js,k)
      if(mag_on)then
      b1(i,js-2:js-1,k) = b1(i,js,k)
      b2(i,js-2:js-1,k) = b2(i,js,k)
      b3(i,js-2:js-1,k) = b3(i,js,k)
      end if
    end do
    end do
  !$omp end do

  case(3) x2_inner_vector ! free ------------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js-2, js-1
      do i = is, ie
      v1(i,j,k) = v1(i,js,k)
      v2(i,j,k) = v2(i,js,k)
      v3(i,j,k) = v3(i,js,k)
      if(mag_on)then
        b1(i,j,k) = b1(i,js,k)
        b2(i,j,k) = b2(i,js,k)
        b3(i,j,k) = b3(i,js,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(9) x2_inner_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = js-2, js-1
      do i = is, ie
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do

  case(10) x2_inner_vector ! Flux -----------------------------------------

  case default x2_inner_vector ! Error ------------------------------------
    print *, "Error from x2 velocity inner boundary condition" ; stop
  end select x2_inner_vector
end if

  ! set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! ========================================================================

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! If physical boundary condition (not MPI)
 if(je==je_global)then
  ! scalar values
  x2_outer_scalar: select case (bc2os)
  case(0) x2_outer_scalar ! periodic --------------------------------------
  if (js==js_global .and. je==je_global) then
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      d(i,je+1:je+2,k) = d(i,js:js+1,k)
      p(i,je+1:je+2,k) = p(i,js:js+1,k)
      if(mag_on)phi(i,je+1:je+2,k) = phi(i,js:js+1,k)
      if(compswitch>=2)spc(1:spn,i,je+1:je+2,k) = spc(1:spn,i,js:js+1,k)
    end do
    end do
  !$omp end do
  endif

  case(1) x2_outer_scalar ! reflective ------------------------------------
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      d(i,je+1,k) = d(i,je,k) ; d(i,je+2,k) = d(i,je-1,k)
      p(i,je+1,k) = p(i,je,k) ; p(i,je+2,k) = p(i,je-1,k)
      if(mag_on)then
      phi(i,je+1,k) = phi(i,je,k) ; phi(i,je+2,k) = phi(i,je-1,k)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,je+1,k) = spc(1:spn,i,je  ,k)
      spc(1:spn,i,je+2,k) = spc(1:spn,i,je-1,k)
      end if
    end do
    end do
  !$omp end do

  case(2:3) x2_outer_scalar ! outgoing/free -------------------------------
  !$omp do private(i,j,k) collapse(3)
    do i = is, ie
    do j = je+1, je+2
      do k = ks, ke
      d(i,j,k) = d(i,je,k)
      p(i,j,k) = p(i,je,k)
      if(mag_on)phi(i,j,k) = phi(i,je,k)
      if(compswitch>=2)spc(1:spn,i,j,k) = spc(1:spn,i,je,k)
      end do
    end do
    end do
  !$omp end do

  case(9) x2_outer_scalar ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = je+1, je+2
      do i = is, ie
      d(i,j,k) = d0(i,j,k)
      p(i,j,k) = p0(i,j,k)
      if(compswitch>=2)spc(1:spn,i,j,k) = spc0(1:spn,i,j,k)
      end do
    end do
    end do
  !$omp end do

  case(10) x2_outer_scalar ! Flux -----------------------------------------

  case default x2_outer_scalar ! Error ------------------------------------
    print *, "Error from x2 scalar outer boundary condition" ; stop
  end select x2_outer_scalar

  ! vector values
  x2_outer_vector: select case (bc2ov)
  case(0) x2_outer_vector ! periodic --------------------------------------
  if (js==js_global .and. je==je_global) then
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      v1(i,je+1:je+2,k) = v1(i,js:js+1,k)
      v2(i,je+1:je+2,k) = v2(i,js:js+1,k)
      v3(i,je+1:je+2,k) = v3(i,js:js+1,k)
      if(mag_on)then
      b1(i,je+1:je+2,k) = b1(i,js:js+1,k)
      b2(i,je+1:je+2,k) = b2(i,js:js+1,k)
      b3(i,je+1:je+2,k) = b3(i,js:js+1,k)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x2_outer_vector ! reflective ------------------------------------
  !$omp do private(i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      v1(i,je+1,k) = v1(i,je,k) ; v1(i,je+2,k) = v1(i,je-1,k)
      v2(i,je+1,k) =-v2(i,je,k) ; v2(i,je+2,k) =-v2(i,je-1,k)
      v3(i,je+1,k) = v3(i,je,k) ; v3(i,je+2,k) = v3(i,je-1,k)
      if(mag_on)then
      b1(i,je+1,k) = b1(i,je,k) ; b1(i,je+2,k) = b1(i,je-1,k)
      b2(i,je+1,k) =-b2(i,je,k) ; b2(i,je+2,k) =-b2(i,je-1,k)
      b3(i,je+1,k) = b3(i,je,k) ; b3(i,je+2,k) = b3(i,je-1,k)
      end if
    end do
    end do
  !$omp end do

  case(2) x2_outer_vector ! outgoing --------------------------------------
  !$omp do private(plug,i,k) collapse(2)
    do k = ks, ke
    do i = is, ie
      plug = 0.5d0+sign(0.5d0,v2(i,je,k))
      v1(i,je+1:je+2,k) = v1(i,je,k)
      v2(i,je+1:je+2,k) = v2(i,je,k)*plug
      v3(i,je+1:je+2,k) = v3(i,je,k)
      if(mag_on)then
      b1(i,je+1:je+2,k) = b1(i,je,k)
      b2(i,je+1:je+2,k) = b2(i,je,k)
      b3(i,je+1:je+2,k) = b3(i,je,k)
      end if
    end do
    end do
  !$omp end do

  case(3) x2_outer_vector ! free ------------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do i = is, ie
    do j = je+1, je+2
      do k = ks, ke
      v1(i,j,k) = v1(i,je,k) ; v2(i,j,k) = v2(i,je,k) ; v3(i,j,k) = v3(i,je,k)
      if(mag_on)then
        b1(i,j,k) = b1(i,je,k) ; b2(i,j,k) = b2(i,je,k) ; b3(i,j,k) = b3(i,je,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(9) x2_outer_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks, ke
    do j = je+1, je+2
      do i = is, ie
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do

  case(10) x2_outer_vector ! Flux -----------------------------------------

  case default x2_outer_vector ! Error ------------------------------------
    print *, "Error from x2 velocity outer boundary condition" ; stop
  end select x2_outer_vector
end if
  ! set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! ========================================================================
end if

! x3-direction ***********************************************************
if(ke>ks)then ! TODO: In MPI, this can be false even when x3 is active
! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! If physical boundary condition (not MPI)
 if(ks==ks_global)then
  ! scalar values
  x3_inner_scalar: select case (bc3is)
  case(0) x3_inner_scalar ! periodic --------------------------------------
  if (ks==ks_global .and. ke==ke_global) then
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ks-2) = d(i,j,ke-1) ; d(i,j,ks-1) = d(i,j,ke)
      p(i,j,ks-2) = p(i,j,ke-1) ; p(i,j,ks-1) = p(i,j,ke)
      if(mag_on)then
      phi(i,j,ks-2) = phi(i,j,ke-1) ; phi(i,j,ks-1) = phi(i,j,ke)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,j,ks-2) = spc(1:spn,i,j,ke-1)
      spc(1:spn,i,j,ks-1) = spc(1:spn,i,j,ke  )
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x3_inner_scalar ! reflective ------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ks-2) = d(i,j,ks+1) ; d(i,j,ks-1) = d(i,j,ks)
      p(i,j,ks-2) = p(i,j,ks+1) ; p(i,j,ks-1) = p(i,j,ks)
      if(mag_on)then
      phi(i,j,ks-2) = phi(i,j,ks+1) ; phi(i,j,ks-1) = phi(i,j,ks)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,j,ks-2) = spc(1:spn,i,j,ks+1)
      spc(1:spn,i,j,ks-1) = spc(1:spn,i,j,ks  )
      end if
    end do
    end do
  !$omp end do

 case(2:3) x3_inner_scalar ! outgoing/free -------------------------------
!$omp do private(i,j) collapse(2)
  do j = js, je
   do i = is, ie
    d(i,j,ks-2:ks-1) = d(i,j,ks)
    p(i,j,ks-2:ks-1) = p(i,j,ks)
    if(mag_on)phi(i,j,ks-2:ks-1) = phi(i,j,ks)
    if(compswitch>=2)then
     spc(1:spn,i,j,ks-2) = spc(1:spn,i,j,ks)
     spc(1:spn,i,j,ks-1) = spc(1:spn,i,j,ks)
    end if
   end do
  end do
!$omp end do

  case(4:5) x3_inner_scalar ! linear --------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ks-1) = d(i,j,ks  ) - (d(i,j,ks+1)-d(i,j,ks))*dx3(ks)/dx3(ks+1)
      d(i,j,ks-2) = d(i,j,ks-1) - (d(i,j,ks)-d(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      p(i,j,ks-1) = p(i,j,ks  ) - (p(i,j,ks+1)-p(i,j,ks))*dx3(ks)/dx3(ks+1)
      p(i,j,ks-2) = p(i,j,ks-1) - (p(i,j,ks)-p(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      if(d(i,j,ks-2)<=0d0)then
      d(i,j,ks-2:ks-1) = d(i,j,ks)
      end if
      if(p(i,j,ks-2)<=0d0)then
      p(i,j,ks-2:ks-1) = p(i,j,ks)
      end if
      if(mag_on)phi(i,j,ks-2:ks-1) = phi(i,j,ks)
      if(compswitch>=2)then
      spc(1:spn,i,j,ks-2) = spc(1:spn,i,j,ks)
      spc(1:spn,i,j,ks-1) = spc(1:spn,i,j,ks)
      end if
    end do
    end do
  !$omp end do

  case(9) x3_inner_scalar ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks-2, ks-1
    do j = js, je
      do i = is, ie
      d(i,j,k) = d0(i,j,k)
      p(i,j,k) = p0(i,j,k)
      if(compswitch>=2) spc(1:spn,i,j,k) = spc0(1:spn,i,j,k)
      end do
    end do
    end do
  !$omp end do

  case(10) x3_inner_scalar ! Flux -----------------------------------------

  case default x3_inner_scalar ! Error ------------------------------------
    print *, "Error from x3 scalar inner boundary condition" ; stop
  end select x3_inner_scalar

  ! vector values
  x3_inner_vector: select case (bc3iv)
  case(0) x3_inner_vector ! periodic --------------------------------------
  if (ks==ks_global .and. ke==ke_global) then
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ks-2:ks-1) = v1(i,j,ke-1:ke)
      v2(i,j,ks-2:ks-1) = v2(i,j,ke-1:ke)
      v3(i,j,ks-2:ks-1) = v3(i,j,ke-1:ke)
      if(mag_on)then
      b1(i,j,ks-2:ks-1) = b1(i,j,ke-1:ke)
      b2(i,j,ks-2:ks-1) = b2(i,j,ke-1:ke)
      b3(i,j,ks-2:ks-1) = b3(i,j,ke-1:ke)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x3_inner_vector ! reflective ------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ks-2) = v1(i,j,ks+1) ; v1(i,j,ks-1) = v1(i,j,ks)
      v2(i,j,ks-2) = v2(i,j,ks+1) ; v2(i,j,ks-1) = v2(i,j,ks)
      v3(i,j,ks-2) =-v3(i,j,ks+1) ; v3(i,j,ks-1) =-v3(i,j,ks)
      if(mag_on)then
      b1(i,j,ks-2) = b1(i,j,ks+1) ; b1(i,j,ks-1) = b1(i,j,ks)
      b2(i,j,ks-2) = b2(i,j,ks+1) ; b2(i,j,ks-1) = b2(i,j,ks)
      b3(i,j,ks-2) =-b3(i,j,ks+1) ; b3(i,j,ks-1) =-b3(i,j,ks)
      end if
    end do
    end do
  !$omp end do

  case(2) x3_inner_vector ! outgoing --------------------------------------
  !$omp do private(plug,i,j) collapse(2)
    do j = js, je
    do i = is, ie
      plug = 0.5d0-sign(0.5d0,v3(i,j,ks))
      v1(i,j,ks-2:ks-1) = v1(i,j,ks)
      v2(i,j,ks-2:ks-1) = v2(i,j,ks)
      v3(i,j,ks-2:ks-1) = v3(i,j,ks)*plug
      if(mag_on)then
      b1(i,j,ks-2:ks-1) = b1(i,j,ks)
      b2(i,j,ks-2:ks-1) = b2(i,j,ks)
      b3(i,j,ks-2:ks-1) = b3(i,j,ks)
      end if
    end do
    end do
  !$omp end do

  case(3) x3_inner_vector ! free ------------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ks-2:ks-1) = v1(i,j,ks)
      v2(i,j,ks-2:ks-1) = v2(i,j,ks)
      v3(i,j,ks-2:ks-1) = v3(i,j,ks)
      if(mag_on)then
      b1(i,j,ks-2:ks-1) = b1(i,j,ks)
      b2(i,j,ks-2:ks-1) = b2(i,j,ks)
      b3(i,j,ks-2:ks-1) = b3(i,j,ks)
      end if
    end do
    end do
  !$omp end do

  case(4) x3_inner_vector ! linear ----------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ks-1) = v1(i,j,ks  ) - (v1(i,j,ks+1)-v1(i,j,ks))*dx3(ks)/dx3(ks+1)
      v1(i,j,ks-2) = v1(i,j,ks-1) - (v1(i,j,ks)-v1(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      v2(i,j,ks-1) = v2(i,j,ks  ) - (v2(i,j,ks+1)-v2(i,j,ks))*dx3(ks)/dx3(ks+1)
      v3(i,j,ks-2) = v3(i,j,ks-1) - (v2(i,j,ks)-v2(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      v3(i,j,ks-1) = v3(i,j,ks  ) - (v3(i,j,ks+1)-v3(i,j,ks))*dx3(ks)/dx3(ks+1)
      v3(i,j,ks-2) = v3(i,j,ks-1) - (v3(i,j,ks)-v3(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      if(mag_on)then
      b1(i,j,ks-1) = b1(i,j,ks  ) - (b1(i,j,ks+1)-b1(i,j,ks))*dx3(ks)/dx3(ks+1)
      b1(i,j,ks-2) = b1(i,j,ks-1) - (b1(i,j,ks)-b1(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      b2(i,j,ks-1) = b2(i,j,ks  ) - (b2(i,j,ks+1)-b2(i,j,ks))*dx3(ks)/dx3(ks+1)
      b3(i,j,ks-2) = b3(i,j,ks-1) - (b2(i,j,ks)-b2(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      b3(i,j,ks-1) = b3(i,j,ks  ) - (b3(i,j,ks+1)-b3(i,j,ks))*dx3(ks)/dx3(ks+1)
      b3(i,j,ks-2) = b3(i,j,ks-1) - (b3(i,j,ks)-b3(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      end if
    end do
    end do
  !$omp end do

  case(5) x3_inner_vector ! linear + outgoing -----------------------------
  !$omp do private(plug,i,j) collapse(2)
    do j = js, je
    do i = is, ie
      plug = 0.5d0-sign(0.5d0,v3(i,j,ks))
      v1(i,j,ks-1) = v1(i,j,ks  ) - (v1(i,j,ks+1)-v1(i,j,ks))*dx3(ks)/dx3(ks+1)
      v1(i,j,ks-2) = v1(i,j,ks-1) - (v1(i,j,ks)-v1(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      v2(i,j,ks-1) = v2(i,j,ks  ) - (v2(i,j,ks+1)-v2(i,j,ks))*dx3(ks)/dx3(ks+1)
      v3(i,j,ks-2) = v3(i,j,ks-1) - (v2(i,j,ks)-v2(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      v3(i,j,ks-1) = v3(i,j,ks  ) - (v3(i,j,ks+1)-v3(i,j,ks))*dx3(ks)/dx3(ks+1)
      v3(i,j,ks-2) = v3(i,j,ks-1) - (v3(i,j,ks)-v3(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      v1(i,j,ks-2:ks-1) = v1(i,j,ks-2:ks-1)
      v2(i,j,ks-2:ks-1) = v2(i,j,ks-2:ks-1)
      v3(i,j,ks-2:ks-1) = v3(i,j,ks-2:ks-1) * plug
      if(mag_on)then
      b1(i,j,ks-1) = b1(i,j,ks  ) - (b1(i,j,ks+1)-b1(i,j,ks))*dx3(ks)/dx3(ks+1)
      b1(i,j,ks-2) = b1(i,j,ks-1) - (b1(i,j,ks)-b1(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      b2(i,j,ks-1) = b2(i,j,ks  ) - (b2(i,j,ks+1)-b2(i,j,ks))*dx3(ks)/dx3(ks+1)
      b3(i,j,ks-2) = b3(i,j,ks-1) - (b2(i,j,ks)-b2(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      b3(i,j,ks-1) = b3(i,j,ks  ) - (b3(i,j,ks+1)-b3(i,j,ks))*dx3(ks)/dx3(ks+1)
      b3(i,j,ks-2) = b3(i,j,ks-1) - (b3(i,j,ks)-b3(i,j,ks-1))*dx3(ks-1)/dx3(ks)
      end if
    end do
    end do
  !$omp end do

  case(9) x3_inner_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ks-2, ks-1
    do j = js, je
      do i = is, ie
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do

  case(10) x3_inner_vector ! Flux -----------------------------------------

  case default x3_inner_vector ! Error ------------------------------------
    print *, "Error from x3 velocity inner boundary condition" ; stop
  end select x3_inner_vector
end if

  ! set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! ========================================================================

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! If physical boundary condition (not MPI)
 if(ke==ke_global)then
  ! scalar values
  x3_outer_scalar: select case (bc3os)
  case(0) x3_outer_scalar ! periodic --------------------------------------
  if (ks==ks_global .and. ke==ke_global) then
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ke+1:ke+2) = d(i,j,ks:ks+1)
      p(i,j,ke+1:ke+2) = p(i,j,ks:ks+1)
      if(mag_on)phi(i,j,ke+1:ke+2) = phi(i,j,ks:ks+1)
      if(compswitch>=2)spc(1:spn,i,j,ke+1:ke+2) = spc(1:spn,i,j,ks:ks+1)
    end do
    end do
  !$omp end do
  endif

  case(1) x3_outer_scalar ! reflective ------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ke+1) = d(i,j,ke) ; d(i,j,ke+2) = d(i,j,ke-1)
      p(i,j,ke+1) = p(i,j,ke) ; p(i,j,ke+2) = p(i,j,ke-1)
      if(mag_on)then
      phi(i,j,ke+1) = phi(i,j,ke) ; phi(i,j,ke+2) = phi(i,j,ke-1)
      end if
      if(compswitch>=2)then
      spc(1:spn,i,j,ke+1) = spc(1:spn,i,j,ke  )
      spc(1:spn,i,j,ke+2) = spc(1:spn,i,j,ke-1)
      end if
    end do
    end do
  !$omp end do

  case(2:3) x3_outer_scalar ! outgoing/free -------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ke+1:ke+2) = d(i,j,ke)
      p(i,j,ke+1:ke+2) = p(i,j,ke)
      if(mag_on)phi(i,j,ke+1:ke+2) = phi(i,j,ke)
      if(compswitch>=2)then
      spc(1:spn,i,j,ke+1) = spc(1:spn,i,j,ke)
      spc(1:spn,i,j,ke+2) = spc(1:spn,i,j,ke)
      end if
    end do
    end do
  !$omp end do

  case(4) x3_outer_scalar ! linear ----------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      d(i,j,ke+1) = d(i,j,ke  ) + (d(i,j,ke)-d(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      d(i,j,ke+2) = d(i,j,ke+1) + (d(i,j,ke+1)-d(i,j,ke))*dx3(ke+1)/dx3(ke)
      p(i,j,ke+1) = p(i,j,ke  ) + (p(i,j,ke)-p(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      p(i,j,ke+2) = p(i,j,ke+1) + (p(i,j,ke+1)-p(i,j,ke))*dx3(ke+1)/dx3(ke)
      if(d(i,j,ke+2)<=0d0)then
      d(i,j,ke+1:ke+2) = d(i,j,ke)
      end if
      if(p(i,j,ke+2)<=0d0)then
      p(i,j,ke+1:ke+2) = p(i,j,ke)
      end if
      if(mag_on)phi(i,j,ke+1:ke+2) = phi(i,j,ke)
      if(compswitch>=2)then
      spc(1:spn,i,j,ke+1) = spc(1:spn,i,j,ke)
      spc(1:spn,i,j,ke+2) = spc(1:spn,i,j,ke)
      end if
    end do
    end do
  !$omp end do

  case(9) x3_outer_scalar ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ke+1, ke+2
    do j = js, je
      do i = is, ie
      d(i,j,k) = d0(i,j,k)
      p(i,j,k) = p0(i,j,k)
      if(compswitch>=2)then
        spc(1:spn,i,j,k) = spc0(1:spn,i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x3_outer_scalar ! Flux -----------------------------------------

  case default x3_outer_scalar ! Error ------------------------------------
    print *, "Error from x3 scalar outer boundary condition" ; stop
  end select x3_outer_scalar

  ! vector values
  x3_outer_vector: select case (bc3ov)
  case(0) x3_outer_vector ! periodic --------------------------------------
  if (ks==ks_global .and. ke==ke_global) then
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ke+1:ke+2) = v1(i,j,ks:ks+1)
      v2(i,j,ke+1:ke+2) = v2(i,j,ks:ks+1)
      v3(i,j,ke+1:ke+2) = v3(i,j,ks:ks+1)
      if(mag_on)then
      b1(i,j,ke+1:ke+2) = b1(i,j,ks:ks+1)
      b2(i,j,ke+1:ke+2) = b2(i,j,ks:ks+1)
      b3(i,j,ke+1:ke+2) = b3(i,j,ks:ks+1)
      end if
    end do
    end do
  !$omp end do
  endif

  case(1) x3_outer_vector ! reflective ------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ke+1) = v1(i,j,ke) ; v1(i,j,ke+2) = v1(i,j,ke-1)
      v2(i,j,ke+1) = v2(i,j,ke) ; v2(i,j,ke+2) = v2(i,j,ke-1)
      v3(i,j,ke+1) =-v3(i,j,ke) ; v3(i,j,ke+2) =-v3(i,j,ke-1)
      if(mag_on)then
      b1(i,j,ke+1) = b1(i,j,ke) ; b1(i,j,ke+2) = b1(i,j,ke-1)
      b2(i,j,ke+1) = b2(i,j,ke) ; b2(i,j,ke+2) = b2(i,j,ke-1)
      b3(i,j,ke+1) =-b3(i,j,ke) ; b3(i,j,ke+2) =-b3(i,j,ke-1)
      end if
    end do
    end do
  !$omp end do

  case(2) x3_outer_vector ! outgoing --------------------------------------
  !$omp do private(plug,i,j) collapse(2)
    do j = js, je
    do i = is, ie
      plug = 0.5d0+sign(0.5d0,v3(i,j,ke))
      v1(i,j,ke+1:ke+2) = v1(i,j,ke)
      v2(i,j,ke+1:ke+2) = v2(i,j,ke)
      v3(i,j,ke+1:ke+2) = v3(i,j,ke)*plug
      if(mag_on)then
      b1(i,j,ke+1:ke+2) = b1(i,j,ke)
      b2(i,j,ke+1:ke+2) = b2(i,j,ke)
      b3(i,j,ke+1:ke+2) = b3(i,j,ke)
      end if
    end do
    end do
  !$omp end do

  case(3) x3_outer_vector ! free ------------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ke+1:ke+2) = v1(i,j,ke)
      v2(i,j,ke+1:ke+2) = v2(i,j,ke)
      v3(i,j,ke+1:ke+2) = v3(i,j,ke)
      if(mag_on)then
      b1(i,j,ke+1:ke+2) = b1(i,j,ke)
      b2(i,j,ke+1:ke+2) = b2(i,j,ke)
      b3(i,j,ke+1:ke+2) = b3(i,j,ke)
      end if
    end do
    end do
  !$omp end do

  case(4) x3_outer_vector ! linear ----------------------------------------
  !$omp do private(i,j) collapse(2)
    do j = js, je
    do i = is, ie
      v1(i,j,ke+1) = v1(i,j,ke  ) + (v1(i,j,ke)-v1(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v1(i,j,ke+2) = v1(i,j,ke+1) + (v1(i,j,ke+1)-v1(i,j,ke))*dx3(ke+1)/dx3(ke)
      v2(i,j,ke+1) = v2(i,j,ke  ) + (v2(i,j,ke)-v2(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v2(i,j,ke+2) = v2(i,j,ke+1) + (v2(i,j,ke+1)-v2(i,j,ke))*dx3(ke+1)/dx3(ke)
      v3(i,j,ke+1) = v3(i,j,ke  ) + (v3(i,j,ke)-v3(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v3(i,j,ke+2) = v3(i,j,ke+1) + (v3(i,j,ke+1)-v3(i,j,ke))*dx3(ke+1)/dx3(ke)
      if(mag_on)then
      b1(i,j,ke+1) = b1(i,j,ke  ) + (b1(i,j,ke)-b1(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b1(i,j,ke+2) = b1(i,j,ke+1) + (b1(i,j,ke+1)-b1(i,j,ke))*dx3(ke+1)/dx3(ke)
      b2(i,j,ke+1) = b2(i,j,ke  ) + (b2(i,j,ke)-b2(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b2(i,j,ke+2) = b2(i,j,ke+1) + (b2(i,j,ke+1)-b2(i,j,ke))*dx3(ke+1)/dx3(ke)
      b3(i,j,ke+1) = b3(i,j,ke  ) + (b3(i,j,ke)-b3(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b3(i,j,ke+2) = b3(i,j,ke+1) + (b3(i,j,ke+1)-b3(i,j,ke))*dx3(ke+1)/dx3(ke)
      end if
    end do
    end do
  !$omp end do

  case(5) x3_outer_vector ! linear + outgoing -----------------------------
  !$omp do private(plug,i,j) collapse(2)
    do j = js, je
    do i = is, ie
      plug = 0.5d0+sign(0.5d0,v3(i,j,ke))
      v1(i,j,ke+1) = v1(i,j,ke  ) + (v1(i,j,ke)-v1(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v1(i,j,ke+2) = v1(i,j,ke+1) + (v1(i,j,ke+1)-v1(i,j,ke))*dx3(ke+1)/dx3(ke)
      v2(i,j,ke+1) = v2(i,j,ke  ) + (v2(i,j,ke)-v2(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v2(i,j,ke+2) = v2(i,j,ke+1) + (v2(i,j,ke+1)-v2(i,j,ke))*dx3(ke+1)/dx3(ke)
      v3(i,j,ke+1) = v3(i,j,ke  ) + (v3(i,j,ke)-v3(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      v3(i,j,ke+2) = v3(i,j,ke+1) + (v3(i,j,ke+1)-v3(i,j,ke))*dx3(ke+1)/dx3(ke)
      v1(i,j,ke+1:ke+2) = v1(i,j,ke+1:ke+2)
      v2(i,j,ke+1:ke+2) = v2(i,j,ke+1:ke+2)
      v3(i,j,ke+1:ke+2) = v3(i,j,ke+1:ke+2) * plug
      if(mag_on)then
      b1(i,j,ke+1) = b1(i,j,ke  ) + (b1(i,j,ke)-b1(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b1(i,j,ke+2) = b1(i,j,ke+1) + (b1(i,j,ke+1)-b1(i,j,ke))*dx3(ke+1)/dx3(ke)
      b2(i,j,ke+1) = b2(i,j,ke  ) + (b2(i,j,ke)-b2(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b2(i,j,ke+2) = b2(i,j,ke+1) + (b2(i,j,ke+1)-b2(i,j,ke))*dx3(ke+1)/dx3(ke)
      b3(i,j,ke+1) = b3(i,j,ke  ) + (b3(i,j,ke)-b3(i,j,ke-1))*dx3(ke)/dx3(ke-1)
      b3(i,j,ke+2) = b3(i,j,ke+1) + (b3(i,j,ke+1)-b3(i,j,ke))*dx3(ke+1)/dx3(ke)
      end if
    end do
    end do
  !$omp end do

  case(9) x3_outer_vector ! Dirichlet -------------------------------------
  !$omp do private(i,j,k) collapse(3)
    do k = ke+1, ke+2
    do j = js, je
      do i = is, ie
      v1(i,j,k) = v10(i,j,k)
      v2(i,j,k) = v20(i,j,k)
      v3(i,j,k) = v30(i,j,k)
      if(mag_on)then
        b1(i,j,k) = b10(i,j,k)
        b2(i,j,k) = b20(i,j,k)
        b3(i,j,k) = b30(i,j,k)
      end if
      end do
    end do
    end do
  !$omp end do

  case(10) x3_outer_vector ! Flux -----------------------------------------

  case default x3_outer_vector ! Error ------------------------------------
    print *, "Error from x3 velocity outer boundary condition" ; stop
  end select x3_outer_vector

end if

  ! set e and ptot =========================================================
  !$omp do private(i,j,k) collapse(3)
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
  !$omp end do
  ! ========================================================================
end if

!$omp end parallel

call stop_clock(wtbnd)

return

end subroutine boundarycondition

end module boundary_mod
