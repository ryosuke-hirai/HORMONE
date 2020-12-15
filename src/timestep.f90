module timestep_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EIGEN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte eigenvalues

 subroutine eigen

  use grid
  use settings,only:eostype
  use physval
  use pressure_mod

  implicit none

  real*8 asq, bb1, bb2, bb3, bbsq
!--------------------------------------------------------------------------
!$omp parallel do private(i,j,k,asq,bbsq,bb1,bb2,bb3)
  do k = ks,ke
   do j = js,je
    do i = is,ie

     select case (eostype)
     case(0:1) ! without recombination
      call eos_p_cf(d(i,j,k), v1(i,j,k), v2(i,j,k), v3(i,j,k), &
                              b1(i,j,k), b2(i,j,k), b3(i,j,k), &
                    e(i,j,k), T (i,j,k), imu(i,j,k),p (i,j,k), cf(i,j,k) )
     case(2) ! with recombination
      call eos_p_cf(d(i,j,k), v1(i,j,k), v2(i,j,k), v3(i,j,k), &
                              b1(i,j,k), b2(i,j,k), b3(i,j,k), &
                    e(i,j,k), T (i,j,k), imu(i,j,k),p (i,j,k), cf(i,j,k), &
                    spc(1,i,j,k), spc(2,i,j,k) )
     end select

     asq = cf(i,j,k)*cf(i,j,k)

     bb1  = b1(i,j,k)**2 / d(i,j,k)
     bb2  = b2(i,j,k)**2 / d(i,j,k)
     bb3  = b3(i,j,k)**2 / d(i,j,k)
     bbsq = bb1 + bb2 + bb3

     cf(i,j,k) = asq + bbsq + &
          sqrt( (asq+bbsq)**2 - 4d0*asq*bb1 )
     cf(i,j,k) = sqrt( 0.5d0*cf(i,j,k) )
    end do
   end do
  end do
!$omp end parallel do
 
  return
 end subroutine eigen


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE TIMESTEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set dt.

 subroutine timestep

  use settings,only:courant,outstyle
  use grid
  use physval
  use constants,only:huge

  implicit none

  real*8,allocatable:: dtdist(:,:,:,:)
  real*8 cfmax
  
!-------------------------------------------------------------------------

  call eigen

  dt = huge
!!$  cfmax = 0.d0
!!$  do k = ks,ke
!!$   do j = js,je
!!$    do i = is,ie
!!$     dt = min( dt , &
!!$                     dxi1(i) / ( abs(v1(i,j,k))+cf(i,j,k) ), &
!!$          g22(i)   * dxi2(j) / ( abs(v2(i,j,k))+cf(i,j,k) ), &
!!$          g33(i,j) * dxi3(k) / ( abs(v3(i,j,k))+cf(i,j,k) )  &
!!$          )
!!$!     cfmax = max(abs(cf(i,j,k)),cfmax)
!!$    end do
!!$   end do
!!$  end do

!!$  do k = ks,ke
!!$   do j = js,je
!!$    do i = is,ie
!!$      dt = min( dt , &
!!$                     dxi1(i) / ( abs(v1(i,j,k))+v3(i,j,k)+cf(i,j,k) ) )
!!$     if(je>1)then
!!$      if(i>sphrn)then
!!$       dt = min( dt , &
!!$           g22(i)   * dxi2(j) / ( abs(v2(i,j,k)+v3(i,j,k))+cf(i,j,k) ) )
!!$      end if
!!$     end if
!!$     if(ke>1)then
!!$      if(i>sphrn)then
!!$       dt = min( dt , &
!!$           g33(i,j) * dxi3(k) / ( abs(v3(i,j,k))+v2(i,j,k)+cf(i,j,k) ) )
!!$      end if
!!$     end if
!!$!     cfmax = max(abs(cf(i,j,k)),cfmax)
!!$    end do
!!$   end do
!!$  end do

  allocate( dtdist(is:ie,js:je,ks:ke,1:3) )
  dtdist = huge
  do k = ks, ke
   do j = js, je
    do i = is, ie
     dtdist(i,j,k,1) = dxi1(i) / ( abs(v1(i,j,k))+abs(v3(i,j,k))+cf(i,j,k) )
     if(je>1)then
      dtdist(i,j,k,2) = g22(i)*dxi2(j) &
                      / ( abs(v2(i,j,k))+abs(v3(i,j,k))+cf(i,j,k) )
     end if
     if(ke>1)then
      dtdist(i,j,k,3) = g33(i,j)*dxi3(k) &
                      / ( abs(v3(i,j,k)+cf(i,j,k) ) )
     end if
    end do
   end do
  end do

! Temporary (for workaround mesh)
  if(sphrn>0.and.crdnt==2)then
   do i = is, is+2
    dtdist(i,js:je,ks:ke,2) = sum( dtdist(i,js:je,ks:ke,2) )
   end do
   do i = is+3, sphrn
    do j = js, je,40
     dtdist(i,j:j+39,ks:ke,2) = sum( dtdist(i,j:j+39,ks:ke,2) )
    end do
   end do
   do i = sphrn+1, sphrn+trnsn1
    do j = js, je,8
     dtdist(i,j:j+7,ks:ke,2) = sum( dtdist(i,j:j+7,ks:ke,2) )    
    end do
   end do
   do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
    do j = js, je,4
     dtdist(i,j:j+3,ks:ke,2) = sum( dtdist(i,j:j+3,ks:ke,2) )
    end do
   end do
   do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
    do j = js, je,2
     dtdist(i,j:j+1,ks:ke,2) = sum( dtdist(i,j:j+1,ks:ke,2) )
    end do
   end do
  end if

  dt = minval( dtdist(is:ie,js:je,ks:ke,1:3) )

  cfmax = maxval(abs(cf))

  dt = courant * dt
!print *,dt,minloc(dtdist(is:ie,js:je,ks:ke,1:3));stop
  if(outstyle==1) dt = min(dt,t_out-time)

  ch = cfmax

 return
 end subroutine timestep

end module timestep_mod
