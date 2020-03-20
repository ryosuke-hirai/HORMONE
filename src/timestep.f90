!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE TIMESTEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set dt

 subroutine timestep

  use settings,only:courant,outstyle
  use grid
  use physval
  use ninewave
  use eigen_mod

  implicit none

  real*8,allocatable:: dtdist(:,:,:,:)
  real*8 cfmax

!-------------------------------------------------------------------------

  call eigen

  dt = 1.d50
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
  dtdist = 1d50
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

  dt = minval( dtdist(is:ie,js:je,ks:ke,1:3) )

  cfmax = maxval(abs(cf))

  dt = courant * dt
!print *,dt,minloc(dtdist(is:ie,js:je,ks:ke,1:3));stop
  if(outstyle==1) dt = min(dt,t_out-time)

  ch = cfmax

 return
 end subroutine timestep
