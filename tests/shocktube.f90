!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SHOCKTUBE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for shock tube problems.

subroutine shocktube

  use grid
  use physval

  implicit none

!--------------------------------------------------------------------

! for shock tube (See Brio & Wu): i direction

!!$  do i = is-2,ie+2
!!$   do j = js-2,je+2
!!$    do k = ks-2,ke+2
!!$     d(i,j,k) = 1.d0
!!$     p(i,j,k) = 1.d0
!!$     v1(i,j,k) = 0.d0
!!$     v2(i,j,k) = 0.d0
!!$     v3(i,j,k) = 0.d0
!!$     b1(i,j,k) = 0.75d0
!!$     b2(i,j,k) = 1.d0
!!$     b3(i,j,k) = 0.d0
!!$    end do
!!$   end do
!!$  end do
!!$
!!$ do i = in/2+1,ie+2
!!$  do j = js-2,je+2
!!$   do k = ks-2,ke+2
!!$    d(i,j,k) = 0.125d0
!!$    p(i,j,k) = 1.d-1
!!$    b2(i,j,k) = -1.d0
!!$   end do
!!$  end do
!!$ end do

! for shock tube (See Brio & Wu): j direction
!!$  do k = ks-1,ke+1
!!$   do j = js-1,je+1
!!$    do i = is-1,ie+1
!!$     d(i,j,k) = 1.d0
!!$     p(i,j,k) = 1.d0
!!$     v1(i,j,k) = 0.d0
!!$     v2(i,j,k) = 0.d0
!!$     v3(i,j,k) = 0.d0
!!$     b1(i,j,k) = 0.d0
!!$     b2(i,j,k) = 0.75d0
!!$     b3(i,j,k) = 1.d0
!!$    end do
!!$   end do
!!$  end do
!!$
!!$  do k = ks-1,ke+1
!!$   do i = is-1,ie+1
!!$    do j = jn/2,je+1
!!$     d(i,j,k) = 0.125d0
!!$     p(i,j,k) = 1.d-1
!!$     b3(i,j,k) = -1.d0
!!$    end do
!!$   end do
!!$  end do

! for shock tube (See Brio & Wu): k direction
!!$  do k = ks-1,ke+1
!!$   do j = js-1,je+1
!!$    do i = is-1,ie+1
!!$     d(i,j,k) = 1.d0
!!$     p(i,j,k) = 1.d0
!!$     v1(i,j,k) = 0.d0
!!$     v2(i,j,k) = 0.d0
!!$     v3(i,j,k) = 0.d0
!!$     b1(i,j,k) = 1.d0
!!$     b2(i,j,k) = 0.d0
!!$     b3(i,j,k) = 0.75d0
!!$    end do
!!$   end do
!!$  end do
!!$
!!$  do k = kn/2,ke+1
!!$   do j = js-1,je+1
!!$    do i = is-1,ie+1
!!$     d(i,j,k) = 0.125d0
!!$     b1(i,j,k) = -1.d0
!!$     p(i,j,k) = 1.d-1
!!$    end do
!!$   end do
!!$  end do

! for shock tube x direction
!!$
!!$k = ks ; j = js
!!$  do i = is,in/2
!!$   d(i,j,k)  = 1.08d0
!!$   p(i,j,k)  = 0.95d0
!!$   v1(i,j,k) = 1.2d0
!!$   v2(i,j,k) = 0.01d0
!!$   v3(i,j,k) = 0.5d0
!!$   b1(i,j,k) = 2.d0/dsqrt(4.d0*pi)
!!$   b2(i,j,k) = 3.6d0/dsqrt(4.d0*pi)
!!$   b3(i,j,k) = 2.d0/dsqrt(4.d0*pi)
!!$  end do
!!$  do i = in/2+1,ie
!!$   d(i,j,k)  = 1.d0
!!$   p(i,j,k)  = 1.d0
!!$   v1(i,j,k) = 0.d0
!!$   v2(i,j,k) = 0.d0
!!$   v3(i,j,k) = 0.d0
!!$   b1(i,j,k) = 2.d0/dsqrt(4.d0*pi)
!!$   b2(i,j,k) = 4.d0/dsqrt(4.d0*pi)
!!$   b3(i,j,k) = 2.d0/dsqrt(4.d0*pi)
!!$  end do

return
end subroutine shocktube
