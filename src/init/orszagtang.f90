module orszagtang_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                    SUBROUTINE ORSZAG-TANG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for the Orszag-Tang problem

subroutine orszagtang

  use grid
  use physval
  use constants,only:pi

  integer:: i,j,k,whichaxis
!--------------------------------------------------------------------

  if(ke==ks)whichaxis=1
  if(je==js)whichaxis=2
  if(ie==is)whichaxis=3
  select case (whichaxis)
  case(1)
! for Orszag-Tang problem : xy surface
   k = ks

   do j = js,je
    do i = is,ie
     d(i,j,k) = gamma**2/(4.d0*pi)
     p(i,j,k) = gamma/(4.d0*pi)
     v1(i,j,k)= -sin(2.d0*pi*x2(j))
     v2(i,j,k)=  sin(2.d0*pi*x1(i))
     v3(i,j,k)= 0d0
     b1(i,j,k)= -sin(2.d0*pi*x2(j))/sqrt(4.d0*pi)
     b2(i,j,k)=  sin(4.d0*pi*x1(i))/sqrt(4.d0*pi)
     b3(i,j,k)= 0d0
    end do
   end do

  case(2)
! for Orszag-Tang problem : xz surface
   j = js

   do k = ks,ke
    do i = is,ie
     d(i,j,k) = gamma**2/(4.d0*pi)
     p(i,j,k) = gamma/(4.d0*pi)
     v1(i,j,k)= -sin(2.d0*pi*x3(k))
     v2(i,j,k)= 0d0
     v3(i,j,k)=  sin(2.d0*pi*x1(i))
     b1(i,j,k)= -sin(2.d0*pi*x3(k))/sqrt(4.d0*pi)
     b2(i,j,k)= 0d0
     b3(i,j,k)=  sin(4.d0*pi*x1(i))/sqrt(4.d0*pi)
    end do
   end do

  case(3)
! for Orszag-Tang problem : yz surface
   i = is

   do k = ks,ke
    do j = js,je
     d(i,j,k) = gamma**2/(4.d0*pi)
     p(i,j,k) = gamma/(4.d0*pi)
     v1(i,j,k)= 0d0
     v2(i,j,k)= -sin(2.d0*pi*x3(k))
     v3(i,j,k)=  sin(2.d0*pi*x2(j))
     b1(i,j,k)= 0d0
     b2(i,j,k)= -sin(2.d0*pi*x3(k))/sqrt(4.d0*pi)
     b3(i,j,k)=  sin(4.d0*pi*x2(j))/sqrt(4.d0*pi)
    end do
   end do

  case default
   stop 'Something is wrong with the grid. Make it 2D.'
  end select

return
end subroutine orszagtang

end module orszagtang_mod
