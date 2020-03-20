!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE KHTESTS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for Kelvini-Helmholtz tests.

subroutine KHtest

  use grid
  use physval

  implicit none
  
  real*8 rdm
  real*8,parameter:: ptb = 1.d-1, frac = 5.d0
  integer,parameter:: l = 5 ! meshes

!--------------------------------------------------------------------

! for Kelvin-Helmholtz test
  do k = ks-1,ke+1
   do j = js-1,je+1
    do i = is-1,ie+1
     d(i,j,k)  = 0.5d0 * ( (1.d0+frac)-(1.d0-frac)*dtanh(x3(k)/(dble(l)*dx3(k))) )
     v1(i,j,k) = 1.d0 * dtanh( x3(k)/(dble(l)*dx3(k)) )
     b1(i,j,k) = 0.d0
     b2(i,j,k) = 0.d0
     b3(i,j,k) = 0.d0

     call random_number(rdm)
     v3(i,j,k) =  ptb * (2.d0*rdm-1.d0)
     v2(i,j,k) = 0.d0
     p(i,j,k)  = 1.d1
    end do
   end do
  end do

  do j = je/2,je+1
   do i = is-1,ie+1
    d(i,j,k) = 5.d0
    p(i,j,k) = 1.d1
    v1(i,j,k) = -1.d0
   end do
  end do

return
end subroutine KHtest
