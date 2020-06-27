!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE PRIMITIVE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert conserved values to primitive values

subroutine primitive

  use grid
  use physval
  use pressure_mod,only:pressure

  implicit none

!--------------------------------------------------------------------
  
!$omp parallel do private(i,j,k)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     d(i,j,k)  = u(i,j,k,1)
     v1(i,j,k) = u(i,j,k,2) / d(i,j,k)
     v2(i,j,k) = u(i,j,k,3) / d(i,j,k)
     v3(i,j,k) = u(i,j,k,4) / d(i,j,k)
     b1(i,j,k) = u(i,j,k,5)
     b2(i,j,k) = u(i,j,k,6)
     b3(i,j,k) = u(i,j,k,7)
     e(i,j,k)  = u(i,j,k,8)
     ! for 9 wave method
     phi(i,j,k)= u(i,j,k,9)
    end do
   end do
  end do
!$omp end parallel do

  call pressure

return
end subroutine primitive
