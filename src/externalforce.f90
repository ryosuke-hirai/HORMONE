module externalforce_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE EXTERNALFORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Add a non-gravitational external force to the source term

subroutine externalforce

 use grid
 use physval

 implicit none

!-----------------------------------------------------------------------------

 

 return
end subroutine externalforce

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE EXTERNALFIELD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To add external gravitational field to the gravitational potential.

subroutine externalfield

 use grid
 use gravmod,only:extgrv,totphi

 integer:: i,j,k

!-----------------------------------------------------------------------------

!$omp parallel do private (i,j,k) collapse(3)
 do k = ks-1, ke+1
  do j = js-1, je+1
   do i = is-1, ie+1
    totphi(i,j,k) = totphi(i,j,k) + extgrv(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

 return
end subroutine externalfield

end module externalforce_mod
