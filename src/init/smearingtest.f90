module smearingtest_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SMEARINGTEST
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up a specific test for smearing

subroutine smearingtest

 use grid
 use physval
 use gravmod

 integer::i,j,k

!-----------------------------------------------------------------------------

 grvphi = -3d-5
 do k = ks, ke
  do j = js, je
   do i = is, ie
     ! Arbitrary arrangement that varies with i, j, k
     d(i,j,k) = 1d1*(real(ie_global-i+1)) + real(j) + real(k)
     p(i,j,k) = 1d1*(real(ie_global-i+1)) + real(j) + real(k)
   end do
  end do
 end do

return
end subroutine smearingtest

end module smearingtest_mod
