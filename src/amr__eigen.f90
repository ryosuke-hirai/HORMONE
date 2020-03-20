module amr_eigen_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE AMR_EIGEN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate eigenvalues for AMR mode

subroutine amr_eigen(lf)

 use grid,only:is,js,ks
 use physval,only:gamma
 use amr_templates
 use amr_module,only:ib,jb,kb

 implicit none

 integer i,j,k
 real*8 asq, bb1, bb2, bb3, bbsq
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

 do k = ks-2, kb+2
  do j = js-2, jb+2
   do i = is-2, ib+2
    asq  = gamma*lf%p(i,j,k) / lf%d(i,j,k)

    bb1  = lf%b1(i,j,k)*lf%b1(i,j,k) / lf%d(i,j,k)
    bb2  = lf%b2(i,j,k)*lf%b2(i,j,k) / lf%d(i,j,k)
    bb3  = lf%b3(i,j,k)*lf%b3(i,j,k) / lf%d(i,j,k)
    bbsq = bb1 + bb2 + bb3

    lf%cf(i,j,k) = asq + bbsq + sqrt( (asq+bbsq)**2.d0 - 4.d0*asq*bb1 )
    lf%cf(i,j,k) = sqrt( 5.d-1 * lf%cf(i,j,k) )
!if(j==1.and.k==1)print *,i,lf%cf(i,j,k),lf%p(i,j,k)
   end do
  end do
 end do

return
end subroutine amr_eigen

end module amr_eigen_mod
