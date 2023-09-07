module conserve_mod
 implicit none
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE CONSERVE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert physical values to conserved values

subroutine conserve

  use settings,only:mag_on
  use grid,only:is,ie,js,je,ks,ke
  use physval,only:u,d,e,v1,v2,v3,b1,b2,b3,phi,&
                   icnt,imo1,imo2,imo3,iene,img1,img2,img3,idcl

  integer::i,j,k

!---------------------------------------------------------------------------

  ! set U
!$omp parallel do private(i,j,k)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     u(i,j,k,icnt) = d(i,j,k)
     u(i,j,k,imo1) = d(i,j,k) * v1(i,j,k)
     u(i,j,k,imo2) = d(i,j,k) * v2(i,j,k)
     u(i,j,k,imo3) = d(i,j,k) * v3(i,j,k)
     u(i,j,k,iene) = e(i,j,k)
     if(.not.mag_on)cycle
     u(i,j,k,img1) = b1(i,j,k)
     u(i,j,k,img2) = b2(i,j,k)
     u(i,j,k,img3) = b3(i,j,k)
     u(i,j,k,idcl) = phi(i,j,k)
    end do
   end do
  end do
!$omp end parallel do

return
end subroutine conserve

end module conserve_mod
