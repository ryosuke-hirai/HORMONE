module eigen_mod

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

end module eigen_mod
