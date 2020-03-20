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
  use physval
  use pressure_mod

  implicit none

  real*8 asq, bb1, bb2, bb3, bbsq
!--------------------------------------------------------------------------
!$omp parallel do private(i,j,k,asq,bbsq,bb1,bb2,bb3)
  do k = ks,ke
   do j = js,je
    do i = is,ie

     call eos_p_cf(d(i,j,k), v1(i,j,k), v2(i,j,k), v3(i,j,k), &
                             b1(i,j,k), b2(i,j,k), b3(i,j,k), &
                   e(i,j,k), T (i,j,k), imu(i,j,k),p (i,j,k), cf(i,j,k) )
     asq = cf(i,j,k)*cf(i,j,k)

     bb1  = b1(i,j,k)*b1(i,j,k) / d(i,j,k)
     bb2  = b2(i,j,k)*b2(i,j,k) / d(i,j,k)
     bb3  = b3(i,j,k)*b3(i,j,k) / d(i,j,k)
     bbsq = bb1 + bb2 + bb3

     cf(i,j,k) = asq + bbsq + &
          sqrt( (asq+bbsq)**2.d0 - 4.d0*asq*bb1 )
     cf(i,j,k) = sqrt( 5.d-1 * cf(i,j,k) )
    end do
   end do
  end do
!$omp end parallel do
 
return
end subroutine eigen


real*8 function soundvel(p,d,b1,b2,b3)

  use physval,only:gamma

  implicit none

  real*8 p,d,b1,b2,b3,c,asq,bb1,bb2,bb3,bbsq

  asq = gamma * p/d

  bb1 = b1*b1/d ; bb2 = b2*b2/d ; bb3 = b3*b3/d
  bbsq = bb1 + bb2 + bb3

  c = asq + bbsq
  soundvel = sqrt( 0.5d0*( c + sqrt( c*c - 4d0*asq*bb1 ) ) )

end function soundvel

end module eigen_mod




module aeigen_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EIGEN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte eigenvalues

subroutine eigen

  use grid
  use physval

  implicit none

  real*8 asq, bb1, bb2, bb3, bbsq
!--------------------------------------------------------------------------
!$omp parallel do private(i,j,k,asq,bbsq,bb1,bb2,bb3)
  do k = ks-2,ke+2
   do j = js-2,je+2
    do i = is-2,ie+2
     asq = gamma*p(i,j,k) / d(i,j,k)

     bb1  = b1(i,j,k)*b1(i,j,k) / d(i,j,k)
     bb2  = b2(i,j,k)*b2(i,j,k) / d(i,j,k)
     bb3  = b3(i,j,k)*b3(i,j,k) / d(i,j,k)
     bbsq = bb1 + bb2 + bb3

     cf(i,j,k) = asq + bbsq + &
          sqrt( (asq+bbsq)**2.d0 - 4.d0*asq*bb1 )
     cf(i,j,k) = sqrt( 5.d-1 * cf(i,j,k) )
    end do
   end do
  end do
!$omp end parallel do
 
return
end subroutine eigen


real*8 function soundvel(p,d,b1,b2,b3)

  use physval,only:gamma

  implicit none

  real*8 p,d,b1,b2,b3,c,asq,bb1,bb2,bb3,bbsq

  asq = gamma * p/d

  bb1 = b1*b1/d ; bb2 = b2*b2/d ; bb3 = b3*b3/d
  bbsq = bb1 + bb2 + bb3

  c = asq + bbsq
  soundvel = sqrt( 0.5d0*( c + sqrt( c*c - 4d0*asq*bb1 ) ) )

end function soundvel

end module aeigen_mod
