!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SHOCKFIND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To identify the location of shocks

subroutine shockfind

  use grid
  use physval
  
 implicit none

 real*8 :: divergence, gradT(1:dim), gradd(1:dim), Mjump

!-----------------------------------------------------------------------------

 shock = 0

 if(crdnt==2.and.ke==1)then
!$omp parallel do private(i,j,k,divergence,gradT,gradd,Mjump)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     divergence = (dx1(i)**2d0*v1(i+1,j,k)-dx1(i+1)**2d0*v1(i-1,j,k)+(dx1(i+1)**2d0-dx1(i)**2d0)*v1(i,j,k))*idx1(i)*idx1(i+1)/sum(dx1(i:i+1)) + 2d0*v1(i,j,k)/x1(i) &
      + ((dx2(j)**2d0*v2(i,j+1,k)-dx2(j+1)**2d0*v2(i,j-1,k)+(dx2(j+1)**2d0-dx2(j)**2d0)*v2(i,j,k))*idx2(j)*idx2(j+1)/sum(dx2(j:j+1))+v2(i,j,k)/tan(x2(j)))/x1(i)
     if(divergence<0d0)then
      gradT(1) = T(i+1,j,k) - T(i-1,j,k)
      gradT(2) = T(i,j+1,k) - T(i,j-1,k)
      gradd(1) = d(i+1,j,k) - d(i-1,j,k)
      gradd(2) = d(i,j+1,k) - d(i,j-1,k)
      if(gradT(1)*gradd(1)>0d0.and.gradT(2)*gradd(2)>0d0)then
       Mjump = max( &
        abs(Mach(i,j,k)-Mach(i-1,j,k)), abs(Mach(i,j,k)-Mach(i+1,j,k)),&
        abs(Mach(i,j,k)-Mach(i,j-1,k)), abs(Mach(i,j,k)-Mach(i,j+1,k)) ) 
       if(Mjump>1.3d0)shock(i,j,k) = 1
      end if
     end if
    end do
   end do
  end do
!$omp end parallel do
 elseif(crdnt==1.and.je==1)then
!$omp parallel do private(i,j,k,divergence,gradT,gradd,Mjump)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     divergence = (dx1(i)**2d0*v1(i+1,j,k)-dx1(i+1)**2d0*v1(i-1,j,k)+(dx1(i+1)**2d0-dx1(i)**2d0)*v1(i,j,k))*idx1(i)*idx1(i+1)/sum(dx1(i:i+1)) + v1(i,j,k)/x1(i) &
      + (dx3(k)**2d0*v3(i,j,k+1)-dx3(k+1)**2d0*v3(i,j,k-1)+(dx3(k+1)**2d0-dx3(k)**2d0)*v3(i,j,k))*idx3(k)*idx3(k+1)/sum(dx3(k:k+1))
     if(divergence<0d0)then
      gradT(1) = T(i+1,j,k) - T(i-1,j,k)
      gradT(2) = T(i,j,k+1) - T(i,j,k-1)
      gradd(1) = d(i+1,j,k) - d(i-1,j,k)
      gradd(2) = d(i,j,k+1) - d(i,j,k-1)
      if(gradT(1)*gradd(1)>0d0.and.gradT(2)*gradd(2)>0d0)then
!!$     Mjump = max( &
!!$      abs(Mach(i,j,k)-Mach(i-1,j,k)), abs(Mach(i,j,k)-Mach(i+1,j,k)),&
!!$      abs(Mach(i,j,k)-Mach(i,j,k-1)), abs(Mach(i,j,k)-Mach(i,j,k+1)) ) 
       Mjump = max( Mach(i+1,j,k),Mach(i-1,j,k), &
        Mach(i,j,k-1),Mach(i,j,k+1), Mach(i,j,k) )
       Mjump = ( sign(0.5d0,Mjump-0.1d0)+0.5d0 )*Mjump &
             / max(0.1d0,min( Mach(i+1,j,k),Mach(i-1,j,k), &
                    Mach(i,j,k-1),Mach(i,j,k+1),Mach(i,j,k)))
       if(Mjump>2.5d0)shock(i,j,k) = 1
      end if
     end if
    end do
   end do
  end do
!$omp end parallel do
 end if
 return
 
contains
 
 real*8 function Mach(ii,jj,kk)
  implicit none
  integer,intent(in):: ii,jj,kk
  Mach = sqrt( v1(ii,jj,kk)*v1(ii,jj,kk)   &
              +v2(ii,jj,kk)*v2(ii,jj,kk)   &
              +v3(ii,jj,kk)*v3(ii,jj,kk) ) &
       / cf(ii,jj,kk)
 end function Mach
 
end subroutine shockfind