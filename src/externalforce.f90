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
 use utils,only:polcar,cylcar,get_vpol,get_vcyl

 integer:: i,j,k
 real(8),dimension(1:3):: atot,ftot,xcar

!-----------------------------------------------------------------------------

 select case(crdnt)
 case(0) ! cartesian coordinates
!$omp parallel do private(i,j,k,atot,ftot) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     atot = frame_acc
     ftot = d(i,j,k) * atot
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atot)*dt
    end do
   end do
  end do
!$omp end parallel do

 case(1) ! cylindrical coordinates
!$omp parallel do private(i,j,k,xcar,atot,ftot) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     xcar = cylcar([x1(i),x2(j),x3(k)])
     call get_vcyl(xcar,frame_acc,atot(1),atot(2),atot(3))
     ftot = d(i,j,k) * atot
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atot)*dt
    end do
   end do
  end do
!$omp end parallel do

 case(2) ! spherical coordinates
!$omp parallel do private(i,j,k,xcar,atot,ftot) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     xcar = polcar([x1(i),x2(j),x3(k)])
     call get_vpol(xcar,x3(k),frame_acc,atot(1),atot(2),atot(3))
     ftot = d(i,j,k) * atot
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atot)*dt
    end do
   end do
  end do
!$omp end parallel do
 end select


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
