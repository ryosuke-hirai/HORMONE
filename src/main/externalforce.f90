module externalforce_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EXTERNALFORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Add a non-gravitational external force to the source term

subroutine externalforce

 use settings,only:frame
 use external_settings,only:include_spinup
 use grid
 use physval
 use utils,only:polcar,cylcar,get_vpol,get_vcyl

 integer:: i,j,k
 real(8),dimension(1:3):: ftot,atotijk
 real(8),allocatable,dimension(:,:,:,:):: atot

!-----------------------------------------------------------------------------

 allocate(atot(1:3,is:ie,js:je,ks:ke))

! Initialize atot
!$omp parallel do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    atot(1:3,i,j,k) = 0d0
   end do
  end do
 end do
!$omp end parallel do

! Get frame force
 if(frame>0)       call frame_force(atot)

! Spin up star
 if(include_spinup)call spinup(atot)


! Add external acceleration to the source term in hydrodynamic equations
 select case(crdnt)
 case(0) ! cartesian coordinates
!$omp parallel do private(i,j,k,ftot,atotijk) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     atotijk = atot(1:3,i,j,k)
     ftot = d(i,j,k) * atotijk
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atotijk)*dt
    end do
   end do
  end do
!$omp end parallel do

 case(1) ! cylindrical coordinates
!$omp parallel do private(i,j,k,ftot,atotijk) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     atotijk =  atot(1:3,i,j,k)
     ftot = d(i,j,k) * atotijk
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atotijk)*dt
    end do
   end do
  end do
!$omp end parallel do

 case(2) ! spherical coordinates
!$omp parallel do private(i,j,k,ftot,atotijk) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     atotijk = atot(1:3,i,j,k)
     ftot = d(i,j,k) * atotijk
     src(i,j,k,imo1) = src(i,j,k,imo1) + ftot(1)
     src(i,j,k,imo2) = src(i,j,k,imo2) + ftot(2)
     src(i,j,k,imo3) = src(i,j,k,imo3) + ftot(3)
     src(i,j,k,iene) = src(i,j,k,iene) &
                     + (ftot(1)*v1(i,j,k)+ftot(2)*v2(i,j,k)+ftot(3)*v3(i,j,k)) &
                     + 0.5d0*dot_product(ftot,atotijk)*dt
    end do
   end do
  end do
!$omp end parallel do
 end select

 deallocate(atot)

end subroutine externalforce

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE FRAME_FORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Add frame force to the total acceleration

subroutine frame_force(atot)

 use grid
 use utils,only:polcar,cylcar,get_vpol,get_vcyl

 real(8),allocatable,dimension(:,:,:,:),intent(inout):: atot
 integer:: i,j,k
 real(8),dimension(1:3):: aframe, xcar

!-----------------------------------------------------------------------------

 select case(crdnt)
 case(0) ! cartesian coordinates
!$omp parallel do private(i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     atot(1:3,i,j,k) = atot(1:3,i,j,k) + frame_acc(1:3)
    end do
   end do
  end do
!$omp end parallel do

 case(1) ! cylindrical coordinates
!$omp parallel do private(i,j,k,xcar,aframe) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     xcar = cylcar([x1(i),x2(j),x3(k)])
     call get_vcyl(xcar,frame_acc,aframe(1),aframe(2),aframe(3))
     atot(1:3,i,j,k) = atot(1:3,i,j,k) + aframe(1:3)
    end do
   end do
  end do
!$omp end parallel do

 case(2) ! spherical coordinates
!$omp parallel do private(i,j,k,aframe) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     call get_vpol(car_x(:,i,j,k),x3(k),frame_acc,aframe(1),aframe(2),aframe(3))
     atot(1:3,i,j,k) = atot(1:3,i,j,k) + aframe(1:3)
    end do
   end do
  end do
!$omp end parallel do
 end select

 return
end subroutine frame_force

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE SPINUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To spin up star in rigid rotation but capped at Keplerian rotation

subroutine spinup(atot)

 use external_settings,only:omegadot,j_max
 use grid,only:crdnt,is,ie,js,je,ks,ke,x1,x2,sinc,spinc_r,spinc_t
 use physval,only:v1,v2,v3
 use gravmod,only:totphi

 real(8),allocatable,dimension(:,:,:,:),intent(inout):: atot
 integer:: i,j,k
 real(8):: v_kep,j_local,v_phi

!-----------------------------------------------------------------------------

 select case(crdnt)
 case(0) ! cartesian coordinates
!$omp parallel do private(i,j,k,j_local,v_kep,v_phi) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     j_local = x1(i)*v2(i,j,k)-x2(j)*v1(i,j,k)
     v_phi = (-v1(i,j,k)*x2(j)+v2(i,j,k)*x1(i))/sqrt(x1(i)**2+x2(j)**2)
     v_kep = sqrt(-totphi(i,j,k))
     if(j_local<=j_max.and.v_phi<=v_kep)then
      atot(1,i,j,k) = atot(1,i,j,k) - omegadot*x2(j)
      atot(2,i,j,k) = atot(2,i,j,k) + omegadot*x1(i)
     end if
    end do
   end do
  end do
!$omp end parallel do

 case(1) ! cylindrical coordinates
!$omp parallel do private(i,j,k,j_local,v_kep) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     j_local = x1(i)*v2(i,j,k)
     v_kep = sqrt(-totphi(i,j,k))
     if(j_local<=j_max.and.v2(i,j,k)<=v_kep)then
      atot(2,i,j,k) = atot(2,i,j,k) + omegadot*spinc_r(i)
     end if
    end do
   end do
  end do
!$omp end parallel do

 case(2) ! spherical coordinates
!$omp parallel do private(i,j,k,j_local,v_kep) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     j_local = x1(i)*sinc(j)*v3(i,j,k)
     v_kep = sqrt(-totphi(i,j,k))
     if(j_local<=j_max.and.v3(i,j,k)<=v_kep)then
      atot(3,i,j,k) = atot(3,i,j,k) + omegadot*spinc_r(i)*spinc_t(j)
     end if
    end do
   end do
  end do
!$omp end parallel do
 end select

end subroutine spinup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EXTERNALFIELD
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
