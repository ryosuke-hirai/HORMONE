!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE SOURCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate source terms

subroutine source

 use grid
 use settings,only:eq_sym
 use physval
 use constants
 use gravmod
 use minmod_mod
 use merger_mod,only:domega_dt,de_dt,sep,spin_coeffr,spin_coefft
 
 implicit none

 real*8 uu(1:3), ptd, phil, phir, ul, ur
 real*8,dimension(1:2):: dx

!----------------------------------------------------------------------------

! To calculate gravitational forces ********************************************
  if(gravswitch==0)then ! gravity off
   grv1 = 0.d0 ; grv2 = 0.d0 ; grv3 = 0.d0

  elseif(gravswitch==1.and.crdnt==2.and.dim<=2)then ! point-source
   k = ks
   do i = is, ie
    mc(i) = mc(i-1) + sum( d(i,js:je,k) * dvol(i,js:je,k) )
   end do

!$omp parallel do private(i,ul)
   do i = is, ie
    ul = 1d0/(xi1(i)*(xi1(i)+xi1(i-1))+xi1(i-1)*xi1(i-1))
!    grv1(i,js,k) = - (mc(i-1)*dxi1(i)+0.25d0*(mc(i)-mc(i-1))*&
!                   (xi1(i)-3d0*xi1(i-1)*xi1(i-1)*xi1(i-1)*ul))&
!                  *3d0*ul * idxi1(i)

!    grv1(i,js,k) = -(mc(i-1)*dxi1(i) + pi/3d0*d(i,js,k)*(3d0*xi1(i-1)**4d0+xi1(i)**4d0-4d0*xi1(i-1)**3d0*xi1(i)))*ul*idxi1(i)*3d0*d(i,js,k)*G
!    grv1(i,js,k) = -(mc(i-1)+d(i,js,k)*4d0/3d0*pi*(x1(i)**3d0-xi1(i-1)**3d0))/(x1(i)*x1(i))*d(i,js,k)*G
    grv1(i,js,k) = -G*d(i,js,k)*mc(i)/x1(i)**2d0

   end do
!$omp end parallel do

!$omp workshare
!  grv1(is:ie,js:je,k) = spread(grv1(is:ie,js,k),2,je) *G* d(is:ie,js:je,k)
!   grv1(is,js:je,k) = 0d0 ; grv2 = 0.d0 ; grv3 = 0.d0
!$omp end workshare

  elseif(gravswitch==2.or.gravswitch==3)then
!$omp parallel do private (i,j,k)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      grv1(i,j,k) = -( (dx1(i  )*idx1(i+1)*grvphi(i+1,j,k)    &
                      - dx1(i+1)*idx1(i  )*grvphi(i-1,j,k) )  &
                       /sum(dx1(i:i+1)) &
                     + (dx1(i+1)-dx1(i))*idx1(i)*idx1(i+1)*grvphi(i,j,k) ) &
                   * d(i,j,k)

      grv2(i,j,k) = -( (dx2(j  )*idx2(j+1)*grvphi(i,j+1,k)   &
                      - dx2(j+1)*idx2(j  )*grvphi(i,j-1,k) ) &
                       /sum(dx2(j:j+1)) &
                     + (dx2(j+1)-dx2(j))*idx2(j)*idx2(j+1)*grvphi(i,j,k) ) &
                   * d(i,j,k)
      grv3(i,j,k) = -( (dx3(k  )*idx3(k+1)*grvphi(i,j,k+1) &
                      - dx3(k+1)*idx3(k  )*grvphi(i,j,k-1) ) &
                       /sum(dx3(k:k+1)) &
                     + (dx3(k+1)-dx3(k))*idx3(k)*idx3(k+1)*grvphi(i,j,k) ) &
                   * d(i,j,k)
     end do
    end do
   end do
!$omp end parallel do
   if(eq_sym.and.crdnt==1)grv3(is:ie,js:je,ks) = 0d0
   if(eq_sym.and.crdnt==2)grv2(is:ie,je,ks:ke) = 0d0
!   if(crdnt==1.or.crdnt==2)grv1(is,js,ks:ke) = 0d0

   if(include_extgrv)call externalfield

   if(ie==1)grv1 = 0d0; if(je==1)grv2 = 0d0; if(ke==1)grv3 = 0d0
  else
   print *,"Error from gravswitch (source.f90)"
   stop
  end if

! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(crdnt==0)then
  do k = ks, ke
   do j = js, je
    do i = is, ie
     src(i,j,k,2) = grv1(i,j,k)
     src(i,j,k,3) = grv2(i,j,k)
     src(i,j,k,4) = grv3(i,j,k)
     src(i,j,k,8) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k) &
                  + grv3(i,j,k)*v3(i,j,k)
    end do
   end do
  end do

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(crdnt==1)then
!$omp parallel do private (i,j,k)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     src(i,j,k,2) = ( d(i,j,k)*v2(i,j,k)*v2(i,j,k) - b2(i,j,k)*b2(i,j,k)   &
                    + ptot(i,j,k) ) * sx1(i) &
                  + grv1(i,j,k)

     src(i,j,k,3) = ( b1(i,j,k)*b2(i,j,k) - d(i,j,k)*v1(i,j,k)*v2(i,j,k) ) &
                    * sx1(i) &
                  + grv2(i,j,k) * sx1(i)

     src(i,j,k,4) = grv3(i,j,k)

     src(i,j,k,6) = ( b1(i,j,k)*v3(i,j,k) - b3(i,j,k)*v1(i,j,k) ) * sx1(i)

     src(i,j,k,8) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k)*sx1(i) &
                  + grv3(i,j,k)*v3(i,j,k)

! injection
!     if(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)<0d0.and.imu(i,j,k)*0.7d0>1d0)then
!      src(i,j,k,8) = src(i,j,k,8) + de_dt*d(i,j,k)
!     end if
     if(-grvphi(i,j,k)>v2(i,j,k)*v2(i,j,k).and.x1(i)*x1(i)+x3(k)*x3(k)<0.25d0*sep*sep)then
      src(i,j,k,3) = src(i,j,k,3) + domega_dt *spin_coeffr(i)*d(i,j,k)
      src(i,j,k,8) = src(i,j,k,8) + v2(i,j,k)*domega_dt*spin_coeffr(i)*d(i,j,k)
     end if

    end do
   end do
  end do
!$omp end parallel do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(crdnt==2)then
  do k = ks, ke
   do j = js, je
    do i = is, ie
     src(i,j,k,2) = ( d(i,j,k)*( v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k) ) &
                    - ( b2(i,j,k)*b2(i,j,k) + b3(i,j,k)*b3(i,j,k) )          &
                    + 2.d0*ptot(i,j,k) ) * sx1(i) &
                  + grv1(i,j,k)

     src(i,j,k,3) = ( b1(i,j,k)*b2(i,j,k) - d(i,j,k)*v1(i,j,k)*v2(i,j,k)     &
                  + ( d(i,j,k)*v3(i,j,k)*v3(i,j,k) - b3(i,j,k)*b3(i,j,k)     &
                    + ptot(i,j,k) ) * scot(j) ) * sx1(i) &
                  + grv2(i,j,k)*sx1(i)

     src(i,j,k,4) = ( b1(i,j,k)*b3(i,j,k) - d(i,j,k)*v1(i,j,k)*v3(i,j,k)     &
                  + ( b2(i,j,k)*b3(i,j,k) - d(i,j,k)*v2(i,j,k)*v3(i,j,k) )   &
                  * scot(j) ) * sx1(i) &
                  + grv3(i,j,k)*sx1(i)*sisin(j)

     src(i,j,k,6) = ( b2(i,j,k)*v1(i,j,k) - b1(i,j,k)*v2(i,j,k) ) * sx1(i)

     src(i,j,k,7) = ( (v2(i,j,k)*b3(i,j,k) - b2(i,j,k)*v3(i,j,k))*scot(j)    &
                    - (b1(i,j,k)*v3(i,j,k)-b3(i,j,k)*v1(i,j,k)) ) * sx1(i)

     src(i,j,k,8) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k)*sx1(i) &
                  + grv3(i,j,k)*v3(i,j,k) * sx1(i)*sisin(j)

     if(-grvphi(i,j,k)>v3(i,j,k)*v3(i,j,k).and.x1(i)<60d0*rsun.and.x1(i)>sep*0.25d0)then
      src(i,j,k,4) = src(i,j,k,4) + domega_dt*spin_coeffr(i)*spin_coefft(j)*d(i,j,k)
      src(i,j,k,8) = src(i,j,k,8) + domega_dt*spin_coeffr(i)*spin_coefft(j)*d(i,j,k)*v3(i,j,k)
     end if


    end do
   end do
  end do
 else
  print *,'Error from source.f90'
  stop
 end if

 if(ie==1)src(is:ie,js:je,ks:ke,2) = 0d0
 if(je<=2.and.crdnt==2)src(is:ie,js:je,ks:ke,3) = 0d0
! if(ke==1)src(is:ie,js:je,ks:ke,4) = 0d0

return
end subroutine source
