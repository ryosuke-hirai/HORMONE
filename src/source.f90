module source_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE SOURCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate source terms

subroutine source

 use grid
 use settings,only:eq_sym,include_extforce,wtime,isrc,mag_on
 use physval
 use constants
 use gravmod
 use externalforce_mod
 use omp_lib
 
!----------------------------------------------------------------------------

 wtime(isrc) = wtime(isrc) - omp_get_wtime()

! To calculate gravitational forces ********************************************
 if(gravswitch==0)then ! gravity off
  grv1 = 0d0 ; grv2 = 0d0 ; grv3 = 0d0

 elseif(gravswitch==1.and.crdnt==2.and.dim<=2)then ! point-source
  k = ks
  if(eq_sym)then
   do i = is, ie
    mc(i) = mc(i-1) + sum( d(i,js:je,k) * dvol(i,js:je,k) )*2d0
   end do
  else
   do i = is, ie
    mc(i) = mc(i-1) + sum( d(i,js:je,k) * dvol(i,js:je,k) )
   end do
  end if


!$omp parallel do private(i)
  do i = is, ie
   grv1(i,js:je,k) = -G*d(i,js:je,k)*mc(i)/x1(i)**2
  end do
!$omp end parallel do
  grv2 = 0d0 ; grv3 = 0d0

 elseif(gravswitch==2.or.gravswitch==3)then
!$omp parallel do private (i,j,k) collapse(3)
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
     if(i<=is+sum(fmr_lvl(1:fmr_max))-1)then
      if(i<=is+fmr_lvl(1)-1)then
       grv2(i,j,k) = 0d0;grv3(i,j,k) = 0d0
      else
       fmr_loop: do n = 2, fmr_max
        if(i<=is+sum(fmr_lvl(1:n))-1)then
         grv2(i,j,k) = grv2(i,j,k)/dble(2**(fmr_max-n+1))
         exit fmr_loop
        end if
       end do fmr_loop
      end if
     end if

    end do
   end do
  end do
!$omp end parallel do

  if(eq_sym.and.crdnt==1)grv3(is:ie,js:je,ks) = 0d0

  if(include_extgrv)call externalfield

  if(ie==1)grv1 = 0d0; if(je==1)grv2 = 0d0!; if(ke==1)grv3 = 0d0
!  if(abs(xi1s)<=tiny)grv1(is,js:je,ks:ke) = 0d0
!  if(crdnt==2)grv2(is:ie,js,ks:ke) = 0d0
!  if(crdnt==2)grv2(is:ie,je,ks:ke) = 0d0
 else
  print *,"Error from gravswitch (source.f90)"
  stop
 end if

 select case(crdnt)
! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(0)
!$omp parallel do private(i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(gravswitch>=1)then
      src(i,j,k,imo1) = grv1(i,j,k)
      src(i,j,k,imo2) = grv2(i,j,k)
      src(i,j,k,imo3) = grv3(i,j,k)
      src(i,j,k,iene) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k) &
                      + grv3(i,j,k)*v3(i,j,k)
     end if
    end do
   end do
  end do
!$omp end parallel do

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(1)
!$omp parallel do private (i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     src(i,j,k,imo1) = ( d(i,j,k)*v2(i,j,k)*v2(i,j,k) + ptot(i,j,k) ) * sx1(i)

     src(i,j,k,imo2) = - d(i,j,k)*v1(i,j,k)*v2(i,j,k) * sx1(i)

     if(mag_on)then
      src(i,j,k,imo1) = src(i,j,k,imo1) - b2(i,j,k)**2*sx1(i)
      src(i,j,k,imo2) = src(i,j,k,imo2) + b1(i,j,k)*b2(i,j,k)*sx1(i)
      src(i,j,k,img2) = ( b1(i,j,k)*v3(i,j,k) - b3(i,j,k)*v1(i,j,k) ) * sx1(i)
     end if
     if(gravswitch>=1)then
      src(i,j,k,imo1) = src(i,j,k,imo1) + grv1(i,j,k)
      src(i,j,k,imo2) = src(i,j,k,imo2) + grv2(i,j,k) * sx1(i)
      src(i,j,k,imo3) = grv3(i,j,k)
      src(i,j,k,iene) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k)*sx1(i) &
                      + grv3(i,j,k)*v3(i,j,k)
     end if
    end do
   end do
  end do
!$omp end parallel do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(2)
!$omp parallel do private (i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     src(i,j,k,imo1) = ( d(i,j,k)*( v2(i,j,k)**2 + v3(i,j,k)**2 ) &
                     + 2d0*ptot(i,j,k) ) * sx1(i)

     if(je>js)then
      src(i,j,k,imo2) = (- d(i,j,k)*v1(i,j,k)*v2(i,j,k) &
                      + ( d(i,j,k)*v3(i,j,k)**2 + ptot(i,j,k) ) * scot(j) ) &
                      * sx1(i)
     else
      src(i,j,k,imo2) = 0d0
     end if

     src(i,j,k,imo3) = ( - d(i,j,k)*v1(i,j,k)*v3(i,j,k)     &
                         - d(i,j,k)*v2(i,j,k)*v3(i,j,k) * scot(j) ) * sx1(i)

     if(mag_on)then
      src(i,j,k,imo1) = src(i,j,k,imo1) - (b2(i,j,k)**2+b3(i,j,k)**2)*sx1(i)
      src(i,j,k,imo2) = src(i,j,k,imo2) &
                      + (b1(i,j,k)*b2(i,j,k)-b3(i,j,k)**2*scot(j))*sx1(i)
      src(i,j,k,imo3) = src(i,j,k,imo3) &
                      + (b1(i,j,k)*b3(i,j,k)+b2(i,j,k)*b3(i,j,k)*scot(j))*sx1(i)
      src(i,j,k,img2) = ( b2(i,j,k)*v1(i,j,k) - b1(i,j,k)*v2(i,j,k) ) * sx1(i)
      src(i,j,k,img3) = ( (v2(i,j,k)*b3(i,j,k) - b2(i,j,k)*v3(i,j,k))*scot(j) &
                      - (b1(i,j,k)*v3(i,j,k)-b3(i,j,k)*v1(i,j,k)) ) * sx1(i)
     end if
     if(gravswitch>=1)then
      src(i,j,k,imo1) = src(i,j,k,imo1) + grv1(i,j,k)
      src(i,j,k,imo2) = src(i,j,k,imo2) + grv2(i,j,k)*sx1(i)
      src(i,j,k,imo3) = src(i,j,k,imo3) + grv3(i,j,k)*sx1(i)*sisin(j)
      src(i,j,k,iene) = grv1(i,j,k)*v1(i,j,k) + grv2(i,j,k)*v2(i,j,k)*sx1(i) &
                      + grv3(i,j,k)*v3(i,j,k) * sx1(i)*sisin(j)

     end if
    end do
   end do
  end do
!$omp end parallel do
 case default
  print *,'Error from source.f90'
  stop
 end select

 if(include_extforce)call externalforce
 
 wtime(isrc) = wtime(isrc) + omp_get_wtime()

 return
end subroutine source

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE PHIDAMP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Damping term for divergence cleaning 9-wave method
! Ref    : Mignone & Tzeferacos 2010, JComP, 229, 2117

subroutine phidamp

 use settings,only:courant
 use grid
 use physval,only:phi

 real(8):: alpha9wave

!-----------------------------------------------------------------------------

! ratio between diffusive and advection timescales (td/ta)
 alpha9wave = 0.1d0
 
!$omp parallel do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    phi(i,j,k) = phi(i,j,k)*exp(-alpha9wave*courant)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine phidamp

end module source_mod
