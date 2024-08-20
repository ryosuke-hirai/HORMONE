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
 use settings,only:eq_sym,include_extforce,mag_on,radswitch,include_sinks
 use physval
 use constants
 use utils,only:masscoordinate
 use gravmod
 use sink_mod,only:sinkfield
 use radiation_mod,only:radiative_force
 use externalforce_mod
 use profiler_mod

 integer:: i,j,k

!----------------------------------------------------------------------------

 call start_clock(wtsrc)

! To calculate gravitational forces ********************************************
 if(gravswitch==0)then ! gravity off
! Don't do anything

 elseif(gravswitch==1.and.crdnt==2)then ! point-source
  call masscoordinate

!$omp parallel do private(i)
  do i = is, ie
   grv1(i,js:je,ks:ke) = -G*d(i,js:je,ks:ke)*mc(i)/x1(i)**2
  end do
!$omp end parallel do

 elseif(gravswitch==2.or.gravswitch==3)then

!$omp parallel do private(i,j,k) collapse(3)
  do k = ks-1, ke+1
   do j = js-1, je+1
    do i = is-1, ie+1
     totphi(i,j,k) = grvphi(i,j,k)
    end do
   end do
  end do
!$omp end parallel do

  if(include_extgrv)call externalfield
  if(include_sinks) call sinkfield

  call get_fieldforce(totphi,d,grv1,grv2,grv3)

  if(eq_sym.and.crdnt==1)grv3(is:ie,js:je,ks) = 0d0

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

 if(radswitch>0)     call radiative_force
 if(include_extforce)call externalforce

 call stop_clock(wtsrc)

 return
end subroutine source

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE GET_FIELDFORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute the force from a field

subroutine get_fieldforce(phi,coeff,frc1,frc2,frc3)

 use grid

 real(8),allocatable,dimension(:,:,:),intent(in):: phi, coeff
 real(8),allocatable,dimension(:,:,:),intent(inout):: frc1,frc2,frc3
 integer:: i,j,k,n

!-----------------------------------------------------------------------------

!$omp parallel do private (i,j,k,n) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     frc1(i,j,k) = -( (dx1(i  )*idx1(i+1)*phi(i+1,j,k)    &
                     - dx1(i+1)*idx1(i  )*phi(i-1,j,k) )  &
                      /sum(dx1(i:i+1)) &
                    + (dx1(i+1)-dx1(i))*idx1(i)*idx1(i+1)*phi(i,j,k) ) &
                  * coeff(i,j,k)

     frc2(i,j,k) = -( (dx2(j  )*idx2(j+1)*phi(i,j+1,k)   &
                     - dx2(j+1)*idx2(j  )*phi(i,j-1,k) ) &
                      /sum(dx2(j:j+1)) &
                    + (dx2(j+1)-dx2(j))*idx2(j)*idx2(j+1)*phi(i,j,k) ) &
                  * coeff(i,j,k)

     frc3(i,j,k) = -( (dx3(k  )*idx3(k+1)*phi(i,j,k+1) &
                     - dx3(k+1)*idx3(k  )*phi(i,j,k-1) ) &
                      /sum(dx3(k:k+1)) &
                    + (dx3(k+1)-dx3(k))*idx3(k)*idx3(k+1)*phi(i,j,k) ) &
                  * coeff(i,j,k)
     if(fmr_max==0)cycle
     if(i<=is_global+sum(fmr_lvl(1:fmr_max))-1)then
      if(i<=is_global+fmr_lvl(1)-1)then
       frc2(i,j,k) = 0d0;frc3(i,j,k) = 0d0
      else
       fmr_loop: do n = 2, fmr_max
        if(i<=is_global+sum(fmr_lvl(1:n))-1)then
         frc2(i,j,k) = frc2(i,j,k)/dble(2**(fmr_max-n+1))
         frc3(i,j,k) = frc3(i,j,k)/dble(2**(fmr_max-n+1))
         exit fmr_loop
        end if
       end do fmr_loop
      end if
     end if

    end do
   end do
  end do
!$omp end parallel do

return
end subroutine get_fieldforce

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE PHIDAMP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Damping term for divergence cleaning 9-wave method
! Ref    : Mignone & Tzeferacos 2010, JComP, 229, 2117

subroutine phidamp

 use settings,only:courant,alpha9wave
 use grid
 use physval,only:phi

 integer:: i,j,k

!-----------------------------------------------------------------------------

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
