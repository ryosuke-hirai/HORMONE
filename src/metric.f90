module metric_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE METRIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate metrics

subroutine metric

 use settings,only:include_extforce
 use constants,only:pi
 use grid

 integer:: i,j,k

!----------------------------------------------------------------------------

 select case(crdnt)
! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(0)
  do i = is-1, ie+1
   detg1(i) = 1d0; idetg1(i) = idxi1(i); g22(i) = 1d0
  end do

  do j = js-1, je+1
   do i = is-1, ie+1
    detg2(i,j) = 1d0; idetg2(i,j) = idxi2(j); g33(i,j) = 1d0
   end do
  end do

  do k = ks-1, ke+1
   do j = js-1, je+1
    do i = is-1, ie+1
     idetg3(i,j,k) = idxi3(k)
     dvol(i,j,k)   = dxi1(i) * dxi2(j) * dxi3(k)
     sa1(i,j,k)   = dxi2(j) * dxi3(k)
     sa2(i,j,k)   = dxi1(i) * dxi3(k)
     sa3(i,j,k)   = dxi1(i) * dxi2(j)
    end do
   end do
  end do

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(1)

  allocate( rdis(-1:gie+2,gks-2:gke+2) )
  allocate( sincyl,coscyl,mold=rdis )
  do i = gis-1, gie+2
   do k = gks-2, gke+2
    rdis(i,k) = sqrt( x1(i)**2+x3(k)**2 )
    sincyl(i,k) = x1(i)/rdis(i,k)
    coscyl(i,k) = x3(k)/rdis(i,k)
   end do
  end do

  do i = is-1, ie+1
   detg1(i) = xi1(i); idetg1(i) = 2d0 / (xi1(i)**2-xi1(i-1)**2)
   sx1(i) = 2d0 / (xi1(i)+xi1(i-1)) ; g22(i) = x1(i)
  end do
  do j = js-1, je+1
   do i = is-1, ie+1
    detg2(i,j) = 1d0; idetg2(i,j) = 1d0 / x1(i) * idxi2(j)
    g33(i,j) = 1d0
   end do
  end do

  do k = ks-1, ke+1
   do j = js-1, je+1
    do i = is-1, ie+1
     idetg3(i,j,k) = idxi3(k)
     dvol(i,j,k) = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j) * dxi3(k)
     if(je==js)dvol(i,j,k) = pi * (xi1(i)**2-xi1(i-1)**2) * dxi3(k)
     sa1(i,j,k) = 0.5d0 * xi1(i) * dxi2(j) * dxi3(k)
     sa2(i,j,k) = 0.5d0 * dxi1(i) * dxi3(k)
     sa3(i,j,k) = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j)
    end do
   end do
  end do

  if(je==js)then
   dvol = dvol * 4d0 
   sa1  = sa1  * 4d0 
   sa3  = sa3  * 4d0 
  end if

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(2)

  allocate( sinc, sini, cosc, cosi, mold=x2 )
  do j = js-2, je+2
   sinc(j)=sin(x2 (j))
   sini(j)=sin(xi2(j))
   cosc(j)=cos(x2 (j))
   cosi(j)=cos(xi2(j))
  end do

  do i = is-1, ie+1
   detg1(i) = xi1(i)**2
   idetg1(i) = 3d0 &
             / ( dxi1(i)*(xi1(i)**2+xi1(i)*xi1(i-1)+xi1(i-1)**2) )
   sx1(i) = 1.5d0 * ( xi1(i)+xi1(i-1) ) &
          / ( xi1(i)**2 + xi1(i)*xi1(i-1) + xi1(i-1)**2 )
   g22(i) = x1(i)
  end do

  do j = js-1, je+1
   do i = is-1, ie+1
    detg2(i,j)  = sini(j)
    idetg2(i,j) = 1d0 / (-cosi(j)+cosi(j-1)) * sx1(i)
    g33(i,j) = x1(i) * sinc(j)
   end do
  end do
  
  do k = ks-1, ke+1
   do j = js-1, je+1
    do i = is-1, ie+1
     idetg3(i,j,k) = 1d0 / sin(x2(j)) * sx1(i) * idxi3(k)
     dvol(i,j,k)   = (xi1(i)**3-xi1(i-1)**3) / 3d0 &
                   * (cosi(j-1)-cosi(j)) * dxi3(k)
     sa1(i,j,k)    = 0.5d0*(xi1(i)**2) &
                   * (cosi(j-1)-cosi(j)) * dxi3(k)
     sa2(i,j,k)    = 0.5d0*sini(j) &
                   * (xi1(i)**2-xi1(i-1)**2) * dxi3(k)
     sa3(i,j,k)    = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j)
    end do
   end do
  end do

  if(ke==ks)then
   dvol = dvol * 4d0
   idetg3 = idetg3 * 0.25d0
   sa1 = sa1 * 4d0 
   sa2 = sa2 * 4d0 
  end if
  if(je==js)then
   dvol = dvol * 2d0
   sa1  = sa1  * 2d0
   sa3  = sa3  * 2d0
  end if

  do j = js-1, je+1
   scot(j)  = ( sini(j) - sini(j-1) ) / ( cosi(j-1) - cosi(j) )
   sisin(j) = ( xi2(j) - xi2(j-1) ) / ( cosi(j-1) - cosi(j) )
  end do

  if(include_extforce)then
   allocate(spinc_r(is:ie),spinc_t(js:je))
   do i = is, ie
    spinc_r(i) = 0.75d0*(xi1(i)**2+xi1(i-1)**2) &
                       *(xi1(i)+xi1(i-1)) &
                       /(xi1(i)**2+xi1(i)*xi1(i-1)+xi1(i-1)**2)
   end do
   do j = js, je
    spinc_t(j) = 0.25d0*(2d0*dxi2(j)-sin(2d0*xi2(j)) &
                         +sin(2d0*xi2(j-1)))/(cosi(j-1)-cosi(j))
   end do
  end if
   
 end select

 return
end subroutine metric

end module metric_mod
