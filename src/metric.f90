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
 use grid

 integer:: i,j,k
 real(8):: x3_mid,fac_j,fac_k

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

  do k = ks_global-1, ke_global+1
   do j = js_global-1, je_global+1
    do i = is_global-1, ie_global+1
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

  x3_mid = 0d0 !0.5d0*(xi3e+xi3s)
  do i = gis_global-1, gie_global+2
   do k = gks_global-2, gke_global+2
    rdis(i,k) = sqrt( x1(i)**2+(x3(k)-x3_mid)**2 )
    sincyl(i,k) = x1(i)/rdis(i,k)
    coscyl(i,k) = (x3(k)-x3_mid)/rdis(i,k)
   end do
  end do

  do i = is_global-1, ie_global+1
   detg1(i) = xi1(i); idetg1(i) = 2d0 / (xi1(i)**2-xi1(i-1)**2)
   sx1(i) = 2d0 / (xi1(i)+xi1(i-1)) ; g22(i) = x1(i)
  end do
  do j = js_global-1, je_global+1
   do i = is_global-1, ie_global+1
    detg2(i,j) = 1d0; idetg2(i,j) = 1d0 / x1(i) * idxi2(j)
    g33(i,j) = 1d0
   end do
  end do

  do k = ks_global-1, ke_global+1
   do j = js_global-1, je_global+1
    do i = is_global-1, ie_global+1
     idetg3(i,j,k) = idxi3(k)
     dvol(i,j,k) = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j) * dxi3(k)
     sa1(i,j,k) = xi1(i) * dxi2(j) * dxi3(k)
     sa2(i,j,k) = dxi1(i) * dxi3(k)
     sa3(i,j,k) = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j)
    end do
   end do
  end do

  if(je_global==js_global)then
   dvol = dvol * 4d0
   sa1  = sa1  * 4d0
   sa3  = sa3  * 4d0
  end if

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(2)

  do i = is_global-1, ie_global+1
   detg1(i) = xi1(i)**2
   idetg1(i) = 3d0 &
             / ( dxi1(i)*(xi1(i)**2+xi1(i)*xi1(i-1)+xi1(i-1)**2) )
   sx1(i) = 1.5d0 * ( xi1(i)+xi1(i-1) ) &
          / ( xi1(i)**2 + xi1(i)*xi1(i-1) + xi1(i-1)**2 )
   g22(i) = abs(x1(i))
  end do

  do j = js_global-1, je_global+1
   do i = is_global-1, ie_global+1
    detg2(i,j)  = sini(j)
    idetg2(i,j) = 1d0 / (-cosi(j)+cosi(j-1)) * sx1(i)
    g33(i,j) = x1(i) * sinc(j)
   end do
  end do

  fac_k=1d0; if(ke==ks)fac_k=4d0
  fac_j=1d0; if(je==js)fac_j=2d0

  do k = ks_global-1, ke_global+1
   do j = js_global-1, je_global+1
    do i = is_global-1, ie_global+1
     idetg3(i,j,k) = idetg2(i,j) * dxi2(j) * idxi3(k) / fac_k
     dvol(i,j,k)   = (xi1(i)**3-xi1(i-1)**3) / 3d0 &
                   * (cosi(j-1)-cosi(j)) * dxi3(k) * fac_j * fac_k
     sa1(i,j,k)    = 0.5d0*xi1(i)**2 &
                   * (cosi(j-1)-cosi(j)) * dxi3(k) * fac_j * fac_k
     sa2(i,j,k)    = 0.5d0*sini(j) &
                   * (xi1(i)**2-xi1(i-1)**2) * dxi3(k) * fac_k
     sa3(i,j,k)    = 0.5d0 * (xi1(i)**2-xi1(i-1)**2) * dxi2(j) * fac_j
    end do
   end do
  end do

  do j = js_global-1, je_global+1
   scot(j)  = ( sini(j) - sini(j-1) ) / ( cosi(j-1) - cosi(j) )
   sisin(j) = ( xi2(j) - xi2(j-1) ) / ( cosi(j-1) - cosi(j) )
  end do

  if(include_extforce)then
   allocate(spinc_r(is_global:ie_global),spinc_t(js_global:je_global))
   do i = is_global, ie_global
    spinc_r(i) = 0.75d0*(xi1(i)**2+xi1(i-1)**2) &
                       *(xi1(i)+xi1(i-1)) &
                       /(xi1(i)**2+xi1(i)*xi1(i-1)+xi1(i-1)**2)
   end do
   do j = js_global, je_global
    spinc_t(j) = 0.25d0*(2d0*dxi2(j)-sin(2d0*xi2(j)) &
                         +sin(2d0*xi2(j-1)))/(cosi(j-1)-cosi(j))
   end do
  end if

 end select

 return
end subroutine metric

end module metric_mod
