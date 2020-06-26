!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE METRIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate metrics

subroutine metric
  
  use grid
  use physval
  use constants

  implicit none
  
!----------------------------------------------------------------------------

! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==0)then
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
     end do
    end do
   end do

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==1)then
   do i = is-1, ie+1
    detg1(i) = xi1(i); idetg1(i) = 2.d0 / (xi1(i)**2-xi1(i-1)**2)
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
      if(je==1)dvol(i,j,k) = pi * (xi1(i)**2-xi1(i-1)**2) * dxi3(k)
     end do
    end do
   end do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==2)then
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
     end do
    end do
   end do

   if(ke==1) dvol = 4d0 * dvol
   if(ke==1) idetg3 = idetg3 * 0.25d0
   if(je==1) dvol = dvol * 2d0

   do j = js-1, je+1
    scot(j)  = ( sini(j) - sini(j-1) ) / ( cosi(j-1) - cosi(j) )
    sisin(j) = ( xi2(j) - xi2(j-1) ) / ( cosi(j-1) - cosi(j) )
   end do

  end if

return
end subroutine metric
