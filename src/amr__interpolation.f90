module amr_interpolation_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE AMR_INTERPOLATION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To interpolate physical values for refined leaf blocks.

subroutine amr_interpolation(lf,olf)

 use settings
 use grid,only:is,js,ks
 use amr_templates
 use amr_module,only:cuts,cib,face,ib,jb,kb
 use pressure_mod

 implicit none

 integer :: lid, i, j, k, ie, je, ke
 type(leaf_contents),intent(inout),dimension(1:cib):: lf
 type(leaf_contents),intent(in) :: olf
 real*8,allocatable,dimension(:,:,:):: uold

!-----------------------------------------------------------------------------

 do lid = 1, cib ! Redefine dx's
  lf(lid)%dx1 = olf%dx1 / dble(min(ib,cuts)); lf(lid)%idx1 = 1.d0 / lf(lid)%dx1
  lf(lid)%dx2 = olf%dx2 / dble(min(jb,cuts)); lf(lid)%idx2 = 1.d0 / lf(lid)%dx2
  lf(lid)%dx3 = olf%dx3 / dble(min(kb,cuts)); lf(lid)%idx3 = 1.d0 / lf(lid)%dx3
 end do

! Set coordinates +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do lid = 1, cib
! x1 direction
  if(mod(lid,cuts)==1)then
   lf(lid)%xi1(is-1) = olf%xi1(is-1)
  else
   lf(lid)%xi1(is-1) = 5.d-1 * ( olf%xi1(is-1) + olf%xi1(ib) )
  end if
  lf(lid)%xi1(is-2) = lf(lid)%xi1(is-1) - lf(lid)%dx1(is-1)
  lf(lid)%x1(is-1)  = lf(lid)%xi1(is-1) - 5.d-1*lf(lid)%dx1(is-1)
  lf(lid)%x1(is-2)  = lf(lid)%x1(is-1)  - lf(lid)%dx1(is-1)
  do i = is, ib+2
   lf(lid)%xi1(i) = lf(lid)%xi1(i-1) + lf(lid)%dx1(i)
   lf(lid)%x1(i)  = lf(lid)%x1(i-1)  + lf(lid)%dx1(i)
  end do
! x2 direction
  if(mod(lid-1,cuts*cuts)<cuts)then
   lf(lid)%xi2(js-1) = olf%xi2(js-1)
  else
   lf(lid)%xi2(js-1) = 5.d-1 * ( olf%xi2(js-1) + olf%xi2(jb) )
  end if
  lf(lid)%xi2(js-2) = lf(lid)%xi2(js-1) - lf(lid)%dx2(js-1)
  lf(lid)%x2(js-1)  = lf(lid)%xi2(js-1) - 5.d-1*lf(lid)%dx2(js-1)
  lf(lid)%x2(js-2)  = lf(lid)%x2(js-1)  - lf(lid)%dx2(js-1)
  do j = js, jb+2
   lf(lid)%xi2(j) = lf(lid)%xi2(j-1) + lf(lid)%dx2(j)
   lf(lid)%x2(j)  = lf(lid)%x2(j-1)  + lf(lid)%dx2(j)
  end do
! x3 direction
  if(lid<=cuts*cuts)then
   lf(lid)%xi3(ks-1) = olf%xi3(ks-1)
  else
   lf(lid)%xi3(ks-1) = 5.d-1 * ( olf%xi3(ks-1) + olf%xi3(kb) )
  end if
  lf(lid)%xi3(ks-2) = lf(lid)%xi3(ks-1) - lf(lid)%dx3(ks-1)
  lf(lid)%x3(ks-1)  = lf(lid)%xi3(ks-1) - 5.d-1*lf(lid)%dx3(ks-1)
  lf(lid)%x3(ks-2)  = lf(lid)%x3(ks-1)  - lf(lid)%dx3(ks-1)
  do k = ks, kb+2
   lf(lid)%xi3(k) = lf(lid)%xi3(k-1) + lf(lid)%dx3(k)
   lf(lid)%x3(k)  = lf(lid)%x3(k-1)  + lf(lid)%dx3(k)
  end do
 end do

! Set metrics +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 do lid = 1, cib
  call amr_block_metric(lf(lid))  
 end do

! Start interpolation of physical values ++++++++++++++++++++++++++++++++++++

 allocate(uold(0:max(ib/cuts,1)+1,0:max(jb/cuts,1)+1,0:max(kb/cuts,1)+1))
 do lid = 1, cib
  i = ib/cuts*mod(lid-1,cuts)            ; ie = i+max(ib/cuts,1) + 1
  j = jb/cuts*mod((lid-1)/cuts,cuts)     ; je = j+max(jb/cuts,1) + 1
  k = kb/cuts*mod((lid-1)/cuts**2,cuts)  ; ke = k+max(kb/cuts,1) + 1

! density
  uold = olf%d(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%d,lf(lid)%dvol,uold)

! velocity
  uold = olf%d(i:ie,j:je,k:ke) * olf%v1(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%v1,lf(lid)%dvol,uold)
  lf(lid)%v1 = lf(lid)%v1 / lf(lid)%d

  uold = olf%d(i:ie,j:je,k:ke) * olf%v2(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%v2,lf(lid)%dvol,uold)
  lf(lid)%v2 = lf(lid)%v2 / lf(lid)%d

  uold = olf%d(i:ie,j:je,k:ke) * olf%v3(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%v3,lf(lid)%dvol,uold)
  lf(lid)%v3 = lf(lid)%v3 / lf(lid)%d
! magnetic field
  uold = olf%b1(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%b1,lf(lid)%dvol,uold)

  uold = olf%b2(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%b2,lf(lid)%dvol,uold)

  uold = olf%b3(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%b3,lf(lid)%dvol,uold)
! internal energy
  uold = olf%e(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%e,lf(lid)%dvol,uold)
! phi (divergence cleaning)
  uold = olf%phi(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%phi,lf(lid)%dvol,uold)
! gravitational potential
  uold = olf%gphi(i:ie,j:je,k:ke)
  call amr_cellintpol(lf(lid)%gphi,lf(lid)%dvol,uold)

 end do
! pressure and total pressure and fluxes
 do lid = 1, cib
  call pressure_block(lf(lid))
  lf(lid)%u(:,:,:,1) = lf(lid)%d
  lf(lid)%u(:,:,:,2) = lf(lid)%d * lf(lid)%v1
  lf(lid)%u(:,:,:,3) = lf(lid)%d * lf(lid)%v2
  lf(lid)%u(:,:,:,4) = lf(lid)%d * lf(lid)%v3
  lf(lid)%u(:,:,:,5) = lf(lid)%b1
  lf(lid)%u(:,:,:,6) = lf(lid)%b2
  lf(lid)%u(:,:,:,7) = lf(lid)%b3
  lf(lid)%u(:,:,:,8) = lf(lid)%e
  lf(lid)%u(:,:,:,9) = lf(lid)%phi
 end do
 deallocate(uold)

contains



!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_CELLINTPOL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To interpolate physical values into a cell.

subroutine amr_cellintpol(lfu,dvol,u)

 use amr_module,only:cib,ib,jb,kb,cuts

 implicit none

 integer i,j,k,lid, il, jl, kl
 real*8,allocatable,dimension(:,:,:),intent(in) :: u, dvol
 real*8,allocatable,dimension(:,:,:),intent(inout):: lfu
 real*8 :: dltp, dltm
 real*8,dimension(3):: dlt

!-----------------------------------------------------------------------------

 if(cuts/=2)then
  print *,"Error: Interpolation scheme for cuts > 2 is needed",cuts
  stop
 end if

 do k = 1, max(kb/2,1)
  do j = 1, max(jb/2,1)
   do i = 1, max(ib/2,1)
    dltp = u(i+1,j,k) - u(i,j,k) ; dltm = u(i,j,k) - u(i-1,j,k)
    dlt(1) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

    dltp = u(i,j+1,k) - u(i,j,k) ; dltm = u(i,j,k) - u(i,j-1,k)
    dlt(2) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

    dltp = u(i,j,k+1) - u(i,j,k) ; dltm = u(i,j,k) - u(i,j,k-1)
    dlt(3) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

    do lid = 1, cib
     il = 2*i-1+mod(lid-1,2) 
     jl = 2*j-1+mod((lid-1)/2,2)
     kl = 2*k-1+mod((lid-1)/2**2,2)
     lfu(il,jl,kl) = u(i,j,k) &
          + dble((-1)**il) * 0.25d0 * dlt(1) * dvol(2*i-1,jl,kl)/dvol(il,jl,kl)&
          + dble((-1)**jl) * 0.25d0 * dlt(2) * dvol(il,2*j-1,kl)/dvol(il,jl,kl)&
          + dble((-1)**kl) * 0.25d0 * dlt(3) * dvol(il,jl,2*k-1)/dvol(il,jl,kl)
    end do
   end do
  end do
 end do

return
end subroutine amr_cellintpol
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

end subroutine amr_interpolation




!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_BLOCK_METRIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate metrics for each leaf block.

subroutine amr_block_metric(lf)

 use settings,only:crdnt
 use grid,only:is,js,ks
 use constants
 use amr_templates
 use amr_module,only:ib,jb,kb

 implicit none

 integer i,j,k
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==0)then
   do i = is-1, ib+1
    lf%detg1(i) = 1.d0; lf%idetg1(i) = lf%idx1(i); lf%g22(i) = 1.d0
   end do

   do j = js-1, jb+1
    do i = is-1, ib+1
     lf%detg2(i,j) = 1.d0; lf%idetg2(i,j) = lf%idx2(j); lf%g33(i,j) = 1.d0
    end do
   end do

   do k = ks-1, kb+1
    do j = js-1, jb+1
     do i = is-1, ib+1
      lf%idetg3(i,j,k) = lf%idx3(k)      
     end do
    end do
   end do

   do k = ks-2, kb+2
    do j = js-2, jb+2
     do i = is-2, ib+2
      lf%dvol(i,j,k)   = lf%dx1(i) * lf%dx2(j) * lf%dx3(k)
     end do
    end do
   end do

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==1)then
   do i = is-1, ib+1
    lf%detg1(i) = lf%xi1(i)
    lf%idetg1(i) = 2.d0 / (lf%xi1(i)**2.-lf%xi1(i-1)**2.)
    lf%sx1(i) = 1.d0 / lf%x1(i); lf%g22(i) = lf%x1(i)
   end do
   do j = js-1, jb+1
    do i = is-1, ib+1
     lf%detg2(i,j) = 1.d0; lf%idetg2(i,j) = 1.d0 / lf%x1(i) * lf%idx2(j)
     lf%g33(i,j) = 1.d0
    end do
   end do

   do k = ks-1, kb+1
    do j = js-1, jb+1
     do i = is-1, ib+1
      lf%idetg3(i,j,k) = lf%idx3(k)      
     end do
    end do
   end do

   do k = ks-2, kb+2
    do j = js-2, jb+2
     do i = is-2, ib+2
      lf%dvol(i,j,k) = 5.d-1 * (lf%xi1(i)**2.d0-lf%xi1(i-1)**2.d0) &
                             * lf%dx2(j) * lf%dx3(k)
      if(jb==1)&
       lf%dvol(i,j,k) = pi * (lf%xi1(i)**2.d0-lf%xi1(i-1)**2.d0) * lf%dx3(k)
     end do
    end do
   end do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==2)then
   do i = is-1, ib+1
    lf%detg1(i) = lf%xi1(i)**2.d0
    lf%idetg1(i) = 3.d0 / ( (lf%xi1(i)-lf%xi1(i-1))*&
                  (lf%xi1(i)**2.+lf%xi1(i)*lf%xi1(i-1)+lf%xi1(i-1)**2.) )
    lf%sx1(i) = 3.d0 / 2.d0 * ( lf%xi1(i)+lf%xi1(i-1) ) &
   / ( lf%xi1(i)*lf%xi1(i) + lf%xi1(i)*lf%xi1(i-1) + lf%xi1(i-1)*lf%xi1(i-1) )
    lf%g22(i) = lf%x1(i)
   end do

   do j = js-1, jb+1
    do i = is-1, ib+1
     lf%detg2(i,j)  = sin(lf%xi2(j))
     lf%idetg2(i,j) = 1.d0 / (-cos(lf%xi2(j))+cos(lf%xi2(j-1))) * lf%sx1(i)
     lf%g33(i,j) = lf%x1(i) * sin(lf%x2(j))
    end do
   end do

   do k = ks-1, kb+1
    do j = js-1, jb+1
     do i = is-1, ib+1
      lf%idetg3(i,j,k) = 1.d0 / sin(lf%x2(j)) * lf%sx1(i) * lf%idx3(k)
     end do
    end do
   end do

   do k = ks-2, kb+2
    do j = js-2, jb+2
     do i = is-2, ib+2
      lf%dvol(i,j,k) = lf%dx1(i) &
           * (3.d0*lf%xi1(i)*(lf%xi1(i)-lf%dx1(i))+lf%dx1(i)*lf%dx1(i))/3.d0 &
           * (cos(lf%xi2(j)-lf%dx2(j))-cos(lf%xi2(j))) * lf%dx3(k)
      if(kb==1)lf%dvol(i,j,k) = lf%dvol(i,j,k) * 4.d0
     end do
    end do
   end do

   do j = js-1, jb+1
    lf%scot(j)  = ( sin(lf%xi2(j)) - sin(lf%xi2(j-1)) ) &
                / ( cos(lf%xi2(j-1)) - cos(lf%xi2(j)) )
    lf%sisin(j) = ( lf%xi2(j) - lf%xi2(j-1) ) &
                / ( cos(lf%xi1(j-1)) - cos(lf%xi2(j)) )
   end do

  end if
 

return
end subroutine amr_block_metric


end module amr_interpolation_mod
