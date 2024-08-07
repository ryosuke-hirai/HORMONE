module gridset_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRIDSET
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set grid

subroutine gridset

 use settings,only:imesh,jmesh,kmesh,eq_sym,start,gravswitch,crdnt
 use grid
 use constants,only:pi
 use utils,only:geometrical_series
 use readbin_mod,only:readgrid

 integer::i,j,k,jetmp,ketmp

!-------------------------------------------------------------------------

 frame_acc = 0d0

 if(gravswitch/=0)then
  deallocate(x1, xi1, dx1, dxi1, idx1, idxi1, &
             x2, xi2, dx2, dxi2, idx2, idxi2, &
             x3, xi3, dx3, dxi3, idx3, idxi3 )
  allocate( &
   x1(gis_global-2:gie_global+2), xi1(gis_global-2:gie_global+2), dx1(gis_global-2:gie_global+2), &
   dxi1(gis_global-2:gie_global+2), idx1(gis_global-2:gie_global+2), idxi1(gis_global-2:gie_global+2), &
   x2(gjs_global-2:gje_global+2), xi2(gjs_global-2:gje_global+2), dx2(gjs_global-2:gje_global+2), &
   dxi2(gjs_global-2:gje_global+2), idx2(gjs_global-2:gje_global+2), idxi2(gjs_global-2:gje_global+2), &
   x3(gks_global-2:gke_global+2), xi3(gks_global-2:gke_global+2), dx3(gks_global-2:gke_global+2), &
   dxi3(gks_global-2:gke_global+2), idx3(gks_global-2:gke_global+2), idxi3(gks_global-2:gke_global+2) &
  )
 end if

 if(start==0)then
! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  coordinate_system: select case (crdnt)
  case(0) coordinate_system

! set x direction
   mesh_car1: select case (imesh)
   case(0) mesh_car1 ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie_global-is_global+1)
   case(1) mesh_car1 ! geometrical series
    call geometrical_series(dxi1,x1min,is_global,ie_global,xi1s,xi1e)
   case(2) mesh_car1 ! user specified mesh
    call other_imesh(dxi1,is_global,ie_global,xi1s,xi1e)
   end select mesh_car1

   do i = is_global-1,ie_global+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is_global-1)  = xi1s - dxi1(is_global-1)*0.5d0
   x1(is_global-2)  = x1(is_global-1) - dx1(is_global-1)
   xi1(is_global-1) = xi1s
   xi1(is_global-2) = xi1s - dxi1(is_global-1)
   do i = is_global,ie_global+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
    dx1 = dxi1 ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
    do i = is_global-1, gis_global-2, -1
     xi1(i) = xi1(i+1) - dxi1(i)
     x1 (i) = x1 (i+1) - dx1 (i)
    end do
    do i = ie_global+1, gie_global+2
     xi1(i) = xi1(i-1) + dxi1(i)
     x1 (i) = x1 (i-1) + dx1 (i)
    end do
   end if

! set y direction
   mesh_car2: select case (jmesh)
   case(0) mesh_car2 ! uniform mesh
    dxi2 = (xi2e-xi2s) / dble(je_global-js_global+1)
   case(1) mesh_car2 ! geometrical series
    call geometrical_series(dxi2,x2min,js_global,je_global,xi2s,xi2e)
   case(2) mesh_car2 ! user specified mesh
    call other_jmesh(dxi2,js_global,je_global,xi2s,xi2e)
   end select mesh_car2

   do j = js_global-1,je_global+2
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js_global-1)  = xi2s - dxi2(js_global-1)*0.5d0
   x2(js_global-2)  = x2(js_global-1) - dx2(js_global-1)
   xi2(js_global-1) = xi2s
   xi2(js_global-2) = xi2s - dxi2(js_global-1)
   do j = js_global,je_global+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
    dx2 = dxi2 ; idx2 = 1d0 / dx2 ; idxi2 = 1d0 / dxi2
    do j = js_global-1, gjs_global-2, -1
     xi2(j) = xi2(j+1) - dxi2(j)
     x2 (j) = x2 (j+1) - dx2 (j)
    end do
    do j = je_global+1, gje_global+2
     xi2(j) = xi2(j-1) + dxi2(j)
     x2 (j) = x2 (j-1) + dx2 (j)
    end do
   end if

! set z direction
   mesh_car3: select case (kmesh)
   case(0) mesh_car3 ! uniform mesh
    dxi3 = (xi3e-xi3s) / dble(ke_global-ks_global+1)
   case(1) mesh_car3 ! geometrical series
    call geometrical_series(dxi3,x3min,ks_global,ke_global,xi3s,xi3e)
   case(2) mesh_car3 ! user specified mesh
    call other_kmesh(dxi3,ks_global,ke_global,xi3s,xi3e)
   end select mesh_car3

   do k = ks_global-1,ke_global+2
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks_global-1)  = xi3s - dxi3(ks_global-1)*0.5d0
   x3(ks_global-2)  = x3(ks_global-1) - dx3(ks_global-1)
   xi3(ks_global-1) = xi3s
   xi3(ks_global-2) = xi3s - dxi3(ks_global-1)
   do k = ks_global,ke_global+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
    dx3 = dxi3 ; idx3 = 1d0 / dx3 ; idxi3 = 1d0 / dxi3
    do k = ks_global-1, gks_global-2, -1
     xi3(k) = xi3(k+1) - dxi3(k)
     x3 (k) = x3 (k+1) - dx3 (k)
    end do
    do k = ke_global+1, gke_global+2
     xi3(k) = xi3(k-1) + dxi3(k)
     x3 (k) = x3 (k-1) + dx3 (k)
    end do
   end if

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  case(1) coordinate_system
   xi2s = 0.d0 ; xi2e = 2.d0*pi
! set r direction
   mesh_cyl1: select case (imesh)
   case(0) mesh_cyl1 ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie_global-is_global+1)
   case(1) mesh_cyl1 ! geometrical series
    call geometrical_series(dxi1,x1min,is_global,ie_global,xi1s,xi1e)
   case(2) mesh_cyl1 ! user specified mesh
    call other_imesh(dxi1,is_global,ie_global,xi1s,xi1e)
   end select mesh_cyl1

   do i = is_global-1,ie_global+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is_global-1)  = xi1s - dxi1(is_global-1)*0.5d0
   x1(is_global-2)  = x1(is_global-1) - dx1(is_global-1)
   xi1(is_global-1) = xi1s
   xi1(is_global-2) = xi1s - dxi1(is_global-1)
   do i = is_global,ie_global+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do


! for volumetric centre
   do i = is_global-1, ie_global+2
    x1(i)  = sqrt( (xi1(i-1)*xi1(i-1)+xi1(i)*xi1(i)) *5d-1 )
    if(i==is_global-1)x1(i)=-x1(i)
    dx1(i) = x1(i) - x1(i-1)
    idx1(i) = 1d0 / dx1(i)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
    do i = ie_global+1, gie_global+2
     dxi1(i) = dxi1(i-1) * dxi1(ie)/dxi1(ie-1)
     xi1(i) = xi1(i-1) + dxi1(i)
     x1(i)  = 0.5d0*(xi1(i)+xi1(i-1))!sqrt( (xi1(i-1)*xi1(i-1)+xi1(i)*xi1(i)) *5d-1 )
     dx1(i) = x1(i) - x1(i-1)
     idx1(i) = 1d0 / dx1(i)
    end do
    dx1(gis_global-2) = dxi1(gis_global-2) ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
   end if

! set theta direction
   if    (je_global==js_global)then
    jetmp=4
   elseif(je_global>=4)then
    jetmp=je_global
   else
    print *,"Error from je",je
   endif
   do j = js_global-1,je_global+2
    dxi2(j) = 2.d0*pi / dble(jetmp)
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js_global-1)  = xi2s - dxi2(js_global-1)*0.5d0
   x2(js_global-2)  = x2(js_global-1) - dx2(js_global-1)
   xi2(js_global-1) = xi2s
   xi2(js_global-2) = xi2s - dxi2(js_global-1)
   do j = js_global,je_global+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do
   xi2e = xi2(je_global)

! set z direction
   mesh_cyl3: select case (kmesh)
   case(0) mesh_cyl3 ! uniform mesh
    dxi3 = (xi3e-xi3s) / dble(ke_global-ks_global+1)
   case(1) mesh_cyl3 ! geometrical series
    call geometrical_series(dxi3,x3min,ks_global,ke_global,xi3s,xi3e)
   case(2) mesh_cyl3 ! user specified mesh
    call other_kmesh(dxi3,ks_global,ke_global,xi3s,xi3e)
   end select mesh_cyl3

   do k = ks_global-1,ke_global+2
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks_global-1)  = xi3s - dxi3(ks_global-1)*0.5d0
   x3(ks_global-2)  = x3(ks_global-1) - dx3(ks_global-1)
   xi3(ks_global-1) = xi3s
   xi3(ks_global-2) = xi3s - dxi3(ks_global-1)
   do k = ks_global,ke_global+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
!    dx3 = dxi3 ; idx3 = 1d0 / dx3 ; idxi3 = 1d0 / dxi3
    do k = ks_global-1, gks_global-2, -1
     dxi3(k) = dxi3(k+1) * dxi3(ks)/dxi3(ks+1)
     idxi3(k)= 1d0 / dxi3(k)
    end do
    do k = ks_global-1, gks_global-1, -1
     dx3(k) = 0.5d0*(dxi3(k-1) + dxi3(k))
     idx3(k) = 1d0/dx3(k)
     xi3(k) = xi3(k+1) - dxi3(k+1)
     x3 (k) = x3 (k+1) - dx3 (k+1)
    end do
    do k = ke_global+1, gke_global+2
     dxi3(k) = dxi3(k-1) * dxi3(ke_global)/dxi3(ke_global-1)
     idxi3(k)= 1d0 / dxi3(k)
     dx3(k) = 0.5d0*(dxi3(k-1) + dxi3(k))
     idx3(k)= 1d0 / dx3(k)
     xi3(k) = xi3(k-1) + dxi3(k)
     x3 (k) = x3 (k-1) + dx3 (k)
    end do
    if(eq_sym)then
     dxi3(gks_global-2:gks_global-1) = dxi3(gks_global+1:gks_global:-1)
     idxi3(gks_global-2:gks_global-1) = 1d0/dxi3(gks_global-2:gks_global-1)
     dx3(gks_global) = 0.5d0*(dxi3(gks_global-1)+dxi3(gks_global))
     dx3(gks_global-1) = 0.5d0*(dxi3(gks_global-2)+dxi3(gks_global-1))
     idx3(gis_global:gks_global-1) = 1d0/dx3(gks_global:gks_global-1)
     xi3(gks_global-2) = xi3(gks_global-1)-dxi3(gks_global-1)
     x3(gks_global-1) = x3(gks_global) - dx3(gks_global)
     x3(gks_global-2) = x3(gks_global-1) - dx3(gks_global-1)
    end if
   end if



! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  case(2) coordinate_system
   xi2s = 0.d0 ; xi2e = pi ; if(je==js.or.eq_sym) xi2e = pi * 0.5d0
   xi3s = -pi  ; xi3e = pi
! set r direction
   mesh_sph1: select case (imesh)
   case(0) mesh_sph1 ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie_global-is_global+1)
   case(1) mesh_sph1 ! geometrical series
    call geometrical_series(dxi1,x1min,is_global,ie_global,xi1s,xi1e)
   case(2) mesh_sph1 ! user specified mesh
    call other_imesh(dxi1,is_global,ie_global,xi1s,xi1e)
   end select mesh_sph1

   do i = is_global-1,ie_global+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is_global-1)  = xi1s - dxi1(is_global-1)*0.5d0
   x1(is_global-2)  = x1(is_global-1) - dx1(is_global-1)
   xi1(is_global-1) = xi1s
   xi1(is_global-2) = xi1s - dxi1(is_global-1)
   do i = is_global,ie_global+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do

! for volumetric centre
   do i = is_global-1, ie_global+2
    x1(i) = 0.75d0*(xi1(i)+xi1(i-1))*(xi1(i)*xi1(i)+xi1(i-1)*xi1(i-1)) &
                  /(xi1(i)*xi1(i)+xi1(i)*xi1(i-1)+xi1(i-1)*xi1(i-1))
    if(i==is_global)x1(i-1)=-x1(i)
    dx1(i) = x1(i) - x1(i-1)
    idx1(i) = 1d0 / dx1(i)
   end do

! extend grid for gravitational potential
   if(gravswitch/=0)then
    do i = ie_global+1, gie_global+2
     dxi1(i) = dxi1(i-1) * dxi1(ie)/dxi1(ie-1)
     idxi1(i) = 1d0/dxi1(i)
     xi1(i) = xi1(i-1) + dxi1(i)
     x1(i) = 0.75d0*(xi1(i)+xi1(i-1))*(xi1(i)*xi1(i)+xi1(i-1)*xi1(i-1)) &
                   /(xi1(i)*xi1(i)+xi1(i)*xi1(i-1)+xi1(i-1)*xi1(i-1))
     dx1(i) = x1(i) - x1(i-1)
     idx1(i) = 1d0 / dx1(i)
    end do
    dx1(gis_global-2) = dxi1(gis_global-2) ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
   end if

!!$   if(gravswitch/=0)then
!!$    dx1 = dxi1 ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
!!$    do i = ie+1, gie+2
!!$     xi1(i) = xi1(i-1) + dxi1(i)
!!$     x1(i)  = sqrt( (xi1(i-1)**3d0+xi1(i)**3d0) *(1d0/3d0) )
!!$     dx1(i) = x1(i) - x1(i-1)
!!$     idx1(i) = 1d0 / dx1(i)
!!$    end do
!!$   end if

! set theta direction
   mesh_sph2: select case (jmesh)
   case(0) mesh_sph2 ! uniform theta
    if    (je_global==js_global)then
     jetmp=1
    elseif(je_global>=js_global+1)then
     jetmp=je_global
    endif
    dxi2(js_global-2:je_global+2) = xi2e / dble(jetmp)
   case(2) mesh_sph2 ! user specified mesh
    call other_jmesh(dxi2,js_global,je_global,xi2s,xi2e)
   end select mesh_sph2

   do j = js_global-1,je_global+2
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js_global-1)  = xi2s - dxi2(js_global-1)*0.5d0
   x2(js_global-2)  = x2(js_global-1) - dx2(js_global-1)
   xi2(js_global-1) = xi2s
   xi2(js_global-2) = xi2s-dxi2(js_global-1)
   do j = js_global,je_global+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do
   xi2e = xi2(je_global)

  allocate( sinc, sini, cosc, cosi, mold=x2 )
  do j = js_global-2, je_global+2
   sinc(j)=real(sin(real(x2 (j),kind=16)),kind=8)
   sini(j)=real(sin(real(xi2(j),kind=16)),kind=8)
   cosc(j)=real(cos(real(x2 (j),kind=16)),kind=8)
   cosi(j)=real(cos(real(xi2(j),kind=16)),kind=8)
  end do

! for volumetric centre
   do j = js_global-1, je_global+2
    x2(j) = ( xi2(j-1)*cosi(j-1)-xi2(j)*cosi(j) &
             +sini(j)-sini(j-1) ) &
          / (cosi(j-1)-cosi(j))
    if(j==js_global)x2(j-1)=-x2(j)
    dx2(j) = x2(j) - x2(j-1)
    idx2(j) = 1d0 / dx2(j)
   end do

! set phi direction
   mesh_sph3: select case (kmesh)
   case(0) mesh_sph3 ! uniform theta
    if    (ke_global==ks_global)then
     ketmp=4
    elseif(ke_global>=ks_global+3)then
     ketmp=ke_global
    else
     print *,"Error from ke",ke_global
    endif
    dxi3(ks_global-2:ke_global+2) = (xi3e-xi3s) / dble(ketmp)
   case(2) mesh_sph3 ! user specified mesh
    call other_kmesh(dxi3,ks_global,ke_global,xi3s,xi3e)
   end select mesh_sph3

   do k = ks_global-1,ke_global+2
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks_global-1)  = xi3s - dxi3(ks_global-1)*0.5d0
   x3(ks_global-2)  = x3(ks_global-1) - dx3(ks_global-1)
   xi3(ks_global-1) = xi3s
   xi3(ks_global-2) = xi3s-dxi3(ks_global-1)
   do k = ks_global,ke_global+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do
   xi3e = xi3(ke_global)

  case default coordinate_system
   print *,"Error from crdnt",crdnt
  end select coordinate_system

 else

! Read from previous gridfile if not first time step
  call readgrid('data/gridfile.bin')

  idxi1(gis-1:gie+2)=1d0/dxi1(gis-1:gie+2)
  idx1 (gis-1:gie+2)=1d0/dx1 (gis-1:gie+2)
  idxi2(gjs-1:gje+2)=1d0/dxi2(gjs-1:gje+2)
  idx2 (gjs-1:gje+2)=1d0/dx2 (gjs-1:gje+2)
  idxi3(gks-1:gke+2)=1d0/dxi3(gks-1:gke+2)
  idx3 (gks-1:gke+2)=1d0/dx3 (gks-1:gke+2)

  xi1e = xi1(ie)
  xi1s = xi1(is-1)

 end if

 return
end subroutine gridset

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OTHER_IMESH
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: User specified i meshes

subroutine other_imesh(dxi1,is,ie,xi1s,xi1e)

 real(8),intent(in):: xi1s,xi1e
 real(8),intent(inout),allocatable:: dxi1(:)
 integer,intent(in):: is,ie

!-----------------------------------------------------------------------------

 integer:: i
 real(8):: xrng, irng, xr, xrnew, xrmax, maxerr, fx, dfx, xmin
 real(8):: radstar

!-----------------------------------------------------------------------------

 xrmax = 1.015d0
 maxerr = 1d-10
 radstar = 5.6d13

 xr = 1.01d0
 xrng = radstar - xi1s ; irng = dble(400 - is + 1)
 xmin = 9d10

 if(xrng/irng<xmin)then
  print *,"Error from geometrical_series ;"
  print *,"xmin should be smaller or uniform mesh should be chosen",xmin
  stop
 end if

 do i = 1, 10000000
  fx = (xr-1d0)*xrng - xmin * (xr**irng-1d0)
  dfx = xrng - irng * xmin * xr**(irng-1d0)

  xrnew = xr - fx/dfx

  if(abs((xrnew-xr)/xr)<maxerr)then
   xr = xrnew ; exit
  end if
  if(xrnew<1d0)xrnew = 2d0

  xr = xrnew
 end do

 if(xr>xrmax)then
  print *,"xmin too small", xmin, xr
  stop
 end if

 dxi1(is) = xmin
 do i = is+1, 400
  dxi1(i) = dxi1(i-1) * xr
 end do
 dxi1(is-1) = dxi1(is) ; dxi1(is-2) = dxi1(is+1)
 dxi1(401:ie+2) = (xi1e-radstar)/dble(ie-400)

 if(xr-1d0<maxerr) dxi1 = (xi1e-xi1s) / irng


return
end subroutine other_imesh

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OTHER_JMESH
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: User specified j meshes

subroutine other_jmesh(dxi2,js,je,xi2s,xi2e)

 real(8),intent(in):: xi2s,xi2e
 real(8),intent(inout),allocatable:: dxi2(:)
 integer,intent(in):: js,je

 integer::j
 real(8),allocatable:: xi2(:)
 real(8):: dcos

!-----------------------------------------------------------------------------

!!$! for equal spacing in cos(theta)
!!$ allocate( xi2(js-2:je+2) )
!!$! dcos = (cos(xi2e)-cos(xi2s))/dble((je-js+1))
!!$ dcos = (cos(xi2e)-cos(xi2s))/dble(80)
!!$
!!$ xi2(js-1) = 1d0
!!$ do j = js, 40
!!$  xi2(j) = xi2(j-1) + dcos*0.125d0
!!$ end do
!!$ do j = 41, 60
!!$  xi2(j) = xi2(j-1) + dcos*0.25d0
!!$ end do
!!$ do j = 61, 80
!!$  xi2(j) = xi2(j-1) + dcos*0.5d0
!!$ end do
!!$ do j = 81, 120
!!$  xi2(j) = xi2(j-1) + dcos
!!$ end do
!!$ do j = 121, 140
!!$  xi2(j) = xi2(j-1) + dcos*0.5d0
!!$ end do
!!$ do j = 141, 160
!!$  xi2(j) = xi2(j-1) + dcos*0.25d0
!!$ end do
!!$ do j = 161, 200
!!$  xi2(j) = xi2(j-1) + dcos*0.125d0
!!$ end do
!!$ xi2(je) = -1d0
!!$
!!$ do j = js-1, je
!!$  xi2(j) = acos(xi2(j))
!!$ end do
!!$ do j = js, je
!!$  dxi2(j) = xi2(j) - xi2(j-1)
!!$ end do
!!$ dxi2(js-1) = dxi2(js)
!!$ dxi2(js-2) = dxi2(js+1)
!!$ dxi2(je+1) = dxi2(je)
!!$ dxi2(je+2) = dxi2(je-1)


 ! for equal spacing in cos(theta)
 allocate( xi2(js-2:je+2) )
 dcos = (cos(xi2e)-cos(xi2s))/dble(je-js+1)

 xi2(js-1) = cos(xi2s)
 do j = js, je
  xi2(j) = xi2(j-1) + dcos
 end do
 xi2(je) = cos(xi2e)
 do j = js-1, je
  xi2(j) = acos(xi2(j))
 end do
 do j = js, je
  dxi2(j) = xi2(j) - xi2(j-1)
 end do
 dxi2(js-1) = dxi2(js)
 dxi2(js-2) = dxi2(js+1)
 dxi2(je+1) = dxi2(je)
 dxi2(je+2) = dxi2(je-1)

return
end subroutine other_jmesh

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OTHER_KMESH
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: User specified j meshes

subroutine other_kmesh(dxi3,ks,ke,xi3s,xi3e)

 use grid,only:dxi1,xi1,is,ie

 real(8),intent(inout):: xi3s,xi3e
 real(8),intent(inout),allocatable:: dxi3(:)
 integer,intent(in):: ks,ke

!-----------------------------------------------------------------------------

 dxi3(ke/2+is:ke+2)  = dxi1(is:ie/2+2)
 dxi3(ke/2:ks-2:-1) = dxi1(is:ie/2+2)
! dxi3(ks:ke+2) = dxi1(is:ke+2)
! dxi3(ks-1) = dxi3(ks)
! dxi3(ks-2) = dxi3(ks+1)
! dxi3(ke+1) = dxi3(ke)
! dxi3(ke+2) = dxi3(ke-1)

! do i = is, ie
!  if(xi1(i)>xi3e)then
 xi3s = -xi1(ke/2)!-xi1(527)
 xi3e =  xi1(ke/2)
!  end if
! end do

return
end subroutine other_kmesh

end module gridset_mod
