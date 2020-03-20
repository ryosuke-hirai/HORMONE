!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRIDSET
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set grid

subroutine gridset

  use settings,only:courant,imesh,jmesh,kmesh,eq_sym,start
  use grid
  use constants,only:pi
  use gravmod
  use merger_mod,only:inifile2

  implicit none

  integer jetmp,ketmp,dummy
  real*8 x1min, x2min, x3min

  namelist /geoscon/ x1min, x2min, x3min

!-------------------------------------------------------------------------
gie=1200!temp
 if(gravswitch/=0)then
  deallocate(x1, xi1, dx1, dxi1, idx1, idxi1, &
             x2, xi2, dx2, dxi2, idx2, idxi2, &
             x3, xi3, dx3, dxi3, idx3, idxi3 )
  allocate( &
   x1(gis-2:gie+2), xi1(gis-2:gie+2), dx1(gis-2:gie+2), dxi1(gis-2:gie+2), &
   idx1(gis-2:gie+2), idxi1(gis-2:gie+2), &
   x2(gjs-2:gje+2), xi2(gjs-2:gje+2), dx2(gjs-2:gje+2), dxi2(gjs-2:gje+2), &
   idx2(gjs-2:gje+2), idxi2(gjs-2:gje+2), &
   x3(gks-2:gke+2), xi3(gks-2:gke+2), dx3(gks-2:gke+2), dxi3(gks-2:gke+2), &
   idx3(gks-2:gke+2), idxi3(gks-2:gke+2) &
  )
 end if
gie=900!temp
 if(start==0)then
! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==0)then

   ! set x direction
   if(imesh==0)then ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie-is+1)
   elseif(imesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi1,x1min,is,ie,xi1s,xi1e)
   elseif(imesh==2)then ! user specified mesh
    call other_imesh(dxi1,is,ie,xi1s,xi1e)
   end if

   do i = is-1,ie+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is-1)  = xi1s - dxi1(is-1)*0.5d0
   x1(is-2)  = x1(is-1) - dx1(is-1)
   xi1(is-1) = xi1s
   xi1(is-2) = xi1s - dxi1(is-1)
   do i = is,ie+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do

   if(gravswitch/=0)then
    dx1 = dxi1 ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
    do i = is-1, gis-2, -1
     xi1(i) = xi1(i+1) - dxi1(i)
     x1 (i) = x1 (i+1) - dx1 (i)
    end do
    do i = ie+1, gie+2
     xi1(i) = xi1(i-1) + dxi1(i)
     x1 (i) = x1 (i-1) + dx1 (i)
    end do
   end if

   ! set y direction
   if(jmesh==0)then ! uniform mesh
    dxi2 = (xi2e-xi2s) / dble(je-js+1)
   elseif(jmesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi2,x2min,js,je,xi2s,xi2e)
   elseif(jmesh==2)then ! user specified mesh
    call other_jmesh(dxi2,js,je,xi2s,xi2e)
   end if

   do j = js-1,je+2
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js-1)  = xi2s - dxi2(js-1)*0.5d0
   x2(js-2)  = x2(js-1) - dx2(js-1)
   xi2(js-1) = xi2s
   xi2(js-2) = xi2s - dxi2(js-1)
   do j = js,je+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do

   if(gravswitch/=0)then
    dx2 = dxi2 ; idx2 = 1d0 / dx2 ; idxi2 = 1d0 / dxi2
    do j = js-1, gjs-2, -1
     xi2(j) = xi2(j+1) - dxi2(j)
     x2 (j) = x2 (j+1) - dx2 (j)
    end do
    do j = je+1, gje+2
     xi2(j) = xi2(j-1) + dxi2(j)
     x2 (j) = x2 (j-1) + dx2 (j)
    end do
   end if

   ! set z direction
   if(kmesh==0)then ! uniform mesh
    dxi3 = (xi3e-xi3s) / dble(ke-ks+1)
   elseif(kmesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi3,x3min,ks,ke,xi3s,xi3e)
   elseif(kmesh==3)then ! user specified mesh
    call other_kmesh(dxi3,ks,ke,xi3s,xi3e)
   end if

   do k = ks-1,ke+2
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks-1)  = xi3s - dxi3(ks-1)*0.5d0
   x3(ks-2)  = x3(ks-1) - dx3(ks-1)
   xi3(ks-1) = xi3s
   xi3(ks-2) = xi3s - dxi3(ks-1)
   do k = ks,ke+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do

   if(gravswitch/=0)then
    dx3 = dxi3 ; idx3 = 1d0 / dx3 ; idxi3 = 1d0 / dxi3
    do k = ks-1, gks-2, -1
     xi3(k) = xi3(k+1) - dxi3(k)
     x3 (k) = x3 (k+1) - dx3 (k)
    end do
    do k = ke+1, gke+2
     xi3(k) = xi3(k-1) + dxi3(k)
     x3 (k) = x3 (k-1) + dx3 (k)
    end do
   end if

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  else if(crdnt==1)then
   xi2s = 0.d0 ; xi2e = 2.d0*pi
   ! set r direction
   if(imesh==0)then ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie-is+1)
   elseif(imesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi1,x1min,is,ie,xi1s,xi1e)
   elseif(imesh==2)then ! user specified mesh
    call other_imesh(dxi1,is,ie,xi1s,xi1e)
   end if

   do i = is-1,ie+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is-1)  = xi1s - dxi1(is-1)*0.5d0
   x1(is-2)  = x1(is-1) - dx1(is-1)
   xi1(is-1) = xi1s
   xi1(is-2) = xi1s - dxi1(is-1)
   do i = is,ie+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do


! temp: for volumetric centre
!!$   do i = is-1, ie+2
!!$    x1(i)  = sqrt( (xi1(i-1)*xi1(i-1)+xi1(i)*xi1(i)) *5d-1 )
!!$    if(i==is-1)x1(i)=-x1(i)
!!$    dx1(i) = x1(i) - x1(i-1)
!!$    idx1(i) = 1d0 / dx1(i)
!!$   end do

   if(gravswitch/=0)then
    do i = ie+1, gie+2
     dxi1(i) = dxi1(i-1) * dxi1(ie)/dxi1(ie-1)
     xi1(i) = xi1(i-1) + dxi1(i)
     x1(i)  = 0.5d0*(xi1(i)+xi1(i-1))!sqrt( (xi1(i-1)*xi1(i-1)+xi1(i)*xi1(i)) *5d-1 )
     dx1(i) = x1(i) - x1(i-1)
     idx1(i) = 1d0 / dx1(i)
    end do
    dx1(gis-2) = dxi1(gis-2) ; idx1 = 1d0 / dx1 ; idxi1 = 1d0 / dxi1
   end if

   ! set theta direction
   if    (je==1)then;jetmp=4
   elseif(je>=4)then;jetmp=je
   else ;print *,"Error from je",je;endif
   do j = js-1,je+2
    dxi2(j) = 2.d0*pi / dble(jetmp)
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js-1)  = xi2s - dxi2(js-1)*0.5d0
   x2(js-2)  = x2(js-1) - dx2(js-1)
   xi2(js-1) = xi2s
   xi2(js-2) = xi2s - dxi2(js-1)
   do j = js,je+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do
   xi2e = xi2(je)

   ! set z direction
   if(kmesh==0)then ! uniform mesh
    dxi3 = (xi3e-xi3s) / dble(ke-ks+1)
   elseif(kmesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi3,x3min,ks,ke,xi3s,xi3e)
   elseif(kmesh==2)then ! user specified mesh
    call other_kmesh(dxi3,ks,ke,xi3s,xi3e)
   end if

   do k = ks-1,ke+2
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks-1)  = xi3s - dxi3(ks-1)*0.5d0
   x3(ks-2)  = x3(ks-1) - dx3(ks-1)
   xi3(ks-1) = xi3s
   xi3(ks-2) = xi3s - dxi3(ks-1)
   do k = ks,ke+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do

   if(gravswitch/=0)then
!    dx3 = dxi3 ; idx3 = 1d0 / dx3 ; idxi3 = 1d0 / dxi3
    do k = ks-1, gks-2, -1
     dxi3(k) = dxi3(k+1) * dxi3(ks)/dxi3(ks+1)
     idxi3(k)= 1d0 / dxi3(k)
    end do
    do k = ks-1, gks-1, -1
     dx3(k) = 0.5d0*(dxi3(k-1) + dxi3(k))
     idx3(k) = 1d0/dx3(k)
     xi3(k) = xi3(k+1) - dxi3(k+1)
     x3 (k) = x3 (k+1) - dx3 (k+1)
    end do
    do k = ke+1, gke+2
     dxi3(k) = dxi3(k-1) * dxi3(ke)/dxi3(ke-1)
     idxi3(k)= 1d0 / dxi3(k)
     dx3(k) = 0.5d0*(dxi3(k-1) + dxi3(k))
     idx3(k)= 1d0 / dx3(k)
     xi3(k) = xi3(k-1) + dxi3(k)
     x3 (k) = x3 (k-1) + dx3 (k)
    end do
    if(eq_sym)then
     dxi3(gks-2:gks-1) = dxi3(gks+1:gks:-1)
     idxi3(gks-2:gks-1) = 1d0/dxi3(gks-2:gks-1)
     dx3(gks) = 0.5d0*(dxi3(gks-1)+dxi3(gks))
     dx3(gks-1) = 0.5d0*(dxi3(gks-2)+dxi3(gks-1))
     idx3(gis:gks-1) = 1d0/dx3(gks:gks-1)
     xi3(gks-2) = xi3(gks-1)-dxi3(gks-1)
     x3(gks-1) = x3(gks) - dx3(gks)
     x3(gks-2) = x3(gks-1) - dx3(gks-1)
    end if
   end if



! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  else if(crdnt==2)then
   xi2s = 0.d0 ; xi2e = pi ; if(je==1.or.eq_sym) xi2e = pi * 0.5d0
   xi3s = 0.d0 ; xi3e = 2.d0*pi
   ! set r direction
   if(imesh==0)then ! uniform mesh
    dxi1 = (xi1e-xi1s) / dble(ie-is+1)
   elseif(imesh==1)then ! geometrical series
    open(unit=1,file='parameters',status='old')
    read(1,NML=geoscon)
    close(1)
    call geometrical_series(dxi1,x1min,is,ie,xi1s,xi1e)
   elseif(imesh==2)then ! user specified mesh
    call other_imesh(dxi1,is,ie,xi1s,xi1e)
   end if

   do i = is-1,ie+2
    dx1(i)  = 0.5d0 * ( dxi1(i-1) + dxi1(i) )
    idxi1(i)= 1.d0 / dxi1(i)
    idx1(i) = 1.d0 / dx1(i)
   end do

   x1(is-1)  = xi1s - dxi1(is-1)*0.5d0
   x1(is-2)  = x1(is-1) - dx1(is-1)
   xi1(is-1) = xi1s
   xi1(is-2) = xi1s - dxi1(is-1)
   do i = is,ie+2
    x1(i)  = x1(i-1)  + dx1(i)
    xi1(i) = xi1(i-1) + dxi1(i)
   end do

!temp
! for volumetric centre
   do i = is-1, is+2
    x1(i) = 0.75d0*(xi1(i)+xi1(i-1))*(xi1(i)*xi1(i)+xi1(i-1)*xi1(i-1)) &
                  /(xi1(i)*xi1(i)+xi1(i)*xi1(i-1)*xi1(i-1)*xi1(i-1))
    if(i==is)x1(i-1)=-x1(i)
    dx1(i) = x1(i) - x1(i-1)
    idx1(i) = 1d0 / dx1(i)
   end do

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
   if(jmesh==0)then
    if    (je==1)then;jetmp=1
    elseif(je>=2)then;jetmp=je;endif
    dxi2(js-2:je+2) = xi2e / dble(jetmp)
   elseif(jmesh==2)then
    call other_jmesh(dxi2,js,je,xi2s,xi2e)
   end if

   do j = js-1,je+2
    dx2(j)  = 0.5d0 * ( dxi2(j-1) + dxi2(j) )
    idxi2(j)= 1.d0 / dxi2(j)
    idx2(j) = 1.d0 / dx2(j)
   end do

   x2(js-1)  = xi2s - dxi2(js-1)*0.5d0
   x2(js-2)  = x2(js-1) - dx2(js-1)
   xi2(js-1) = xi2s
   xi2(js-2) = xi2s-dxi2(js-1)
   do j = js,je+2
    x2(j)  = x2(j-1)  + dx2(j)
    xi2(j) = xi2(j-1) + dxi2(j)
   end do
   xi2e = xi2(je)

!temp
! for volumetric centre
   do j = js-1, je+2
    x2(j) = ( xi2(j-1)*cos(xi2(j-1))-xi2(j)*cos(xi2(j)) &
             +sin(xi2(j))-sin(xi2(j-1)) ) &
          / (cos(xi2(j-1))-cos(xi2(j)))
    if(j==js)x2(j-1)=-x2(j)
    dx2(j) = x2(j) - x2(j-1)
    idx2(j) = 1d0 / dx2(j)
   end do

   ! set phi direction
   if    (ke==1)then;ketmp=4
   elseif(ke>=4)then;ketmp=je
   else ;print *,"Error from ke",ke;endif
   do k = ks-1,ke+2
    dxi3(k) = xi3e / dble(ketmp)
    dx3(k)  = 0.5d0 * ( dxi3(k-1) + dxi3(k) )
    idxi3(k)= 1.d0 / dxi3(k)
    idx3(k) = 1.d0 / dx3(k)
   end do

   x3(ks-1)  = - dxi3(ks-1)*0.5d0
   x3(ks-2)  = x3(ks-1) - dx3(ks-1)
   xi3(ks-1) = 0.d0
   xi3(ks-2) = -dxi3(ks-1)
   do k = ks,ke+2
    x3(k)  = x3(k-1)  + dx3(k)
    xi3(k) = xi3(k-1) + dxi3(k)
   end do
   xi3e = xi3(ke)

!hirai
!!$    open(55,file='gridhirai.bin',form='unformatted')
!!$    do i = is-1,ie
!!$     read(55)dummy,xi1(i),dxi1(i+1)
!!$    end do
!!$    close(55)
!!$    dxi1(is-1) = dxi1(is) ; dxi1(is-2) = dxi1(is+1)
!!$    xi1(is-2) = xi1(is-1) - dxi1(is-1)
!!$    dxi1(ie+1) = dxi1(ie) ; dxi1(ie+2) = dxi1(ie)
!!$    dx1(is-2) = dxi1(is-2) ; x1(is-2) = xi1(is-2) - 0.5d0*dxi1(is-2)
!!$    do i = is-1,ie+2
!!$     idxi1(i)  = 1.0d0/dxi1(i)
!!$     x1(i)     = xi1(i) - 0.5d0*dxi1(i)
!!$     dx1(i)    = x1(i) - x1(i-1)
!!$     idx1(i)   = 1.0d0/dx1(i)
!!$    end do

  else
   print *,"Error from crdnt",crdnt
  end if

 elseif(start>int(inifile2))then

  gie=1200

  open(unit=41,file='data/gridfile2.bin',status='old',form='unformatted')
  read(41)x1(gis-2:gie+2),xi1(gis-2:gie+2),dxi1(gis-2:gie+2),dx1(gis-2:gie+2), &
          x2(gjs-2:gje+2),xi2(gjs-2:gje+2),dxi2(gjs-2:gje+2),dx2(gjs-2:gje+2), &
          x3(gks-2:gke+2),xi3(gks-2:gke+2),dxi3(gks-2:gke+2),dx3(gks-2:gke+2)
  close(41)

  idxi1(gis-1:gie+2)=1d0/dxi1(gis-1:gie+2)
  idx1 (gis-1:gie+2)=1d0/dx1 (gis-1:gie+2)
  idxi2(gjs-1:gje+2)=1d0/dxi2(gjs-1:gje+2)
  idx2 (gjs-1:gje+2)=1d0/dx2 (gjs-1:gje+2)
  idxi3(gks-1:gke+2)=1d0/dxi3(gks-1:gke+2)
  idx3 (gks-1:gke+2)=1d0/dx3 (gks-1:gke+2)

 elseif(start>0)then
  gie=900;ie=gie
  open(unit=41,file='data/gridfile_old.bin',status='old',form='unformatted')
  read(41)x1(gis-2:gie+2),xi1(gis-2:gie+2),dxi1(gis-2:gie+2),dx1(gis-2:gie+2), &
          x2(gjs-2:gje+2),xi2(gjs-2:gje+2),dxi2(gjs-2:gje+2),dx2(gjs-2:gje+2), &
          x3(gks-2:gke+2),xi3(gks-2:gke+2),dxi3(gks-2:gke+2),dx3(gks-2:gke+2)
  close(41)

  idxi1(gis-1:gie+2)=1d0/dxi1(gis-1:gie+2)
  idx1 (gis-1:gie+2)=1d0/dx1 (gis-1:gie+2)
  idxi2(gjs-1:gje+2)=1d0/dxi2(gjs-1:gje+2)
  idx2 (gjs-1:gje+2)=1d0/dx2 (gjs-1:gje+2)
  idxi3(gks-1:gke+2)=1d0/dxi3(gks-1:gke+2)
  idx3 (gks-1:gke+2)=1d0/dx3 (gks-1:gke+2)

 end if

return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GEOMETRICAL_SERIES
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate dxi's in a geometrical series.

subroutine geometrical_series(dxi,xmin,is,ie,xis,xie)

 implicit none

 integer,intent(in):: is,ie
 real*8,intent(in):: xis,xie,xmin
 real*8,intent(inout),allocatable:: dxi(:)
 integer i
 real*8 xrng, irng, xr, xrnew, xrmax, err, maxerr, fx, dfx

!-----------------------------------------------------------------------------

 xrmax = 1.015d0
 maxerr = 1d-10

 xr = 1.01d0
 xrng = xie - xis ; irng = dble(ie - is + 1)

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

 dxi(is) = xmin
 do i = is+1, ie
  dxi(i) = dxi(i-1) * xr
 end do
 dxi(is-1) = dxi(is) ; dxi(is-2) = dxi(is+1)
 dxi(ie+1) = dxi(ie)*xr ; dxi(ie+2) = dxi(ie)*xr*xr

 if(xr-1d0<maxerr) dxi = (xie-xis) / irng

return
end subroutine geometrical_series

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OTHER_IMESH
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: User specified i meshes

subroutine other_imesh(dxi1,is,ie,xi1s,xi1e)

 implicit none

 real*8,intent(in):: xi1s,xi1e
 real*8,intent(inout),allocatable:: dxi1(:)
 integer,intent(in):: is,ie

!-----------------------------------------------------------------------------

 

return
end subroutine other_imesh

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE OTHER_JMESH
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: User specified j meshes

subroutine other_jmesh(dxi2,js,je,xi2s,xi2e)

 implicit none

 real*8,intent(in):: xi2s,xi2e
 real*8,intent(inout),allocatable:: dxi2(:)
 integer,intent(in):: js,je

 integer j
 real*8,allocatable:: xi2(:)
 real*8 dcos, cosnow

!-----------------------------------------------------------------------------

! for equal spacing in cos(theta)
 allocate( xi2(js-2:je+2) )
 dcos = (cos(xi2e)-cos(xi2s))/dble(je-js+1)

 xi2(js-1) = 1d0
 do j = js, je
  xi2(j) = xi2(j-1) + dcos
 end do
 xi2(je) = 0d0
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

 use grid,only:dxi1,xi1

 implicit none

 real*8,intent(inout):: xi3s,xi3e
 real*8,intent(inout),allocatable:: dxi3(:)
 integer,intent(in):: ks,ke
 integer i

!-----------------------------------------------------------------------------

 dxi3(ke/2+is:ke+2)  = dxi1(is:ke/2+2)
 dxi3(ke/2:ks-2:-1) = dxi1(is:ke/2+2)
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


end subroutine gridset


