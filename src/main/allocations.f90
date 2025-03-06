module allocation_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE ALLOCATIONS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate various variables.

subroutine allocations

 use settings
 use derived_types,only:null_sink
 use grid
 use physval
 use gravmod
 use sink_mod
 use dirichlet_mod

 integer::i,j,k,n

!-----------------------------------------------------------------------------

! 1 dimensional arrays =======================================================
! grid-related variables
 allocate(x1(gis_global-2:gie_global+2))
 allocate(xi1,dx1,dxi1,idx1,idxi1,mold=x1)

 allocate(x2(gjs_global-2:gje_global+2))
 allocate(xi2,dx2,dxi2,idx2,idxi2,mold=x2)

 allocate(x3(gks_global-2:gke_global+2))
 allocate(xi3,dx3,dxi3,idx3,idxi3,mold=x3)

!  metric-related variables
 allocate(detg1(is_global-2:ie_global+2))
 allocate(idetg1,sx1,g22,mold=detg1)

 allocate(scot(js_global-2:je_global+2))
 allocate(sisin,mold=scot)

 allocate(detg2(is_global-2:ie_global+2,js_global-2:je_global+2))
 allocate(idetg2,g33,mold=detg2)

 allocate(dvol(is_global-2:ie_global+2,js_global-2:je_global+2,ks_global-2:ke_global+2))
 allocate(idetg3,sa1,sa2,sa3,Imom,mold=dvol)
 allocate(car_x(1:3,is-1:ie+1,js-1:je+1,ks-1:ke+1))

! 3 dimensional arrays =======================================================
! physical variables
!  Strictly non-zero quantities
 allocate(d(is-2:ie+2,js-2:je+2,ks-2:ke+2))
 allocate(p,e,T,ptot,cs,eint,erad,imu,mold=d)

!  Initially zero quantities
 allocate(phi(is-2:ie+2,js-2:je+2,ks-2:ke+2))
 allocate(v1,v2,v3,b1,b2,b3,grv1,grv2,grv3,mold=phi)
 allocate(shock(is-2:ie+2,js-2:je+2,ks-2:ke+2))

! Parallel first touch for OpenMP optimization on NUMA cores
!$omp parallel
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js, je
   do i = is, ie
! Strictly non-zero quantities
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
! Initially zero quantities
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
    grv1(i,j,k) = 0d0; grv2(i,j,k) = 0d0; grv3(i,j,k) = 0d0
    phi(i,j,k) = 0d0; erad(i,j,k) = 0d0
    shock(i,j,k) = 0
   end do
  end do
 end do
!$omp end do

! Boundary values
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js, je
   do i = is-2, is-1
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js, je
   do i = ie+1, ie+2
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js-2, js-1
   do i = is, ie
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = je+1, je+2
   do i = is, ie
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ks-2, ks-1
  do j = js, je
   do i = is, ie
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp do private(i,j,k) collapse(3) schedule(static)
 do k = ke+1, ke+2
  do j = js, je
   do i = is, ie
    d(i,j,k) = 1d0; p(i,j,k) = 1d0; e(i,j,k) = 1d0
    ptot(i,j,k) = 1d0; cs(i,j,k) = 1d0; eint(i,j,k) = 1d0
    T(i,j,k) = 1d3; imu(i,j,k) = 1d0/muconst
    v1(i,j,k) = 0d0; v2(i,j,k) = 0d0; v3(i,j,k) = 0d0
    b1(i,j,k) = 0d0; b2(i,j,k) = 0d0; b3(i,j,k) = 0d0
   end do
  end do
 end do
!$omp end do
!$omp end parallel

! 4 dimensional arrays =======================================================
! gradients
 allocate(dd(is-2:ie+2,js-2:je+2,ks-2:ke+2,1:3))
 allocate(de,dphi,dm1,dm2,dm3,db1,db2,db3,dmu,mold=dd)
 if(radswitch>0)allocate(der,mold=dd)
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js, je
   do i = is-1, ie+1
    dd(i,j,k,1) = 0d0 ; dphi(i,j,k,1) = 0d0
    dm1(i,j,k,1) = 0d0 ; dm2(i,j,k,1) = 0d0 ; dm3(i,j,k,1) = 0d0
    db1(i,j,k,1) = 0d0 ; db2(i,j,k,1) = 0d0 ; db3(i,j,k,1) = 0d0
    dmu(i,j,k,1) = 0d0
    if(radswitch>0)der(i,j,k,1) = 0d0
   end do
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js-1, je+1
   do i = is, ie
    dd(i,j,k,2) = 0d0 ; dphi(i,j,k,2) = 0d0
    dm1(i,j,k,2) = 0d0 ; dm2(i,j,k,2) = 0d0 ; dm3(i,j,k,2) = 0d0
    db1(i,j,k,2) = 0d0 ; db2(i,j,k,2) = 0d0 ; db3(i,j,k,2) = 0d0
    dmu(i,j,k,2) = 0d0
    if(radswitch>0)der(i,j,k,2) = 0d0
   end do
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks-1, ke+1
  do j = js, je
   do i = is, ie
    dd(i,j,k,3) = 0d0 ; dphi(i,j,k,3) = 0d0
    dm1(i,j,k,3) = 0d0 ; dm2(i,j,k,3) = 0d0 ; dm3(i,j,k,3) = 0d0
    db1(i,j,k,3) = 0d0 ; db2(i,j,k,3) = 0d0 ; db3(i,j,k,3) = 0d0
    dmu(i,j,k,3) = 0d0
    if(radswitch>0)der(i,j,k,3) = 0d0
   end do
  end do
 end do
!$omp end parallel do

! conserved quantities and flux
 allocate(u(is-2:ie+2,js-2:je+2,ks-2:ke+2,1:ufnmax))
 allocate(flux1,flux2,flux3,mold=u)
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js, je
   do i = is-1, ie
    flux1(i,j,k,1:ufnmax) = 0d0
   end do
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks, ke
  do j = js-1, je
   do i = is, ie
    flux2(i,j,k,1:ufnmax) = 0d0
   end do
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
 do k = ks-1, ke
  do j = js, je
   do i = is, ie
    flux3(i,j,k,1:ufnmax) = 0d0
   end do
  end do
 end do
!$omp end parallel do
 allocate(src(is:ie,js:je,ks:ke,1:ufnmax))
 allocate(uorg,mold=src)
!$omp parallel do private(i,j,k,n) collapse(4) schedule(static)
 do n = 1, ufnmax
  do k = ks, ke
   do j = js, je
    do i = is, ie
     uorg(i,j,k,n) = 0d0
     u(i,j,k,n) = 0d0
     src(i,j,k,n) = 0d0
    end do
   end do
  end do
 end do
!$omp end parallel do

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! gravity-related variables
 if(gravswitch>=1)then
  allocate(grvphi(gis-2:gie+2,gjs-2:gje+2,gks-2:gke+2))
  allocate(totphi,mold=grvphi)
  allocate(hgsrc(gis:gie,gjs:gje,gks:gke))
  allocate(grvpsi,gsrc,lapphi,mold=hgsrc)
  allocate(grvphiorg(gis:gie,gjs:gje,gks:gke,1:2))
! Parallel first touch for OpenMP optimization on NUMA cores
!$omp parallel do private(i,j,k) collapse(3) schedule(static)
  do k = gks, gke
   do j = gjs, gje
    do i = gis, gie
     grvphi(i,j,k) = 0d0
     grvpsi(i,j,k) = 0d0
     grvphiorg(i,j,k,1:2) = 0d0
     hgsrc(i,j,k) = 1d0
    end do
   end do
  end do
!$omp end parallel do

!  for gravbound
  allocate(phiio(gie+1:gie+2,gjs-2:gje+2), phiii(gis-2:gis-1,gjs-2:gje+2), &
           phi1o(gie+1:gie+2,gks-2:gke+2), phi3i(gis-2:gie+2,gks-2:gks-1), &
           phi3o(gis-2:gie+2,gke+1:gke+2), mc(is_global-1:ie_global+2) )
 end if

! allocate Dirichlet variables if Dirichlet boundary is applied
 if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or. &
    bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9.or. &
    is_test)then
  allocate(d0(is-2:ie+2,js-2:je+2,ks-2:ke+2)); d0=0d0
  allocate(p0,v10,v20,v30,b10,b20,b30,source=d0)
 end if

! allocate chemical composition if compswitch/=0
 if(compswitch>0)then
  allocate( &
   spc(1:spn,is-2:ie+2,js-2:je+2,ks-2:ke+2), spcorg(1:spn,is:ie,js:je,ks:ke), &
   dspc  (1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
   spcflx(1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
   species(1:spn) )
  if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or. &
     bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9)then
   allocate( spc0,source=spc )
  end if

  species(:) = 'aaa'
!$omp parallel do private(i,j,k,n) collapse(4) schedule(static)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     do n = 1, spn
      spc(n,i,j,k) = 0d0
     end do
    end do
   end do
  end do
!$omp end parallel do
 end if

! allocate external gravitational field if necessary
 if(include_extgrv)allocate(extgrv,mold=grvphi)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! allocate sink particle-related variables if necessary
 if(include_sinks)then
  allocate(snkphi,mold=d)
  allocate(sink(1:nsink))
  do n = 1, nsink
   call null_sink(sink(n))
  end do
 end if

 return
end subroutine allocations

end module allocation_mod
