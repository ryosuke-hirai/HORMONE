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
 use grid
 use physval
 use gravmod
 use dirichlet_mod

!-----------------------------------------------------------------------------

 lmax = (gie-gis+1)*(gje-gjs+1)*(gke-gks+1)

! 1 dimensional arrays
! grid-related variables
 allocate(x1(gis-2:gie+2)); x1=0d0
 allocate(xi1,dx1,dxi1,idx1,idxi1,source=x1)

 allocate(x2(gjs-2:gje+2)); x2=0d0
 allocate(xi2,dx2,dxi2,idx2,idxi2,source=x2)

 allocate(x3(gks-2:gke+2)); x3=0d0
 allocate(xi3,dx3,dxi3,idx3,idxi3,source=x3)

!  metric-related variables
 allocate(detg1(is-2:ie+2)); detg1=0d0
 allocate(idetg1,sx1,g22,source=detg1)

 allocate(scot(js-2:je+2)); scot=0d0
 allocate(sisin,source=scot)

 allocate(detg2(is-2:ie+2,js-2:je+2)); detg2=0d0
 allocate(idetg2,g33,source=detg2)

 allocate(dvol(is-2:ie+2,js-2:je+2,ks-2:ke+2)); dvol=0d0
 allocate(idetg3,sa1,sa2,sa3,source=dvol)

! 3 dimensional arrays
! physical variables
!  Strictly non-zero quantities
 allocate(d(is-2:ie+2,js-2:je+2,ks-2:ke+2)); d = 1d0
 allocate(p,e,T,ptot,cs,eint,imu,source=d)

!  Initially zero quantities
 allocate(phi(is-2:ie+2,js-2:je+2,ks-2:ke+2)); phi = 0d0
 allocate(v1,v2,v3,b1,b2,b3,grv1,grv2,grv3,source=phi)
 allocate(shock(is-2:ie+2,js-2:je+2,ks-2:ke+2)); shock = 0

! 4 dimensional  arrays
! gradients
 allocate(dd(is-2:ie+2,js-2:je+2,ks-2:ke+2,1:3)); dd = 0d0
 allocate(de,dphi,dm1,dm2,dm3,db1,db2,db3,dmu,source=dd)

! conserved quantities and flux
 allocate(u(is-2:ie+2,js-2:je+2,ks-2:ke+2,1:9)); u = 0d0
 allocate(flux1,flux2,flux3,source=u)
 allocate(src(is:ie,js:je,ks:ke,1:9)); src = 0d0
 allocate(uorg,source=src)

! gravity-related variables
 if(gravswitch>=1)then
  allocate(grvphi(gis-2:gie+2,gjs-2:gje+2,gks-2:gke+2)); grvphi=0d0
  allocate(grvphiold,source=grvphi)
  allocate(hgsrc(gis:gie,gjs:gje,gks:gke))
!  for MICCG method
  allocate(modlimax(1:lmax)); modlimax=0d0
  allocate(a1(0:lmax)); a1=0d0
  allocate(a2,a3,preca,precb,precc,precd,prece,source=a1)
!  for gravbound
  allocate(phiio(gie+1:gie+2,gjs-2:gje+2), phiii(gis-2:gis-1,gjs-2:gje+2), &
           phi1o(gie+1:gie+2,gks-2:gke+2), phi3i(gis-2:gie+2,gks-2:gks-1), &
           phi3o(gis-2:gie+2,gke+1:gke+2), mc(is-1:ie+2) )
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
 end if

! allocate external gravitational field if necessary
 if(include_extgrv)allocate(extgrv,source=grvphi)
 
 T = 1d3  ! initial guess for temperature
 
 return
end subroutine allocations

end module allocation_mod
