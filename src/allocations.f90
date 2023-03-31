!\ \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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

 implicit none

!-----------------------------------------------------------------------------

 in = ie + 2 ; jn = je + 2 ; kn = ke + 2
 gin = gie + 2 ; gjn = gje + 2 ; gkn = gke + 2
 lmax = (gie-gis+1)*(gje-gjs+1)*(gke-gks+1)

 allocate( &
  x1(-1:gin),xi1(-1:gin),dx1(-1:gin),dxi1(-1:gin),idx1(-1:gin),idxi1(-1:gin), &
  x2(-1:gjn),xi2(-1:gjn),dx2(-1:gjn),dxi2(-1:gjn),idx2(-1:gjn),idxi2(-1:gjn), &
  x3(-1:gkn),xi3(-1:gkn),dx3(-1:gkn),dxi3(-1:gkn),idx3(-1:gkn),idxi3(-1:gkn), &

  d (-1:in,-1:jn,-1:kn), p (-1:in,-1:jn,-1:kn), e (-1:in,-1:jn,-1:kn), &
  v1(-1:in,-1:jn,-1:kn), v2(-1:in,-1:jn,-1:kn), v3(-1:in,-1:jn,-1:kn), &
  b1(-1:in,-1:jn,-1:kn), b2(-1:in,-1:jn,-1:kn), b3(-1:in,-1:jn,-1:kn), &
  ptot(-1:in,-1:jn,-1:kn), cs(-1:in,-1:jn,-1:kn), phi(-1:in,-1:jn,-1:kn), &
  T(-1:in,-1:jn,-1:kn), eint(-1:in,-1:jn,-1:kn), imu (-1:in,-1:jn,-1:kn), &
  shock(-1:in,-1:jn,-1:kn),&
  
  dd (-1:in,-1:jn,-1:kn,1:3), de(-1:in,-1:jn,-1:kn,1:3), dphi(-1:in,-1:jn,-1:kn,1:3),&
  dm1(-1:in,-1:jn,-1:kn,1:3), dm2(-1:in,-1:jn,-1:kn,1:3), dm3(-1:in,-1:jn,-1:kn,1:3),&
  db1(-1:in,-1:jn,-1:kn,1:3), db2(-1:in,-1:jn,-1:kn,1:3), db3(-1:in,-1:jn,-1:kn,1:3),&
  dmu(-1:in,-1:jn,-1:kn,1:3), &
! temporary
  dw1(-1:in,-1:jn,-1:kn,1:3), dw2(-1:in,-1:jn,-1:kn,1:3), dw3(-1:in,-1:jn,-1:kn,1:3),&
  dw4(-1:in,-1:jn,-1:kn,1:3), dw5(-1:in,-1:jn,-1:kn,1:3), dw6(-1:in,-1:jn,-1:kn,1:3),&
  dw7(-1:in,-1:jn,-1:kn,1:3), dp(-1:in,-1:jn,-1:kn,1:3),&
! end temporary
  u(-1:in,-1:jn,-1:kn,1:9), uorg(is:ie,js:je,ks:ke,1:9), &
  flux1(-1:in,-1:jn,-1:kn,1:9), &
  flux2(-1:in,-1:jn,-1:kn,1:9), &
  flux3(-1:in,-1:jn,-1:kn,1:9), &
  src  (is:ie,js:je,ks:ke,1:9), &
  grv1(-1:in,-1:jn,-1:kn), grv2(-1:in,-1:jn,-1:kn),grv3(-1:in,-1:jn,-1:kn), &

  detg1(-1:in), idetg1(-1:in), sx1(-1:in), g22(-1:in), &
  scot(-1:jn), sisin(-1:jn), &
  detg2(-1:in,-1:jn), idetg2(-1:in,-1:jn), g33(-1:in,-1:jn), &
  idetg3(-1:in,-1:jn,-1:kn), dvol(-1:in,-1:jn,-1:kn), &
  sa1(-1:in,-1:jn,-1:kn), sa2(-1:in,-1:jn,-1:kn), sa3(-1:in,-1:jn,-1:kn)  &
 )

! allocate gravity related quantities if gravswitch>=1
 if(gravswitch>=1)then
  allocate( &
   grvphi   (gis-2:gin,gjs-2:gjn,gks-2:gkn), &
   grvphiold(gis-2:gin,gjs-2:gjn,gks-2:gkn), &
   hgsrc(gis:gie,gjs:gje,gks:gke), &

   modlimax(1:lmax), &
   a1(0:lmax), a2(0:lmax), a3(0:lmax), &
   preca(0:lmax), precb(0:lmax), precc(0:lmax), precd(0:lmax), prece(0:lmax), &
 
   phiio(gie+1:gin,gjs-2:gjn), phiii(gis-2:gis-1,gjs-2:gjn), &
   phi1o(gie+1:gin,gks-2:gkn), phi3i(gis-2:gin,gks-2:gks-1), &
   phi3o(gis-2:gin,gke+1:gkn), &
   mc(is-1:in) &
  )
 end if

! allocate Dirichlet variables if Dirichlet boundary is applied
 if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or. &
    bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9)then
  allocate( &
   d0 (-1:in,-1:jn,-1:kn), p0 (-1:in,-1:jn,-1:kn), &
   b10(-1:in,-1:jn,-1:kn), b20(-1:in,-1:jn,-1:kn), b30(-1:in,-1:jn,-1:kn), &
   v10(-1:in,-1:jn,-1:kn), v20(-1:in,-1:jn,-1:kn), v30(-1:in,-1:jn,-1:kn)  &
  )
 end if

! allocate chemical composition if compswitch/=0
 if(compswitch>0)then
  allocate( &
   spc(1:spn,-1:in,-1:jn,-1:kn), spcorg(1:spn,is:ie,js:je,ks:ke), &
   dspc  (1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
   spcflx(1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
   species(1:spn) )
  if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or. &
     bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9)then
   allocate( spc0(1:spn,-1:in,-1:jn,-1:kn) )
  end if
 end if

! allocate external gravitational field if necessary
 if(include_extgrv)then
  allocate( extgrv(gis-2:gin,gjs-2:gjn,gks-2:gkn) )
 end if
 
 d = 1d0  ! set for numerical reasons
 T = 1d3  ! initial guess for temperature
 cs = 1d0
 b1 = 0d0; b2 = 0d0; b3 = 0d0
 
 return
end subroutine allocations
