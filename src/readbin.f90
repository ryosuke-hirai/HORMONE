module readbin_mod
 implicit none
 public:: allocate_readgrid,readbin
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READBIN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate variables and read gridfile

 subroutine allocate_readgrid

  use settings
  use grid
  use physval
  use constants
  use gravmod
  use ionization_mod,only:ionization_setup

  implicit none

 !-----------------------------------------------------------------------------

! read parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                     sphrn, trnsn1, trnsn2, trnsn3
  namelist /out_con/ outstyle, endstyle, tnlim, t_end, dt_out, tn_out, &
                     dt_unit, sigfig, write_other_vel, write_shock
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn, include_cooling
  namelist /simucon/ crdnt,courant, rktype, start, mag_on, flux_limiter
  namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                     bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
  namelist /gravcon/ gravswitch, grverr, cgerr, HGfac, hgcfl, &
                     include_extgrv, gis, gie, gjs, gje, gks, gke
  
  open(unit=1,file='../parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=out_con)
  read(1,NML=eos_con)
  read(1,NML=simucon)
  read(1,NML=bouncon)
  read(1,NML=gravcon)
  close(1)

! allocate variables
  call checksetup
  call allocations

!!$  allocate( &
!!$   x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
!!$   x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
!!$   x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2), &
!!$   d(is:ie,js:je,ks:ke), e(is:ie,js:je,ks:ke), p(is:ie,js:je,ks:ke), &
!!$   eint(is:ie,js:je,ks:ke), T(is:ie,js:je,ks:ke), &
!!$   imu(is:ie,js:je,ks:ke), &
!!$   v1(is:ie,js:je,ks:ke), v2(is:ie,js:je,ks:ke), v3(is:ie,js:je,ks:ke), &
!!$   b1(is:ie,js:je,ks:ke), b2(is:ie,js:je,ks:ke), b3(is:ie,js:je,ks:ke), &
!!$   phi(is:ie,js:je,ks:ke), grvphi(is-1:ie+1,js:je,ks-1:ke+1), &
!!$   grvphiold(is:ie,js:je,ks:ke), dvol(is-1:ie+1,js-1:je+1,ks-1:ke+1), &
!!$   spc(1:spn,is:ie,js:je,ks:ke) )

! read coordinate data
  open(unit=1,file='gridfile.bin',form='unformatted',status='old')
  read(1) x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
          x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
          x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2)
  close(1)
  
! calculate volume element  
 do k = ks, ke
  do j = js-1, je+1
   do i = is-1, ie+1
    dvol(i,j,k)   = (xi1(i)**3d0-xi1(i-1)**3d0) / 3.d0 &
                  * (cos(xi2(j-1))-cos(xi2(j))) * dxi3(k)
!    dvol(i,j,k) = 5.d-1 * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi2(j) * dxi3(k)
!    if(je==1)dvol(i,j,k) = pi * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi3(k)
   end do
  end do
 end do
 if(ke==1) dvol = 4d0 * dvol
 T = 1d3
 gamma = 5d0/3d0
 fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
 fac_pgas = kbol/amu ! frequently used factor for Pgas

 if(eostype==2)then
  call ionization_setup
 end if

 return
end subroutine allocate_readgrid

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READBIN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read physval's from a binary file

subroutine readbin(filename)

 use settings
 use grid
 use physval
 use pressure_mod
 use composition_mod
 use gravmod

 implicit none

 character(len=*),intent(in):: filename
!-----------------------------------------------------------------------------

 open(unit=20,file=filename,status='old',form='unformatted')
 read(20) tn,time
 read(20) d (is:ie,js:je,ks:ke), &
          v1(is:ie,js:je,ks:ke), &
          v2(is:ie,js:je,ks:ke), &
          v3(is:ie,js:je,ks:ke), &
          e (is:ie,js:je,ks:ke)
 if(gravswitch>=2)read(20)grvphi(gis:gie,gjs:gje,gks:gke)
 if(gravswitch==3)then
  read(20)grvphiold(gis:gie,gjs:gje,gks:gke), &
          dt_old
 end if
 if(compswitch>=2)read(20)spc(1:spn,is:ie,js:je,ks:ke)
 if(mag_on)then
  read(20) b1(is:ie,js:je,ks:ke), &
           b2(is:ie,js:je,ks:ke), &
           b3(is:ie,js:je,ks:ke), &
           phi(is:ie,js:je,ks:ke)
 end if
 close(20)

 call meanmolweight
 call pressure

return
end subroutine readbin

end module readbin_mod
