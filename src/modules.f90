!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                  MODULES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Modules

module funcs
 implicit none
 
 contains

  real*8 function pw(p,a)
   implicit none
   integer:: i
   integer,intent(in):: p
   real*8,intent(inout):: a

   pw = a
   do i = 2, p
    pw = pw * a
   end do
  end function pw
 
end module funcs



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module settings

  implicit none

! boundary conditions
  integer bc1is, bc1os, bc2is, bc2os, bc3is, bc3os
  integer bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov
  logical eq_sym
! numerical setups
  integer rktype, crdnt, tnlim, gravswitch, start, tn_out, outstyle, endstyle
  integer eostype, spn, compswitch
  real*8 courant, t_end, dt_out
  real*8 grverr, cgerr, eoserr
  integer imesh, jmesh, kmesh
! switches
  logical include_extgrv, include_particles

end module settings



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module grid

  use settings,only:crdnt

  implicit none

! number of grids
  integer is, ie, js, je, ks, ke, in, jn, kn
  integer i,j,k,n,tn,dim
  integer rungen, ufn
  integer musize,sphrn,trnsn1,trnsn2,trnsn3
! grid center = x, grid interface = xi
  real*8,allocatable,dimension(:):: x1, xi1, dx1, dxi1, idx1, idxi1
  real*8,allocatable,dimension(:):: x2, xi2, dx2, dxi2, idx2, idxi2
  real*8,allocatable,dimension(:):: x3, xi3, dx3, dxi3, idx3, idxi3
  real*8 time, dt, t_out
  real*8 xi1s, xi1e, xi2s, xi2e, xi3s, xi3e
  real*8,allocatable,dimension(:):: sinc, sini, cosc, cosi

end module grid


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module constants

  implicit none

! all in cgs units
  real*8,parameter:: pi = acos(-1d0)
  real*8,parameter:: G  = 6.67428d-8
  real*8,parameter:: msun = 1.989d33, rsun = 6.963d10, &
                     msolar = msun, rsolar = rsun
  real*8,parameter:: kbol = 1.38064852d-16, amu = 1.6605402d-24
  real*8,parameter:: arad = 7.5646d-15, year = 24d0*36d2*365.25d0
  real*8:: fac_pgas, fac_egas

end module constants


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module physval
  
  use grid,only:in,jn,kn

  implicit none

  real*8,allocatable,dimension(:,:,:):: d, p, e, v1, v2, v3, b1, b2, b3, ptot
  real*8,allocatable,dimension(:,:,:):: T, eint, imu
  real*8,allocatable,dimension(:,:,:,:):: dd, de, dm1, dm2, dm3 
  real*8,allocatable,dimension(:,:,:,:):: db1, db2, db3, dphi, dmu
  real*8,allocatable,dimension(:,:,:):: cf, phi
  real*8,allocatable,dimension(:,:,:,:):: u, flux1, flux2, flux3, uorg, src
  real*8,allocatable,dimension(:,:):: mudata
  real*8,allocatable,dimension(:,:,:,:):: spc, spcorg
  real*8,allocatable,dimension(:,:,:,:,:):: dspc, spcflx

  real*8 gamma, muconst

  real*8,allocatable,dimension(:):: detg1, idetg1, sx1, g22, scot, sisin
  real*8,allocatable,dimension(:,:):: detg2, idetg2, g33
  real*8,allocatable,dimension(:,:,:):: idetg3, dvol

end module physval


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module ninewave

  implicit none

  real*8 ch

end module ninewave


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module gravmod

  use grid,only:ie,je,ke,in,jn,kn
  use settings,only:gravswitch,grverr,cgerr,include_extgrv

  implicit none

  integer l, ll, lmax
  integer gis, gie, gjs, gje, gks, gke, gin, gjn, gkn
  integer,parameter:: llmax = 1000
  integer,allocatable,dimension(:):: modlimax
  real*8,allocatable,dimension(:,:,:):: grvphi, grvphiold, grv1, grv2, grv3, phicg
  real*8,allocatable,dimension(:):: a1,a2,a3, preca, precb, precc, precd, prece
  real*8,allocatable,dimension(:):: x,y,z,r,aw, mc
  real*8,allocatable,dimension(:,:):: Pl
  real*8,allocatable,dimension(:,:,:):: Plc
  real*8,allocatable,dimension(:,:):: phiio, phiii, phi1o, phi3i, phi3o
  real*8 ,dimension(0:llmax):: ml
  real*8,allocatable,dimension(:,:):: rdis, sincyl, coscyl
  real*8 dt_old, HGfac, hgcfl, l2norm, grvtime
  real*8,allocatable,dimension(:):: hg11,hg12,hg21,hg22,hg31,hg32
  real*8,allocatable,dimension(:,:):: lag
  real*8,allocatable,dimension(:,:,:):: hg123,orgdis, extgrv, hgsrc

!experimental
  real*8 h
  real*8,allocatable,dimension(:,:,:,:):: lag11,lag12,lag21,lag22,lag31,lag32

end module gravmod


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module dirichlet

  use grid,only:in,jn,kn

  implicit none

  real*8,allocatable,dimension(:,:,:):: d0,p0,b10,b20,b30,v10,v20,v30

end module dirichlet


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module amr_templates

  implicit none

  type amr_header
   integer :: parnt, level
   integer,allocatable,dimension(:):: child, neigh, ldif
   integer :: ref_deref
  end type amr_header

  type leaf_contents
   real*8,allocatable,dimension(:):: x1, xi1, dx1, idx1!dxi1, idxi1
   real*8,allocatable,dimension(:):: x2, xi2, dx2, idx2!dxi2, idxi2
   real*8,allocatable,dimension(:):: x3, xi3, dx3, idx3!dxi3, idxi3
   real*8,allocatable,dimension(:):: detg1, idetg1, g22, sx1, scot, sisin
   real*8,allocatable,dimension(:,:):: detg2, idetg2, g33
   real*8,allocatable,dimension(:,:,:):: dvol, idetg3
   real*8,allocatable,dimension(:,:,:):: d, p, e, v1, v2, v3, b1, b2, b3
   real*8,allocatable,dimension(:,:,:):: ptot, cf, phi, gphi
   real*8,allocatable,dimension(:,:,:,:):: u, uorg, flux1, flux2, flux3, src
   real*8,allocatable,dimension(:,:,:,:):: tflx1, tflx2, tflx3
   real*8,allocatable,dimension(:,:,:,:):: midu1, midu2, midu3
  end type leaf_contents

end module amr_templates


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module amr_module

  use amr_templates
  use grid,only: dim

  implicit none

  integer maxamr, cuts, cib, face
  integer ib, jb, kb, bkn, totbloks, lfn, totleafs, lid, curlvl, lvl
  integer,allocatable,dimension(:):: regrid
  type(amr_header),allocatable,dimension(:):: bk, obk
  type(leaf_contents),allocatable,dimension(:):: lf, olf
  integer,dimension(8,6)::localneigh
  real*8,allocatable,dimension(:):: schedule
  real*8 dt_local
  logical sync

end module amr_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module ejectamod

  implicit none

  integer count, ejdatnum, tstartn
  real*8,allocatable,dimension(:):: t_ej, d_ej, p_ej, e_ej, v_ej
  real*8,allocatable,dimension(:,:,:):: nsdis, nsdfr, nssin, nscos
  real*8 tstart, pmax, psmass, ejectadistance, sep
  character*40 ejtbinfile

end module ejectamod

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module particle_mod

 use settings,only:include_particles

 implicit none
 
 integer np, maxptc, jmax, npl
 integer,allocatable,dimension(:,:):: ptci
 integer,allocatable,dimension(:):: ptc_in
 real*8,allocatable:: ptcx(:,:), cosiptc(:), sincptc(:), coscptc(:)
 
end module particle_mod

