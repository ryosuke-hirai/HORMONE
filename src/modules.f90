!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                  MODULES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Modules

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module settings

 implicit none

! boundary conditions
 integer:: bc1is, bc1os, bc2is, bc2os, bc3is, bc3os
 integer:: bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov
 logical:: eq_sym, dirichlet_on, fluxbound_on
! numerical setups
 integer:: rktype, crdnt, tnlim, start, tn_out, tn_evo, outstyle, endstyle
 integer:: gravswitch, compswitch, radswitch
 integer:: eostype, spn, sigfig, outres, gbtype, grktype, maxptc
 integer:: grvsrctype, opacitytype, lambdatype
 real(8):: courant, t_end, dt_out, dt_unit_in_sec
 character(len=5):: dt_unit
 real(8):: grverr, cgerr, eoserr, HGfac, hgcfl
 integer:: imesh, jmesh, kmesh
! switches
 logical:: include_extgrv, include_particles, include_cooling, mag_on
 logical:: include_extforce, include_sinks, is_test
 logical:: write_other_vel, write_shock, grav_init_other, write_evo
 logical:: write_other_slice, write_temp
 character(len=30):: flux_limiter, simtype
 character(len=50):: parafile,extrasfile

end module settings


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module grid

 use settings,only:crdnt

  implicit none

! number of grids
  integer:: is, ie, js, je, ks, ke
  integer:: gis, gie, gjs, gje, gks, gke

  integer:: is_global, ie_global, js_global, je_global, ks_global, ke_global
  integer:: gis_global, gie_global, gjs_global, gje_global, gks_global, gke_global

  integer:: tn,dim
  integer:: rungen
  integer:: musize,fmr_max,fmr_lvl(0:20)
! grid center = x, grid interface = xi
  real(8),allocatable,dimension(:):: x1, xi1, dx1, dxi1, idx1, idxi1
  real(8),allocatable,dimension(:):: x2, xi2, dx2, dxi2, idx2, idxi2
  real(8),allocatable,dimension(:):: x3, xi3, dx3, dxi3, idx3, idxi3
  real(8):: time, dt, t_out
  real(8):: xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, x1min, x2min, x3min
  real(8),allocatable,dimension(:):: sinc, sini, cosc, cosi
  real(8),allocatable,dimension(:,:,:):: dvol, idetg3
  real(8),allocatable,dimension(:):: detg1, idetg1, sx1, g22, scot, sisin
  real(8),allocatable,dimension(:,:):: detg2, idetg2, g33
  real(8),allocatable,dimension(:,:,:):: sa1, sa2, sa3, Imom
  real(8),allocatable,dimension(:,:):: rdis, sincyl, coscyl
  real(8),allocatable,dimension(:,:,:,:):: car_x
  real(8),allocatable,dimension(:):: spinc_r,spinc_t

end module grid


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module constants

  implicit none

! all in cgs units
  real(8),parameter:: pi = acos(-1d0)
  real(8),parameter:: G  = 6.67428d-8, clight = 2.99792458d10
  real(8),parameter:: msun = 1.989d33, rsun = 6.963d10
  real(8),parameter:: msolar = msun, rsolar = rsun, au = 1.49598073d13
  real(8),parameter:: kbol = 1.38064852d-16, amu = 1.6605402d-24
  real(8),parameter:: arad = 7.5646d-15, sigma = 5.67051d-5
  real(8),parameter:: hplanck = 6.6260755d-27, m_p = 1.6726231d-24
  real(8),parameter:: m_n = 1.6749286e-24, m_e = 9.1093897e-28
  real(8),parameter:: N_A = 6.0221367e23
  real(8),parameter:: year = 3600d0*24d0*365.25d0
  real(8),parameter:: huge = 1d99, tiny = 1d-99
  real(8):: fac_pgas, fac_egas

end module constants


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module physval

  implicit none

  integer:: icnt,imo1,imo2,imo3,iene,img1,img2,img3,i9wv,ufnmax
  real(8),allocatable,dimension(:,:,:):: d, p, e, v1, v2, v3, b1, b2, b3, ptot
  real(8),allocatable,dimension(:,:,:):: T, eint, imu
  real(8),allocatable,dimension(:,:,:,:):: dd, de, dm1, dm2, dm3
  real(8),allocatable,dimension(:,:,:,:):: db1, db2, db3, dphi, dmu
  real(8),allocatable,dimension(:,:,:):: cs, phi, grv1, grv2, grv3
  real(8),allocatable,dimension(:,:,:,:):: u, flux1, flux2, flux3, uorg, src
  real(8),allocatable,dimension(:,:):: mudata
  real(8),allocatable,dimension(:,:,:,:):: spc, spcorg
  real(8),allocatable,dimension(:,:,:,:,:):: dspc, spcflx
  character(len=10),allocatable:: species(:)
  real(8),allocatable,dimension(:,:,:):: d0,p0,b10,b20,b30,v10,v20,v30
  real(8),allocatable,dimension(:,:,:,:):: spc0

  real(8):: gamma, muconst, ch

  integer,allocatable,dimension(:,:,:):: shock

end module physval


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module gravmod

  use settings,only:gravswitch,grverr,cgerr,include_extgrv,hgcfl,HGfac,&
                    grav_init_other,gbtype

  implicit none

  integer,parameter:: llmax = 1000
  real(8),allocatable,dimension(:,:,:),target:: grvphi, grvphiold, grvphidot, totphi
  real(8),allocatable,dimension(:,:,:,:):: grvphiorg
  real(8),allocatable,dimension(:,:):: Pl
  real(8),allocatable,dimension(:,:,:):: Plc
  real(8),allocatable,dimension(:,:):: phiio, phiii, phi1o, phi3i, phi3o
  real(8),target:: dt_old, grvtime, dtgrav, cgrav2, cgrav
  real(8),allocatable,dimension(:):: hg11,hg12,hg21,hg22,hg31,hg32, mc
  real(8),allocatable,dimension(:,:):: lag
  real(8),allocatable,dimension(:,:,:):: hg123,orgdis, extgrv, hgsrc, gsrc
  real(8),allocatable,dimension(:,:,:,:):: lag11,lag12,lag21,lag22,lag31,lag32
  real(8):: coremass,hg_dx

end module gravmod

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module ejectamod

  implicit none

  integer:: count, ejdatnum, tstartn, compsize, j_ejmax
  real(8),allocatable,dimension(:):: t_ej, d_ej, p_ej, e_ej, v_ej, m_ej
  real(8),allocatable,dimension(:,:):: comp_ej
  real(8),allocatable,dimension(:,:,:):: nsdis, nsdfr, nssin, nscos
  real(8):: tstart, pmax, psmass, ejectadistance, sep
  character(len=40):: ejtbinfile

end module ejectamod

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
