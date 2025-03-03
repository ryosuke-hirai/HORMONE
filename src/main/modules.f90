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
 integer:: gravswitch, compswitch, radswitch, frame
 integer:: eostype, spn, sigfig, outres, gbtype, grktype, maxptc
 integer:: grvsrctype, opacitytype, lambdatype
 real(8):: courant, t_end, dt_out, dt_unit_in_sec, alpha9wave
 character(len=5):: dt_unit
 real(8):: grverr, cgerr, eoserr, HGfac, hgcfl, alphagrv
 integer:: imesh, jmesh, kmesh
 integer:: nsink, maxtngrv
! test tolerance
 real(8):: test_tol, Mach_tol
! switches
 logical:: solve_i, solve_j, solve_k
 logical:: include_extgrv, include_particles, include_cooling, mag_on
 logical:: include_extforce, include_sinks, include_accretion, is_test
 logical:: write_other_vel, write_shock, write_evo, write_other_slice
 logical:: write_temp, write_mc, output_ascii
 logical:: grav_init_other, grav_init_relax
 logical:: in_loop
 character(len=30):: flux_limiter, simtype
 character(len=50):: parafile,extrasfile, outdir

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
  real(8):: frame_acc(1:3)

end module grid


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module constants

  implicit none

! all in cgs units

! Exact values
  real(8),parameter:: pi = acos(-1d0), clight = 2.99792458d10
  real(8),parameter:: kbol = 1.380649d-16, hplanck = 6.62607015d-27
  real(8),parameter:: sigma=2d0*pi**5*kbol**4/(15d0*clight**2*hplanck**3)
  real(8),parameter:: arad = 4d0*sigma/clight
! Measured values
  real(8),parameter:: G = 6.67430d-8
  real(8),parameter:: Navo = 6.02214076d23, amu = 1d0/Navo
  real(8),parameter:: m_p = 1.6726231d-24, m_n = 1.6749286d-24
  real(8),parameter:: m_e = 9.1093897d-28, Rgas = kbol*Navo
  real(8),parameter:: year = 3600d0*24d0*365.25d0
! IAU 2015 Resolution B3
  real(8),parameter:: msun = 1.3271244d26/G, rsun = 6.957d10
  real(8),parameter:: au = 1.49597870700d13
! Others
  real(8),parameter:: huge = 1d99, tiny = 1d-99
  real(8):: fac_egas

end module constants


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module physval

  implicit none

  integer:: icnt,imo1,imo2,imo3,iene,img1,img2,img3,i9wv,ufnmax
  real(8),allocatable,dimension(:,:,:):: d, p, e, v1, v2, v3, b1, b2, b3, ptot
  real(8),allocatable,dimension(:,:,:):: T, eint, erad, imu
  real(8),allocatable,dimension(:,:,:,:):: dd, de, dm1, dm2, dm3
  real(8),allocatable,dimension(:,:,:,:):: db1, db2, db3, dphi, dmu, der
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
  real(8),allocatable,dimension(:,:,:):: grvphi,grvpsi, totphi,lapphi
  real(8),allocatable,dimension(:,:,:,:):: grvphiorg, lap_coeff
  real(8),allocatable,dimension(:,:):: Pl
  real(8),allocatable,dimension(:,:,:):: Plc
  real(8),allocatable,dimension(:,:):: phiio, phiii, phi1o, phi3i, phi3o
  real(8),target:: grvtime, dtgrav, cgrav2, cgrav, cgrav_old
  real(8),allocatable,dimension(:):: mc
  real(8),allocatable,dimension(:,:,:):: orgdis, extgrv, hgsrc, gsrc
  real(8):: coremass,dtg_unit

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
module derived_types

 implicit none

 type sink_prop
  sequence
  integer:: i,j,k,pad ! pad is added to align memory for 64-bit machines
  real(8):: mass, softfac, lsoft, laccr, locres, dt, mdot, racc, facc, jet_ang
  real(8),dimension(1:3):: x,v,a,xpol,Jspin,jdot,jet_dir
 end type sink_prop

contains

 subroutine null_sink(sink)
  type(sink_prop),intent(out):: sink
  sink%i=0
  sink%j=0
  sink%k=0
  sink%pad=0
  sink%mass=0d0
  sink%softfac=0d0
  sink%lsoft=0d0
  sink%laccr=0d0
  sink%locres=0d0
  sink%dt=0d0
  sink%mdot=0d0
  sink%racc=0d0
  sink%facc=0d0
  sink%jet_ang=0d0
  sink%x=0d0
  sink%v=0d0
  sink%a=0d0
  sink%xpol=0d0
  sink%Jspin=0d0
  sink%jdot=0d0
  sink%jet_dir=[0d0,0d0,1d0]
 end subroutine null_sink

end module derived_types
