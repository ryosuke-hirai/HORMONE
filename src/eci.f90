!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE ECI
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for head-on collision of stars.

subroutine eci

 use grid
 use settings,only:dt_out
 use physval
 use constants
 use utils,only:intpol
 use ejectamod
 use gravmod,only:extgrv,grvtime,include_extgrv,coremass
 use pressure_mod

 implicit none

 real(8):: mass, radius, dbg, corerad, c1, c2, c3, rhoh, rhoph, Mdot, vwind, mej
 real(8):: mnow, mold, rnow, shellv, shelld, dr, rsoft, msoft, shellp
 integer:: nn, ii, lines, rows
 character(10000):: dumc
 character(30),allocatable:: header(:),dum(:)
 real(8),allocatable,dimension(:,:):: dat, comp
 real(8),allocatable,dimension(:):: m, r, pres,rho,ene, XX, YY, ZZ, mumu, Temp
 logical setcore

!-----------------------------------------------------------------------------

 dbg = 0d0
 Mdot = 1d-6*msun/year
 vwind = 1d7
 d(is:ie,js:je,ks:ke) = 0d0
 spc(1,is-2:ie+2,js-2:je+2,ks-2:ke+2) = 0.73d0
 spc(2,is-2:ie+2,js-2:je+2,ks-2:ke+2) = 0.25d0
 spc(3,is-2:ie+2,js-2:je+2,ks-2:ke+2) = 0.02d0
! soft = 3d0*dx1(is)

! reading data from datafile ! -----------------------------------------------
 open(unit=40,file='16earlyRSG.data',status='old')
 read(40,'()')
 read(40,'()')
 read(40,*) lines, lines
 read(40,'()')
 read(40,'()')
 read(40,'(a)') dumc

! counting rows
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
101 do i = 1, 500
  if(dum(i)=='aaa')then
   rows = i - 1
   exit
  end if
 end do

 allocate(header(1:rows),dat(1:lines,1:rows))
 header(1:rows) = dum(1:rows)
 deallocate(dum)

 do i = 1, lines
  read(40,*) dat(lines-i+1,1:rows) 
 end do

 allocate(m(1:lines),r(1:lines),pres(1:lines),rho(1:lines),ene(1:lines), Temp(1:lines),&
          XX(1:lines),YY(1:lines),ZZ(1:lines),mumu(1:lines),comp(1:8,1:lines))
 do i = 1, rows
  if(trim(header(i))=='mass') m(1:lines) = dat(1:lines,i) * msun
  if(trim(header(i))=='density') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='energy') ene(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='temperature') Temp(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='mu') mumu(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius_cm') r(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='x_mass_fraction_H') XX(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='y_mass_fraction_He') YY(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='z_mass_fraction_metals') ZZ(1:lines) = dat(1:lines,i)
!  if(trim(header(i))=='total_energy') mumu(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='h1'  ) comp(1,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='he4' ) comp(2,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='he3' ) comp(3,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='c12' ) comp(4,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='n14' ) comp(5,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='o16' ) comp(6,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='ne20') comp(7,1:lines) = dat(1:lines,i)
  if(trim(header(i))=='mg24') comp(8,1:lines) = dat(1:lines,i)
 end do
!!$mass=1d3;radius=1d3
!!$ do i = 1, lines
!!$  rsoft=0d0
!!$!  call eos_p_cf(rho(i),rsoft,rsoft,rsoft,rsoft,rsoft,rsoft,ene(i)*rho(i),mass,shellv,coremass,corerad,XX(i),YY(i))
!!$!  call get_erec_imurec(rho(i),Temp(i),XX(i),YY(i),coremass,corerad)
!!$  mumu(i) = eos_e(rho(i),pres(i),mass,shellv,XX(i),YY(i))
!!$!  shelld = eos_p(rho(i),mumu(i),radius,shellv,XX(i),YY(i))
!!$  call eos_p_cf(rho(i),rsoft,rsoft,rsoft,rsoft,rsoft,rsoft,mumu(i),mass,shellv,shelld,corerad,XX(i),YY(i))
!!$!  mumu(i)=get_erec(rho(i),Temp(i),XX(i),YY(i))
!!$  print*,m(i)/msun,shelld,pres(i)
!!$!  print*,m(i)/msun,mumu(i)/rho(i),ene(i)
!!$!  print*,Temp(i),1d0/mumu(i),corerad
!!$ end do
!!$stop
! start setting initial conditions ! -----------------------------------------
! find He core surface
 do i = 1, lines
  if(YY(i)<0.4d0)then
   coremass = m(i)
   corerad  = r(i)
   rsoft = corerad*2d0
   exit
  end if
 end do
 mass = m(lines)
 radius = r(lines)
! coremass = 0d0 !coremass

! set chemical distribution
 j = js; setcore=.false.
 do k = ks, ke
  do i = is, ie
   if(rdis(i,k)<radius)then
    do n = 1, lines-1
     if(rdis(i,k)>r(n).and.rdis(i,k)<=r(n+1))then
      spc(2,i,j,k) = intpol(r(n:n+1),comp(2,n:n+1),rdis(i,k))
      spc(3,i,j,k) = intpol(r(n:n+1),comp(3,n:n+1),rdis(i,k))
      spc(4,i,j,k) = intpol(r(n:n+1),comp(4,n:n+1),rdis(i,k))
      spc(5,i,j,k) = intpol(r(n:n+1),comp(5,n:n+1),rdis(i,k))
      spc(6,i,j,k) = intpol(r(n:n+1),comp(6,n:n+1),rdis(i,k))
      spc(7,i,j,k) = intpol(r(n:n+1),comp(7,n:n+1),rdis(i,k))
      spc(8,i,j,k) = intpol(r(n:n+1),comp(8,n:n+1),rdis(i,k))
      spc(1,i,j,k) = 1d0-sum(spc(2:8,i,j,k))
      exit
     end if
    end do
   end if
  end do
 end do

 rnow = dxi1(is)*1.1
 dr = dxi1(is)*1.1

 ! find boundary value between softened core and envelope
 do n = 1, lines-1
  if(r(n)<rsoft.and.r(n+1)>=rsoft)then
   rhoh = intpol(r(n:n+1),rho(n:n+1),rsoft)
   rhoph = (rho(n+1)-rho(n))/(r(n+1)-r(n))
   shellp = intpol(r(n:n+1),pres(n:n+1),rsoft)
   msoft = intpol(r(n:n+1),m(n:n+1),rsoft)-coremass
   exit
  end if
 end do
 c1 = 2d0/(rsoft*rsoft)*rhoph - 10d0/rsoft**3d0*rhoh &
    + 7.5d0/(pi*rsoft**6d0)*msoft
 c2 = rhoph*0.5d0/rsoft-1.5d0*c1*rsoft
 c3 = rhoh-0.5d0*rsoft*rhoph+0.5d0*c1*rsoft**3d0
 shellv = 0d0 ; mold = 0d0 ; shellp = 0d0
 do while(rnow<=radius)
  ! first set distribution in softened core
  if(rnow<=rsoft)then
   mnow = 4d0*pi*rnow**3d0*(c1*rnow**3d0/6d0+0.2d0*c2*rnow**2+c3/3d0)
!   shellp = shellp - G*mnow/rnow/rnow*shelld*dr
  ! then set distribution in star
  else
   do n = 1, lines-1
    if(rnow>r(n).and.rnow<=r(n+1))then
     mnow = intpol(r(n:n+1),m(n:n+1),rnow)-coremass
!     shellp = intpol(rnow,r(n:n+1),pres(n:n+1))
     exit
    end if
   end do
  end if

  do k = ks, ke
   do i = is, ie
    if(rdis(i,k)<=rnow.and.d(i,j,k)==0d0)then
     shellv = shellv + dvol(i,j,k)
    end if
   end do
  end do
  shelld = (mnow-mold)/shellv
  if(shellv==0d0)then
   rnow = rnow + dr
   cycle
  end if
  if(rnow<=rsoft)then
   shellp = shellp - G*(mnow/(rnow*rnow)+coremass*rnow/rsoft**3d0)*shelld*dr
  else
   shellp = shellp - G*(mnow+coremass)/(rnow*rnow)*shelld*dr   
  end if
  do k = ks, ke
   do i = is, ie
    if(rdis(i,k)<=rnow.and.d(i,j,k)==0d0)then
     d(i,j,k) = shelld
     p(i,j,k) = shellp
    end if
   end do
  end do
  shellv = 0d0
  rnow = rnow + dr
  mold = mnow
 end do

 p(is:ie,js:je,ks:ke) = p(is:ie,js:je,ks:ke) - minval(p(is:ie,js:je,ks:ke))+pres(lines)
 
 shellv = 0d0
 do k = ks, ke
  do i = is, ie
   if(rdis(i,k)<=radius+dr.and.d(i,j,k)==0d0)then
    shellv = shellv + dvol(i,j,k)
   end if
  end do
 end do
 do k = ks, ke
  do i = is, ie
   if(rdis(i,k)<radius+dr.and.d(i,j,k)==0d0)then
    d(i,j,k) = (mass-coremass-mnow)/shellv
    p(i,j,k) = pres(lines)
   end if
  end do
 end do

 do k = ks, ke
  do i = is, ie
   if(d(i,j,k)==0d0)then
    d(i,j,k) = Mdot/(4d0*pi*rdis(i,k)*rdis(i,k)*vwind)
    p(i,j,k) = pres(lines)*1d-2
    v1(i,j,k) = vwind*sincyl(i,k)
    v3(i,j,k) = vwind*coscyl(i,k)
   end if
  end do
 end do

! set external gravitational field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 j = js

if(include_extgrv)then
 do k = gks-2, gke+2
  do i = is-2, gie+2
   if(rdis(i,k)<rsoft)then
    extgrv(i,j,k) = G*coremass/(2d0*rsoft)*((rdis(i,k)/rsoft)**2-3d0)
   else
    extgrv(i,j,k) = -G*coremass/rdis(i,k)
   end if
  end do
 end do
end if



open(unit=9191,file='data/extgrv.bin',status='replace',form='unformatted')
write(9191)coremass
write(9191)extgrv(gis-2:gie+2,gjs:gje,gks-2:gke+2)
close(9191)



! set ejecta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 time = tstart / ejectadistance * (sep-radius*1.01d0)
 t_out = (int(time / dt_out) + 1 )* dt_out

 j = js

 do k = ks,ke
  do i = is,ie

   if(nsdis(i,j,k)<=(sep-radius*1.01d0))then
    do nn = 1, count-1
     if(time*nsdfr(i,j,k)>=t_ej(nn).and.time*nsdfr(i,j,k)<t_ej(nn+1))then
      d(i,j,k)  = (d_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  d_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * nsdfr(i,j,k)**3.d0
!!$      p(i,j,k)  = (p_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
!!$                -  p_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
!!$                / (t_ej(nn+1)-t_ej(nn)) &
!!$                * nsdfr(i,j,k)**(3.d0)
      v1(i,j,k) = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * ( nssin(i,j,k))
      v3(i,j,k) = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * (-nscos(i,j,k))
      p(i,j,k) = 0.5d0*d(i,j,k)*(v1(i,j,k)**2+v3(i,j,k)**2)*1d-2
!!$      if(p_ej(nn)>pmax*1.d-3)then
!!$       p(i,j,k) = pmax !for shock condition
!!$!       e(i,j,k) = p(i,j,k) / (gamma-1.d0)
!!$      end if
      mej = (m_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
          -  m_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
          / (t_ej(nn+1)-t_ej(nn))
      do ii = 1, compsize
       if(mej>comp_ej(0,ii))then
        spc(2:8,i,j,k) = comp_ej(2:8,ii)
        spc(1,i,j,k) = 1d0-sum(spc(2:8,i,j,k))
        exit
       end if
      end do

     end if

    end do
   end if

  end do
 end do

 grvtime = time

 return
end subroutine eci

