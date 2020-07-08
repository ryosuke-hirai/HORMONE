!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE FAILEDSN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial conditions for a failed supernova

subroutine failedSN

 use physval
 use grid
 use constants
 use pressure_mod
 use gravmod,only:mc,grvphi,grvphiold,gis,gie,gjs,gje,gks,gke,lag11,lag12,lag21,lag31,lag32,h,dt_old
 use merger_mod

 implicit none

 integer lines, rows, iinj
 character*10000 dumc
 character*24,allocatable:: header(:),dum(:)
 real*8,allocatable,dimension(:,:):: dat
 real*8,allocatable,dimension(:):: m, r, pres,rho,ene, XX, YY, ZZ,gpot,temp
 real*8 Eexp, PNSmass, Ebind, mass, radius
 real*8 ecell, mnow, mold, rnow, dr, shell, shelld, dbg, fac, dnow, dold
!-----------------------------------------------------------------------------

 
! setting parameters for explosion ! -----------------------------------------
 !Eexp = 1.d51 ! erg
 PNSmass = 0d0 * msun
 dbg = 1d-14
 d(is:ie,js:je,ks:ke) = 0d0
 v1 = 0d0 ; v2 = 0d0 ; v3 = 0d0
 b1 = 0d0 ; b2 = 0d0 ; b3 = 0d0
fac = 1.21d0

! reading data from datafile ! -----------------------------------------------
 open(unit=40,file='../../100Msunprofile.data',status='old')
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

 allocate(m(0:lines),r(0:lines),pres(1:lines),rho(0:lines),ene(1:lines), &
          XX(1:lines),YY(1:lines),ZZ(1:lines),gpot(0:lines),temp(0:lines))
 do i = 1, rows
  if(trim(header(i))=='mass') m(1:lines) = dat(1:lines,i) * msun
  if(trim(header(i))=='density') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='energy') ene(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius') r(1:lines) = dat(1:lines,i) * rsun
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='temperature') temp(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='x_mass_fraction_H') XX(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='y_mass_fraction_He') YY(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='z_mass_fraction_metals') ZZ(1:lines) = dat(1:lines,i)
 end do

! start setting initial conditions ! -----------------------------------------

! reading composition data
 musize = lines
 allocate( mudata(1:lines,0:1) )
 mudata(1:lines,0) = m(1:lines)
 mudata(1:lines,1) = 4d0/(6d0*XX(1:lines)+YY(1:lines)+2d0)

 Ebind = 0d0
 do i = 2, lines
  Ebind = Ebind - G*m(i-1)*(m(i)-m(i-1)) / r(i-1) + ene(i)*(m(i)-m(i-1))
 end do
 mass = m(lines)
 radius = r(lines)

 Eexp = -Ebind * 0.5d0 !Eexp - Ebind
!print *,Eexp,Ebind
!stop
 mc(0) = 0d0
 r(0) = 0d0
 m(0) = 0d0
 rho(0) = rho(1)

! calculate gravitational potential
 do n = 0, lines-1
  rnow = 0d0;mnow = 0d0
  do i = n+1, lines
   mnow = mnow + (rho(i-1)*r(i)-rho(i)*r(i-1))/(r(i)-r(i-1))*0.5d0*(r(i)*r(i)-r(i-1)*r(i-1))+(rho(i)-rho(i-1))*(r(i)*r(i)+r(i)*r(i-1)+r(i-1)*r(i-1))/3d0!r(i)*(r(i)-r(i-1))
  end do
  gpot(n) = -G*m(n)/(r(n)+1d-99)-4d0*pi*G*mnow
 end do

 k = ks
 do i = is, ie
  if(xi1(i)>=radius)then
   mc(i:ie) = mass
   exit
  end if
  do n = 0, lines-1
   if(r(n+1)>xi1(i).and.r(n)<=xi1(i))then
    mc(i) = intpol(xi1(i),r(n:n+1),m(n:n+1))
    exit
   end if
  end do
 end do

 do j = js, je
  do i = is, ie
   d(i,j,k) = (mc(i)-mc(i-1))/sum(2d0*dvol(i,js:je,k))
   if(x1(i)<r(1))then
    p(i,j,k) = pres(1)
    spc(1,i,j,k) = XX(1)
    spc(2,i,j,k) = YY(1)
    spc(3,i,j,k) = ZZ(1)
   elseif(x1(i)<radius)then
    do n = 0, lines-1
     if(r(n+1)>x1(i).and.r(n)<=x1(i))then
      p(i,j,k) = intpol(x1(i),r(n:n+1),pres(n:n+1))
      spc(2,i,j,k) = intpol(x1(i),r(n:n+1),YY(n:n+1))
      spc(3,i,j,k) = intpol(x1(i),r(n:n+1),ZZ(n:n+1))
      spc(1,i,j,k) = 1d0-spc(2,i,j,k)-spc(3,i,j,k)
      exit
     end if
    end do
   else
    d(i,j,k) = dbg
    p(i,j,k) = G*mass*d(i,j,k)/x1(i)
    spc(1,i,j,k) = XX(lines)
    spc(2,i,j,k) = YY(lines)
    spc(3,i,j,k) = ZZ(lines)
   end if
  end do
 end do

 call meanmolweight
 k = ks
 do j = js, je
  do i = is, ie
   e(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
  end do
 end do

 call pressure

 k = ks
 Iinertia = 0d0
 do j = js, je
  do i = is, ie
   if(d(i,j,k)>dbg)then
    Iinertia = Iinertia + d(i,j,k)*dvol(i,j,k)*(x1(i)*sin(x2(j)))**2d0*2d0
   end if
  end do
 end do
! domega_dt = Jinject / Iinertia / Tinject

 do j = gjs-1, gje+1
  do i = gis-1, gie+1
   if(x1(i)<r(lines-1))then
    do n = 0, lines-1
     if(r(n+1)>x1(i).and.r(n)<=x1(i))then
      grvphi(i,j,k) = intpol(x1(i),r(n:n+1),gpot(n:n+1))
     end if
    end do
   else
    grvphi(i,j,k) = -G*mass/x1(i)
   end if
  end do
 end do

 grvphiold = grvphi


 v1=0d0;v2=0d0;v3=0d0
 b1=0d0;b2=0d0;b3=0d0

return

contains

 real*8 function intpol(x,rr,u)
  implicit none
  real*8 x
  real*8,dimension(1:2):: rr, u

  intpol = ( u(1)*(rr(2)-x) + u(2)*(x-rr(1)) ) / (rr(2)-rr(1))
 end function intpol

end subroutine failedSN
