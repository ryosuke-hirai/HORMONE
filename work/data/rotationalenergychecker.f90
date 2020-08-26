program rotationcheck

! purpose: To calculate rotational contribution to binding energy

  implicit none

  integer n,i,j,k,tn,vn, is,ie,js,je,ks,ke,imesh,jmesh,kmesh, &
          thetan,sphrn,trnsn1,trnsn2,trnsn3, spn, eostype, compswitch
  real(8),parameter:: pi = acos(-1d0), msun = 1.98855d33, gamma = 5d0/3d0,&
                     G = 6.67408d-8, rsun = 6.963d10, &
                     amu = 1.6605402d-24, kbol = 1.380658d-16, arad = 7.5646d-15
  real(8) x,y,dt_old,geta,time,mc, phih, eoserr, muconst
  real(8) xi1s,xi1e,xi2s,xi2e,xi3s,xi3e
  real(8),allocatable,dimension(:,:,:):: d,v1,v2,v3,b1,b2,b3,e,p,phi,grvphi,&
                                        grvphiold, dvol, lapphi,eint,imu, s, T
  real(8),allocatable,dimension(:,:,:,:):: spc
  real(8),allocatable,dimension(:):: x1, x2, x3, xi1, xi2, xi3,&
                                    dxi1, dxi2, dxi3, dx1, dx2, dx3
  real(8),allocatable:: lag(:,:), spin_coeffr(:), spin_coefft(:)

  real(8) mass,bmass,rmass,res,res2, vsq, Egrav, Ekin, com, Erot, corr
  real(8) iniEbind,inimass,de_dt,domega_dt,theta,fac_egas,fac_pgas, pnow
  real(8) Eenv, Erotenv, mcore, menv
  character*50 binfile

  integer bndcll, kizami
  integer,parameter:: ie1d=100
  real(8),allocatable,dimension(:):: d1d, e1d, p1d, m1d, s1d, j1d, T1d, kin1d
  real(8),allocatable:: spc1d(:,:),spcsum(:)
  real(8):: vsum,ssum,msum,esum,jsum,Tsum,imusum
  logical outflag
  logical,allocatable:: bnd(:,:,:)

! PREREQUISITES #############################################################

! read parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, sphrn,trnsn1,trnsn2,trnsn3
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn

  open(unit=1,file='../parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=eos_con)
  close(1)
ie=1200
! allocate variables
  allocate( &
   x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
   x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
   x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2), &
   d(is:ie,js:je,ks:ke), e(is:ie,js:je,ks:ke), p(is:ie,js:je,ks:ke), &
   eint(is:ie,js:je,ks:ke), s(is:ie,js:je,ks:ke), T(is:ie,js:je,ks:ke), &
   imu(is:ie,js:je,ks:ke), bnd(is:ie,js:je,ks:ke),&
   v1(is:ie,js:je,ks:ke), v2(is:ie,js:je,ks:ke), v3(is:ie,js:je,ks:ke), &
   b1(is:ie,js:je,ks:ke), b2(is:ie,js:je,ks:ke), b3(is:ie,js:je,ks:ke), &
   phi(is:ie,js:je,ks:ke), grvphi(is-1:ie+1,js:je,ks-1:ke+1), &
   grvphiold(is:ie,js:je,ks:ke), dvol(is-1:ie+1,js-1:je+1,ks-1:ke+1), &
   spc(1:spn,is:ie,js:je,ks:ke) )

! read coordinate data
  open(unit=1,file='gridfile2.bin',form='unformatted',status='old')
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

 allocate(spin_coeffr(is:ie),spin_coefft(js:je))
 do i = is, ie
  spin_coeffr(i) = 0.75d0*(xi1(i)*xi1(i)+xi1(i-1)*xi1(i-1))*(xi1(i)+xi1(i-1))/(xi1(i)*xi1(i)+xi1(i)*xi1(i-1)+xi1(i-1)*xi1(i-1))
 end do
 do j = js, je
  spin_coefft(j) = 0.25d0*(2d0*dxi2(j)-sin(2d0*xi2(j))+sin(2d0*xi2(j-1)))/(cos(xi2(j-1))-cos(xi2(j)))
 end do

 T = 1d3
 fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
 fac_pgas = kbol/amu ! frequently used factor for Pgas

! START ANALYSIS ############################################################

 n = 1772000

 write(binfile,'(a3,i11.11,a5)')'bin',n,'s.dat'
 open(unit=20,file=binfile,status='old',form='unformatted')
 read(20) tn,time,iniEbind,inimass,de_dt,domega_dt
 read(20) d (is:ie,js:je,ks:ke), &
      v1(is:ie,js:je,ks:ke), &
      v2(is:ie,js:je,ks:ke), &
      v3(is:ie,js:je,ks:ke), &
      e (is:ie,js:je,ks:ke), &
      grvphi(is:ie,js:je,ks:ke)
 read(20) grvphiold(is:ie,js:je,ks:ke), &
          dt_old
 if(compswitch>=2)then
  read(20)spc(1:spn,is:ie,js:je,ks:ke)
 end if
 close(20)

! Calculate other physical quantities
 bndcll = 0
 bnd = .false.
 do k = ks,ke
  do j = js,je
   do i = is,ie
    imu(i,j,k) = 0.25d0*(6d0*spc(1,i,j,k)+spc(2,i,j,k)+2d0)
    vsq = v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k)
    eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq
    corr = 1d99
    do while(abs(corr)>eoserr*T(i,j,k))
     corr = (eint(i,j,k) - ( arad*T(i,j,k)*T(i,j,k)*T(i,j,k) &
                           + d(i,j,k)*fac_egas*imu(i,j,k) ) *T(i,j,k) )&
          / ( -4d0*arad*T(i,j,k)*T(i,j,k)*T(i,j,k) &
              -d(i,j,k)*fac_egas*imu(i,j,k) )
     T(i,j,k) = T(i,j,k) - corr
    end do
    p(i,j,k) = ( d(i,j,k)*fac_pgas*imu(i,j,k) & ! gas pressure
               + arad*T(i,j,k)**3d0/3d0 ) *T(i,j,k) ! radiation pressure
    s(i,j,k) = imu(i,j,k)/amu*log(T(i,j,k)**1.5d0/d(i,j,k)) &
             + 4d0*arad/3d0*T(i,j,k)**3d0/d(i,j,k)/kbol ! in erg/g/K
    if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
     bndcll=bndcll+1
     bnd(i,j,k) = .true.
    end if
   end do
  end do
 end do

 Eenv = 0d0; mcore = 0d0

 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(bnd(i,j,k).and.spc(1,i,j,k)<=0.25d0)then
     Eenv = Eenv + d(i,j,k)*v3(i,j,k)*dvol(i,j,k)*spin_coeffr(i)*spin_coefft(j)*2d0
     mcore = mcore + d(i,j,k)*dvol(i,j,k)*2d0
    end if
   end do
  end do
 end do

 print *,Eenv,mcore/msun

end program rotationcheck
