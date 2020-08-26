program onedaverage

! purpose: To calculate latitudinal dependences

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
  real(8),allocatable:: lag(:,:)

  real(8) mass,bmass,rmass,res,res2, vsq, Egrav, Ekin, com, Erot, corr
  real(8) iniEbind,inimass,de_dt,domega_dt,theta,fac_egas,fac_pgas, pnow
  character*50 binfile, filename

  integer bndcll, kizami
  integer,parameter:: ie1d=300
  real(8),allocatable,dimension(:):: d1d, e1d, p1d, m1d, s1d, j1d, T1d
  real(8),allocatable:: spc1d(:,:),spcsum(:)
  real(8):: vsum,ssum,msum,esum,jsum,Tsum,imusum
  logical outflag
  logical,allocatable,dimension(:,:,:):: ub

! PREREQUISITES #############################################################

! read parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, sphrn,trnsn1,trnsn2,trnsn3
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn

  open(unit=1,file='../parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=eos_con)
  close(1)
ie = 1200

! allocate variables
  allocate( &
   x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
   x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
   x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2), &
   d(is:ie,js:je,ks:ke), e(is:ie,js:je,ks:ke), p(is:ie,js:je,ks:ke), &
   eint(is:ie,js:je,ks:ke), s(is:ie,js:je,ks:ke), T(is:ie,js:je,ks:ke), &
   imu(is:ie,js:je,ks:ke), &
   v1(is:ie,js:je,ks:ke), v2(is:ie,js:je,ks:ke), v3(is:ie,js:je,ks:ke), &
   b1(is:ie,js:je,ks:ke), b2(is:ie,js:je,ks:ke), b3(is:ie,js:je,ks:ke), &
   phi(is:ie,js:je,ks:ke), grvphi(is-1:ie+1,js:je,ks-1:ke+1), &
   grvphiold(is:ie,js:je,ks:ke), dvol(is-1:ie+1,js-1:je+1,ks-1:ke+1), &
   spc(1:spn,is:ie,js:je,ks:ke), ub(is:ie,js:je,ks:ke) )

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
 T = 1d3
 fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
 fac_pgas = kbol/amu ! frequently used factor for Pgas

! START ANALYSIS ############################################################

 do n = 812000, 812000, 100000
  
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

  bndcll = 0

  ub=.false.
  do k = ks,ke
   do j = js,je
    do i = is,ie
     imu(i,j,k) = 0.25d0*(6d0*spc(1,i,j,k)+spc(2,i,j,k)+2d0)
     vsq = v1(i,j,k)*v1(i,j,k) + v2(i,j,k)*v2(i,j,k) + v3(i,j,k)*v3(i,j,k)
     eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*vsq
     if(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)>0d0)then
      ub(i,j,k)=.true.
     end if
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
     if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<0d0)bndcll=bndcll+1
    end do
   end do
 end do

! write(filename,'(i)')n/100000
 open(unit=50,file='mass_dist.dat',status='replace')
 do j = js, je
  write(50,'(4(1PE13.5e2))')x2(j)/pi*180d0,&
   sum(d(is:ie,j,ks)*dvol(is:ie,j,ks),MASK=ub(is:ie,j,ks))/msun*200d0*2d0,&
   sum(e(is:ie,j,ks)*dvol(is:ie,j,ks),MASK=ub(is:ie,j,ks))*200d0*2d0,&
   sum(d(is:ie,j,ks)*v1(is:ie,j,ks)*dvol(is:ie,j,ks),MASK=ub(is:ie,j,ks))*200d0*2d0
 end do

end do

!!$ kizami = bndcll/ie1d
!!$ allocate( d1d(1:ie1d), e1d(1:ie1d), p1d(1:ie1d), m1d(0:ie1d), s1d(1:ie1d), &
!!$           spc1d(1:spn,1:ie1d), spcsum(1:spn), j1d(1:ie1d), T1d(1:ie1d) )
!!$ m1d = 0d0
!!$! pnow = 1d99!isobar averaged
!!$ pnow = 0d0 ! entropy sorted
!!$ do n = 1, ie1d
!!$  vsum=0d0;ssum=0d0;msum=0d0;esum=0d0;spcsum=0d0;jsum=0d0
!!$  if(n==ie1d)kizami = bndcll-(ie1d-1)*kizami
!!$  do vn = 1, kizami
!!$   pnow = minval(s(is:ie,js:je,ks),s(is:ie,js:je,ks)>pnow.and.e(is:ie,js:je,ks)+d(is:ie,js:je,ks)*grvphi(is:ie,js:je,ks)<0d0)! entropy sorted
!!$!   pnow = maxval(p(is:ie,js:je,ks),p(is:ie,js:je,ks)<pnow.and.e(is:ie,js:je,ks)+d(is:ie,js:je,ks)*grvphi(is:ie,js:je,ks)<0d0)! isobar averaged
!!$   outflag = .false.
!!$   do k = ks, ke
!!$    do j = js, je
!!$     do i = is, ie
!!$!      if(p(i,j,k)==pnow)then ! isobar averaged
!!$      if(s(i,j,k)==pnow)then ! entropy sorted
!!$       vsum = vsum + dvol(i,j,k)*2d0
!!$       msum = msum + d(i,j,k)*dvol(i,j,k)*2d0
!!$       esum = esum + eint(i,j,k)*dvol(i,j,k)*2d0
!!$       ssum = ssum + s(i,j,k)*d(i,j,k)*dvol(i,j,k)*2d0
!!$       jsum = jsum + d(i,j,k)*v3(i,j,k)*x1(i)*sin(x2(j))*dvol(i,j,k)*2d0
!!$       spcsum(1:spn) = spcsum(1:spn) + spc(1:spn,i,j,k)*d(i,j,k)*dvol(i,j,k)*2d0
!!$       outflag = .true.
!!$      end if;if(outflag)exit
!!$     end do;if(outflag)exit
!!$    end do;if(outflag)exit
!!$   end do
!!$  end do
!!$
!!$  d1d(n) = msum/vsum
!!$  e1d(n) = esum/msum
!!$  s1d(n) = ssum/msum
!!$  spc1d(1:spn,n) = spcsum(1:spn)/sum(spcsum(1:spn))
!!$  j1d(n) = jsum/msum
!!$  m1d(n) = m1d(n-1) + msum
!!$
!!$  imusum = 0.25d0*(6d0*spc1d(1,n)+spc1d(2,n)+2d0)
!!$
!!$  corr = 1d99;Tsum=1d3
!!$  do while(abs(corr)>eoserr*Tsum)
!!$   corr = (esum/vsum - ( arad*Tsum*Tsum*Tsum &
!!$                    + d1d(n)*fac_egas*imusum ) *Tsum )&
!!$        / ( -4d0*arad*Tsum*Tsum*Tsum &
!!$            -d1d(n)*fac_egas*imusum )
!!$   Tsum = Tsum - corr
!!$  end do
!!$
!!$  T1d(n) = Tsum
!!$  s1d(n) = imusum/amu*log(Tsum**1.5d0/d1d(n)) &
!!$         + 4d0*arad/3d0*Tsum*Tsum*Tsum/d1d(n)/kbol
!!$
!!$!  write(*,'(i5,13(1PE13.5e2))') n,m1d(n)/msun,d1d(n),e1d(n),s1d(n),spc1d(1:spn,n),j1d(n)
!!$ end do
!!$
!!$! open(unit=20,file='isobar_averaged_dist.dat',status='replace')
!!$! open(unit=30,file='isobar_averaged_comp.dat',status='replace')
!!$! open(unit=40,file='isobar_averaged_angm.dat',status='replace')
!!$ open(unit=20,file='entropy_sorted_dist.dat',status='replace')
!!$ open(unit=30,file='entropy_sorted_comp.dat',status='replace')
!!$ open(unit=40,file='entropy_sorted_angm.dat',status='replace')
!!$ write(20,'(i5)') ie1d-1
!!$ write(30,'(2i5)') ie1d-1, spn
!!$ write(40,'(i5)') ie1d-1
!!$ do n = ie1d-1, 1, -1
!!$  write(20,'(3(1PE14.6e2))') 1d0-m1d(n)/m1d(ie1d), d1d(n), T1d(n)
!!$  write(30,'(9(1PE14.6e2))') 1d0-m1d(n)/m1d(ie1d), spc1d(1:spn,n)
!!$  write(40,'(2(1PE14.6e2))') 1d0-m1d(n)/m1d(ie1d), j1d(n)
!!$ end do
!!$ close(20)
!!$ close(30)
!!$ close(40)
!!$
!!$ rmass = m1d(ie1d)
!!$print *,rmass/msun
!!$
!!$ n=1
!!$spc1d=1d0;m1d(n)=0d0
!!$ do i = is, ie
!!$  m1d(n) = m1d(n) + 2d0*sum(d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)
!!$  d1d(n) = sum(d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)/sum(dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)
!!$  e1d(n) = sum(e(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)/sum(d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)
!!$  s1d(n) = sum(s(i,js:je,ks)*d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)/sum(d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)
!!$  j1d(n) = x1(i)*sum(v3(i,js:je,ks)*sin(x2(js:je))*d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)/sum(d(i,js:je,ks)*dvol(i,js:je,ks),MASK=e(i,js:je,ks)+d(i,js:je,ks)*grvphi(i,js:je,ks)<0d0)
!!$!  write(*,'(i5,13(1PE13.5e2))') i,m1d(n)/msun,d1d(n),e1d(n),s1d(n),spc1d(1:spn,n),j1d(n)
!!$  if(m1d(n)>rmass)exit
!!$ end do

end program onedaverage
