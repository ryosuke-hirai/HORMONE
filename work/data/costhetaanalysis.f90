program analysis

! purpose: To analyze the data for axisymmetric spherical coordinates.

  implicit none

  integer n,i,j,k,tn,vn, is,ie,js,je,ks,ke,imesh,jmesh,kmesh, thetan,sphrn,trnsn1,trnsn2,trnsn3
  integer,parameter:: nof = 1000000, bins=20
  real*8,parameter:: pi = acos(-1d0), msun = 1.98855d33, gamma = 5d0/3d0,&
                     G = 6.67408d-8, rsun = 6.963d10
  real*8 x,y,dt_old,geta,time,mc, phih
  real*8 xi1s,xi1e,xi2s,xi2e,xi3s,xi3e
  real*8,allocatable,dimension(:,:,:):: d,v1,v2,v3,b1,b2,b3,e,p,phi,grvphi,&
                                        grvphiold, hg123, dvol, lapphi
  real*8,allocatable,dimension(:):: x1, x2, x3, xi1, xi2, xi3,&
                                    dxi1, dxi2, dxi3, dx1, dx2, dx3, &
                                    hg11, hg12, hg21, hg22, hg31, hg32
  real*8,allocatable:: lag(:,:)

  real*8 mass,bmass,rmass,res,res2, vsq, Eint, Egrav, Ekin, com, Erot
  real*8 iniEbind,inimass,de_dt,domega_dt,theta
  character*50 binfile
  real*8,dimension(1:bins):: fv, fke, fmom
  real*8,dimension(0:5,1:bins):: fm

! PREREQUISITES #############################################################

! read parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, sphrn,trnsn1,trnsn2,trnsn3
  open(unit=1,file='../parameters',status='old')
  read(1,NML=gridcon)
  close(1)
! allocate variables
  allocate( &
   x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
   x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
   x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2), &
   d (is:ie,js:je,ks:ke), e (is:ie,js:je,ks:ke), p (is:ie,js:je,ks:ke), &
   v1(is:ie,js:je,ks:ke), v2(is:ie,js:je,ks:ke), v3(is:ie,js:je,ks:ke), &
   b1(is:ie,js:je,ks:ke), b2(is:ie,js:je,ks:ke), b3(is:ie,js:je,ks:ke), &
   phi(is:ie,js:je,ks:ke), grvphi(is-1:ie+1,js:je,ks-1:ke+1), &
   lapphi(is:ie,js:je,ks:ke),&
   grvphiold(is:ie,js:je,ks:ke), dvol(is-1:ie+1,js-1:je+1,ks-1:ke+1), &
   hg11(0:ie+1),hg12(0:ie+1),hg21(0:ie+1),hg22(0:ie+1),&
   hg31(0:ke+1),hg32(0:ke+1),hg123(0:ie+1,1:1,0:ke+1), lag(-1:1,0:ie+1) )

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
!stop

! START ANALYSIS ############################################################

! open(unit=10,file='timeevo.dat',status='replace')

! write(10,'(11a14)') "time", "mass", "ber_mass","bound_mass", "Ekin", "Eint", "Egrav", "Poisson_err","Poisson_err2", "geta", "centroid"

 open(unit=30,file='costheta_distribution.dat',status='replace')
 write(30,'(a6,8a13)')'mu','m(v>0km/s)','m(v>100km/s)','m(v>200km/s)','m(v>300km/s)','m(v>400km/s)','m(v>500km/s)','KE','momentum'

 do n = 1000000, nof,1000
  write(binfile,'(a3,i11.11,a5)')'bin',n,'s.dat'
  open(unit=20,file=binfile,status='old',form='unformatted',err=500)
  read(20) tn,time,iniEbind,inimass,de_dt,domega_dt
  read(20) d (is:ie,js:je,ks:ke), &
           v1(is:ie,js:je,ks:ke), &
           v2(is:ie,js:je,ks:ke), &
           v3(is:ie,js:je,ks:ke), &
!           b1(is:ie,js:je,ks:ke), &
!           b2(is:ie,js:je,ks:ke), &
!           b3(is:ie,js:je,ks:ke), &
           e (is:ie,js:je,ks:ke), &
!           phi(is:ie,js:je,ks:ke), &
           grvphi(is:ie,js:je,ks:ke)
  read(20) grvphiold(is:ie,js:je,ks:ke), &
           dt_old
  read(20) spc(1:spn,is:ie,js:je,ks:ke)
  close(20)

  fm = 0d0; fke=0d0; fv=0d0; fmom = 0d0
  do k = ks, ke
   do j = js, je
    do i = is, ie
     theta = cos(x2(j))!x3(k)/sqrt(x1(i)*x1(i)+x3(k)*x3(k))
     thetan = int(theta*bins)+1

  !   if(d(i,j,k)*grvphi(i,j,k)+e(i,j,k)>0d0)then
      if    (v1(i,j,k)>5d7)then
       fm(0:5,thetan) = fm(0:5,thetan) + d(i,j,k)*dvol(i,j,k)
      elseif(v1(i,j,k)>4d7)then
       fm(0:4,thetan) = fm(0:4,thetan) + d(i,j,k)*dvol(i,j,k)
      elseif(v1(i,j,k)>3d7)then
       fm(0:3,thetan) = fm(0:3,thetan) + d(i,j,k)*dvol(i,j,k)
      elseif(v1(i,j,k)>2d7)then
       fm(0:2,thetan) = fm(0:2,thetan) + d(i,j,k)*dvol(i,j,k)
      elseif(v1(i,j,k)>1d7)then
       fm(0:1,thetan) = fm(0:1,thetan) + d(i,j,k)*dvol(i,j,k)
!      elseif(v1(i,j,k)>0d7)then
      else
       fm(0,thetan) = fm(0,thetan) + d(i,j,k)*dvol(i,j,k)
      end if
      fmom(thetan)=fmom(thetan) + d(i,j,k)*v1(i,j,k)*dvol(i,j,k)
      fke(thetan)= fke(thetan)+ 0.5d0*d(i,j,k)*(v1(i,j,k)*v1(i,j,k)+v2(i,j,k)*v2(i,j,k)+v3(i,j,k)*v3(i,j,k))*dvol(i,j,k)
  !   end if

    end do
   end do
  end do

  write(30,'(a1,i)')'#',n-1000000
  do i = 1, bins
   write(30,'(f6.3,8(1PE13.5e2))')dble(i)/dble(bins)-0.5d0/dble(bins),fm(0:5,i)*2d0/msun,fke(i),fmom(i)
  end do
  write(30,'()')
!  print *,sum(fm(0,1:bins))*2d0/msun

!  write(10,'(11(1pe14.6e2))') time, mass/msun, bmass/msun,rmass/msun, Ekin, Erot, Eint, Egrav,com/msun!, res,res2, geta, com
500 cycle
 end do

 close(10)

end program analysis
