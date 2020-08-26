program analysis

! purpose: To analyze the data for axisymmetric cylindrical coordinates.

  implicit none

  integer n,i,j,k,tn, is,ie,js,je,ks,ke,imesh,jmesh,kmesh, thetan, inow, know,jnow
  integer,parameter:: nof = 1000
  real*8,parameter:: pi = acos(-1d0), msun = 1.98855d33, gamma = 5d0/3d0,&
                     G = 6.67408d-8, rsun = 6.963d10
  real*8 x,y,dt_old,geta,time,mc, phih
  real*8 xi1s,xi1e,xi2s,xi2e,xi3s,xi3e
  real*8,allocatable,dimension(:,:,:):: d,v1,v2,v3,b1,b2,b3,e,p,phi,grvphi,&
                                        grvphiold, hg123, dvol, lapphi, rdis
  real*8,allocatable,dimension(:):: x1, x2, x3, xi1, xi2, xi3,&
                                    dxi1, dxi2, dxi3, dx1, dx2, dx3, &
                                    hg11, hg12, hg21, hg22, hg31, hg32
  real*8,allocatable:: lag(:,:)

  real*8 mass,bmass,rmass,res,res2, vsq, Eint, Egrav, Ekin, com, rnow, rold
  real*8 iniEbind,inimass,de_dt,domega_dt,theta, shellv, shellm, shellvel, shellber, shellvth
  character*50 binfile, outfile
  real*8,dimension(1:18):: fm, fke,fv, fmom
  logical out


! PREREQUISITES #############################################################

! read parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh
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
   lapphi(is:ie,js:je,ks:ke), rdis(is:ie,js:je,ks:ke),&
   grvphiold(is:ie,js:je,ks:ke), dvol(is-1:ie+1,js-1:je+1,ks-1:ke+1), &
   hg11(0:ie+1),hg12(0:ie+1),hg21(0:ie+1),hg22(0:ie+1),&
   hg31(0:ke+1),hg32(0:ke+1),hg123(0:ie+1,1:1,0:ke+1), lag(-1:1,0:ie+1) )

! read coordinate data
  open(unit=1,file='gridfile.bin',form='unformatted',status='old')
  read(1) x1(is-2:ie+2),xi1(is-2:ie+2),dxi1(is-2:ie+2),dx1(is-2:ie+2), &
          x2(js-2:je+2),xi2(js-2:je+2),dxi2(js-2:je+2),dx2(js-2:je+2), &
          x3(ks-2:ke+2),xi3(ks-2:ke+2),dxi3(ks-2:ke+2),dx3(ks-2:ke+2)
  close(1)

  do k = ks, ke
   do i = is, ie
    rdis(i,js,k) = sqrt(x1(i)*x1(i)+x3(k)*x3(k))
   end do
  end do

! calculate volume element  
 do k = ks-1, ke+1
  do j = js-1, je+1
   do i = is-1, ie+1
    dvol(i,j,k) = 5.d-1 * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi2(j) * dxi3(k)
    if(je==1)dvol(i,j,k) = pi * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi3(k)
   end do
  end do
 end do

! calculate Laplacian
 do i = is-1, ie+1
  hg11(i) = 1d0/dx1(i+1)/sum(dx1(i:i+1))
  hg12(i) = 1d0/dx1(i  )/sum(dx1(i:i+1))
 end do
 do i = is, ie+1
  hg21(i) = 1d0/(dx1(i+1)*dx1(i+1)) ! 1/h**2
  hg22(i) = sqrt(x1(i)*x1(i)+dx1(i+1)*dx1(i+1)) ! x(h)
! for Lagrange interpolation (second order)
  lag(-1,i) = (hg22(i)-x1(i))*(hg22(i)-x1(i+1))/dx1(i)/sum(dx1(i:i+1))
  lag( 0,i) = -(hg22(i)-x1(i-1))*(hg22(i)-x1(i+1))/dx1(i)/dx1(i+1)
  lag( 1,i) = (hg22(i)-x1(i-1))*(hg22(i)-x1(i))/sum(dx1(i:i+1))/dx1(i+1)
 end do
 do k = ks-1, ke+1
  hg31(k) = 1d0/(dx3(k+1)*sum(dx3(k:k+1))) 
  hg32(k) = 1d0/(dx3(k  )*sum(dx3(k:k+1))) 
 end do
 do k = ks-1, ke+1
  do i = is-1, ie+1
   hg123(i,js,k) = 1d0/dx1(i+1)*( 1d0/dx1(i)+1d0/dx1(i+1) ) + 1d0/(dx3(k)*dx3(k+1))
  end do
 end do

! START ANALYSIS ############################################################

do thetan = 0,85,5

 write(outfile,'("1Daveraged_",i2.2,"-",i2.2,"deg.data")')thetan,thetan+5

 open(unit=10,file=outfile,status='replace')

 do n = 900000, 2950000, 10000
  write(binfile,'(a3,i11.11,a5)')'bin',n,'s.dat'
  open(unit=20,file=binfile,status='old',form='unformatted')
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
  close(20)

  grvphi(is-1,js:je,ks:ke)= grvphi(is,js:je,ks:ke)

  mass = sum( d(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke) )

  bmass = 0d0
  i=is;j=js;k=ks
!!$  do k = ks, ke
!!$!  do i = is, ie
!!$   bmass = bmass + 4d0/3d0*pi*(xi3(k)**3d0-xi3(k-1)**3d0)*d(is,js,k)/msun
!!$!   bmass = bmass + 4d0/3d0*pi*(xi1(i)**3d0-xi1(i-1)**3d0)*d(i,js,ks)/msun
!!$   write(10,'(i8,5(1PE13.5e2))') n, xi3(k), d(i,j,k),v3(i,j,k), grvphi(i,j,k)+(e(i,j,k)+p(i,j,k))/d(i,j,k), bmass
!!$!   write(10,'(i8,5(1PE13.5e2))') n, xi1(i), d(i,j,k),v1(i,j,k), grvphi(i,j,k)+(e(i,j,k)+p(i,j,k))/d(i,j,k), bmass
!!$
!!$  end do
!!$  write(10,'()')

!!$  com = 2d12;j=js
!!$  do
!!$   rnow = 1d99
!!$   do k = ks, ke
!!$    do i = is, ie
!!$     if(rnow>rdis(i,j,k))then
!!$      if(rdis(i,j,k)>com)then
!!$       if(atan(x3(k)/x1(i))>89d0/180d0*pi.and.atan(x3(k)/x1(i))<=90d0/180d0*pi)then
!!$        rnow = rdis(i,j,k)
!!$        inow = i
!!$        know = k
!!$       end if
!!$      end if
!!$     end if
!!$    end do
!!$   end do
!!$   i=inow;k=know
!!$   write(10,'(i8,2i5,6(1PE13.5e2))') n,i,k,rdis(i,j,k), d(i,j,k),x1(i)*v1(i,j,k)+x3(k)*v3(i,j,k), grvphi(i,j,k)+(e(i,j,k)+p(i,j,k))/d(i,j,k), bmass
!!$   com = rdis(i,j,k)
!!$   if(com>x3(ke))exit
!!$  end do
!!$  write(10,'()')
!!$ end do

  j=js; jnow = 9
  rold = 0d0
  do
   shellv = 0d0;shellm = 0d0;shellvel = 0d0; shellber = 0d0; shellvth = 0d0
   rnow = xi1(jnow)
   do k = ks, ke
    do i = is, ie
     if(rnow>=rdis(i,j,k))then
      if(rdis(i,j,k)>rold)then
       if(atan(x3(k)/x1(i))>dble(85-thetan)/180d0*pi.and.atan(x3(k)/x1(i))<=dble(90-thetan)/180d0*pi)then
        shellv = shellv + dvol(i,j,k)
        shellm = shellm + d(i,j,k)*dvol(i,j,k)
        shellvel = shellvel + d(i,j,k)*dvol(i,j,k)*(x1(i)*v1(i,j,k)+x3(k)*v3(i,j,k))/rdis(i,j,k)
        shellvth = shellvth + d(i,j,k)*dvol(i,j,k)*(-x3(k)*v1(i,j,k)+x1(i)*v3(i,j,k))/rdis(i,j,k)
        shellber = shellber + dvol(i,j,k)*(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)+p(i,j,k))
       end if
      end if
     end if
    end do
   end do
   write(10,'(i8,i5,5(1PE13.5e2))') n,jnow,0.5d0*(rnow+rold), shellm/shellv,shellvel/shellm, shellvth/shellm, shellber/shellm
   jnow = jnow + 3
   rold = rnow
   if(jnow >= ie)exit
  end do
  write(10,'()')
 end do

 close(10)
end do

end program analysis
