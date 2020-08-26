program analysis

! purpose: To analyze the data for axisymmetric cylindrical coordinates.

  implicit none

  integer n,i,j,k,tn, is,ie,js,je,ks,ke,imesh,jmesh,kmesh, thetan
  integer,parameter:: nof = 2250000
  real*8,parameter:: pi = acos(-1d0), msun = 1.98855d33, gamma = 5d0/3d0,&
                     G = 6.67408d-8, rsun = 6.963d10
  real*8 x,y,dt_old,geta,time,mc, phih
  real*8 xi1s,xi1e,xi2s,xi2e,xi3s,xi3e
  real*8,allocatable,dimension(:,:,:):: d,v1,v2,v3,b1,b2,b3,e,p,phi,grvphi,&
                                        grvphiold, hg123, dvol, lapphi, divv
  real*8,allocatable,dimension(:):: x1, x2, x3, xi1, xi2, xi3,&
                                    dxi1, dxi2, dxi3, dx1, dx2, dx3, &
                                    hg11, hg12, hg21, hg22, hg31, hg32, idx1,idx3
  real*8,allocatable:: lag(:,:)

  real*8 mass,bmass,rmass,res,res2, vsq, Eint, Egrav, Ekin, com, Erot
  real*8 iniEbind,inimass,de_dt,domega_dt,theta
  character*50 binfile, pltfile
  real*8,dimension(1:18):: fm, fke,fv, fmom

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
   idx1(is-2:ie+2),idx3(ks-2:ke+2), &
   d (is:ie,js:je,ks:ke), e (is:ie,js:je,ks:ke), p (is:ie,js:je,ks:ke), &
   v1(is:ie,js:je,ks:ke), v2(is:ie,js:je,ks:ke), v3(is:ie,js:je,ks:ke), &
   b1(is:ie,js:je,ks:ke), b2(is:ie,js:je,ks:ke), b3(is:ie,js:je,ks:ke), &
   phi(is:ie,js:je,ks:ke), grvphi(is-1:ie+1,js:je,ks-1:ke+1), &
   lapphi(is:ie,js:je,ks:ke),divv(is:ie,js:je,ks:ke),&
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

 do n = 900000, nof,1000
  write(binfile,'(a3,i11.11,a5)')'bin',n,'s.dat'
  write(pltfile,'(a3,i11.11,a5)')'divv',n,'s.dat'
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

  open(unit=30,file=pltfile,status='replace')
  write(30,*)
  write(30,*)

  do k = ks, ke,2
   do j = js, je
    do i = is, ie,2
     divv(i,j,k) = v1(i,j,k)/x1(i)+&
          dx1(i)*idx1(i+1)/sum(dx1(i:i+1))*v1(i+1,j,k)-dx1(i+1)*idx1(i)/sum(dx1(i:i+1))*v1(i-1,j,k)+(dx1(i+1)-dx1(i))*idx1(i)*idx1(i+1)*v1(i,j,k) &
        + dx3(k)*idx3(k+1)/sum(dx3(k:k+1))*v3(i,j,k+1)-dx3(k+1)*idx3(k)/sum(dx3(k:k+1))*v3(i,j,k-1)+(dx3(k+1)-dx3(k))*idx3(k)*idx3(k+1)*v3(i,j,k)

     write(30,'(1PE13.5e2)') divv(i,j,k)
    end do
    if(ie/=1)write(30,*)
   end do
  end do

  close(30)


 end do

end program analysis
