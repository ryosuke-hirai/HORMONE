program trajectory

! purpose: To trace the trajectory of selected particles

  implicit none

  integer,parameter:: nof=2660000, nop=50
  integer maxptc, n, i,j,k, np, npl
  integer,dimension(1:nop):: selp, selp0
  integer,allocatable,dimension(:,:):: ptci
  real*8,allocatable,dimension(:,:):: ptcx
  real*8,parameter:: msun = 1.989d33, rsun = 6.963d10, &
                     msolar = msun, rsolar = rsun, pi = acos(-1d0), &
                     sampler = 50d0*rsun
  real*8 radius, theta, dtheta, drnow, thnow, xnow, znow
  character*40 ptcfile, format1
  logical include_particles

! PREREQUISITES #############################################################

! read parameters
  namelist /partcon/ include_particles, maxptc
  open(unit=1,file='../parameters',status='old')
  read(1,NML=partcon)
  close(1)

  allocate( ptcx(0:2,0:maxptc), ptci(0:2,1:maxptc) )
  dtheta = 0.5d0*pi/dble(nop)
  ptcx(0:2,0) = 1d0/0d0
  write(format1,'(a,i,a)')'(i10,',nop,'(i3,2(1PE13.5e2)))'

  open(unit=20,file='particle_trajectories.dat',status='replace')

  do n = 1000000, nof, 1000
   write(ptcfile,'(a3,i11.11,a5)') 'bpt',n,'s.dat'
   open(unit=10,file=ptcfile,status='old',form='unformatted')
   read(10)np,npl
   read(10)ptci(0:2,1:np), ptcx(0:2,1:np)
   close(10)

! select particles on a certain surface initially (set selp(1:nop))
   if(n==1000000)then
    do k = 1, nop
     drnow = 1d99
     theta = dble(k)*dtheta
     xnow = sampler*cos(theta)
     znow = sampler*sin(theta)
     do i = 1, np
      radius = sqrt(ptcx(1,i)*ptcx(1,i)+ptcx(2,i)*ptcx(2,i))
      if(radius>=sampler-0.9d0.and.radius<sampler*1.1d0)then
       thnow  = atan(ptcx(2,i)/ptcx(1,i))
       if(thnow<theta.and.thnow>=theta-dtheta)then
        if((xnow-ptcx(1,i))**2d0+(znow-ptcx(2,i))**2d0<drnow)then
         selp0(k) = i
         drnow = (ptcx(1,selp(k))-ptcx(1,i))**2d0+(ptcx(2,selp(k))-ptcx(2,i))**2d0
        end if
       end if
      end if
     end do
    end do
   end if

! write down trajectory

   selp(1:nop) = 0
   do k = 1, nop
    do i = 1, np
     if(ptci(0,i)==selp0(k))then
      selp(k) = i
     end if
    end do
   end do

   write(20,format1) n, ((ptci(2,k),ptcx(1:2,selp(k))),k=1,nop)


  end do

  close(20)


end program trajectory
