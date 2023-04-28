module particle_mod

 use settings,only:include_particles,maxptc

 implicit none
 
 integer:: np, jmax, npl
 integer,allocatable,dimension(:,:):: ptci
 integer,allocatable,dimension(:):: ptc_in
 real(8),allocatable:: ptcx(:,:), cosiptc(:), sincptc(:), coscptc(:)
 
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE PARTICLES_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place the particles at initial position.

subroutine particles_setup

 use grid
 use physval
 use constants
 use ejectamod

 implicit none

!-----------------------------------------------------------------------------

! ptcx(0,n): mass of nth particle
! ptcx(1:dim,n): position vector of nth particle
! ptci(0,n): label of particle
! ptci(1,n): origin of particle (1 if star, 2 if ejecta)
! ptcl(2,n): status of particle (0 if bound, 1 if unbound)
! np: current number of particles
! npl: label for all particles that have been created
 np = 0

! for axisymmetrical cylindrical coordinates
 if(crdnt==1.and.dim==2.and.je==1)then
  j = js
!  rhoc = maxval(d(is:ie,js:je,ks:ke),mask=v1(is:ie,js:je,ks:ke)==0d0)
  allocate( ptcx(0:dim,1:maxptc), ptci(0:2,1:maxptc) )

  ! One particle in each cell for stellar matter
  do k = ks, ke
   do i = is, ie
!    if(d(i,j,k)>=rhoc/4d0.and.v1(i,j,k)==0d0)then
    if(v2(i,j,k)>0d0)then
     np = np+1
     if(np>maxptc)then
      print *,"number of particles exceeded maximum array size",np
      stop
     end if

     ptci(0,np) = np
     ptci(1,np) = 1 !star
     ptci(2,np) = 0 !bound
     ptcx(0,np) = d(i,j,k)*dvol(i,j,k)
     ptcx(1,np) = x1(i)
     ptcx(2,np) = x3(k)

    end if
   end do
  end do

  npl = np

! for axisymmetrical spherical coordinates
 elseif(crdnt==2.and.dim==2.and.ke==1)then
  k = ks
  allocate( ptcx(0:dim,1:maxptc), ptci(0:2,1:maxptc) )

  ! One particle in each cell for stellar matter
  do j = js, je
   do i = is, ie
    if(v3(i,j,k)>0d0)then
     np = np+1
     if(np>maxptc)then
      print *,"number of particles exceeded maximum array size",np
      stop
     end if

     ptci(0,np) = np
     ptci(1,np) = 1 ! star
     ptci(2,np) = 0 ! bound
     ptcx(0,np) = d(i,j,k)*dvol(i,j,k)
     ptcx(1,np) = x1(i)*sin(x2(j))
     ptcx(2,np) = x1(i)*cos(x2(j))

    end if
   end do
  end do

 else
  print *,'particles_setup not ready for this coordinate system',crdnt
  stop
 end if

return
end subroutine particles_setup


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE PARTICLES
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To follow trajectories of tracer particles.

subroutine particles

 use grid
 use physval
 use constants
 use ejectamod

 implicit none

 logical:: extflag
 real(8):: v1tmp1, v1tmp2, v2tmp1, v2tmp2, ptclr, ptclt

!-----------------------------------------------------------------------------

! for axisymmetrical cylindrical coordinates
 if(crdnt==1.and.dim==2.and.je==1)then
  
  n = 1 ; j = js

  do
   extflag = .false.
   do k = ks, ke+1
    if(ptcx(2,n)>=x3(k-1))then;if(ptcx(2,n)<x3(k))then
     do i = is, ie+1
      if(ptcx(1,n)>=x1(i-1))then;if(ptcx(1,n)<x1(i))then
       v1tmp1=(v1(i-1,j,k-1)*(x1(i)-ptcx(1,n))+v1(i,j,k-1)*(ptcx(1,n)-x1(i-1)))&
             * idx1(i)
       v2tmp1=(v3(i-1,j,k-1)*(x1(i)-ptcx(1,n))+v3(i,j,k-1)*(ptcx(1,n)-x1(i-1)))&
             * idx1(i)
       v1tmp2=(v1(i-1,j,k  )*(x1(i)-ptcx(1,n))+v1(i,j,k  )*(ptcx(1,n)-x1(i-1)))&
             * idx1(i)
       v2tmp2=(v3(i-1,j,k  )*(x1(i)-ptcx(1,n))+v3(i,j,k  )*(ptcx(1,n)-x1(i-1)))&
             * idx1(i)

       v1tmp1 = (v1tmp1*(x3(k)-ptcx(2,n))+v1tmp2*(ptcx(2,n)-x3(k-1))) * idx3(k)
       v2tmp1 = (v2tmp1*(x3(k)-ptcx(2,n))+v2tmp2*(ptcx(2,n)-x3(k-1))) * idx3(k)

       ptcx(1,n) = ptcx(1,n) + v1tmp1*dt
       ptcx(2,n) = ptcx(2,n) + v2tmp1*dt

       extflag = .true.
       exit

      end if;end if
     end do
     if(extflag)exit
    end if;end if
   end do

   ! to remove particles that have exceeded boundaries
   if(ptcx(1,n)>x1(ie).or.ptcx(2,n)>x3(ke))then
    ptcx(0:2,n) = ptcx(0:2,np)
    ptci(0:2,n) = ptci(0:2,np)
    np = np-1
    cycle
   end if

   if(n>=np)exit

   n = n+1

  end do

! for axisymmetrical spherical coordinates
 elseif(crdnt==2.and.dim==2.and.ke==1)then
! 1: equatorial direction
! 2: axial direction

  n = 1 ; k = ks

  do
   extflag = .false.
   ptclr = sqrt( ptcx(1,n)*ptcx(1,n)+ptcx(2,n)*ptcx(2,n) )
   ptclt = atan( ptcx(1,n)/ptcx(2,n) )-(sign(0.5d0,ptcx(2,n))-0.5d0)*pi

   do j = js, je+1
    if(ptclt>=x2(j-1))then;if(ptclt<x2(j))then
     do i = is, ie+1
      if(ptclr>=x1(i-1))then;if(ptclr<x1(i))then
       
       v1tmp1 =( v1(i-1,j-1,k)*(x1(i)-ptclr)+v1(i,j-1,k)*(ptclr-x1(i-1)) )&
              * idx1(i)
       v1tmp2 =( v1(i-1,j  ,k)*(x1(i)-ptclr)+v1(i,j  ,k)*(ptclr-x1(i-1)) )&
              * idx1(i)
       v2tmp1 =( v2(i-1,j-1,k)*(x1(i)-ptclr)+v2(i,j-1,k)*(ptclr-x1(i-1)) )&
              * idx1(i)
       v2tmp2 =( v2(i-1,j  ,k)*(x1(i)-ptclr)+v2(i,j  ,k)*(ptclr-x1(i-1)) )&
              * idx1(i)

       v1tmp1 = (v1tmp1*(x2(j)-ptclt)+v1tmp2*(ptclt-x2(j-1)))&
              * idx2(j)
       v2tmp1 = (v2tmp1*(x2(j)-ptclt)+v2tmp2*(ptclt-x2(j-1)))&
              * idx2(j)

       ptcx(1,n) = ptcx(1,n) + &
                  (v1tmp1*sin(ptclt)+v2tmp1*cos(ptclt)) * dt
       ptcx(2,n) = ptcx(2,n) + &
                  (v1tmp1*cos(ptclt)+v2tmp1*sin(ptclt)) * dt

       extflag = .true.
       exit

      end if;end if
     end do
     if(extflag)exit
    end if;end if
   end do

   ! to remove particles that have exceeded boundaries
   if(ptcx(1,n)**2d0+ptcx(2,n)**2d0>xi1e*xi1e.or.ptcx(1,n)**2d0+ptcx(2,n)**2d0<xi1s*xi1s)then
    ptcx(0:2,n) = ptcx(0:2,np)
    ptci(0:2,n) = ptci(0:2,np)
    np = np-1
    cycle
   end if

   if(n>=np)exit

   n=n+1

  end do

 else
  print *,'particles not ready for this coordinate system',crdnt
  stop
 end if

return
end subroutine particles


end module particle_mod
