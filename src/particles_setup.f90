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
 use particle_mod
 use ejectamod

 implicit none

 integer nn
 real*8 rptc, rhoc

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
