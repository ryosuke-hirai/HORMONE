module radstar_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE RADSTAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up a polytropic model of a star with radiation hydro

subroutine radstar

 use settings,only:eostype
 use constants,only:msun,rsun,arad,G
 use grid
 use physval
 use star_mod
 use pressure_mod,only:eos_e

 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8)::mass,radius,mcore,rsoft,imu_const,dbg,pbg
 integer::i,j,k

!-----------------------------------------------------------------------------

 mass = 20d0*msun
 radius = 4d0*rsun
 mcore = 0.0d0
 rsoft = 0.0d0
 imu_const = 1.0d0

 call isentropic_star(mass,radius,mcore,rsoft,imu_const,m,r,rho,pres)

 ! Place the star at the origin
 call set_star_sph_grid(r,m,pres)

! Attach a wind-like atmosphere
! (Hardwire values to avoid compiler-dependent atmospheres)
 dbg = 1d-5
 pbg = 2d0*G*mass/(3d0*radius**2*0.3d0)

 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(d(i,j,k)<0d0)then
     d(i,j,k) = dbg*(radius/x1(i))**2
     p(i,j,k) = pbg*(radius/x1(i))**2
    end if
   end do
  end do
 end do

! Set radiation pressure
 eostype=1
 do k = ks, ke
  do j = js, je
   do i = is, ie
    eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    erad(i,j,k) = arad*T(i,j,k)**4
    p(i,j,k) = p(i,j,k) - erad(i,j,k)/3d0
   end do
  end do
 end do
 eostype=0

return
end subroutine radstar

end module radstar_mod
