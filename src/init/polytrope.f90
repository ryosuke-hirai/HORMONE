module polytrope_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE POLYTROPE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up a polytropic model of a star

subroutine polytrope

 use constants,only:msun,rsun,G
 use grid
 use physval
 use star_mod


 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8)::mass,radius,mcore,rsoft,imu_const,dbg
 integer::i,j,k

!-----------------------------------------------------------------------------

 mass = msun
 radius = rsun
 mcore = 0.0d0
 rsoft = 0.0d0
 imu_const = 1.0d0

 call isentropic_star(mass,radius,mcore,rsoft,imu_const,m,r,rho,pres)

 ! Place the star at the origin
 call set_star_sph_grid(r,m,pres)

 ! Attach a wind-like atmosphere
 dbg = rho(size(rho)-1)

 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(d(i,j,k)<0d0)then
     d(i,j,k) = dbg*(radius/x1(i))**2
     p(i,j,k) = G*mass*d(i,j,k)/x1(i)
    end if
   end do
  end do
 end do

return
end subroutine polytrope

end module polytrope_mod
