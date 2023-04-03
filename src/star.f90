module star_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE REPLACE_CORE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To replace core with a constant entropy profile + point particle

subroutine replace_core(rcore,r,m,rho,pres,comp,comp_list)

 use constants,only:pi
 use composition_mod,only:get_imu

 real*8,intent(inout)::rcore
 real*8,allocatable,dimension(:),intent(inout):: r, m, rho, pres
 real*8,allocatable,dimension(:,:),intent(inout):: comp
 character(len=10),allocatable,intent(in):: comp_list(:)
 integer i,j,ih1,ihe4,ierr
 real*8:: mcore,mpt,imuh
 real*8,allocatable,dimension(:):: softr,softrho,softp
 
!-----------------------------------------------------------------------------

! get indices for hydrogen and helium
 do i = 1, size(comp_list)
  if(trim(comp_list(i))=='h1')ih1=i
  if(trim(comp_list(i))=='he4')ihe4=i
 end do

! find index for rcore
 do i = 1, size(pres)-1
  if(r(i)>rcore)exit
 end do
 
 rcore = r(i)
 mcore = m(i)
 imuh = get_imu((/comp(ih1,i),comp(ihe4,i)/))
 allocate(softr(0:i+1),softrho(0:i),softp(0:i+1))
 softr(0:i+1) = r(0:i+1)
 softrho(0:i) = rho(0:i)
 softp(0:i+1) = pres(0:i+1)

 call get_softened_profile(softr,mpt,mcore,imuh,softrho,softp,ierr)

 rho(0:i) = softrho(0:i)
 pres(0:i) = softp(0:i)
 
! recalculate mass coordinate
 m(0) = mpt
 do j = 1, i
  m(j) = m(j-1) + 4d0*pi/3d0*(r(j)**3-r(j-1)**3)*rho(j)
 end do
! Set everything inside rcore to uniform composition
 do j = 1, size(comp_list)
  comp(j,0:i-1) = comp(j,i)
 end do

return
end subroutine replace_core

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE SET_STAR_SPH_GRID
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place star at the origin of a spherical coordinate grid

subroutine set_star_sph_grid(r,m,rho,pres,comp,comp_list)

 use settings,only:spn,compswitch,eq_sym
 use constants,only:G,pi
 use grid
 use physval
 use gravmod,only:gravswitch,grvphi,grvphiold,mc
 use utils,only:intpol
 
 real*8,allocatable,dimension(:),intent(in):: r,m,rho,pres
 real*8,allocatable,dimension(:,:),intent(in),optional:: comp
 character(len=10),allocatable,intent(in),optional:: comp_list(:)
 real*8,allocatable,dimension(:)::gpot
 integer lines,nn,sn
 real*8:: Eexp,PNSmass, mass, radius
 real*8:: mnow,rnow,volfac
!-----------------------------------------------------------------------------

 lines = size(r)-1
 mass = m(lines)
 radius = r(lines)
 
! calculate gravitational potential
 allocate(gpot(0:lines))
 do n = 0, lines-1
  rnow = 0d0;mnow = 0d0
  do i = n+1, lines
   mnow = mnow &
        + (rho(i-1)*r(i)-rho(i)*r(i-1))/(r(i)-r(i-1))*0.5d0*(r(i)**2-r(i-1)**2)&
        + (rho(i)-rho(i-1))*(r(i)**2+r(i)*r(i-1)+r(i-1)**2)/3d0
  end do
  gpot(n) = -G*m(n)/(r(n)+1d-99)-4d0*pi*G*mnow
 end do

 k = ks
 do i = is, ie
  if(xi1(i)>=radius)then
   mc(i:ie) = mass
   exit
  end if
  do n = 0, lines-1
   if(r(n+1)>xi1(i).and.r(n)<=xi1(i))then
    mc(i) = intpol(r(n:n+1),m(n:n+1),xi1(i))
    exit
   end if
  end do
 end do

 if(eq_sym)then
  volfac=2d0
 else
  volfac=1d0
 end if
 
 do k = ks, ke
  do j = js, je
   do i = is, ie
    d(i,j,k) = (mc(i)-mc(i-1))/(volfac*sum(dvol(i,js:je,ks:ke)))
    if(x1(i)<r(1))then
     p(i,j,k) = pres(1)
     if(compswitch==2)spc(1:spn,i,j,k) = comp(1:spn,1)
    elseif(x1(i)<radius)then
     do n = 0, lines-1
      if(r(n+1)>x1(i).and.r(n)<=x1(i))then
       p(i,j,k) = intpol(r(n:n+1),pres(n:n+1),x1(i))
       if(compswitch==2)then
        do nn = 1, spn-1
         do sn = 1, size(comp_list)-1
          if(trim(comp_list(sn))==trim(species(nn)))then
           spc(nn,i,j,k) = intpol(r(n:n+1),comp(sn,n:n+1),x1(i))
           exit
          end if
         end do
        end do
        spc(spn,i,j,k) = 1d0-sum(spc(1:spn-1,i,j,k)) !dump the rest into others
       end if
       exit
      end if
     end do
    else
     d(i,j,k) = -1d0 ! Set negative for easy identification
    end if
   end do
  end do
 end do
 
!!$ p(ie+1,js:je,ks:ke) = 1d-99
!!$ do i = ie, is, -1
!!$  p(i,js:je,ks:ke) = p(i+1,js:je,ks:ke) + G*mc(i)*max(d(i,js,ks),rho(lines)*1d-5)/xi1(i)**2*dx1(i+1)
!!$ end do

 do j = gjs-1, gje+1
  do i = gis-1, gie+1
   if(x1(i)<r(lines-1))then
    do n = 0, lines-1
     if(r(n+1)>x1(i).and.r(n)<=x1(i))then
      grvphi(i,j,k) = intpol(r(n:n+1),gpot(n:n+1),x1(i))
     end if
    end do
   else
    grvphi(i,j,k) = -G*mass/x1(i)
   end if
  end do
 end do

 if(gravswitch==3)grvphiold = grvphi

return
end subroutine set_star_sph_grid

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE SET_STAR_CYL_GRID
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place star on a cylindrical coordinate grid

subroutine set_star_cyl_grid

 use grid

 

!-----------------------------------------------------------------------------

 

return
end subroutine set_star_cyl_grid
 

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE ONE_SHOT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculate a hydrostatic structure for a given entropy

subroutine one_shot(Sc,imu,r,mcore,msoft,rho,p,mass)

 use constants,only:G,pi
 use pressure_mod,only:get_d_from_ps
 use utils,only:softened_acc
 
 real*8,intent(in)::Sc,mcore,msoft
 real*8,intent(inout)::imu
 real*8,allocatable,dimension(:),intent(in)::r
 real*8,allocatable,dimension(:),intent(inout)::rho,p
 real*8,intent(out)::mass

 integer i,j,Nmax
 real*8::hsoft
 real*8,allocatable,dimension(:)::dr,vol

!-----------------------------------------------------------------------------

 Nmax=size(rho)-1
 allocate(dr(1:Nmax+1),vol(1:Nmax+1))
 do i = 1, Nmax+1
  dr(i) = r(i)-r(i-1)
  vol(i) = 4d0*pi/3d0*(r(i)**3-r(i-1)**3)
 end do
 hsoft=r(Nmax)

 mass=msoft

 do i = Nmax, 1, -1
  p(i-1)=(dr(i)*dr(i+1)*sum(dr(i:i+1))&
        *rho(i)*G*(mass/r(i)**2+mcore*softened_acc(r(i),hsoft))&
        +dr(i)**2*p(i+1) &
        +(dr(i+1)**2-dr(i)**2)*p(i))/dr(i+1)**2
  rho(i-1) = get_d_from_ps(p(i-1),Sc,imu)
  mass=mass-rho(i)*vol(i)
  if(mass<0d0)return
 end do

return
end subroutine one_shot

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                     SUBROUTINE GET_SOFTENED_PROFILE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Returns softened core profile with fixed entropy

subroutine get_softened_profile(r,mpt,mh,imuh,rho,p,ierr)

 use settings,only:eostype
 use pressure_mod,only:entropy_from_dp

 real*8,allocatable,dimension(:),intent(in)::r
 real*8,intent(in)::mh,imuh
 real*8,intent(inout)::mpt
 real*8,allocatable,dimension(:),intent(inout)::rho,p
 integer,intent(out)::ierr

 integer Nmax,i
 real*8::Sc,mass,mold,msoft,fac,Sedge,T,imu,eostype0

!-----------------------------------------------------------------------------

! Instructions

! input variables should be given in the following format
 
! r(0:Nmax+1): Array of radial grid. Should be set so that r(0)=0 and r(Nmax)=hsoft
! mpt: Core particle mass
! mh: Mass coordinate at hsoft
! imuh: 1/mu at hsoft (mu is mean molecular weight)
! rho(0:Nmax): Give rho(Nmax)=(rho at hsoft) as input. Outputs density profile.
! p(0:Nmax+1): Give p(Nmax:Nmax+1)=(p at r(Nmax:Nmax+1)) as input. Outputs pressure profile.

! ierr: Is set to 1 when we cannot find the solution.

! This module does not work with non-ideal EoSs
 eostype0 = eostype
 if(eostype>=2) eostype = 1
 
 ierr=0
 mpt=mh*0.5d0 ! initial guess for point particle mass
 msoft=mh-mpt
 Nmax=size(rho)-1
 imu = imuh
 T = 1d3
 Sedge=entropy_from_dp(rho(Nmax),p(Nmax),T,imu)

! Start shooting method
 fac=0.05d0
 mass=msoft
 Sc=Sedge
 do i = 1, 500
  mold=mass
  call one_shot(Sc,imu,r,mpt,msoft,rho,p,mass)
  if(mass<0d0)then
   mpt=mpt*(1d0-fac)
   msoft=mh-mpt
  elseif(mass/msoft<1d-10)then
   exit
  else
   mpt=mpt*(1d0+fac)
   msoft=mh-mpt
  end if
  if(mold*mass<0d0)fac=fac*0.5d0
 end do

 if(i>500) ierr = 1

 eostype = eostype0 ! set back to original eostype

return
end subroutine get_softened_profile


end module star_mod
