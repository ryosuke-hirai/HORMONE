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

 use settings,only:compswitch
 use constants,only:pi
 use composition_mod,only:get_imu
 use physval,only:muconst

 real(8),intent(inout)::rcore
 real(8),allocatable,dimension(:),intent(inout):: r, m, rho, pres
 real(8),allocatable,dimension(:,:),intent(inout):: comp
 character(len=10),allocatable,intent(in):: comp_list(:)
 integer::i,j,ih1,ihe4,ierr
 real(8):: mcore,mpt,imuh
 real(8),allocatable,dimension(:):: softr,softrho,softp

!-----------------------------------------------------------------------------

 if(rcore<=0d0)return

! get indices for hydrogen and helium
 if(compswitch==2)then
  do i = 1, size(comp_list)
   if(trim(comp_list(i))=='h1')ih1=i
   if(trim(comp_list(i))=='he4')ihe4=i
  end do
 end if

! find index for rcore
 do i = 1, size(pres)-1
  if(r(i)>rcore)exit
 end do

 rcore = r(i)
 mcore = m(i)
 if(compswitch==2)then
  imuh = get_imu([comp(ih1,i),comp(ihe4,i)])
 else
  imuh = 1d0/muconst
 end if
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
 if(compswitch==2)then
  do j = 1, size(comp_list)
   comp(j,0:i-1) = comp(j,i)
  end do
 end if

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
 use gravmod,only:grvphi,mc
 use utils,only:intpol
 use mpi_domain,only:sum_global_array

 real(8),allocatable,dimension(:),intent(in):: r,m,rho,pres
 real(8),allocatable,dimension(:,:),intent(in),optional:: comp
 character(len=10),allocatable,intent(in),optional:: comp_list(:)
 real(8),allocatable,dimension(:)::gpot,Vshell
 integer::i,j,k,n,lines,nn,sn
 real(8):: mass, radius
 real(8):: mnow,volfac
 real(8):: term1, term2, term3, term4

!-----------------------------------------------------------------------------

 lines = size(r)-1
 mass = m(lines)
 radius = r(lines)

! calculate gravitational potential
 allocate(gpot(0:lines))
 do n = 0, lines-1
  ! To reduce roundoff error
  !   1) split summation into 4 terms
  term1 = 0d0; term2 = 0d0; term3 = 0d0; term4 = 0d0
  !   2) sum in reverse, since outer shells have less mass,
  !      so we want to add those together first
  do i = lines, n+1, -1
   term1 = term1 + rho(i-1)*r(i)  *(r(i)+r(i-1))
   term2 = term2 + rho(i)  *r(i-1)*(r(i)+r(i-1))
   term3 = term3 + rho(i)  *(r(i)**2+r(i)*r(i-1)+r(i-1)**2)
   term4 = term4 + rho(i-1)*(r(i)**2+r(i)*r(i-1)+r(i-1)**2)
  end do
  mnow = 0.5d0*(term1-term2)+(term3-term4)/3d0
  gpot(n) = -G*m(n)/(r(n)+1d-99)-4d0*pi*G*mnow
 end do

 do i = max(is_global,is-1), ie
  if(xi1(i)>=radius)then
   mc(i:ie) = mass
   exit
  end if
  do n = 0, lines-1
   if(r(n+1)>xi1(i).and.r(n)<=xi1(i))then
    ! Interpolate as m \propto r^3 (maybe not the best way?)
    mc(i) = intpol((r(n:n+1)/xi1(i))**3,m(n:n+1),1d0)
    exit
   end if
  end do
 end do

 allocate(Vshell(is:ie))
 do i = is_global, ie_global
  if (is <= i .and. i <= ie) then
   Vshell(i) = sum(dvol(i,js_global:je_global,ks_global:ke_global))
  end if
 enddo

 if(eq_sym)then
  volfac=2d0
 else
  volfac=1d0
 end if

!$omp parallel do collapse(3) private(i,j,k,n,nn,sn)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    d(i,j,k) = (mc(i)-mc(i-1))/(volfac*Vshell(i))
    if(x1(i)<r(1))then
     p(i,j,k) = pres(1)
     if(compswitch==2)then
      do nn = 1, spn-1
       do sn = 1, size(comp_list)-1
        if(trim(comp_list(sn))==trim(species(nn)))then
         spc(nn,i,j,k) = comp(sn,1)
         exit
        end if
       end do
      end do
      spc(spn,i,j,k) = 1d0-sum(spc(1:spn-1,i,j,k)) !dump the rest into others
     end if
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
!$omp end parallel do

!$omp parallel do private(i,j,k,n) collapse(3)
 do k = gks-2, gke+2
  do j = gjs-2, gje+2
   do i = gis-2, gie+2
    if(x1(i)<r(lines-1))then
     do n = 0, lines-1
      if(r(n+1)>abs(x1(i)).and.r(n)<=abs(x1(i)))then
       grvphi(i,j,k) = intpol(r(n:n+1),gpot(n:n+1),abs(x1(i)))
      end if
     end do
    else
     grvphi(i,j,k) = -G*mass/x1(i)
    end if
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine set_star_sph_grid

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE SET_STAR_CYL_GRID
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place star on a cylindrical coordinate grid

subroutine set_star_cyl_grid(r,m,pres,comp,comp_list)

 use settings,only:compswitch,spn
 use constants,only:G
 use grid
 use physval
 use utils,only:intpol

 real(8),allocatable,dimension(:),intent(in):: r,m,pres
 real(8),allocatable,dimension(:,:),intent(in),optional:: comp
 character(len=10),allocatable,intent(in),optional:: comp_list(:)
 real(8)::dr,mnow,rnow,mold,shell,shelld,mass,radius
 integer:: i,j,k,n,nn,sn
 real(8),allocatable:: rdis0(:,:)

!-----------------------------------------------------------------------------

 allocate(rdis0(is:ie,ks:ke))
 do i = is, ie
  do k = ks, ke
   rdis0(i,k) = sqrt( x1(i)**2+x3(k)**2 )
  end do
 end do

 mass=m(size(m)-1)
 radius=r(size(r)-1)
 dr = dx1(is)*1.d0
 rnow = dx1(is) ; mold = 0d0
 d = -1d0
 do
  shell = 0d0
  ! Choose cells in current mass shell
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(rdis0(i,k)<rnow.and.d(i,j,k)<=1d-99)then
      shell = shell + dvol(i,j,k)
     end if
    end do
   end do
  end do
  ! Evaluate the mass coordinate
  do n = 0, size(m)-2
   if(rnow>r(n).and.rnow<=r(n+1))then
    mnow = intpol(r(n:n+1),m(n:n+1),rnow)
    shelld = (mnow-mold)/shell
    exit
   end if
  end do
  ! Set density for chosen cells
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(rdis0(i,k)<=rnow.and.d(i,j,k)<=1d-99)then
      d(i,j,k) = shelld
      do n = 0, size(m)-2
       if(rdis0(i,k)>r(n).and.rdis0(i,k)<=r(n+1))then
        p(i,j,k) = intpol(r(n:n+1),pres(n:n+1),rdis0(i,k))

        if(compswitch==2)then
         do nn = 1, spn-1
          do sn = 1, size(comp_list)-1
           if(trim(comp_list(sn))==trim(species(nn)))then
            spc(nn,i,j,k) = intpol(r(n:n+1),comp(sn,n:n+1),rdis0(i,k))
            exit
           end if
          end do
         end do
         spc(spn,i,j,k) = 1d0-sum(spc(1:spn-1,i,j,k)) !dump the rest into others
        end if

        exit
       end if
      end do
     end if
    end do
   end do
  end do
  mold = mnow ; rnow = rnow + dr
  if(rnow>radius)then
   shell = 0d0
   do k = ks, ke
    do j = js, je
     do i = is, ie
      if(rdis0(i,k)<=rnow.and.d(i,j,k)<=1d-99)then
       shell = shell + dvol(i,js,k)
      end if
     end do
    end do
   end do
   do k = ks, ke
    do j = js, je
     do i = is, ie
      if(rdis0(i,k)<=rnow.and.d(i,j,k)<=1d-99)then
       d(i,j,k) = (mass-mold)/shell
       p(i,j,k) = G*mass*(mass-mold)/shell/rdis0(i,k)
       if(compswitch==2)then
        do nn = 1, spn-1
         do sn = 1, size(comp_list)-1
          if(trim(comp_list(sn))==trim(species(nn)))then
           spc(nn,i,j,k) = comp(sn,size(r)-1)
           exit
          end if
         end do
        end do
        spc(spn,i,j,k) = 1d0-sum(spc(1:spn-1,i,j,k)) !dump the rest into others
       end if
      end if
     end do
    end do
   end do
   exit
  end if
 end do

return
end subroutine set_star_cyl_grid


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                       SUBROUTINE ONE_SHOT_INWARDS
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculate a hydrostatic structure for a given entropy

subroutine one_shot_inwards(Sc,imu,r,mcore,msoft,rho,p,mass)

 use constants,only:G,pi
 use pressure_mod,only:get_d_from_ps
 use utils,only:softened_acc

 real(8),intent(in)::Sc,mcore,msoft
 real(8),intent(inout)::imu
 real(8),allocatable,dimension(:),intent(in)::r
 real(8),allocatable,dimension(:),intent(inout)::rho,p
 real(8),intent(out)::mass

 integer i,Nmax
 real(8)::hsoft
 real(8),allocatable,dimension(:)::dr,vol

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
end subroutine one_shot_inwards

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                       SUBROUTINE ISENTROPIC_STAR1
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculate a hydrostatic structure for a given entropy

subroutine isentropic_star1(Sc,imu,m,rsoft,r,rho,p)

 use constants,only:G,pi
 use pressure_mod,only:get_d_from_ps
 use utils,only:softened_acc

 real(8),intent(in)::Sc,rsoft
 real(8),intent(inout)::imu
 real(8),allocatable,dimension(:),intent(in)::m
 real(8),allocatable,dimension(:),intent(inout)::r,rho,p
 integer i,Nmax,which
 real(8):: fac
 real(8),parameter:: err=1d-8

!-----------------------------------------------------------------------------

 Nmax=size(rho)-1
 p = -1d0
 fac = 0.5d0

 p(0:1) = G*(m(Nmax)-m(0))**2/r(Nmax)**4/(8d0*pi) ! Initial guess
 which = 0

 convergence_loop:do
  rho(0:1) = get_d_from_ps(p(0),Sc,imu)
  r(1) = ((m(1)-m(0))/(4d0*pi/3d0*rho(1)))**(1d0/3d0)
  shot_loop:do i = 2, Nmax
   p(i) = p(i-1) - (m(i)-m(i-2))*G/(8d0*pi*r(i-1)**2) &
                 *((m(i-1)-m(0))/r(i-1)**2+m(0)*softened_acc(r(i-1),rsoft))
   if(p(i)<=0d0)exit shot_loop
   rho(i) = get_d_from_ps(p(i),Sc,imu)
   r(i) = (r(i-1)**3+3d0*(m(i)-m(i-1))/(4d0*pi*rho(i)))**(1d0/3d0)
  end do shot_loop

  if(p(Nmax)/(p(0)*1d-8)-1d0>err)then
   p(Nmax) = -1d0
   p(0:1)=p(0)*(1d0-fac)
   if(which>0)fac=fac*0.5d0
   which=-1
  elseif(p(Nmax)/(p(0)*1d-8)-1d0<-err)then
   p(0:1)=p(0)*(1d0+fac)
   if(which<0)fac=fac*0.5d0
   which=1
  else
   exit convergence_loop
  end if

 end do convergence_loop

! print'(4(1PE13.5e2))',r(i),m(i),rho(i),p(i)

return
end subroutine isentropic_star1

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                     SUBROUTINE GET_SOFTENED_PROFILE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Returns softened core profile with fixed entropy

subroutine get_softened_profile(r,mpt,mh,imuh,rho,p,ierr)

 use settings,only:eostype
 use pressure_mod,only:entropy_from_dp

 real(8),allocatable,dimension(:),intent(in)::r
 real(8),intent(in)::mh,imuh
 real(8),intent(inout)::mpt
 real(8),allocatable,dimension(:),intent(inout)::rho,p
 integer,intent(out)::ierr

 integer Nmax,i,eostype0
 real(8)::Sc,mass,mold,msoft,fac,Sedge,T,imu

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
  call one_shot_inwards(Sc,imu,r,mpt,msoft,rho,p,mass)
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


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE ISENTROPIC_STAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To make polytropic (isentropic) stars with arbitrary cores

subroutine isentropic_star(mass,radius,mcore,rsoft,imu,m,r,rho,p)

 use constants,only:G,pi
 use pressure_mod,only:eostype,entropy_from_dp
 use utils,only:geometrical_series
 use mpi_utils,only:myrank,allreduce_mpi

 real(8),intent(in)::mass,radius,mcore,rsoft
 real(8),intent(inout):: imu
 real(8),allocatable,dimension(:),intent(inout)::m,r,rho,p
 real(8):: Sc, fac, T, err, dmmin
 real(8),allocatable::dm(:)
 integer:: i,Nmax,eostype0,which

!-----------------------------------------------------------------------------

 Nmax = 3000
 allocate( m(0:Nmax),r(0:Nmax),rho(0:Nmax),p(0:Nmax) )
 m = 0d0
 r = 0d0
 rho = 0d0
 p = 0d0

 ! Run on master task and copy to other tasks
 if (myrank==0) then
 ! Can only do uniform composition for now
 ! This module does not work with non-ideal EoSs
  eostype0 = eostype
  if(eostype>=2) eostype = 1

  allocate( dm(-1:Nmax+2) )

  m(0) = mcore
  dmmin = (mass-mcore)/dble(Nmax)*1d-3
  call geometrical_series(dm,dmmin,1,Nmax/2,0d0,(mass-mcore)/2d0)
  do i = 1, Nmax/2
   m(i) = m(i-1) + dm(i)
  end do
  do i = Nmax/2+1, Nmax
   m(i) = m(i-1) + dm(Nmax-i+1)
  end do

  fac = 0.05d0
  err = 1d-5
  which=0

 ! Initial guess
  p(0) = G*(mass-mcore)**2/radius**4/(8d0*pi)
  rho(0) = (mass-mcore)/(4d0*pi/3d0*radius**3)
  T = 1d3
  Sc = entropy_from_dp(rho(0),p(0),T,imu)

  do
   r(Nmax) = radius
   call isentropic_star1(Sc,imu,m,rsoft,r,rho,p)

   if(r(Nmax)/radius-1d0>err)then
    Sc=Sc*(1d0-fac)
    if(which>0)fac=fac*0.5d0
    which=-1
   elseif(r(Nmax)/radius-1d0<-err)then
    Sc=Sc*(1d0+fac)
    if(which<0)fac=fac*0.5d0
    which=1
   else
    exit
   end if

  end do

  eostype = eostype0 ! set back to original eostype

 endif

 ! Broadcast to all tasks using allreduce
 call allreduce_mpi('sum', m)
 call allreduce_mpi('sum', r)
 call allreduce_mpi('sum', rho)
 call allreduce_mpi('sum', p)

return
end subroutine isentropic_star


end module star_mod
