module sn2022jli_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SN2022JLI
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a stellar model at the centre of a spherical grid and a
!          point particle colliding into it. The frame is fixed to the centre
!          of the star and we also inject heat to the envelope due to ejecta-
!          companion interaction.

subroutine sn2022jli

 use settings,only:compswitch,spn,extrasfile,eostype,eq_sym
 use constants,only:G,msun,rsun,pi,tiny
 use grid
 use physval
 use input_mod
 use star_mod,only:isentropic_star,replace_core,set_star_sph_grid
 use gravmod,only:mc ,grvphi,totphi
 use sink_mod,only:sink,sinkfield
 use utils,only:softened_pot,polcar,get_vpol,gravpot1d
 use composition_mod,only:get_imu
 use pressure_mod,only:eos_e,eos_p,get_e_from_ds,entropy_from_dT
 use mpi_utils,only:allreduce_mpi

 character(len=100):: mesafile
 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 character(len=10)::spc_list(1:1000)
 integer:: i,j,k,istat,nn,sn,ih1,ihe4,which
 real(8)::rcore,mcore,dbg,mass,spc_bg(1:spn),radius,imu_const,gradphi
 real(8)::nsmass,nssoft,asep,dis,Porb,mprog,orbv,ecc
 real(8),allocatable:: comptmp(:),p1d(:)
 logical::isentropic
 real(8):: Eexp,Ebind,Ebind0,entr,entr0,mheat,Omega,Eheat,fac,dfac,TT,dnow,phinow(3),newphi
 real(8),parameter:: err=1d-8

!-----------------------------------------------------------------------------

 namelist /jli_con/ mesafile,spc_list,rcore,isentropic,&
                    nsmass,nssoft,dis,mprog,Porb,Eexp

 spc_list='aaa'
! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('sn2022jli',extrasfile)
 read(nn,NML=jli_con,iostat=istat)
 if(istat/=0)call error_nml('sn2022jli',extrasfile)
 close(nn)

 rcore = rcore*rsun
 nsmass = nsmass*msun
 mprog  = mprog *msun
 nssoft = nssoft*rsun

! Re-count spn based on spc_list and reallocate relevant arrays
 if(compswitch==2)then
  spn=-1
  do nn = 1, 1000
   spn = spn + 1
   if(spc_list(nn)=='aaa')exit
  end do
  deallocate(spc,spcorg,dspc,spcflx,species)
  allocate(spc    (1:spn,is-2:ie+2,js-2:je+2,ks-2:ke+2), &
           spcorg (1:spn,is:ie,js:je,ks:ke), &
           dspc   (1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
           spcflx (1:spn,is-1:ie+1,js-1:je+1,ks-1:ke+1,1:3), &
           species(1:spn) )
  species(1:spn) = spc_list(1:spn)
 end if

! Read MESA file
 call read_mesa(mesafile,r,m,rho,pres,comp=comp,comp_list=comp_list)

 mass = m(size(m)-1)
 dbg = rho(size(rho)-1)*1d0
! Use outermost composition as the ambient gas composition
 if(compswitch==2)then
  do nn = 1, spn-1
   do sn = 1, size(comp_list)
    if(trim(comp_list(sn))==trim(species(nn)))then
     spc_bg(nn) = comp(sn,size(rho)-1)
     exit
    end if
   end do
  end do
  spc_bg(spn) = 1d0-sum(spc_bg(1:spn-1))
 end if

! Replace the core with a point particle + softened gas
 call replace_core(rcore,r,m,rho,pres,comp,comp_list)

! Set external gravity
 sink(1)%mass = m(0)
 sink(1)%lsoft = rcore
 sink(1)%laccr = 0d0
 sink(1)%softfac = 3d0
 sink(1)%x(1:3) = 0d0
 sink(1)%v(1:3) = 0d0
 sink(1)%mdot = 0d0
 sink(1)%jdot = 0d0
 sink(1)%Jspin = 0d0

 mcore = m(0)
 mass  = m(size(m)-1)
 radius= r(size(r)-1)

! Isentropic envelope
 if(isentropic)then
! get indices for hydrogen and helium
  do i = 1, size(comp_list)
   if(trim(comp_list(i))=='h1')ih1=i
   if(trim(comp_list(i))=='he4')ihe4=i
  end do
! remember central composition
  allocate(comptmp(1:size(comp_list)))
  comptmp = comp(1:size(comp_list),1)
  imu_const = get_imu((/comp(ih1,1),comp(ihe4,1)/))

  deallocate(m,r,rho,pres,comp)
  call isentropic_star(mass,radius,mcore,rcore,imu_const,m,r,rho,pres)
! re-insert composition
  allocate(comp(1:size(comp_list),0:size(m)-1))
  do i = 0, size(m)-1
   comp(:,i) = comptmp
  end do
 end if

! Place the star at the origin
 m = m-mcore
 call set_star_sph_grid(r,m,rho,pres,comp,comp_list)

! Attach a uniform density atmosphere
 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(d(i,j,k)<0d0)then
     d(i,j,k) = dbg*(radius/x1(i))**2
     p(i,j,k) = G*mass*d(i,j,k)/x1(i)
     if(compswitch==2)spc(1:spn,i,j,k) = spc_bg(1:spn)
    end if
   end do
  end do
 end do

! Re-solve hydrostatic equilibrium
 call gravpot1d
 totphi = grvphi
 call sinkfield
 allocate(p1d(is_global-1:ie_global+2))
 p1d(ie_global+1) = G*mc(ie_global)/x1(ie_global+1)*dbg*(radius/x1(ie_global+1))**2
 p1d(ie_global) = G*mc(ie_global)/x1(ie_global)*dbg*(radius/x1(ie_global))**2
 if(ie==ie_global)phinow(1:3) = totphi(ie:ie+2,js,ks)

 do i = ie_global-1, is_global, -1
  phinow(2:3) = phinow(1:2)
  dnow = 0d0 ; newphi = 0d0
  if(i  >=is.and.i  <=ie)newphi = totphi(i,js,ks)
  if(i+1>=is.and.i+1<=ie)dnow   = d(i+1,js,ks)

  call allreduce_mpi('max',dnow)
  call allreduce_mpi('min',newphi)
  phinow(1) = newphi
  gradphi = dx1(i+1)**2*phinow(3)-dx1(i+2)**2*phinow(1)+(dx1(i+2)**2-dx1(i+1)**2)*phinow(2)
  p1d(i) = (dnow*gradphi+dx1(i+1)**2*p1d(i+2)+(dx1(i+2)**2-dx1(i+1)**2)*p1d(i+1))/dx1(i+2)**2
  if(i>=is.and.i<=ie)p(i,js:je,ks:ke) = p1d(i)
 end do

! Remember core mass
 if (is==is_global) mc(is-1) = mcore

! place colliding object
 sink(2)%mass = nsmass
 sink(2)%lsoft = nssoft
 sink(2)%softfac = 3d0
 sink(2)%laccr = nssoft
 sink(2)%mdot = 0d0
 sink(2)%racc = 0d0
 sink(2)%facc = 0d0
 sink(2)%jdot = 0d0
 sink(2)%Jspin = 0d0
 sink(2)%jet_ang = 45d0
 sink(2)%jet_dir = [sin(0.25d0*pi),0d0,cos(0.25d0*pi)]

! For SN2022jli
 asep = (G*(mass+nsmass)*(Porb*3600d0*24d0/2d0/pi)**2)**(1d0/3d0)
 orbv = sqrt(G*(mass+nsmass)/asep)
 ecc = 1d0-dis*rsun/asep

 sink(2)%x(1)=-asep*(1d0-ecc)
 sink(2)%x(2:3)=0d0
 sink(2)%v(1)=0d0
 sink(2)%v(2)=orbv*sqrt((1d0+ecc)/(1d0-ecc))*mass/(mass+nsmass)
 sink(2)%v(3)=0d0

 sink(1)%v(2)=-orbv*sqrt((1d0+ecc)/(1d0-ecc))*nsmass/(mass+nsmass)

 sink(2)%v = sink(2)%v - sink(1)%v
 sink(1)%v = 0d0

! Try to make NS atmosphere hydrostatic
!$omp parallel do private(i,j,k) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    select case(eostype)
    case(0,1)
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2)
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
    eint(i,j,k) = eint(i,j,k) &
                - G*nsmass*d(i,j,k)&
                  *softened_pot(norm2(car_x(:,i,j,k)-sink(2)%x),sink(2)%lsoft)
    select case(eostype)
    case(0,1)
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))
    case(2)
     p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k),&
                      spc(1,i,j,k),spc(2,i,j,k))
    end select
   end do
  end do
 end do
!$omp end parallel do

 Ebind0 = 0d0
 fac = 1d0; if(eq_sym)fac = 2d0
 !$omp parallel do private(i,j,k,entr) collapse(3) reduction(+:Ebind0)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(x1(i)<=radius)then
     Ebind0 = Ebind0 + (eint(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)*fac
    end if
   end do
  end do
 end do
!$omp end parallel do
 call allreduce_mpi('sum',Ebind0)

 Omega = (1d0-sqrt(1d0-(radius/(asep*(1d0-ecc)))**2))/2d0
 mheat = 3d0*msun*Omega/2d0
 Eheat = Eexp*Omega/12d0

 entr0 = 5d0
 dfac = 0.2d0
 which=0

! Find the right entropy normalization through the bisection method
 do

  Ebind = 0d0
!$omp parallel do private(i,j,k,entr,TT) collapse(3) reduction(+:Ebind)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(x1(i)<=radius)then

      entr = entropy_from_dT(d(i,j,k),T(i,j,k),imu(i,j,k),spc(1,i,j,k),spc(2,i,j,k))

      entr = entr + entr0*min(1d0,mheat/((mass-mc(i))+tiny))
      eint(i,j,k) = get_e_from_ds(d(i,j,k),entr,imu(i,j,k),spc(1,i,j,k),spc(2,i,j,k))
      TT = T(i,j,k)
      select case(eostype)
      case(0,1)
       p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),TT,imu(i,j,k))
      case(2)
       p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),TT,imu(i,j,k),&
                       spc(1,i,j,k),spc(2,i,j,k))
      end select
      Ebind = Ebind + (eint(i,j,k)+0.5d0*d(i,j,k)*grvphi(i,j,k))*dvol(i,j,k)*fac
     end if
    end do
   end do
  end do
!$omp end parallel do
  call allreduce_mpi('sum',Ebind)

! bisection method
  if((Ebind-Ebind0-Eheat)/abs(Ebind0)>err)then
   entr0 = entr0*(1d0-dfac)
   if(which>0)dfac=dfac*0.5d0
   which=-1
  elseif((Ebind-Ebind0-Eheat)/abs(Ebind0)<-err)then
   entr0 = entr0*(1d0+dfac)
   if(which<0)dfac=dfac*0.5d0
   which=1
  else
   exit
  end if

 end do

return
end subroutine sn2022jli

end module sn2022jli_mod
