module commonenvelope_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE COMMONENVELOPE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a stellar model at the centre of a spherical grid and a
!          point particle orbiting it. The frame is fixed to the centre
!          of the star.

subroutine commonenvelope

 use settings,only:compswitch,spn,extrasfile,eostype
 use constants,only:G,msun,rsun,pi
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
 integer:: i,j,k,istat,nn,sn,ih1,ihe4
 real(8)::rcore,mcore,dbg,mass,spc_bg(1:spn),radius,imu_const,gradphi
 real(8)::compmass,compsoft,comprad,asep,dis,Porb,orbv,ecc
 real(8),allocatable:: comptmp(:),p1d(:)
 logical::isentropic
 real(8):: dnow,phinow(3),newphi,comp_jet_ang

!-----------------------------------------------------------------------------

 namelist /cee_con/ mesafile,spc_list,mass,radius,mcore,rcore,isentropic,&
                    compmass,compsoft,comprad,dis,Porb,ecc,comp_jet_ang

 spc_list='aaa'
! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('commonenvelope',extrasfile)
 read(nn,NML=cee_con,iostat=istat)
 if(istat/=0)call error_nml('commonenvelope',extrasfile)
 close(nn)

 rcore = rcore*rsun
 compmass = compmass*msun
 compsoft = compsoft*rsun
 comprad  = comprad *1d5

 if(trim(mesafile)/='')then
! If using MESA file

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
  dbg = rho(size(rho)-1)*1d-3
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

 else
! Create isentropic star if no MESA file is provided
  mass=mass*msun
  radius=radius*rsun
  mcore = mcore*msun
  imu_const = 1d0/muconst
  call isentropic_star(mass,radius,mcore,rcore,imu_const,m,r,rho,pres)
  dbg = rho(size(rho)-1)*1d-3

 end if

! Set external gravity
 sink(1)%mass = mcore
 sink(1)%lsoft = rcore
 sink(1)%laccr = 0d0
 sink(1)%softfac = 3d0
 sink(1)%x(1:3) = 0d0
 sink(1)%v(1:3) = 0d0
 sink(1)%mdot = 0d0
 sink(1)%jdot = 0d0
 sink(1)%Jspin = 0d0

! Place the star at the origin
 m = m-mcore
 call set_star_sph_grid(r,m,pres,comp,comp_list)

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
 sink(2)%mass = compmass
 sink(2)%lsoft = compsoft
 sink(2)%softfac = 3d0
 sink(2)%laccr = compsoft
 sink(2)%mdot = 0d0
 sink(2)%racc = comprad
 sink(2)%facc = 0d0
 sink(2)%jdot = 0d0
 sink(2)%Jspin = 0d0
 sink(2)%jet_ang = comp_jet_ang
 sink(2)%jet_dir = [0d0,0d0,1d0]

! Set orbit
 if(dis>0d0)then
  asep = dis
 else
  asep = (G*(mass+compmass)*(Porb*3600d0*24d0/2d0/pi)**2)**(1d0/3d0)  
 end if
 orbv = sqrt(G*(mass+compmass)/asep)

 sink(2)%x(1)=-asep*(1d0-ecc)
 sink(2)%x(2:3)=0d0
 sink(2)%v(1)=0d0
 sink(2)%v(2)=orbv*sqrt((1d0+ecc)/(1d0-ecc))*mass/(mass+compmass)
 sink(2)%v(3)=0d0

 sink(1)%v(2)=-orbv*sqrt((1d0+ecc)/(1d0-ecc))*compmass/(mass+compmass)

 sink(2)%v = sink(2)%v - sink(1)%v
 sink(1)%v = 0d0

! Try to make companion atmosphere hydrostatic
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
                - G*compmass*d(i,j,k)&
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

return
end subroutine commonenvelope

end module commonenvelope_mod
