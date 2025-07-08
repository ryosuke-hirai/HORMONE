module stellarcollision_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE STELLARCOLLISION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a stellar model at the centre of a spherical grid and a
!          point particle colliding into it.

subroutine stellarcollision

 use settings,only:compswitch,spn,extrasfile
 use constants,only:G,msun,rsun
 use grid
 use physval
 use input_mod
 use star_mod,only:isentropic_star,replace_core,set_star_sph_grid
 use gravmod,only:mc ,grvphi,totphi
 use sink_mod,only:sink,sinkfield
 use utils,only:softened_pot,polcar
 use composition_mod,only:get_imu
 use pressure_mod,only:eos_e,eos_p

 character(len=100):: mesafile
 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 character(len=10)::spc_list(1:1000)
 integer:: i,j,k,istat,nn,sn,ih1,ihe4
 real(8)::rcore,mcore,dbg,mass,spc_bg(1:spn),radius,imu_const,gradphi
 real(8)::nsmass,kickvel,kicktheta,kickphi,nssoft,asep,mprog,orbv
 real(8),allocatable:: comptmp(:)
 logical::isentropic

!-----------------------------------------------------------------------------

 namelist /col_con/ mesafile,spc_list,rcore,isentropic,&
                    nsmass,kickvel,kicktheta,kickphi,nssoft,asep,mprog

 spc_list='aaa'
! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('stellarcollision',extrasfile)
 read(nn,NML=col_con,iostat=istat)
 if(istat/=0)call error_nml('stellarcollision',extrasfile)
 close(nn)

 rcore = rcore*rsun
 nsmass = nsmass*msun
 mprog  = mprog *msun
 kickvel= kickvel*1d5
 nssoft = nssoft*rsun
 asep   = asep*rsun

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
 dbg = rho(size(rho)-1)*1d-5
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
 sink(1)%softfac = 3d0
 sink(1)%x(1:3) = 0d0
 sink(1)%v(1:3) = 0d0

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
 totphi = grvphi
 call sinkfield
 p(ie+1:ie+2,js:je,ks:ke) = G*mass/x1(ie+1)*d(ie,js,ks)
 do i = ie, is, -1
  gradphi = dx1(i+1)**2*totphi(i+2,js,ks)-dx1(i+2)**2*totphi(i,js,ks)+(dx1(i+2)**2-dx1(i+1)**2)*totphi(i+1,js,ks)
  p(i,js,ks) = (d(i+1,js,ks)*gradphi+dx1(i+1)**2*p(i+2,js,ks)+(dx1(i+2)**2-dx1(i+1)**2)*p(i+1,js,ks))/dx1(i+2)**2
  p(i,js:je,ks:ke) = p(i,js,ks)
 end do

! Remember core mass
 if (is==is_global) mc(is-1) = mcore

! place colliding object
 sink(2)%mass = nsmass
 sink(2)%x(1) = -asep
 sink(2)%x(2) = 0d0
 sink(2)%x(3) = 0d0
 orbv = sqrt(G*(mass+mprog)/asep)
 sink(2)%v(1) = kickvel*sin(kicktheta)*cos(kickphi)
 sink(2)%v(2) = kickvel*cos(kicktheta)+orbv
 sink(2)%v(3) = kickvel*sin(kicktheta)*sin(kickphi)
 sink(2)%v(2) = norm2(sink(2)%v(2:3))
 sink(2)%v(3) = 0d0
 sink(2)%lsoft = nssoft
 sink(2)%softfac = 3d0

return
end subroutine stellarcollision

end module stellarcollision_mod
