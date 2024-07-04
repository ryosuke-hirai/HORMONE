module redsupergiant_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE REDSUPERGIANT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a red supergiant model at the centre of a spherical grid.

subroutine redsupergiant

 use settings,only:compswitch,spn,extrasfile,include_sinks
 use constants,only:G
 use grid
 use physval
 use input_mod
 use star_mod,only:isentropic_star,replace_core,set_star_sph_grid
 use gravmod,only:extgrv,mc
 use sink_mod,only:sink,sinkfield
 use utils,only:softened_pot
 use output_mod,only:write_extgrv
 use composition_mod,only:get_imu

 character(len=100):: mesafile
 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 character(len=10)::spc_list(1:1000)
 integer:: i,j,k,istat,nn,sn,ih1,ihe4
 real(8)::rcore,mcore,dbg,mass,spc_bg(1:spn),radius,imu_const
 real(8),allocatable:: comptmp(:)
 logical::isentropic,core_is_sink

!-----------------------------------------------------------------------------

 namelist /rsg_con/ mesafile,spc_list,rcore,isentropic,core_is_sink

 spc_list='aaa'
! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('rsg',extrasfile)
 read(nn,NML=rsg_con,iostat=istat)
 if(istat/=0)call error_nml('rsg',extrasfile)
 close(nn)

! Re-count spn based on spc_list and reallocate relevant arrays
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

! Read MESA file
 call read_mesa(mesafile,r,m,rho,pres,comp=comp,comp_list=comp_list)

 mass = m(size(m)-1)
 dbg = rho(size(rho)-1)*1d-5
! Use outermost composition as the ambient gas composition
 do nn = 1, spn-1
  do sn = 1, size(comp_list)
   if(trim(comp_list(sn))==trim(species(nn)))then
    spc_bg(nn) = comp(sn,size(rho)-1)
    exit
   end if
  end do
 end do
 spc_bg(spn) = 1d0-sum(spc_bg(1:spn-1))

! Replace the core with a point particle + softened gas
 call replace_core(rcore,r,m,rho,pres,comp,comp_list)

! Set external gravity
 if(core_is_sink)then
  if(.not.include_sinks)stop 'Set include_sinks=.true. if core_is_sink=.true.'
  sink(1)%mass = m(0)
  sink(1)%lsoft = rcore
  sink(1)%softfac = 3d0
  sink(1)%x(1:3) = 0d0
  sink(1)%v(1:3) = 0d0
  call sinkfield
 else
  do i = is, ie+2
   extgrv(i,js-2:je+2,ks-2:ke+2) = G*m(0)*softened_pot(x1(i),rcore)
  end do
  extgrv(is-1,js-2:je+2,ks-2:ke+2) =  extgrv(is,js-2:je+2,ks-2:ke+2)
 end if


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
     d(i,j,k) = dbg
     p(i,j,k) = G*mass*d(i,j,k)/x1(i)
     if(compswitch==2)spc(1:spn,i,j,k) = spc_bg(1:spn)
    end if
   end do
  end do
 end do

! Remember core mass
 if (is==is_global) mc(is-1) = mcore

 call write_extgrv

return
end subroutine redsupergiant

end module redsupergiant_mod
