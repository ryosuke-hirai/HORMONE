module tools_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                               SUBROUTINE TOOLS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Setting frequently used variables

subroutine tools

 use settings,only:crdnt,gravswitch,radswitch,include_cooling,eostype
 use grid,only:coscyl,gis,gie,gks,gke,cosc,is,ie,js,je,ks,ke,fmr_max,fmr_lvl
 use physval,only:gamma,imu,muconst
 use gravmod,only:llmax,Plc,Pl,dtg_unit
 use constants,only:Cv,Rgas
 use timestep_mod,only:dtgrav_cell
 use ionization_mod,only:ionization_setup
 use cooling_mod,only:cooling_setup
 use gravity_mod,only:gravsetup
 use radiation_mod,only:radiation_setup
 use smear_mod,only:smear_setup,block_j,block_k,lijk_from_id,nsmear
 use mpi_utils,only:allreduce_mpi
 use mpi_domain,only:partially_my_domain

 integer:: i,j,k,n,l,ll,jb,kb
 real(8),allocatable,dimension(:,:,:):: dtg

!-----------------------------------------------------------------------------

! ************************* Legendre polynomials *****************************

 select case(crdnt)
! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(1)
  allocate( Plc(0:llmax,gis:gie+2,gks-2:gke+2) )
  do k = gks-2, gke+2
   do i = gis, gie+2
    Plc(0,i,k) = 1d0
    Plc(1,i,k) = coscyl(i,k)
    Plc(2,i,k) = (3d0*coscyl(i,k)*Plc(1,i,k) - Plc(0,i,k)) * 0.5d0

    do ll = 3, llmax
     Plc(ll,i,k) = (dble(2*ll-1)*coscyl(i,k)*Plc(ll-1,i,k) &
                 -  dble(ll-1)              *Plc(ll-2,i,k)) /dble(ll)
    end do
   end do
  end do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 case(2)
  allocate( Pl(0:llmax,js-1:je+1) )
  do j=js-1,je+1
   Pl(0,j) = 1d0
   Pl(1,j) = cosc(j)
   Pl(2,j) = (3d0*cosc(j)*Pl(1,j) - Pl(0,j)) * 0.5d0

   do ll=3,llmax
    Pl(ll,j)   = (dble(2*ll-1)*cosc(j)*Pl(ll-1,j) &
               -  dble(ll-1)          *Pl(ll-2,j)) /dble(ll)
   end do
  end do

 end select

! reading composition data from datafile ! --------------------------------
!!$ open(unit=40,file='17lateRSG.data',status='old')
!!$ read(40,'()')
!!$ read(40,'()')
!!$ read(40,*) lines, lines
!!$ read(40,'()')
!!$ read(40,'()')
!!$ read(40,'(a)') dumc
!!$
!!$! counting rows
!!$ allocate(dum(500)) ; dum = 'aaa'
!!$ read(dumc,*,end=101) dum
!!$101 do i = 1, 500
!!$  if(dum(i)=='aaa')then
!!$   rows = i - 1
!!$   exit
!!$  end if
!!$ end do
!!$
!!$ allocate(header(1:rows),dat(1:lines,1:rows))
!!$ header(1:rows) = dum(1:rows)
!!$ deallocate(dum)
!!$
!!$ do i = 1, lines
!!$  read(40,*) dat(lines-i+1,1:rows)
!!$ end do
!!$
!!$ allocate(comp_ej(0:8,1:lines))
!!$ do i = 1, rows
!!$  if(trim(header(i))=='mass') comp_ej(0,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='h1'  ) comp_ej(1,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='he4' ) comp_ej(2,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='he3' ) comp_ej(3,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='c12' ) comp_ej(4,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='n14' ) comp_ej(5,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='o16' ) comp_ej(6,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='ne20') comp_ej(7,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='mg24') comp_ej(8,1:lines) = dat(1:lines,i)
!!$ end do
!!$
!!$! reduce unnecessarily high resolution
!!$ i = 2; compsize = lines
!!$ do
!!$  if(abs(comp_ej(1,i)-comp_ej(1,i+1))<1d-3*comp_ej(1,i).and.&
!!$     abs(comp_ej(2,i)-comp_ej(2,i+1))<1d-3*comp_ej(2,i).and.&
!!$     abs(comp_ej(3,i)-comp_ej(3,i+1))<1d-3*comp_ej(3,i).and.&
!!$     abs(comp_ej(4,i)-comp_ej(4,i+1))<1d-3*comp_ej(4,i).and.&
!!$     abs(comp_ej(5,i)-comp_ej(5,i+1))<1d-3*comp_ej(5,i).and.&
!!$     abs(comp_ej(6,i)-comp_ej(6,i+1))<1d-3*comp_ej(6,i).and.&
!!$     abs(comp_ej(7,i)-comp_ej(7,i+1))<1d-3*comp_ej(7,i).and.&
!!$     abs(comp_ej(8,i)-comp_ej(8,i+1))<1d-3*comp_ej(8,i))then
!!$   comp_ej(0:8,i:compsize-1) = comp_ej(0:8,i+1:compsize)
!!$   compsize = compsize-1
!!$   if(compsize<=i+1)exit
!!$   cycle
!!$  end if
!!$  i = i+1
!!$ end do
!!$ deallocate(dat)
!!$ allocate(dat(0:8,1:compsize))
!!$ dat(0:8,1:compsize) = comp_ej(0:8,1:compsize)
!!$ deallocate(comp_ej)
!!$ allocate(comp_ej(0:8,1:compsize))
!!$ comp_ej(0:8,1:compsize) = dat(0:8,compsize:1:-1)

! EoS parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Cv = Rgas/(gamma-1d0) ! frequently used factor for egas
 imu(is-2:ie+2,js-2:je+2,ks-2:ke+2) = 1d0/muconst

 if(eostype==2) call ionization_setup

! Set cooling parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(include_cooling) call cooling_setup

! Set gravity parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(gravswitch>0) call gravsetup

! Set radiation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(radswitch>0) call radiation_setup

! Set smearing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(fmr_max>0) call smear_setup

! Set baseline dtgrav %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(gravswitch==3)then
  allocate(dtg(is:ie,js:je,ks:ke))
!$omp parallel do private(i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     call dtgrav_cell(i,j,k,dtg,1d0)
    end do
   end do
  end do
!$omp end parallel do

! Use longer time step if using nested grids
  if(fmr_max>0.and.crdnt==2)then
!$omp parallel do private(i,j,k,l,jb,kb) schedule(dynamic)
   do n = 1, nsmear
    l = lijk_from_id(0,n)
    if(fmr_lvl(l)==0)cycle
    i = lijk_from_id(1,n)
    j = lijk_from_id(2,n)
    k = lijk_from_id(3,n)
    jb = block_j(l)
    kb = block_k(l)
    if(partially_my_domain(i,j,k,1,jb,kb)) &
                call dtgrav_cell(i,j,k,dtg,1d0,jb=jb,kb=kb)
   end do
!$omp end parallel do
  end if
  dtg_unit = minval(dtg)
  call allreduce_mpi('min',dtg_unit)
 end if

return
end subroutine tools

end module tools_mod
