module eruption_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE ERUPTION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To simulate eruptions of stars (Owocki+2019)

subroutine eruption

 use settings,only:compswitch,spn,extrasfile,is_test,eostype
 use constants,only:G,rsun
 use grid
 use physval
 use input_mod
 use star_mod
 use pressure_mod,only:eos_e,pressure
 use composition_mod,only:meanmolweight
 use gravmod,only:mc

 character(len=100):: mesafile
 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 character(len=10)::spc_list(1:1000)
 integer::i,j,k,nn,sn,istat,i_inj_in,i_inj_out
 real(8)::dbg,mass,radius,spc_bg(1:spn),Eexp,Ebind,ene,Temp,&
          inj_in,inj_out,imu_local,ecell

!-----------------------------------------------------------------------------

 namelist /erupcon/ mesafile,spc_list,Eexp,inj_in,inj_out

 spc_list='aaa'
 if(is_test) extrasfile='../para/extras_eruption'

! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('eruption',extrasfile)
 read(nn,NML=erupcon,iostat=istat)
 if(istat/=0)call error_nml('eruption',extrasfile)
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

 ! Not all of the species in the list may be in the mesa file
 ! so they need to initialised to zero to avoid errors
 spc = 0d0
 spc_bg = 0d0

 select case (compswitch)
 case(1)
  call read_mesa(mesafile,r,m,rho,pres,mu=mudata)
 case(2)
  call read_mesa(mesafile,r,m,rho,pres,comp=comp,comp_list=comp_list)
 case default
  call read_mesa(mesafile,r,m,rho,pres)
 end select

 mass = m(size(m)-1)
 radius = r(size(r)-1)
 dbg = rho(size(rho)-1)*1d-3
 musize = size(m)-1

! Compute binding energy
 Ebind = 0d0
 do i = 2, size(rho)-1
  ! Get local internal energy
  select case (eostype)
  case(1)
   imu_local = 1d0/mudata(i,1)
   ene = eos_e(rho(i),pres(i),Temp,imu_local)
  case(2)
   ene = eos_e(rho(i),pres(i),Temp,imu_local,X=comp(1,i),Y=comp(2,i))
  case default
   ene = eos_e(rho(i),pres(i),Temp,imu_local)
  end select
  Ebind = Ebind - G*m(i-1)*(m(i)-m(i-1)) / r(i-1) + ene*(m(i)-m(i-1))
 end do

! Set explosion energy to fraction of binding energy
 Eexp = - Ebind * Eexp

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

! Place the star at the origin
 call set_star_sph_grid(r,m,pres,comp,comp_list)

! Attach a wind-like atmosphere
 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(d(i,j,k)<0d0)then
     d(i,j,k) = dbg*(radius/x1(i))**2
     if(compswitch==2)spc(1:spn,i,j,k) = spc_bg(1:spn)
    end if
   end do
  end do
 end do

! Re-calculate hydrostatic equilibrium
 call meanmolweight
 do k = ks, ke
  do j = js, je
   do i = ie, is, -1
    p(i,j,k) = p(i+1,j,k)+G*mc(i)*d(i,j,k)/xi1(i)**2*dx1(i)
    select case (eostype)
    case(0,1)
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    case(2)
     eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k),&
                         spc(1,i,j,k),spc(2,i,j,k))
    end select
   end do
  end do
 end do

! Find inner and outer boundary
 i_inj_in = is-1
 do i = is, ie
  if(xi1(i)>inj_in*rsun.and.i_inj_in<is)then
   i_inj_in = i
  end if
  if(xi1(i)>inj_out*rsun.or.xi1(i)>r(size(r)-1))then
   i_inj_out = i
   exit
  end if
 end do

 ecell = Eexp / (mc(i_inj_out)-mc(i_inj_in-1))

 do k = ks, ke
  do j = js, je
   do i = i_inj_in, i_inj_out
    eint(i,j,k) = eint(i,j,k) + ecell*d(i,j,k)
   end do
  end do
 end do
 e = eint

 call pressure

return
end subroutine eruption

end module eruption_mod
