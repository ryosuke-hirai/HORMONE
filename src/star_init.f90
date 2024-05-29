module star_init_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE STAR_INIT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a red supergiant model at the centre of a spherical grid.

subroutine star_init

 use settings,only:compswitch,spn,extrasfile,is_test
 use constants,only:G
 use grid
 use physval
 use input_mod
 use star_mod

 character(len=100):: mesafile
 real(8),allocatable,dimension(:):: r,m,rho,pres
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 character(len=10)::spc_list(1:1000)
 integer::i,j,k,nn,sn,istat
 real(8)::dbg,mass,radius,spc_bg(1:spn)

!-----------------------------------------------------------------------------

 namelist /starcon/ mesafile,spc_list

 spc_list='aaa'
 if(is_test) extrasfile='../para/extras_star'

! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('star',extrasfile)
 read(nn,NML=starcon,iostat=istat)
 if(istat/=0)call error_nml('star',extrasfile)
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

 call read_mesa(mesafile,r,m,rho,pres,comp,comp_list)

 mass = m(size(m)-1)
 radius = r(size(r)-1)
 dbg = rho(size(rho)-1)*1d-0

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
 call set_star_sph_grid(r,m,rho,pres,comp,comp_list)

! Attach a wind-like atmosphere
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

return
end subroutine star_init

end module star_init_mod
