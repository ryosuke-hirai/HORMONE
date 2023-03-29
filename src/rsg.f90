!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE REDSUPERGIANT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To place a red supergiant model at the centre of a spherical grid.

subroutine redsupergiant

 use settings,only:compswitch,spn
 use constants,only:G
 use grid
 use physval
 use input_mod
 use star_mod
 use gravmod,only:grvphi,extgrv
 use utils,only:softened_pot

 implicit none

 character*100:: mesafile
 real*8,allocatable,dimension(:):: r,m,rho,pres
 real*8,allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 integer nel,nn,sn
 real*8::rcore,dbg,mass,spc_bg(1:spn)
 
!-----------------------------------------------------------------------------

! Specify input file, elements you want to track, and a softening length
 mesafile='16lateRSG.data'
 species(1:spn) = [character(len=10):: &
                  'h1','he4','he3','c12','n14','o16','fe56','others']
 rcore = 1.5e12
 
 call read_mesa(mesafile,r,m,rho,pres,comp,comp_list)

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
 do i = is, ie+2
  extgrv(i,js-2:je+2,ks-2:ke+2) = G*m(0)*softened_pot(x1(i),rcore)
 end do
 extgrv(is-1,js-2:je+2,ks-2:ke+2) =  extgrv(is,js-2:je+2,ks-2:ke+2)
 
! Place the star at the origin
 m = m-m(0)
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
 
return
end subroutine redsupergiant
