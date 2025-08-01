module windtunnel_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE WINDTUNNEL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up wind tunnel simulation

subroutine windtunnel

 use constants,only:msun,rsun,Rgas,arad
 use settings,only:extrasfile
 use grid
 use physval
 use input_mod
 use star_mod
 use fluxbound_mod,only:vwind,dwind,pwind,Twind

 real(8),allocatable,dimension(:):: r,m,rho,pres
 character(len=100)::star_type,mesafile
 real(8)::mass,radius,poly_n,gamma0,imu_const
 real(8),allocatable,dimension(:,:):: comp
 character(len=10),allocatable:: comp_list(:)
 integer::i,j,k,nn,istat

!-----------------------------------------------------------------------------

 namelist /wtnlcon/ star_type,mass,radius,poly_n,mesafile,vwind,dwind,Twind

! Specify input file, elements you want to track, and a softening length
 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('windtunnel',extrasfile)
 read(nn,NML=wtnlcon,iostat=istat)
 if(istat/=0)call error_nml('windtunnel',extrasfile)
 close(nn)

 mass = mass*msun
 radius = radius*rsun

! Create polytrope
 gamma0=gamma
 gamma=1d0+1d0/poly_n
 imu_const=1d0/muconst
 call isentropic_star(mass,radius,0d0,0d0,imu_const,m,r,rho,pres)
 gamma=gamma0

! Set passive scalars
 allocate(comp_list(1:2),comp(1:2,0:size(m)))
 comp_list(1) = 'planet'
 comp_list(2) = 'wind'
 comp(1,:) = 1d0
 comp(2,:) = 0d0
 species=comp_list

! Place the star at the origin
 call set_star_cyl_grid(r,m,pres,comp,comp_list)

 pwind = (Rgas/muconst*dwind + arad*Twind**3/3d0*0d0)*Twind
! Embed the star in a uniform wind tunnel
 do k = ks, ke
  do j = js, je
   do i = is, ie
    if(d(i,j,k)<=0d0)then
     if(x3(k)<-1.05d0*radius)then
      d(i,j,k) = dwind*min(-x3(k)/radius-1d0,1d0)
      p(i,j,k) = pwind*min(-x3(k)/radius-1d0,1d0)
      v3(i,j,k) = vwind
     else
      d(i,j,k) = dwind*1d-2
      p(i,j,k) = pwind*1d-2
      v3(i,j,k) = 0d0
     end if
      v1(i,j,k) = 0d0
      v2(i,j,k) = 0d0
      spc(1,i,j,k) = 0d0
      spc(2,i,j,k) = 1d0
    end if
   end do
  end do
 end do

return
end subroutine windtunnel

end module windtunnel_mod
