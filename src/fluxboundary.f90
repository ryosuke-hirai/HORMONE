module fluxbound_mod
 implicit none

 real(8):: dwind,vwind,Twind,pwind,ewind
 logical:: read_wind = .false.

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE FLUXBOUNDARY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set flux boundary

subroutine fluxboundary

 use settings

!-----------------------------------------------------------------------------

!!$ if(bc1is==10)call fluxboundary1i
!!$ if(bc1os==10)call fluxboundary1o
!!$ if(bc2is==10)call fluxboundary2i
!!$ if(bc2os==10)call fluxboundary2o
 if(bc3is==10)call fluxboundary3i
!!$ if(bc3os==10)call fluxboundary3o

return
end subroutine fluxboundary


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE FLUXBOUNDARY3I
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set flux for x3 inner boundary

subroutine fluxboundary3i

 use settings,only:compswitch,spn,extrasfile
 use constants,only:arad,fac_pgas,fac_egas
 use grid
 use physval,only:flux3,spcflx,muconst,imo3,iene,icnt,spc,dspc
 use input_mod,only:error_extras,error_nml

 integer:: i,j,k,n,istat
 real(8):: poly_n,mass,radius
 character(len=100):: mesafile,star_type

!-----------------------------------------------------------------------------
 
! Read wind parameters
 namelist /wtnlcon/ star_type,mass,radius,poly_n,mesafile,vwind,dwind,Twind

 if(.not.read_wind)then
! Specify input file, elements you want to track, and a softening length
  open(newunit=n,file=extrasfile,status='old',iostat=istat)
  if(istat/=0)call error_extras('windtunnel',extrasfile)
  read(n,NML=wtnlcon,iostat=istat)
  if(istat/=0)call error_nml('windtunnel',extrasfile)
  close(n)
  ewind = fac_egas/muconst*dwind*Twind + 0.5d0*dwind*vwind**2
  pwind = fac_pgas/muconst*dwind*Twind
  read_wind = .true.
 end if

 k = ks-1
!$omp parallel do private(i,j,n) collapse(2)
 do j = js, je
  do i = is, ie
   flux3(i,j,k,:) = 0d0
   flux3(i,j,k,icnt) = dwind*vwind
   flux3(i,j,k,imo3) = dwind*vwind**2 + pwind
   flux3(i,j,k,iene) = (ewind+pwind) * vwind
   if(compswitch>=2)then
    spcflx(1,i,j,k,3) = 0d0
    spcflx(2,i,j,k,3) = flux3(i,j,k,1)
   end if
  end do
 end do
!$omp end parallel do

return
end subroutine fluxboundary3i

end module fluxbound_mod
