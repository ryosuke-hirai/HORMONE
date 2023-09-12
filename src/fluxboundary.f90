module fluxbound_mod
 implicit none

 real(8):: dwind,vwind,Twind,pwind

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

 use settings,only:compswitch,spn
 use constants,only:arad,fac_pgas,fac_egas
 use grid
 use physval,only:flux3,spcflx,muconst,imo3,iene,icnt,spc,dspc

 integer:: i,j,k,n
 real(8)::fix,signdflx,ewind,dx(1:2)
 real(8),dimension(1:spn):: spcl,spcr

!-----------------------------------------------------------------------------

! wind parameters
! dwind = 0.04d-2 ! g/cc
! vwind = 230d5 ! cm/s
! Twind = 1.2d6 ! K
 ewind = fac_egas/muconst*dwind*Twind
 pwind = (fac_pgas/muconst*dwind + arad*Twind**3/3d0*0d0)*Twind

 k = ks-1
 dx(1) = xi3(k)-x3(k) ; dx(2) = x3(k+1)-xi3(k)
 do i = is, ie
  do j = js, je
   flux3(i,j,k,:) = 0d0
   flux3(i,j,k,icnt) = dwind*vwind
   flux3(i,j,k,imo3) = dwind*vwind**2 + pwind
   flux3(i,j,k,iene) = (ewind+0.5d0*dwind*vwind**2+pwind) * vwind
   if(compswitch>=2)then
!!$    do n = 1, spn
!!$     spcl(n) = spc(n,i  ,j,k) + dx(1) * dspc(n,i  ,j,k,1)
!!$     spcr(n) = spc(n,i+1,j,k) - dx(2) * dspc(n,i+1,j,k,1)
!!$    end do
!!$    signdflx = sign(0.5d0,flux3(i,j,k,1))
!!$    fix = 1d0/( sum(spcl(1:spn))*(0.5d0+signdflx) &
!!$              + sum(spcr(1:spn))*(0.5d0-signdflx) )
!!$    do n = 1, spn
!!$     spcflx(n,i,j,k,1) = fix * flux3(i,j,k,1) &
!!$      * ( spcl(n)*(0.5d0+signdflx) + spcr(n)*(0.5d0-signdflx) )
!!$    end do
    do n = 1, spn
     spcflx(n,i,j,k,3) = flux3(i,j,k,1) * spc(n,i,j,ks)
    end do
   end if
  end do
 end do

return
end subroutine fluxboundary3i

end module fluxbound_mod
