module fluxbound_mod
 implicit none

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

 if(bc1is==10)call fluxboundary1i
 if(bc1os==10)call fluxboundary1o
 if(bc2is==10)call fluxboundary2i
 if(bc2os==10)call fluxboundary2o
 if(bc3is==10)call fluxboundary3i
 if(bc3os==10)call fluxboundary3o

return
end subroutine fluxboundary


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE FLUXBOUNDARY1I
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set flux for x1 inner boundary

subroutine fluxboundary1i

 use grid
 use physval,only:flux1

 real(8)::fix,signdflx

!-----------------------------------------------------------------------------

 i = is-1
 do k = ks, ke
  do j = js, je
   flux1(i,j,k,1) = 0d0
   if(compswitch>=2)then
    do n = 1, spn
     spcl(n) = spc(n,i  ,j,k) + dx(1) * dspc(n,i  ,j,k,1)
     spcr(n) = spc(n,i+1,j,k) - dx(2) * dspc(n,i+1,j,k,1)
    end do
    signdflx = sign(0.5d0,flux1(i,j,k,1))
    fix = 1d0/( sum(spcl(1:spn))*(0.5d0+signdflx) + &
     sum(spcr(1:spn))*(0.5d0-signdflx) )
    do n = 1, spn
     spcflx(n,i,j,k,1) = fix * flux1(i,j,k,1) &
      * ( spcl(n)*(0.5d0+signdflx) + spcr(n)*(0.5d0-signdflx) )
    end do
   end if
  end do
 end do

return
end subroutine fluxboundary1i

end module fluxbound_mod
