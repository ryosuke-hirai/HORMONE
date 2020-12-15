!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE INITIALCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine initialcondition

 use settings,only:start,eostype
 use grid
 use physval
 use pressure_mod

 implicit none

!----------------------------------------------------------------------------

 if(start==0)then

  call agndisk

! to set internal energy consistently ------------------------------!
  call meanmolweight                                                !
  select case (eostype)                                             !
  case(0:1) ! without recombination                                 !
   do k=ks,ke ; do j=js,je ; do i=is,ie                             !
    e(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k)) &       !
  + d(i,j,k)*0.5d0*( v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2 ) & !
  +          0.5d0*( b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2 )   !
   end do; end do ; end do                                          !
  case(2) ! with recombination                                      !
   do k=ks,ke ; do j=js,je ; do i=is,ie                             !
    e(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k), &       !
                     spc(1,i,j,k),spc(2,i,j,k)) &                   !
  + d(i,j,k)*0.5d0*( v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2 ) & !
  +          0.5d0*( b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2 )   !
   end do; end do ; end do                                          !
  end select                                                        !
! ----------------------------------------------------------------- !

 else

  call restart
  call meanmolweight

 end if

return
end subroutine initialcondition
