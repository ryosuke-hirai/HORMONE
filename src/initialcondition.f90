!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE INITIALCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine initialcondition

 use settings,only:start,eostype,simtype
 use grid
 use physval
 use pressure_mod
 use composition_mod
 
 implicit none

!----------------------------------------------------------------------------

 if(start==0)then

  select case(simtype)
  case('eostest')
   call eostest
  case('orszagtang_xy','orszagtang_yz','orszagtang_xz')
   call orszagtang
  case('sodshock_x')
   gamma = 1.4d0
   call shocktube(1,    1d0,  1d0,0d0,0d0,0d0,0d0,0d0,0d0,&
                    0.125d0,0.1d0,0d0,0d0,0d0,0d0,0d0,0d0)
  case('sodshock_y')
   gamma = 1.4d0
   call shocktube(2,    1d0,  1d0,0d0,0d0,0d0,0d0,0d0,0d0,&
                    0.125d0,0.1d0,0d0,0d0,0d0,0d0,0d0,0d0)
  case('sodshock_z')
   gamma = 1.4d0
   call shocktube(3,    1d0,  1d0,0d0,0d0,0d0,0d0,0d0,0d0,&
                    0.125d0,0.1d0,0d0,0d0,0d0,0d0,0d0,0d0)
  case('briowushock_x')
   gamma = 2d0
   call shocktube(1,    1d0,  1d0,0d0,0d0,0d0,0.75d0, 1d0,0d0,&
    & 0.125d0,0.1d0,0d0,0d0,0d0,0.75d0,-1d0,0d0)
  case('briowushock_y')
   gamma = 2d0
   call shocktube(2,    1d0,  1d0,0d0,0d0,0d0,0d0,0.75d0, 1d0,&
                    0.125d0,0.1d0,0d0,0d0,0d0,0d0,0.75d0,-1d0)
  case('briowushock_z')
   gamma = 2d0
   call shocktube(3,    1d0,  1d0,0d0,0d0,0d0, 1d0,0d0,0.75d0,&
                    0.125d0,0.1d0,0d0,0d0,0d0,-1d0,0d0,0.75d0)
  case('rsg_sph')
   call redsupergiant
  case default
   print*,'Chosen simtype not available, simtype = ',simtype
   stop
  end select

! to set internal energy consistently ------------------------------!
  call meanmolweight                                                !
  select case (eostype)                                             !
  case(0:1) ! without recombination                                 !
   do k=ks,ke ; do j=js,je ; do i=is,ie                             !
    eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))      !
    e(i,j,k) = eint(i,j,k) &                                        !
  + d(i,j,k)*0.5d0*( v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2 ) & !
  +          0.5d0*( b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2 )   !
   end do; end do ; end do                                          !
  case(2) ! with recombination                                      !
   do k=ks,ke ; do j=js,je ; do i=is,ie                             !
    eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k), &    !
                        spc(1,i,j,k),spc(2,i,j,k))                  !
    e(i,j,k) = eint(i,j,k) &                                        !
  + d(i,j,k)*0.5d0*( v1(i,j,k)**2 + v2(i,j,k)**2 + v3(i,j,k)**2 ) & !
  +          0.5d0*( b1(i,j,k)**2 + b2(i,j,k)**2 + b3(i,j,k)**2 )   !
   end do; end do ; end do                                          !
  end select                                                        !
! ----------------------------------------------------------------- !

 else

  call restart

 end if

return
end subroutine initialcondition
