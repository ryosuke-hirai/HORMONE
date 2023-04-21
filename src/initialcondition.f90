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

 integer:: dir
!----------------------------------------------------------------------------

 if(start==0)then

  select case(simtype)
  case('eostest')
   call eostest
   
  case('orszagtang_xy','orszagtang_yz','orszagtang_xz')
   call orszagtang

  case('sodshock_x','sodshock_y','sodshock_z')
   gamma = 1.4d0
   if(simtype(10:10)=='x')dir=1
   if(simtype(10:10)=='y')dir=2
   if(simtype(10:10)=='z')dir=3
   call shocktube(dir,    1d0,  1d0,0d0,0d0,0d0,0d0,0d0,0d0,&
                      0.125d0,0.1d0,0d0,0d0,0d0,0d0,0d0,0d0)

  case('briowushock_x','briowushock_y','briowushock_z')
   gamma = 2d0
   if(simtype(13:13)=='x')dir=1
   if(simtype(13:13)=='y')dir=2
   if(simtype(13:13)=='z')dir=3
   call shocktube(dir,    1d0,  1d0,0d0,0d0,0d0,0.75d0, 1d0,0d0,&
                      0.125d0,0.1d0,0d0,0d0,0d0,0.75d0,-1d0,0d0)

  case('sedov')
   gamma=1.4d0
   call sedov(1d0,1d0)

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
