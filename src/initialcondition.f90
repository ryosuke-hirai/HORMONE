module initialcondition_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE INITIALCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine initialcondition

 use settings,only:start,eostype,simtype
 use grid
 use physval
 use pressure_mod
 use composition_mod
 use restart_mod

 use eostest_mod
 use shocktube_mod
 use sedov_mod
 use orszagtang_mod
 use KHinstability_mod
 use rad_box_mod
 use star_init_mod
 use redsupergiant_mod
 use agndisk_mod
 use modify_mod
 use windtunnel_mod
 use stellarcollision_mod
 use iotest_mod

 integer:: i,j,k

!----------------------------------------------------------------------------

 if(start==0)then

  select case(simtype)
  case('eostest')
   call eostest

  case('sodshock_x','sodshock_y','sodshock_z',&
       'briowushock_x','briowushock_y','briowushock_z',&
       'other_shocktube_x','other_shocktube_y','other_shocktube_z')
   call shocktube

  case('sedov_default','sedov_other')
   call sedov

  case('orszagtang_xy','orszagtang_yz','orszagtang_xz')
   call orszagtang

  case('KHinstability')
   call KHinstability

  case('rad_box')
   call rad_box

  case('star_sph')
   call star_init

  case('rsg_sph')
   call redsupergiant

  case('agndisk')
   call agndisk

  case('windtunnel')
   call windtunnel

  case('stellarcollision')
   call stellarcollision

  case('modify')
   call modify

  case('iotest')
    call iotest

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

end module initialcondition_mod
