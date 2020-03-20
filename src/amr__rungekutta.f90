module amr_rungekutta_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_RUNGEKUTTA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate numflux using TVD Runge-Kutta method in AMR mode.

subroutine amr_rungekutta(lf)

 use settings,only:rktype
 use grid,only:is,js,ks,rungen, ufn
 use amr_templates
 use amr_module,only:ib,jb,kb,dt_local,sync,bkn
 use pressure_mod

 implicit none

 integer i,j,k,l,m
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

! Third order TVD Runge-Kutta ################################################
 if(rktype==3)then
! for tflx numbering
  if(.not.sync)then
   l = 1 ; m = 2
  else
   l = 3 ; m = 4
  end if

  if(rungen==1)then

   lf%uorg = lf%u
   if(.not.sync)then
    lf%tflx1 = 0.d0 ; lf%tflx2 = 0.d0 ; lf%tflx3 = 0.d0
   end if

   do ufn = 1, 9
    
    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = lf%uorg(i,j,k,ufn) &
                       + dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

   end do

! Keep temporary surface flux for later restriction
   do k = ks, kb
    do j = js, jb
     lf%tflx1(l,j,k,1:9) = lf%flux1(is-1,j,k,1:9) / 6.d0
     lf%tflx1(m,j,k,1:9) = lf%flux1(ib  ,j,k,1:9) / 6.d0
    end do
   end do
   do k = ks, kb
    do i = is, ib
     lf%tflx2(i,l,k,1:9) = lf%flux2(i,js-1,k,1:9) / 6.d0
     lf%tflx2(i,m,k,1:9) = lf%flux2(i,jb  ,k,1:9) / 6.d0
    end do
   end do
   do j = js, jb
    do i = is, ib
     lf%tflx3(i,j,l,1:9) = lf%flux3(i,j,ks-1,1:9) / 6.d0
     lf%tflx3(i,j,m,1:9) = lf%flux3(i,j,kb  ,1:9) / 6.d0
    end do
   end do

  elseif(rungen==2)then
   do ufn = 1, 9
    
    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = 7.5d-1 * lf%uorg(i,j,k,ufn) + 2.5d-1 * lf%u(i,j,k,ufn)&
                       + 2.5d-1 * dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

   end do

! Keep temporary surface flux for later restriction
   do k = ks, kb
    do j = js, jb
     lf%tflx1(l,j,k,1:9) = lf%tflx1(l,j,k,1:9) + lf%flux1(is-1,j,k,1:9)/6.d0
     lf%tflx1(m,j,k,1:9) = lf%tflx1(m,j,k,1:9) + lf%flux1(ib  ,j,k,1:9)/6.d0
    end do
   end do
   do k = ks, kb
    do i = is, ib
     lf%tflx2(i,l,k,1:9) = lf%tflx2(i,l,k,1:9) + lf%flux2(i,js-1,k,1:9)/6.d0
     lf%tflx2(i,m,k,1:9) = lf%tflx2(i,m,k,1:9) + lf%flux2(i,jb  ,k,1:9)/6.d0
    end do
   end do
   do j = js, jb
    do i = is, ib
     lf%tflx3(i,j,l,1:9) = lf%tflx3(i,j,l,1:9) + lf%flux3(i,j,ks-1,1:9)/6.d0
     lf%tflx3(i,j,m,1:9) = lf%tflx3(i,j,m,1:9) + lf%flux3(i,j,kb  ,1:9)/6.d0
    end do
   end do


  else

! To conserve intermediate values for boundary conditions
   lf%midu1(-2:0,:,:,1:9) = lf%u(is-1:is+1,js-1:jb+1,ks-1:kb+1,1:9)
   lf%midu1( 1:3,:,:,1:9) = lf%u(ib-1:ib+1,js-1:jb+1,ks-1:kb+1,1:9)
   lf%midu2(:,-2:0,:,1:9) = lf%u(is-1:ib+1,js-1:js+1,ks-1:kb+1,1:9)
   lf%midu2(:, 1:3,:,1:9) = lf%u(is-1:ib+1,jb-1:jb+1,ks-1:kb+1,1:9)
   lf%midu3(:,:,-2:0,1:9) = lf%u(is-1:ib+1,js-1:jb+1,ks-1:ks+1,1:9)
   lf%midu3(:,:, 1:3,1:9) = lf%u(is-1:ib+1,js-1:jb+1,kb-1:kb+1,1:9)
    
   do ufn = 1, 9

    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = 1.d0/3.d0 * lf%uorg(i,j,k,ufn) &
                       + 2.d0/3.d0 * lf%u(i,j,k,ufn) + 2.d0/3.d0 * dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

   end do

! Keep temporary surface flux for later restriction
   do k = ks, kb
    do j = js, jb
     lf%tflx1(l,j,k,1:9) = lf%tflx1(l,j,k,1:9) + lf%flux1(is-1,j,k,1:9)*2/3.d0
     lf%tflx1(m,j,k,1:9) = lf%tflx1(m,j,k,1:9) + lf%flux1(ib  ,j,k,1:9)*2/3.d0
    end do
   end do
   do k = ks, kb
    do i = is, ib
     lf%tflx2(i,l,k,1:9) = lf%tflx2(i,l,k,1:9) + lf%flux2(i,js-1,k,1:9)*2/3.d0
     lf%tflx2(i,m,k,1:9) = lf%tflx2(i,m,k,1:9) + lf%flux2(i,jb  ,k,1:9)*2/3.d0
    end do
   end do
   do j = js, jb
    do i = is, ib
     lf%tflx3(i,j,l,1:9) = lf%tflx3(i,j,l,1:9) + lf%flux3(i,j,ks-1,1:9)*2/3.d0
     lf%tflx3(i,j,m,1:9) = lf%tflx3(i,j,m,1:9) + lf%flux3(i,j,kb  ,1:9)*2/3.d0
    end do
   end do

  end if

! Second order TVD Runge-Kutta ###############################################
 elseif(rktype==2)then

  if(rungen==1)then

   lf%uorg = lf%u
   if(.not.sync)then
    lf%tflx1 = 0.d0 ; lf%tflx2 = 0.d0 ; lf%tflx3 = 0.d0
   end if

   do ufn = 1, 9
    
    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = lf%uorg(i,j,k,ufn) &
                       + dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

! Keep temporary surface flux for later restriction
    do k = ks, kb
     do j = js, jb
      lf%tflx1(1,j,k,ufn) = lf%flux1(is-1,j,k,ufn) * 5.d-1
      lf%tflx1(2,j,k,ufn) = lf%flux1(ib  ,j,k,ufn) * 5.d-1
      lf%midu1(-2:0,j,k,ufn) = 5.d-1*&
                          (lf%u(is-1:is+1,j,k,ufn)+lf%uorg(is-1:is+1,j,k,ufn))
      lf%midu1( 1:3,j,k,ufn) = 5.d-1*&
                          (lf%u(ib-1:ib+1,j,k,ufn)+lf%uorg(ib-1:ib+1,j,k,ufn))
     end do
    end do
    do k = ks, kb
     do i = is, ib
      lf%tflx2(i,1,k,ufn) = lf%flux2(i,js-1,k,ufn) * 5.d-1
      lf%tflx2(i,2,k,ufn) = lf%flux2(i,jb  ,k,ufn) * 5.d-1
      lf%midu2(i,-2:0,k,ufn) = 5.d-1*&
                          (lf%u(i,js-1:js+1,k,ufn)+lf%uorg(i,js-1:js+1,k,ufn))
      lf%midu2(i, 1:3,k,ufn) = 5.d-1*&
                          (lf%u(i,jb-1:jb+1,k,ufn)+lf%uorg(i,jb-1:jb+1,k,ufn))
     end do
    end do
    do j = js, jb
     do i = is, ib
      lf%tflx3(i,j,1,ufn) = lf%flux3(i,j,ks-1,ufn) * 5.d-1
      lf%tflx3(i,j,2,ufn) = lf%flux3(i,j,kb  ,ufn) * 5.d-1
      lf%midu3(i,j,-2:0,ufn) = 5.d-1*&
                          (lf%u(i,j,ks-1:ks+1,ufn)+lf%uorg(i,j,ks-1:ks+1,ufn))
      lf%midu3(i,j, 1:3,ufn) = 5.d-1*&
                          (lf%u(i,j,kb-1:kb+1,ufn)+lf%uorg(i,j,kb-1:kb+1,ufn))
     end do
    end do


   end do

  else
   do ufn = 1, 9
    
    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = 5.d-1 * lf%uorg(i,j,k,ufn) + 5.d-1 * lf%u(i,j,k,ufn) &
                       + 5.d-1 * dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

! Keep temporary surface flux for later restriction
    do k = ks, kb
     do j = js, jb
      lf%tflx1(1,j,k,ufn) = lf%tflx1(1,j,k,ufn) + lf%flux1(is-1,j,k,ufn)/2.d0
      lf%tflx1(2,j,k,ufn) = lf%tflx1(2,j,k,ufn) + lf%flux1(ib  ,j,k,ufn)/2.d0
     end do
    end do
    do k = ks, kb
     do i = is, ib
      lf%tflx2(i,1,k,ufn) = lf%tflx2(i,1,k,ufn) + lf%flux2(i,js-1,k,ufn)/2.d0
      lf%tflx2(i,2,k,ufn) = lf%tflx2(i,2,k,ufn) + lf%flux2(i,jb  ,k,ufn)/2.d0
     end do
    end do
    do j = js, jb
     do i = is, ib
      lf%tflx3(i,j,1,ufn) = lf%tflx3(i,j,1,ufn) + lf%flux3(i,j,ks-1,ufn)/2.d0
      lf%tflx3(i,j,2,ufn) = lf%tflx3(i,j,2,ufn) + lf%flux3(i,j,kb  ,ufn)/2.d0
     end do
    end do

   end do
  end if

! Euler method ##############################################################
 elseif(rktype==1)then

   lf%uorg = lf%u
   do ufn = 1, 9
    
    do k = ks, kb
     do j = js, jb
      do i = is, ib
       lf%u(i,j,k,ufn) = lf%uorg(i,j,k,ufn) &
                       + dt_local * &
                       ( lf%idetg1(i) * & 
                         ( lf%detg1(i-1)*lf%flux1(i-1,j,k,ufn)   &
                         - lf%detg1(i  )*lf%flux1(i  ,j,k,ufn) ) &

                       + lf%idetg2(i,j) * &
                         ( lf%detg2(i,j-1)*lf%flux2(i,j-1,k,ufn)   &
                         - lf%detg2(i,j  )*lf%flux2(i,j  ,k,ufn) ) &

                       + lf%idetg3(i,j,k) * &
                         ( lf%flux3(i,j,k-1,ufn)-lf%flux3(i,j,k,ufn) ) &

                       + lf%src(i,j,k,ufn) )
      end do
     end do
    end do

! Keep temporary surface flux for later restriction
    do k = ks, kb
     do j = js, jb
      lf%tflx1(1,j,k,ufn) = lf%flux1(is-1,j,k,ufn)
      lf%tflx1(2,j,k,ufn) = lf%flux1(ib  ,j,k,ufn)
     end do
    end do
    do k = ks, kb
     do i = is, ib
      lf%tflx2(i,1,k,ufn) = lf%flux2(i,js-1,k,ufn)
      lf%tflx2(i,2,k,ufn) = lf%flux2(i,jb  ,k,ufn)
     end do
    end do
    do j = js, jb
     do i = is, ib
      lf%tflx3(i,j,1,ufn) = lf%flux3(i,j,ks-1,ufn)
      lf%tflx3(i,j,2,ufn) = lf%flux3(i,j,kb  ,ufn)
     end do
    end do

   end do

 else
  print *,"Error from rktype in AMR mode",rktype
  stop
 end if

 call amr_primitive(lf)

return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE AMR_PRIMITIVE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert conserved values to primitive values in AMR mode

subroutine amr_primitive(lf)

 use grid,only:is,js,ks
 use amr_templates
 use amr_module,only:ib,jb,kb

 implicit none

 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

 do k = ks, kb
  do j = js, jb
   do i = is, ib
    lf% d(i,j,k) = lf%u(i,j,k,1)
    lf%v1(i,j,k) = lf%u(i,j,k,2) / lf%d(i,j,k)
    lf%v2(i,j,k) = lf%u(i,j,k,3) / lf%d(i,j,k)
    lf%v3(i,j,k) = lf%u(i,j,k,4) / lf%d(i,j,k)
    lf%b1(i,j,k) = lf%u(i,j,k,5)
    lf%b2(i,j,k) = lf%u(i,j,k,6)
    lf%b3(i,j,k) = lf%u(i,j,k,7)
    lf% e(i,j,k) = lf%u(i,j,k,8)
    ! for 9 wave method
    lf%phi(i,j,k) = lf%u(i,j,k,9)
   end do
  end do
 end do

!!$ do k = ks, kb
!!$  do j = js, jb
!!$   do i = is, ib
!!$    if(lf%d(i,j,k)<=0.d0)lf%d(i,j,k) = minval(lf%d(is-1:ib+1,js-1:jb+1,ks-1:kb+1),MASK=lf%d>0.d0)
!!$    if(lf%e(i,j,k)<=0.d0)lf%e(i,j,k) = minval(lf%e(is-1:ib+1,js-1:jb+1,ks-1:kb+1),MASK=lf%e>0.d0)
!!$   end do
!!$  end do
!!$ end do

 call pressure_block(lf)

return
end subroutine amr_primitive

end subroutine amr_rungekutta

end module amr_rungekutta_mod
