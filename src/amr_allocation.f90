module amr_allocation

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_ALLOCATEBLOCK
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate memory for refined blocks.

! CAUTION: Always deallocate before usage.

subroutine amr_allocateblock

 use amr_module

 implicit none

!-----------------------------------------------------------------------------

 allocate  ( bk(totbloks) )

 do bkn = 1, totbloks
  allocate( bk(bkn)%child(cib), bk(bkn)%neigh(face), bk(bkn)%ldif(face) )
  bk(bkn)%child = 0
  bk(bkn)%neigh = 0
  bk(bkn)%ldif  = 0
  bk(bkn)%parnt = 0
  bk(bkn)%level = 0
  bk(bkn)%ref_deref = 0
 end do

return
end subroutine amr_allocateblock


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE AMR_ALLOCATEOLDBLOCK
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate memory for refined blocks.

! CAUTION: Always deallocate before usage.

subroutine amr_allocateoldblock

 use amr_module

 implicit none

!-----------------------------------------------------------------------------

 allocate  ( obk(totbloks) )

 do bkn = 1, totbloks
  allocate( obk(bkn)%child(cib), obk(bkn)%neigh(face), obk(bkn)%ldif(face) )
  obk(bkn)%child = 0
  obk(bkn)%neigh = 0
  obk(bkn)%ldif  = 0
  obk(bkn)%parnt = 0
  obk(bkn)%level = 0
  obk(bkn)%ref_deref = 0
 end do

return
end subroutine amr_allocateoldblock


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_ALLOCATELEAFS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate memory for leafs.

! CAUTION: Always deallocate before usage.

subroutine amr_allocateleafs

 use settings,only:crdnt
 use amr_module

 implicit none

!-----------------------------------------------------------------------------

 allocate  ( lf(totbloks)  )

 do bkn = 1, totbloks
  allocate( &
   lf(bkn)%x1(-1:ib+2)  , lf(bkn)%xi1(-1:ib+2) , lf(bkn)%dx1(-1:ib+2),   &
!   lf(bkn)%dxi1(-1:ib+2),
   lf(bkn)%idx1(-1:ib+2),&! lf(bkn)%idxi1(-1:ib+2), &
   lf(bkn)%x2(-1:jb+2)  , lf(bkn)%xi2(-1:jb+2) , lf(bkn)%dx2(-1:jb+2),   &
!   lf(bkn)%dxi2(-1:jb+2),
   lf(bkn)%idx2(-1:jb+2),&! lf(bkn)%idxi2(-1:jb+2), &
   lf(bkn)%x3(-1:kb+2)  , lf(bkn)%xi3(-1:kb+2) , lf(bkn)%dx3(-1:kb+2),   &
!   lf(bkn)%dxi3(-1:kb+2),
   lf(bkn)%idx3(-1:kb+2),&! lf(bkn)%idxi3(-1:kb+2), &
   
   lf(bkn)%d(-1:ib+2,-1:jb+2,-1:kb+2), lf(bkn)%p(-1:ib+2,-1:jb+2,-1:kb+2), &
   lf(bkn)%e(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   lf(bkn)%v1(-1:ib+2,-1:jb+2,-1:kb+2),lf(bkn)%v2(-1:ib+2,-1:jb+2,-1:kb+2),&
   lf(bkn)%v3(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   lf(bkn)%b1(-1:ib+2,-1:jb+2,-1:kb+2),lf(bkn)%b2(-1:ib+2,-1:jb+2,-1:kb+2),&
   lf(bkn)%b3(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   lf(bkn)%ptot(-1:ib+2,-1:jb+2,-1:kb+2),lf(bkn)%cf(-1:ib+2,-1:jb+2,-1:kb+2), &
   lf(bkn)% phi(-1:ib+2,-1:jb+2,-1:kb+2), &
   lf(bkn)%gphi(-1:ib+2,-1:jb+2,-1:kb+2), &
   lf(bkn)%dvol(-1:ib+2,-1:jb+2,-1:kb+2), &

   lf(bkn)%detg1(0:ib+1), lf(bkn)%idetg1(0:ib+1),&
   lf(bkn)%detg2(0:ib+1,0:jb+1), lf(bkn)%idetg2(0:ib+1,0:jb+1),&
   lf(bkn)%idetg3(0:ib+1,0:jb+1,0:kb+1),&

   lf(bkn)%sx1(0:ib+1),&

   lf(bkn)%g22(0:ib+1), lf(bkn)%g33(0:ib+1,0:jb+1),&
   
   lf(bkn)%u(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   lf(bkn)%uorg(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   lf(bkn)%flux1(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   lf(bkn)%flux2(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   lf(bkn)%flux3(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   lf(bkn)%src(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   
! tflx is for flux restriction at coarser-cell boundaries
! 1,3 => for coarse inner face flux ; 2,4 => for coarse outer face flux
! 1,2 => for first step ; 3,4 => for second step
   lf(bkn)%tflx1(1:4,1:jb,1:kb,1:9), &
   lf(bkn)%tflx2(1:ib,1:4,1:kb,1:9), &
   lf(bkn)%tflx3(1:ib,1:jb,1:4,1:9), &

! midu is for coarse boundary condition time interpolation
! -2:0 => for inner boundaries
!  1:3 => for outer boundaries
   lf(bkn)%midu1(-2:3,0:jb+1,0:kb+1,1:9),&
   lf(bkn)%midu2(0:ib+1,-2:3,0:kb+1,1:9),&
   lf(bkn)%midu3(0:ib+1,0:jb+1,-2:3,1:9) )

  if(crdnt==2) allocate( lf(bkn)%scot(0:jb+1), lf(bkn)%sisin(0:jb+1) )
  if(crdnt/=2) allocate( lf(bkn)%scot(1:1)   , lf(bkn)%sisin(1:1)    )

  lf(bkn)%d = 2.d0
  lf(bkn)%u(:,:,:,1) = 1.d5 ; lf(bkn)%u(:,:,:,2:9) = 0.d0
  lf(bkn)%v1 = 0.d0 ; lf(bkn)%v2 = 0.d0 ; lf(bkn)%v3 = 0.d0
  lf(bkn)%b1 = 0.d0 ; lf(bkn)%b2 = 0.d0 ; lf(bkn)%b3 = 0.d0
  lf(bkn)%phi= 0.d0 ; lf(bkn)%gphi = 0.d0
  lf(bkn)%uorg = 0.d0
  lf(bkn)%midu1 = 0.d0 ; lf(bkn)%midu2 = 0.d0 ; lf(bkn)%midu3 = 0.d0

 end do

return
end subroutine amr_allocateleafs



!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE AMR_ALLOCATEOLDLEAFS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To allocate memory for leafs.

! CAUTION: Always deallocate before usage.

subroutine amr_allocateoldleafs

 use settings,only:crdnt
 use amr_module

 implicit none

!-----------------------------------------------------------------------------

 allocate  ( olf(totbloks)  )

 do bkn = 1, totbloks
  allocate( &
   olf(bkn)%x1(-1:ib+2)  , olf(bkn)%xi1(-1:ib+2) , olf(bkn)%dx1(-1:ib+2),   &
!   olf(bkn)%dxi1(-1:ib+2),
   olf(bkn)%idx1(-1:ib+2),&! olf(bkn)%idxi1(-1:ib+2), &
   olf(bkn)%x2(-1:jb+2)  , olf(bkn)%xi2(-1:jb+2) , olf(bkn)%dx2(-1:jb+2),   &
!   olf(bkn)%dxi2(-1:jb+2),
   olf(bkn)%idx2(-1:jb+2),&! olf(bkn)%idxi2(-1:jb+2), &
   olf(bkn)%x3(-1:kb+2)  , olf(bkn)%xi3(-1:kb+2) , olf(bkn)%dx3(-1:kb+2),   &
!   olf(bkn)%dxi3(-1:kb+2),
   olf(bkn)%idx3(-1:kb+2),&! olf(bkn)%idxi3(-1:kb+2), &
   
   olf(bkn)%d(-1:ib+2,-1:jb+2,-1:kb+2), olf(bkn)%p(-1:ib+2,-1:jb+2,-1:kb+2), &
   olf(bkn)%e(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   olf(bkn)%v1(-1:ib+2,-1:jb+2,-1:kb+2),olf(bkn)%v2(-1:ib+2,-1:jb+2,-1:kb+2),&
   olf(bkn)%v3(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   olf(bkn)%b1(-1:ib+2,-1:jb+2,-1:kb+2),olf(bkn)%b2(-1:ib+2,-1:jb+2,-1:kb+2),&
   olf(bkn)%b3(-1:ib+2,-1:jb+2,-1:kb+2), &
   
   olf(bkn)%ptot(-1:ib+2,-1:jb+2,-1:kb+2),olf(bkn)%cf(-1:ib+2,-1:jb+2,-1:kb+2), &
   olf(bkn)%phi(-1:ib+2,-1:jb+2,-1:kb+2), &
   olf(bkn)%gphi(-1:ib+2,-1:jb+2,-1:kb+2), &
   olf(bkn)%dvol(-1:ib+2,-1:jb+2,-1:kb+2),&

   olf(bkn)%detg1(0:ib+1), olf(bkn)%idetg1(0:ib+1),&
   olf(bkn)%detg2(0:ib+1,0:jb+1), olf(bkn)%idetg2(0:ib+1,0:jb+1),&
   olf(bkn)%idetg3(0:ib+1,0:jb+1,0:kb+1),&

   olf(bkn)%sx1(0:ib+1),&

   olf(bkn)%g22(0:ib+1), olf(bkn)%g33(0:ib+1,0:jb+1),&
   
   olf(bkn)%u(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   olf(bkn)%uorg(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   olf(bkn)%flux1(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   olf(bkn)%flux2(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   olf(bkn)%flux3(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   olf(bkn)%src(-1:ib+2,-1:jb+2,-1:kb+2,1:9), &
   
   olf(bkn)%tflx1(1:4,1:jb,1:kb,1:9),&
   olf(bkn)%tflx2(1:ib,1:4,1:kb,1:9),&
   olf(bkn)%tflx3(1:ib,1:jb,1:4,1:9), &

   olf(bkn)%midu1(-2:3,1:jb,1:kb,1:9),&
   olf(bkn)%midu2(1:ib,-2:3,1:kb,1:9),&
   olf(bkn)%midu3(1:ib,1:jb,-2:3,1:9) )

  if(crdnt==2) allocate( olf(bkn)%scot(0:jb+1), olf(bkn)%sisin(0:jb+1) )

  olf(bkn)%d = 3.d0

 end do

return
end subroutine amr_allocateoldleafs

end module amr_allocation
