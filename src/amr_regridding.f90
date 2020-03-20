module amr_regridding_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_REGRIDDING
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set subgrids for Adaptive Mesh Refinement.

subroutine amr_regridding

 use amr_module
 use amr_interpolation_mod
 use amr_allocation
 use amr_boundary_mod
use grid,only:rungen
 implicit none

 integer old_totbloks, lastb, level

!-----------------------------------------------------------------------------

 old_totbloks = totbloks

 call amr_allocateoldblock
 call amr_allocateoldleafs

! Flagging refinement blocks
 do level = 0, curlvl
  do bkn = 1, old_totbloks
   if(bk(bkn)%level==level)then
    call amr_boundary
    call amr_criterion(lf(bkn),bk(bkn))
   end if
  end do
 end do
!print *,bk%ref_deref;stop
 obk = bk

! do not derefine unless all blocks in the parent block is flagged -1
! parents of refining blocks should be flagged 1
 do bkn = 1, old_totbloks
  if(bk(bkn)%ref_deref==-1)then
   do lid = 1, cib
    if(bk(bk(bk(bkn)%parnt)%child(lid))%ref_deref/=-1)then
     bk(bkn)%ref_deref=0 ; exit
    end if
   end do
  elseif(bk(bkn)%ref_deref==1.and.bk(bkn)%level>0)then
   bk(bk(bkn)%parnt)%ref_deref=1
  end if
 end do

! flag grandparents of leaf blocks
 do bkn = 1, old_totbloks
  if(bk(bkn)%level>=2)then
   if(bk(bkn)%child(1)==0)then
    bk(bk(bk(bkn)%parnt)%parnt)%ref_deref=1
   end if
  end if
 end do

! refining all neighbouring cells to refining blocks
 do bkn = 1, old_totbloks
  if(obk(bkn)%ref_deref==1)then
   do lid = 1, face
    if(bk(bkn)%neigh(lid)>0)then
     bk(bk(bkn)%neigh(lid))%ref_deref = 1
    end if
   end do
  end if
 end do

 obk = bk

!!$! reflagging so that level differences don't exceed unity
!!$ do bkn = 1, old_totbloks
!!$  do lid = 1, face
!!$   if(bk(bkn)%neigh(lid)>0)then
!!$    if(bk(bkn)%level+bk(bkn)%ref_deref-bk(bk(bkn)%neigh(lid))%level-bk(bk(bkn)%neigh(lid))%ref_deref>1)then
!!$     bk(bk(bkn)%neigh(lid))%ref_deref = bk(bk(bkn)%neigh(lid))%ref_deref + 1
!!$    end if
!!$   end if
!!$  end do
!!$ end do

! reflagging so that all levels at least have 2 blocks as margin
!!$ do bkn = 1, old_totbloks
!!$  do lid = 1, face
!!$   if(obk(bkn)%neigh(lid)>0)then
!!$    if(obk(obk(bkn)%neigh(lid))%ref_deref==1.or.obk(bkn)%ldif(lid)==1)then
!!$     do lfn = 1, face
!!$      if(obk(bkn)%ldif(lfn)==-1)bk(bk(bkn)%neigh(lfn))%ref_deref = 1     
!!$     end do
!!$     exit
!!$    end if
!!$   end if
!!$  end do
!!$ end do

 obk = bk

! derefine only if derefinement was flagged for all blocks in a parent block
 do bkn = 1, totbloks
  if(obk(bkn)%ref_deref==-1)then
! do not refine if parent or its neighbour is flagged to refine
   if(obk(obk(bkn)%parnt)%ref_deref==1)then
    bk(bkn)%ref_deref = 0 ; cycle
!!$   elseif(maxval(obk(obk(bkn)%parnt)%ldif)==1)then
!!$    bk(bkn)%ref_deref = 0 ; cycle
   end if
! do not derefine if all neighbours are higher levels
   do lfn = 1, face
    if(obk(obk(bkn)%parnt)%neigh(lfn)>=0)then
     if(obk(obk(bkn)%parnt)%ldif(lfn) + &
        obk(obk(obk(bkn)%parnt)%neigh(lfn))%ref_deref==2)then
      bk(bkn)%ref_deref = 0
      exit
     end if
    end if
   end do
!!$   do lid = 1, cib
!!$    if(obk(obk(obk(bkn)%parnt)%child(lid))%ref_deref/=-1)then
!!$     bk(bkn)%ref_deref = 0 ; exit
!!$    end if
!!$   end do
  end if
  if(bk(bkn)%level==maxamr.and.bk(bkn)%ref_deref==1)bk(bkn)%ref_deref = 0
  if(bk(bkn)%child(1)/=0)bk(bkn)%ref_deref = 0
 end do


! counting number of new total blocks
 do bkn = 1, old_totbloks
  if(bk(bkn)%ref_deref==1)then
   totbloks = totbloks + cib
  elseif(bk(bkn)%ref_deref==-1)then
   totbloks = totbloks - 1
  end if
 end do

 obk = bk

 deallocate(bk)
 call amr_allocateblock

 allocate(regrid(-face:old_totbloks)); regrid = 0
 do lid = -face, -1
  regrid(lid) = lid
 end do

! Start renumbering blocks
!  regrid : converts old ID to new ID
 bkn = 1 ; lastb = 1
 do
  if(obk(lastb)%ref_deref==0)then
   bk(bkn)%level = obk(lastb)%level
   regrid(lastb) = bkn
   bkn = bkn + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==1)then
   bk(bkn)%level = obk(lastb)%level
   regrid(lastb) = bkn
   do lid = 1, cib
    bk(bkn+lid)%parnt = bkn
    bk(bkn+lid)%level = obk(lastb)%level + 1
   end do
   bkn = bkn + cib + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==-1)then
   regrid(lastb) = 0
   lastb = lastb + cib
  else
   print *,"Error from amr renumbering",obk%ref_deref
   stop
  end if
  if(lastb>old_totbloks)exit
 end do
! End of renumbering

! Renumbering amr header values for each block
 bkn = 1 ; lastb = 1
 do
  if(obk(lastb)%ref_deref==0)then ! for unchanged blocks-------------------
   bk(bkn)%parnt = regrid(obk(lastb)%parnt)
   do lid = 1, cib
    bk(bkn)%child(lid) = regrid(obk(lastb)%child(lid))
   end do
   do lid = 1, face
    bk(bkn)%neigh(lid) = regrid(obk(lastb)%neigh(lid))
    if(regrid(obk(lastb)%neigh(lid))==0)then ! if neighbour is derefined
     bk(bkn)%neigh(lid) = regrid(obk(obk(lastb)%neigh(lid))%parnt)
    end if
    if(obk(lastb)%neigh(lid)<0)then ! if neighbour is boundary
     bk(bkn)%neigh(lid) = obk(lastb)%neigh(lid)
     bk(bkn)%ldif(lid)  = -2
    else
     bk(bkn)%ldif(lid)  = obk(lastb)%ldif(lid) &
                          + obk(obk(lastb)%neigh(lid))%ref_deref
    end if
   end do
   bkn   = bkn  + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==1)then ! for refined blocks-----------------
!  Refined block becomes a parent
   bk(bkn)%parnt = regrid(obk(lastb)%parnt)
   do lid = 1, cib
    bk(bkn)%child(lid) = bkn + lid!regrid(obk(lastb)%child(lfn))
   end do
   do lid = 1, face
    bk(bkn)%neigh(lid) = regrid(obk(lastb)%neigh(lid))
    if(bk(bkn)%neigh(lid)==0)then ! if neighbour is derefined
     bk(bkn)%neigh(lid) = regrid(obk(obk(lastb)%neigh(lid))%parnt)
     bk(bkn)%ldif(lid)  = obk(lastb)%ldif(lid) - 1!&
!                          + obk(obk(lastb)%neigh(lid))%ref_deref
    end if
    if(obk(lastb)%neigh(lid)<0)then ! if neighbour is boundary
     bk(bkn)%neigh(lid) = obk(lastb)%neigh(lid)
     bk(bkn)%ldif(lid)  = -2
    else
     bk(bkn)%ldif(lid)  = obk(lastb)%ldif(lid) &
                          + obk(obk(lastb)%neigh(lid))%ref_deref
    end if
   end do
   bkn = bkn + cib + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==-1)then ! for derefined blocks--------------
   bk(bkn-1)%child = 0
   lastb = lastb + cib
  else
   print *,"Error from amr renumbering 2"
   stop
  end if
  if(lastb>old_totbloks)exit
 end do

! Set neighbours for new blocks+++++++++++++++++++++++++++++++++++++++++++++
 bkn = 1 ; lastb = 1
 do
  if(obk(lastb)%ref_deref==0)then
   bkn  = bkn  + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==1)then
   do lid = 1, cib
    bk(bkn+lid)%parnt = bkn
    do lfn = 1, face
     if(bk(bkn)%neigh(lfn)<0)then
      bk(bkn+lid)%neigh(lfn) = bk(bkn)%neigh(lfn)
      bk(bkn+lid)%ldif(lfn)  = -2
     elseif(bk(bkn)%neigh(lfn)>0)then
      bk(bkn+lid)%neigh(lfn) = bk(bk(bkn)%neigh(lfn))%&
                       child(lid+(-1)**((lid-1)/2**((lfn-1)/2))*2**((lfn-1)/2))
      bk(bkn+lid)% ldif(lfn) = 0
     end if
     if(bk(bkn+lid)%neigh(lfn)==0)then
      bk(bkn+lid)%neigh(lfn) = bk(bkn)%neigh(lfn)
      bk(bkn+lid)% ldif(lfn) = -1
     end if
     if(localneigh(lid,lfn)/=0)then
      bk(bkn+lid)%neigh(lfn) = bkn + localneigh(lid,lfn)
      bk(bkn+lid)% ldif(lfn) = 0
     end if
     if(bk(bkn+lid)%ldif(lfn)==0.or.bk(bkn+lid)%ldif(lfn)==1)then
      bk(bk(bkn+lid)%neigh(lfn))%neigh(lfn-(-1)**lfn) = bkn + lid
      bk(bk(bkn+lid)%neigh(lfn))%ldif (lfn-(-1)**lfn) = 0
     end if
    end do
   end do
   bkn = bkn + cib + 1
   lastb = lastb + 1
  elseif(obk(lastb)%ref_deref==-1)then
   lastb = lastb + cib
  end if
  if(lastb>old_totbloks)exit
 end do

! Calculating ldif for each block++++++++++++++++++++++++++++++++++++++++++++

 do bkn = 1, totbloks
  do lfn = 1, face
   if(bk(bkn)%neigh(lfn)<0)then
    bk(bkn)%ldif(lfn) = -2
   elseif(bk(bk(bkn)%neigh(lfn))%child(1)/=0)then
    bk(bkn)%ldif(lfn) = 1
   elseif(bk(bk(bkn)%neigh(lfn))%level-bk(bkn)%level==-1)then
    bk(bkn)%ldif(lfn) = - 1
   else
    bk(bkn)%ldif(lfn) = 0
   end if
  end do
 end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$do bkn=1,totbloks
!!$print *,bkn,bk(bkn)%ldif(1:2),bk(bkn)%level;enddo
! Reordering leaf data

 olf = lf

 deallocate(lf)
 call amr_allocateleafs

do bkn = 1, old_totbloks
! Just inherit old leaf for unrefined blocks
 if(obk(bkn)%ref_deref==0)then
  lf(regrid(bkn)) = olf(bkn)

! Refine blocks into smaller leaves
 elseif(obk(bkn)%ref_deref==1)then
  lfn = regrid(bkn)
  lf(lfn) = olf(bkn)
  call amr_interpolation(lf(lfn+1:lfn+cib),olf(bkn))

! Deallocate derefined blocks
 elseif(obk(bkn)%ref_deref==-1)then
  cycle

 else
  print *,"Error from ref_deref of block number = ",bkn
  stop
 end if
end do

deallocate(obk,olf,regrid)

curlvl=maxval(bk%level)

bk%ref_deref = 0

return
end subroutine amr_regridding

end module amr_regridding_mod
