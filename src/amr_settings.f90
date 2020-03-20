!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE AMR_SETTINGS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set variables for AMR.(Should be after gridset.f)

subroutine amr_settings

  use settings,only:start
  use grid
  use physval
  use gravmod
  use amr_module
  use amr_allocation
  use amr_interpolation_mod
  use amr_regridding_mod
  use pressure_mod

  implicit none

  integer bni, bnj, bnk, i1,i2,j1,j2,k1,k2
  character*20 startfile

!-----------------------------------------------------------------------------

  cib  = cuts**dim ! Number of cells per subgrid block for AMR mode.
  face = dim*2     ! Number of faces per cell for AMR mode.

! set matrix for local neighbour numbering
 localneigh = 0
 localneigh(1,2) = 2 ;  localneigh(1,4) = 3 ;  localneigh(1,6) = 5
 localneigh(2,1) = 1 ;  localneigh(2,4) = 4 ;  localneigh(2,6) = 6
 localneigh(3,2) = 4 ;  localneigh(3,3) = 1 ;  localneigh(3,6) = 7
 localneigh(4,1) = 3 ;  localneigh(4,3) = 2 ;  localneigh(4,6) = 8
 localneigh(5,2) = 6 ;  localneigh(5,4) = 7 ;  localneigh(5,5) = 1
 localneigh(6,1) = 5 ;  localneigh(6,4) = 8 ;  localneigh(6,5) = 2
 localneigh(7,2) = 8 ;  localneigh(7,3) = 5 ;  localneigh(7,5) = 3
 localneigh(8,1) = 7 ;  localneigh(8,3) = 6 ;  localneigh(8,5) = 4

if(start==0)then
 
 call initialcondition
 call boundarycondition
 call pressure
 call gravsetup

 ! Construct blocks on original grid
 bni = ie/ib ; bnj = je/jb ; bnk = ke/kb
 totbloks = bni * bnj * bnk

 call amr_allocateblock
 call amr_allocateleafs

 do bkn = 1, totbloks
  bk(bkn)%parnt = 0 ;  bk(bkn)%level = 0
  bk(bkn)%child = 0 ;  bk(bkn)%ldif  = 0
  do lfn = 1, face
   bk(bkn)%neigh(lfn) = bkn +(-1)**lfn * bni**((lfn+1)/4) * bnj**(lfn/5)
  end do
  if(mod(bkn,bni)==1.or.bni==1)    bk(bkn)%neigh(1) = -1 ! ii boundary
  if(mod(bkn,bni)==1.or.bni==1)    bk(bkn)% ldif(1) = -2 ! ii boundary
  if(mod(bkn,bni)==0)              bk(bkn)%neigh(2) = -2 ! io boundary
  if(mod(bkn,bni)==0)              bk(bkn)% ldif(2) = -2 ! io boundary
  if(dim>=2)then
   if(mod(bkn,bni*bnj)<=bni)        bk(bkn)%neigh(3) = -3 ! ji boundary
   if(mod(bkn,bni*bnj)<=bni)        bk(bkn)% ldif(3) = -2 ! ji boundary
   if(mod(bkn,bni*bnj)>bni*(bnj-1).or.mod(bkn,bni*bnj)==0)then
                                    bk(bkn)%neigh(4) = -4 ! jo boundary
                                    bk(bkn)% ldif(4) = -2 ! jo boundary
   end if
   if(dim==3)then
    if(bkn<=bni*bnj)                 bk(bkn)%neigh(5) = -5 ! ki boundary
    if(bkn<=bni*bnj)                 bk(bkn)% ldif(5) = -2 ! ki boundary
    if(bkn>totbloks-bni*bnj)         bk(bkn)%neigh(6) = -6 ! ko boundary
    if(bkn>totbloks-bni*bnj)         bk(bkn)% ldif(6) = -2 ! ko boundary
   end if
  end if

  lf(bkn)%xi1(is-1) = xi1s + (xi1e-xi1s)/dble(bni) * dble(mod(bkn,bni)-1)
  lf(bkn)%xi1(ib  ) = xi1s + (xi1e-xi1s)/dble(bni) * dble(mod(bkn,bni))
  if(mod(bkn,bni)==0)then
   lf(bkn)%xi1(is-1) = xi1s + (xi1e-xi1s)/dble(bni) * dble(bni-1)
   lf(bkn)%xi1(ib  ) = xi1s + (xi1e-xi1s)/dble(bni) * dble(bni)
  end if
  lf(bkn)%xi2(js-1) = xi2s+(xi2e-xi2s)/dble(bnj)*dble((mod(bkn,bni*bnj)-1)/bni)
  lf(bkn)%xi2(jb  ) = (xi2e-xi2s)/dble(bnj) * dble((mod(bkn,bni*bnj)-1)/bni+1)
  if(mod(bkn,bni*bnj)==0)then
   lf(bkn)%xi2(js-1) = xi2s + (xi2e-xi2s)/dble(bnj)*dble((bni*bnj-1)/bni)
   lf(bkn)%xi2(jb  ) = xi2s + (xi2e-xi2s)/dble(bnj)*dble((bni*bnj-1)/bni+1)
  end if
  lf(bkn)%xi3(ks-1) = xi3s + (xi3e-xi3s)/dble(bnk) * dble((bkn-1)/(bni*bnj))
  lf(bkn)%xi3(kb  ) = xi3s + (xi3e-xi3s)/dble(bnk) * dble((bkn-1)/(bni*bnj)+1)

  lf(bkn)%dx1 = ( lf(bkn)%xi1(ib) - lf(bkn)%xi1(is-1) ) / dble(ib)
  lf(bkn)%dx2 = ( lf(bkn)%xi2(jb) - lf(bkn)%xi2(js-1) ) / dble(jb)
  lf(bkn)%dx3 = ( lf(bkn)%xi3(kb) - lf(bkn)%xi3(ks-1) ) / dble(kb)
  lf(bkn)%idx1 = 1.d0 / lf(bkn)%dx1
  lf(bkn)%idx2 = 1.d0 / lf(bkn)%dx2
  lf(bkn)%idx3 = 1.d0 / lf(bkn)%dx3

  lf(bkn)%xi1(is-2) = lf(bkn)%xi1(is-1) - lf(bkn)%dx1(is-1)
  lf(bkn)% x1(is-2) = lf(bkn)%xi1(is-2) - 5.d-1*lf(bkn)%dx1(is-2)
  do i = is-1, ib+2
   lf(bkn)%xi1(i) = lf(bkn)%xi1(i-1) + lf(bkn)%dx1(i)
   lf(bkn)% x1(i) = lf(bkn)%xi1(i)   - 5.d-1*lf(bkn)%dx1(i)
  end do
  lf(bkn)%xi2(js-2) = lf(bkn)%xi2(js-1) - lf(bkn)%dx2(js-1)
  lf(bkn)% x2(js-2) = lf(bkn)%xi2(js-2) - 5.d-1*lf(bkn)%dx2(js-2)
  do j = js-1, jb+2
   lf(bkn)%xi2(j) = lf(bkn)%xi2(j-1) + lf(bkn)%dx2(j)
   lf(bkn)% x2(j) = lf(bkn)%xi2(j)   - 5.d-1*lf(bkn)%dx2(j)
  end do
  lf(bkn)%xi3(ks-2) = lf(bkn)%xi3(ks-1) - lf(bkn)%dx3(ks-1)
  lf(bkn)% x3(ks-2) = lf(bkn)%xi3(ks-2) - 5.d-1*lf(bkn)%dx3(ks-2)
  do k = ks-1, kb+2
   lf(bkn)%xi3(k) = lf(bkn)%xi3(k-1) + lf(bkn)%dx3(k)
   lf(bkn)% x3(k) = lf(bkn)%xi3(k)   - 5.d-1*lf(bkn)%dx3(k)
  end do

  call amr_block_metric(lf(bkn))

 end do

! insert initial condition into AMR blocks
 do k = 1, bnk
  do j = 1, bnj
   do i = 1, bni
    i1 = (i-1)*ib - 1 ; i2 = i1 + ib + 1
    j1 = (j-1)*jb - 1 ; j2 = j1 + jb + 1
    k1 = (k-1)*kb - 1 ; k2 = k1 + kb + 1
    bkn = i + (j-1)*bni + (k-1)*bni*bnj
    lf(bkn)% d(is-2:ib+2,js-2:jb+2,ks-2:kb+2) =  d(i1:i2,j1:j2,k1:k2)
    lf(bkn)% p(is-2:ib+2,js-2:jb+2,ks-2:kb+2) =  p(i1:i2,j1:j2,k1:k2)
    lf(bkn)% e(is-2:ib+2,js-2:jb+2,ks-2:kb+2) =  e(i1:i2,j1:j2,k1:k2)
    lf(bkn)%v1(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = v1(i1:i2,j1:j2,k1:k2)
    lf(bkn)%v2(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = v2(i1:i2,j1:j2,k1:k2)
    lf(bkn)%v3(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = v3(i1:i2,j1:j2,k1:k2)
    lf(bkn)%b1(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = b1(i1:i2,j1:j2,k1:k2)
    lf(bkn)%b2(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = b2(i1:i2,j1:j2,k1:k2)
    lf(bkn)%b3(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = b3(i1:i2,j1:j2,k1:k2)
    lf(bkn)%ptot(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = ptot(i1:i2,j1:j2,k1:k2)
    lf(bkn)%phi(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = phi(i1:i2,j1:j2,k1:k2)
    lf(bkn)%cf(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = cf(i1:i2,j1:j2,k1:k2)
    lf(bkn)%gphi(is-2:ib+2,js-2:jb+2,ks-2:kb+2) = grvphi(i1:i2,j1:j2,k1:k2)
   end do
  end do
 end do

 do bkn = 1, totbloks
  lf(bkn)%u(:,:,:,1) = lf(bkn)%d
  lf(bkn)%u(:,:,:,2) = lf(bkn)%d*lf(bkn)%v1
  lf(bkn)%u(:,:,:,3) = lf(bkn)%d*lf(bkn)%v2
  lf(bkn)%u(:,:,:,4) = lf(bkn)%d*lf(bkn)%v3
  lf(bkn)%u(:,:,:,5) = lf(bkn)%d*lf(bkn)%b1
  lf(bkn)%u(:,:,:,6) = lf(bkn)%d*lf(bkn)%b2
  lf(bkn)%u(:,:,:,7) = lf(bkn)%d*lf(bkn)%b3
  lf(bkn)%u(:,:,:,8) = lf(bkn)%e
  lf(bkn)%u(:,:,:,9) = lf(bkn)%phi
  lf(bkn)%uorg = lf(bkn)%u
 end do
 
 rungen=2 ; sync = .true. ! for boundary reasons
 do
  totleafs = totbloks
  call amr_regridding
  if(totbloks<=totleafs)exit
 end do

 elseif(start>0)then ! for restarting

  write(startfile,'(a4,i5.5,a4)')'bina',start,'.dat'
  open(unit=10,file=startfile,status='old',form='unformatted')

  read(10)tn,time
  read(10)totbloks,curlvl
  call amr_allocateblock
  call amr_allocateleafs

  do lfn = 1, totbloks
   read(10)bk(lfn)%parnt,bk(lfn)%level,bk(lfn)%neigh,bk(lfn)%child,bk(lfn)%ldif
   read(10) lf(lfn)%x1, lf(lfn)%xi1, lf(lfn)%dx1, lf(lfn)%idx1
   read(10) lf(lfn)%x2, lf(lfn)%xi2, lf(lfn)%dx2, lf(lfn)%idx2
   read(10) lf(lfn)%x3, lf(lfn)%xi3, lf(lfn)%dx3, lf(lfn)%idx3
   read(10) lf(lfn)%d, lf(lfn)%p, lf(lfn)%e, lf(lfn)%ptot, lf(lfn)%phi
   read(10) lf(lfn)%v1, lf(lfn)%v2, lf(lfn)%v3
   read(10) lf(lfn)%b1, lf(lfn)%b2, lf(lfn)%b3
   read(10) lf(lfn)%cf, lf(lfn)%gphi, lf(lfn)%dvol

   read(10) lf(lfn)%detg1, lf(lfn)%idetg1, lf(lfn)%detg2, lf(lfn)%idetg2
   read(10) lf(lfn)%idetg3, lf(lfn)%sx1
   read(10) lf(lfn)%g22, lf(lfn)%g33
   read(10) lf(lfn)%scot, lf(lfn)%sisin
  end do

  close(10)

 else
  print *,'Error from start'
  stop
 end if

return
end subroutine amr_settings
