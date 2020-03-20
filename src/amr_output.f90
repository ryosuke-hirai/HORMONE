!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE AMR_OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output data in AMR mode.

subroutine amr_output

 use grid
 use amr_module
use amr_boundary_mod !temp

 implicit none

 character*20 pltfile, binfile, blkfile
 character*30 form1

!-----------------------------------------------------------------------------

!pltfile---------------------------------------------------------------------
 write(pltfile,'(a4,i5.5,a4)')'plta',tn,'.dat'


 open(unit=30,file = pltfile, status='replace')

 write(30,'(2(a,1PE12.4e2,2X))') '#time= ',time, "dt = ",dt
 write(30,'(4(a5),13(a13))')'# bkn','i','j','k','x1(5)','x2(6)','x3(7)','d(8)',&
      'e(9)','p(10)','v1(11)','v2(12)','v3(13)','b1(14)','b2(15)','b3(16)',&
      'phi(17)'

 do lfn = 1, totbloks
  if(bk(lfn)%child(1)==0)then
!   bkn=lfn!temp
!   call amr_boundary!temp
   do k = ks,kb
    do j = js,jb
     do i = is,ib
      write(30,'(4(i5),13(1x,1PE12.4e2))') &
           lfn,i,j,k, lf(lfn)%x1(i), lf(lfn)%x2(j), lf(lfn)%x3(k),&
           lf(lfn)% d(i,j,k), lf(lfn)% e(i,j,k), lf(lfn)% p(i,j,k), &
           lf(lfn)%v1(i,j,k), lf(lfn)%v2(i,j,k), lf(lfn)%v3(i,j,k), &
           lf(lfn)%b1(i,j,k), lf(lfn)%b2(i,j,k), lf(lfn)%b3(i,j,k), &
           lf(lfn)%gphi(i,j,k)!,lf(lfn)%dvol(i,j,k)
     end do
     if(ib/=1)write(30,*)
    end do
    !   if(ie==1)write(30,*)
   end do

  end if
 end do

 close(30)

!binfile---------------------------------------------------------------------
 write(binfile,'(a4,i5.5,a4)')'bina',tn,'.dat'

 open(unit=10,file=binfile,status='replace',form='unformatted')

 write(10)tn,time
 write(10)totbloks,curlvl
 do lfn = 1, totbloks
  write(10)bk(lfn)%parnt,bk(lfn)%level,bk(lfn)%neigh,bk(lfn)%child,bk(lfn)%ldif
  write(10) lf(lfn)%x1, lf(lfn)%xi1, lf(lfn)%dx1, lf(lfn)%idx1
  write(10) lf(lfn)%x2, lf(lfn)%xi2, lf(lfn)%dx2, lf(lfn)%idx2
  write(10) lf(lfn)%x3, lf(lfn)%xi3, lf(lfn)%dx3, lf(lfn)%idx3
  write(10) lf(lfn)%d, lf(lfn)%p, lf(lfn)%e, lf(lfn)%ptot, lf(lfn)%phi
  write(10) lf(lfn)%v1, lf(lfn)%v2, lf(lfn)%v3
  write(10) lf(lfn)%b1, lf(lfn)%b2, lf(lfn)%b3
  write(10) lf(lfn)%cf, lf(lfn)%gphi, lf(lfn)%dvol

  write(10) lf(lfn)%detg1, lf(lfn)%idetg1, lf(lfn)%detg2, lf(lfn)%idetg2
  write(10) lf(lfn)%idetg3, lf(lfn)%sx1
  write(10) lf(lfn)%g22, lf(lfn)%g33
  write(10) lf(lfn)%scot, lf(lfn)%sisin
 end do


 close(10)

!blkfile---------------------------------------------------------------------
 write(blkfile,'(a3,i5.5,a4)')'blk',tn,'.dat'

 open(unit=20,file=blkfile,status='replace')

 write(20,'(a,i3,a,i4)')'# current AMR level=',curlvl,'total blocks=',totbloks
 write(form1,'(a,i2,a)')"(a5,2a7,2X,3(2X,a",cib*5,"))"
 write(20,form1)'# bkn','level','parnt','neigh','children','ldif'

 write(form1,'(a,i1,a)')"(i5,2i7,2X,3(2X,",cib,"i5))"
 do lfn = 1, totbloks
  write(20,form1) lfn, bk(lfn)%level, bk(lfn)%parnt, &
                      bk(lfn)%neigh, bk(lfn)%child, bk(lfn)%ldif
 end do

 close(20)

return
end subroutine amr_output
