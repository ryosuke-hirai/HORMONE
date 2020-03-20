module amr_source_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE AMR_SOURCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate source terms for a block in AMR mode.

subroutine amr_source(lf)

 use settings,only:crdnt
 use grid,only:is,js,ks
 use settings,only:gravswitch
 use amr_templates
 use amr_module,only:ib,jb,kb
 use minmod_mod

 implicit none

 integer i,j,k
 real*8 u1, u2, u3, u4, phil, phir, ul, ur
 type(leaf_contents),intent(inout):: lf
 real*8,dimension(is:ib,js:jb,ks:kb):: grv1, grv2, grv3

!-----------------------------------------------------------------------------

! To calculate gravitational forces *******************************************
!!$ if(gravswitch==0)then
!!$  grv1 = 0.d0 ; grv2 = 0.d0 ; grv3 = 0.d0
!!$ else
!!$!$omp parallel do private(i,j,k,ul,ur,u1,u2,u3,u4,phil,phir)
!!$  do k = ks,kb
!!$   do j = js,jb
!!$    do i = is, ib
!!$     u3=lf%gphi(i,j,k)
!!$     u1=lf%gphi(i-2,j,k); u2=lf%gphi(i-1,j,k); u4=lf%gphi(i+1,j,k)
!!$     call minmod(ul,u1,u2,u3); call minmod(ur,u2,u3,u4)
!!$     ul = u2 + 5.d-1*ul ; ur = u3 - 5.d-1*ur ; phil = 5.d-1 * (ul+ur)
!!$     u1 = lf%gphi(i+2,j,k)
!!$     call minmod(ul,u2,u3,u4); call minmod(ur,u3,u4,u1)
!!$     ul = u3 + 5.d-1*ul ; ur = u4 - 5.d-1*ur ; phir = 5.d-1 * (ul+ur)
!!$
!!$     grv1(i,j,k) = - (phir-phil) * lf%idx1(i) * lf%d(i,j,k)
!!$
!!$     u1=lf%gphi(i,j-2,k); u2=lf%gphi(i,j-1,k); u4=lf%gphi(i,j+1,k)
!!$     call minmod(ul,u1,u2,u3); call minmod(ur,u2,u3,u4)
!!$     ul = u2 + 5.d-1*ul ; ur = u3 - 5.d-1*ur ; phil = 5.d-1 * (ul+ur)
!!$     u1 = lf%gphi(i,j+2,k)
!!$     call minmod(ul,u2,u3,u4); call minmod(ur,u3,u4,u1)
!!$     ul = u3 + 5.d-1*ul ; ur = u4 - 5.d-1*ur ; phir = 5.d-1 * (ul+ur)
!!$
!!$     grv2(i,j,k) = - (phir-phil) * lf%idx2(j) * lf%d(i,j,k)
!!$
!!$     u1=lf%gphi(i,j,k-2); u2=lf%gphi(i,j,k-1); u4=lf%gphi(i,j,k+1)
!!$     call minmod(ul,u1,u2,u3); call minmod(ur,u2,u3,u4)
!!$     ul = u2 + 5.d-1*ul ; ur = u3 - 5.d-1*ur ; phil = 5.d-1 * (ul+ur)
!!$     u1 = lf%gphi(i,j,k+2)
!!$     call minmod(ul,u2,u3,u4); call minmod(ur,u3,u4,u1)
!!$     ul = u3 + 5.d-1*ul ; ur = u4 - 5.d-1*ur ; phir = 5.d-1 * (ul+ur)
!!$
!!$     grv3(i,j,k) = - (phir-phil) * lf%idx3(k) * lf%d(i,j,k)
!!$     !print *,grv1(i,j,k)
!!$    end do
!!$   end do
!!$  end do
!!$!$omp end parallel do
!!$ end if
!!$
!!$
!!$
!!$! Calculate source terms ****************************************************
!!$
!!$! Cartesian >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$ if(crdnt==0)then
!!$  do k = ks, kb
!!$   do j = js, jb
!!$    do i = is, ib
!!$     lf%src(i,j,k,1) = 0.d0
!!$     lf%src(i,j,k,2) = grv1(i,j,k)
!!$     lf%src(i,j,k,3) = grv2(i,j,k)
!!$     lf%src(i,j,k,4) = grv3(i,j,k)
!!$     lf%src(i,j,k,5:9) = 0.d0
!!$    end do
!!$   end do
!!$  end do
!!$
!!$! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$ elseif(crdnt==1)then
!!$  do k = ks, kb
!!$   do j = js, jb
!!$    do i = is, ib
!!$     lf%src(i,j,k,1) = 0.d0
!!$
!!$     lf%src(i,j,k,2) = ( lf%d(i,j,k) * lf%v2(i,j,k)*lf%v2(i,j,k) &
!!$                                     - lf%b2(i,j,k)*lf%b2(i,j,k) &
!!$                       + lf%ptot(i,j,k) ) * lf%sx1(i) &
!!$                     + grv1(i,j,k)
!!$
!!$     lf%src(i,j,k,3) = (-lf%d(i,j,k) * lf%v1(i,j,k)*lf%v2(i,j,k)   &
!!$                                     + lf%b1(i,j,k)*lf%b2(i,j,k) ) &
!!$                       * lf%sx1(i) &
!!$                     + grv2(i,j,k)
!!$
!!$     lf%src(i,j,k,4:5) = 0.d0
!!$
!!$     lf%src(i,j,k,6) = ( lf%b1(i,j,k)*lf%v3(i,j,k) &
!!$                       - lf%b3(i,j,k)*lf%v1(i,j,k) ) * lf%sx1(i) &
!!$                     + grv3(i,j,k)
!!$
!!$     lf%src(i,j,k,7:9) = 0.d0
!!$
!!$    end do
!!$   end do
!!$  end do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$ elseif(crdnt==2)then
!!$  do k = ks, kb
!!$   do j = js, jb
!!$    do i = is, ib
!!$     lf%src(i,j,k,1) = 0.d0
!!$
!!$     lf%src(i,j,k,2) = ( lf%d(i,j,k)*( lf%v2(i,j,k)*lf%v2(i,j,k)   &
!!$                                     + lf%v3(i,j,k)*lf%v3(i,j,k) ) &
!!$                     - ( lf%b2(i,j,k)*lf%b2(i,j,k) + lf%b3(i,j,k)*lf%b3(i,j,k))&
!!$                       + 2.d0*lf%ptot(i,j,k) ) * lf%sx1(i) &
!!$                     + grv1(i,j,k)
!!$
!!$     lf%src(i,j,k,3) = ( lf%b1(i,j,k)*lf%b2(i,j,k) &
!!$                       - lf%d(i,j,k)*lf%v1(i,j,k)*lf%v2(i,j,k) &
!!$                     + ( lf%d(i,j,k)*lf%v3(i,j,k)*lf%v3(i,j,k) &
!!$                       - lf%b3(i,j,k)*lf%b3(i,j,k) + lf%ptot(i,j,k) ) &
!!$                       * lf%scot(j) ) * lf%sx1(i) &
!!$                     + grv2(i,j,k)*lf%sx1(i)
!!$
!!$     lf%src(i,j,k,4) = ( lf%b1(i,j,k)*lf%b3(i,j,k) &
!!$                       - lf%d(i,j,k)*lf%v1(i,j,k)*lf%v3(i,j,k) &
!!$                       + ( lf%b2(i,j,k)*lf%b3(i,j,k) &
!!$                         - lf%d(i,j,k)*lf%v2(i,j,k)*lf%v3(i,j,k) ) &
!!$                         * lf%scot(j) ) * lf%sx1(i) &
!!$                     + grv3(i,j,k)*lf%sx1(i)*lf%sisin(j)
!!$
!!$     lf%src(i,j,k,5) = 0.d0
!!$
!!$     lf%src(i,j,k,6) = ( lf%b2(i,j,k)*lf%v1(i,j,k) &
!!$                       - lf%b1(i,j,k)*lf%v2(i,j,k) ) * lf%sx1(i)
!!$
!!$     lf%src(i,j,k,7) = ( ( lf%v2(i,j,k)*lf%b3(i,j,k) &
!!$                         - lf%b2(i,j,k)*lf%v3(i,j,k) ) * lf%scot(j) &
!!$                     - ( lf%b1(i,j,k)*lf%v3(i,j,k) &
!!$                       - lf%b3(i,j,k)*lf%v1(i,j,k)) ) * lf%sx1(i)
!!$
!!$     lf%src(i,j,k,8) = grv1(i,j,k)*lf%v1(i,j,k) &
!!$                     + grv2(i,j,k)*lf%v2(i,j,k) * lf%sx1(i) &
!!$                     + grv3(i,j,k)*lf%v3(i,j,k) * lf%sx1(i)*lf%sisin(j)
!!$
!!$     lf%src(i,j,k,9) = 0.d0
!!$
!!$    end do
!!$   end do
!!$  end do
!!$ else
!!$  print *,'Error from amr_source.f90, crdnt'
!!$  stop
!!$ end if


return
end subroutine amr_source

end module amr_source_mod
