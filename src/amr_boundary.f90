module amr_bound_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_BOUND_SAME
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundaries for blocks with same level neighbours.

subroutine amr_bound_same(n,lfu,neiu)

 use grid,only:is,js,ks
 use amr_module,only:ib,jb,kb

 implicit none

 integer,intent(in):: n
 real*8,allocatable,intent(inout),dimension(:,:,:):: lfu
 real*8,allocatable,intent(in   ),dimension(:,:,:):: neiu

!-----------------------------------------------------------------------------

! x1 direction
 if((n+1)/2==1)then
  if(mod(n,2)==1)then ! inner boundary
   lfu(is-2:is-1,js-1:jb+1,ks-1:kb+1) = neiu(ib-1:ib,js-1:jb+1,ks-1:kb+1)
  else                ! outer boundary
   lfu(ib+1:ib+2,js-1:jb+1,ks-1:kb+1) = neiu(is:is+1,js-1:jb+1,ks-1:kb+1)
  end if
! x2 direction
 elseif((n+1)/2==2)then
  if(mod(n,2)==1)then ! inner boundary
   lfu(is-1:ib+1,js-2:js-1,ks-1:kb+1) = neiu(is-1:ib+1,jb-1:jb,ks-1:kb+1)
  else                ! outer boundary
   lfu(is-1:ib+1,jb+1:jb+2,ks-1:kb+1) = neiu(is-1:ib+1,js:js+1,ks-1:kb+1)
  end if
! x3 direction
 elseif((n+1)/2==3)then
  if(mod(n,2)==1)then ! inner boundary
   lfu(is-1:ib+1,js-1:jb+1,ks-2:ks-1) = neiu(is-1:ib+1,js-1:jb+1,kb-1:kb)
  else                ! outer boundary
   lfu(is-1:ib+1,js-1:jb+1,kb+1:kb+2) = neiu(is-1:ib+1,js-1:jb+1,ks:ks+1)
  end if

 else
  print *,"Error from boundary, amr_bound_same"
  stop
 end if

end subroutine amr_bound_same

!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$!
!!$!                        SUBROUTINE AMR_BOUND_FINE
!!$!
!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$
!!$! PURPOSE: To set boundary condition for finer neighbours.
!!$
!!$subroutine amr_bound_fine(n,lf,blf)
!!$
!!$ use grid,only:is,js,ks
!!$ use physval,only:gamma
!!$ use amr_module,only:ib,jb,kb,cib,cuts
!!$ use amr_templates
!!$
!!$ implicit none
!!$
!!$ integer,intent(in):: n
!!$ integer i, j, k, lid,i1,i2,j1,j2,k1,k2
!!$ type(leaf_contents),intent(in),dimension(1:cib/cuts):: blf
!!$ type(leaf_contents),intent(inout):: lf
!!$ real*8 vsq, bsq
!!$
!!$!-----------------------------------------------------------------------------
!!$
!!$!x1 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$ if((n+1)/2==1)then
!!$
!!$   do k = ks, kb
!!$    do j = js, jb
!!$     lid = 1 + (j-1)*jb/cuts + (k-1)*kb/cuts
!!$     j1 = mod(j,max(1,jb/cuts))
!!$     j1 = cuts*j1 + (1-sign(1,j1-1))/2*jb - cuts + 1
!!$     j1 = (1+sign(1,j1))/2 *  j1         + (1-sign(1,j1))/2
!!$     j2 = (1+sign(1,j1))/2 * (j1+cuts-1) + (1-sign(1,j1))/2
!!$     k1 = mod(k,max(1,kb/cuts))
!!$     k1 = cuts*k1 + (1-sign(1,k1-1))/2*kb - cuts + 1
!!$     k1 = (1+sign(1,k1))/2 *  k1         + (1-sign(1,k1))/2
!!$     k2 = (1+sign(1,k1))/2 * (k1+cuts-1) + (1-sign(1,k1))/2
!!$print *,j1,j2,k1,k2
!!$     if(mod(n,2)==1)then ! inner boundary
!!$      do i = 2,1,-1
!!$       i1 = ib-i*cuts+1 ; i2 = i1-1 + cuts
!!$
!!$       lf% d(is-i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k)
!!$
!!$       lf%v1(is-i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k) / lf%d(is-i,j,k)
!!$
!!$       lf%v2(is-i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k) / lf%d(is-i,j,k)
!!$
!!$       lf%v3(is-i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k) / lf%d(is-i,j,k)
!!$
!!$       vsq = lf%v1(is-i,j,k)*lf%v1(is-i,j,k) &
!!$           + lf%v2(is-i,j,k)*lf%v2(is-i,j,k) &
!!$           + lf%v3(is-i,j,k)*lf%v3(is-i,j,k)
!!$
!!$       lf%b1(is-i,j,k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k)
!!$
!!$       lf%b2(is-i,j,k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k)
!!$
!!$       lf%b3(is-i,j,k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k)
!!$
!!$       bsq = lf%b1(is-i,j,k)*lf%b1(is-i,j,k) &
!!$           + lf%b2(is-i,j,k)*lf%b2(is-i,j,k) &
!!$           + lf%b3(is-i,j,k)*lf%b3(is-i,j,k)
!!$
!!$       lf% e(is-i,j,k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(is-i,j,k)
!!$
!!$       lf%phi(is-i,j,k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(is-i,j,k)
!!$
!!$       lf% p(is-i,j,k) = (gamma-1.d0) &
!!$                       * ( lf%e(is-i,j,k) - 5.d-1*lf%d(is-i,j,k)*vsq-5.d-1*bsq )
!!$       lf%ptot(is-i,j,k) = lf%p(is-i,j,k) + 5.d-1*bsq
!!$      end do
!!$     else ! outer boundary
!!$      do i = 1,2
!!$       i1 = 1 + (i-1)*cuts ; i2 = i1 + cuts-1
!!$
!!$       lf% d(ib+i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k)
!!$
!!$       lf%v1(ib+i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k) / lf%d(ib+i,j,k)
!!$
!!$       lf%v2(ib+i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k) / lf%d(ib+i,j,k)
!!$
!!$       lf%v3(ib+i,j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k) / lf%d(ib+i,j,k)
!!$
!!$       vsq = lf%v1(ib+i,j,k)*lf%v1(ib+i,j,k) &
!!$           + lf%v2(ib+i,j,k)*lf%v2(ib+i,j,k) &
!!$           + lf%v3(ib+i,j,k)*lf%v3(ib+i,j,k)
!!$
!!$       lf%b1(ib+i,j,k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k)
!!$
!!$       lf%b2(ib+i,j,k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k)
!!$
!!$       lf%b3(ib+i,j,k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k)
!!$
!!$       bsq = lf%b1(ib+i,j,k)*lf%b1(ib+i,j,k) &
!!$           + lf%b2(ib+i,j,k)*lf%b2(ib+i,j,k) &
!!$           + lf%b3(ib+i,j,k)*lf%b3(ib+i,j,k)
!!$
!!$       lf% e(ib+i,j,k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(ib+i,j,k)
!!$
!!$       lf%phi(ib+i,j,k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(ib+i,j,k)
!!$
!!$       lf% p(ib+i,j,k) = (gamma-1.d0) &
!!$                       * ( lf%e(ib+i,j,k) - 5.d-1*lf%d(ib+i,j,k)*vsq-5.d-1*bsq )
!!$       lf%ptot(ib+i,j,k) = lf%p(ib+i,j,k) + 5.d-1*bsq
!!$      end do
!!$     end if
!!$    end do
!!$   end do
!!$
!!$! x2 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  elseif((n+1)/2==2)then
!!$
!!$   do k = ks, kb
!!$    do i = is, ib
!!$     lid = 1 + (i-1)*ib/cuts + (k-1)*kb/cuts
!!$     i1 = mod(i,max(1,ib/cuts))
!!$     i1 = cuts*i1 + (1-sign(1,i1-1))/2*ib - cuts + 1
!!$     i1 = (1+sign(1,i1))/2 *  i1         + (1-sign(1,i1))/2
!!$     i2 = (1+sign(1,i1))/2 * (i1+cuts-1) + (1-sign(1,i1))/2
!!$     k1 = mod(k,max(1,kb/cuts))
!!$     k1 = cuts*k1 + (1-sign(1,k1-1))/2*kb - cuts + 1
!!$     k1 = (1+sign(1,k1))/2 *  k1         + (1-sign(1,k1))/2
!!$     k2 = (1+sign(1,k1))/2 * (k1+cuts-1) + (1-sign(1,k1))/2
!!$
!!$     if(mod(n,2)==1)then ! inner boundary
!!$      do j = 2,1,-1
!!$       j1 = jb-j*cuts+1 ; j2 = j1-1 + cuts
!!$
!!$       lf% d(i,js-j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k)
!!$
!!$       lf%v1(i,js-j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k) / lf%d(i,js-j,k)
!!$
!!$       lf%v2(i,js-j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k) / lf%d(i,js-j,k)
!!$
!!$       lf%v3(i,js-j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k) / lf%d(i,js-j,k)
!!$
!!$       vsq = lf%v1(i,js-j,k)*lf%v1(i,js-j,k) &
!!$           + lf%v2(i,js-j,k)*lf%v2(i,js-j,k) &
!!$           + lf%v3(i,js-j,k)*lf%v3(i,js-j,k)
!!$
!!$       lf%b1(i,js-j,k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k)
!!$
!!$       lf%b2(i,js-j,k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k)
!!$
!!$       lf%b3(i,js-j,k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k)
!!$
!!$       bsq = lf%b1(i,js-j,k)*lf%b1(i,js-j,k) &
!!$           + lf%b2(i,js-j,k)*lf%b2(i,js-j,k) &
!!$           + lf%b3(i,js-j,k)*lf%b3(i,js-j,k)
!!$
!!$       lf% e(i,js-j,k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,js-j,k)
!!$
!!$       lf%phi(i,js-j,k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(i,js-j,k)
!!$
!!$       lf% p(i,js-j,k) = (gamma-1.d0) &
!!$                       * ( lf%e(i,js-j,k) - 5.d-1*lf%d(i,js-j,k)*vsq-5.d-1*bsq )
!!$       lf%ptot(i,js-j,k) = lf%p(i,js-j,k) + 5.d-1*bsq
!!$      end do
!!$     else ! outer boundary
!!$      do j = 1,2
!!$       j1 = 1 + (j-1)*cuts ; j2 = j1 + cuts-1
!!$
!!$       lf% d(i,jb+j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k)
!!$
!!$       lf%v1(i,jb+j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k) / lf%d(i,jb+j,k)
!!$
!!$       lf%v2(i,jb+j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k) / lf%d(i,jb+j,k)
!!$
!!$       lf%v3(i,jb+j,k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k) / lf%d(i,jb+j,k)
!!$
!!$       vsq = lf%v1(i,jb+j,k)*lf%v1(i,jb+j,k) &
!!$           + lf%v2(i,jb+j,k)*lf%v2(i,jb+j,k) &
!!$           + lf%v3(i,jb+j,k)*lf%v3(i,jb+j,k)
!!$
!!$       lf%b1(i,jb+j,k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k)
!!$
!!$       lf%b2(i,jb+j,k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k)
!!$
!!$       lf%b3(i,jb+j,k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k)
!!$
!!$       bsq = lf%b1(i,jb+j,k)*lf%b1(i,jb+j,k) &
!!$           + lf%b2(i,jb+j,k)*lf%b2(i,jb+j,k) &
!!$           + lf%b3(i,jb+j,k)*lf%b3(i,jb+j,k)
!!$
!!$       lf% e(i,jb+j,k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,jb+j,k)
!!$
!!$       lf%phi(i,jb+j,k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(i,jb+j,k)
!!$
!!$       lf% p(i,jb+j,k) = (gamma-1.d0) &
!!$                       * ( lf%e(i,jb+j,k) - 5.d-1*lf%d(i,jb+j,k)*vsq-5.d-1*bsq )
!!$       lf%ptot(i,jb+j,k) = lf%p(i,jb+j,k) + 5.d-1*bsq
!!$      end do
!!$     end if
!!$    end do
!!$   end do
!!$
!!$! x3 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$  elseif((n+1)/2==3)then
!!$
!!$   do j = js, jb
!!$    do i = is, ib
!!$     lid = 1 + (i-1)*ib/cuts + (j-1)*jb/cuts
!!$     i1 = mod(i,max(1,ib/cuts))
!!$     i1 = cuts*i1 + (1-sign(1,i1-1))/2*ib - cuts + 1
!!$     i1 = (1+sign(1,i1))/2 *  i1         + (1-sign(1,i1))/2
!!$     i2 = (1+sign(1,i1))/2 * (i1+cuts-1) + (1-sign(1,i1))/2
!!$     j1 = mod(j,max(1,jb/cuts))
!!$     j1 = cuts*j1 + (1-sign(1,j1-1))/2*jb - cuts + 1
!!$     j1 = (1+sign(1,j1))/2 *  j1         + (1-sign(1,j1))/2
!!$     j2 = (1+sign(1,j1))/2 * (j1+cuts-1) + (1-sign(1,j1))/2
!!$
!!$     if(mod(n,2)==1)then ! inner boundary
!!$      do k = 2,1,-1
!!$       k1 = kb-k*cuts+1 ; k2 = k1-1 + cuts
!!$
!!$       lf% d(i,j,ks-k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k)
!!$
!!$       lf%v1(i,j,ks-k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k) / lf%d(i,j,ks-k)
!!$
!!$       lf%v2(i,j,ks-k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k) / lf%d(i,j,ks-k)
!!$
!!$       lf%v3(i,j,ks-k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k) / lf%d(i,j,ks-k)
!!$
!!$       vsq = lf%v1(i,j,ks-k)*lf%v1(i,j,ks-k) &
!!$           + lf%v2(i,j,ks-k)*lf%v2(i,j,ks-k) &
!!$           + lf%v3(i,j,ks-k)*lf%v3(i,j,ks-k)
!!$
!!$       lf%b1(i,j,ks-k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k)
!!$
!!$       lf%b2(i,j,ks-k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k)
!!$
!!$       lf%b3(i,j,ks-k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k)
!!$
!!$       bsq = lf%b1(i,j,ks-k)*lf%b1(i,j,ks-k) &
!!$           + lf%b2(i,j,ks-k)*lf%b2(i,j,ks-k) &
!!$           + lf%b3(i,j,ks-k)*lf%b3(i,j,ks-k)
!!$
!!$       lf% e(i,j,ks-k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,ks-k)
!!$
!!$       lf%phi(i,j,ks-k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(i,j,ks-k)
!!$
!!$       lf% p(i,j,ks-k) = (gamma-1.d0) &
!!$                       * ( lf%e(i,j,ks-k) - 5.d-1*lf%d(i,j,ks-k)*vsq-5.d-1*bsq )
!!$       lf%ptot(i,j,ks-k) = lf%p(i,j,ks-k) + 5.d-1*bsq
!!$      end do
!!$     else ! outer boundary
!!$      do k = 1,2
!!$       k1 = 1 + (k-1)*cuts ; k2 = k1 + cuts-1
!!$
!!$       lf% d(i,j,kb+k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k)
!!$
!!$       lf%v1(i,j,kb+k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k) / lf%d(i,j,kb+k)
!!$
!!$       lf%v2(i,j,kb+k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k) / lf%d(i,j,kb+k)
!!$
!!$       lf%v3(i,j,kb+k) = sum( blf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k) / lf%d(i,j,kb+k)
!!$
!!$       vsq = lf%v1(i,j,kb+k)*lf%v1(i,j,kb+k) &
!!$           + lf%v2(i,j,kb+k)*lf%v2(i,j,kb+k) &
!!$           + lf%v3(i,j,kb+k)*lf%v3(i,j,kb+k)
!!$
!!$       lf%b1(i,j,kb+k) = sum( blf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k)
!!$
!!$       lf%b2(i,j,kb+k) = sum( blf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k)
!!$
!!$       lf%b3(i,j,kb+k) = sum( blf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
!!$                            * blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k)
!!$
!!$       bsq = lf%b1(i,j,kb+k)*lf%b1(i,j,kb+k) &
!!$           + lf%b2(i,j,kb+k)*lf%b2(i,j,kb+k) &
!!$           + lf%b3(i,j,kb+k)*lf%b3(i,j,kb+k)
!!$
!!$       lf% e(i,j,kb+k) = sum( blf(lid)%   e(i1:i2,j1:j2,k1:k2) &
!!$                           *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                       / lf%dvol(i,j,kb+k)
!!$
!!$       lf%phi(i,j,kb+k) = sum( blf(lid)% phi(i1:i2,j1:j2,k1:k2) &
!!$                            *  blf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
!!$                        / lf%dvol(i,j,kb+k)
!!$
!!$       lf% p(i,j,kb+k) = (gamma-1.d0) &
!!$                       * ( lf%e(i,j,kb+k) - 5.d-1*lf%d(i,j,kb+k)*vsq-5.d-1*bsq )
!!$       lf%ptot(i,j,kb+k) = lf%p(i,j,kb+k) + 5.d-1*bsq
!!$      end do
!!$     end if
!!$    end do
!!$   end do
!!$
!!$  else
!!$   print *,"Error from amr_bound_fine",n
!!$   stop
!!$  end if
!!$
!!$return
!!$end subroutine amr_bound_fine


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_BOUND_COARSE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set guard cells for coarser neighbours.

subroutine amr_bound_coarse(n,lid,lf,nlf)

 use grid,only:is,js,ks,rungen
 use physval,only:gamma
 use amr_templates
 use amr_module,only:ib,jb,kb,cib,face,cuts,sync,bkn

 implicit none

 integer,intent(in):: n,lid
 integer i1,i2,j1,j2,k1,k2
 real*8,allocatable,dimension(:,:,:):: vsq, bsq
 type(leaf_contents),intent(in):: nlf
 type(leaf_contents),intent(inout):: lf

 real*8,allocatable,dimension(:,:,:,:):: nu

!-----------------------------------------------------------------------------

 if(cuts/=2)then
  print *,"Error: Boundary scheme for cuts > 2 is needed",cuts
  stop
 end if

 allocate( bsq(-1:ib+2,-1:jb+2,-1:kb+2),vsq(-1:ib+2,-1:jb+2,-1:kb+2) )

!x1 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
allocate( nu(-1:1,0:max(jb/cuts+1,2),0:max(kb/cuts+1,2),1:9) )

if((n+1)/2==1)then
 j1 = jb/cuts * (mod(lid-1,cuts*cuts)/cuts) ; j2 = j1 + max(jb/cuts+1,2)
 k1 = kb/cuts * ((lid - 1)/(cuts*cuts))     ; k2 = k1 + max(kb/cuts+1,2) 

 if(mod(n,2)==1)then ! inner boundary
  i1 = is-2 ; i2 = is-1

  if(sync)then
   if    (rungen==1)then ; nu = nlf%midu1(1:3,j1:j2,k1:k2,1:9)
   elseif(rungen==2)then ; nu = nlf%u(ib-1:ib+1,j1:j2,k1:k2,1:9)
   else
    nu = (nlf%u(ib-1:ib+1,j1:j2,k1:k2,1:9)+nlf%midu1(1:3,j1:j2,k1:k2,1:9))*5.d-1
   end if
  else
   if    (rungen==1)then ; nu = nlf%uorg(ib-1:ib+1,j1:j2,k1:k2,1:9)
   elseif(rungen==2)then ; nu = nlf%midu1(1:3,j1:j2,k1:k2,1:9)
   else
    nu = (nlf%uorg(ib-1:ib+1,j1:j2,k1:k2,1:9)+nlf%midu1(1:3,j1:j2,k1:k2,1:9))*5.d-1
   end if
  end if

  call amr_bound_intpol_x1(i1,nu,lf%u,lf%dvol)
  lf%  d(i1:i2,:,:) = lf%u(i1:i2,:,:,1)
  lf% v1(i1:i2,:,:) = lf%u(i1:i2,:,:,2) / lf%u(i1:i2,:,:,1)
  lf% v2(i1:i2,:,:) = lf%u(i1:i2,:,:,3) / lf%u(i1:i2,:,:,1)
  lf% v3(i1:i2,:,:) = lf%u(i1:i2,:,:,4) / lf%u(i1:i2,:,:,1)
  lf% b1(i1:i2,:,:) = lf%u(i1:i2,:,:,5)
  lf% b2(i1:i2,:,:) = lf%u(i1:i2,:,:,6)
  lf% b3(i1:i2,:,:) = lf%u(i1:i2,:,:,7)
  lf%  e(i1:i2,:,:) = lf%u(i1:i2,:,:,8)
  lf%phi(i1:i2,:,:) = lf%u(i1:i2,:,:,9)

! pressure and total pressure
  vsq(i1:i2,:,:) = lf%v1(i1:i2,:,:)*lf%v1(i1:i2,:,:) &
                 + lf%v2(i1:i2,:,:)*lf%v2(i1:i2,:,:) &
                 + lf%v3(i1:i2,:,:)*lf%v3(i1:i2,:,:)
  bsq(i1:i2,:,:) = lf%b1(i1:i2,:,:)*lf%b1(i1:i2,:,:) &
                 + lf%b2(i1:i2,:,:)*lf%b2(i1:i2,:,:) &
                 + lf%b3(i1:i2,:,:)*lf%b3(i1:i2,:,:)
  lf%p(i1:i2,:,:) = (gamma-1.d0) &
                  * (lf%e(i1:i2,:,:) - 0.5d0*lf%d(i1:i2,:,:) * vsq(i1:i2,:,:) &
                                     - 0.5d0*                  bsq(i1:i2,:,:) )
  lf%ptot(i1:i2,:,:) = lf%p(i1:i2,:,:) + 5.d-1*bsq(i1:i2,:,:)

 else ! outer boundary
  i1 = ib+1 ; i2 = ib+2

  if(sync)then
   if    (rungen==1)then ; nu = nlf%midu1(-2:0,j1:j2,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%u(is-1:is+1,j1:j2,k1:k2,:)
   else
    nu = (nlf%u(is-1:is+1,j1:j2,k1:k2,:)+nlf%midu1(-2:0,j1:j2,k1:k2,:))*5.d-1
   end if
  else
   if    (rungen==1)then ; nu = nlf%uorg(is-1:is+1,j1:j2,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%midu1(-2:0,j1:j2,k1:k2,:)
   else
    nu = (nlf%uorg(is-1:is+1,j1:j2,k1:k2,:)+nlf%midu1(-2:0,j1:j2,k1:k2,:))*5.d-1
   end if
  end if

  call amr_bound_intpol_x1(i1,nu,lf%u,lf%dvol)
  lf%  d(i1:i2,:,:) = lf%u(i1:i2,:,:,1)
  lf% v1(i1:i2,:,:) = lf%u(i1:i2,:,:,2) / lf%u(i1:i2,:,:,1)
  lf% v2(i1:i2,:,:) = lf%u(i1:i2,:,:,3) / lf%u(i1:i2,:,:,1)
  lf% v3(i1:i2,:,:) = lf%u(i1:i2,:,:,4) / lf%u(i1:i2,:,:,1)
  lf% b1(i1:i2,:,:) = lf%u(i1:i2,:,:,5)
  lf% b2(i1:i2,:,:) = lf%u(i1:i2,:,:,6)
  lf% b3(i1:i2,:,:) = lf%u(i1:i2,:,:,7)
  lf%  e(i1:i2,:,:) = lf%u(i1:i2,:,:,8)
  lf%phi(i1:i2,:,:) = lf%u(i1:i2,:,:,9)

! pressure and total pressure
  vsq(i1:i2,:,:) = lf%v1(i1:i2,:,:)*lf%v1(i1:i2,:,:) &
                 + lf%v2(i1:i2,:,:)*lf%v2(i1:i2,:,:) &
                 + lf%v3(i1:i2,:,:)*lf%v3(i1:i2,:,:)
  bsq(i1:i2,:,:) = lf%b1(i1:i2,:,:)*lf%b1(i1:i2,:,:) &
                 + lf%b2(i1:i2,:,:)*lf%b2(i1:i2,:,:) &
                 + lf%b3(i1:i2,:,:)*lf%b3(i1:i2,:,:)
  lf%p(i1:i2,:,:) = (gamma-1.d0) &
                  * (lf%e(i1:i2,:,:) - 0.5d0*lf%d(i1:i2,:,:) * vsq(i1:i2,:,:) &
                                     - 0.5d0*                  bsq(i1:i2,:,:) )
  lf%ptot(i1:i2,:,:) = lf%p(i1:i2,:,:) + 5.d-1*bsq(i1:i2,:,:)

 end if

! x2 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
deallocate(nu)
allocate( nu(0:max(ib/cuts+1,2),-1:1,0:max(kb/cuts+1,2),1:9) )

elseif((n+1)/2==2)then
 i1 = ib/cuts * mod(lid+1,cuts)       ; i2 = i1 + max(ib/cuts+1,2)
 k1 = kb/cuts * (lid - 1)/(cuts*cuts) ; k2 = k1 + max(kb/cuts+1,2)

 if(mod(n,2)==1)then ! inner boundary
  j1 = js-2 ; j2 = js-1
  if(sync)then
   if(rungen==1)then     ; nu = nlf%midu2(i1:i2,1:3,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%u(i1:i2,jb-1:jb+1,k1:k2,:)
   else
    nu = (nlf%u(i1:i2,jb-1:jb+1,k1:k2,:)+nlf%midu2(i1:i2,1:3,k1:k2,:))*5.d-1
   end if
  else
   if(rungen==1)then     ; nu = nlf%uorg(i1:i2,jb-1:jb+1,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%midu2(i1:i2,1:3,k1:k2,:)
   else
    nu = (nlf%uorg(i1:i2,jb-1:jb+1,k1:k2,:)+nlf%midu2(i1:i2,1:3,k1:k2,:))*5.d-1
   end if
  end if
  call amr_bound_intpol_x2(j1,nu,lf%u,lf%dvol)
  lf% d(:,j1:j2,:) = lf%u(:,j1:j2,:,1)
  lf%v1(:,j1:j2,:) = lf%u(:,j1:j2,:,2) / lf%u(:,j1:j2,:,1)
  lf%v2(:,j1:j2,:) = lf%u(:,j1:j2,:,3) / lf%u(:,j1:j2,:,1)
  lf%v3(:,j1:j2,:) = lf%u(:,j1:j2,:,4) / lf%u(:,j1:j2,:,1)
  lf%b1(:,j1:j2,:) = lf%u(:,j1:j2,:,5)
  lf%b2(:,j1:j2,:) = lf%u(:,j1:j2,:,6)
  lf%b3(:,j1:j2,:) = lf%u(:,j1:j2,:,7)
  lf% e(:,j1:j2,:) = lf%u(:,j1:j2,:,8)
  lf%phi(:,j1:j2,:) = lf%u(:,j1:j2,:,9)

! pressure and total pressure
  vsq(:,j1:j2,:) = lf%v1(:,j1:j2,:)*lf%v1(:,j1:j2,:) &
                 + lf%v2(:,j1:j2,:)*lf%v2(:,j1:j2,:) &
                 + lf%v3(:,j1:j2,:)*lf%v3(:,j1:j2,:)
  bsq(:,j1:j2,:) = lf%b1(:,j1:j2,:)*lf%b1(:,j1:j2,:) &
                 + lf%b2(:,j1:j2,:)*lf%b2(:,j1:j2,:) &
                 + lf%b3(:,j1:j2,:)*lf%b3(:,j1:j2,:)
  lf%p(:,j1:j2,:) = (gamma-1.d0) &
                  * (lf%e(:,j1:j2,:) - 0.5d0*lf%d(:,j1:j2,:) * vsq(:,j1:j2,:) &
                                     - 0.5d0*                  bsq(:,j1:j2,:) )
  lf%ptot(:,j1:j2,:) = lf%p(:,j1:j2,:) + 5.d-1 * bsq(:,j1:j2,:)

  else ! outer boundary
  j1 = jb+1 ; j2 = jb+2
  if(sync)then
   if(rungen==1)then     ; nu = nlf%midu2(i1:i2,-2:0,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%u(i1:i2,js-1:js+1,k1:k2,:)
   else
    nu = (nlf%u(i1:i2,js-1:js+1,k1:k2,:)+nlf%midu2(i1:i2,-2:0,k1:k2,:))*5.d-1
   end if
  else
   if(rungen==1)then     ; nu = nlf%uorg(i1:i2,js-1:js+1,k1:k2,:)
   elseif(rungen==2)then ; nu = nlf%midu2(i1:i2,-2:0,k1:k2,:)
   else
    nu = (nlf%uorg(i1:i2,js-1:js+1,k1:k2,:)+nlf%midu2(i1:i2,-2:0,k1:k2,:))*5.d-1
   end if
  end if
  call amr_bound_intpol_x2(j1,nu,lf%u,lf%dvol)
  lf%  d(:,j1:j2,:) = lf%u(:,j1:j2,:,1)
  lf% v1(:,j1:j2,:) = lf%u(:,j1:j2,:,2) / lf%u(:,j1:j2,:,1)
  lf% v2(:,j1:j2,:) = lf%u(:,j1:j2,:,3) / lf%u(:,j1:j2,:,1)
  lf% v3(:,j1:j2,:) = lf%u(:,j1:j2,:,4) / lf%u(:,j1:j2,:,1)
  lf% b1(:,j1:j2,:) = lf%u(:,j1:j2,:,5)
  lf% b2(:,j1:j2,:) = lf%u(:,j1:j2,:,6)
  lf% b3(:,j1:j2,:) = lf%u(:,j1:j2,:,7)
  lf%  e(:,j1:j2,:) = lf%u(:,j1:j2,:,8)
  lf%phi(:,j1:j2,:) = lf%u(:,j1:j2,:,9)

! pressure and total pressure
  vsq(:,j1:j2,:) = lf%v1(:,j1:j2,:)*lf%v1(:,j1:j2,:) &
                 + lf%v2(:,j1:j2,:)*lf%v2(:,j1:j2,:) &
                 + lf%v3(:,j1:j2,:)*lf%v3(:,j1:j2,:)
  bsq(:,j1:j2,:) = lf%b1(:,j1:j2,:)*lf%b1(:,j1:j2,:) &
                 + lf%b2(:,j1:j2,:)*lf%b2(:,j1:j2,:) &
                 + lf%b3(:,j1:j2,:)*lf%b3(:,j1:j2,:)
  lf%p(:,j1:j2,:) = (gamma-1.d0) &
                  * (lf%e(:,j1:j2,:) - 0.5d0*lf%d(:,j1:j2,:) * vsq(:,j1:j2,:) &
                                     - 0.5d0*                  bsq(:,j1:j2,:) )
  lf%ptot(:,j1:j2,:) = lf%p(:,j1:j2,:) + 5.d-1 * bsq(:,j1:j2,:)
 end if

! x3 direction ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
deallocate(nu)
allocate( nu(0:max(ib/cuts+1,2),0:max(jb/cuts+1,2),-1:1,1:9) )

elseif((n+1)/2==3)then
 i1 = ib/cuts * mod(lid+1,cuts)           ; i2 = i1 + max(ib/cuts+1,2)
 j1 = jb/cuts * mod(lid-1,cuts*cuts)/cuts ; j2 = j1 + max(jb/cuts+1,2)

 if(mod(n,2)==1)then ! inner boundary
  k1 = ks-2 ; k2 = ks-1
  if(sync)then
   if(rungen==1)then     ; nu = nlf%midu3(i1:i2,j1:j2,1:3,:)
   elseif(rungen==2)then ; nu = nlf%u(i1:i2,j1:j2,kb-1:kb+1,:)
   else
    nu = (nlf%u(i1:i2,j1:j2,kb-1:kb+1,:)+nlf%midu3(i1:i2,j1:j2,1:3,:))*5.d-1
   end if
  else
   if(rungen==1)then     ; nu = nlf%uorg(i1:i2,j1:j2,kb-1:kb+1,:)
   elseif(rungen==2)then ; nu = nlf%midu3(i1:i2,j1:j2,1:3,:)
   else
    nu = (nlf%uorg(i1:i2,j1:j2,kb-1:kb+1,:)+nlf%midu3(i1:i2,j1:j2,1:3,:))*5.d-1
   end if
  end if
  call amr_bound_intpol_x3(k1,nu,lf%u,lf%dvol)
  lf%  d(:,:,k1:k2) = lf%u(:,:,k1:k2,1)
  lf% v1(:,:,k1:k2) = lf%u(:,:,k1:k2,2) / lf%u(:,:,k1:k2,1)
  lf% v2(:,:,k1:k2) = lf%u(:,:,k1:k2,3) / lf%u(:,:,k1:k2,1)
  lf% v3(:,:,k1:k2) = lf%u(:,:,k1:k2,4) / lf%u(:,:,k1:k2,1)
  lf% b1(:,:,k1:k2) = lf%u(:,:,k1:k2,5)
  lf% b2(:,:,k1:k2) = lf%u(:,:,k1:k2,6)
  lf% b3(:,:,k1:k2) = lf%u(:,:,k1:k2,7)
  lf%  e(:,:,k1:k2) = lf%u(:,:,k1:k2,8)
  lf%phi(:,:,k1:k2) = lf%u(:,:,k1:k2,9)

! pressure and total pressure
  vsq(:,:,k1:k2) = lf%v1(:,:,k1:k2)*lf%v1(:,:,k1:k2) &
                 + lf%v2(:,:,k1:k2)*lf%v2(:,:,k1:k2) &
                 + lf%v3(:,:,k1:k2)*lf%v3(:,:,k1:k2)
  bsq(:,:,k1:k2) = lf%b1(:,:,k1:k2)*lf%b1(:,:,k1:k2) &
                 + lf%b2(:,:,k1:k2)*lf%b2(:,:,k1:k2) &
                 + lf%b3(:,:,k1:k2)*lf%b3(:,:,k1:k2)
  lf%p(:,:,k1:k2) = (gamma-1.d0) &
                  * (lf%e(:,:,k1:k2) - 0.5d0*lf%d(:,:,k1:k2) * vsq(:,:,k1:k2) &
                                     - 0.5d0*                  bsq(:,:,k1:k2) )
  lf%ptot(:,:,k1:k2) = lf%p(:,:,k1:k2) + 5.d-1 * bsq(:,:,k1:k2)

 else ! outer boundary
  k1 = kb+1 ; k2 = kb+2
  if(sync)then
   if(rungen==1)then     ; nu = nlf%midu3(i1:i2,j1:j2,-2:0,:)
   elseif(rungen==2)then ; nu = nlf%u(i1:i2,j1:j2,ks-1:ks+1,:)
   else
    nu = (nlf%u(i1:i2,j1:j2,ks-1:ks+1,:)+nlf%midu3(i1:i2,j1:j2,-2:0,:))*5.d-1
   end if
  else
   if(rungen==1)then     ; nu = nlf%uorg(i1:i2,j1:j2,ks-1:ks+1,:)
   elseif(rungen==2)then ; nu = nlf%midu3(i1:i2,j1:j2,-2:0,:)
   else
    nu = (nlf%uorg(i1:i2,j1:j2,ks-1:ks+1,:)+nlf%midu3(i1:i2,j1:j2,-2:0,:))*5.d-1
   end if
  end if
  call amr_bound_intpol_x3(k1,nu,lf%u,lf%dvol)
  lf%  d(:,:,k1:k2) = lf%u(:,:,k1:k2,1)
  lf% v1(:,:,k1:k2) = lf%u(:,:,k1:k2,2) / lf%u(:,:,k1:k2,1)
  lf% v2(:,:,k1:k2) = lf%u(:,:,k1:k2,3) / lf%u(:,:,k1:k2,1)
  lf% v3(:,:,k1:k2) = lf%u(:,:,k1:k2,4) / lf%u(:,:,k1:k2,1)
  lf% b1(:,:,k1:k2) = lf%u(:,:,k1:k2,5)
  lf% b2(:,:,k1:k2) = lf%u(:,:,k1:k2,6)
  lf% b3(:,:,k1:k2) = lf%u(:,:,k1:k2,7)
  lf%  e(:,:,k1:k2) = lf%u(:,:,k1:k2,8)
  lf%phi(:,:,k1:k2) = lf%u(:,:,k1:k2,9)

! pressure and total pressure
  vsq(:,:,k1:k2) = lf%v1(:,:,k1:k2)*lf%v1(:,:,k1:k2) &
                 + lf%v2(:,:,k1:k2)*lf%v2(:,:,k1:k2) &
                 + lf%v3(:,:,k1:k2)*lf%v3(:,:,k1:k2)
  bsq(:,:,k1:k2) = lf%b1(:,:,k1:k2)*lf%b1(:,:,k1:k2) &
                 + lf%b2(:,:,k1:k2)*lf%b2(:,:,k1:k2) &
                 + lf%b3(:,:,k1:k2)*lf%b3(:,:,k1:k2)
  lf%p(:,:,k1:k2) = (gamma-1.d0) &
                  * (lf%e(:,:,k1:k2) - 0.5d0*lf%d(:,:,k1:k2) * vsq(:,:,k1:k2) &
                                     - 0.5d0*                  bsq(:,:,k1:k2) )
  lf%ptot(:,:,k1:k2) = lf%p(:,:,k1:k2) + 5.d-1 * bsq(:,:,k1:k2)

 end if
else
 print *,"Error from amr_bound_coarse",n
 stop
end if


return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE AMR_BOUND_INTPOL_X1
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To interpolate boundary values for coarser neighbours in x direction.

subroutine amr_bound_intpol_x1(is,nu,u,dvol)

 use amr_module,only:jb,kb,cuts,cib

 implicit none

 integer j,k,id,il,jl,kl,n
 integer,intent(in):: is ! -1 -> inner, ib+1 -> outer
 real*8,allocatable,intent(in),dimension(:,:,:,:):: nu
 real*8,allocatable,intent(in),dimension(:,:,:):: dvol
 real*8,allocatable,intent(inout),dimension(:,:,:,:):: u
 real*8 :: dltp, dltm
 real*8,dimension(1:3):: dlt

!-----------------------------------------------------------------------------

do n = 1,9
 do k = 1, max(kb/cuts,1)
  do j = 1, max(jb/cuts,1)
   dltp = nu(1  ,j,k,n) - nu(0,j,k,n) ; dltm = nu(0,j,k,n) - nu( -1,j,k,n)
   dlt(1) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   dltp = nu(0,j+1,k,n) - nu(0,j,k,n) ; dltm = nu(0,j,k,n) - nu(0,j-1,k,n)
   dlt(2) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   dltp = nu(0,j,k+1,n) - nu(0,j,k,n) ; dltm = nu(0,j,k,n) - nu(0,j,k-1,n)
   dlt(3) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   do id = 1, cib
    il = is   +mod(id-1,2) 
    jl = 2*j-1+mod((id-1)/2,2)
    kl = 2*k-1+mod((id-1)/2**2,2)

    u(il,jl,kl,n) = nu(0,j,k,n) &
                + (-1)**il * 0.25d0 * dlt(1) * dvol(is   ,jl,kl)/dvol(il,jl,kl)&
                + (-1)**jl * 0.25d0 * dlt(2) * dvol(il,2*j-1,kl)/dvol(il,jl,kl)&
                + (-1)**kl * 0.25d0 * dlt(3) * dvol(il,jl,2*k-1)/dvol(il,jl,kl)
   end do

  end do
 end do
end do

return
end subroutine amr_bound_intpol_x1

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE AMR_BOUND_INTPOL_X2
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To interpolate boundary values for coarser neighbours in y direction.

subroutine amr_bound_intpol_x2(js,nu,u,dvol)

 use amr_module,only:ib,kb,cuts,cib

 implicit none

 integer i,k,id,il,jl,kl,n
 integer,intent(in):: js ! js-2 -> inner, jb+1 -> outer
 real*8,intent(in),dimension(1:ib/cuts,-1:1,1:kb/cuts,9):: nu
 real*8,allocatable,intent(in),dimension(:,:,:):: dvol
 real*8,allocatable,intent(out),dimension(:,:,:,:):: u
 real*8 :: dltp, dltm
 real*8,dimension(1:3):: dlt

!-----------------------------------------------------------------------------

do n = 1,9
 do k = 1, max(kb/2,1)
  do i = 1, max(ib/2,1)
   dltp = nu(i+1,0,k,n) - nu(i,0,k,n) ; dltm = nu(i,0,k,n) - nu(i-1,0,k,n)
   dlt(1) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   dltp = nu(i,0+1,k,n) - nu(i,0,k,n) ; dltm = nu(i,0,k,n) - nu(i,0-1,k,n)
   dlt(2) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))
 
   dltp = nu(i,0,k+1,n) - nu(i,0,k,n) ; dltm = nu(i,0,k,n) - nu(i,0,k-1,n)
   dlt(3) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   do id = 1, cib
    il = 1*i-1+mod(id-1,2) 
    jl = js   +mod((id-1)/2,2)
    kl = 2*k-1+mod((id-1)/2**2,2)
    u(il,jl,kl,n) = nu(i,0,k,n) &
                + (-1)**il * 0.25d0 * dlt(1) * dvol(2*i-1,jl,kl)/dvol(il,jl,kl)&
                + (-1)**jl * 0.25d0 * dlt(2) * dvol(il,js   ,kl)/dvol(il,jl,kl)&
                + (-1)**kl * 0.25d0 * dlt(3) * dvol(il,jl,2*k-1)/dvol(il,jl,kl)
   end do

  end do
 end do
end do

return
end subroutine amr_bound_intpol_x2

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE AMR_BOUND_INTPOL_X3
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To interpolate boundary values for coarser neighbours in y direction.

subroutine amr_bound_intpol_x3(ks,nu,u,dvol)

 use amr_module,only:ib,jb,cuts,cib

 implicit none

 integer i,j,id,il,jl,kl,n
 integer,intent(in):: ks ! js-2 -> inner, jb+1 -> outer
 real*8,intent(in),dimension(1:ib/cuts,1:jb/cuts,-1:1,9):: nu
 real*8,allocatable,intent(in),dimension(:,:,:):: dvol
 real*8,allocatable,intent(out),dimension(:,:,:,:):: u
 real*8 :: dltp, dltm
 real*8,dimension(1:3):: dlt

!-----------------------------------------------------------------------------

do n = 1,9
 do j = 1, max(jb/2,1)
  do i = 1, max(ib/2,1)
   dltp = nu(i+1,j,0,n) - nu(i,j,0,n) ; dltm = nu(i,j,0,n) - nu(i-1,j,0,n)
   dlt(1) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   dltp = nu(i,j+1,0,n) - nu(i,j,0,n) ; dltm = nu(i,j,0,n) - nu(i,j-1,0,n)
   dlt(2) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))
 
   dltp = nu(i,j,0+1,n) - nu(i,j,0,n) ; dltm = nu(i,j,0,n) - nu(i,j,0-1,n)
   dlt(3) = (sign(0.5d0,dltp)+sign(0.5d0,dltm))*min(abs(dltp),abs(dltm))

   do id = 1, cib
    il = 1*i-1+mod(id-1,2) 
    jl = 2*j-1+mod((id-1)/2,2)
    kl = ks   +mod((id-1)/2**2,2)
    u(il,jl,kl,n) = nu(i,j,0,n) &
                + (-1)**il * 0.25d0 * dlt(1) * dvol(2*i-1,jl,kl)/dvol(il,jl,kl)&
                + (-1)**jl * 0.25d0 * dlt(2) * dvol(il,2*j-1,kl)/dvol(il,jl,kl)&
                + (-1)**kl * 0.25d0 * dlt(3) * dvol(il,jl,ks   )/dvol(il,jl,kl)
   end do

  end do
 end do
end do

return
end subroutine amr_bound_intpol_x3

end subroutine amr_bound_coarse

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_BOUND_BOUND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for AMR blocks next to the domain boundary

subroutine amr_bound_bound(n,lf)

 use grid,only:is,js,ks
 use settings
 use amr_templates
 use amr_module,only:ib,jb,kb

 implicit none

 integer i,j,k
 integer,intent(in):: n
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

! Note:
!  0: Periodic b.c.
!  1: Reflective b.c.
!  2: Outgoing b.c.
!  3: Dirichlet b.c. (boundary values should be given elsewhere!)

! x1-direction ***************************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(n==1)then
! scalar values
  if(bc1is==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x1is"
   stop
  elseif(bc1is==1)then ! reflective
   do k = ks, kb
    do j = js, jb
     lf%d(is-2,j,k) = lf%d(is+1,j,k) ; lf%d(is-1,j,k) = lf%d(is,j,k)
     lf%ptot(is-2,j,k) = lf%ptot(is+1,j,k) ; lf%ptot(is-1,j,k) = lf%ptot(is,j,k)
     lf%e(is-2,j,k) = lf%e(is+1,j,k) ; lf%e(is-1,j,k) = lf%e(is,j,k)
     lf%b1(is-2,j,k) = lf%b1(is+1,j,k) ; lf%b1(is-1,j,k) = lf%b1(is,j,k)
     lf%b2(is-2,j,k) = lf%b2(is+1,j,k) ; lf%b2(is-1,j,k) = lf%b2(is,j,k)
     lf%b3(is-2,j,k) = lf%b3(is+1,j,k) ; lf%b3(is-1,j,k) = lf%b3(is,j,k)
    end do
   end do
  elseif(bc1is==2)then ! outgoing
   do k = ks, kb
    do j = js, jb
     lf%d(is-2,j,k) = lf%d(is,j,k) ; lf%d(is-1,j,k) = lf%d(is,j,k)
     lf%ptot(is-2,j,k) = lf%ptot(is,j,k) ; lf%ptot(is-1,j,k) = lf%ptot(is,j,k)
     lf%e(is-2,j,k) = lf%e(is,j,k) ; lf%e(is-1,j,k) = lf%e(is,j,k)
     lf%b1(is-2,j,k) = lf%b1(is,j,k) ; lf%b1(is-1,j,k) = lf%b1(is,j,k)
     lf%b2(is-2,j,k) = lf%b2(is,j,k) ; lf%b2(is-1,j,k) = lf%b2(is,j,k)
     lf%b3(is-2,j,k) = lf%b3(is,j,k) ; lf%b3(is-1,j,k) = lf%b3(is,j,k)
    end do
   end do
  elseif(bc1is==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x1is"
   stop
  else
   print *,"Error from AMR x1 scalar inner boundary condition",bc1is
   stop
  end if

! vector values
  if(bc1iv==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x1iv"
   stop
  elseif(bc1iv==1)then ! reflective
   do k = ks, kb
    do j = js, jb
     lf%v1(is-1,j,k) = -lf%v1(is,j,k) ; lf%v1(is-2,j,k) = -lf%v1(is+1,j,k)
     lf%v2(is-1,j,k) =  lf%v2(is,j,k) ; lf%v2(is-2,j,k) =  lf%v2(is+1,j,k)
     lf%v3(is-1,j,k) =  lf%v3(is,j,k) ; lf%v3(is-2,j,k) =  lf%v3(is+1,j,k)
    end do
   end do
  elseif(bc1iv==2)then ! outgoing
   do k = ks, kb
    do j = js, jb
     lf%v1(is-1,j,k) = lf%v1(is,j,k) ; lf%v1(is-2,j,k) = lf%v1(is,j,k)
     lf%v2(is-1,j,k) = lf%v2(is,j,k) ; lf%v2(is-2,j,k) = lf%v2(is,j,k)
     lf%v3(is-1,j,k) = lf%v3(is,j,k) ; lf%v3(is-2,j,k) = lf%v3(is,j,k)
    end do
   end do
  elseif(bc1iv==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x1iv"
   stop
  else
   print *, "Error from AMR x1 vector inner boundary condition",bc1iv
   stop
  end if

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(n==2)then
! scalar values
  if(bc1os==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x1os"
   stop
  elseif(bc1os==1)then ! reflective
   do k = ks, kb
    do j = js, jb
     lf%d(ib+1,j,k) = lf%d(ib,j,k) ; lf%d(ib+2,j,k) = lf%d(ib-1,j,k)
     lf%ptot(ib+1,j,k) = lf%ptot(ib,j,k) ; lf%ptot(ib+2,j,k) = lf%ptot(ib-1,j,k)
     lf%e(ib+1,j,k) = lf%e(ib,j,k) ; lf%e(ib+2,j,k) = lf%e(ib-1,j,k)
     lf%b1(ib+1,j,k) = lf%b1(ib,j,k) ; lf%b1(ib+2,j,k) = lf%b1(ib-1,j,k)
     lf%b2(ib+1,j,k) = lf%b2(ib,j,k) ; lf%b2(ib+2,j,k) = lf%b2(ib-1,j,k)
     lf%b3(ib+1,j,k) = lf%b3(ib,j,k) ; lf%b3(ib+2,j,k) = lf%b3(ib-1,j,k)
    end do
   end do
  elseif(bc1os==2)then ! outgoing
   do k = ks, kb
    do j = js, jb
     lf%d(ib+1,j,k) = lf%d(ib,j,k) ; lf%d(ib+2,j,k) = lf%d(ib,j,k)
     lf%ptot(ib+1,j,k) = lf%ptot(ib,j,k) ; lf%ptot(ib+2,j,k) = lf%ptot(ib,j,k)
     lf%e(ib+1,j,k) = lf%e(ib,j,k) ; lf%e(ib+2,j,k) = lf%e(ib,j,k)
     lf%b1(ib+1,j,k) = lf%b1(ib,j,k) ; lf%b1(ib+2,j,k) = lf%b1(ib,j,k)
     lf%b2(ib+1,j,k) = lf%b2(ib,j,k) ; lf%b2(ib+2,j,k) = lf%b2(ib,j,k)
     lf%b3(ib+1,j,k) = lf%b3(ib,j,k) ; lf%b3(ib+2,j,k) = lf%b3(ib,j,k)
    end do
   end do
  elseif(bc1os==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x1os"
   stop
  else
   print *, "Error from AMR x1 scalar outer boundary condition",bc1os
   stop
  end if

! vector values
  if(bc1ov==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x1ov"
   stop
  elseif(bc1ov==1)then ! reflective
   do k = ks, kb
    do j = js, jb
     lf%v1(ib+1,j,k) = -lf%v1(ib,j,k) ; lf%v1(ib+2,j,k) = -lf%v1(ib-1,j,k)
     lf%v2(ib+1,j,k) =  lf%v2(ib,j,k) ; lf%v2(ib+2,j,k) =  lf%v2(ib-1,j,k)
     lf%v3(ib+1,j,k) =  lf%v3(ib,j,k) ; lf%v3(ib+2,j,k) =  lf%v3(ib-1,j,k)
    end do
   end do
  elseif(bc1ov==2)then ! outgoing
   do k = ks, kb
    do j = js, jb
     lf%v1(ib+1,j,k) = lf%v1(ib,j,k) ; lf%v1(ib+2,j,k) = lf%v1(ib,j,k)
     lf%v2(ib+1,j,k) = lf%v2(ib,j,k) ; lf%v2(ib+2,j,k) = lf%v2(ib,j,k)
     lf%v3(ib+1,j,k) = lf%v3(ib,j,k) ; lf%v3(ib+2,j,k) = lf%v3(ib,j,k)
    end do
   end do
  elseif(bc1ov==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x1ov"
   stop
  else
   print *, "Error from AMR x1 vector outer boundary condition",bc1ov
   stop
  end if

! x2-direction ***************************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(n==3)then
! scalar values
  if(bc2is==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x2is"
   stop
  elseif(bc2is==1)then ! reflective
   do k = ks, kb
    do i = is, ib
     lf%d(i,js-1,k) = lf%d(i,js,k) ; lf%d(i,js-2,k) = lf%d(i,js+1,k)
     lf%ptot(i,js-1,k) = lf%ptot(i,js,k) ; lf%ptot(i,js-2,k) = lf%ptot(i,js+1,k)
     lf%e(i,js-1,k) = lf%e(i,js,k) ; lf%e(i,js-2,k) = lf%e(i,js+1,k)
     lf%b1(i,js-1,k) = lf%b1(i,js,k) ; lf%b1(i,js-2,k) = lf%b1(i,js+1,k)
     lf%b2(i,js-1,k) = lf%b2(i,js,k) ; lf%b2(i,js-2,k) = lf%b2(i,js+1,k)
     lf%b3(i,js-1,k) = lf%b3(i,js,k) ; lf%b3(i,js-2,k) = lf%b3(i,js+1,k)
    end do
   end do
  elseif(bc2is==2)then ! outgoing
   do k = ks, kb
    do i = is, ib
     lf%d(i,js-1,k) = lf%d(i,js,k) ; lf%d(i,js-2,k) = lf%d(i,js,k)
     lf%ptot(i,js-1,k) = lf%ptot(i,js,k) ; lf%ptot(i,js-2,k) = lf%ptot(i,js,k)
     lf%e(i,js-1,k) = lf%e(i,js,k) ; lf%e(i,js-2,k) = lf%e(i,js,k)
     lf%b1(i,js-1,k) = lf%b1(i,js,k) ; lf%b1(i,js-2,k) = lf%b1(i,js,k)
     lf%b2(i,js-1,k) = lf%b2(i,js,k) ; lf%b2(i,js-2,k) = lf%b2(i,js,k)
     lf%b3(i,js-1,k) = lf%b3(i,js,k) ; lf%b3(i,js-2,k) = lf%b3(i,js,k)  
    end do
   end do
  elseif(bc2is==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x2is"
   stop
  else
   print *, "Error from AMR x2 scalar inner boundary condition",bc2is
   stop
  end if

! vector values
  if(bc2iv==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x2iv"
   stop
  elseif(bc2iv==1)then ! reflective
   do k = ks, kb
    do i = is, ib
     lf%v1(i,js-1,k) =  lf%v1(i,js,k) ; lf%v1(i,js-2,k) =  lf%v1(i,js+1,k)
     lf%v2(i,js-1,k) = -lf%v2(i,js,k) ; lf%v2(i,js-2,k) = -lf%v2(i,js+1,k)
     lf%v3(i,js-1,k) =  lf%v3(i,js,k) ; lf%v3(i,js-2,k) =  lf%v3(i,js+1,k)
    end do
   end do
  elseif(bc2iv==2)then ! outgoing
   do k = ks, kb
    do i = is, ib
     lf%v1(i,js-1,k) = lf%v1(i,js,k) ; lf%v1(i,js-2,k) = lf%v1(i,js,k)
     lf%v2(i,js-1,k) = lf%v2(i,js,k) ; lf%v2(i,js-2,k) = lf%v2(i,js,k)
     lf%v3(i,js-1,k) = lf%v3(i,js,k) ; lf%v3(i,js-2,k) = lf%v3(i,js,k)
    end do
   end do
  elseif(bc2iv==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x2iv"
   stop
  else
   print *, "Error from AMR x2 vector inner boundary condition",bc2iv
   stop
  end if

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(n==4)then
!scalar values
  if(bc2os==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x2os"
   stop
  elseif(bc2os==1)then ! reflective
   do k = ks, kb
    do i = is, ib
     lf%d(i,jb+1,k) = lf%d(i,jb,k) ; lf%d(i,jb+2,k) = lf%d(i,jb-1,k)
     lf%ptot(i,jb+1,k) = lf%ptot(i,jb,k) ; lf%ptot(i,jb+2,k) = lf%ptot(i,jb-1,k)
     lf%e(i,jb+1,k) = lf%e(i,jb,k) ; lf%e(i,jb+2,k) = lf%e(i,jb-1,k)
     lf%b1(i,jb+1,k) = lf%b1(i,jb,k) ; lf%b1(i,jb+2,k) = lf%b1(i,jb-1,k)
     lf%b2(i,jb+1,k) = lf%b2(i,jb,k) ; lf%b2(i,jb+2,k) = lf%b2(i,jb-1,k)
     lf%b3(i,jb+1,k) = lf%b3(i,jb,k) ; lf%b3(i,jb+2,k) = lf%b3(i,jb-1,k)
    end do
   end do
  elseif(bc2os==2)then ! outgoing
   do k = ks, kb
    do i = is, ib
     lf%d(i,jb+1,k) = lf%d(i,jb,k) ; lf%d(i,jb+2,k) = lf%d(i,jb,k)
     lf%ptot(i,jb+1,k) = lf%ptot(i,jb,k) ; lf%ptot(i,jb+2,k) = lf%ptot(i,jb,k)
     lf%e(i,jb+1,k) = lf%e(i,jb,k) ; lf%e(i,jb+2,k) = lf%e(i,jb,k)
     lf%b1(i,jb+1,k) = lf%b1(i,jb,k) ; lf%b1(i,jb+2,k) = lf%b1(i,jb,k)
     lf%b2(i,jb+1,k) = lf%b2(i,jb,k) ; lf%b2(i,jb+2,k) = lf%b2(i,jb,k)
     lf%b3(i,jb+1,k) = lf%b3(i,jb,k) ; lf%b3(i,jb+2,k) = lf%b3(i,jb,k)  
    end do
   end do
  elseif(bc2os==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x2os"
   stop
  else
   print *, "Error from AMR x2 scalar outer boundary condition",bc2os
   stop
  end if

! vector values
  if(bc2ov==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x2ov"
   stop
  elseif(bc2ov==1)then ! reflective
   do k = ks, kb
    do i = is, ib
     lf%v1(i,jb+1,k) =  lf%v1(i,jb,k) ; lf%v1(i,jb+2,k) =  lf%v1(i,jb-1,k)
     lf%v2(i,jb+1,k) = -lf%v2(i,jb,k) ; lf%v2(i,jb+2,k) = -lf%v2(i,jb-1,k)
     lf%v3(i,jb+1,k) =  lf%v3(i,jb,k) ; lf%v3(i,jb+2,k) =  lf%v3(i,jb-1,k)
    end do
   end do
  elseif(bc2ov==2)then ! outgoing
   do k = ks, kb
    do i = is, ib
     lf%v1(i,jb+1,k) = lf%v1(i,jb,k) ; lf%v1(i,jb+2,k) = lf%v1(i,jb,k)
     lf%v2(i,jb+1,k) = lf%v2(i,jb,k) ; lf%v2(i,jb+2,k) = lf%v2(i,jb,k)
     lf%v3(i,jb+1,k) = lf%v3(i,jb,k) ; lf%v3(i,jb+2,k) = lf%v3(i,jb,k)
    end do
   end do
  elseif(bc2ov==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x2ov"
   stop
  else
   print *, "Error from AMR x2 vector outer boundary condition",bc2ov
   stop
  end if

! x3-direction **************************************************************

! >>> inner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(n==5)then
! scalar values
  if(bc3is==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x3is"
   stop
  elseif(bc3is==1)then ! reflective
   do j = js, jb
    do i = is, ib
     lf%d(i,j,ks-1) = lf%d(i,j,ks) ; lf%d(i,j,ks-2) = lf%d(i,j,ks+1)
     lf%ptot(i,j,ks-1) = lf%ptot(i,j,ks) ; lf%ptot(i,j,ks-2) = lf%ptot(i,j,ks+1)
     lf%e(i,j,ks-1) = lf%e(i,j,ks) ; lf%e(i,j,ks-2) = lf%e(i,j,ks+1)
     lf%b1(i,j,ks-1) = lf%b1(i,j,ks) ; lf%b1(i,j,ks-2) = lf%b1(i,j,ks+1)
     lf%b2(i,j,ks-1) = lf%b2(i,j,ks) ; lf%b2(i,j,ks-2) = lf%b2(i,j,ks+1)
     lf%b3(i,j,ks-1) = lf%b3(i,j,ks) ; lf%b3(i,j,ks-2) = lf%b3(i,j,ks+1)
    end do
   end do
  elseif(bc3is==2)then ! outgoing
   do j = js, jb
    do i = is, ib
     lf%d(i,j,ks-1) = lf%d(i,j,ks) ; lf%d(i,j,ks-2) = lf%d(i,j,ks)
     lf%ptot(i,j,ks-1) = lf%ptot(i,j,ks) ; lf%ptot(i,j,ks-2) = lf%ptot(i,j,ks)
     lf%e(i,j,ks-1) = lf%e(i,j,ks) ; lf%e(i,j,ks-2) = lf%e(i,j,ks)
     lf%b1(i,j,ks-1) = lf%b1(i,j,ks) ; lf%b1(i,j,ks-2) = lf%b1(i,j,ks)
     lf%b2(i,j,ks-1) = lf%b2(i,j,ks) ; lf%b2(i,j,ks-2) = lf%b2(i,j,ks)
     lf%b3(i,j,ks-1) = lf%b3(i,j,ks) ; lf%b3(i,j,ks-2) = lf%b3(i,j,ks)  
    end do
   end do
  elseif(bc3is==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x3is"
   stop
  else
   print *, "Error from AMR x3 scalar inner boundary condition",bc3is
   stop
  end if

! vector values
  if(bc3iv==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x3iv"
   stop
  elseif(bc3iv==1)then ! reflective
   do j = js, jb
    do i = is, ib
     lf%v1(i,j,ks-1) =  lf%v1(i,j,ks) ; lf%v1(i,j,ks-2) =  lf%v1(i,j,ks+1)
     lf%v2(i,j,ks-1) =  lf%v2(i,j,ks) ; lf%v2(i,j,ks-2) =  lf%v2(i,j,ks+1)
     lf%v3(i,j,ks-1) = -lf%v3(i,j,ks) ; lf%v3(i,j,ks-2) = -lf%v3(i,j,ks+1)
    end do
   end do
  elseif(bc3iv==2)then ! outgoing
   do j = js, jb
    do i = is, ib
     lf%v1(i,j,ks-1) = lf%v1(i,j,ks) ; lf%v1(i,j,ks-2) = lf%v1(i,j,ks)
     lf%v2(i,j,ks-1) = lf%v2(i,j,ks) ; lf%v2(i,j,ks-2) = lf%v2(i,j,ks)
     lf%v3(i,j,ks-1) = lf%v3(i,j,ks) ; lf%v3(i,j,ks-2) = lf%v3(i,j,ks)
    end do
   end do
  elseif(bc3iv==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x3iv"
   stop
  else
   print *, "Error from AMR x3 vector inner boundary condition",bc3iv
   stop
  end if

! >>> outer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 elseif(n==6)then
! scalar values
  if(bc3os==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x3os"
   stop
  elseif(bc3os==1)then ! reflective
   do j = js, jb
    do i = is, ib
     lf%d(i,j,kb+1) = lf%d(i,j,kb) ; lf%d(i,j,kb+2) = lf%d(i,j,kb-1)
     lf%ptot(i,j,kb+1) = lf%ptot(i,j,kb) ; lf%ptot(i,j,kb+2) = lf%ptot(i,j,kb-1)
     lf%e(i,j,kb+1) = lf%e(i,j,kb) ; lf%e(i,j,kb+2) = lf%e(i,j,kb-1)
     lf%b1(i,j,kb+1) = lf%b1(i,j,kb) ; lf%b1(i,j,kb+2) = lf%b1(i,j,kb-1)
     lf%b2(i,j,kb+1) = lf%b2(i,j,kb) ; lf%b2(i,j,kb+2) = lf%b2(i,j,kb-1)
     lf%b3(i,j,kb+1) = lf%b3(i,j,kb) ; lf%b3(i,j,kb+2) = lf%b3(i,j,kb-1)
    end do
   end do
  elseif(bc3os==2)then ! outgoing
   do j = js, jb
    do i = is, ib
     lf%d(i,j,kb+1) = lf%d(i,j,kb) ; lf%d(i,j,kb+2) = lf%d(i,j,kb)
     lf%ptot(i,j,kb+1) = lf%ptot(i,j,kb) ; lf%ptot(i,j,kb+2) = lf%ptot(i,j,kb)
     lf%e(i,j,kb+1) = lf%e(i,j,kb) ; lf%e(i,j,kb+2) = lf%e(i,j,kb)
     lf%b1(i,j,kb+1) = lf%b1(i,j,kb) ; lf%b1(i,j,kb+2) = lf%b1(i,j,kb)
     lf%b2(i,j,kb+1) = lf%b2(i,j,kb) ; lf%b2(i,j,kb+2) = lf%b2(i,j,kb)
     lf%b3(i,j,kb+1) = lf%b3(i,j,kb) ; lf%b3(i,j,kb+2) = lf%b3(i,j,kb)  
    end do
   end do
  elseif(bc3os==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x3os"
   stop
  else
   print *, "Error from AMR x3 scalar outer boundary condition",bc3os
   stop
  end if

! vector values
  if(bc3ov==0)then ! periodic
   print *,"Periodic boundary for AMR is under construction : x3ov"
   stop
  elseif(bc3ov==1)then ! reflective
   do j = js, jb
    do i = is, ib
     lf%v1(i,j,kb+1) =  lf%v1(i,j,kb) ; lf%v1(i,j,kb+2) =  lf%v1(i,j,kb-1)
     lf%v2(i,j,kb+1) =  lf%v2(i,j,kb) ; lf%v2(i,j,kb+2) =  lf%v2(i,j,kb-1)
     lf%v3(i,j,kb+1) = -lf%v3(i,j,kb) ; lf%v3(i,j,kb+2) = -lf%v3(i,j,kb-1)
    end do
   end do
  elseif(bc3ov==2)then ! outgoing
   do j = js, jb
    do i = is, ib
     lf%v1(i,j,kb+1) = lf%v1(i,j,kb) ; lf%v1(i,j,kb+2) = lf%v1(i,j,kb)
     lf%v2(i,j,kb+1) = lf%v2(i,j,kb) ; lf%v2(i,j,kb+2) = lf%v2(i,j,kb)
     lf%v3(i,j,kb+1) = lf%v3(i,j,kb) ; lf%v3(i,j,kb+2) = lf%v3(i,j,kb)
    end do
   end do
  elseif(bc3ov==3)then ! Dirichlet
   print *,"Dirichlet boundary for AMR is under construction : x3ov"
   stop
  else
   print *, "Error from AMR x3 vector outer boundary condition",bc3ov
   stop
  end if

 else
  print *,"Error in number of faces",n
  stop
 end if

return
end subroutine amr_bound_bound

end module amr_bound_mod


module amr_boundary_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE AMR_BOUNDARY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for AMR blocks.
! CAUTION: needs bkn

subroutine amr_boundary

 use amr_module
 use amr_bound_mod

 implicit none

 integer n

!-----------------------------------------------------------------------------

 do n = 1, face
  if    (bk(bkn)%ldif(n)>= 0)then ! If neighbour is same level
   call amr_bound_same(n,lf(bkn)% d,lf(bk(bkn)%neigh(n))% d)
   call amr_bound_same(n,lf(bkn)%v1,lf(bk(bkn)%neigh(n))%v1)
   call amr_bound_same(n,lf(bkn)%v2,lf(bk(bkn)%neigh(n))%v2)
   call amr_bound_same(n,lf(bkn)%v3,lf(bk(bkn)%neigh(n))%v3)
   call amr_bound_same(n,lf(bkn)%b1,lf(bk(bkn)%neigh(n))%b1)
   call amr_bound_same(n,lf(bkn)%b2,lf(bk(bkn)%neigh(n))%b2)
   call amr_bound_same(n,lf(bkn)%b3,lf(bk(bkn)%neigh(n))%b3)
   call amr_bound_same(n,lf(bkn)% e,lf(bk(bkn)%neigh(n))% e)
   call amr_bound_same(n,lf(bkn)% p,lf(bk(bkn)%neigh(n))% p)
   call amr_bound_same(n,lf(bkn)%ptot,lf(bk(bkn)%neigh(n))%ptot)
   call amr_bound_same(n,lf(bkn)% phi,lf(bk(bkn)%neigh(n))% phi)

  elseif(bk(bkn)%ldif(n)==-1)then ! If neighbour is coarser
   do lid = 1, cib
    if(bk(bk(bkn)%parnt)%child(lid)==bkn)then
     call amr_bound_coarse(n,lid,lf(bkn),lf(bk(bkn)%neigh(n)))
     exit
    end if
   end do
 
  elseif(bk(bkn)%ldif(n)==-2)then ! If neighbour is boundary
   call amr_bound_bound(n,lf(bkn))

  end if
 end do

return
end subroutine amr_boundary

end module amr_boundary_mod
