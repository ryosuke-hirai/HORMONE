!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE AMR_CRITERION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set criterion for refinement in AMR mode.

subroutine amr_criterion(lf,bk)

  use grid,only:is,js,ks
  use amr_templates
  use amr_module,only:ib,jb,kb,maxamr

  implicit none

  integer i,j,k, i1,i2,j1,j2,k1,k2
  type(amr_header),intent(inout)::bk
  type(leaf_contents),intent(in)::lf
  real*8 error
  real*8,parameter:: epsl = 0.01d0

!-----------------------------------------------------------------------------

! ref_deref =  1 : refine
! ref_deref = -1 : derefine

  i1 = is-1 ; i2 = ib+1 ; j1 = js-1 ; j2 = jb+1 ; k1 = ks-1 ; k2 = kb+1
  if(jb==1)then ; j1 = 1 ; j2 = 1 ; endif
  if(kb==1)then ; k1 = 1 ; k2 = 1 ; endif

  bk%ref_deref = 0
  error = 0.d0

  do k = k1, k2
   do j = j1, j2
    do i = i1, i2
     error = max(abs(lf%d(i+1,j,k)-2.d0*lf%d(i,j,k)+lf%d(i-1,j,k)) &
           / ( abs(lf%d(i+1,j,k)-lf%d(i,j,k)) + abs(lf%d(i,j,k)-lf%d(i-1,j,k)) &
           + epsl*(abs(lf%d(i+1,j,k))+2.d0*abs(lf%d(i,j,k))+abs(lf%d(i-1,j,k)))),error)
error = max(error,max(abs(lf%e(i+1,j,k)-2.d0*lf%e(i,j,k)+lf%e(i-1,j,k)) &
           / ( abs(lf%e(i+1,j,k)-lf%e(i,j,k)) + abs(lf%e(i,j,k)-lf%e(i-1,j,k)) &
           + epsl*(abs(lf%e(i+1,j,k))+2.d0*abs(lf%e(i,j,k))+abs(lf%e(i-1,j,k)))),error))
    end do
   end do
  end do

if(error>0.8d0)then
 bk%ref_deref = 1
elseif(error<0.5d0.and.bk%level>0.and.maxval(bk%ldif)<=0)then
 bk%ref_deref = -1
end if

return
end subroutine amr_criterion
