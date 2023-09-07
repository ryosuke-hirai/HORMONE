module composition_mod
 implicit none

 public:: meanmolweight, get_imu
 
 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE MEANMOLWEIGHT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Evolve chemical composition distribution.

subroutine meanmolweight

 use settings
 use grid
 use physval
 use gravmod,only:mc
 use utils,only:masscoordinate

 integer:: i,j,k

!-----------------------------------------------------------------------------

 composition_type: select case (compswitch)
 case(0) composition_type ! for uniform mu
  return

 case(1) composition_type ! for fixed mu distribution

  call masscoordinate

  do i = is, ie
   if(mc(i)<mudata(1,0))then
    imu(i,js,ks) = 1d0/mudata(1,1)
   elseif(mc(i)<mudata(musize,0))then
    do j = 1, musize-1
     if(mudata(j+1,0)>mc(i).and.mudata(j,0)<=mc(i))then
      imu(i,js,ks) = ( (mudata(j+1,0)-mc(i)) * mudata(j,1)     &
                     + (mc(i)-mudata(  j,0)) * mudata(j+1,1) ) &
                   / (mudata(j+1,0)-mudata(j,0))
      imu(i,js,ks) = 1d0/imu(i,js,ks)
     end if
    end do
   else
    imu(i,js,ks) = 1d0/muconst
   end if
  end do
  imu(is-2:is-1,js,ks) = imu(is,js,ks)
  imu(ie+1:ie+2,js,ks) = imu(ie,js,ks)

 case(2) composition_type ! for composition advection
!$omp parallel
!$omp do private(i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     imu(i,j,k) = get_imu(spc(1:2,i,j,k))!0.25d0*(6d0*spc(1,i,j,k)+spc(2,i,j,k)+2d0)
    end do
   end do
  end do
!$omp end do
!$omp do private(i,j,k)
  do k = ks, ke
   do j = js-2, js-1
    do i = is, ie
     imu(i,j,k) = imu(i,js,k)
    end do
   end do
   do j = js, je
    do i = is-2, is-1
     imu(i,j,k) = imu(is,j,k)
    end do
    do i = ie+1, ie+2
     imu(i,j,k) = imu(ie,j,k)
    end do
   end do
   do j = je+1, je+2
    do i = is, ie
     imu(i,j,k) = imu(i,je,k)
    end do
   end do
  end do
!$omp end do nowait
!$omp do private(i,j) collapse(2)
  do j = js, je
   do i = is, ie
    imu(i,j,ks-2:ks-1) = imu(i,j,ks)
    imu(i,j,ke+1:ke+2) = imu(i,j,ke)
   end do
  end do
!$omp end do
!$omp end parallel
 case default composition_type
  print *, 'Error in compswitch',compswitch
  stop
 end select composition_type

return
end subroutine meanmolweight

function get_imu(spc) result(imu)
 implicit none
 real(8),intent(in)::spc(1:2)
 real(8):: imu

 imu = 0.25d0*(6d0*spc(1)+spc(2)+2d0)
 
end function get_imu

end module composition_mod
