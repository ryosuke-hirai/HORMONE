module amr_timestep_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE AMR_TIMESTEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set timestep for AMR mode

subroutine amr_timestep

 use settings,only:courant
 use grid,only:dt,is,js,ks
 use ninewave
 use amr_module
 use amr_eigen_mod

 implicit none

 integer i,j,k
 real*8 lf_dt, cfmax

!-----------------------------------------------------------------------------

 dt = 1.d50
 cfmax = 0.d0

 do lfn = 1, totbloks
  call amr_eigen(lf(lfn))
  do k = ks, kb
   do j = js, jb
    do i = is, ib
     lf_dt = min( &
             lf(lfn)%dx1(i) &
           / ( abs(lf(lfn)%v1(i,j,k))+lf(lfn)%cf(i,j,k) ), &
             lf(lfn)%g22(i)   * lf(lfn)%dx2(j) &
           / ( abs(lf(lfn)%v2(i,j,k))+lf(lfn)%cf(i,j,k) ), &
             lf(lfn)%g33(i,j) * lf(lfn)%dx3(k) &
           / ( abs(lf(lfn)%v3(i,j,k))+lf(lfn)%cf(i,j,k) ) )
! compare it with dt's of all AMR levels
     dt = min( lf_dt * dble(2**bk(lfn)%level) , dt )

     cfmax = max(abs(lf(lfn)%cf(i,j,k)),cfmax)
    end do
   end do
  end do
 end do

 dt = courant * dt ! timestep for lowest level

 ch = cfmax

return
end subroutine amr_timestep

end module amr_timestep_mod
