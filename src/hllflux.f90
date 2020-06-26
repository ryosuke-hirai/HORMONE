!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE HLLFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte hllflux

subroutine hllflux(fflux,fl,fr,ul,ur,cfl,cfr,vl,vr)

  use grid
  use physval

  implicit none

  real*8,intent(in) :: fl, fr, ul, ur, cfl, cfr, vl, vr
  real*8 sl, sr
  real*8,intent(inout)::fflux

!--------------------------------------------------------------------

! find left and right waves
  sr  = max(vr+cfr,vl+cfl,0.d0)
  sl  = min(vr-cfr,vl-cfl,0.d0)

  fflux = (sr*fl - sl*fr + sr*sl*(ur-ul)) / (sr-sl)

return
end subroutine hllflux
