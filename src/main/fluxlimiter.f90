module fluxlimiter_mod

 implicit none

 public:: fluxlimiter
 private:: minmod,modified_mc

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE FLUXLIMITER
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To apply flux limiter for Godunov flux

subroutine fluxlimiter(uu,dx,x,xi,du)

 use settings,only:flux_limiter

 real(8),intent(in):: uu(1:3),dx(1:2),x(1:3),xi(1:2)
 real(8),intent(out):: du

!-----------------------------------------------------------------------------

 select case(flux_limiter)
 case('minmod')
  call minmod(du,uu,dx)
 case('modified_mc')
  call modified_mc(du,uu,x,xi)
 case default
  print*, "Error in flux_limiter",flux_limiter
  stop
 end select

return
end subroutine fluxlimiter

subroutine minmod(mm,u,dx)

  real(8),intent(in ):: u(1:3),dx(1:2)
  real(8),intent(out):: mm
  real(8):: x,y,z

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Generalized minmod
  x = 1.4d0*(u(3) - u(2))*dx(2)
  y = 1.4d0*(u(2) - u(1))*dx(1)
  z = (u(3)-u(1))*dx(1)*dx(2)/sum(dx)

  mm = sign(1d0,x) * max(0.d0,min(abs(x),sign(1d0,x)*y,sign(1d0,x)*z))

return
end subroutine minmod

subroutine modified_mc(mm,u,x,xi)

 real(8),intent(in ):: u(1:3),x(1:3),xi(1:2)
 real(8),intent(out):: mm
 real(8):: cf,cb,top,bot,sig

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  cf = (x(3)-x(2))/(xi(2)-x(2))
  cb = (x(2)-x(1))/(x(2)-xi(1))
  top = (u(2)-u(1))*(x(3)-x(2))/(x(2)-x(1))
  bot = u(3)-u(2)
  sig = sign(1d0,bot)

  mm =  sig*max(0d0,min(sig*0.5d0*(bot+top),sig*cf*bot,sig*cb*top))/(x(3)-x(2))

return
end subroutine modified_mc

end module fluxlimiter_mod
