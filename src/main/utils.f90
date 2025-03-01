module utils
 implicit none

contains

! In-house sine function to avoid machine-dependent behaviour
 function sin0(x) result(s)
  use constants,only:pi
  real(8), intent(in) :: x
  real(8):: s, x_reduced, cosx
  integer:: sign_factor

! Reduce x to [0, 2π] using periodicity
  x_reduced = modulo(x,2d0*pi)

! Reduce to [0, π] and keep track of sign
  sign_factor = 1
  if (x_reduced > pi) then
   x_reduced = 2d0*pi - x_reduced
   sign_factor = -1
  end if

! Reduce to [0, π/2] using symmetry sin(π - x) = sin(x)
  if (x_reduced > 0.5d0*pi) then
   x_reduced = pi - x_reduced
  end if

! Reduce to [0, π/4] using sin(x) = cos(π/2 - x) and cos(x) = sin(π/2 - x)
  if (x_reduced > 0.25d0*pi) then
   cosx = sine_taylor(0.5d0*pi-x_reduced)
   s = dble(sign_factor) * sqrt(1d0-cosx*cosx)
  else
   s = dble(sign_factor) * sine_taylor(x_reduced)
  end if
 end function sin0

 function cos0(x) result(c)
  use constants,only: pi
  real(8),intent(in):: x
  real(8):: c
  c = sin0(0.5d0*pi-x)
 end function cos0

 pure function sine_taylor(x) result(s)
  real(8), intent(in) :: x
  real(8) :: s, term
  real(8),parameter:: tolerance=1d-15
  integer :: n

! Taylor series initialization
  s = 0d0
  term = x
  n = 1

! Compute sine via Taylor series: sin(x) = x - x^3/3! + x^5/5! - ...
  do while (abs(term) > tolerance)
   s = s + term
   n = n + 2
   term = -term*x*x / dble(n*(n-1))
  end do
 end function sine_taylor
 
! convert polar to cartesian coordinates
 pure function polcar(xp) result(x)
  implicit none
  real(8),intent(in):: xp(1:3)
  real(8):: x(1:3)
  real(16) :: xp1,xp2,xp3

  ! Convert to quad precision
  xp1 = real(xp(1),kind=16)
  xp2 = real(xp(2),kind=16)
  xp3 = real(xp(3),kind=16)

  x(1) = xp(1)*real(sin(xp2)*cos(xp3),kind=8)
  x(2) = xp(1)*real(sin(xp2)*sin(xp3),kind=8)
  x(3) = xp(1)*real(cos(xp2),kind=8)

 end function polcar

! convert cartesian to polar coordinates
 pure function carpol(x) result(xp)
  implicit none
  real(8),intent(in):: x(1:3)
  real(8):: xp(1:3), arg
  xp(1) = norm2(x)
  arg = x(3)/max(xp(1),tiny(1d0)) ! avoid dividing by 0
  arg = sign(min(abs(arg),1d0),arg) ! avoid arguments for acos being >1 or <-1
  xp(2) = acos(arg)
  xp(3) = atan2(x(2),x(1))
 end function carpol

! convert cylindrical to cartesian coordinates
 pure function cylcar(xp) result(x)
  implicit none
  real(8),intent(in):: xp(1:3)
  real(8):: x(1:3)
  real(16) :: xp2

  ! Convert to quad precision
  xp2 = real(xp(2),kind=16)

  x(1) = xp(1)*real(cos(xp2),kind=8)
  x(2) = xp(1)*real(sin(xp2),kind=8)
  x(3) = xp(3)

 end function cylcar

! convert cartesian to cylindrical coordinates
 pure function carcyl(x) result(xp)
  implicit none
  real(8),intent(in):: x(1:3)
  real(8):: xp(1:3)
  xp(1) = norm2(x(1:2))
  xp(2) = atan2(x(2),x(1))
  xp(3) = x(3)
 end function carcyl

 pure function softened_pot(r,hsoft) result(phi)
! Softened potential from Price & Monaghan 2007 Appendix A
  implicit none
  real(8),intent(in):: r,hsoft
  real(8):: phi,h,q
  h = 0.5d0*max(hsoft,r*0.5d0,1d-99)
  q=r/h
  if(q>=2d0)then
   phi = -1d0/r
  elseif(q>=1d0)then
   phi = (4d0/3d0*q**2-q**3+0.3d0*q**4-q**5/30d0-1.6d0+1d0/(15d0*q))/h
  elseif(q>=0d0)then
   phi = (2d0/3d0*q**2-0.3d0*q**4+0.1d0*q**5-1.4d0)/h
  end if

 end function softened_pot

 pure function softened_acc(r,hsoft) result(acc)
! Softened acceleration from Price & Monaghan 2007 Appendix A
  implicit none
  real(8),intent(in):: r,hsoft
  real(8)::acc,h,q
  h = 0.5d0*max(hsoft,r*0.5d0,1d-99)
  q = r/h
  if(q<1d0)then
   acc=(4d0/3d0*q-1.2d0*q**3+0.5d0*q**4)/h**2
  elseif(q<2d0)then
   acc=(8d0/3d0*q-3d0*q**2+1.2d0*q**3-q**4/6d0-1d0/(15d0*q**2))/h**2
  else
   acc=1d0/r**2
  end if

 end function softened_acc


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VCAR
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the Cartesian vector components from polar coordinates

subroutine get_vcar(xcar,x3,v1,v2,v3,vcar)
 implicit none

 real(8),intent(in):: xcar(1:3),x3,v1,v2,v3
 real(8),intent(out)::vcar(1:3)
 real(8),dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = get_uvec2(uvec1, x3)
 uvec3 = cross(uvec1,uvec2)

 vcar = v1*uvec1 + v2*uvec2 + v3*uvec3

return
end subroutine get_vcar

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VPOL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the Polar vector components from Cartesian coordinates

subroutine get_vpol(xcar,x3,vcar,v1,v2,v3)
 implicit none

 real(8),intent(in):: xcar(1:3),x3,vcar(1:3)
 real(8),intent(out)::v1,v2,v3
 real(8),dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = get_uvec2(uvec1, x3)
 uvec3 = cross(uvec1,uvec2)

 v1 = dot_product(vcar,uvec1)
 v2 = dot_product(vcar,uvec2)
 v3 = dot_product(vcar,uvec3)

return
end subroutine get_vpol

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VCYL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate Cylindrical vector components from Cartesian coordinates

subroutine get_vcyl(xcar,vcar,v1,v2,v3)

 real(8),intent(in):: xcar(1:3),vcar(1:3)
 real(8),intent(out)::v1,v2,v3
 real(8),dimension(1:3)::uvec1,uvec2,uvec3
 real(8):: r

!-----------------------------------------------------------------------------

 r = norm2(xcar(1:2))
 uvec1 = [ xcar(1)/r,xcar(2)/r,0d0]
 uvec2 = [-xcar(2)/r,xcar(1)/r,0d0]
 uvec3 = [0d0,0d0,1d0]

 v1 = dot_product(vcar,uvec1)
 v2 = dot_product(vcar,uvec2)
 v3 = vcar(3)

return
end subroutine get_vcyl

! This function is equivalent to rotz(roty(rotz(uvec1,-x3),0.5*pi),x3)
! but expanded out, to reduce roundoff error
 pure function get_uvec2(x,theta) result(xp)
  implicit none
  real(8),intent(in):: x(1:3), theta
  real(8):: xp(1:3)
  xp(1) = sin(theta)**2*x(1) - 0.5d0*sin(2d0*theta)*x(2) + cos(theta)*x(3)
  xp(2) = -0.5d0*sin(2d0*theta)*x(1) + cos(theta)**2*x(2) + sin(theta)*x(3)
  xp(3) = -cos(theta)*x(1) - sin(theta)*x(2)
 end function get_uvec2

! rotate about the x axis
 pure function rotx(x,theta) result(xp)
  implicit none
  real(8),intent(in):: x(1:3), theta
  real(8):: xp(1:3)
  xp(1) = x(1)
  xp(2) = cos(theta)*x(2) - sin(theta)*x(3)
  xp(3) = sin(theta)*x(2) + cos(theta)*x(3)
 end function rotx

! rotate about the y axis
 pure function roty(x,theta) result(xp)
  implicit none
  real(8),intent(in):: x(1:3), theta
  real(8):: xp(1:3)
  xp(1) = cos(theta)*x(1) + sin(theta)*x(3)
  xp(2) = x(2)
  xp(3) =-sin(theta)*x(1) + cos(theta)*x(3)
 end function roty

! rotate about the z axis
 pure function rotz(x,theta) result(xp)
  implicit none
  real(8),intent(in):: x(1:3), theta
  real(8):: xp(1:3)
  xp(1) = cos(theta)*x(1) - sin(theta)*x(2)
  xp(2) = sin(theta)*x(1) + cos(theta)*x(2)
  xp(3) = x(3)
 end function rotz

! cross product
 pure function cross(x,y) result(z)
  implicit none
  real(8),intent(in):: x(1:3),y(1:3)
  real(8):: z(1:3)
  z(1) = x(2)*y(3)-x(3)*y(2)
  z(2) = x(3)*y(1)-x(1)*y(3)
  z(3) = x(1)*y(2)-x(2)*y(1)
 end function cross

 ! linear interpolation
 pure function intpol(x,y,z) result(val)
  implicit none
  real(8),intent(in):: x(1:2), y(1:2)
  real(8),intent(in):: z
  real(8):: val
  val = ((x(2)-z)*y(1)+(z-x(1))*y(2))/(x(2)-x(1))
 end function intpol

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE GEOMETRICAL_SERIES
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate dxi's in a geometrical series.

subroutine geometrical_series(dxi,xmin,is,ie,xis,xie)

 implicit none

 integer,intent(in):: is,ie
 real(8),intent(in):: xis,xie,xmin
 real(8),intent(inout),allocatable:: dxi(:)
 integer:: i
 real(8):: xrng, irng, xr, xrnew, xrmax, maxerr, fx, dfx

!-----------------------------------------------------------------------------

 xrmax = 1.05d0
 maxerr = 1d-10

 xr = 1.01d0
 xrng = xie - xis ; irng = dble(ie - is + 1)

 if(xrng/irng<xmin)then
  print *,"Error from geometrical_series ;"
  print *,"xmin should be smaller or uniform mesh should be chosen",xmin
  stop
 end if

 do i = 1, 10000000
  fx = (xr-1d0)*xrng - xmin * (xr**irng-1d0)
  dfx = xrng - irng * xmin * xr**(irng-1d0)

  xrnew = xr - fx/dfx

  if(abs((xrnew-xr)/xr)<maxerr)then
   xr = xrnew ; exit
  end if
  if(xrnew<1d0)xrnew = 2d0

  xr = xrnew
 end do

 if(xr>xrmax)then
  print *,"xmin too small", xmin, xr
  stop
 end if

 dxi(is) = xmin
 do i = is+1, ie
  dxi(i) = dxi(i-1) * xr
 end do
 dxi(is-1) = dxi(is) ; dxi(is-2) = dxi(is+1)
 dxi(ie+1) = dxi(ie)*xr ; dxi(ie+2) = dxi(ie)*xr*xr

 if(xr-1d0<maxerr) dxi = (xie-xis) / irng

return
end subroutine geometrical_series


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GRAVPOT1D
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate 1D gravitational potential

subroutine gravpot1d
 use grid
 use constants,only:pi,G
 use physval,only:d
 use gravmod,only:grvphi,mc
 use mpi_utils,only:allreduce_mpi
 real(8) :: num_n, denom_n, ishell
 integer:: i,n,j,k

 ! Loop over each shell
 do i = is_global, ie_global-1

  ! Add up contributions from all exterior shells
  ishell = 0d0
  do n = i+1, ie_global

    num_n = 0d0
    denom_n = 0d0

    ! Loop over each cell in the shell, adding to counters if the cell is in the current mpi task
    if (is<=n .and. n<=ie) then
      do j = js, je
        do k = ks, ke
          num_n = num_n + d(n,j,k)*dvol(n,j,k)
          denom_n = denom_n + dvol(n,j,k)
        end do
      end do
    end if

    ! Add up counters across tasks
    call allreduce_mpi('sum',num_n)
    call allreduce_mpi('sum',denom_n)

    ! Add the contribution from shell n
    ! (This should now be the same value on each MPI task...)
    ishell = ishell - num_n/denom_n * x1(n)*dxi1(n)

  end do

  ! If the cell(s) is in the current MPI task, update grvphi.
  if (is<=i .and. i<=ie) then
    grvphi(i,js-2:je+2,ks-2:ke+2) = G*(-mc(i)/x1(i)+4d0*pi*ishell)
  endif

 end do

 ! Outermost shell
 if (ie==ie_global) then
    grvphi(ie,js-2:je+2,ks-2:ke+2) = -G*mc(ie)/x1(ie)
 end if

 return
end subroutine gravpot1d

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE MASSCOORDINATE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get mass coordinates

subroutine masscoordinate

 use settings,only:eq_sym
 use grid,only:is_global,ie_global,is,ie,js,je,ks,ke,dvol
 use physval,only:d
 use gravmod,only:mc
 use mpi_utils,only:allreduce_mpi

 implicit none

 integer:: i,j,k
 real(8):: fac, shellmass

!-----------------------------------------------------------------------------

 if(eq_sym)then
  fac=2d0
 else
  fac=1d0
 end if

 do i = is_global, ie_global
  shellmass = 0d0
  if (is<=i .and. i<=ie) then
!$omp parallel do private(j,k) reduction(+:shellmass) collapse(2)
   do k = ks, ke
    do j = js, je
     shellmass = shellmass + d(i,j,k)*dvol(i,j,k)
    end do
   end do
!$omp end parallel do
  end if
  call allreduce_mpi('sum',shellmass)
  mc(i) = mc(i-1) + fac*shellmass
 end do

return
end subroutine masscoordinate

logical function isequal(a,b)
  implicit none
  real(8),intent(in):: a,b
  isequal = abs(a-b) < epsilon(real(0.,kind=8))
end function isequal

end module utils
