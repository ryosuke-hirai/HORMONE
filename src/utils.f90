module utils
 implicit none

contains

! convert polar to cartesian coordinates
 pure function polcar(xp) result(x)
  implicit none
  real(8),intent(in):: xp(1:3)
  real(8):: x(1:3)
  x(1) = xp(1)*sin(xp(2))*cos(xp(3))
  x(2) = xp(1)*sin(xp(2))*sin(xp(3))
  x(3) = xp(1)*cos(xp(2))
 end function polcar

! convert cartesian to polar coordinates
 pure function carpol(x) result(xp)
  implicit none
  real(8),intent(in):: x(1:3)
  real(8):: xp(1:3)
  real(8),parameter:: tiny=1d-99
  xp(1) = norm2(x)
  xp(2) = acos(x(3)/max(xp(1),tiny))
  xp(3) = atan2(x(2),x(1))
 end function carpol


 pure function softened_pot(r,hsoft) result(phi)
! Softened potential from Price & Monaghan 2007 Appendix A
  implicit none
  real(8),intent(in):: r,hsoft
  real(8):: phi,h,q
  h = 0.5d0*max(hsoft,1d-99)
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
  h = 0.5d0*max(hsoft,1d-99)
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

 use constants,only:pi

 implicit none

 real(8),intent(in):: xcar(1:3),x3,v1,v2,v3
 real(8),intent(out)::vcar(1:3)
 real(8),dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = rotz(roty(rotz(uvec1,-x3),0.5*pi),x3)
 uvec3 = cross(uvec1,uvec2)

 vcar = v1*uvec1 + v2*uvec2 + v3*uvec3

return
end subroutine get_vcar


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_VPOL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the Polar vector components from Cartesian coordinates

subroutine get_vpol(xcar,x3,vcar,v1,v2,v3)

 use constants,only:pi

 implicit none

 real(8),intent(in):: xcar(1:3),x3,vcar(1:3)
 real(8),intent(out)::v1,v2,v3
 real(8),dimension(1:3)::uvec1,uvec2,uvec3

!-----------------------------------------------------------------------------

 uvec1 = xcar/norm2(xcar)
 uvec2 = rotz(roty(rotz(uvec1,-x3),0.5*pi),x3)
 uvec3 = cross(uvec1,uvec2)

 v1 = dot_product(vcar,uvec1)
 v2 = dot_product(vcar,uvec2)
 v3 = dot_product(vcar,uvec3)

return
end subroutine get_vpol

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
    grvphi(i,js:je,ks:ke) = G*(-mc(i)/x1(i)+4d0*pi*ishell)
  endif

 end do

 ! Outermost shell
 if (ie==ie_global) then
    grvphi(ie,js:je,ks:ke) = -G*mc(ie)/x1(ie)
 end if

 return
end subroutine gravpot1d

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE MASSCOORDINATE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get mass coordinates

subroutine masscoordinate

 use settings,only:eq_sym
 use grid,only:is,ie,js,je,ks,ke,dvol
 use physval,only:d
 use gravmod,only:mc

 implicit none

 integer:: i,j,k
 real(8):: fac, shellmass

!-----------------------------------------------------------------------------

 if(eq_sym)then
  fac=2d0
 else
  fac=1d0
 end if

!$omp parallel
 do i = is, ie
  shellmass = 0d0
!$omp do private(j,k) reduction(+:shellmass)
  do k = ks, ke
   do j = js, je
    shellmass = shellmass + d(i,j,k)*dvol(i,j,k)
   end do
  end do
!$omp end do
!$omp single
  mc(i) = mc(i-1) + fac*shellmass
!$omp end single
 end do
!$omp end parallel

return
end subroutine masscoordinate

logical function isequal(a,b)
  implicit none
  real(8),intent(in):: a,b
  isequal = abs(a-b) < epsilon(real(0.,kind=8))
end function isequal

end module utils
