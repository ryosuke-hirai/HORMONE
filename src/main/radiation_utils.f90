module radiation_utils
 implicit none

 real(8),allocatable,dimension(:,:,:,:),public:: geo

contains

! Get flux limiter for flux-limited diffusion
function lambda(R)
 use settings,only:lambdatype
 real(8),intent(in)::R
 real(8):: lambda
 select case(lambdatype)
 case(0) ! Eddington approximation
  lambda = 1d0/3d0
 case(1) ! Levermore & Pomraning 1981
  lambda = lambda_LP81(R)
 case(2) ! Minerbo 1978
  lambda = lambda_M78(R)
 case(3) ! Kley 1989
  lambda = lambda_K89(R)
 case default
  print*,'Error from choice of radiation flux limiter'
  print*,'lambdatype=',lambdatype
  stop
 end select
end function lambda

! Levermore & Pomraning 1981
pure function lambda_LP81(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 real(8),parameter:: eps=1d-4,third=1d0/3d0
 if(R<eps)then
  lambda = ((1d0/tanh(eps)-1d0/eps)/eps-third)*R/eps + third
 else
  lambda = (1d0/tanh(R)-1d0/R)/R
 end if
end function lambda_LP81

! Minerbo 1978
pure function lambda_M78(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 if(R<=1.5d0)then
  lambda = 2d0/(3d0+sqrt(9d0+12d0*R**2))
 else
  lambda = 1d0/(1d0+R+sqrt(1d0+2d0*R))
 end if
end function lambda_M78

! Kley 1989
pure function lambda_K89(R) result(lambda)
 real(8),intent(in):: R
 real(8):: lambda
 if(R<=2d0)then
  lambda = 2d0/(3d0+sqrt(9d0+10d0*R**2))
 else
  lambda = 10d0/(10d0*R+9d0+sqrt(180d0*R+81d0))
 end if
end function lambda_K89

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE GET_GEO
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute geometrical coefficients for diffusion equation

subroutine get_geo

 use grid,only:is,ie,js,je,ks,ke,sa1,sa2,sa3,dx1,dx2,dx3,g22,g33

 integer:: i,j,k

!-----------------------------------------------------------------------------

 allocate(geo(1:3,is-1:ie,js-1:je,ks-1:ke))

!$omp parallel do private(i,j,k) collapse(3)
 do k = ks-1, ke
  do j = js-1, je
   do i = is-1, ie
    geo(1,i,j,k) = sa1(i,j,k) /  dx1(i+1)
    geo(2,i,j,k) = sa2(i,j,k) / (dx2(j+1)*g22(i))
    geo(3,i,j,k) = sa3(i,j,k) / (dx3(k+1)*g33(i,j))
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_geo

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SOLVE_QUARTIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To solve quartic equation required for rad_heat_cool

subroutine solve_quartic(c1,c2,x)

 use constants,only:huge

 integer:: n
 real(8),intent(in):: c1,c2
 real(8),intent(inout):: x
 real(8):: corr, xdot, dt
 real(8),parameter:: W4err = 1d-2, raderr = 1d-15

!-----------------------------------------------------------------------------

 corr=huge ; xdot=0d0 ; dt=0.1d0
 do n = 1, 500
  corr = (x**4+c1*x+c2) / (4d0*x**3+c1)
  if(abs(corr)>W4err*x)then
   x = x + xdot*dt
   xdot = (1d0-2d0*dt)*xdot - dt*corr
  else
   x = x-corr
   xdot = 0d0
  end if
  if(abs(corr)<raderr*x)exit
  if(n>50)dt=0.5d0
 end do

 if(n>500)then
  print*,'Error in solve_quartic in radiation.f90'
  stop
 end if

return
end subroutine solve_quartic


function update_Tgas(X,Z,d,erad,T,dt) result(Tnew)
! PURPOSE: Update Tgas based on linear approximation
 use constants,only:a=>arad,c=>clight,Cv
 use opacity_mod,only:kappa_p

 real(8),intent(in):: X,Z,d,erad,T,dt
 real(8):: Tnew,kappap

!-----------------------------------------------------------------------------

 kappap = kappa_p(X,Z,d,T)

 Tnew = (kappap*c*(3d0*a*T**4+erad)*dt+Cv*T)/(Cv+4d0*kappap*c*a*T**3*dt)

end function update_Tgas

end module radiation_utils
