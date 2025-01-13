module utils_analysis
  implicit none

 contains

! electron scattering opacity
 pure function kap_es(x)
  real(8),intent(in):: x
  real(8):: kap_es
  kap_es = 0.2d0*(1d0+x)
 end function kap_es

! simplified opacity function
 elemental function kap(X,T,nu)
  real(8),intent(in):: X,T
  real(8),intent(in),optional:: nu
  real(8):: x_ion
  real(8),parameter:: Trec = 5d3
  real(8),parameter:: opacity_floor = 1d-2 ! based on Bersten+2011
  real(8):: kap
  if(present(nu))then
! frequency-dependent opacity (same as grey opacity for now)
   x_ion = 1d0/(1d0+(max(T,1d1)/Trec)**(-11))
   kap = max(kap_es(X)*x_ion,opacity_floor)
  else
   ! grey opacity
   x_ion = 1d0/(1d0+(max(T,1d1)/Trec)**(-11))
   kap = max(kap_es(X)*x_ion,opacity_floor)
  end if
 end function kap

 ! Planck function per frequency
 elemental function planck_nu(nu,T) result(B)
  use constants,only:k=>kbol,c=>clight,h=>hplanck
  real(8),intent(in):: nu,T
  real(8):: B
  B = 2d0*h*nu**3/c**2/(exp(h*nu/(k*T))-1d0)
 end function planck_nu

 ! Planck function per wavelength
 elemental function planck_lambda(lambda,T) result(B)
  use constants,only:k=>kbol,c=>clight,h=>hplanck
  real(8),intent(in):: lambda,T
  real(8):: B
  B = 2d0*h*c**2/lambda**5/(exp(h*c/(lambda*k*max(T,1d2)))-1d0)
 end function planck_lambda

 subroutine fit_blackbody(lambda, spec, T_eff)
  real(8),allocatable,intent(in) :: lambda(:),spec(:)
  real(8),intent(out) :: T_eff
  ! Variables for fitting
  real(8) :: T_min, T_max, T_mid
  real(8) :: chi2_min, chi2_mid, chi2_max
  real(8) :: step, chi2, scale_factor
  integer :: i, j

  ! Initialize temperature range (adjust as needed)
  T_min = 1000.0d0      ! Minimum temperature to search (Kelvin)
  T_max = 500000.0d0     ! Maximum temperature to search (Kelvin)
  step = 100.0d0         ! Initial step size (Kelvin)

  ! Minimize chi-squared
  chi2_min = huge(0d0)  ! Large initial value
  T_eff = T_min

  do while (step > 0.1d0)
   do i = 0, ceiling((T_max - T_min) / step)
    T_mid = T_min + i * step
    chi2 = 0.0d0

! Compute scale factor for normalization
    scale_factor = sum(spec)/sum(planck_lambda(lambda,T_mid))

! Compute chi-squared with normalization
    do j = 1, size(lambda)
     chi2 = chi2 + (spec(j) - scale_factor*planck_lambda(lambda(j), T_mid))**2
    end do

    if (chi2 < chi2_min) then
     chi2_min = chi2
     T_eff = T_mid
    end if

   end do

! Refine the search range
   T_min = max(T_min, T_eff - step)
   T_max = min(T_max, T_eff + step)
   step = step / 2.0d0
  end do

 end subroutine fit_blackbody

 function get_bbradius(Teff,Lum) result(Reff)
  use constants,only:pi,sigma
  real(8),intent(in):: Teff,Lum
  real(8):: Reff
  Reff = sqrt(Lum/(4d0*pi*sigma*Teff**4))
 end function get_bbradius

end module utils_analysis
