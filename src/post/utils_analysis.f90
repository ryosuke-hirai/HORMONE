module utils_analysis
  implicit none

 contains

! electron scattering opacity
 pure function kap_es(x)
  real(8),intent(in):: x
  real(8):: kap_es
  kap_es = 0.2d0*(1d0+x)
 end function kap_es

! electron scattering opacity
 function kap(X,T)
  real(8),intent(in):: X,T
!  real(8):: x_ion
!  real(8),parameter:: Trec = 1.2d4
  real(8):: kap
!  x_ion = 1d0/(1d0+(T/Trec)**(-11))
  kap = kap_es(X)!*x_ion+1d-4
 end function kap

end module utils_analysis
