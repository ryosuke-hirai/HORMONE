module opacity_mod
 implicit none

 real(8),public:: c_kappa_p, c_kappa_r, c_kappa_f, opacity_floor

contains

! Thomson scattering opacity
 pure elemental function kap_thomson(X)
  real(8),intent(in):: X
  real(8):: kap_thomson
  kap_thomson = 0.2d0*(1d0+X)
 end function kap_thomson

! electron scattering opacity (Klein-Nishina correction)
 pure elemental function kap_es(X,d,T)
  real(8),intent(in):: X,d,T
  real(8):: kap_es
  kap_es = kap_thomson(X) / (1d0+2.7d11*d/T**2)/(1d0+(T/4.5d8)**0.86d0)
 end function kap_es

! Kramers opacity (free-free + bound-free + bound-bound)
 pure elemental function kap_kramers(X,Z,d,T)
  real(8),intent(in):: X,Z,d,T
  real(8):: kap_kramers
  kap_kramers = 4d25*(1d0+X)*(Z+1d-3)*d/T**3.5d0
 end function kap_kramers

! fit to H- opacity
 pure elemental function kap_negH(Z,d,T)
  real(8),intent(in):: Z,d,T
  real(8):: kap_negH
!  kap_negH = 1.1d-25*sqrt(Z*d)*T**7.7d0 ! Metzger & Pejcha 2017
  kap_negH = 1.1d-40*sqrt(Z)*d**0.2d0*T**11 ! Hirai+2025
 end function kap_negH

! approximate molecular opacity
 pure elemental function kap_mol(Z)
  real(8),intent(in):: Z
  real(8):: kap_mol
  kap_mol = 0.1d0*Z
 end function kap_mol

! approximate radiative opacity
 pure elemental function kap_rad(X,Z,d,T)
  real(8),intent(in):: X,Z,d,T
  real(8):: kap_rad
  kap_rad = kap_mol(Z) &
          + 1d0 / ( 1d0/kap_negH(Z,d,T) &
                  + 1d0/(kap_es(X,d,T)+kap_kramers(X,Z,d,T)+kap_hline(d,T)) )
 end function kap_rad

! approximate conductive opacity
 pure elemental function kap_cond(Z,d,T)
  real(8),intent(in):: Z,d,T
  real(8):: kap_cond
  kap_cond = 2.6d-7*Z*(T/d)**2*(1d0+(d/2d6)**(2d0/3d0))
 end function kap_cond

! (TEST FEATURE, DO NOT USE) attempt at adding hydrogen line effect
 pure elemental function kap_hline(d,T)
  real(8),intent(in):: d,T
  real(8):: kap_hline, logT_h0
  logT_h0 = log10(d)/40+4.+11/40.
  kap_hline=1d2*d**0.2d0/(1d0+(2d0*(log10(T)-logT_h0)/0.15d0)**4)
 end function kap_hline

! approximate opacity function
 elemental function kap_approx(X,Z,d,T) result(kap)
  real(8),intent(in):: X,Z,d,T
  real(8):: kap
! grey opacity
  kap = max(1d0/(1d0/kap_rad(X,Z,d,T)+1d0/kap_cond(Z,d,T)),opacity_floor)
 end function kap_approx

! Get Rosseland mean opacity
 function kappa_r(X,Z,d,T) result(kappa)
  use settings,only:opacitytype
  real(8),intent(in)::X,Z,d,T
  real(8):: kappa

  select case(opacitytype)
  case(0) ! fixed opacity
   kappa = c_kappa_r
  case(1) ! approximate analytical formula
   kappa = kap_approx(X,Z,d,T)
  end select

 end function kappa_r

! Get Planck mean opacity
 function kappa_p(X,Z,d,T) result(kappa)
  use settings,only:opacitytype
  real(8),intent(in)::X,Z,d,T
  real(8):: kappa

  select case(opacitytype)
  case(0) ! fixed opacity
   kappa = c_kappa_p
  case(1) ! approximate analytical formula
   kappa = kap_approx(X,Z,d,T)
  end select

 end function kappa_p

! Get Flux mean opacity
 function kappa_f(X,Z,d,T) result(kappa)
  use settings,only:opacitytype
  real(8),intent(in)::X,Z,d,T
  real(8):: kappa

  select case(opacitytype)
  case(0) ! fixed opacity
   kappa = c_kappa_f
  case(1) ! approximate analytical formula
   kappa = kap_approx(X,Z,d,T)
  end select

 end function kappa_f

end module opacity_mod
