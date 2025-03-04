module opacity_mod
 implicit none

 real(8),public:: c_kappa_p, c_kappa_r, c_kappa_f

contains

! Get Rosseland mean opacity
function kappa_r(d,T) result(kappa)
 use settings,only:opacitytype
 real(8),intent(in)::d,T
 real(8):: kappa

 select case(opacitytype)
 case(0)
  kappa = c_kappa_r
 end select

 ! TEMPORARY: Suppress warnings for function under construction
 if (.false.) then
  if (d>0.d0 .and. T>0.d0) continue
 end if
end function kappa_r

! Get Planck mean opacity
function kappa_p(d,T) result(kappa)
 use settings,only:opacitytype
 real(8),intent(in)::d,T
 real(8):: kappa

 select case(opacitytype)
 case(0)
  kappa = c_kappa_p
 end select

 ! TEMPORARY: Suppress warnings for function under construction
 if (.false.) then
  if (d>0.d0 .and. T>0.d0) continue
 end if
end function kappa_p

! Get Flux mean opacity
function kappa_f(d,T) result(kappa)
 use settings,only:opacitytype
 real(8),intent(in)::d,T
 real(8):: kappa

 select case(opacitytype)
 case(0)
  kappa = c_kappa_f
 end select

 ! TEMPORARY: Suppress warnings for function under construction
 if (.false.) then
  if (d>0.d0 .and. T>0.d0) continue
 end if
end function kappa_f

end module opacity_mod
