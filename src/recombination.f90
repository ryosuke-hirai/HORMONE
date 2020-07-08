module recombination_mod
 implicit none

 real*8,allocatable,public,dimension(:):: eion
 real*8,allocatable,private,dimension(:):: logeion, &
                                   arec, brec, crec, drec, arec1c, brec1c
 real*8,private:: frec, edge, tanh_c, dtanh_c
 real*8,parameter,private:: tanh_edge = 3.6467d0

 interface rapid_tanh
    module procedure rapid_tanhs,rapid_tanhv
 end interface rapid_tanh
 interface rapid_dtanh
    module procedure rapid_dtanhs,rapid_dtanhv
 end interface rapid_dtanh
 
contains
! **************************************************************************
 function rapid_tanhs(x)
  implicit none
  real*8,intent(in)::x
  real*8:: a,b,x2,rapid_tanhs
  
  if(abs(x)>=tanh_edge)then
   rapid_tanhs = sign(1d0,x)-tanh_c/x
  else
   x2 = x*x
   a = ((x2+105d0)*x2+945d0)*x
   b = ((x2+28d0)*x2+63d0)*15d0
   rapid_tanhs = a/b
  end if
 end function rapid_tanhs

 function rapid_tanhv(x)
  implicit none
  real*8,intent(in)::x(:)
  real*8,dimension(size(x)):: rapid_tanhv
  real*8:: a,b,x2
  integer:: i
  
  do i = 1, size(x)
   if(abs(x(i))>=tanh_edge)then
    rapid_tanhv(i) = sign(1d0,x(i))-tanh_c/x(i)
   else
    x2 = x(i)*x(i)
    a = ((x2+105d0)*x2+945d0)*x(i)
    b = ((x2+28d0)*x2+63d0)*15d0
    rapid_tanhv(i) = a/b
   end if
  end do
 end function rapid_tanhv
! **************************************************************************
 function rapid_dtanhs(x)
  implicit none
  real*8,intent(in)::x
  real*8:: a,b,x2,rapid_dtanhs
  
  if(abs(x)>=tanh_edge)then
   rapid_dtanhs = dtanh_c/x**2
  else
   x2=x**2
   a = (((x2-21d0)*x2+420d0)*x2-6615d0)*x2+59535d0
   b = (x2+28d0)*x2+63d0
   rapid_dtanhs = a/(15d0*b**2)
  end if
 end function rapid_dtanhs
 
 function rapid_dtanhv(x)
  implicit none
  real*8,intent(in)::x(:)
  real*8,dimension(size(x)):: rapid_dtanhv
  real*8:: a,b,x2
  integer:: i
  
  do i = 1, size(x)
   if(abs(x(i))>=tanh_edge)then
    rapid_dtanhv(i) = dtanh_c/(x(i)*x(i))
   else
    x2=x(i)**2
    a = (((x2-21d0)*x2+420d0)*x2-6615d0)*x2+59535d0
    b = (x2+28d0)*x2+63d0
    rapid_dtanhv(i) = a/(15d0*b**2)
   end if
  end do
 end function rapid_dtanhv
! **************************************************************************
 subroutine recombination_setup
  use constants,only:amu,kbol
  implicit none
  real*8:: x
  allocate( eion(1:4), arec(2:4), brec(2:4), crec(1:4), drec(1:4),&
            arec1c(1:2), brec1c(1:2), logeion(1:4) )

  eion(1) = 4.36d12   ! H2   [erg/mol]
  eion(2) = 1.312d13  ! HI   [erg/mol]
  eion(3) = 2.3723d13 ! HeI  [erg/mol]
  eion(4) = 5.2505d13 ! HeII [erg/mol]
  logeion(1:4) = log10(eion(1:4)*amu/kbol)

! These fitting parameters are tuned for eosDT in the range X=0.6-0.8
  frec = 0.005d0
  arec(2:4) = (/ 0.821d0, 0.829d0, 0.846d0 /)
  brec(2:4) = (/ 0.055d0, 0.055d0, 0.055d0 /)
  crec(1:4) = (/ 0.02d0, 0.025d0, 0.015d0, 0.015d0 /)
  drec(1:4) = (/ 0.05d0, 0.05d0, 0.05d0, 0.05d0 /)

  arec1c(1:2) = (/ log10(3500d0)/logeion(1), 0.753d0 /)
  brec1c(1:2) = (/ 0d0, 0.055d0 /)
  edge = -9.4d0

! Parameter for rapid tanh function
  x=tanh_edge
  tanh_c = x*(1d0-((x*x+105d0)*x*x+945d0)*x/(((x*x+28d0)*x*x+63d0)*15d0))
  dtanh_c= x*x*((((x*x-21d0)*x*x+420d0)*x*x-6615d0)*x*x+59535d0)&
               / (15d0*((x*x+28d0)*x*x+63d0)**2d0)

  return
 end subroutine recombination_setup
! ***************************************************************************
 real*8 function arec1(x)
! molecular hydrogen fit coefficients
  implicit none
  real*8,intent(in):: x ! logd

  if(x<edge)then
   arec1 = arec1c(1)
  else
   arec1 = arec1c(2)
  end if
 end function arec1
! ***************************************************************************
 real*8 function brec1(x)
! molecular hydrogen fit coefficients
  implicit none
  real*8,intent(in):: x

  if(x<edge)then
   brec1 = brec1c(1)
  else
   brec1 = brec1c(2)
  end if
 end function brec1
! *************************************************************************** 
 subroutine get_erec_imurec(logd,T,X,Y,erec,imurec,derecdT,dimurecdT)
  implicit none
  real*8,intent(in):: logd,T,X,Y
  real*8,intent(out):: erec,derecdT,imurec,dimurecdT
  real*8:: logQ, logT, Yfac
  real*8,dimension(1:4):: e, Ttra, width, xi, arg, zi

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.

  e(1) = eion(1)*X*0.5d0
  e(2) = eion(2)*X
  e(3) = eion(3)*Y*0.25d0
  e(4) = eion(4)*Y*0.25d0

  logT = log10(T)
  logQ = max(-14d0,logd)-2d0*logT+12d0

  Yfac=(1d0-frec*Y)

! calculate transition temperature for given logQ
  Ttra (1)   = arec1(logd)*logeion(1)+brec1(logd)*logQ
  Ttra (2:4) = Yfac*arec(2:4)*logeion(2:4)+brec(2:4)*logQ
  width(1:4) = Ttra(1:4)*crec(1:4)*(1d0+drec(1:4)*logQ)
  arg(1:4) = (logT-Ttra(1:4))/width(1:4)

! calculate ionisation fraction and derivatives
  xi(1:4) = 0.5d0*(rapid_tanh(arg(1:4))+1d0)
  zi(1) = ( width(1)*(1d0+2d0*brec1(logd)) &
            + 2d0*crec(1)*(logT-Ttra(1))&
              *(brec1(logd)*(1d0+drec(1)*logQ)+drec(1)*Ttra(1)) )&
          / (2d0*T*width(1)*width(1)) * rapid_dtanh(arg(1))
  zi(2:4) = ( width(2:4)*(1d0+2d0*brec(2:4)) &
            + 2d0*crec(2:4)*(logT-Ttra(2:4))&
              *(brec(2:4)*(1d0+drec(2:4)*logQ)+drec(2:4)*Ttra(2:4)) )&
          / (2d0*T*width(2:4)*width(2:4)) * rapid_dtanh(arg(2:4))

  erec = sum(e(1:4)*xi(1:4))
  imurec = (0.5d0*xi(1)+xi(2))*X+0.25d0*(xi(3)+xi(4)-1d0)*Y+0.5d0

  derecdT = sum(e(1:4)*zi(1:4))
  dimurecdT = (0.5d0*zi(1)+zi(2))*X+0.25d0*(zi(3)+zi(4))*Y

  return
 end subroutine get_erec_imurec

! ***************************************************************************
 
 subroutine get_imurec(logd,T,X,Y,imurec,dimurecdT)
  ! Calculate the mean molecular weight for partially ionised plasma
  implicit none
  real*8,intent(in):: logd,T,X,Y
  real*8,intent(out):: imurec,dimurecdT
  real*8:: logQ, logT, Yfac
  real*8,dimension(1:4):: Ttra, width, xi, zi, arg

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.

  logT = log10(T)
  logQ = max(-14d0,logd)-2d0*logT+12d0

  Yfac=(1d0-frec*Y)

  Ttra (1)   = arec1(logd)*logeion(1)+brec1(logd)*logQ
  Ttra (2:4) = Yfac*arec(2:4)*logeion(2:4)+brec(2:4)*logQ
  width(1:4) = Ttra(1:4)*crec(1:4)*(1d0+drec(1:4)*logQ)
  arg(1:4) = (logT-Ttra(1:4))/width(1:4)

  xi(1:4) = 0.5d0*(rapid_tanh(arg(1:4))+1d0)
  zi(1) = ( width(1)*(1d0+2d0*brec1(logd)) &
            + 2d0*crec(1)*(logT-Ttra(1))&
              *(brec1(logd)*(1d0+drec(1)*logQ)+drec(1)*Ttra(1)) )&
          / (2d0*T*width(1)*width(1)) * rapid_dtanh(arg(1))
  zi(2:4) = ( width(2:4)*(1d0+2d0*brec(2:4)) &
            + 2d0*crec(2:4)*(logT-Ttra(2:4))&
              *(brec(2:4)*(1d0+drec(2:4)*logQ)+drec(2:4)*Ttra(2:4)) )&
          / (2d0*T*width(2:4)*width(2:4)) * rapid_dtanh(arg(2:4))

  imurec = (0.5d0*xi(1)+xi(2))*X+0.25d0*(xi(3)+xi(4)-1d0)*Y+0.5d0

  dimurecdT = (0.5d0*zi(1)+zi(2))*X+0.25d0*(zi(3)+zi(4))*Y
  
  return
 end subroutine get_imurec

! ***************************************************************************
 real*8 function get_erec(logd,T,X,Y)
  implicit none
  real*8,intent(in):: logd,T,X,Y
  real*8:: logQ, logT, Yfac
  real*8,dimension(1:4):: e, Ttra, width, xi, arg, efac

! CAUTION: This is only a poor man's way of implementing recombination energy.
!          It only should be used for -3.5<logQ<-6 where logQ=logrho-2logT+12.
  
  e(1) = eion(1)*X*0.5d0
  e(2) = eion(2)*X
  e(3) = eion(3)*Y*0.25d0
  e(4) = eion(4)*Y*0.25d0

  logT = log10(T)
  logQ = max(-14d0,logd)-2d0*logT+12d0

  Yfac=(1d0-frec*Y)

  Ttra (1)   = arec1(logd)*logeion(1)+brec1(logd)*logQ
  Ttra (2:4) = Yfac*arec(2:4)*logeion(2:4)+brec(2:4)*logQ
  width(1:4) = Ttra(1:4)*crec(1:4)*(1d0+drec(1:4)*logQ)
  arg(1:4) =(logT-Ttra(1:4))/width(1:4)

  xi(1:4) = 0.5d0*(rapid_tanh(arg(1:4))+1d0)

  get_erec = sum(e(1:4)*xi(1:4))

  return
 end function get_erec

end module recombination_mod
