module derived_types

 implicit none

 type sink_prop
  sequence
  integer:: i,j,k,pad ! pad is added to align memory for 64-bit machines
  real(8):: mass, softfac, lsoft, laccr, locres, dt, mdot, racc, facc, jet_ang
  real(8),dimension(1:3):: x,v,a,xpol,Jspin,jdot,jet_dir
 end type sink_prop

contains

 subroutine null_sink(sink)
  type(sink_prop),intent(out):: sink
  sink%i=0
  sink%j=0
  sink%k=0
  sink%pad=0
  sink%mass=0d0
  sink%softfac=0d0
  sink%lsoft=0d0
  sink%laccr=0d0
  sink%locres=0d0
  sink%dt=0d0
  sink%mdot=0d0
  sink%racc=0d0
  sink%facc=0d0
  sink%jet_ang=0d0
  sink%x=0d0
  sink%v=0d0
  sink%a=0d0
  sink%xpol=0d0
  sink%Jspin=0d0
  sink%jdot=0d0
  sink%jet_dir=[0d0,0d0,1d0]
 end subroutine null_sink

end module derived_types
