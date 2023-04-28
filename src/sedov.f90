module sedov_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE SEDOV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Set initial conditions for Sedov-Taylor blast wave

subroutine sedov

 use settings,only:simtype,extrasfile
 use grid
 use physval
 use pressure_mod,only:eos_p
 use input_mod,only:error_extras,error_nml

 real(8):: damb, Eexp, ein, pin, pamb, Tin, imuconst
 integer:: i_inj,ui,istat

!-----------------------------------------------------------------------------

 namelist /sedocon/ gamma,damb,Eexp

 select case(simtype)
 case('sedov_default')
  open(newunit=ui,file='../para/extras_sedov',status='old',iostat=istat)
 case('sedov_other')
  open(newunit=ui,file=extrasfile,status='old',iostat=istat)
  if(istat/=0)call error_extras('sedov',extrasfile)
 end select
 read(ui,NML=sedocon,iostat=istat)
 if(istat/=0)call error_nml('sedov',extrasfile)
 close(ui)

 i_inj = 10
 ein = Eexp/sum(dvol(is:i_inj,js:je,ks:ke))
 Tin = 1d3
 imuconst = 1d0/muconst
 pin = eos_p(damb,ein,Tin,imuconst)
 pamb = 1d-3*Eexp/sum(dvol(is:ie,js:je,ks:ke)) ! Low energy

 d(is:ie,js:je,ks:ke) = damb
 p(is:ie,js:je,ks:ke) = pamb
 p(is:i_inj,js:je,ks:ke) = pin
 v1(is:ie,js:je,ks:ke) = 0d0
 v2(is:ie,js:je,ks:ke) = 0d0
 v3(is:ie,js:je,ks:ke) = 0d0

return
end subroutine sedov

end module sedov_mod
