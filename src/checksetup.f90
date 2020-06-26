!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE CHECKSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
subroutine checksetup
  
! purpose: To check whether the setups are consistent.

  use settings
  use grid
  use physval
  use gravmod
  
  implicit none

!-----------------------------------------------------------------------------

! Reading dimension of grid
  if(ie==1.and.je==1.and.ke==1) print *,"Error from gridset"
  if(ie==1.and.je==1.and.ke==1) stop
  if(ie/=1.and.je==1.and.ke==1) dim=1
  if(ie==1.and.je/=1.and.ke==1) dim=1
  if(ie==1.and.je==1.and.ke/=1) dim=1
  if(ie/=1.and.je/=1.and.ke==1) dim=2
  if(ie/=1.and.je==1.and.ke/=1) dim=2
  if(ie==1.and.je/=1.and.ke/=1) dim=2
  if(ie/=1.and.je/=1.and.ke/=1) dim=3
  courant = courant / dble(dim) ! To set a consistent courant number.
  hgcfl   = hgcfl   / dble(dim) ! To set a consistent hgcfl number.

! To warn the periodic boundary condition
  if(bc1is*bc1os==0.and.bc1is+bc1os/=0)then
   print *,"Error from boundary condition 'bc1is,bc1os'"
   print *,"Both should be periodic"
   stop
  elseif(bc2is*bc2os==0.and.bc2is+bc2os/=0)then
   print *,"Error from boundary condition 'bc2is,bc2os'"
   print *,"Both should be periodic"
   stop
  elseif(bc3is*bc3os==0.and.bc3is+bc3os/=0)then
   print *,"Error from boundary condition 'bc3is,bc3os'"
   print *,"Both should be periodic"
   stop
  elseif(bc1iv*bc1ov==0.and.bc1iv+bc1ov/=0)then
   print *,"Error from boundary condition 'bc1iv,bc1ov'"
   print *,"Both should be periodic"
   stop
  elseif(bc2iv*bc2ov==0.and.bc2iv+bc2ov/=0)then
   print *,"Error from boundary condition 'bc2iv,bc2ov'"
   print *,"Both should be periodic"
   stop
  elseif(bc3iv*bc3ov==0.and.bc3iv+bc3ov/=0)then
   print *,"Error from boundary condition 'bc3iv,bc3ov'"
   print *,"Both should be periodic"
   stop
  end if

! Setting boundary conditions for 1D and 2D simulations
  if(is==ie)then
   bc1is=2 ; bc1os=2 ; bc1iv=2 ; bc1ov=2
  end if

  if(js==je)then
   bc2is=3 ; bc2os=3 ; bc2iv=3 ; bc2ov=3
  end if

  if(ks==ke)then
   bc3is=3 ; bc3os=3 ; bc3iv=3 ; bc3ov=3
  end if

! Setting boundary condition for equatorial symmetry
  if(eq_sym.and.crdnt==1)then
   bc3is=1
  elseif(eq_sym.and.crdnt==2)then
   bc3os=1
  end if

! Set flag if any boundary is set to Dirichlet boundary condition
  dirichlet_on = .false.
  if(bc1is==9.or.bc1os==9.or.bc2is==9.or.bc2os==9.or.bc3is==9.or.bc3os==9.or.&
     bc1iv==9.or.bc1ov==9.or.bc2iv==9.or.bc2ov==9.or.bc3iv==9.or.bc3ov==9)then
   dirichlet_on = .true.
  end if
  
! check CFL condition
  if(courant>1.d0)then
   print *,"Error from courant number, courant = ",courant
   print *,"courant number should be <= 1.0"
   stop
  end if

! check gravity
  if(dim==1.and.gravswitch==2)then
   print *,"Self-gravity calculations are not necessary for 1D calculations"
   stop
  end if
  if(crdnt/=2.and.gravswitch==1)then
   print *,"Point source gravity routine only for spherical coordinates"
   stop
  end if

  if(gravswitch/=0.and.imesh/=0)then
!   gis = is ; gie = ie ; gjs = js ; gje = je ; gks = ks ; gke = ke
  elseif(gravswitch/=0.and.crdnt==1)then
   gis = is ; gjs = js ; gje = je
  elseif(gravswitch/=0.and.crdnt==2)then
   gis = is ; gjs = js ; gje = je ; gks = ks ; gke = ke
  end if

  gis = min(gis,is) ; gjs = min(gjs,js) ; gks = min(gks,ks)
  gie = max(gie,ie) ; gje = max(gje,je) ; gke = max(gke,ke)

!  if(gravswitch==3)courant = courant / HGfac
  t_out = dt_out

! Spherical composition only for spherical coordinates
  if(compswitch==1.and.crdnt/=2)then
   print *,"compswitch is not consistent with coordinate", compswitch
   stop
  end if

! Cooling only for ideal gas EoS
  if(include_cooling.and.eostype>0)then
   stop "Cooling is only for optically thin equation of state"
  end if

! Recombination only if composition is provided
  if(compswitch<=1.and.eostype==2)then
   print *,'Provide composition to use recombination module'
  end if
  
return
end subroutine checksetup
