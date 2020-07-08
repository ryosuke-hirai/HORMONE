
module merger_mod

 use constants
 
 implicit none
 real*8,parameter:: sep = 20d0*rsun
 real*8,allocatable:: spin_coeffr(:), spin_coefft(:)
 real*8 Mtot, M1, Einject, Jinject, Tinject
 real*8 Iinertia, Jinitial, Einitial
 real*8 domega_dt, de_dt, iniEtot, iniJtot, curEtot, curJtot, inifile, inifile2

 logical:: heatdone
 integer:: iedge
 real*8:: v2old, Erot, edge, din, dout
 real*8,allocatable:: heatV(:), Edistorg(:), Edist(:)

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE MERGER_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: For initial setup of merger simulation.

 subroutine merger_setup

  use grid
  use physval
  use gravmod
  use pressure_mod

  implicit none

!-----------------------------------------------------------------------------

  allocate(Edist(js:je),Edistorg(js:je),heatV(js:je))
  
  if(time==inifile)then
   open(unit=60,file='data/angmom.dat',status='replace')

   call conserve
   call pressure
   call output

   k = ks

   Iinertia = 0d0;Jinitial=0d0;Einitial=0d0;Mtot=0d0
   do j = js, je
    do i = is, ie
     if(x1(i)<0.5d0*sep)then
      Iinertia = Iinertia + d(i,j,k)*(0.4d0*pi*(xi1(i)**5-xi1(i-1)**5)*((cosi(j)**3-cosi(j-1)**3)/3d0-cosi(j)+cosi(j-1)))*2d0
      Jinitial = Jinitial + d(i,j,k)*v3(i,j,k)*x1(i)*sinc(j)*dvol(i,j,k)*2d0
      Einitial = Einitial + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))*dvol(i,j,k)*2d0
      Mtot = Mtot + d(i,j,k)*dvol(i,j,k)*2d0
     end if
    end do
   end do

   iniEtot = 0d0; iniJtot = 0d0 ; Edistorg=0d0
   do j = js, je
    do i = is, ie
     iniEtot = iniEtot + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))*dvol(i,j,k)*2d0
     iniJtot = iniJtot + d(i,j,k)*v3(i,j,k)*x1(i)*sinc(j)*dvol(i,j,k)*2d0
     Edistorg(j) = Edistorg(j) + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))*dvol(i,j,k)*2d0
    end do
   end do

   open(unit=100,file='data/iniEJ.dat',status='replace',form='unformatted')
   write(100)iniEtot,iniJtot,Iinertia,Jinitial,Einitial,Mtot,Edistorg
   close(100)
 
  else
   open(unit=60,file='data/angmom.dat',status='old',position='append')
   open(unit=100,file='data/iniEJ.dat',status='old',form='unformatted')
   read(100)iniEtot,iniJtot,Iinertia,Jinitial,Einitial,Mtot,Edistorg
   close(100)
  end if

  
  M1 = 0.5d0*Mtot
  Einject = -G*M1*(Mtot-M1)/(2d0*sep)
  Jinject = M1*(Mtot-M1)/Mtot*sqrt(G*Mtot*sep)
  Tinject = 2d0*pi*sqrt(sep**3d0/G/Mtot)

  domega_dt = Jinject / Iinertia / Tinject*0d0

 return
 end subroutine merger_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE MERGER
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: For various procedures related to stellar mergers.

subroutine merger

 use grid
 use physval
 use gravmod

 implicit none

!-----------------------------------------------------------------------------

 curEtot = 0d0 ; curJtot = 0d0
 do k = ks, ke
  do j = js, je
   do i = is, ie
    !      if(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)<0d0)then
    curEtot = curEtot + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))&
     *dvol(i,j,k)*2d0
    curJtot = curJtot + d(i,j,k)*v3(i,j,k)*x1(i)*sinc(j)*dvol(i,j,k)*2d0
    !      end if
   end do
  end do
 end do
 if(curJtot-iniJtot>Jinject-Jinitial)then
  domega_dt = 0d0
 end if

 if(domega_dt==0d0.and.(.not.heatdone))then
  heatV = 0d0;Erot=0d0;curEtot = 0d0
  din = 0d0;dout = 0d0
  do i = is, ie
   if(x1(i)>0.5d0*sep.and.din==0d0)then
    din = d(i,je,ks)
   elseif(x1(i)>1d0*sep.and.dout==0d0)then
    dout = d(i,je,ks)
    exit
   end if
  end do
  Edist = 0d0;k=ks
  do j = js, je
   do i = is, ie
    Edist(j) = Edist(j) + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))*dvol(i,j,k)*2d0
   end do
  end do
  do k = ks, ke
   do j = js, je
    do i = is, ie
     !        if(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)<0d0)then
     !        if(x1(i)>sep*0.5d0.and.x1(i)<sep*1d0)then
     if(d(i,j,k)<=din.and.d(i,j,k)>=dout)then
      heatV(j) = heatV(j) + d(i,j,k)*dvol(i,j,k)*2d0          
     end if
     !        Erot=Erot+0.5d0*d(i,j,k)*v2(i,j,k)*v2(i,j,k)*dvol(i,j,k)*2d0
     curEtot = curEtot + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))&
      *dvol(i,j,k)*2d0
    end do
   end do
  end do
  open(unit=1010,file='injectedenergy.dat',status='replace')
  write(1010,'(6(1PE16.8e2))')time,Einject,Einitial,curEtot,iniEtot,Einject-Einitial-curEtot+iniEtot
  do k = ks, ke
   do j = js, je
    do i = is, ie
     !        if(grvphi(i,j,k)*d(i,j,k)+e(i,j,k)<0d0)then
     !        if(x1(i)>sep*0.5d0.and.x1(i)<sep*1d0)then
     if(d(i,j,k)<=din.and.d(i,j,k)>=dout)then
      !         e(i,j,k) = e(i,j,k) + (Einject-Einitial-curEtot+iniEtot)*d(i,j,k)/heatV
      e(i,j,k) = e(i,j,k)+((Einject-Einitial)/dble(je)-Edist(j)+Edistorg(j))&
       * d(i,j,k)/heatV(j)
     end if
    end do
   end do
  end do
  call conserve
  heatdone = .true.
  close(1010)
 end if

return
end subroutine merger

end module merger_mod
