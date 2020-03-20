!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!                                                                           !
!                                                                           !
!                                                                           !
!                             PROGRAM HORMONE                               !
!                                                                           !
!                                                                           !
!                                                                           !
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!
!                                                                           !
!    High ORder Magnetohydrodynamic cOde with Numerous Enhancements         !
!                                                                           !
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\!

! PURPOSE: Main program to solve ideal MHD problems

!############################################################################

program hormone

  use settings
  use grid
  use physval
  use constants
  use pressure_mod
  use particle_mod
  use merger_mod
  use funcs

  use gravmod

  implicit none

  logical:: heatdone
  integer:: iedge
  real*8:: v2old, Erot, edge, din, dout
  real*8,allocatable:: heatV(:), Edistorg(:), Edist(:)

! set parameters
  namelist /gridcon/ xi1s, xi1e, xi2s, xi2e, xi3s, xi3e, &
                     is, ie, js, je, ks, ke, imesh, jmesh, kmesh, &
                     sphrn, trnsn1, trnsn2, trnsn3
  namelist /timecon/ outstyle, endstyle, tnlim, t_end, dt_out, tn_out
  namelist /eos_con/ eostype, eoserr, compswitch, muconst, spn
  namelist /simucon/ crdnt,courant, rktype, start
  namelist /bouncon/ bc1is, bc1os, bc2is, bc2os, bc3is, bc3os, &
                     bc1iv, bc1ov, bc2iv, bc2ov, bc3iv, bc3ov, eq_sym
  namelist /gravcon/ gravswitch, grverr, cgerr, HGfac, hgcfl, &
                     include_extgrv, gis, gie, gjs, gje, gks, gke
  namelist /partcon/ include_particles, maxptc


!############################## start program ################################

! Initial setups -------------------------------------------------------------

  time = 0.d0; tn = 0
  heatdone = .true.
  inifile = 7d5
  inifile2= 7d5

  ! Reading parameters
  open(unit=1,file='parameters',status='old')
  read(1,NML=gridcon)
  read(1,NML=timecon)
  read(1,NML=eos_con)
  read(1,NML=simucon)
  read(1,NML=bouncon)
  read(1,NML=gravcon)
  read(1,NML=partcon)
  close(1)

  call checksetup
  call allocations
  call gridset
  call tools
  call metric
  call gravsetup

  call initialcondition
  call boundarycondition
  call meanmolweight
  call conserve
  call pressure
  
  if(time==inifile)then
   open(unit=60,file='data/angmom.dat',status='replace')
  else
   open(unit=60,file='data/angmom.dat',status='old',position='append')
  end if
  allocate(Edist(js:je),Edistorg(js:je),heatV(js:je))

! Initial output
  if(bc3os==9)call dirichletbound ! for ejecta
  if(include_particles.and.tn==0)call particles_setup
  if(include_particles.and.time==inifile)call particles_setup
  call boundarycondition
  call timestep
  if(gravswitch==3.and.tn==0)dt_old=dt / (courant*HGfac) * hgcfl

  if(time==inifile)then

   call conserve
   call pressure
   call output

   k = ks

   Iinertia = 0d0;Jinitial=0d0;Einitial=0d0;Mtot=0d0
   do j = js, je
    do i = is, ie
     if(x1(i)<0.5d0*sep)then
      Iinertia = Iinertia + d(i,j,k)*(0.4d0*pi*(pw(5,xi1(i))-pw(5,xi1(i-1)))*((pw(3,cosi(j))-pw(3,cosi(j-1)))/3d0-cosi(j)+cosi(j-1)))*2d0
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
   open(unit=100,file='data/iniEJ.dat',status='old',form='unformatted')
   read(100)iniEtot,iniJtot,Iinertia,Jinitial,Einitial,Mtot,Edistorg
   close(100)
  end if

  M1 = 0.5d0*Mtot
  Einject = -G*M1*(Mtot-M1)/(2d0*sep)
  Jinject = M1*(Mtot-M1)/Mtot*sqrt(G*Mtot*sep)
  Tinject = 2d0*pi*sqrt(sep**3d0/G/Mtot)

!!$  edge = 60d0*rsun
!!$  do i = is, ie
!!$   if(x1(i)>edge)then
!!$    iedge = i
!!$    exit
!!$   end if
!!$  end do
!!$  do while(curJtot<iniJtot-Jinitial+Jinject)
!!$   iedge = iedge - 1
!!$   edge = x1(iedge)
!!$   do j = js, je
!!$    do i = is, ie
!!$     if(x1(i)>edge.and.x1(i)<70d0*rsun.and.-grvphi(i,j,k)>v3(i,j,k)*v3(i,j,k))then
!!$      v3(i,j,k) = sqrt(-grvphi(i,je,k))*sinc(j)
!!$      e(i,j,k) = eint(i,j,k)+0.5d0*d(i,j,k)*(v1(i,j,k)**2d0+v2(i,j,k)**2d0+v3(i,j,k)**2d0)
!!$     end if
!!$    end do
!!$   end do
!!$   curJtot = 0d0;curEtot = 0d0
!!$   do j = js, je
!!$    do i = is, ie
!!$     curEtot = curEtot + (0.5d0*grvphi(i,j,k)*d(i,j,k)+e(i,j,k))&
!!$                        *dvol(i,j,k)*2d0
!!$     curJtot = curJtot + d(i,j,k)*v3(i,j,k)*x1(i)*sinc(j)*dvol(i,j,k)*2d0
!!$    end do
!!$   end do
!!$  end do

  domega_dt = Jinject / Iinertia / Tinject*0d0

! Start integration ----------------------------------------------------------
  if(tnlim/=0)then
   do while(time < t_end.and.tn<=tnlim)

    call timestep

    if(tn>0)call gravity
    if(bc3os==9)call dirichletbound ! for ejecta

    do rungen = 1, rktype
     call boundarycondition
     call numflux
     call source
     call rungekutta
    end do


    if(include_particles)call particles

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

    time = time + dt ; tn = tn + 1
! Output sequence ---------------------- !
    if(outstyle==2)then                  !
     if(tn/=0.and.mod(tn,tn_out)==0)then !
      call boundarycondition             !
      call output                        !
     end if                              !
    elseif(outstyle==1)then              !
     if(time>=t_out)then                 !
      call boundarycondition             !
      call output                        !
      t_out = t_out + dt_out             !
     end if                              !
    else                                 !
     print *,'outstyle out of range'     !
     stop                                !
    end if                               !
! -------------------------------------- !

    if(endstyle==1)then
     if(time>=t_end)exit
    elseif(endstyle==2)then
     if(tn>=tnlim)exit
    end if

   end do
  end if
! End integration ------------------------------------------------------------

  call output ! To see final state

!------------------------------- end program ---------------------------------

  print *, 'Calculation complete! tn = ',tn

end program hormone
