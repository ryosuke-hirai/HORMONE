module output_mod
 implicit none

 integer:: ievo
 public:: output,set_file_name,write_extgrv,evo_output
 private:: write_grid,write_bin,write_plt,get_header,add_column, &
           write_val

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE OUTPUT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output snapshot data

subroutine output

  use settings,only:is_test
  use grid,only:tn

!----------------------------------------------------------------------------

  if(is_test) return ! do not bother outputting if it is a test

  if(tn==0)call write_grid

  call write_bin
  
  call write_plt

!!$  if(include_particles)then
!!$!ptcfile----------------------------------------------------------------
!!$  if(outstyle==1)then
!!$   write(ptcfile,'(a8,i11.11,a5)')'data/ptc',nint(time),'s.dat'
!!$   if(time>2147483647d0)then
!!$    write(ptcfile,'(a8,i9.9,a7)')'data/ptc',nint(time*0.01d0),'00s.dat'
!!$   end if
!!$  elseif(outstyle==2)then
!!$   write(ptcfile,'(a8,i8.8,a4)')'data/ptc',tn,'.dat'
!!$  end if
!!$
!!$  if(crdnt==1.and.je==1)then
!!$   do n = 1, np
!!$    extflag = .false.
!!$    do k = ks, ke
!!$     if(ptcx(2,n)>=xi3(k-1))then;if(ptcx(2,n)<xi3(k))then
!!$      do i = is, ie
!!$       if(ptcx(1,n)>=xi1(i-1))then;if(ptcx(1,n)<xi1(i))then
!!$        if(e(i,js,k)+grvphi(i,js,k)*d(i,js,k)>0d0)then
!!$         ptci(2,n) = 1
!!$        else
!!$         ptci(2,n) = 0
!!$        end if
!!$        extflag=.true.
!!$        exit
!!$       end if;end if
!!$      end do
!!$      if(extflag)exit
!!$     end if;end if
!!$    end do
!!$   end do
!!$  elseif(crdnt==2.and.ke==1)then
!!$   do n = 1, np
!!$    extflag = .false.
!!$    pr = sqrt(ptcx(1,n)**2d0+ptcx(2,n)**2d0)
!!$    pt = acos(ptcx(2,n)/pr)
!!$    do j = js, je
!!$     if(pt>=xi2(j-1))then;if(pt<xi2(j))then
!!$      do i = is, ie
!!$       if(pr>=xi1(i-1))then;if(pr<xi1(i))then
!!$        if(e(i,j,ks)+grvphi(i,j,ks)*d(i,j,ks)>0d0)then
!!$         ptci(2,n) = 1
!!$        else
!!$         ptci(2,n) = 0
!!$        end if
!!$        extflag=.true.
!!$        exit
!!$       end if;end if
!!$      end do
!!$      if(extflag)exit
!!$     end if;end if
!!$    end do
!!$   end do
!!$  end if
!!$
!!$  open(unit=70,file = ptcfile, status='replace')
!!$
!!$  write(70,'(a,i7,a,1PE12.4e2,2(a5,i9))')&
!!$       '#tn =',tn,'  time= ',time,'np= ',np,'npl=',npl
!!$  write(70,'(a7,2a3,3a15)')'label','ej','ub','mass','x1','x3'
!!$
!!$  do i = 1, np
!!$   write(70,'(i7,2i3,3(1PE15.7e2))')ptci(0:2,i),ptcx(0:2,i)
!!$  end do
!!$
!!$  close(70)
!!$
!!$!bptfile----------------------------------------------------------------
!!$  if(outstyle==1)then
!!$   write(bptfile,'(a8,i11.11,a5)')'data/bpt',nint(time),'s.dat'
!!$   if(time>2147483647d0)then
!!$    write(bptfile,'(a8,i9.9,a7)')'data/bpt',nint(time*0.01d0),'00s.dat'
!!$   end if
!!$  elseif(outstyle==2)then
!!$   write(bptfile,'(a8,i8.8,a4)')'data/bpt',tn,'.dat'
!!$  end if
!!$
!!$  open(unit=80,file = bptfile, status='replace',form='unformatted')
!!$
!!$  write(80)np,npl
!!$  write(80)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)
!!$
!!$  close(80)
!!$
!!$  end if

!remmassfile---------------------------------------------------------------
!!$  if(tn==0)then
!!$   open(unit=60,file='data/remmass.dat',status='replace')
!!$   inimass = 0d0
!!$   do k = ks, ke ; do j = js, je ; do i = is, ie
!!$    if(e(i,j,k)+d(i,j,k)*grvphi(i,j,k)<=0d0)then
!!$     inimass = inimass + dvol(i,j,k)*d(i,j,k) + coremass
!!$    end if
!!$   end do ; end do ; end do
!!$   write(60,'(a,1PE16.8e2,a15,1PE16.8e2)') &
!!$    'initial_env_mass = ',inimass/msun,'core_mass = ',coremass/msun
!!$   write(60,'(7a16)')'time','remmass','belmass','v_kick','v_kickb'
!!$  else
!!$   open(unit=60,file='data/remmass.dat',status='old',position='append')
!!$  end if
!!$
!!$  remmass = coremass ; belmass = coremass
!!$  centre_of_mass = 0d0 ; centre_of_massb = 0d0
!!$  vel = 0d0 ; velb = 0d0
!!$  do k = ks, ke
!!$   do j = js, je
!!$    do i = is, ie
!!$     if(e(i,j,k)+d(i,j,k)*(grvphi(i,j,k)+extgrv(i,j,k))<=0d0)then
!!$      remmass = remmass + d(i,j,k)*dvol(i,j,k)
!!$      vel = vel + d(i,j,k)*dvol(i,j,k)*v3(i,j,k)
!!$     end if
!!$     if(v1(i,j,k)**2d0+v3(i,j,k)**2d0+(grvphi(i,j,k)+extgrv(i,j,k))<=0d0)then
!!$      belmass = belmass + d(i,j,k)*dvol(i,j,k)
!!$      velb = velb + d(i,j,k)*dvol(i,j,k)*v3(i,j,k)
!!$     end if
!!$    end do
!!$   end do
!!$  end do
!!$  write(60,'(7(1PE16.8e2))') time, remmass/msun, belmass/msun, vel/remmass, velb/belmass
!!$  call flush(60)

return
end subroutine output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE OPEN_EVOFILE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Start or open evofile

subroutine open_evofile

 use settings,only:sigfig,gravswitch,mag_on,crdnt
 use grid,only:tn,dim
 
 character(50):: forma

!-----------------------------------------------------------------------------

 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 if(tn==0)then
  open(newunit=ievo,file='data/evo.dat',status='replace')
  write(ievo,'(a10)',advance='no')'tn'
  write(ievo,forma,advance="no")'time'
  write(ievo,forma,advance="no")'tot_mass'
  write(ievo,forma,advance="no")'tot_e'
  write(ievo,forma,advance="no")'tot_eint'
  write(ievo,forma,advance="no")'tot_ekin'
  if(gravswitch>=1)write(ievo,forma,advance="no")'tot_egrv'
  if(mag_on)write(ievo,forma,advance="no")'tot_emag'
  if(dim>=2.and.crdnt>=1)write(ievo,forma,advance="no")'tot_angmom'
  write(ievo,'()')
 else
  open(newunit=ievo,file='data/evo.dat',status='old',position='append')
 end if

return
end subroutine open_evofile

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE EVO_OUTPUT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output quantities as a function of time

subroutine evo_output

 use settings
 use grid
 use physval
 use gravmod

 implicit none

 character*50:: forme
 real*8:: Mtot, Etot, Eitot, Ektot, Egtot, Ebtot, Jtot

!-----------------------------------------------------------------------------
 
 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers

 Mtot = sum(d(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Etot = sum(e(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Eitot = sum(eint(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))
 Ektot = 0.5d0*sum(d(is:ie,js:je,ks:ke)&
                   * ( v1(is:ie,js:je,ks:ke)**2 &
                     + v2(is:ie,js:je,ks:ke)**2 &
                     + v3(is:ie,js:je,ks:ke)**2 ) &
                   * dvol(is:ie,js:je,ks:ke) )
 if(mag_on)then
  Ebtot = 0.5d0*sum( ( b1(is:ie,js:je,ks:ke)**2 &
                     + b2(is:ie,js:je,ks:ke)**2 &
                     + b3(is:ie,js:je,ks:ke)**2 ) &
                   * dvol(is:ie,js:je,ks:ke) )
 end if
 if(gravswitch>=1)then
  if(include_extgrv)then
! grvphi is self-gravity so is double counted, whereas extgrv is not
   Egtot = sum(d(is:ie,js:je,ks:ke) &
              *(0.5d0*grvphi(is:ie,js:je,ks:ke)+extgrv(is:ie,js:je,ks:ke)) &
              *dvol(is:ie,js:je,ks:ke))
  else
   Egtot = 0.5d0*sum(d(is:ie,js:je,ks:ke)*grvphi(is:ie,js:je,ks:ke) &
                    *dvol(is:ie,js:je,ks:ke))
  end if
 end if

 if(dim>=2.and.crdnt>=1)then
  Jtot = 0d0
  select case(crdnt)
  case(1)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      Jtot = Jtot + d(i,j,k)*x1(i)*v2(i,j,k)*dvol(i,j,k)
     end do
    end do
   end do
  case(2)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      Jtot = Jtot + d(i,j,k)*x1(i)*sinc(j)*v3(i,j,k)*dvol(i,j,k)
     end do
    end do
   end do
  end select
 end if

 write(ievo,'(i10)',advance='no')tn
 call write_anyval(ievo,forme,time)
 call write_anyval(ievo,forme,Mtot)
 call write_anyval(ievo,forme,Etot)
 call write_anyval(ievo,forme,Eitot)
 call write_anyval(ievo,forme,Ektot)
 if(gravswitch>=1)call write_anyval(ievo,forme,Egtot)
 if(mag_on)call write_anyval(ievo,forme,Ebtot)
 if(dim>=2.and.crdnt>=1)call write_anyval(ievo,forme,Jtot)
 write(ievo,'()')
 flush(ievo)
 
return
end subroutine evo_output

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE WRITE_GRID
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output gridfile.bin

subroutine write_grid

 use settings
 use grid

 implicit none

 character*50:: formhead,formval,formnum
 integer:: ui

!-----------------------------------------------------------------------------

!binary gridfile--------------------------------------------------------------
 open(newunit=ui,file='data/gridfile.bin',status='replace',form='unformatted')

 write(ui)x1(gis-2:gie+2),xi1(gis-2:gie+2),dx1(gis-2:gie+2),dxi1(gis-2:gie+2),&
          x2(gjs-2:gje+2),xi2(gjs-2:gje+2),dx2(gjs-2:gje+2),dxi2(gjs-2:gje+2),&
          x3(gks-2:gke+2),xi3(gks-2:gke+2),dx3(gks-2:gke+2),dxi3(gks-2:gke+2)
 close(ui)

!gridfile---------------------------------------------------------------------

 open(newunit=ui,file='data/gridfile.dat',status='replace')
 write(ui,'()')
 write(formhead,'("(",i1,"a5,",i1,"a",i2,")")')dim,dim+1,sigfig+8
 write(formval ,'("(",i1,"i5,",i1,"(1x,1PE",i2,".",i2,"e2))")')&
     dim,dim+1,sigfig+7,sigfig-1

 select case (dim)
! 1D outputs
 case(1)
  write(formnum,'("(",a4,"i4,2i",i2,")")')'"#",',sigfig+8
  write(ui,formnum)1,2,3
  
  if(ie>1)then
   write(ui,formhead)'  i','x1','dvol'
   j=js;k=ks
   do i = is, ie
    write(ui,formval)i,x1(i),dvol(i,j,k)    
   end do
  elseif(je>1)then
   write(ui,formhead)'  j','x2','dvol'
   i=is;k=ks
   do j = js, je
    write(ui,formval)j,x2(j),dvol(i,j,k)    
   end do
  elseif(ke>1)then
   write(ui,formhead)'  k','x3','dvol'
   i=is;j=js
   do k = ks, ke
    write(ui,formval)k,x3(k),dvol(i,j,k)    
   end do
  end if

! 2D outputs
 case(2)
  write(formnum,'("(",a4,"i4,i5,3i",i2,")")')'"#",',sigfig+8
  write(ui,formnum)1,2,3,4,5
  
  if(ke==1)then! For 2D Cartesian, polar coordinates or axisymmetrical spherical
   write(ui,formhead)'  i','j','x1','x2','dvol'
   k=ks
! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=js-1
    do i = is, ie
     write(ui,formval)i,j,x1(i),xi2s,dvol(i,j,k)
    end do
    write(ui,'()')
   end if

   do j = js, je
    do i = is, ie
     write(ui,formval)i,j,x1(i),x2(j),dvol(i,j,k)
    end do
    write(ui,'()')
   end do

! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=je+1
    do i = is, ie
     write(ui,formval)i,j,x1(i),xi2e,dvol(i,j,k)
    end do
    write(ui,'()')
   end if
   
  elseif(je==1)then! mainly for 2D Cartesian or axisymmetrical cylindrical
   write(ui,formhead)'  i','k','x1','x3','dvol'
   j=js
   do k = ks, ke, outres
! writing inner boundary for polar coordinates
    if(crdnt==1.or.crdnt==2)&
     write(ui,formval)is-1,k,xi1(is-1),x3(k),dvol(is,j,k)
    do i = is, ie, outres
     write(ui,formval)i,k,xi1(i),x3(k),dvol(i,j,k)
    end do
    write(ui,'()')
   end do   
   
  elseif(ie==1)then! For 2D Cartesian
!CAUTION: Not designed for cylindrical or spherical yet
   write(ui,formhead)'  j','k','x2','x3','dvol'
   i=is
   do k = ks, ke
    do j = js, je
     write(ui,formval)j,k,x2(j),x3(k),dvol(i,j,k)
    end do
    write(ui,'()')
   end do
  end if

 case(3)
  write(formnum,'("(",a4,"i4,2i5,4i",i2,")")')'"#",',sigfig+8
  write(ui,formnum)1,2,3,4,5,6,7
  
  write(ui,formhead)'  i','j','k','x1','x2','x3','dvol'
  if(crdnt==2)then
   do j = je, je
    do i = is, ie
     write(ui,formval)i,j,k,x1(i),x2(j),xi3(ks-1),dvol(i,j,k)
    end do
    write(ui,'()')
   end do
  end if
  
  do k = ks, ke
   do j = je, je
    do i = is, ie
     write(ui,formval)i,j,k,x1(i),x2(j),x3(k),dvol(i,j,k)
    end do
    write(ui,'()')
   end do
  end do

  if(crdnt==2)then
   do j = je, je
    do i = is, ie
     write(ui,formval)i,j,k,x1(i),x2(j),xi3(ke),dvol(i,j,k)
    end do
    write(ui,'()')
   end do
  end if
  
 case default
  stop 'Something wrong with dimension'
 end select

 close(ui)
 
return
end subroutine write_grid


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE WRITE_BIN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output binary full dump file

subroutine write_bin

 use settings
 use grid,only:is,ie,js,je,ks,ke,gis,gie,gjs,gje,gks,gke,time,tn
 use physval
 use gravmod,only:grvphi,grvphiold,dt_old

 implicit none

 character*50:: binfile
 integer:: un
 
!-----------------------------------------------------------------------------

 call set_file_name('bin',tn,time,binfile)
 open(newunit=un,file=binfile,status='replace',form='unformatted')

 write(un)tn,time
 write(un) d (is:ie,js:je,ks:ke), &
           v1(is:ie,js:je,ks:ke), &
           v2(is:ie,js:je,ks:ke), &
           v3(is:ie,js:je,ks:ke), &
           e (is:ie,js:je,ks:ke)
 if(gravswitch>=2)write(un)grvphi(gis:gie,gjs:gje,gks:gke)
 if(gravswitch==3)write(un)grvphiold(gis:gie,gjs:gje,gks:gke),dt_old
 if(compswitch>=2)write(un)spc(1:spn,is:ie,js:je,ks:ke),species(1:spn)
 if(mag_on)then
  write(un) b1(is:ie,js:je,ks:ke), &
            b2(is:ie,js:je,ks:ke), &
            b3(is:ie,js:je,ks:ke), &
            phi(is:ie,js:je,ks:ke)
 end if

 close(un)
 
return
end subroutine write_bin

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE WRITE_PLT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output ascii file for plotting

subroutine write_plt

 use settings
 use grid,only:n,i,j,k,is,ie,js,je,ks,ke,time,tn,dim
 use utils,only:gravpot1d
 use shockfind_mod,only:shockfind

 implicit none

 character*50:: pltfile
 character*20:: header(50)='aaa',forma,forme,formi
 integer:: ui,columns

!-----------------------------------------------------------------------------

! Calculate gravitational potential here if spherical
 if(gravswitch==1) call gravpot1d
! Calculate shock position if required
 if(write_shock)call shockfind
 
! Set format
 write(forma,'("(a",i2,")")')sigfig+8 ! for strings
 write(forme,'("(1x,1PE",i2,".",i2,"e2)")')sigfig+7,sigfig-1 ! for real numbers
 write(formi,'("(i",i2,")")')sigfig+8 ! for integers
! ui = 30
 
! Open file
 call set_file_name('plt',tn,time,pltfile)
 open(newunit=ui,file = pltfile, status='replace')

! Write time and time step
 write(ui,'(a,i7,a,1PE12.4e2,a)')&
  '#tn =',tn,'  time= ',time/dt_unit_in_sec,dt_unit

! Decide what quantities to write
 call get_header(header,columns)
 do n = 1, columns
  write(ui,formi,advance='no') n+2*dim+1
 end do
 write(ui,'()')
 do n = 1, columns
  write(ui,forma,advance="no") trim(adjustl(header(n)))
 end do
 write(ui,'()')

! Write file
 select case (dim)
! 1D outputs
 case(1)
 
  if(ie>1)then
   j=js;k=ks
   do i = is, ie
    call write_val(ui,i,j,k,forme,header)
   end do
  elseif(je>1)then
   i=is;k=ks
   do j = js, je
    call write_val(ui,i,j,k,forme,header)
   end do
  elseif(ke>1)then
   i=is;j=js
   do k = ks, ke
    call write_val(ui,i,j,k,forme,header)
   end do
  end if

! 2D outputs
 case(2)
  
  if(ke==1)then! For 2D Cartesian, polar coordinates or axisymmetrical spherical
   k=ks
! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=js
    do i = is, ie
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end if

   do j = js, je
    do i = is, ie
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end do

! output coordinate axis if cylindrical or spherical coordinates
   if(crdnt==1.or.crdnt==2)then
    j=je
    do i = is, ie
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end if
   
  elseif(je==1)then! mainly for 2D Cartesian or axisymmetrical cylindrical
   j=js
   do k = ks, ke, outres
! writing inner boundary for polar coordinates
    if(crdnt==1.or.crdnt==2)call write_val(ui,is,j,k,forme,header)
    do i = is, ie, outres
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end do   
   
  elseif(ie==1)then! For 2D Cartesian
!CAUTION: Not designed for cylindrical or spherical yet
   i=is
   do k = ks, ke
    do j = js, je
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end do
  end if

 case(3)
  if(crdnt==2)then
   do j = je, je
    do i = is, ie
     call write_val(ui,i,j,ks,forme,header)
    end do
    write(ui,'()')
   end do
  end if
  
  do k = ks, ke
   do j = je, je
    do i = is, ie
     call write_val(ui,i,j,k,forme,header)
    end do
    write(ui,'()')
   end do
  end do

  if(crdnt==2)then
   do j = je, je
    do i = is, ie
     call write_val(ui,i,j,ke,forme,header)
    end do
    write(ui,'()')
   end do
  end if
  
 case default
  stop 'Something wrong with dimension'
 end select


 close(ui)
 

return
end subroutine write_plt

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE WRITE_EXTGRV
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write extgrv file

subroutine write_extgrv

 use grid
 use gravmod,only:extgrv,mc

 implicit none

 integer:: ui

!-----------------------------------------------------------------------------

 open(newunit=ui,file='data/extgrv.bin',status='replace',form='unformatted')
 write(ui)mc(is-1)
 write(ui)extgrv(gis-2:gie+2,gjs:gje,gks-2:gke+2)
 close(ui)

return
end subroutine write_extgrv

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                       SUBROUTINE SET_FILE_NAME
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To create a file name with a given prefix

subroutine set_file_name(prefix,tn,time,filename)

 use settings,only:dt_unit_in_sec,dt_unit,outstyle

 implicit none

 character(*),intent(in):: prefix
 character(*),intent(out)::filename
 integer,intent(in):: tn
 real*8,intent(in):: time
 
!-----------------------------------------------------------------------------
 select case(outstyle)
 case(1)
  write(filename,'(a5,a,i11.11,a,a4)')&
   'data/',trim(prefix),nint(time/dt_unit_in_sec),trim(dt_unit),'.dat'
  if(time/dt_unit_in_sec>2147483647d0)then ! if integer is overflowed
   write(filename,'(a5,a,i9.9,"00",a,a4)')&
    'data/',trim(prefix),nint(time/dt_unit_in_sec*0.01d0),trim(dt_unit),'.dat'
  end if
 case(2)
  write(filename,'(a5,a,i8.8,a4)')'data/',trim(prefix),tn,'.dat'   
 end select
 
return
end subroutine set_file_name


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE GET_HEADER
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To determine what to output

subroutine get_header(header,columns)

 use settings
 use physval,only:species
 use grid,only:dim,ie,je,ke

 implicit none

 character(*),intent(inout):: header(50)
 integer,intent(out):: columns
 integer n

!-----------------------------------------------------------------------------
! Output density/internal energy/pressure by default
 columns = 0
 call add_column('d',columns,header)
 call add_column('e',columns,header)
 call add_column('p',columns,header)  
 
! Decide which velocity components to output
 select case(dim)
 case(1)
  if(ie>1)then
   call add_column('v1',columns,header)
  elseif(je>1)then
   call add_column('v2',columns,header)
  elseif(ke>1)then
   call add_column('v3',columns,header)
  end if
 case(2)
  if(ke==1)then
   call add_column('v1',columns,header)
   call add_column('v2',columns,header)
   if(write_other_vel)call add_column('v3',columns,header)
  elseif(je==1)then
   call add_column('v1',columns,header)
   if(write_other_vel)call add_column('v2',columns,header)
   call add_column('v3',columns,header)
  elseif(ie==1)then
   if(write_other_vel)call add_column('v1',columns,header)
   call add_column('v2',columns,header)
   call add_column('v3',columns,header)
  end if
 case(3)
  call add_column('v1',columns,header)
  call add_column('v2',columns,header)
  call add_column('v3',columns,header)
 case default
  stop 'Something wrong with dimensions'
 end select

! Output magnetic field if mag_on
 if(mag_on)then
  call add_column('b1',columns,header)
  call add_column('b2',columns,header)
  call add_column('b3',columns,header)
 end if

! Output temperature if EoS is non-ideal gas
 if(eostype>=1)then
  call add_column('T',columns,header)
 end if

! Output gravitational potential if gravswitch>=1
 if(gravswitch>=1)then
  call add_column('phi',columns,header)
  if(include_extgrv)then
   call add_column('extphi',columns,header)
  end if
 end if

! Output mean molecular weight if compswitch>=1
 if(compswitch>=1)then
  call add_column('mu',columns,header)
  if(compswitch>=2)then
   do n = 1, spn
    call add_column(species(n),columns,header)
   end do
  end if
 end if

! Output shock position if write_shock=.true.
 if(write_shock)call add_column('shock',columns,header)

return
end subroutine get_header

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE ADD_COLUMN
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To add a column of output variable into the header array

pure subroutine add_column(string,columns,header)

 implicit none

 character(*),intent(inout):: header(50)
 character(*),intent(in):: string
 integer,intent(inout):: columns

!-----------------------------------------------------------------------------

 columns = columns + 1
 header(columns) = string

return
end subroutine add_column

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE WRITE_VAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write a single row of values

subroutine write_val(ui,i,j,k,forme,header)

 use settings,only:spn
 use physval
 use gravmod,only:grvphi,extgrv
 
 implicit none

 integer,intent(in):: ui,i,j,k
 character(*),intent(in):: forme, header(:)
 integer:: n, nn

!-----------------------------------------------------------------------------

 do n = 1, 50
  select case(header(n))
  case('d')!density
   call write_anyval(ui,forme,d(i,j,k))
  case('e')!internal energy
   call write_anyval(ui,forme,eint(i,j,k))
  case('p')!pressure
   call write_anyval(ui,forme,p(i,j,k))
  case('v1')!velocity 1
   call write_anyval(ui,forme,v1(i,j,k))
  case('v2')!velocity 2
   call write_anyval(ui,forme,v2(i,j,k))
  case('v3')!velocity 3
   call write_anyval(ui,forme,v3(i,j,k))
  case('b1')!magnetic field 1
   call write_anyval(ui,forme,b1(i,j,k))
  case('b2')!magnetic field 2
   call write_anyval(ui,forme,b2(i,j,k))
  case('b3')!magnetic field 3
   call write_anyval(ui,forme,b3(i,j,k))
  case('T')!temperature
   call write_anyval(ui,forme,T(i,j,k))
  case('phi')!gravitational potential
   call write_anyval(ui,forme,grvphi(i,j,k))
  case('extphi')!external gravitational potential
   call write_anyval(ui,forme,extgrv(i,j,k))
  case('mu')!mean molecular weight
   call write_anyval(ui,forme,1d0/imu(i,j,k))
  case('shock')!shock position
   call write_anyval(ui,forme,dble(shock(i,j,k)))
  case('aaa')!end of line
   write(ui,'()')
   exit
  case default!chemical composition
   do nn = 1, spn
    if(header(n)==species(nn))then
     call write_anyval(ui,forme,spc(nn,i,j,k))
    end if
   end do
  end select
 end do

return
end subroutine write_val

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE WRITE_ANYVAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To write one value

subroutine write_anyval(ui,forme,val)

 integer,intent(in):: ui
 character(*),intent(in):: forme
 real*8,intent(in):: val
!-----------------------------------------------------------------------------

 write(ui,forme,advance='no')val

return
end subroutine write_anyval

end module output_mod
