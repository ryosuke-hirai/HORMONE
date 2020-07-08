!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE RESTART
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To restart calculation from arbitrary timestep.
!          from INITIALCONDITION

subroutine restart

 use settings!,only:start,outstyle,gravswitch,dt_out,include_particles
 use grid
 use physval
 use gravmod
 use pressure_mod
 use ejtfilemod,only:inimass
 use particle_mod
 
 implicit none

 character*30 startfile, bptfile

!-----------------------------------------------------------------------------

 if(outstyle==1)then
  write(startfile,'(a8,i11.11,a5)')'data/bin',start,'s.dat'
 elseif(outstyle==2)then
  write(startfile,'(a8,i8.8,a4)')'data/bin',start,'.dat'
 end if

 open(unit=11,file=startfile,status='old',form='unformatted')

 read(11)tn,time,mc(is-1)
 read(11) d (is:ie,js:je,ks:ke), &
          v1(is:ie,js:je,ks:ke), &
          v2(is:ie,js:je,ks:ke), &
          v3(is:ie,js:je,ks:ke), &
          e (is:ie,js:je,ks:ke)
 if(gravswitch>=2)read(11)grvphi(gis:gie,gjs:gje,gks:gke)
 if(gravswitch==3)then
  read(11)grvphiold(gis:gie,gjs:gje,gks:gke), &
          dt_old
 end if
 if(compswitch>=2)read(11)spc(1:spn,is:ie,js:je,ks:ke)
 if(mag_on)then
  read(11) b1(is:ie,js:je,ks:ke), &
           b2(is:ie,js:je,ks:ke), &
           b3(is:ie,js:je,ks:ke), &
           phi(is:ie,js:je,ks:ke)
 end if

 close(11)


 call dirichletbound

!bptfile----------------------------------------------------------------
 if(include_particles)then
  if(outstyle==1)then
   write(bptfile,'(a8,i11.11,a5)')'data/bpt',start,'s.dat'
  elseif(outstyle==2)then
   write(bptfile,'(a8,i11.11,a4)')'data/bpt',start,'.dat'
  end if

  open(unit=81,file = bptfile, status='old',form='unformatted')
  allocate( ptcx(0:dim,1:maxptc), ptci(0:2,1:maxptc) )

  read(81)np,npl
  read(81)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

  close(81)
 end if

 t_out = time + dt_out
 if(gravswitch==3)grvtime = time

!extgrvfile-------------------------------------------------------------
 if(include_extgrv)then
  extgrv = 0d0
  open(unit=9191,file='data/extgrv.bin',status='old',form='unformatted')
  read(9191)coremass
  read(9191)extgrv(gis-2:gie+2,gjs:gje,gks-2:gke+2)
  close(9191)
 end if

end subroutine restart
