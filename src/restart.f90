!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE RESTART
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To restart calculation from arbitrary timestep.
!          from INITIALCONDITION

subroutine restart

 use settings
 use grid
 use physval
 use gravmod
 use pressure_mod
 use particle_mod
 use dirichlet_mod
 use output_mod,only:set_file_name
 
 implicit none

 real*8:: starttime
 character*30 startfile, bptfile

!-----------------------------------------------------------------------------

 starttime = dble(start)*dt_unit_in_sec
 call set_file_name('bin',start,starttime,startfile)
 
 open(unit=11,file=startfile,status='old',form='unformatted')

 read(11)tn,time
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


 if(dirichlet_on)call dirichletbound

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
