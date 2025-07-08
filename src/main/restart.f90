module restart_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE RESTART
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
 use readbin_mod,only:readbin,read_extgrv

 real(8):: starttime
 character(len=30):: startfile, bptfile
 integer:: ui

!-----------------------------------------------------------------------------

 starttime = dble(start)*dt_unit_in_sec
 call set_file_name('bin',start,starttime,startfile)

 call readbin(startfile)


 if(dirichlet_on)call dirichletbound

!bptfile----------------------------------------------------------------
 if(include_particles)then
  if(outstyle==1)then
   write(bptfile,'(a8,i11.11,a5)')'data/bpt',start,'s.dat'
  elseif(outstyle==2)then
   write(bptfile,'(a8,i11.11,a4)')'data/bpt',start,'.dat'
  end if

  open(newunit=ui,file = bptfile, status='old',form='unformatted')
  allocate( ptcx(0:dim,1:maxptc), ptci(0:2,1:maxptc) )

  read(ui)np,npl
  read(ui)ptci(0:2,1:np),ptcx(0:2,1:np),ptc_in(1:jmax)

  close(ui)
 end if

 t_out = time + dt_out
 if(gravswitch==3)grvtime = time

!extgrvfile-------------------------------------------------------------
 if(include_extgrv)call read_extgrv('data/extgrv.bin')

end subroutine restart

end module restart_mod
