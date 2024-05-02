module readbin_mod
 implicit none

 public:: readbin,readgrid,read_extgrv

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE READBIN
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read physval's from a binary file

subroutine readbin(filename)

 use settings
 use grid
 use physval
 use io
 use pressure_mod
 use composition_mod
 use gravmod
 use sink_mod,only:nsink,sink

 character(len=*),intent(in):: filename
 integer:: un

!-----------------------------------------------------------------------------

 call open_file_read(filename, un)

 call read_dummy_recordmarker(un)
 call read_var(un, tn)
 call read_var(un, time)
 call read_dummy_recordmarker(un)

 call read_dummy_recordmarker(un)
 call read_var(un, d, is, ie, js, je, ks, ke)
 call read_var(un, v1, is, ie, js, je, ks, ke)
 call read_var(un, v2, is, ie, js, je, ks, ke)
 call read_var(un, v3, is, ie, js, je, ks, ke)
 call read_var(un, e, is, ie, js, je, ks, ke)
 call read_dummy_recordmarker(un)

 if(gravswitch>=2) then
   call read_dummy_recordmarker(un)
   call read_var(un, grvphi, gis, gie, gjs, gje, gks, gke)
   call read_dummy_recordmarker(un)
 endif

 if(gravswitch==3) then
   call read_dummy_recordmarker(un)
   call read_var(un, grvphidot, gis, gie, gjs, gje, gks, gke)
   call read_var(un, dt_old)
   call read_dummy_recordmarker(un)
 endif

 if(compswitch>=2) then
   call read_dummy_recordmarker(un)
   call read_var(un, spc, 1, spn, is, ie, js, je, ks, ke)
   call read_var(un, species, 1, spn)
   call read_dummy_recordmarker(un)
 endif

 if(mag_on)then
  call read_dummy_recordmarker(un)
  call read_var(un, b1, is, ie, js, je, ks, ke)
  call read_var(un, b2, is, ie, js, je, ks, ke)
  call read_var(un, b3, is, ie, js, je, ks, ke)
  call read_var(un, phi, is, ie, js, je, ks, ke)
  call read_dummy_recordmarker(un)
 end if

 if(include_sinks) then
   call read_dummy_recordmarker(un)
   call read_var(un, sink, 1, nsink)
   call read_dummy_recordmarker(un)
 endif

 call close_file(un)

 call meanmolweight
 call pressure

return
end subroutine readbin

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE READGRID
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read gridfile.bin

subroutine readgrid(filename)

 use grid

 character(len=*),intent(in)::filename
 integer:: ui
!-----------------------------------------------------------------------------

 open(newunit=ui,file=filename,status='old',form='unformatted')
 read(ui)x1(gis-2:gie+2),xi1(gis-2:gie+2),dx1(gis-2:gie+2),dxi1(gis-2:gie+2), &
         x2(gjs-2:gje+2),xi2(gjs-2:gje+2),dx2(gjs-2:gje+2),dxi2(gjs-2:gje+2), &
         x3(gks-2:gke+2),xi3(gks-2:gke+2),dx3(gks-2:gke+2),dxi3(gks-2:gke+2)
 close(ui)

 return
end subroutine readgrid

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE READ_EXTGRV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read gridfile.bin

subroutine read_extgrv(filename)
 use grid, only:is,gis,gie,gjs,gje,gks,gke,is_global
 use gravmod, only:coremass,extgrv,mc
 use io, only:open_file_read,read_dummy_recordmarker,read_var,close_file
 character(len=*),intent(in)::filename
 integer:: ui

 extgrv = 0d0

 call open_file_read(filename, ui)

 call read_dummy_recordmarker(ui)
 call read_var(ui, coremass)
 call read_dummy_recordmarker(ui)

 call read_dummy_recordmarker(ui)
 call read_var(ui, extgrv, gis, gie, gjs, gje, gks, gke, grav=.true.)
 call read_dummy_recordmarker(ui)

 call close_file(ui)

 if (is==is_global) mc(is-1) = coremass

 return
end subroutine read_extgrv

end module readbin_mod
