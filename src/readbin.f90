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
 use pressure_mod
 use composition_mod
 use gravmod
 use sink_mod,only:nsink,sink

 character(len=*),intent(in):: filename
 integer:: un,istat

!-----------------------------------------------------------------------------

 open(newunit=un,file=filename,status='old',form='unformatted',iostat=istat)
 if(istat/=0)then
  print*,'Binary dump file not found'
  print'(3a)','File name = "',trim(filename),'"'
  stop
 end if
 read(un) tn,time
 read(un) d (is:ie,js:je,ks:ke), &
          v1(is:ie,js:je,ks:ke), &
          v2(is:ie,js:je,ks:ke), &
          v3(is:ie,js:je,ks:ke), &
          e (is:ie,js:je,ks:ke)
 if(gravswitch>=2)read(un)grvphi(gis:gie,gjs:gje,gks:gke)
 if(gravswitch==3)read(un)grvphidot(gis:gie,gjs:gje,gks:gke),dt_old
 if(compswitch>=2)read(un)spc(1:spn,is:ie,js:je,ks:ke),species(1:spn)
 if(mag_on)then
  read(un) b1(is:ie,js:je,ks:ke), &
           b2(is:ie,js:je,ks:ke), &
           b3(is:ie,js:je,ks:ke), &
           phi(is:ie,js:je,ks:ke)
 end if
 if(include_sinks)read(un)sink(1:nsink)
 close(un)

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

 use grid,only:is,gis,gie,gjs,gje,gks,gke
 use gravmod,only:coremass,extgrv,mc

 character(len=*),intent(in)::filename
 integer:: ui
!-----------------------------------------------------------------------------

 open(newunit=ui,file=filename,status='old',form='unformatted')
 extgrv = 0d0
 read(ui)coremass
 read(ui)extgrv(gis-2:gie+2,gjs-2:gje+2,gks-2:gke+2)
 close(ui)

 mc(is-1) = coremass
 
 return
end subroutine read_extgrv

end module readbin_mod
