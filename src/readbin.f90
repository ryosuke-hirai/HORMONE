module readbin_mod
 implicit none

 public:: readbin,readgrid,read_extgrv

contains

!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$!
!!$!                     SUBROUTINE ALLOCATE_READGRID
!!$!
!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$
!!$! PURPOSE: To allocate variables and read gridfile
!!$
!!$ subroutine allocate_readgrid
!!$
!!$  use settings
!!$  use grid
!!$  use physval
!!$  use constants
!!$  use gravmod
!!$  use ionization_mod,only:ionization_setup
!!$  use setup_mod,only:read_parameters
!!$
!!$!-----------------------------------------------------------------------------
!!$
!!$! read parameters
!!$  call read_parameters('../parameters')
!!$
!!$! allocate variables
!!$  call checksetup
!!$  call allocations
!!$
!!$! read coordinate data
!!$  call readgrid('gridfile.bin')
!!$
!!$! calculate volume element
!!$ do k = ks, ke
!!$  do j = js-1, je+1
!!$   do i = is-1, ie+1
!!$    dvol(i,j,k)   = (xi1(i)**3d0-xi1(i-1)**3d0) / 3.d0 &
!!$                  * (cos(xi2(j-1))-cos(xi2(j))) * dxi3(k)
!!$!    dvol(i,j,k) = 5.d-1 * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi2(j) * dxi3(k)
!!$!    if(je==1)dvol(i,j,k) = pi * (xi1(i)**2.d0-xi1(i-1)**2.d0) * dxi3(k)
!!$   end do
!!$  end do
!!$ end do
!!$ if(ke==1) dvol = 4d0 * dvol
!!$ T = 1d3
!!$ gamma = 5d0/3d0
!!$ fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
!!$ fac_pgas = kbol/amu ! frequently used factor for Pgas
!!$
!!$ if(eostype==2)then
!!$  call ionization_setup
!!$ end if
!!$
!!$ return
!!$end subroutine allocate_readgrid

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
 ! NOTE: record marker may not always be integer. Could depend on compiler etc
 integer(kind=4) :: record_marker
 logical :: legacy = .true.

!-----------------------------------------------------------------------------

 open(newunit=un,file=filename,status='old',form='unformatted',iostat=istat, access='stream')
 if(istat/=0)then
  print*,'Binary dump file not found'
  print'(3a)','File name = "',trim(filename),'"'
  stop
 end if

 if (legacy) read(un) record_marker
 read(un) tn
 read(un) time
 if (legacy) read(un) record_marker

 if (legacy) read(un) record_marker
 read(un) d (is:ie,js:je,ks:ke)
 read(un) v1(is:ie,js:je,ks:ke)
 read(un) v2(is:ie,js:je,ks:ke)
 read(un) v3(is:ie,js:je,ks:ke)
 read(un) e (is:ie,js:je,ks:ke)
 if (legacy) read(un) record_marker

 if(gravswitch>=2) then
   if (legacy) read(un) record_marker
   read(un) grvphi(gis:gie,gjs:gje,gks:gke)
   if (legacy) read(un) record_marker
 endif

 if(gravswitch==3) then
   if (legacy) read(un) record_marker
   read(un) grvphidot(gis:gie,gjs:gje,gks:gke)
   read(un) dt_old
   if (legacy) read(un) record_marker
 endif

 if(compswitch>=2) then
   if (legacy) read(un) record_marker
   read(un) spc(1:spn,is:ie,js:je,ks:ke)
   read(un) species(1:spn)
   if (legacy) read(un) record_marker
 endif

 if(mag_on)then
  if (legacy) read(un) record_marker
  read(un) b1(is:ie,js:je,ks:ke)
  read(un) b2(is:ie,js:je,ks:ke)
  read(un) b3(is:ie,js:je,ks:ke)
  read(un) phi(is:ie,js:je,ks:ke)
  if (legacy) read(un) record_marker
 end if

 if(include_sinks) then
   if (legacy) read(un) record_marker
   read(un) sink(1:nsink)
   if (legacy) read(un) record_marker
 endif

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
