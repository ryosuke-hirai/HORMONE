module modify_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE MODIFY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To modify an existing binfile to use as an initial condition

subroutine modify

!-----------------------------------------------------------------------------

 call extend2Dto3D

return
end subroutine modify

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EXTEND2DTO3D
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To extend from 2D to 3D

subroutine extend2Dto3D

 use settings,only:gravswitch,compswitch,mag_on,start,include_extgrv,&
                   dt_unit,spn,dt_out
 use grid,only:is,ie,js,je,ks,ke,gis,gie,gjs,gje,gks,gke,time,t_out,tn
 use physval,only:d,e,v1,v2,v3,b1,b2,b3,phi,spc
 use readbin_mod,only:readbin,read_extgrv
 use pressure_mod,only:pressure
 use composition_mod,only:meanmolweight
 use gravmod,only:grvphi,grvpsi,extgrv,grvtime
 use output_mod,only:write_extgrv

 integer:: holdke,holdgke
 character(len=30)::startfile

!-----------------------------------------------------------------------------

 start = 560000
 write(startfile,'(a,i11.11,a,a)')'data/bin',start,trim(dt_unit),'_old.dat'

 holdke = ke; holdgke = gke
 ke = ks ; gke = gks
 call readbin(startfile)
 if(include_extgrv)call read_extgrv('data/extgrv_old.bin')
 ke = holdke ; gke = holdgke

 d (is:ie,js:je,ks:ke) = spread(d (is:ie,js:je,ks),3,ke-ks+1)
 e (is:ie,js:je,ks:ke) = spread(e (is:ie,js:je,ks),3,ke-ks+1)
 v1(is:ie,js:je,ks:ke) = spread(v1(is:ie,js:je,ks),3,ke-ks+1)
 v2(is:ie,js:je,ks:ke) = spread(v2(is:ie,js:je,ks),3,ke-ks+1)
 v3(is:ie,js:je,ks:ke) = spread(v3(is:ie,js:je,ks),3,ke-ks+1)

 if(gravswitch>=2)then
  grvphi(gis:gie,gjs:gje,gks:gke) &
                            = spread(grvphi(gis:gie,gjs:gje,gks),3,gke-gks+1)
  if(gravswitch==3)then
   grvpsi(gis:gie,gjs:gje,gks:gke) &
                         = spread(grvpsi(gis:gie,gjs:gje,gks),3,gke-gks+1)
  end if
  if(include_extgrv)then
   extgrv(gis-2:gie+2,gjs-2:gje+2,gks-2:gke+2) &
                    = spread(extgrv(gis-2:gie+2,gjs-2:gje+2,gks),3,gke-gks+5)
  end if
 end if

 if(compswitch>=2)then
  spc(1:spn,is:ie,js:je,ks:ke) = spread(spc(1:spn,is:ie,js:je,ks),4,ke-ks+1)
 end if

 if(mag_on)then
  b1 (is:ie,js:je,ks:ke) = spread(b1 (is:ie,js:je,ks),3,ke-ks+1)
  b2 (is:ie,js:je,ks:ke) = spread(b2 (is:ie,js:je,ks),3,ke-ks+1)
  b3 (is:ie,js:je,ks:ke) = spread(b3 (is:ie,js:je,ks),3,ke-ks+1)
  phi(is:ie,js:je,ks:ke) = spread(phi(is:ie,js:je,ks),3,ke-ks+1)
 end if

 call meanmolweight
 call pressure

 tn = 0
 t_out = time + dt_out
 if(gravswitch==3)grvtime = time

 call write_extgrv

return
end subroutine extend2Dto3D

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE RADIFY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To turn a hydro simulation dump into a radiation hydro dump

subroutine radify

 use settings,only:extrasfile,eostype,radswitch,eq_sym
 use constants,only:arad
 use grid,only:is,ie,js,je,ks,ke,time,dvol
 use physval,only:d,p,T,imu,eint,erad,e
 use pressure_mod,only:eos_e
 use readbin_mod,only:readbin
 use input_mod,only:error_extras,error_nml
 use output_mod,only:write_bin,write_ascii

 character(len=100):: infile,outfile
 integer:: i,j,k,nn,istat

!-----------------------------------------------------------------------------

 namelist /rdfycon/ infile,outfile

 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('radify',extrasfile)
 read(nn,NML=rdfycon,iostat=istat)
 if(istat/=0)call error_nml('radify',extrasfile)

 eostype=1
 radswitch=0
 call readbin(infile)

! Set radiation pressure
 do k = ks, ke
  do j = js, je
   do i = is, ie
    eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
    erad(i,j,k) = arad*T(i,j,k)**4
    e(i,j,k) = e(i,j,k) - erad(i,j,k)
   end do
  end do
 end do

 eostype=0
 radswitch=1
 time = 3600d3

!!$! Add thermal bomb
!!$ Eheat = 1d50
!!$ vol   = sum(dvol(is:is+30,js:je,ks:ke))
!!$ if(eq_sym)vol=vol*2d0
!!$ do k = ks, ke
!!$  do j = js, je
!!$   do i = is, is+30
!!$!    eint(i,j,k) = eint(i,j,k) + Eheat/vol
!!$    !    e   (i,j,k) = e   (i,j,k) + Eheat/vol
!!$    erad(i,j,k) = erad(i,j,k) + Eheat/vol
!!$   end do
!!$  end do
!!$ end do

 call write_bin(outfile)
 call write_ascii('plt')

 print*,'File converted to a radiation hydrodynamics dump.'
 print*,'Make sure to update the parameters file to switch on radiation.'
 print*,'e.g.'
 print*,'- eostype=0   in &eos_con'
 print*,'- radswitch=1 in &rad_con'
 stop

 return
end subroutine radify

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE BLOWUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To blow up a star

subroutine blowup

 use settings,only:start,dt_unit,eq_sym,dt_out,gravswitch,extrasfile
 use constants,only:rsun
 use grid
 use physval
 use readbin_mod,only:readbin,read_extgrv
 use pressure_mod,only:eos_p
 use source_mod,only:get_totphi
 use gravmod,only:grvtime,grvphi,totphi
 use input_mod,only:error_extras,error_nml
 use output_mod,only:write_bin,write_ascii

 integer:: i,j,k,iinj,nn,istat
 real(8):: Ebind,Eexp,Rinj,rad,Mheat
 character(len=100):: infile,outfile

!-----------------------------------------------------------------------------

 namelist /blwpcon/ infile,outfile,Eexp,Rinj

 open(newunit=nn,file=extrasfile,status='old',iostat=istat)
 if(istat/=0)call error_extras('blowup',extrasfile)
 read(nn,NML=blwpcon,iostat=istat)
 if(istat/=0)call error_nml('blowup',extrasfile)

 call readbin(infile)
 call get_totphi

 Rinj = Rinj * rsun

 Ebind = 0d0
 do k = ks, ke
  do j = js, je
   do i = is, ie
    Ebind = Ebind + (e(i,j,k)+0.5d0*totphi(i,j,k)*d(i,j,k))*dvol(i,j,k)
   end do
  end do
 end do

 t_out = time + dt_out
 if(gravswitch==3)grvtime = time

 Eexp = Eexp*abs(Ebind)

 do i = is, ie
  if(xi1(i)>=Rinj)then
   Rinj = xi1(i)
   iinj = i
   Mheat = sum(d(is:i,js:je,ks:ke)*dvol(is:i,js:je,ks:ke))
   if(eq_sym)Mheat = 2d0*Mheat
   exit
  end if
 end do

 do k = ks, ke
  do j = js, je
   do i = is, iinj
    erad(i,j,k) = erad(i,j,k) + Eexp/Mheat*d(i,j,k)
    p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k),spc(1,i,j,k),spc(2,i,j,k))
    !    e(i,j,k) = e(i,j,k) + Eexp/Mheat*d(i,j,k)
   end do
  end do
 end do

 call write_bin(outfile)

 print*,'Energy injected into specified dump.'
 print*,'Eexp=',Eexp
 stop

 return
end subroutine blowup

end module modify_mod
