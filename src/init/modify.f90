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
!                              SUBROUTINE BLOWUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To blow up a star

subroutine blowup

 use settings,only:start,dt_unit,eq_sym,dt_out,gravswitch
 use grid
 use physval
 use readbin_mod,only:readbin,read_extgrv
 use pressure_mod,only:eos_p
 use gravmod,only:grvtime,grvphi,extgrv

 integer:: i,j,k,iinj
 real(8):: Ebind,Eexp,Rinj,rad,Mheat
 character(len=30)::startfile

!-----------------------------------------------------------------------------

 start = 20000
 write(startfile,'(a,i11.11,a,a)')'data/bin',start,trim(dt_unit),'_old.dat'
 call readbin(startfile)
 call read_extgrv('data/extgrv.bin')

 start = 0

 Ebind = 0d0
 do k = ks, ke
  do j = js, je
   do i = is, ie
    Ebind = Ebind + (e(i,j,k)+(0.5d0*grvphi(i,j,k)+extgrv(i,j,k))*d(i,j,k))*dvol(i,j,k)
   end do
  end do
 end do

 t_out = time + dt_out
 if(gravswitch==3)grvtime = time

 rad = 7.4d13
 Rinj = rad*1d0/15d0
 Eexp = 0.5d0*abs(Ebind)

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
    eint(i,j,k) = eint(i,j,k) + Eexp/Mheat*d(i,j,k)
    p(i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k),spc(1,i,j,k),spc(2,i,j,k))
   end do
  end do
 end do

 return
end subroutine blowup

end module modify_mod
