module tools_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE TOOLS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Setting frequently used variables

subroutine tools

 use settings
 use grid
 use physval
 use gravmod
 use constants
 use particle_mod
 use pressure_mod
 use ejectamod
 use ionization_mod,only:ionization_setup
 use cooling_mod,only:cooling_setup
 
 integer:: ll

!-----------------------------------------------------------------------------

! ************************ Trigonometric Function ****************************

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(crdnt==1)then
!   allocate( rdis(-1:in,-1:kn), sincyl(-1:in,-1:kn),coscyl(-1:in,-1:kn) )
  allocate( rdis(-1:gin,gks-2:gkn), &
            sincyl(-1:gin,gks-2:gkn),coscyl(-1:gin,gks-2:gkn) )
  do i = gis-1, gie+2
   do k = gks-2, gke+2
    rdis(i,k) = sqrt( x1(i)*x1(i)+x3(k)*x3(k) )
    sincyl(i,k) = x1(i)/rdis(i,k)
    coscyl(i,k) = x3(k)/rdis(i,k)
   end do
  end do
 end if

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(crdnt==2)then
  allocate( sinc(-1:jn), sini(-1:jn), cosc(-1:jn), cosi(-1:jn) )
  do j = js-2, je+2
   sinc(j)=sin(x2 (j))
   sini(j)=sin(xi2(j))
   cosc(j)=cos(x2 (j))
   cosi(j)=cos(xi2(j))
  end do
 end if

! ************************* Legendre polynomials *****************************

! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(crdnt==1)then
  allocate( Plc(0:llmax,0:gin,gks-2:gkn) )
  do k = gks-2, gke+2
   do i = gis, gie+2
    Plc(0,i,k) = 1d0
    Plc(1,i,k) = coscyl(i,k)
    Plc(2,i,k) = (3d0*coscyl(i,k)*Plc(1,i,k) - Plc(0,i,k)) * 0.5d0

    do ll = 3, llmax
     Plc(ll,i,k) = (dble(2*ll-1)*coscyl(i,k)*Plc(ll-1,i,k) &
                 -  dble(ll-1)              *Plc(ll-2,i,k)) /dble(ll)
    end do
   end do
  end do
 end if

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(crdnt==2)then
  allocate( Pl(0:llmax,0:jn) )
  do j=js-1,je+1
   Pl(0,j) = 1d0
   Pl(1,j) = cosc(j)
   Pl(2,j) = (3d0*cosc(j)*Pl(1,j) - Pl(0,j)) * 0.5d0

   do ll=3,llmax 
    Pl(ll,j)   = (dble(2*ll-1)*cosc(j)*Pl(ll-1,j) &
               -  dble(ll-1)          *Pl(ll-2,j)) /dble(ll)
   end do
  end do
 end if

! reading composition data from datafile ! --------------------------------
!!$ open(unit=40,file='17lateRSG.data',status='old')
!!$ read(40,'()')
!!$ read(40,'()')
!!$ read(40,*) lines, lines
!!$ read(40,'()')
!!$ read(40,'()')
!!$ read(40,'(a)') dumc
!!$
!!$! counting rows
!!$ allocate(dum(500)) ; dum = 'aaa'
!!$ read(dumc,*,end=101) dum
!!$101 do i = 1, 500
!!$  if(dum(i)=='aaa')then
!!$   rows = i - 1
!!$   exit
!!$  end if
!!$ end do
!!$
!!$ allocate(header(1:rows),dat(1:lines,1:rows))
!!$ header(1:rows) = dum(1:rows)
!!$ deallocate(dum)
!!$
!!$ do i = 1, lines
!!$  read(40,*) dat(lines-i+1,1:rows) 
!!$ end do
!!$
!!$ allocate(comp_ej(0:8,1:lines))
!!$ do i = 1, rows
!!$  if(trim(header(i))=='mass') comp_ej(0,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='h1'  ) comp_ej(1,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='he4' ) comp_ej(2,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='he3' ) comp_ej(3,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='c12' ) comp_ej(4,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='n14' ) comp_ej(5,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='o16' ) comp_ej(6,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='ne20') comp_ej(7,1:lines) = dat(1:lines,i)
!!$  if(trim(header(i))=='mg24') comp_ej(8,1:lines) = dat(1:lines,i)
!!$ end do
!!$
!!$! reduce unnecessarily high resolution
!!$ i = 2; compsize = lines
!!$ do
!!$  if(abs(comp_ej(1,i)-comp_ej(1,i+1))<1d-3*comp_ej(1,i).and.&
!!$     abs(comp_ej(2,i)-comp_ej(2,i+1))<1d-3*comp_ej(2,i).and.&
!!$     abs(comp_ej(3,i)-comp_ej(3,i+1))<1d-3*comp_ej(3,i).and.&
!!$     abs(comp_ej(4,i)-comp_ej(4,i+1))<1d-3*comp_ej(4,i).and.&
!!$     abs(comp_ej(5,i)-comp_ej(5,i+1))<1d-3*comp_ej(5,i).and.&
!!$     abs(comp_ej(6,i)-comp_ej(6,i+1))<1d-3*comp_ej(6,i).and.&
!!$     abs(comp_ej(7,i)-comp_ej(7,i+1))<1d-3*comp_ej(7,i).and.&
!!$     abs(comp_ej(8,i)-comp_ej(8,i+1))<1d-3*comp_ej(8,i))then
!!$   comp_ej(0:8,i:compsize-1) = comp_ej(0:8,i+1:compsize)
!!$   compsize = compsize-1
!!$   if(compsize<=i+1)exit
!!$   cycle
!!$  end if
!!$  i = i+1
!!$ end do
!!$ deallocate(dat)
!!$ allocate(dat(0:8,1:compsize))
!!$ dat(0:8,1:compsize) = comp_ej(0:8,1:compsize)
!!$ deallocate(comp_ej)
!!$ allocate(comp_ej(0:8,1:compsize))
!!$ comp_ej(0:8,1:compsize) = dat(0:8,compsize:1:-1)

! set default values for composition
 if(compswitch>=1)then
  species(1) = 'h1'
  species(2) = 'he4'
  species(3) = 'metal'
  spc(1,:,:,:) = 0.7d0
  spc(2,:,:,:) = 0.28d0
  spc(3,:,:,:) = 1d0-spc(1,:,:,:)-spc(3,:,:,:)
 end if

! set external gravitational field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$if(include_extgrv)then
!!$ coremass = 8d0
!!$ do k = ks-2, ke+2
!!$  do i = is-2, ie+2
!!$   extgrv(i,js,k) = -coremass*msun*G/sqrt(x1(i)*x1(i)+x3(k)*x3(k)+(3d0*dx1(is))**2d0)!-1.4d0*msun*G/nsdis(i,j,k)
!!$  end do
!!$ end do
!!$ extgrv(is-1,js:je,ks:ke) = extgrv(is,js:je,ks:ke)
!!$ extgrv(is-2,js:je,ks:ke) = extgrv(is+1,js:je,ks:ke)
!!$end if

!!$ grvtime = time

! EoS parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 gamma = 5.d0/3.d0 ! for EoS
 fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
 fac_pgas = kbol/amu ! frequently used factor for Pgas
 imu(is-2:ie+2,js-2:je+2,ks-2:ke+2) = 1d0/muconst

 if(eostype==2)then
  call ionization_setup
 end if

! Set cooling parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(include_cooling) call cooling_setup

return
end subroutine tools

end module tools_mod
