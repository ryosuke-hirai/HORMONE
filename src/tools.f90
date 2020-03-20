!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                              SUBROUTINE TOOLS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Setting frequently used variables

subroutine tools

 use grid
 use physval
 use gravmod
 use constants
 use particle_mod
 use pressure_mod
 use merger_mod
 
 implicit none

 real*8 coremass

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

! set external gravitational field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(include_extgrv)then
 coremass = 8d0
 do k = ks-2, ke+2
  do i = is-2, ie+2
   extgrv(i,js,k) = -coremass*msun*G/sqrt(x1(i)*x1(i)+x3(k)*x3(k)+(3d0*dx1(is))**2d0)!-1.4d0*msun*G/nsdis(i,j,k)
  end do
 end do
 extgrv(is-1,js:je,ks:ke) = extgrv(is,js:je,ks:ke)
 extgrv(is-2,js:je,ks:ke) = extgrv(is+1,js:je,ks:ke)
end if

 grvtime = time

! EoS parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 gamma = 5.d0/3.d0 ! for EoS
 fac_egas = kbol/((gamma-1d0)*amu) ! frequently used factor for egas
 fac_pgas = kbol/amu ! frequently used factor for Pgas
 imu(is:ie,js:je,ks:ke) = 1d0/muconst

! for merger product %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! domega_dt = Jinject / Iinertia / Tinject
! de_dt = 5d50 /1d2/msun/ 3d10 *0d0

return
end subroutine tools
