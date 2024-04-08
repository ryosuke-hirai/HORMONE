module smear_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To smear out quantities for nested-grid-like feature

subroutine smear

 use grid,only:is,js,je,ks,ke,fmr_max,fmr_lvl,dim,crdnt
 use physval
 use composition_mod

 implicit none

 integer:: i,j,k,n,jb,kb

!-----------------------------------------------------------------------------

 if(fmr_max==0.or.dim==1)return

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.


 if(crdnt==2)then

!$omp parallel
  do n = 1, fmr_max
   if(fmr_lvl(n)==0)cycle
   if(n==1)then
    jb=je;kb=ke
   else
    jb=min(2**(fmr_max-n+1),je) ; kb=min(2**(fmr_max-n+1),ke)
   end if
!$omp do private(i,j,k) collapse(3)
   do k = ks, ke, kb
    do j = js, je, jb
     do i = is+sum(fmr_lvl(0:n-1)), is+sum(fmr_lvl(0:n))-1
      call angular_smear(i,j,j+jb-1,k,k+kb-1)
     end do
    end do
   end do
!$omp end do
  end do
!$omp end parallel
  
 end if
 
 return
end subroutine smear


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE ANGULAR_SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out quantities over several cells in the angular direction

 subroutine angular_smear(i,js,je,ks,ke)

  use settings,only:spn,compswitch
  use grid,only:x1,x2,x3,dvol,sinc
  use physval,only:u,spc,v1,v2,v3,eint,icnt,iene,imo1,imo2,imo3,ufnmax
  use utils
  use gravmod,only:totphi,extgrv,gravswitch,include_extgrv

  implicit none

  integer,intent(in)::i,js,je,ks,ke
  integer:: n,j,k
  real(8):: mtot, etot, Jtot, Itot, vol, eave
  real(8),dimension(1:3):: vcar, xcar, momtot, vave
  
!-----------------------------------------------------------------------------

  vol = sum(dvol(i,js:je,ks:ke))
  mtot = sum( u(i,js:je,ks:ke,icnt)*dvol(i,js:je,ks:ke) )

  if(compswitch>=2)then
   do n = 1, spn
    spc(n,i,js:je,ks:ke) = sum( u(i,js:je,ks:ke,icnt)*spc(n,i,js:je,ks:ke) &
                               *dvol(i,js:je,ks:ke) ) &
                         / mtot
   end do
  end if
  
!!$  do n = 1, ufnmax
!!$   u(i,js:je,ks:ke,n) = sum(u(i,js:je,ks:ke,n)*dvol(i,js:je,ks:ke))/vol
!!$  end do
!!$  return
  
  momtot=0d0;etot=0d0
  do j = js, je
   do k = ks, ke
    xcar = polcar([x1(i),x2(j),x3(k)])
    call get_vcar(xcar,x3(k),u(i,j,k,imo1),u(i,j,k,imo2),u(i,j,k,imo3),vcar)
    momtot = momtot + vcar*dvol(i,j,k)! add up momenta
    etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
    if(gravswitch>0)then
     etot = etot + u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)! and gravitational ene
    end if
   end do
  end do
  vave = momtot/mtot ! get average cartesian velocity
  u(i,js:je,ks:ke,icnt) = mtot/vol ! density

  do j = js, je
   do k = ks, ke
    xcar = polcar([x1(i),x2(j),x3(k)])
    call get_vpol(xcar,x3(k),vave,v1(i,j,k),v2(i,j,k),v3(i,j,k))
    u(i,j,k,2) = v1(i,j,k)*u(i,j,k,icnt)
    u(i,j,k,3) = v2(i,j,k)*u(i,j,k,icnt)
    u(i,j,k,4) = v3(i,j,k)*u(i,j,k,icnt)
    if(gravswitch>0)then
     etot = etot - u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
    end if
   end do
  end do

  u(i,js:je,ks:ke,iene) = etot / sum( dvol(i,js:je,ks:ke) )
!  u(i,js:je,ks:ke,8) = sum( u(i,js:je,ks:ke,8)*dvol(i,js:je,ks:ke) )&
!                      / sum( dvol(i,js:je,ks:ke) )

!!$  Jtot = 0d0; Itot=0d0; momtot=0d0; vol = sum(dvol(i,js:je,ks:ke))
!!$  do k = ks, ke
!!$   do j = js, je
!!$    Jtot = Jtot + u(i,j,k,imo3)*x1(i)*sinc(j)*dvol(i,j,k)! add up angular momenta
!!$    Itot = Itot + (x1(i)*sinc(j))**2*dvol(i,j,k)! add up moment of inertia
!!$    momtot(1) = momtot(1) + u(i,j,k,imo1)*dvol(i,j,k)! add up radial momenta
!!$    momtot(2) = momtot(2) + u(i,j,k,imo2)*dvol(i,j,k)! add up polar momenta
!!$    etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
!!$    eave = eave + eint(i,j,k)*dvol(i,j,k)
!!$    if(gravswitch>0)then
!!$     etot = etot + u(i,j,k,icnt)*grvphi(i,j,k)*dvol(i,j,k)! and gravitational ene
!!$    end if
!!$    if(include_extgrv)then
!!$     etot = etot + u(i,j,k,icnt)*extgrv(i,j,k)*dvol(i,j,k)! and external gravity
!!$    end if
!!$   end do
!!$  end do
!!$  mtot = sum( u(i,js:je,ks:ke,icnt)*dvol(i,js:je,ks:ke) )
!!$  u(i,js:je,ks:ke,icnt) = mtot/vol ! density
!!$  Itot = Itot * u(i,js,ks,icnt)
!!$  vave(1:2) = momtot(1:2)/mtot
!!$  vave(3) = Jtot/Itot ! get average angular velocity
!!$  eave = eave / vol
!!$
!!$  do k = ks, ke
!!$   do j = js, je
!!$    v3(i,j,k) = vave(3)*x1(i)*sinc(j)
!!$    u(i,j,k,imo3) = u(i,j,k,icnt) * v3(i,j,k)
!!$    if(gravswitch>0)then
!!$     etot = etot - u(i,j,k,icnt)*grvphi(i,j,k)*dvol(i,j,k)
!!$    end if
!!$    if(include_extgrv)then
!!$     etot = etot - u(i,j,k,icnt)*extgrv(i,j,k)*dvol(i,j,k)
!!$    end if
!!$   end do
!!$  end do
!!$
!!$  vave(1:2) = vave(1:2) / norm2(vave(1:2))
!!$  etot = etot - eave*vol - 0.5d0*Itot*vave(3)**2
!!$  if(etot<0d0)then
!!$   if(abs(etot/(eave*vol))>1d-2)then
!!$    print*,'Error in smear second',i,js,je,ks,ke,sum(u(i,js:je,ks:ke,iene)*dvol(i,js:je,ks:ke)),etot,eave*vol,0.5d0*Itot*vave(3)**2
!!$    stop
!!$   else
!!$    etot = 0d0
!!$   end if
!!$  end if
!!$  vave(1:2) = vave(1:2) * sqrt(2d0/u(i,js,ks,icnt)*etot/vol)
!!$  do k = ks, ke
!!$   do j = js, je
!!$    v1(i,j,k) = vave(1)
!!$    v2(i,j,k) = vave(2)
!!$    u(i,j,k,imo1) = v1(i,j,k)*u(i,j,k,icnt)
!!$    u(i,j,k,imo2) = v2(i,j,k)*u(i,j,k,icnt)
!!$   end do
!!$  end do
!!$
!!$  do k = ks, ke
!!$   do j = js, je
!!$    u(i,j,k,iene) = eave + 0.5d0*(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k))*u(i,j,k,icnt)
!!$   end do
!!$  end do
  
 return
end subroutine angular_smear


end module smear_mod
