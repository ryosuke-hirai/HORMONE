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

 use grid,only:is_global,js_global,je_global,ks_global,ke_global,&
               fmr_max,fmr_lvl,dim,crdnt
 use physval
 use composition_mod
 use profiler_mod

 implicit none

 integer:: i,j,k,n,jb,kb

!-----------------------------------------------------------------------------

 if(fmr_max==0.or.dim==1)return

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.

 call start_clock(wtsmr)

 if(crdnt==2)then

  do n = 1, fmr_max
   if(fmr_lvl(n)==0)cycle
   if(n==1)then
    jb=je_global-js_global; kb=ke_global-ks_global
   else
    jb=min(2**(fmr_max-n+1),je_global)-1 ; kb=min(2**(fmr_max-n+1),ke_global)-1
   end if
   do k = ks_global, ke_global, kb+1
    do j = js_global, je_global, jb+1
     do i = is_global+sum(fmr_lvl(0:n-1)), is_global+sum(fmr_lvl(0:n))-1
      call angular_smear(i,j,j+jb,k,k+kb)
     end do
    end do
   end do
  end do

 end if

 call stop_clock(wtsmr)

 return
end subroutine smear


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE ANGULAR_SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out quantities over several cells in the angular direction

 subroutine angular_smear(i,js_,je_,ks_,ke_)

  use settings,only:spn,compswitch
  use grid,only:x1,x2,x3,dvol,is,ie,js,je,ks,ke,car_x
  use physval,only:u,spc,v1,v2,v3,icnt,iene,imo1,imo2,imo3
  use utils
  use gravmod,only:totphi,gravswitch
  use mpi_utils,only:allreduce_mpi
  use mpi_domain,only:sum_global_array,partially_my_domain

  implicit none

  integer,intent(in)::i,js_,je_,ks_,ke_
  integer:: n,j,k
  integer:: jl,jr,kl,kr
  real(8):: mtot, etot, spctot, vol
  real(8),dimension(1:3)::  vave
  logical:: overlap
  real(8),dimension(1:3):: momtot, compen, tempsum, element, vcar, momtotlocal, compenglobal

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)
  overlap = partially_my_domain(i,js_,ks_,0,je_-js_+1,ke_-ks_+1)

  vol = sum(dvol(i,js_:je_,ks_:ke_))
  mtot = sum_global_array(u,i,i,js_,je_,ks_,ke_,icnt,weight=dvol)

  if(compswitch>=2)then
   do n = 1, spn
    spctot = sum_global_array(u,i,i,js_,je_,ks_,ke_,icnt, &
                              l_weight2=n, weight=dvol, weight2=spc )
    if (overlap) spc(n,i,jl:jr,kl:kr) = spctot / mtot
   end do
  end if

  momtot=0d0;etot=0d0;compen=0d0;momtotlocal=0d0;compenglobal=0d0
  if (overlap) then
!!$omp parallel firstprivate(compen,momtotlocal)
!!$omp do private(j,k,vcar,tempsum,element) &
!!$omp reduction(+:etot) collapse(2)
   do k = kl, kr
    do j = jl, jr
     call get_vcar(car_x(:,i,j,k),x3(k),&
                   u(i,j,k,imo1),u(i,j,k,imo2),u(i,j,k,imo3),vcar)
     element = real(vcar*dvol(i,j,k),kind=8)-compen
     tempsum = momtotlocal + element
     compen = get_compensation(element,tempsum,momtotlocal)
     momtotlocal = tempsum ! add up momenta using the Kahan method
!!$     momtot = momtot + real(vcar*dvol(i,j,k),kind=8)
     etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
     if(gravswitch>0)& ! and gravitational energy
      etot = etot + u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
    end do
   end do
!!$omp end do
!!$omp critical
   tempsum = momtot
   element = momtotlocal - compenglobal
   momtot = tempsum + element
   compenglobal = get_compensation(element,momtot,tempsum)
!!$omp end critical
!!$omp end parallel
  end if

  call allreduce_mpi('sum',momtot)
  vave = real(momtot,kind=8)/mtot ! get average cartesian velocity

  if (overlap) u(i,jl:jr,kl:kr,icnt) = mtot/vol ! density

  if (overlap) then
!!$omp parallel do private(j,k) collapse(2) reduction(+:etot)
   do k = kl, kr
    do j = jl, jr
     call get_vpol(car_x(:,i,j,k),x3(k),vave,v1(i,j,k),v2(i,j,k),v3(i,j,k))
     u(i,j,k,2) = v1(i,j,k)*u(i,j,k,icnt)
     u(i,j,k,3) = v2(i,j,k)*u(i,j,k,icnt)
     u(i,j,k,4) = v3(i,j,k)*u(i,j,k,icnt)
     if(gravswitch>0)&
      etot = etot - u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
    end do
   end do
!!$omp end parallel do
  endif
  call allreduce_mpi('sum',etot)
  if (overlap) u(i,jl:jr,kl:kr,iene) = etot / vol
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

pure elemental function get_compensation(y,t,s) result(c)
! Get compensation term for the Kahan-Babuska-Neumaier method
 real(8),intent(in):: y,t,s
 real(8):: c

 if(abs(s)>=abs(y))then
  c = (s-t)+y
 else
  c = (y-t)+s
 end if

end function get_compensation

end module smear_mod
