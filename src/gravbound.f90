

module gravbound_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVBOUND
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for miccg with companion star

 subroutine gravbound

  use settings,only:bc3is,crdnt,grverr
  use grid,only:is,ie,js,je,ks,ke,gis,gie,gjs,gje,gks,gke,xi1s,x1,tn,rdis
  use constants,only:G
  use gravmod,only:llmax,grvphi,phiii,phiio,phi1o,phi3i,phi3o,mc,Pl,Plc

  integer:: i,j,k,ll, error1
  real(8):: dphiii, dphiio, dphi1o, dphi3i, dphi3o
  real(8),dimension(0:llmax):: ml
  logical,dimension(0:llmax):: got

!------------------------------------------------------------------------------

! cylindrical (axial symmetry) ################################################
  if(crdnt==1.and.je==1.and.bc3is/=1)then

   got = .false.

   error1=0
   do i = gie+1,gie+2
    do k=gks,gke
     phi1o(i,k) = 0d0
     do ll=0,llmax
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphi1o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**ll/rdis(i,k)
      if(ll/=0.and.abs(dphi1o) < grverr*abs(phi1o(i,k))) exit
      phi1o(i,k) = phi1o(i,k) + dphi1o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi1o(i,k) = G * phi1o(i,k)! + pripot1o(i,k)
    end do !k-loop

    if(error1==1)then
     write(6,*)"Error from gravbound 1: Number of terms in multipole expansion &
         &is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do k = gks, gke
    do i = gie+1, gie+2
     grvphi(i,js,k) = phi1o(i,k)
    end do
   end do

   error1=0
   do k = gks-2, gks-1
    do i = gis, gie
     phi3i(i,k) = 0d0
     do ll=0,llmax
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphi3i = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**ll/rdis(i,k)
     if(i==gis)print'(4i5,5(1PE15.7e2))',tn,i,k,ll,phi3i(i,k),dphi3i,(rdis(ie,ke)/rdis(i,k))**ll,ml(ll),Plc(ll,i,k)
      if(ll/=0.and.abs(dphi3i) < grverr*abs(phi3i(i,k))) exit
      phi3i(i,k) = phi3i(i,k) + dphi3i
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3i(i,k) = G*phi3i(i,k)! + pripot3i(i)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3i: Number of terms in multipole expansion &
        &  is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   error1=0
   do k = gke+1,gke+2
    do i=gis,gie
     phi3o(i,k) = 0d0
     do ll=0,llmax
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphi3o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**ll/rdis(i,k)
      if(ll/=0.and.abs(dphi3o) < grverr*abs(phi3o(i,k))) exit
      phi3o(i,k) = phi3o(i,k) + dphi3o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3o(i,k) = G*phi3o(i,k)! + pripot3o(i)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3o: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = gis, gie
    do k = gke+1, gke+2
     grvphi(i,js,k) = phi3o(i,k)     
    end do
    do k = gks-2, gks-1
     grvphi(i,js,k) = phi3i(i,k)
    end do
   end do

! cylindrical (equatorial+axial symmetry) #####################################
  elseif(crdnt==1.and.je==1.and.bc3is==1)then

   got = .false.

   error1=0
   do i = gie+1, gie+2
    do k=gks,gke
     phi1o(i,k) = 0d0
     do ll=0,llmax,2
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphi1o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**ll/rdis(i,k)
      if(ll/=0.and.abs(dphi1o) < grverr*abs(phi1o(i,k))) exit
      phi1o(i,k) = phi1o(i,k) + dphi1o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi1o(i,k) = 2d0*G* phi1o(i,k)! + pripot1o(i,k)
    end do !k-loop

    if(error1==1)then
     write(6,*)"Error from gravbound 1: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = ie+1, ie+2
    do k = ks, ke
     grvphi(i,js,k) = phi1o(i,k)
    end do
   end do

   error1=0
   do k = gke+1, gke+2
    do i=gis,gie
     phi3o(i,k) = 0d0
     do ll=0,llmax,2
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphi3o = -ml(ll)*Plc(ll,i,k)*(rdis(ie,ke)/rdis(i,k))**ll/rdis(i,k)
      if(ll/=0.and.abs(dphi3o) < grverr*abs(phi3o(i,k))) exit
      phi3o(i,k) = phi3o(i,k) + dphi3o
     end do !ll-loop
     if(ll>=llmax) error1=1
     phi3o(i,k) = 2d0*G*phi3o(i,k)! + pripot3o(i,k)
    end do !i-loop

    if(error1==1)then
     write(6,*)"Error from gravbound3o: Number of terms in multipole expansion &
         & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do i = gis, gie
    do k = gke+1, gke+2
     grvphi(i,js,k) = phi3o(i,k)
    end do
   end do
   grvphi(is:ie,js,ks-1) = grvphi(is:ie,js,ks)
   grvphi(is:ie,js,ks-2) = grvphi(is:ie,js,ks+1)
   phi3i = 0d0
   grvphi(is-1,js,ks:ke) = grvphi(is,js,ks:ke)
   grvphi(is-2,js,ks:ke) = grvphi(is+1,js,ks:ke)

! spherical coordinates (axial symmetry) #######################################
  elseif(crdnt==2.and.ke==1)then

   got = .false.

   error1=0
   do i = gie+1, gie+2
    do j = gjs,gje
     phiio(i,j) = 0d0
     do ll=0,llmax
      if(.not.got(ll))call multipole(ll,ml(ll),got(ll))
      dphiio = -G*Pl(ll,j)*ml(ll)/x1(i)
      if(ll/=0.and.abs(dphiio) < grverr*abs(phiio(i,j))) exit
      phiio(i,j) = phiio(i,j) + dphiio
     end do !ll-loop
     phiio(i,j) = phiio(i,j)! - G*mc(is-1)/x1(i)
     if(ll>=llmax) error1=1
    end do

    if(error1==1)then
     write(6,*)"Error from gravbound i: Number of terms in multipole expansion &
        & is not enough. Computation stopped. tn=",tn,"i=",i
     stop
    end if
   end do

   do j = gjs, gje
    do i = gie+1, gie+2
     grvphi(i,j,ks) = phiio(i,j)
    end do
   end do

   if(xi1s>0d0)then
    error1=0
    do i = gis-2, gis-1
     do j = gjs, gje
      phiii(i,j) = 0d0
      do ll=0,llmax
       call multipoleinner(ll,ml(ll))
       dphiii = -G*Pl(ll,j)*ml(ll)/x1(is-1)
       if(ll/=0.and.abs(dphiii) < grverr*abs(phiii(i,j))) exit
       phiii(i,j) = phiii(i,j) + dphiii
      end do !ll-loop
      phiii(i,j) = phiii(i,j) - G*mc(is-1)/x1(is-1)
      if(ll>=llmax) error1=1
     end do

     if(error1==1)then
      write(6,*)"Error from gravbound i: Number of terms in inner multipole &
           & expansion is not enough. Computation stopped. tn=",tn,"i=",i
      stop
     end if
    end do
   end if


   if(mc(is-1)==0.d0.and.xi1s/=0.d0)error1=1
   if(error1==1)then
    write(6,*)"Error from gravbound i: inner boundary is wrong",xi1s,mc(is-1)
    stop
   end if

  end if

  return
 end subroutine gravbound

!------------------------------------------------------------------------------
 subroutine multipole(ll,ml,got)

  use settings,only:eq_sym
  use grid,only:is,ie,js,je,ks,ke,x1,dvol,rdis,crdnt
  use gravmod,only:Pl,Plc,gsrc

  integer,intent(in):: ll
  real(8),intent(out):: ml
  logical,intent(out):: got
  integer:: i, j, k

  ml = 0d0
! Cylindrical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==1.and.je==1)then
   j = js
!$omp parallel do private (i,k) collapse(2) reduction(+:ml)
   do k = ks,ke
    do i = is,ie
     ml = ml + gsrc(i,j,k)*(rdis(i,k)/rdis(ie,ke))**ll &
             * Plc(ll,i,k)*dvol(i,j,k)
     if(eq_sym)then
      ml = ml + gsrc(i,j,k)*(-rdis(i,k)/rdis(ie,ke))**ll &
              * Plc(ll,i,k)*dvol(i,j,k)
     end if
    end do
   end do
!$omp end parallel do

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elseif(crdnt==2.and.ke==1)then
   k = ks
!$omp parallel do private (i,j) collapse(2) reduction(+:ml)
   do j = js,je
    do i = is,ie
     ml = ml + gsrc(i,j,k)*(x1(i)/x1(ie+1))**ll &
             * Pl(ll,j) * dvol(i,j,k)
     if(eq_sym)then
      ml = ml + gsrc(i,j,k)*(-x1(i)/x1(ie+1))**ll &
              * Pl(ll,j) * dvol(i,j,k)
     end if
    end do
   end do
!$omp end parallel do

  end if

  got = .true.

  return
  end subroutine multipole


!------------------------------------------------------------------------------
 subroutine multipoleinner(ll,ml)

   use grid,only:crdnt,is,ie,js,je,ks,ke,x1,dvol
   use gravmod,only:Pl,gsrc

   integer,intent(in):: ll
   real(8),intent(out):: ml
   integer:: i, j, k

!------------------------------------------------------------------------------

! Spherical >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(crdnt==2.and.ke==1)then
   k = ks
   ml = 0d0
!$omp parallel do private(i,j) collapse(2) reduction(+:ml)
   do j = js,je
    do i = is,ie
     ml = ml + gsrc(i,j,k)*(x1(is-1)/x1(i))**(ll+1) &
             * Pl(ll,j) * dvol(i,j,k)
    end do
   end do
!$omp end parallel do
  end if

  return
  end subroutine multipoleinner

end module gravbound_mod
