module analysis_mod
 implicit none

 real(8),allocatable,dimension(:,:):: dat,comp
 real(8),allocatable,dimension(:):: m, r, rho, pres
 integer lines
 
 contains
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LUMINOSITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate luminosity assuming black body radiation

subroutine get_luminosity(angle,lum)

 use grid
 use physval
 use constants
 use utils
 use utils_analysis

 implicit none

 real(8),intent(in):: angle
 real(8),intent(out):: lum
 real(8):: tau, tau0, drmin, dtheta, dA, rho, XX, TT, dzabs, lumt, dzabs0
 real(8),dimension(1:3):: dz, x, xp, dz0
 real(8),parameter:: tauerr=1d-3
 real(8),allocatable:: dr(:), rr(:)
 integer:: ii, jj, kk, iin, jjn, kkn

!-----------------------------------------------------------------------------

 tau0 = 1d0 ! optical depth at photosphere
 iin = 300 ! radial resolution
 jjn = 10 ! angular resolution
 kkn = 100 ! depth resolution
 allocate(dr(-1:iin+2),rr(-1:iin+2))
 drmin = dxi1(is)*0.2d0
 call geometrical_series(dr,drmin,1,iin,0d0,xi1(ie))
 rr(0) = 0d0
 do ii = 1, iin
  rr(ii) = rr(ii-1)+dr(ii)
 end do
 dtheta = 2d0*pi/dble(jjn)
 dz0 = (/0d0,0d0,-xi1(ie)/dble(kkn)/)
 dz0 = rotx(dz0,angle)
 dzabs0 = norm2(dz0)
 lum = 0d0
 lumt = 0d0

!$omp parallel do default(none) reduction(+:lumt) &
!$omp private(ii,jj,kk,dA,x,xp,tau,rho,XX,TT,dz,dzabs) &
!$omp shared(dtheta,tau0,angle,iin,jjn,kkn,xi1,ie,dz0,dzabs0,rr)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   tau = 0d0;dz=dz0;dzabs=dzabs0
   x(1) = 0.5d0*(rr(ii-1)+rr(ii))*cos(jj*dtheta)
   x(2) = 0.5d0*(rr(ii-1)+rr(ii))*sin(jj*dtheta)
   x(3) = sqrt(xi1(ie)**2-x(1)**2-x(2)**2)
   dA = 0.5d0*(rr(ii)**2-rr(ii-1)**2)*dtheta
   x = rotx(x,angle)
   xp = carpol(x)
! find photosphere
   do while (dot_product(x,dz0)<xi1(ie)*dzabs0)
    x = x + dz
    xp = carpol(x)
    call get_local_val(xp,rho,XX,TT)
    tau = tau + rho*kap_es(XX)*dzabs

    if(tau-tau0>tauerr*tau0)then ! tan ran over tau0
     tau = tau - rho*kap_es(XX)*dzabs
     x = x - dz
     xp = carpol(x)
     dz = dz*0.1d0
     dzabs = norm2(dz)!dzabs*0.1d0
     cycle
    elseif(tau>tau0.and.abs(tau-tau0)<tauerr*tau0)then ! tau=tau0
! add black body flux
     call get_local_val(xp,rho,XX,TT)
     lumt = lumt+dA*sigma*TT**4
!     if(TT>0d0)print'(2i5,3(1PE13.5e2))',ii,jj,x(3),TT,tau
     exit
    end if
   end do

  end do
 end do
!$omp end parallel do

 lum = lumt

 return
end subroutine get_luminosity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LUMINOSITY2
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate luminosity assuming black body radiation

subroutine get_luminosity2(angle,lum)

 use grid
 use physval
 use constants,only:pi,sigma
 use utils,only:carpol,rotx,geometrical_series
 use utils_analysis

 implicit none

 real(8),intent(in):: angle
 real(8),intent(out):: lum
 real(8):: tau, drmin, dtheta, dA, rho, XX, TT, dzabs, lumt, dtau
 real(8),dimension(1:3):: x, xp, dz0
 real(8),allocatable:: dr(:), rr(:)
 integer:: ii, jj, iin, jjn
 real(8):: source, intens

!-----------------------------------------------------------------------------

 iin = 600 ! radial resolution
 jjn = 100 ! angular resolution
 allocate(dr(-1:iin+2),rr(-1:iin+2))
 drmin = dxi1(is)*0.2d0
 call geometrical_series(dr,drmin,1,iin,0d0,xi1(ie))
 rr(0) = 0d0
 do ii = 1, iin
  rr(ii) = rr(ii-1)+dr(ii)
 end do
 dtheta = 2d0*pi/dble(jjn)
 dz0 = (/0d0,0d0,1d0/)
 dz0 = rotx(dz0,angle)
 lum = 0d0
 lumt = 0d0

!$omp parallel do reduction(+:lumt) collapse(2)  &
!$omp private(ii,jj,dA,x,xp,tau,rho,XX,TT,dzabs,dtau,source,intens)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   x(1) = 0.5d0*(rr(ii-1)+rr(ii))*cos(jj*dtheta)
   x(2) = 0.5d0*(rr(ii-1)+rr(ii))*sin(jj*dtheta)
   x(3) = -sqrt(xi1(ie)**2-x(1)**2-x(2)**2)*0.97d0
   dA = 0.5d0*(rr(ii)**2-rr(ii-1)**2)*dtheta
   x = rotx(x,angle)

! Shoot ray
   intens=0d0;tau=0d0;dzabs=0d0;dtau=1d99
   ray_loop:do while (dot_product(x,x)<xi1(ie)**2)

    x = x + dz0*dzabs
    xp = carpol(x)

    if(xp(1)<norm2(x))then
     print*,'error'
     print*,ii,jj,xp(1),norm2(x)
     stop
    end if

    call get_local_val(xp,rho,XX,TT,dzabs)
    dtau = rho*kap(XX,TT)*dzabs
    tau = tau + dtau
    source = sigma*TT**4/pi
    intens = (intens-source)*exp(-dtau)+source

   end do ray_loop
   lumt = lumt+dA*intens
  end do
 end do
!$omp end parallel do

 lum = lumt*4d0*pi

 return
end subroutine get_luminosity2


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LUMINOSITY3
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate luminosity assuming black body radiation
! v3:: Compute multiple frequencies so that we can get an effective temperature

subroutine get_luminosity3(angle,lmda,lum,spec)

 use settings,only:eq_sym
 use grid
 use physval
 use constants,only:pi,sigma
 use utils,only:carpol,roty,rotz,geometrical_series
 use utils_analysis

 implicit none

 real(8),intent(in):: angle(1:2)
 real(8),allocatable,intent(in):: lmda(:)
 real(8),intent(out):: lum
 real(8),allocatable,intent(inout):: spec(:)
 real(8):: drmin, dtheta, dA, rho, XX, TT, dzabs, lumt
 real(8),dimension(1:3):: x, xp, dz0
 real(8),allocatable,dimension(:):: dr, rr, int_nu, S_nu, dtau, tau
 integer:: ii, jj, iin, jjn
 real(8):: source, intens, dtau_grey
 logical:: half

!-----------------------------------------------------------------------------

 iin = 600 ! radial resolution
 jjn = 100 ! angular resolution
 allocate(dr(-1:iin+2),rr(-1:iin+2))
 drmin = dxi1(is)
 call geometrical_series(dr,drmin,1,iin,0d0,xi1(ie))
 rr(0) = 0d0
 do ii = 1, iin
  rr(ii) = rr(ii-1)+dr(ii)
 end do
 dtheta = 2d0*pi/dble(jjn)
 dz0 = [0d0,0d0,1d0]
 dz0 = rotz(roty(dz0,angle(1)),angle(2))
 lum = 0d0
 lumt = 0d0
 spec = 0d0
 half = .false.
 if(eq_sym.and.abs(angle(1)*2d0/pi-1d0)<1d-10)then
  half=.true.
  jjn=jjn/2
 end if

 allocate(int_nu,S_nu,tau,dtau,mold=lmda)

!$omp parallel do reduction(+:lumt,spec) collapse(2)  &
!$omp private(ii,jj,dA,x,xp,tau,rho,XX,TT,dzabs,dtau,source,intens,int_nu,S_nu,dtau_grey)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   x(1) =-0.5d0*(rr(ii-1)+rr(ii))*sin(jj*dtheta)
   x(2) = 0.5d0*(rr(ii-1)+rr(ii))*cos(jj*dtheta)
   x(3) = -sqrt(xi1(ie)**2-x(1)**2-x(2)**2)*0.97d0
   dA = 0.5d0*(rr(ii)**2-rr(ii-1)**2)*dtheta
   x = rotz(roty(x,angle(1)),angle(2))

! Shoot ray
   dzabs=0d0
   intens=0d0
   tau=0d0
   dtau=1d99
   int_nu=0d0
   dtau_grey=1d99
   ray_loop:do while (dot_product(x,x)<xi1(ie)**2)

    x = x + dz0*dzabs
    xp = carpol(x)

    if(xp(1)<norm2(x))then
     ! for some reason the Intel compiler causes this to happen
     print*,'error'
     print*,ii,jj,xp(1),norm2(x)
     stop
    end if

    call get_local_val(xp,rho,XX,TT,dzabs)
! Assuming grey opacities
    dtau = rho*kap(XX,TT,lmda)*dzabs
    dtau_grey = rho*kap(XX,TT)*dzabs
    tau = tau + dtau
    source = sigma*TT**4/pi
    intens = (intens-source)*exp(-dtau_grey)+source
    S_nu = planck_lambda(lmda,TT)
    int_nu = (int_nu-S_nu)*exp(-dtau)+S_nu

   end do ray_loop
   lumt = lumt+dA*intens
   spec = spec+dA*int_nu
  end do
 end do
!$omp end parallel do

 lum = lumt*4d0*pi
 if(half)then
  lum = lum*2d0
  spec = spec*2d0
 end if

 return
end subroutine get_luminosity3

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LUMINOSITY4
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate luminosity assuming black body radiation
! v4:: Faster integration by shooting rays from the observer and truncating when
!      it gets optically thick

subroutine get_luminosity4(angle,lmda,lum,spec)

 use settings,only:eq_sym
 use grid
 use physval
 use constants,only:pi,sigma
 use utils,only:carpol,roty,rotz,geometrical_series
 use utils_analysis

 implicit none

 real(8),intent(in):: angle(1:2)
 real(8),allocatable,intent(in):: lmda(:)
 real(8),intent(out):: lum
 real(8),allocatable,intent(inout):: spec(:)
 real(8):: drmin, dtheta, dA, rho, XX, TT, dzabs, lumt
 real(8),dimension(1:3):: x, xp, dz0
 real(8),allocatable,dimension(:):: dr, rr, int_nu, S_nu, dtau, tau, dint_nu
 integer:: ii, jj, iin, jjn
 real(8):: source, intens, dtau_grey, tau_grey, dintens
 logical:: half

!-----------------------------------------------------------------------------

 iin = 600 ! radial resolution
 jjn = 100 ! angular resolution
 allocate(dr(-1:iin+2),rr(-1:iin+2))
 drmin = dxi1(is)
 call geometrical_series(dr,drmin,1,iin,0d0,xi1(ie))
 rr(0) = 0d0
 do ii = 1, iin
  rr(ii) = rr(ii-1)+dr(ii)
 end do
 dtheta = 2d0*pi/dble(jjn)
 dz0 = [0d0,0d0,1d0]
 dz0 = rotz(roty(dz0,angle(1)),angle(2))
 lum = 0d0
 lumt = 0d0
 spec = 0d0
 half = .false.
 if(eq_sym.and.abs(angle(1)*2d0/pi-1d0)<1d-10)then
  half=.true.
  jjn=jjn/2
 end if

 allocate(int_nu,S_nu,tau,dtau,mold=lmda)

!$omp parallel do reduction(+:lumt,spec) collapse(2)  &
!$omp private(ii,jj,dA,x,xp,rho,XX,TT,dzabs,dtau,tau,source,intens,dintens,&
!$omp         int_nu,dint_nu,S_nu,dtau_grey,tau_grey)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   x(1) =-0.5d0*(rr(ii-1)+rr(ii))*sin(jj*dtheta)
   x(2) = 0.5d0*(rr(ii-1)+rr(ii))*cos(jj*dtheta)
   x(3) = sqrt(xi1(ie)**2-x(1)**2-x(2)**2)*0.97d0
   dA = 0.5d0*(rr(ii)**2-rr(ii-1)**2)*dtheta
   x = rotz(roty(x,angle(1)),angle(2))

! Shoot ray
   dzabs=0d0
   intens=0d0
   tau=0d0
   dtau=1d99
   int_nu=0d0
   tau_grey=0d0
   dtau_grey=1d99
   ray_loop:do while (dot_product(x,x)<xi1(ie)**2)

    x = x - dz0*dzabs
    xp = carpol(x)

    if(xp(1)<norm2(x))then
     ! for some reason the Intel compiler causes this to happen
     print*,'error'
     print*,ii,jj,xp(1),norm2(x)
     stop
    end if

    call get_local_val(xp,rho,XX,TT,dzabs)
! Assuming grey opacities
    dtau = rho*kap(XX,TT,lmda)*dzabs
    dtau_grey = rho*kap(XX,TT)*dzabs
    source = sigma*TT**4/pi
    dintens = source*(1d0-exp(-dtau_grey))*exp(-tau_grey)
    intens = intens + dintens
    S_nu = planck_lambda(lmda,TT)
    dint_nu = S_nu*(1d0-exp(-dtau))*exp(-tau)
    int_nu = int_nu + dint_nu
    if(maxval(abs(dint_nu/int_nu))<1d-8)exit ray_loop
    tau = tau + dtau
    tau_grey = tau_grey + dtau_grey
   end do ray_loop

   lumt = lumt+dA*intens
   spec = spec+dA*int_nu
  end do
 end do
!$omp end parallel do

 lum = lumt*4d0*pi
 if(half)then
  lum = lum*2d0
  spec = spec*2d0
 end if

 return
end subroutine get_luminosity4


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GET_LOCAL_VAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get local values of density, temperature and hydrogen fraction

 subroutine get_local_val(xp,rho,XX,TT,dz)

  use grid
  use settings
  use physval
  use constants,only:pi
  use utils,only:intpol
  
  implicit none

  real(8),intent(in):: xp(1:3)
  real(8),intent(out):: rho
  real(8),intent(out),optional:: XX, TT
  real(8),intent(inout),optional:: dz
  real(8):: xp1, xp2, xp3, rhot(1:2), XXt(1:2), TTt(1:2)
  real(8):: uu,vv,ww
  integer:: i,j,k,ii,jj,kk

!-----------------------------------------------------------------------------

! For 2D spherical coordinates
  if(crdnt==2.and..not.solve_k)then
   if(xp(1)>=x1(ie))then
    rho = 0d0
    if(present(XX))XX  = 0.7d0
    if(present(TT))TT  = 1d-30
    return
   end if
   kk = ks
   xp1 = xp(1)
   xp2 = min(xp(2),pi-xp(2))
   do jj = js, je
    if(x2(jj)>=xp2)exit
   end do
   if(jj==js)xp2=x2(js)
   do ii = is, ie
    if(x1(ii)>=xp1)then
     if(ii==is)xp1 = x1(ii)
     rhot(1) = intpol(x1(ii-1:ii),d(ii-1:ii,jj-1,kk),xp1)
     rhot(2) = intpol(x1(ii-1:ii),d(ii-1:ii,jj  ,kk),xp1)
     rho = intpol(x2(jj-1:jj),rhot(1:2),xp2)
     if(present(XX))then
      XXt(1) = intpol(x1(ii-1:ii),spc(1,ii-1:ii,jj-1,kk),xp1)
      XXt(2) = intpol(x1(ii-1:ii),spc(1,ii-1:ii,jj  ,kk),xp1)
      XX = intpol(x2(jj-1:jj),XXt(1:2),xp2)
     end if
     if(present(TT))then
      TTt(1) = intpol(x1(ii-1:ii),T(ii-1:ii,jj-1,kk),xp1)
      TTt(2) = intpol(x1(ii-1:ii),T(ii-1:ii,jj  ,kk),xp1)
      TT = intpol(x2(jj-1:jj),TTt(1:2),xp2)
     end if
     exit
    end if
   end do
  elseif(crdnt==2.and.solve_i.and.solve_j.and.solve_k)then
! For 3D spherical coordinates
   if(xp(1)>=x1(ie))then ! Ignore when outside computational domain
    rho = 0d0
    if(present(XX))XX  = 0.7d0
    if(present(TT))TT  = 1d-30
    if(present(dz))dz  = dxi1(ie)
    return
   end if

   xp1 = xp(1)
   xp2 = xp(2)
   if(eq_sym)xp2 = min(xp(2),pi-xp(2))
   xp3 = xp(3)

   do ii = is, ie
    if(x1(ii)>=xp1)exit
   end do
   do jj = js, je
    if(x2(jj)>=xp2)exit
   end do
   if(jj==js)xp2=x2(js)
   do kk = ks, ke
    if(x3(kk)>=xp3)exit
   end do

   uu = (xp1-x1(ii-1))/(x1(ii)-x1(ii-1))
   vv = (xp2-x2(jj-1))/(x2(jj)-x2(jj-1))
   ww = (xp3-x3(kk-1))/(x3(kk)-x3(kk-1))

   rho = 0d0
   XX = 0.7d0
   TT = 0d0
   do k = 0, 1
    do j = 0, 1
     do i = 0, 1
      rho = rho + d(ii+i-1,jj+j-1,kk+k-1)*(1d0-uu)**(1-i)*uu**i &
                                         *(1d0-vv)**(1-j)*vv**j &
                                         *(1d0-ww)**(1-k)*ww**k
      if(present(XX).and.spn>0)then
       XX = XX + spc(1,ii+i-1,jj+j-1,kk+k-1)*(1d0-uu)**(1-i)*uu**i &
                                            *(1d0-vv)**(1-j)*vv**j &
                                            *(1d0-ww)**(1-k)*ww**k
      end if
      if(present(TT).and.eostype>0)then
       TT = TT + T(ii+i-1,jj+j-1,kk+k-1)*(1d0-uu)**(1-i)*uu**i &
                                        *(1d0-vv)**(1-j)*vv**j &
                                        *(1d0-ww)**(1-k)*ww**k
      end if
     end do
    end do
   end do

   if(present(dz))dz=0.3d0*min(dxi1(ii),x1(ii)*dxi2(jj),x1(ii)*dxi3(kk))

  end if


 return
 end subroutine get_local_val

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GET_LOCAL_VAL3D
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get local values of whatever quantity is input

!!$ subroutine get_local_val3D(xp,data,val)
!!$
!!$  use grid
!!$  use constants,only:pi
!!$  use utils,only:intpol
!!$  
!!$  implicit none
!!$
!!$  real(8),intent(in):: xp(1:3)
!!$  real(8),allocatable,intent(in):: data(:,:,:)
!!$  real(8),intent(out):: val
!!$  integer:: ii,jj,kk
!!$
!!$!-----------------------------------------------------------------------------
!!$  
!!$  if(ke>1)then
!!$   if(xp(1)>=x1(ie))then
!!$    val = 0d0
!!$    return
!!$   end if
!!$   phi_loop: do kk = ks, ke
!!$    if(xp(3)>xi3(kk-1).and.xp(3)<=xi3(kk))then
!!$     do jj = js, je
!!$      if(xp(2)>xi2(jj-1).and.xp(2)<=xi2(jj))then
!!$       do ii = is, ie
!!$        if(xp(1)>xi1(ii-1).and.xp(1)<=xi1(ii))then
!!$         val = data(ii,jj,kk)
!!$         exit phi_loop
!!$        end if
!!$       end do
!!$      end if
!!$     end do
!!$    end if
!!$   end do phi_loop
!!$
!!$  end if
!!$
!!$ return
!!$end subroutine get_local_val3D


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                           SUBROUTINE GET_COM
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get centre of mass position and velocity

 subroutine get_com(xcom,vcom,spin,totmass)

  use grid
  use physval
  use utils
  
  implicit none

  integer:: i,j,k
  real(8),dimension(1:3),intent(out):: xcom, vcom, spin
  real(8),intent(out):: totmass
  real(8),dimension(1:3):: xcar, xpol, vcar, angm
  real(8):: dm
!-----------------------------------------------------------------------------

  xcom = 0d0
  totmass = 0d0
  vcom = 0d0
  angm = 0d0
!$omp parallel do default(none) &
!$omp private(i,j,k,xpol,xcar,vcar,dm)&
!$omp reduction (+:totmass,xcom,vcom,angm)&
!$omp shared (v1,v2,v3,d,dvol,x1,x2,x3,ks,ke,js,je,is,ie)
  do k = ks, ke
   xpol(3) = x3(k)
   do j = js, je
    xpol(2) = x2(j)
    do i = is, ie
     xpol(1) = x1(i)
     xcar = polcar(xpol)
     call get_vcar(xcar,x3(k),v1(i,j,k),v2(i,j,k),v3(i,j,k),vcar)
     dm = d(i,j,k)*dvol(i,j,k)
     totmass = totmass + dm
     xcom = xcom + dm*xcar
     vcom = vcom + dm*vcar
     angm = angm + dm*cross(xcar,vcar)
    end do
   end do
  end do
!$omp end parallel do

  xcom = xcom/totmass
  xcom(3) = 0d0
  vcom = vcom/totmass
  vcom(3) = 0d0
  totmass = totmass*2d0

  angm(1:2) = 0d0
  angm(3) = 2d0*angm(3)
  spin = angm - totmass*cross(xcom,vcom)

 return
 end subroutine get_com

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GET_ACCRETION
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculated accreted amount of mass onto the NS

!!$ subroutine get_accretion(xcom,vcom,accmass,accJ,mstar,Jspin,Mej)
!!$
!!$  use grid
!!$  use physval
!!$  use constants,only:G,pi
!!$  use gravmod,only:nsmass,nspos,nsvel,nssoft,grvphi
!!$  use utils
!!$
!!$  implicit none
!!$
!!$  real(8),dimension(1:3),intent(in):: xcom,vcom
!!$  real(8),intent(out):: accmass,accJ,mstar,Mej,Jspin
!!$  real(8):: nspot
!!$  real(8),dimension(1:3):: angm, xcar,xpol,vcar,spin
!!$  logical:: boundtoNS,boundtostar
!!$
!!$!-----------------------------------------------------------------------------
!!$
!!$  accmass = 0d0
!!$  angm = 0d0
!!$  mstar = 0d0
!!$  Mej = 0d0
!!$!$omp parallel do private(i,j,k,xpol,xcar,vcar,nspot,boundtoNS,boundtostar) &
!!$!$omp reduction (+:accmass,angm,mstar,Mej)
!!$  do k = ks, ke
!!$   xpol(3) = x3(k)
!!$   do j = js, je
!!$    xpol(2) = x2(j)
!!$    do i = is, ie
!!$     boundtoNS=.false.;boundtostar=.false.
!!$     xpol(1) = x1(i)
!!$     xcar = polcar(xpol)
!!$     call get_vcar(xcar,x3(k),v1(i,j,k),v2(i,j,k),v3(i,j,k),vcar)
!!$     nspot = G*nsmass*softened_potential(norm2(xcar-nspos),nssoft)
!!$     if(eint(i,j,k)/d(i,j,k)+0.5d0*dot_product(vcar-nsvel,vcar-nsvel)+nspot<=0d0)then ! bound to NS
!!$      boundtoNS=.true.
!!$!      accmass = accmass + d(i,j,k)*dvol(i,j,k)
!!$!      angm = angm + d(i,j,k)*dvol(i,j,k)*cross(xcar-nspos,vcar-nsvel)
!!$      
!!$     end if
!!$     if(eint(i,j,k)/d(i,j,k)+0.5d0*dot_product(vcar-vcom,vcar-vcom)+grvphi(i,j,k)<=0d0)then ! bound to star
!!$      if(boundtoNS)then
!!$       if(norm2(xcar-xcom) < norm2(xcar-nspos))then !closer to star
!!$        boundtostar=.true.
!!$        boundtoNS=.false.
!!$       end if
!!$      else
!!$       boundtostar=.true.       
!!$      end if
!!$     end if
!!$     if(boundtoNS)then
!!$      accmass = accmass + d(i,j,k)*dvol(i,j,k)
!!$      angm = angm + d(i,j,k)*dvol(i,j,k)*cross(xcar-nspos,vcar-nsvel)
!!$     elseif(boundtostar)then
!!$      mstar = mstar + d(i,j,k)*dvol(i,j,k)
!!$      spin = spin + d(i,j,k)*dvol(i,j,k)*cross(xcar-xcom,vcar-vcom)
!!$     else
!!$      Mej = Mej + d(i,j,k)*dvol(i,j,k)
!!$     end if
!!$    end do
!!$   end do
!!$  end do
!!$!$omp end parallel do
!!$
!!$  accmass = 2d0*accmass
!!$  accJ = 2d0*angm(3)
!!$  mstar = 2d0*mstar
!!$  Jspin = 2d0*spin(3)
!!$  Mej = 2d0*Mej
!!$
!!$ return
!!$ end subroutine get_accretion


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE GET_ISOCONTOUR
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get isocontour surface data

!!$ subroutine get_isocontour(xc,levels,data,contour)
!!$
!!$  use grid
!!$  use constants,only:pi
!!$  use utils,only:carpol
!!$
!!$  implicit none
!!$
!!$  real(8),intent(in):: xc(1:3)
!!$  real(8),allocatable,intent(in):: data(:,:,:), levels(:)
!!$  real(8),allocatable,intent(inout):: contour(:,:,:,:)
!!$  integer:: isize(4), ii, jj, kk, ll, iie, jje, kke, lle
!!$  real(8):: dr(0:3), r(1:3), theta, phip, rpol(1:3), rho, rho_old
!!$
!!$!-----------------------------------------------------------------------------
!!$
!!$  isize = shape(contour)
!!$  jje = isize(3)-1
!!$  kke = isize(4)-1
!!$  lle = size(levels)
!!$  dr(0) = minval(dxi1(is:ie))
!!$
!!$  contour = xi1(ie)*3d0 ! outside the computational domain
!!$
!!$  do kk = 0, kke
!!$   phip   = dble(kk)/dble(kke)*2d0*pi-pi
!!$   do jj = 0, jje
!!$    r = xc
!!$    theta = dble(jj)/dble(jje)*pi*0.5d0
!!$    if(jj==0)theta = 1d-6
!!$
!!$    dr(1) = dr(0)*sin(theta)*cos(phip)
!!$    dr(2) = dr(0)*sin(theta)*sin(phip)
!!$    dr(3) = dr(0)*cos(theta)
!!$    rho=1d99
!!$    do while (norm2(r)<xi1(ie))
!!$     r=r+dr(1:3)
!!$     rpol = carpol(r)
!!$     rho_old = rho
!!$     call get_local_val3D(rpol,data,rho)
!!$     do ll = 1, lle
!!$      if(rho<=levels(ll).and.rho_old>levels(ll))then
!!$!       if(contour(1,ll,jj,kk)>xi1(ie))
!!$       contour(1:3,ll,jj,kk) = r-dr(1:3)*0.5d0
!!$      end if
!!$     end do
!!$     if(contour(1,lle,jj,kk)<=xi1(ie))exit
!!$    end do
!!$   end do
!!$  end do
!!$  
!!$
!!$ return
!!$end subroutine get_isocontour

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LOCAL_GRAV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculate simplified gravitational acceleration on NS

subroutine get_local_grav(x,y,acc)

 use constants,only:G
 use utils,only:intpol

 implicit none

 real(8),intent(in):: x(1:3),y(1:3)
 real(8),intent(out):: acc(1:3)
 real(8):: relx(1:3), mass, dis
 integer i
 
!-----------------------------------------------------------------------------

! x is the companion star
! y is the neutron star
 relx = y-x
 dis = norm2(relx)

 mass = 0d0
 if(dis>=r(lines))then
  mass = m(lines)
 else
  do i = 1, lines-1
   if(dis>=r(i).and.dis<r(i+1))then
    mass = intpol(r(i:i+1),m(i:i+1),dis)
    exit
   end if
  end do
 end if
 
 acc = -G*mass*relx/dis**3
 

return
end subroutine get_local_grav

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE GET_ANAL_DRAG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculate the analytic drag

!!$subroutine get_anal_drag(vrel,rrel,drag,cs,dragHL)
!!$
!!$ use constants,only:G,pi
!!$ use utils,only:intpol
!!$ use gravmod,only:nsmass
!!$
!!$ implicit none
!!$
!!$ real(8),intent(in):: vrel, rrel
!!$ real(8),intent(out):: drag(1:2), cs, dragHL
!!$ integer i, j
!!$ real(8):: mass, dens, p, Rbondi, rhomax,rhomin, Rrange
!!$ 
!!$!-----------------------------------------------------------------------------
!!$
!!$ cs=1d-99
!!$ drag=0d0
!!$ dragHL=0d0
!!$! if(rrel>=r(lines))then
!!$!  drag = 0d0
!!$! else
!!$!  do i = 1, lines-1
!!$!   if(rrel>=r(i).and.rrel<r(i+1))then
!!$!    mass = intpol(r(i:i+1),m(i:i+1),rrel)
!!$!    dens = intpol(r(i:i+1),rho(i:i+1),rrel)
!!$!    p = intpol(r(i:i+1),pres(i:i+1),rrel)
!!$!    Rbondi = 2d0*G*nsmass/vrel**2
!!$!    cs = sqrt(5d0/3d0*p/dens)
!!$!!    drag = 4d0*pi*G**2*nsmass/vrel**2*dens
!!$!    do j = 1, i
!!$!     if(rrel-Rbondi>=r(j).and.rrel-Rbondi<r(j+1))then
!!$!      dens = intpol(r(j:j+1),rho(j:j+1),rrel-Rbondi)
!!$!      exit
!!$!     end if
!!$!    end do
!!$!    drag = pi*(Rbondi*vrel)**2*dens/nsmass
!!$!    exit
!!$!   end if
!!$!  end do
!!$! end if
!!$
!!$ Rbondi = 2d0*G*nsmass/vrel**2
!!$ if(rrel>=r(lines))then
!!$  drag = 0d0
!!$ else
!!$  do i = 1, lines-1
!!$   if(rrel>=r(i).and.rrel<r(i+1))then
!!$    mass = intpol(r(i:i+1),m(i:i+1),rrel)
!!$    dens = intpol(r(i:i+1),rho(i:i+1),rrel)
!!$    p = intpol(r(i:i+1),pres(i:i+1),rrel)
!!$    cs = sqrt(5d0/3d0*p/dens)
!!$    dragHL = 4d0*pi*G**2*nsmass/vrel**2*dens
!!$    do j = i, 1, -1
!!$     if(m(j)/r(j)**2>nsmass/(rrel-r(j))**2.or.r(j)<0.75d0*rrel)then
!!$!     if(rrel-r(j)>2d0*Rbondi)then
!!$      rhomax = rho(j)
!!$      Rrange = rrel-r(j)
!!$      exit
!!$     end if
!!$    end do
!!$!     if(rrel-Rbondi>=r(j).and.rrel-Rbondi<r(j+1))then
!!$!      rhomax = intpol(r(j:j+1),rho(j:j+1),rrel-Rbondi)
!!$    do j = i,lines-1
!!$     if(rrel+Rrange>=r(lines))then
!!$      rhomin = 0d0
!!$      exit
!!$     elseif(rrel+Rrange>=r(j).and.rrel+Rrange<r(j+1))then
!!$      rhomin = intpol(r(j:j+1),rho(j:j+1),rrel+Rrange)
!!$      exit
!!$     end if
!!$    end do
!!$    dens = 0.5d0*(rhomax+rhomin)
!!$    Rrange = Rrange/(G*nsmass/vrel**2)
!!$    print'(4(1PE13.5e2))',rrel/(G*nsmass/vrel**2),Rrange,rhomax,rhomin
!!$    drag(1) = 0.5d0*pi*log(Rrange**2+1d0)*(Rbondi*vrel)**2*dens/nsmass
!!$    drag(2) = pi/8d0*(Rrange**2-log(Rrange**2+1d0))*(rhomax-rhomin)*(Rbondi*vrel)**2/nsmass/Rrange
!!$    exit
!!$   end if
!!$  end do
!!$ end if
!!$
!!$return
!!$end subroutine get_anal_drag
! convert cartesian to polar coordinates
end module analysis_mod
