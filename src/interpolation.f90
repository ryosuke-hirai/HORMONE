!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE INTERPOLATION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Calculates interpolated values at cell boundaries based on MUSCL 
!          interpolation.

subroutine interpolation

 use settings,only:compswitch,spn,eostype
 use grid
 use physval
 use ninewave
 use fluxlimiter
 use pressure_mod

 implicit none

 real*8 dl, dr, ptl, ptr, el, er, m1l, m1r, m2l, m2r, m3l, m3r
 real*8 b1l, b1r, b2l, b2r, b3l, b3r, phil, phir, eintl, eintr,imul,imur
 real*8 uu(1:3), du, Xl, Xr, Yl, Yr
 real*8 dx(1:2), x(1:3), xi(1:2)

!-----------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Notations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! slopes : dd, de, dm1, dm2, dm3, db1, db2, db3
! *l, *r are cell boundary values at left and right looking from x1(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! slope1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,m1l,m1r,m2l,m2r,m3l,m3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,phil,phir,uu,du,dx,eintl,eintr,imul,imur,n,x,xi,&
!$omp Xl,Xr,Yl,Yr)

 do k = ks,ke
  do j = js,je
   do i = is-1, ie+1
    dx(1:2) = idx1(i:i+1)
    x(1:3) = x1(i-1:i+1) ; xi(1:2) = xi1(i-1:i)

    uu(1:3) = d(i-1:i+1,j,k)
!    call minmod(du,uu,dx) ; dd(i,j,k,1) = du
    call modified_mc(du,uu,x,xi) ; dd(i,j,k,1) = du
    dl = uu(2) - (x1(i)-xi1(i-1))*du ; dr = uu(2) + (xi1(i)-x1(i))*du
    
    uu(1:3) = e(i-1:i+1,j,k)
!    call minmod(du,uu,dx) ; de(i,j,k,1) = du
    call modified_mc(du,uu,x,xi) ; de(i,j,k,1) = du
    el = uu(2) - (x1(i)-xi1(i-1))*du ; er = uu(2) + (xi1(i)-x1(i))*du

    uu(1:3) = d(i-1:i+1,j,k)*v1(i-1:i+1,j,k)
!    call minmod(du,uu,dx) ; dm1(i,j,k,1) = du
    call modified_mc(du,uu,x,xi) ; dm1(i,j,k,1) = du
    m1l = uu(2) - (x1(i)-xi1(i-1))*du ; m1r = uu(2) + (xi1(i)-x1(i))*du

    uu(1:3) = d(i-1:i+1,j,k)*v2(i-1:i+1,j,k)
!    call minmod(du,uu,dx) ; dm2(i,j,k,1) = du
    call modified_mc(du,uu,x,xi) ; dm2(i,j,k,1) = du
    m2l = uu(2) - (x1(i)-xi1(i-1))*du ; m2r = uu(2) + (xi1(i)-x1(i))*du

    uu(1:3) = d(i-1:i+1,j,k)*v3(i-1:i+1,j,k)
!    call minmod(du,uu,dx) ; dm3(i,j,k,1) = du
    call modified_mc(du,uu,x,xi) ; dm3(i,j,k,1) = du
    m3l = uu(2) - (x1(i)-xi1(i-1))*du ; m3r = uu(2) + (xi1(i)-x1(i))*du

!!$    uu(1:3) = b1(i-1:i+1,j,k)
!!$    call minmod(du,uu,dx) ; db1(i,j,k,1) = du
!!$    b1l = uu(2) - (x1(i)-xi1(i-1))*du ; b1r = uu(2) + (xi1(i)-x1(i))*du
!!$
!!$    uu(1:3) = b2(i-1:i+1,j,k)
!!$    call minmod(du,uu,dx) ; db2(i,j,k,1) = du
!!$    b2l = uu(2) - (x1(i)-xi1(i-1))*du ; b2r = uu(2) + (xi1(i)-x1(i))*du
!!$
!!$    uu(1:3) = b3(i-1:i+1,j,k)
!!$    call minmod(du,uu,dx) ; db3(i,j,k,1) = du
!!$    b3l = uu(2) - (x1(i)-xi1(i-1))*du ; b3r = uu(2) + (xi1(i)-x1(i))*du

! check stability at cell boundary ---------------------------------------- !
    ! first calculate mean molecular weight at the surface
    select case (compswitch)
    case(0) ! uniform composition
     imul = 1d0/muconst ; imur = 1d0/muconst ; dmu(i,j,k,1) = 0d0
    case(1:2) ! nonuniform composition
     uu(1:3) = 1d0/imu(i-1:i+1,j,k)
     call minmod(du,uu,dx) ; dmu(i,j,k,1) = du
     imul = uu(2) - (x1(i)-xi1(i-1))*du ; imur = uu(2) + (xi1(i)-x1(i))*du
     imul = 1d0/imul ; imur = 1d0/imur
     if(compswitch==2)then
      do n = 1, spn
       uu(1:3) = spc(n,i-1:i+1,j,k)
       call minmod(du,uu,dx) ; dspc(n,i,j,k,1) = du
      end do
     end if
    case default
     print *, "Error in compswitch",compswitch
     stop
    end select
    ! then calculate the energies at boundaries
    eintl = el - 0.5d0*(m1l*m1l+m2l*m2l+m3l*m3l)/dl
    eintr = er - 0.5d0*(m1r*m1r+m2r*m2r+m3r*m3r)/dr
    if( eintl>=maxval(eint(i-1:i,j,k)).or.eintl<=minval(eint(i-1:i,j,k)).or.&
        eintr>=maxval(eint(i:i+1,j,k)).or.eintr<=minval(eint(i:i+1,j,k)))then
     dd (i,j,k,1) = 0d0 ; de (i,j,k,1) = 0d0
     dm1(i,j,k,1) = 0d0 ; dm2(i,j,k,1) = 0d0 ; dm3(i,j,k,1) = 0d0
     db1(i,j,k,1) = 0d0 ; db2(i,j,k,1) = 0d0 ; db3(i,j,k,1) = 0d0
    else
     select case (eostype)
     case(0:1) ! without recombination
      ptl = eos_p(dl,eintl,T(i,j,k),imul) ; ptr = eos_p(dr,eintr,T(i,j,k),imur)
     case(2) ! with recombination
      Xl = spc(1,i,j,k)-(x1(i)-xi1(i-1))* dspc(1,i,j,k,1)
      Yl = spc(2,i,j,k)-(x1(i)-xi1(i-1))* dspc(2,i,j,k,1)
      Xr = spc(1,i,j,k)+(xi1(i)-x1(i  ))* dspc(1,i,j,k,1)
      Yr = spc(2,i,j,k)+(xi1(i)-x1(i  ))* dspc(2,i,j,k,1)
      ptl = eos_p(dl,eintl,T(i,j,k),imul,Xl,Yl)
      ptr = eos_p(dr,eintr,T(i,j,k),imur,Xr,Yr)
     end select
     if( ptl>=maxval(ptot(i-1:i,j,k)).or.ptl<=minval(ptot(i-1:i,j,k)).or.&
         ptr>=maxval(ptot(i:i+1,j,k)).or.ptr<=minval(ptot(i:i+1,j,k)))then
      dd (i,j,k,1) = 0d0 ; de (i,j,k,1) = 0d0
      dm1(i,j,k,1) = 0d0 ; dm2(i,j,k,1) = 0d0 ; dm3(i,j,k,1) = 0d0
      db1(i,j,k,1) = 0d0 ; db2(i,j,k,1) = 0d0 ; db3(i,j,k,1) = 0d0
     end if
    end if
! ------------------------------------------------------------------------- !

    uu(1:3) = phi(i-1:i+1,j,k)
    call minmod(du,uu,dx) ; dphi(i,j,k,1) = du

   end do
  end do
 end do
!$omp end parallel do

! slope2 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if(je/=1)then
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,m1l,m1r,m2l,m2r,m3l,m3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,phil,phir,uu,du,dx,eintl,eintr,imul,imur,n,x,xi,&
!$omp Xl,Xr,Yl,Yr)
 do k = ks,ke
  do j = js-1,je+1
   do i = is, ie
    dx(1:2) = idx2(j:j+1)
    x(1:3) = x2(j-1:j+1) ; xi(1:2) = xi2(j-1:j)

    uu(1:3) = d(i,j-1:j+1,k)
!    call minmod(du,uu,dx) ; dd(i,j,k,2) = du
    call modified_mc(du,uu,x,xi) ; dd(i,j,k,2) = du
    dl = uu(2) - (x2(j)-xi2(j-1))*du ; dr = uu(2) + (xi2(j)-x2(j))*du

    uu(1:3) = e(i,j-1:j+1,k)
!    call minmod(du,uu,dx) ; de(i,j,k,2) = du
    call modified_mc(du,uu,x,xi) ; de(i,j,k,2) = du
    el = uu(2) - (x2(j)-xi2(j-1))*du ; er = uu(2) + (xi2(j)-x2(j))*du

    uu(1:3) = d(i,j-1:j+1,k)*v1(i,j-1:j+1,k)
!    call minmod(du,uu,dx) ; dm1(i,j,k,2) = du
    call modified_mc(du,uu,x,xi) ; dm1(i,j,k,2) = du
    m1l = uu(2) - (x2(j)-xi2(j-1))*du ; m1r = uu(2) + (xi2(j)-x2(j))*du

    uu(1:3) = d(i,j-1:j+1,k)*v2(i,j-1:j+1,k)
!    call minmod(du,uu,dx) ; dm2(i,j,k,2) = du
    call modified_mc(du,uu,x,xi) ; dm2(i,j,k,2) = du
    m2l = uu(2) - (x2(j)-xi2(j-1))*du ; m2r = uu(2) + (xi2(j)-x2(j))*du

    uu(1:3) = d(i,j-1:j+1,k)*v3(i,j-1:j+1,k)
!    call minmod(du,uu,dx) ; dm3(i,j,k,2) = du
    call modified_mc(du,uu,x,xi) ; dm3(i,j,k,2) = du
    m3l = uu(2) - (x2(j)-xi2(j-1))*du ; m3r = uu(2) + (xi2(j)-x2(j))*du

!!$    uu(1:3) = b1(i,j-1:j+1,k)
!!$    call minmod(du,uu,dx) ; db1(i,j,k,2) = du
!!$    b1l = uu(2) - (x2(j)-xi2(j-1))*du ; b1r = uu(2) + (xi2(j)-x2(j))*du
!!$
!!$    uu(1:3) = b2(i,j-1:j+1,k)
!!$    call minmod(du,uu,dx) ; db2(i,j,k,2) = du
!!$    b2l = uu(2) - (x2(j)-xi2(j-1))*du ; b2r = uu(2) + (xi2(j)-x2(j))*du
!!$
!!$    uu(1:3) = b3(i,j-1:j+1,k)
!!$    call minmod(du,uu,dx) ; db3(i,j,k,2) = du
!!$    b3l = uu(2) - (x2(j)-xi2(j-1))*du ; b3r = uu(2) + (xi2(j)-x2(j))*du

! check stability at cell boundary ---------------------------------------- !
    ! first calculate mean molecular weight at boundaries
    select case (compswitch)
    case(0) ! uniform composition
     imul = 1d0/muconst ; imur = 1d0/muconst ; dmu(i,j,k,2) = 0d0
    case(1:2) ! nonuniform composition
     uu(1:3) = 1d0/imu(i,j-1:j+1,k)
     call minmod(du,uu,dx) ; dmu(i,j,k,2) = du
     imul = uu(2) - (x2(j)-xi2(j-1))*du ; imur = uu(2) + (xi2(j)-x2(j))*du
     imul = 1d0/imul ; imur = 1d0/imur
     if(compswitch==2)then
      do n = 1, spn
       uu(1:3) = spc(n,i,j-1:j+1,k)
       call minmod(du,uu,dx) ; dspc(n,i,j,k,2) = du
      end do
     end if
    case default
     print *, "Error in compswitch",compswitch
     stop
    end select
    ! then calculate the energies at boundaries
    eintl = el - 0.5d0*(m1l*m1l+m2l*m2l+m3l*m3l)/dl
!               - 0.5d0*(b1l*b1l+b2l*b2l+b3l*b3l)
    eintr = er - 0.5d0*(m1r*m1r+m2r*m2r+m3r*m3r)/dr
!               - 0.5d0*(b1r*b1r+b2r*b2r+b3r*b3r)
    if( eintl>=maxval(eint(i,j-1:j,k)).or.eintl<=minval(eint(i,j-1:j,k)).or.&
        eintr>=maxval(eint(i,j:j+1,k)).or.eintr<=minval(eint(i,j:j+1,k)))then
     dd (i,j,k,2) = 0d0 ; de (i,j,k,2) = 0d0
     dm1(i,j,k,2) = 0d0 ; dm2(i,j,k,2) = 0d0 ; dm3(i,j,k,2) = 0d0
     db1(i,j,k,2) = 0d0 ; db2(i,j,k,2) = 0d0 ; db3(i,j,k,2) = 0d0
    else
     select case (eostype)
     case(0:1) ! without recombination
      ptl = eos_p(dl,eintl,T(i,j,k),imul) ; ptr = eos_p(dr,eintr,T(i,j,k),imur)
     case(2) ! with recombination
      Xl = spc(1,i,j,k)-(x2(j)-xi2(j-1))* dspc(1,i,j,k,2)
      Yl = spc(2,i,j,k)-(x2(j)-xi2(j-1))* dspc(2,i,j,k,2)
      Xr = spc(1,i,j,k)+(xi2(j)-x2(j  ))* dspc(1,i,j,k,2)
      Yr = spc(2,i,j,k)+(xi2(j)-x2(j  ))* dspc(2,i,j,k,2)
      ptl = eos_p(dl,eintl,T(i,j,k),imul,Xl,Yl)
      ptr = eos_p(dr,eintr,T(i,j,k),imur,Xr,Yr)
     end select
     if( ptl>=maxval(ptot(i,j-1:j,k)).or.ptl<=minval(ptot(i,j-1:j,k)).or.&
         ptr>=maxval(ptot(i,j:j+1,k)).or.ptr<=minval(ptot(i,j:j+1,k)))then
      dd (i,j,k,2) = 0d0 ; de (i,j,k,2) = 0d0
      dm1(i,j,k,2) = 0d0 ; dm2(i,j,k,2) = 0d0 ; dm3(i,j,k,2) = 0d0
      db1(i,j,k,2) = 0d0 ; db2(i,j,k,2) = 0d0 ; db3(i,j,k,2) = 0d0
     end if
    end if
! -------------------------------------------------------------------------- !

    uu(1:3) = phi(i,j-1:j+1,k)
    call minmod(du,uu,dx) ; dphi(i,j,k,2) = du

   end do
  end do
 end do
!$omp end parallel do
end if
! slope3 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if(ke/=1)then
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,m1l,m1r,m2l,m2r,m3l,m3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,phil,phir,uu,du,dx,eintl,eintr,imul,imur,n,x,xi,&
!$omp Xl,Xr,Yl,Yr)

 do k = ks-1,ke+1
  do j = js,je
   do i = is, ie
    dx(1:2) = idx3(k:k+1)
    x(1:3) = x3(k-1:k+1) ; xi(1:2) = xi3(k-1:k)

    uu(1:3) = d(i,j,k-1:k+1)
!    call minmod(du,uu,dx) ; dd(i,j,k,3) = du
    call modified_mc(du,uu,x,xi) ; dd(i,j,k,3) = du
    dl = uu(2) - (x3(k)-xi3(k-1))*du ; dr = uu(2) + (xi3(k)-x3(k))*du

    uu(1:3) = e(i,j,k-1:k+1)
    !    call minmod(du,uu,dx) ; de(i,j,k,3) = du
    call modified_mc(du,uu,x,xi) ; de(i,j,k,3) = du
    el = uu(2) - (x3(k)-xi3(k-1))*du ; er = uu(2) + (xi3(k)-x3(k))*du

    uu(1:3) = d(i,j,k-1:k+1)*v1(i,j,k-1:k+1)
!    call minmod(du,uu,dx) ; dm1(i,j,k,3) = du
    call modified_mc(du,uu,x,xi) ; dm1(i,j,k,3) = du
    m1l = uu(2) - (x3(k)-xi3(k-1))*du ; m1r = uu(2) + (xi3(k)-x3(k))*du

    uu(1:3) = d(i,j,k-1:k+1)*v2(i,j,k-1:k+1)
!    call minmod(du,uu,dx) ; dm2(i,j,k,3) = du
    call modified_mc(du,uu,x,xi) ; dm2(i,j,k,3) = du
    m2l = uu(2) - (x3(k)-xi3(k-1))*du ; m2r = uu(2) + (xi3(k)-x3(k))*du

    uu(1:3) = d(i,j,k-1:k+1)*v3(i,j,k-1:k+1)
!    call minmod(du,uu,dx) ; dm3(i,j,k,3) = du
    call modified_mc(du,uu,x,xi) ; dm3(i,j,k,3) = du
    m3l = uu(2) - (x3(k)-xi3(k-1))*du ; m3r = uu(2) + (xi3(k)-x3(k))*du

!!$    uu(1:3) = b1(i,j,k-1:k+1)
!!$    call minmod(du,uu,dx) ; db1(i,j,k,3) = du
!!$    b1l = uu(2) - (x3(k)-xi3(k-1))*du ; b1r = uu(2) + (xi3(k)-x3(k))*du
!!$
!!$    uu(1:3) = b2(i,j,k-1:k+1)
!!$    call minmod(du,uu,dx) ; db2(i,j,k,3) = du
!!$    b2l = uu(2) - (x3(k)-xi3(k-1))*du ; b2r = uu(2) + (xi3(k)-x3(k))*du
!!$
!!$    uu(1:3) = b3(i,j,k-1:k+1)
!!$    call minmod(du,uu,dx) ; db3(i,j,k,3) = du
!!$    b3l = uu(2) - (x3(k)-xi3(k-1))*du ; b3r = uu(2) + (xi3(k)-x3(k))*du

! check stability at cell boundary ---------------------------------------- !
    ! first calculate mean molecular weight at boundaries
    select case (compswitch)
    case(0) ! uniform composition
     imul = 1d0/muconst ; imur = 1d0/muconst ; dmu(i,j,k,3) = 0d0
    case(1:2) ! nonuniform composition
     uu(1:3) = 1d0/imu(i,j,k-1:k+1)
     call minmod(du,uu,dx) ; dmu(i,j,k,3) = du
     imul = uu(2) - (x3(k)-xi3(k-1))*du ; imur = uu(2) + (xi3(k)-x3(k))*du
     imul = 1d0/imul ; imur = 1d0/imur
     if(compswitch==2)then
      do n = 1, spn
       uu(1:3) = spc(n,i,j,k-1:k+1)
       call minmod(du,uu,dx) ; dspc(n,i,j,k,3) = du
      end do
     end if
    case default
     print *, "Error in compswitch",compswitch
     stop
    end select
    ! then calculate the condition at boundaries
    eintl = el - 0.5d0*(m1l*m1l+m2l*m2l+m3l*m3l)/dl
!               - 0.5d0*(b1l*b1l+b2l*b2l+b3l*b3l)
    eintr = er - 0.5d0*(m1r*m1r+m2r*m2r+m3r*m3r)/dr
!               - 0.5d0*(b1r*b1r+b2r*b2r+b3r*b3r)
    if( eintl>=maxval(eint(i,j,k-1:k)).or.eintl<=minval(eint(i,j,k-1:k)).or.&
        eintr>=maxval(eint(i,j,k:k+1)).or.eintr<=minval(eint(i,j,k:k+1)))then
     dd (i,j,k,3) = 0d0 ; de (i,j,k,3) = 0d0
     dm1(i,j,k,3) = 0d0 ; dm2(i,j,k,3) = 0d0 ; dm3(i,j,k,3) = 0d0
     db1(i,j,k,3) = 0d0 ; db2(i,j,k,3) = 0d0 ; db3(i,j,k,3) = 0d0
    else
     select case (eostype)
     case(0:1) ! without recombination
      ptl = eos_p(dl,eintl,T(i,j,k),imul) ; ptr = eos_p(dr,eintr,T(i,j,k),imur)
     case(2) ! with recombination
      Xl = spc(1,i,j,k)-(x3(k)-xi3(k-1))* dspc(1,i,j,k,3)
      Yl = spc(2,i,j,k)-(x3(k)-xi3(k-1))* dspc(2,i,j,k,3)
      Xr = spc(1,i,j,k)+(xi3(k)-x3(k  ))* dspc(1,i,j,k,3)
      Yr = spc(2,i,j,k)+(xi3(k)-x3(k  ))* dspc(2,i,j,k,3)
      ptl = eos_p(dl,eintl,T(i,j,k),imul,Xl,Yl)
      ptr = eos_p(dr,eintr,T(i,j,k),imur,Xr,Yr)
     end select
     if( ptl>=maxval(ptot(i,j,k-1:k)).or.ptl<=minval(ptot(i,j,k-1:k)).or.&
         ptr>=maxval(ptot(i,j,k:k+1)).or.ptr<=minval(ptot(i,j,k:k+1)))then
      dd (i,j,k,3) = 0d0 ; de (i,j,k,3) = 0d0
      dm1(i,j,k,3) = 0d0 ; dm2(i,j,k,3) = 0d0 ; dm3(i,j,k,3) = 0d0
      db1(i,j,k,3) = 0d0 ; db2(i,j,k,3) = 0d0 ; db3(i,j,k,3) = 0d0
     end if
    end if
! -------------------------------------------------------------------------- !

    uu(1:3) = phi(i,j,k-1:k+1)
    call minmod(du,uu,dx) ; dphi(i,j,k,3) = du

   end do
  end do
 end do
!$omp end parallel do
end if

return
end subroutine interpolation
