!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE AMR_NUMFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate numerical flux for each cell in AMR blocks.

subroutine amr_numflux(lf)

 use grid,only:is,js,ks
 use hlldflux_mod
 use amr_templates
 use amr_module,only:ib,jb,kb
 use amr_eigen_mod
 use pressure_mod

 implicit none

 integer i,j,k
 real*8 cfl, cfr, v1l, v1r, dl, dr, ptl, ptr, el, er
 real*8 b1l, b1r, b2l, b2r, b3l, b3r, v2l, v2r, v3l, v3r, phil, phir
 real*8,dimension(1:9)::tmpflux
 real*8,dimension(4):: u
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

 call amr_conserve(lf)
 call pressure_block(lf)
 call amr_eigen(lf)

! use MUSCL interpolation (with minmod limiter)

! flux1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,tmpflux,u)
 do k = ks, kb
  do j = js, jb
   do i = is-1, ib
    u(1:4) = lf%ptot(i-1:i+2,j,k)
    call amr_minmod(ptl,u(1:3)) ; call amr_minmod(ptr,u(2:4))
    ptl = u(2) + 5.d-1*ptl ; ptr = u(3) - 5.d-1*ptr

    u(1:4) = lf%d(i-1:i+2,j,k)
    call amr_minmod(dl,u(1:3))  ; call amr_minmod(dr,u(2:4))
    dl = u(2) + 5.d-1*dl ; dr = u(3) - 5.d-1*dr

    u(1:4) = lf%e(i-1:i+2,j,k)
    call amr_minmod(el,u(1:3))  ; call amr_minmod(er,u(2:4))
    el = u(2) + 5.d-1*el ; er = u(3) - 5.d-1*er

    u(1:4) = lf%v1(i-1:i+2,j,k)
    call amr_minmod(v1l,u(1:3))  ; call amr_minmod(v1r,u(2:4))
    v1l = u(2) + 5.d-1*v1l ; v1r = u(3) - 5.d-1*v1r

    u(1:4) = lf%v2(i-1:i+2,j,k)
    call amr_minmod(v2l,u(1:3))  ; call amr_minmod(v2r,u(2:4))
    v2l = u(2) + 5.d-1*v2l ; v2r = u(3) - 5.d-1*v2r

    u(1:4) = lf%v3(i-1:i+2,j,k)
    call amr_minmod(v3l,u(1:3))  ; call amr_minmod(v3r,u(2:4))
    v3l = u(2) + 5.d-1*v3l ; v3r = u(3) - 5.d-1*v3r

    u(1:4) = lf%b1(i-1:i+2,j,k)
    call amr_minmod(b1l,u(1:3))  ; call amr_minmod(b1r,u(2:4))
    b1l = u(2) + 5.d-1*b1l ; b1r = u(3) - 5.d-1*b1r

    u(1:4) = lf%b2(i-1:i+2,j,k)
    call amr_minmod(b2l,u(1:3))  ; call amr_minmod(b2r,u(2:4))
    b2l = u(2) + 5.d-1*b2l ; b2r = u(3) - 5.d-1*b2r

    u(1:4) = lf%b3(i-1:i+2,j,k)
    call amr_minmod(b3l,u(1:3))  ; call amr_minmod(b3r,u(2:4))
    b3l = u(2) + 5.d-1*b3l ; b3r = u(3) - 5.d-1*b3r

    u(1:4) = lf%cf(i-1:i+2,j,k)
    call amr_minmod(cfl,u(1:3))  ; call amr_minmod(cfr,u(2:4))
    cfl = u(2) + 5.d-1*cfl ; cfr = u(3) - 5.d-1*cfr

    u(1:4) = lf%phi(i-1:i+2,j,k)
    call amr_minmod(phil,u(1:3))  ; call amr_minmod(phir,u(2:4))
    phil = u(2) + 5.d-1*phil ; phir = u(3) - 5.d-1*phir
!print *,i,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
!          el,er,ptl,ptr
    call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

    lf%flux1(i,j,k,1:9) = tmpflux
!print *,i,lf%flux1(i,j,k,1),lf%e(i,j,k)
   end do
  end do
 end do

!$omp end parallel do

! flux2 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if(jb>1)then
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,tmpflux,u)
 do k = ks, kb
  do j = js-1, jb
   do i = is, ib
    u(1:4) = lf%ptot(i,j-1:j+2,k)
    call amr_minmod(ptl,u(1:3)) ; call amr_minmod(ptr,u(2:4))
    ptl = u(2) + 5.d-1*ptl ; ptr = u(3) - 5.d-1*ptr

    u(1:4) = lf%d(i,j-1:j+2,k)
    call amr_minmod(dl,u(1:3))  ; call amr_minmod(dr,u(2:4))
    dl = u(2) + 5.d-1*dl ; dr = u(3) - 5.d-1*dr

    u(1:4) = lf%e(i,j-1:j+2,k)
    call amr_minmod(el,u(1:3))  ; call amr_minmod(er,u(2:4))
    el = u(2) + 5.d-1*el ; er = u(3) - 5.d-1*er

    u(1:4) = lf%v2(i,j-1:j+2,k)
    call amr_minmod(v1l,u(1:3))  ; call amr_minmod(v1r,u(2:4))
    v1l = u(2) + 5.d-1*v1l ; v1r = u(3) - 5.d-1*v1r

    u(1:4) = lf%v3(i,j-1:j+2,k)
    call amr_minmod(v2l,u(1:3))  ; call amr_minmod(v2r,u(2:4))
    v2l = u(2) + 5.d-1*v2l ; v2r = u(3) - 5.d-1*v2r

    u(1:4) = lf%v1(i,j-1:j+2,k)
    call amr_minmod(v3l,u(1:3))  ; call amr_minmod(v3r,u(2:4))
    v3l = u(2) + 5.d-1*v3l ; v3r = u(3) - 5.d-1*v3r

    u(1:4) = lf%b2(i,j-1:j+2,k)
    call amr_minmod(b1l,u(1:3))  ; call amr_minmod(b1r,u(2:4))
    b1l = u(2) + 5.d-1*b1l ; b1r = u(3) - 5.d-1*b1r

    u(1:4) = lf%b3(i,j-1:j+2,k)
    call amr_minmod(b2l,u(1:3))  ; call amr_minmod(b2r,u(2:4))
    b2l = u(2) + 5.d-1*b2l ; b2r = u(3) - 5.d-1*b2r

    u(1:4) = lf%b1(i,j-1:j+2,k)
    call amr_minmod(b3l,u(1:3))  ; call amr_minmod(b3r,u(2:4))
    b3l = u(2) + 5.d-1*b3l ; b3r = u(3) - 5.d-1*b3r

    u(1:4) = lf%cf(i,j-1:j+2,k)
    call amr_minmod(cfl,u(1:3))  ; call amr_minmod(cfr,u(2:4))
    cfl = u(2) + 5.d-1*cfl ; cfr = u(3) - 5.d-1*cfr

    u(1:4) = lf%phi(i,j-1:j+2,k)
    call amr_minmod(phil,u(1:3))  ; call amr_minmod(phir,u(2:4))
    phil = u(2) + 5.d-1*phil ; phir = u(3) - 5.d-1*phir

    call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

    lf%flux2(i,j,k,1) = tmpflux(1)
    lf%flux2(i,j,k,2) = tmpflux(4)
    lf%flux2(i,j,k,3) = tmpflux(2)
    lf%flux2(i,j,k,4) = tmpflux(3)
    lf%flux2(i,j,k,5) = tmpflux(7)
    lf%flux2(i,j,k,6) = tmpflux(5)
    lf%flux2(i,j,k,7) = tmpflux(6)
    lf%flux2(i,j,k,8) = tmpflux(8)
    lf%flux2(i,j,k,9) = tmpflux(9)

   end do
  end do
 end do
!$omp end parallel do
end if

! flux3 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if(kb>1)then
!$omp parallel do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,tmpflux,u)
 do k = ks-1, kb
  do j = js, jb
   do i = is, ib
    u(1:4) = lf%ptot(i,j,k-1:k+2)
    call amr_minmod(ptl,u(1:3)) ; call amr_minmod(ptr,u(2:4))
    ptl = u(2) + 5.d-1*ptl ; ptr = u(3) - 5.d-1*ptr

    u(1:4) = lf%d(i,j,k-1:k+2)
    call amr_minmod(dl,u(1:3))  ; call amr_minmod(dr,u(2:4))
    dl = u(2) + 5.d-1*dl ; dr = u(3) - 5.d-1*dr

    u(1:4) = lf%e(i,j,k-1:k+2)
    call amr_minmod(el,u(1:3))  ; call amr_minmod(er,u(2:4))
    el = u(2) + 5.d-1*el ; er = u(3) - 5.d-1*er

    u(1:4) = lf%v3(i,j,k-1:k+2)
    call amr_minmod(v1l,u(1:3))  ; call amr_minmod(v1r,u(2:4))
    v1l = u(2) + 5.d-1*v1l ; v1r = u(3) - 5.d-1*v1r

    u(1:4) = lf%v1(i,j,k-1:k+2)
    call amr_minmod(v2l,u(1:3))  ; call amr_minmod(v2r,u(2:4))
    v2l = u(2) + 5.d-1*v2l ; v2r = u(3) - 5.d-1*v2r

    u(1:4) = lf%v2(i,j,k-1:k+2)
    call amr_minmod(v3l,u(1:3))  ; call amr_minmod(v3r,u(2:4))
    v3l = u(2) + 5.d-1*v3l ; v3r = u(3) - 5.d-1*v3r

    u(1:4) = lf%b3(i,j,k-1:k+2)
    call amr_minmod(b1l,u(1:3))  ; call amr_minmod(b1r,u(2:4))
    b1l = u(2) + 5.d-1*b1l ; b1r = u(3) - 5.d-1*b1r

    u(1:4) = lf%b1(i,j,k-1:k+2)
    call amr_minmod(b2l,u(1:3))  ; call amr_minmod(b2r,u(2:4))
    b2l = u(2) + 5.d-1*b2l ; b2r = u(3) - 5.d-1*b2r

    u(1:4) = lf%b2(i,j,k-1:k+2)
    call amr_minmod(b3l,u(1:3))  ; call amr_minmod(b3r,u(2:4))
    b3l = u(2) + 5.d-1*b3l ; b3r = u(3) - 5.d-1*b3r

    u(1:4) = lf%cf(i,j,k-1:k+2)
    call amr_minmod(cfl,u(1:3))  ; call amr_minmod(cfr,u(2:4))
    cfl = u(2) + 5.d-1*cfl ; cfr = u(3) - 5.d-1*cfr

    u(1:4) = lf%phi(i,j,k-1:k+2)
    call amr_minmod(phil,u(1:3))  ; call amr_minmod(phir,u(2:4))
    phil = u(2) + 5.d-1*phil ; phir = u(3) - 5.d-1*phir

    call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

    lf%flux3(i,j,k,1) = tmpflux(1)
    lf%flux3(i,j,k,2) = tmpflux(3)
    lf%flux3(i,j,k,3) = tmpflux(4)
    lf%flux3(i,j,k,4) = tmpflux(2)
    lf%flux3(i,j,k,5) = tmpflux(6)
    lf%flux3(i,j,k,6) = tmpflux(7)
    lf%flux3(i,j,k,7) = tmpflux(5)
    lf%flux3(i,j,k,8) = tmpflux(8)
    lf%flux3(i,j,k,9) = tmpflux(9)

   end do
  end do
 end do
!$omp end parallel do
end if

if(ib==1)lf%flux1=0.d0; if(jb==1)lf%flux2=0.d0; if(kb==1)lf%flux3=0.d0

return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE AMR_CONSERVE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert physical values to conserved values

subroutine amr_conserve(lf)

 use grid,only:is,js,ks
 use amr_templates
 use amr_module,only:ib,jb,kb

 implicit none

 integer i,j,k
 type(leaf_contents),intent(inout):: lf

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k)
 do k = ks-1, kb+1
  do j = js-1, jb+1
   do i = is-1, ib+1
    lf%u(i,j,k,1) = lf%d(i,j,k)
    lf%u(i,j,k,2) = lf%d(i,j,k) * lf%v1(i,j,k)
    lf%u(i,j,k,3) = lf%d(i,j,k) * lf%v2(i,j,k)
    lf%u(i,j,k,4) = lf%d(i,j,k) * lf%v3(i,j,k)
    lf%u(i,j,k,5) = lf%b1(i,j,k)
    lf%u(i,j,k,6) = lf%b2(i,j,k)
    lf%u(i,j,k,7) = lf%b3(i,j,k)
    lf%u(i,j,k,8) = lf%e(i,j,k)
    lf%u(i,j,k,9) = lf%phi(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine amr_conserve

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine amr_minmod(mm,u)

  implicit none

  real*8,intent(in),dimension(3):: u
  real*8,intent(out):: mm
  real*8 x,y

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  x = u(3) - u(2) ; y = u(2) - u(1)

  mm = sign(1.d0,x) * max(0.d0,min(abs(x),sign(1.d0,x)*y))

return
end subroutine amr_minmod

end subroutine amr_numflux
