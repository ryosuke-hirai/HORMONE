!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE NUMFLUX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calcaulte numflux

subroutine numflux

  use settings,only:compswitch,spn,eq_sym,eostype,mag_on
  use grid
  use physval
  use hllflux_mod
  use pressure_mod

  implicit none

  real*8:: cfl, cfr, v1l, v1r, dl, dr, ptl, ptr, el, er, Tl, Tr, imul, imur, &
           b1l=0., b1r=0., b2l=0., b2r=0., b3l=0., b3r=0., phil=0., phir=0., &
           v2l, v2r, v3l, v3r, fix
  real*8,dimension(1:9)::tmpflux
  real*8,dimension(1:2):: dx
  real*8,dimension(1:spn):: spcl,spcr
  real*8 signdflx,ul,ur,fl,fr,rinji, rotfac

!--------------------------------------------------------------------

! calculate flux

! use MUSCL interpolation
  if(eostype==1)call meanmolweight
  call interpolation

! flux1 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!$omp parallel
!$omp do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,ufn,tmpflux,dx,Tl,Tr,&
!$omp imul,imur,fix,spcl,spcr,signdflx,n,ul,ur,fl,fr,rinji)
  do k = ks,ke
   do j = js,je
    do i = is-1, ie
     dx(1) = xi1(i)-x1(i) ; dx(2) = x1(i+1)-xi1(i)

     dl = d(i  ,j,k) + dx(1) * dd(i  ,j,k,1)
     dr = d(i+1,j,k) - dx(2) * dd(i+1,j,k,1)

     el = e(i  ,j,k) + dx(1) * de(i  ,j,k,1)
     er = e(i+1,j,k) - dx(2) * de(i+1,j,k,1)

     v1l = ( d(i  ,j,k)*v1(i  ,j,k) + dx(1)*dm1(i  ,j,k,1) ) / dl
     v1r = ( d(i+1,j,k)*v1(i+1,j,k) - dx(2)*dm1(i+1,j,k,1) ) / dr

     v2l = ( d(i  ,j,k)*v2(i  ,j,k) + dx(1)*dm2(i  ,j,k,1) ) / dl
     v2r = ( d(i+1,j,k)*v2(i+1,j,k) - dx(2)*dm2(i+1,j,k,1) ) / dr

     v3l = ( d(i  ,j,k)*v3(i  ,j,k) + dx(1)*dm3(i  ,j,k,1) ) / dl
     v3r = ( d(i+1,j,k)*v3(i+1,j,k) - dx(2)*dm3(i+1,j,k,1) ) / dr

     if(mag_on)then
      b1l = b1(i  ,j,k) + dx(1) * db1(i  ,j,k,1)
      b1r = b1(i+1,j,k) - dx(2) * db1(i+1,j,k,1)

      b2l = b2(i  ,j,k) + dx(1) * db2(i  ,j,k,1)
      b2r = b2(i+1,j,k) - dx(2) * db2(i+1,j,k,1)

      b3l = b3(i  ,j,k) + dx(1) * db3(i  ,j,k,1)
      b3r = b3(i+1,j,k) - dx(2) * db3(i+1,j,k,1)

      phil = phi(i  ,j,k) + dx(1) * dphi(i  ,j,k,1)
      phir = phi(i+1,j,k) - dx(2) * dphi(i+1,j,k,1)
     else
      b1l = 0d0 ; b1r = 0d0 ; b2l = 0d0 ; b2r = 0d0 ; b3l = 0d0 ; b3r = 0d0
      phil = 0d0 ; phir = 0d0
     end if

     Tl = T(i,j,k) ; Tr = T(i+1,j,k)

     select case (eostype)
     case(0:1) ! without recombination
      imul = 1d0/imu(i  ,j,k) + dx(1) * dmu(i  ,j,k,1)
      imur = 1d0/imu(i+1,j,k) - dx(2) * dmu(i+1,j,k,1)
      imul = 1d0/imul ; imur = 1d0/imur
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl)
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr)
     case(2) ! with recombination
      spcl(1:2) = spc(1:2,i  ,j,k) + dx(1) * dspc(1:2,i  ,j,k,1)
      spcr(1:2) = spc(1:2,i+1,j,k) - dx(2) * dspc(1:2,i+1,j,k,1)
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl, &
                    spcl(1),spcl(2))
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr, &
                    spcr(1),spcr(2))
     end select

     call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

     if(maxval(shock(i:i+1,j-1:j+1,k-1:k+1))==1)then ! if supersonic
      if(max(abs(v3l),abs(v3r))>1d-1*max(abs(v2l),abs(v2r),abs(v1l),abs(v1r)))then ! and aligned

       fl = dl*v1l;fr = dr*v1r ; ul = dl ; ur = dr
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(1) = rinji

       fl = dl*v1l*v1l+ptl-b1l*b1l;fr = dr*v1r*v1r+ptr-b1r*b1r 
       ul = dl*v1l ; ur = dr*v1r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(2) = rinji

       fl = dl*v1l*v2l-b1l*b2l;fr = dr*v1r*v2r-b1r*b2r 
       ul = dl*v2l ; ur = dr*v2r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(3) = rinji

       fl = dl*v1l*v3l-b1l*b3l;fr = dr*v1r*v3r-b1r*b3r 
       ul = dl*v3l ; ur = dr*v3r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(4) = rinji

       fl = (el+ptl)*v1l-b1l*(v1l*b1l+v2l*b2l+v3l*b3l)
       fr = (er+ptr)*v1r-b1r*(v1r*b1r+v2r*b2r+v3r*b3r)
       ul = el ; ur = er
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(8) = rinji

       if(mag_on)then
        fl = phil;fr = phir ; ul = b1l ; ur = b1r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(5) = rinji

        fl = b2l*v1l-b1l*v2l;fr = b2r*v1r-b1r*v2r
        ul = b2l ; ur = b2r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(6) = rinji

        fl = b3l*v1l-b1l*v3l;fr = b3r*v1r-b1r*v3r
        ul = b3l ; ur = b3r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(7) = rinji

        fl = ch*ch*b1l ; fr = ch*ch*b1r ; ul = phil;ur=phir
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(9) = rinji
       end if
       
      end if
     end if


     do ufn = 1,9
      flux1(i,j,k,ufn) = tmpflux(ufn)
     end do

     if(compswitch>=2)then
      do n = 1, spn
       spcl(n) = spc(n,i  ,j,k) + dx(1) * dspc(n,i  ,j,k,1)
       spcr(n) = spc(n,i+1,j,k) - dx(2) * dspc(n,i+1,j,k,1)
      end do
      signdflx = sign(0.5d0,flux1(i,j,k,1))
      fix = 1d0/( sum(spcl(1:spn))*(0.5d0+signdflx) + &
                  sum(spcr(1:spn))*(0.5d0-signdflx) )
      do n = 1, spn
       spcflx(n,i,j,k,1) = fix * flux1(i,j,k,1) &
              * ( spcl(n)*(0.5d0+signdflx) + spcr(n)*(0.5d0-signdflx) )
      end do
     end if

    end do
   end do
  end do
!$omp end do

! flux2 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if(je/=1)then
!$omp do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,ufn,tmpflux,dx,Tl,Tr,&
!$omp imul,imur,fix,spcl,spcr,signdflx,n)
  do k = ks,ke
   do j = js-1,je
    do i = is, ie

     dx(1) = xi2(j)-x2(j) ; dx(2) = x2(j+1)-xi2(j)

     dl = d(i,j  ,k) + dx(1) * dd(i,j  ,k,2)
     dr = d(i,j+1,k) - dx(2) * dd(i,j+1,k,2)

     el = e(i,j  ,k) + dx(1) * de(i,j  ,k,2)
     er = e(i,j+1,k) - dx(2) * de(i,j+1,k,2)

     v1l = ( d(i,j  ,k)*v2(i,j  ,k) + dx(1)*dm2(i,j  ,k,2) ) / dl
     v1r = ( d(i,j+1,k)*v2(i,j+1,k) - dx(2)*dm2(i,j+1,k,2) ) / dr

     v2l = ( d(i,j  ,k)*v3(i,j  ,k) + dx(1)*dm3(i,j  ,k,2) ) / dl
     v2r = ( d(i,j+1,k)*v3(i,j+1,k) - dx(2)*dm3(i,j+1,k,2) ) / dr

     v3l = ( d(i,j  ,k)*v1(i,j  ,k) + dx(1)*dm1(i,j  ,k,2) ) / dl
     v3r = ( d(i,j+1,k)*v1(i,j+1,k) - dx(2)*dm1(i,j+1,k,2) ) / dr

     if(mag_on)then
      b1l = b1(i,j  ,k) + dx(1) * db1(i,j  ,k,2)
      b1r = b1(i,j+1,k) - dx(2) * db1(i,j+1,k,2)

      b2l = b2(i,j  ,k) + dx(1) * db2(i,j  ,k,2)
      b2r = b2(i,j+1,k) - dx(2) * db2(i,j+1,k,2)

      b3l = b3(i,j  ,k) + dx(1) * db3(i,j  ,k,2)
      b3r = b3(i,j+1,k) - dx(2) * db3(i,j+1,k,2)

      phil = phi(i,j  ,k) + dx(1) * dphi(i,j  ,k,2)
      phir = phi(i,j+1,k) - dx(2) * dphi(i,j+1,k,2)
     else
      b1l = 0d0 ; b1r = 0d0 ; b2l = 0d0 ; b2r = 0d0 ; b3l = 0d0 ; b3r = 0d0
      phil = 0d0 ; phir = 0d0
     end if

     Tl = T(i,j,k) ; Tr = T(i,j+1,k)

     select case (eostype)
     case(0:1) ! without recombination
      imul = 1d0/imu(i,j  ,k) + dx(1) * dmu(i,j  ,k,2)
      imur = 1d0/imu(i,j+1,k) - dx(2) * dmu(i,j+1,k,2)
      imul = 1d0/imul ; imur = 1d0/imur
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl)
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr)
     case(2) ! with recombination
      spcl(1:2) = spc(1:2,i,j  ,k) + dx(1) * dspc(1:2,i,j  ,k,2)
      spcr(1:2) = spc(1:2,i,j+1,k) - dx(2) * dspc(1:2,i,j+1,k,2)
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl, &
                    spcl(1),spcl(2))
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr, &
                    spcr(1),spcr(2))
     end select
     
     call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)


     if(maxval(shock(i-1:i+1,j:j+1,k-1:k+1))==1)then ! if shock
      if(max(abs(v3l),abs(v3r))>1d-1*max(abs(v1l),abs(v1r)))then ! and aligned

       fl = dl*v1l;fr = dr*v1r ; ul = dl ; ur = dr
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(1) = rinji

       fl = dl*v1l*v1l+ptl-b1l*b1l;fr = dr*v1r*v1r+ptr-b1r*b1r 
       ul = dl*v1l ; ur = dr*v1r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(2) = rinji

       fl = dl*v1l*v2l-b1l*b2l;fr = dr*v1r*v2r-b1r*b2r 
       ul = dl*v2l ; ur = dr*v2r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(3) = rinji

       fl = dl*v1l*v3l-b1l*b3l;fr = dr*v1r*v3r-b1r*b3r 
       ul = dl*v3l ; ur = dr*v3r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(4) = rinji

       fl = (el+ptl)*v1l-b1l*(v1l*b1l+v2l*b2l+v3l*b3l)
       fr = (er+ptr)*v1r-b1r*(v1r*b1r+v2r*b2r+v3r*b3r)
       ul = el ; ur = er
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(8) = rinji

       if(mag_on)then
        fl = phil;fr = phir ; ul = b1l ; ur = b1r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(5) = rinji

        fl = b2l*v1l-b1l*v2l;fr = b2r*v1r-b1r*v2r
        ul = b2l ; ur = b2r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(6) = rinji
 
        fl = b3l*v1l-b1l*v3l;fr = b3r*v1r-b1r*v3r
        ul = b3l ; ur = b3r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(7) = rinji

        fl = ch*ch*b1l ; fr = ch*ch*b1r ; ul = phil; ur=phir
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(9) = rinji
       end if

      end if
     end if

     flux2(i,j,k,1) = tmpflux(1)
     flux2(i,j,k,2) = tmpflux(4)
     flux2(i,j,k,3) = tmpflux(2)
     flux2(i,j,k,4) = tmpflux(3)
     flux2(i,j,k,5) = tmpflux(7)
     flux2(i,j,k,6) = tmpflux(5)
     flux2(i,j,k,7) = tmpflux(6)
     flux2(i,j,k,8) = tmpflux(8)
     flux2(i,j,k,9) = tmpflux(9)

     if(compswitch>=2)then
      do n = 1, spn
       spcl(n) = spc(n,i,j  ,k) + dx(1) * dspc(n,i,j  ,k,2)
       spcr(n) = spc(n,i,j+1,k) - dx(2) * dspc(n,i,j+1,k,2)
      end do
      signdflx = sign(0.5d0,flux2(i,j,k,1))
      fix = 1d0/( sum(spcl(1:spn))*(0.5d0+signdflx) + &
                  sum(spcr(1:spn))*(0.5d0-signdflx) )
      do n = 1, spn
       spcflx(n,i,j,k,2) = fix * flux2(i,j,k,1) &
              * ( spcl(n)*(0.5d0+signdflx) + spcr(n)*(0.5d0-signdflx) )
      end do
     end if

    end do
   end do
  end do
!$omp end do
 end if

! flux3 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 if(ke/=1)then
!$omp do private(i,j,k,ptl,ptr,dl,dr,el,er,v1l,v1r,v2l,v2r,v3l,v3r,&
!$omp b1l,b1r,b2l,b2r,b3l,b3r,cfl,cfr,phil,phir,ufn,tmpflux,dx,Tl,Tr,&
!$omp imul,imur,fix,spcl,spcr,signdflx,n,ul,ur,fl,fr,rinji)
  do k = ks-1,ke
   do j = js,je
    do i = is, ie

     dx(1) = xi3(k)-x3(k) ; dx(2) = x3(k+1)-xi3(k)

     dl = d(i,j,k  ) + dx(1) * dd(i,j,k  ,3)
     dr = d(i,j,k+1) - dx(2) * dd(i,j,k+1,3)

     el = e(i,j,k  ) + dx(1) * de(i,j,k  ,3)
     er = e(i,j,k+1) - dx(2) * de(i,j,k+1,3)

     v1l = ( d(i,j,k  )*v3(i,j,k  ) + dx(1)*dm3(i,j,k  ,3) ) / dl
     v1r = ( d(i,j,k+1)*v3(i,j,k+1) - dx(2)*dm3(i,j,k+1,3) ) / dr

     v2l = ( d(i,j,k  )*v1(i,j,k  ) + dx(1)*dm1(i,j,k  ,3) ) / dl
     v2r = ( d(i,j,k+1)*v1(i,j,k+1) - dx(2)*dm1(i,j,k+1,3) ) / dr

     v3l = ( d(i,j,k  )*v2(i,j,k  ) + dx(1)*dm2(i,j,k  ,3) ) / dl
     v3r = ( d(i,j,k+1)*v2(i,j,k+1) - dx(2)*dm2(i,j,k+1,3) ) / dr

     if(mag_on)then
      b1l = b1(i,j,k  ) + dx(1) * db1(i,j,k  ,3)
      b1r = b1(i,j,k+1) - dx(2) * db1(i,j,k+1,3)

      b2l = b2(i,j,k  ) + dx(1) * db2(i,j,k  ,3)
      b2r = b2(i,j,k+1) - dx(2) * db2(i,j,k+1,3)

      b3l = b3(i,j,k  ) + dx(1) * db3(i,j,k  ,3)
      b3r = b3(i,j,k+1) - dx(2) * db3(i,j,k+1,3)

      phil = phi(i,j,k  ) + dx(1) * dphi(i,j,k  ,3)
      phir = phi(i,j,k+1) - dx(2) * dphi(i,j,k+1,3)
     else
      b1l = 0d0 ; b1r = 0d0 ; b2l = 0d0 ; b2r = 0d0 ; b3l = 0d0 ; b3r = 0d0
      phil = 0d0 ; phir = 0d0
     end if

     Tl = T(i,j,k) ; Tr = T(i,j,k+1)

     select case (eostype)
     case(0:1) ! without recombination
      imul = 1d0/imu(i,j,k  ) + dx(1) * dmu(i,j,k  ,3)
      imur = 1d0/imu(i,j,k+1) - dx(2) * dmu(i,j,k+1,3)
      imul = 1d0/imul ; imur = 1d0/imur
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl)
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr)
     case(2) ! with recombination
      spcl(1:2) = spc(1:2,i,j,k  ) + dx(1) * dspc(1:2,i,j,k  ,3)
      spcr(1:2) = spc(1:2,i,j,k+1) - dx(2) * dspc(1:2,i,j,k+1,3)
      call eos_p_cf(dl,v1l,v2l,v3l,b1l,b2l,b3l,el,Tl,imul,ptl,cfl, &
                    spcl(1),spcl(2))
      call eos_p_cf(dr,v1r,v2r,v3r,b1r,b2r,b3r,er,Tr,imur,ptr,cfr, &
                    spcr(1),spcr(2))
     end select

     call hlldflux(tmpflux,cfl,cfr,v1l,v1r,v2l,v2r,v3l,v3r,dl,dr, &
          el,er,ptl,ptr,b1l,b1r,b2l,b2r,b3l,b3r,phil,phir)

     if(maxval(shock(i-1:i+1,j,k:k+1))==1)then ! if supersonic
      if(max(abs(v2l),abs(v2r))>1d-1*max(abs(v1l),abs(v1r),abs(v3l),abs(v3r)))then ! and aligned

       fl = dl*v1l;fr = dr*v1r ; ul = dl ; ur = dr
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(1) = rinji

       fl = dl*v1l*v1l+ptl-b1l*b1l;fr = dr*v1r*v1r+ptr-b1r*b1r 
       ul = dl*v1l ; ur = dr*v1r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(2) = rinji

       fl = dl*v1l*v2l-b1l*b2l;fr = dr*v1r*v2r-b1r*b2r 
       ul = dl*v2l ; ur = dr*v2r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(3) = rinji

       fl = dl*v1l*v3l-b1l*b3l;fr = dr*v1r*v3r-b1r*b3r 
       ul = dl*v3l ; ur = dr*v3r
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(4) = rinji

       fl = (el+ptl)*v1l-b1l*(v1l*b1l+v2l*b2l+v3l*b3l)
       fr = (er+ptr)*v1r-b1r*(v1r*b1r+v2r*b2r+v3r*b3r)
       ul = el ; ur = er
       call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
       tmpflux(8) = rinji

       if(mag_on)then
        fl = phil;fr = phir ; ul = b1l ; ur = b1r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(5) = rinji

        fl = b2l*v1l-b1l*v2l;fr = b2r*v1r-b1r*v2r
        ul = b2l ; ur = b2r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(6) = rinji

        fl = b3l*v1l-b1l*v3l;fr = b3r*v1r-b1r*v3r
        ul = b3l ; ur = b3r
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(7) = rinji

        fl = ch*ch*b1l ; fr = ch*ch*b1r ; ul = phil;ur=phir
        call hllflux(rinji,fl,fr,ul,ur,cfl,cfr,v1l,v1r)
        tmpflux(9) = rinji
       end if

      end if
     end if

     flux3(i,j,k,1) = tmpflux(1)
     flux3(i,j,k,2) = tmpflux(3)
     flux3(i,j,k,3) = tmpflux(4)
     flux3(i,j,k,4) = tmpflux(2)
     flux3(i,j,k,5) = tmpflux(6)
     flux3(i,j,k,6) = tmpflux(7)
     flux3(i,j,k,7) = tmpflux(5)
     flux3(i,j,k,8) = tmpflux(8)
     flux3(i,j,k,9) = tmpflux(9)

     if(compswitch>=2)then
      do n = 1, spn
       spcl(n) = spc(n,i,j,k  ) + dx(1) * dspc(n,i,j,k  ,3)
       spcr(n) = spc(n,i,j,k+1) - dx(2) * dspc(n,i,j,k+1,3)
      end do
      signdflx = sign(0.5d0,flux3(i,j,k,1))
      fix = 1d0/( sum(spcl(1:spn))*(0.5d0+signdflx) + &
                  sum(spcr(1:spn))*(0.5d0-signdflx) )
      do n = 1, spn
       spcflx(n,i,j,k,3) = fix * flux3(i,j,k,1) &
              * ( spcl(n)*(0.5d0+signdflx) + spcr(n)*(0.5d0-signdflx) )
      end do
     end if

   end do
  end do
 end do
!$omp end do
end if

if(ie==1)flux1=0d0
if(ie==1.and.compswitch>=2)spcflx(1:spn,is-1:ie,js:je,ks:ke,1)=0d0
if(je<=2.and.crdnt==2)flux2=0d0
if(je==1.and.compswitch>=2)spcflx(1:spn,is:ie,js-1:je,ks:ke,2)=0d0
if(ke==1)flux3=0d0
if(ke==1.and.compswitch>=2)spcflx(1:spn,is:ie,js:je,ks-1:ke,3)=0d0

if(eq_sym.and.crdnt==1)flux3(is:ie,js:je,ks-1,1:9)=0d0

!$omp end parallel
return

end subroutine numflux
