!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                    　  SUBROUTINE RUNGEKUTTA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate numflux using Runge-Kutta method

subroutine rungekutta

  use settings,only:rktype,compswitch,spn
  use grid
  use physval
  use merger_mod!temp
  use funcs! temp

  implicit none

  real*8 v0x,v0z,v0p

!-----------------------------------------------------------------------

  if(rktype==3)then

   if(rungen==1)then
!$omp parallel
!$omp do private (ufn,i,j,k)
    do ufn = 1,9

     do k = ks,ke
      do j = js,je
       do i = is,ie
        uorg(i,j,k,ufn) = u(i,j,k,ufn)
        u(i,j,k,ufn) = uorg(i,j,k,ufn) + dt * &
             ( idetg1(i) * &
               (detg1(i-1  )*flux1(i-1,j,k,ufn)-detg1(i  )*flux1(i,j,k,ufn)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*flux2(i,j-1,k,ufn)-detg2(i,j)*flux2(i,j,k,ufn)) &
             + idetg3(i,j,k) * (flux3(i,j,k-1,ufn)-flux3(i,j,k,ufn)) &
             + src(i,j,k,ufn) )
       end do
      end do
     end do

    end do
!$omp end do
    if(compswitch>=2)then
!$omp do private (i,j,k,n)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         spcorg(n,i,j,k) = spc(n,i,j,k)*uorg(i,j,k,1)
         spc(n,i,j,k) = (spcorg(n,i,j,k) + dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,1)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel
   elseif(rungen==2)then
!$omp parallel
!$omp do private(ufn,i,j,k)
    do ufn = 1,9
     
     do k = ks,ke
      do j = js,je
       do i = is,ie
        u(i,j,k,ufn) = 0.75d0 * uorg(i,j,k,ufn) + 0.25d0 * u(i,j,k,ufn) + &
             0.25d0 * dt * &
             ( idetg1(i) * &
               (detg1(i-1  )*flux1(i-1,j,k,ufn)-detg1(i  )*flux1(i,j,k,ufn)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*flux2(i,j-1,k,ufn)-detg2(i,j)*flux2(i,j,k,ufn)) &
             + idetg3(i,j,k) * (flux3(i,j,k-1,ufn)-flux3(i,j,k,ufn)) &
             + src(i,j,k,ufn) )
       end do
      end do
     end do

    end do
!$omp end do

    if(compswitch>=2)then
!$omp do private(i,j,k,n)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         spc(n,i,j,k) =( 0.75d0*spcorg(n,i,j,k) &
                      + 0.25d0*spc(n,i,j,k)*d(i,j,k) &
              + 0.25d0 * dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,1)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel
   else
!$omp parallel
!$omp do private(ufn,i,j,k)
    do ufn = 1,9

     do k = ks,ke
      do j = js,je
       do i = is,ie
        u(i,j,k,ufn) = 1.d0/3.d0 * uorg(i,j,k,ufn) + 2.d0/3.d0 * u(i,j,k,ufn) +&
             2.d0/3.d0 * dt * &
             ( idetg1(i) * &
               (detg1(i-1  )*flux1(i-1,j,k,ufn)-detg1(i  )*flux1(i,j,k,ufn)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*flux2(i,j-1,k,ufn)-detg2(i,j)*flux2(i,j,k,ufn)) &
             + idetg3(i,j,k) * (flux3(i,j,k-1,ufn)-flux3(i,j,k,ufn)) &
             + src(i,j,k,ufn) )
       end do
      end do
     end do

    end do
!$omp end do

    if(compswitch>=2)then
!$omp do private(i,j,k,n)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         spc(n,i,j,k) = (1d0/3d0*spcorg(n,i,j,k) &
                      + 2d0/3d0*spc(n,i,j,k)*d(i,j,k) &
              + 2d0/3d0 * dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,1)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel

   end if

  else if(rktype==2)then

   if(rungen==1)then
!$omp parallel
!$omp do private(ufn,i,j,k)
    do ufn = 1,9

     do k = ks,ke
      do j = js,je
       do i = is,ie
        uorg(i,j,k,ufn) = u(i,j,k,ufn)
        u(i,j,k,ufn) = uorg(i,j,k,ufn) + dt * &
             ( idetg1(i) * &
               (detg1(i-1  )*flux1(i-1,j,k,ufn)-detg1(i  )*flux1(i,j,k,ufn)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*flux2(i,j-1,k,ufn)-detg2(i,j)*flux2(i,j,k,ufn)) &
             + idetg3(i,j,k) * (flux3(i,j,k-1,ufn)-flux3(i,j,k,ufn)) &
             + src(i,j,k,ufn) )
       end do
      end do
     end do

    end do
!$omp end do
    if(compswitch>=2)then
!$omp do private(i,j,k,n)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         spcorg(n,i,j,k) = spc(n,i,j,k)
         spc(n,i,j,k) = (spcorg(n,i,j,k)*uorg(i,j,k,1) + dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,1)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel

   else
!$omp parallel
!$omp do private(ufn,i,j,k)
    do ufn = 1,9
     
     do k = ks,ke
      do j = js,je
       do i = is,ie
        u(i,j,k,ufn) = 5.d-1 * uorg(i,j,k,ufn) + 5.d-1 * u(i,j,k,ufn) + &
             5.d-1 * dt * &
             ( idetg1(i) * &
               (detg1(i-1  )*flux1(i-1,j,k,ufn)-detg1(i  )*flux1(i,j,k,ufn)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*flux2(i,j-1,k,ufn)-detg2(i,j)*flux2(i,j,k,ufn)) &
             + idetg3(i,j,k) * (flux3(i,j,k-1,ufn)-flux3(i,j,k,ufn)) &
             + src(i,j,k,ufn) )
       end do
      end do
     end do

    end do
!$omp end do
    if(compswitch>=2)then
!$omp do private(i,j,k,n)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         spc(n,i,j,k) = (0.5d0*spcorg(n,i,j,k)*uorg(i,j,k,1) &
                      + 0.5d0*spc(n,i,j,k)*d(i,j,k) &
              + 0.5d0 * dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,1)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel
   end if

  elseif(rktype==1)then

   call euler

  else
   print *,"Error from rktype",rktype
   stop
  end if

!temp
!  u(is,js:je,ks:ke,2) = 0.2d0*(u(is+2,js:je,ks:ke,2)/u(is+2,js:je,ks:ke,1))*u(is,js:je,ks:ke,1)
!  u(is+1,js:je,ks:ke,2) = 0.6d0*(u(is+2,js:je,ks:ke,2)/u(is+2,js:je,ks:ke,1))*u(is+1,js:je,ks:ke,1)

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.
  if(sphrn>0.and.crdnt==2)then
   ! for composition
   if(compswitch>=2)then
    do i = is, is+2
     do n = 1, spn
      spc(n,i,js:je,ks:ke) = sum( u(i,js:je,ks:ke,1)*spc(n,i,js:je,ks:ke) &
                                 *dvol(i,js:je,ks:ke) ) &
                           / sum( u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke) )
     end do
    end do
    do i = is+3, sphrn
     do n = 1, spn
      do j = js, je,40
       spc(n,i,j:j+39,ks:ke) = sum( u(i,j:j+39,ks:ke,1)*spc(n,i,j:j+39,ks:ke) &
                                  *dvol(i,j:j+39,ks:ke) ) &
                             / sum( u(i,j:j+39,ks:ke,1)*dvol(i,j:j+39,ks:ke) )
      end do
     end do
    end do
    do i = sphrn+1, sphrn+trnsn1
     do n = 1, spn
      do j = js, je,8
       spc(n,i,j:j+7,ks:ke) = sum( u(i,j:j+7,ks:ke,1)*spc(n,i,j:j+7,ks:ke) &
                                  *dvol(i,j:j+7,ks:ke) ) &
                            / sum( u(i,j:j+7,ks:ke,1)*dvol(i,j:j+7,ks:ke) )
      end do
     end do
    end do
    do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
     do n = 1, spn
      do j = js, je,4
       spc(n,i,j:j+3,ks:ke) = sum( u(i,j:j+3,ks:ke,1)*spc(n,i,j:j+3,ks:ke) &
                                  *dvol(i,j:j+3,ks:ke) ) &
                            / sum( u(i,j:j+3,ks:ke,1)*dvol(i,j:j+3,ks:ke) )
      end do
     end do
    end do
    do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
     do n = 1, spn
      do j = js, je,2
       spc(n,i,j:j+1,ks:ke) = sum( u(i,j:j+1,ks:ke,1)*spc(n,i,j:j+1,ks:ke) &
                                  *dvol(i,j:j+1,ks:ke) ) &
                            / sum( u(i,j:j+1,ks:ke,1)*dvol(i,j:j+1,ks:ke) )
      end do
     end do
    end do
   end if
  end if
  ! for angular momentum
  k = ks
  do i = is, is+2
   v0p=0d0
   do j = js, je
    v0p = v0p + spin_coefft(j)*u(i,j,k,4)*dvol(i,j,k)
   end do
   v0p=v0p/sum(spin_coefft(js:je)*dvol(i,js:je,k)*spin_coefft(js:je))
   u(i,js:je,k,4) = spin_coefft(js:je)*v0p
  end do
  do i = is+3, sphrn
   do n = js, je, 40
    v0p = 0d0
    do j = n, n+39
     v0p = v0p + spin_coefft(j)*u(i,j,k,4)*dvol(i,j,k)
    end do
    v0p=v0p/sum(spin_coefft(n:n+39)*dvol(i,n:n+39,k)*spin_coefft(n:n+39))
    u(i,n:n+39,k,4) = spin_coefft(n:n+39) * v0p
   end do
  end do
  do i = sphrn+1, sphrn+trnsn1
   do n = js, je, 8
    v0p = 0d0
    do j = n, n+7
     v0p = v0p + spin_coefft(j)*u(i,j,k,4)*dvol(i,j,k)
    end do
    v0p=v0p/sum(spin_coefft(n:n+7)*dvol(i,n:n+7,k)*spin_coefft(n:n+7))
    u(i,n:n+7,k,4) = spin_coefft(n:n+7) * v0p
   end do
  end do
  do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
   do n = js, je, 4
    v0p = 0d0
    do j = n, n+3
     v0p = v0p + spin_coefft(j)*u(i,j,k,4)*dvol(i,j,k)
    end do
    v0p=v0p/sum(spin_coefft(n:n+3)*dvol(i,n:n+3,k)*spin_coefft(n:n+3))
    u(i,n:n+3,k,4) = spin_coefft(n:n+3) * v0p
   end do
  end do
  do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
   do n = js, je, 2
    v0p = 0d0
    do j = n, n+1
     v0p = v0p + spin_coefft(j)*u(i,j,k,4)*dvol(i,j,k)
    end do
    v0p=v0p/sum(spin_coefft(n:n+1)*dvol(i,n:n+1,k)*spin_coefft(n:n+1))
    u(i,n:n+1,k,4) = spin_coefft(n:n+1) * v0p
   end do
  end do

  ! for scalar quantities
  do ufn = 1, 9
   if(ufn==4)cycle
   do i = is, is+2
    u(i,js:je,ks:ke,ufn) = sum( u(i,js:je,ks:ke,ufn)*dvol(i,js:je,ks:ke) ) &
                         / sum( dvol(i,js:je,ks:ke) )
   end do
   do i = is+3, sphrn
    do j = js, je, 40
     u(i,j:j+39,ks:ke,ufn) = sum( u(i,j:j+39,ks:ke,ufn)*dvol(i,j:j+39,ks:ke) ) &
                           / sum( dvol(i,j:j+39,ks:ke) )
      
    end do
   end do
   do i = sphrn+1, sphrn+trnsn1
    do j = js, je,8
     u(i,j:j+7,ks:ke,ufn) = sum( u(i,j:j+7,ks:ke,ufn)*dvol(i,j:j+7,ks:ke) ) &
                          / sum( dvol(i,j:j+7,ks:ke) )
    end do
   end do
   do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
    do j = js, je,4
     u(i,j:j+3,ks:ke,ufn) = sum( u(i,j:j+3,ks:ke,ufn)*dvol(i,j:j+3,ks:ke) ) &
                          / sum( dvol(i,j:j+3,ks:ke) )
    end do
   end do
   do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
    do j = js, je,2
     u(i,j:j+1,ks:ke,ufn) = sum( u(i,j:j+1,ks:ke,ufn)*dvol(i,j:j+1,ks:ke) ) &
                          / sum( dvol(i,j:j+1,ks:ke) )
    end do
   end do
  end do
   ! for vector quantities
!!$   k = ks
!!$   do i = is, is+2
!!$    v0x=0d0;v0z=0d0;v0p=0d0
!!$    do j = js, je
!!$     v0x = v0x + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$              ( u(i,j,k,2)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j)))&
!!$              + u(i,j,k,3)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j))) )
!!$     v0z = v0z + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$              ( u(i,j,k,2)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j)))&
!!$              - u(i,j,k,3)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j))) )
!!$     v0p = v0p + u(i,j,k,4)/x1(i)*spin_coeffr(i)*spin_coefft(j)*dvol(i,j,k)
!!$    end do
!!$    v0x = v0x/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$    v0z = v0z/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$    v0p = v0p*16d0/(xi1(i)**4d0-xi1(i-1)**4d0)/(2d0*xi2(je)-2d0*xi2(js-1)-sin(2d0*xi2(je))+sin(2d0*xi2(js-1)))
!!$    do j = js, je
!!$     u(i,j,k,2) = u(i,j,k,1)*(v0x*sinc(j)+v0z*cosc(j))
!!$     u(i,j,k,3) = u(i,j,k,1)*(v0x*cosc(j)-v0z*sinc(j))
!!$!     u(i,j,k,4) = v0p * x1(i) *sinc(j)
!!$    end do
!!$   end do
!!$   do i = is+3, sphrn
!!$    do n = js,je,40
!!$     v0x=0d0;v0z=0d0
!!$     do j = n, n+39
!!$      v0x = v0x + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j)))&
!!$               + u(i,j,k,3)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j))) )
!!$      v0z = v0z + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j)))&
!!$               - u(i,j,k,3)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j))) )
!!$     end do
!!$     v0x = v0x/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     v0z = v0z/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     do j = n, n+7
!!$      u(i,j,k,2) = u(i,j,k,1)*(v0x*sinc(j)+v0z*cosc(j))
!!$      u(i,j,k,3) = u(i,j,k,1)*(v0x*cosc(j)-v0z*sinc(j))
!!$     end do
!!$    end do
!!$   end do
!!$   do i = sphrn+1, sphrn+trnsn1
!!$    do n = js,je,8
!!$     v0x=0d0;v0z=0d0
!!$     do j = n, n+7
!!$      v0x = v0x + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j)))&
!!$               + u(i,j,k,3)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j))) )
!!$      v0z = v0z + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j)))&
!!$               - u(i,j,k,3)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j))) )
!!$     end do
!!$     v0x = v0x/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     v0z = v0z/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     do j = n, n+7
!!$      u(i,j,k,2) = u(i,j,k,1)*(v0x*sinc(j)+v0z*cosc(j))
!!$      u(i,j,k,3) = u(i,j,k,1)*(v0x*cosc(j)-v0z*sinc(j))
!!$     end do
!!$    end do
!!$   end do
!!$   do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
!!$    do n = js,je,4
!!$     v0x=0d0;v0z=0d0
!!$     do j = n, n+3
!!$      v0x = v0x + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j)))&
!!$               + u(i,j,k,3)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j))) )
!!$      v0z = v0z + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j)))&
!!$               - u(i,j,k,3)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j))) )
!!$     end do
!!$     v0x = v0x/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     v0z = v0z/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     do j = n, n+3
!!$      u(i,j,k,2) = u(i,j,k,1)*(v0x*sinc(j)+v0z*cosc(j))
!!$      u(i,j,k,3) = u(i,j,k,1)*(v0x*cosc(j)-v0z*sinc(j))
!!$     end do
!!$    end do
!!$   end do
!!$   do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
!!$    do n = js,je,2
!!$     v0x=0d0;v0z=0d0
!!$     do j = n, n+1
!!$      v0x = v0x + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j)))&
!!$               + u(i,j,k,3)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j))) )
!!$      v0z = v0z + (pw(3,xi1(i))-pw(3,xi1(i-1)))/12d0*&
!!$               ( u(i,j,k,2)*(cos(2d0*xi2(j-1))-cos(2d0*xi2(j)))&
!!$               - u(i,j,k,3)*(2d0*dxi2(j)+sin(2d0*xi2(j-1))-sin(2d0*xi2(j))) )
!!$     end do
!!$     v0x = v0x/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     v0z = v0z/sum(u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke))
!!$     do j = n, n+1
!!$      u(i,j,k,2) = u(i,j,k,1)*(v0x*sinc(j)+v0z*cosc(j))
!!$      u(i,j,k,3) = u(i,j,k,1)*(v0x*cosc(j)-v0z*sinc(j))
!!$     end do
!!$    end do
!!$   end do

  call primitive

return
end subroutine rungekutta
