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
  use composition_mod

  implicit none

  real*8 v0x,v0z,v0p

!-----------------------------------------------------------------------

  runge_type: select case (rktype)
  case(3) runge_type
   
   rk3_number: select case (rungen)
   case(1) rk3_number
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
   case(2) rk3_number
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
   case default rk3_number
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

   end select rk3_number

  case(2) runge_type

   rk2_number: select case (rungen)
   case(1) rk2_number
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

   case default rk2_number
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
   end select rk2_number

  case(1) runge_type

   call euler

  case default runge_type
   print *,"Error from rktype",rktype
   stop
  end select runge_type


! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.
  if(sphrn>0.and.crdnt==2)then
   do ufn = 1, 9
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

  call primitive

 
return
end subroutine rungekutta
