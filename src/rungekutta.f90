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
  use smear_mod

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
        u(i,j,k,ufn) = 0.5d0 * uorg(i,j,k,ufn) + 0.5d0 * u(i,j,k,ufn) + &
             0.5d0 * dt * &
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


  call smear

  call primitive

 
return
end subroutine rungekutta
