module rungekutta_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                    　  SUBROUTINE RUNGEKUTTA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate numflux using Runge-Kutta method

 subroutine rungekutta

  use settings,only:rktype,compswitch,spn,wtime,irng
  use grid
  use physval
  use composition_mod
  use smear_mod
  use omp_lib

  real(8):: faco,fact,facn

!-----------------------------------------------------------------------

  wtime(irng) = wtime(irng) - omp_get_wtime()

  runge_type: select case (rktype)
  case(3) runge_type
   
   rk3_number: select case (rungen)
   case(1) rk3_number
    faco = 1d0 ; fact = 1d0 ; facn = 0d0
   case(2) rk3_number
    faco = 0.75d0 ; fact = 0.25d0 ; facn = fact
   case(3) rk3_number
    faco = 1d0/3d0 ; fact = 2d0/3d0 ; facn = fact
   end select rk3_number

  case(2) runge_type
   rk2_number: select case (rungen)
   case(1) rk2_number
    faco = 1d0 ; fact = 1d0 ; facn = 0d0
   case(2) rk2_number
    faco = 0.5d0 ; fact = 0.5d0 ; facn = fact
   end select rk2_number
   
  case(1) runge_type
   faco = 1d0 ; fact = 1d0 ; facn = 0d0

  case default
   print*,"Error from rktype: rktype =",rktype
   stop
   
  end select runge_type
    
!$omp parallel
!$omp do private (ufn,i,j,k) collapse(4)
    do ufn = 1,ufnmax

     do k = ks,ke
      do j = js,je
       do i = is,ie
        if(rungen==1)uorg(i,j,k,ufn) = u(i,j,k,ufn)
        u(i,j,k,ufn) = faco*uorg(i,j,k,ufn) &
                     + facn*u(i,j,k,ufn) &
             + fact*dt * &
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
!$omp do private (i,j,k,n) collapse(4)
     do k = ks, ke
      do j = js, je
       do i = is, ie
        do n = 1, spn
         if(rungen==1)spcorg(n,i,j,k) = spc(n,i,j,k)*uorg(i,j,k,icnt)
         spc(n,i,j,k) = (faco*spcorg(n,i,j,k) + facn*spcorg(n,i,j,k)+fact*dt * &
              ( idetg1(i) * &
                ( detg1(i-1)*spcflx(n,i-1,j,k,1)   &
                - detg1(i  )*spcflx(n,i  ,j,k,1) ) &
              + idetg2(i,j) * &
                ( detg2(i,j-1)*spcflx(n,i,j-1,k,2)   &
                - detg2(i,j  )*spcflx(n,i,j  ,k,2) ) &
              + idetg3(i,j,k) * &
                ( spcflx(n,i,j,k-1,3) &
                - spcflx(n,i,j,k  ,3) ) ) ) / u(i,j,k,icnt)
        end do
       end do
      end do
     end do
!$omp end do
    end if
!$omp end parallel

  call smear
  call primitive
  
  wtime(irng) = wtime(irng) + omp_get_wtime()

  return
 end subroutine rungekutta

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE EULER
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate numflux using Euler method

 subroutine euler

  use grid
  use physval

!--------------------------------------------------------------------

!$omp parallel do private(i,j,k,ufn) collapse(4)
  do ufn = 1,ufnmax
   do k = ks,ke
    do j = js,je
     do i = is,ie
      u(i,j,k,ufn) = u(i,j,k,ufn) - dt * &
           ( idetg1(i) * &
             (detg1(i  )*flux1(i,j,k,ufn)-detg1(i-1  )*flux1(i-1,j,k,ufn)) &
           + idetg2(i,j) * &
             (detg2(i,j)*flux2(i,j,k,ufn)-detg2(i,j-1)*flux2(i,j-1,k,ufn)) &
           + idetg3(i,j,k) * (flux3(i,j,k,ufn)-flux3(i,j,k-1,ufn)) &
           + src(i,j,k,ufn) )
     end do
    end do
   end do
  end do
!$omp end parallel do

  return
 end subroutine euler

 !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE PRIMITIVE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert conserved values to primitive values

 subroutine primitive

  use settings,only:mag_on
  use grid
  use physval
  use pressure_mod,only:pressure
  use composition_mod,only:meanmolweight
  
!--------------------------------------------------------------------
  
!$omp parallel do private(i,j,k) collapse(3)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     d(i,j,k)  = u(i,j,k,icnt)
     v1(i,j,k) = u(i,j,k,imo1) / d(i,j,k)
     v2(i,j,k) = u(i,j,k,imo2) / d(i,j,k)
     v3(i,j,k) = u(i,j,k,imo3) / d(i,j,k)
     e(i,j,k)  = u(i,j,k,iene)
     if(.not.mag_on)cycle
     b1(i,j,k) = u(i,j,k,img1)
     b2(i,j,k) = u(i,j,k,img2)
     b3(i,j,k) = u(i,j,k,img3)
     ! for 9 wave method
     phi(i,j,k)= u(i,j,k,idcl)
    end do
   end do
  end do
!$omp end parallel do

  call meanmolweight
  call pressure

  return
 end subroutine primitive

 
end module rungekutta_mod

