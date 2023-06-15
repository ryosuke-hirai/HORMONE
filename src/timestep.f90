module timestep_mod
 implicit none
 private:: off
 public:: timestep
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE TIMESTEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set dt.

 subroutine timestep

  use settings,only:courant,outstyle,eostype,mag_on
  use grid
  use physval
  use constants,only:huge
  use pressure_mod,only:eos_p_cs,get_cf

  real(8),allocatable:: dtdist(:,:,:,:), dti(:,:,:)
  real(8):: cfmax, cf1,cf2,cf3
  integer:: ierr
  
!-------------------------------------------------------------------------


  allocate( dtdist(is:ie,js:je,ks:ke,1:3), dti(is:ie,js:je,ks:ke) )
!!$  dtdist = huge
!!$!$omp parallel do private(i,j,k)
!!$  do k = ks, ke
!!$   do j = js, je
!!$    do i = is, ie
!!$     if(ie>1)then
!!$      dtdist(i,j,k,1) = dxi1(i) / ( abs(v1(i,j,k))+abs(v3(i,j,k))+cs(i,j,k) )
!!$     end if
!!$     if(je>1)then
!!$      dtdist(i,j,k,2) = g22(i)*dxi2(j) &
!!$                      / ( abs(v2(i,j,k))+abs(v3(i,j,k))+cs(i,j,k) )
!!$     end if
!!$     if(ke>1)then
!!$      dtdist(i,j,k,3) = g33(i,j)*dxi3(k) &
!!$                      / ( abs(v3(i,j,k)+cs(i,j,k) ) )
!!$     end if
!!$    end do
!!$   end do
!!$  end do
!!$!$omp end parallel do

  dti = huge
  cfmax = 0d0
!$omp parallel do private(i,j,k,cf1,cf2,cf3) reduction (max:cfmax) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     select case(eostype)
     case(0:1)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), ierr=ierr )
     case(2)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
     end select

     if(mag_on)then
      cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
      cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))
      cf3 = get_cf(d(i,j,k),cs(i,j,k),b3(i,j,k),b2(i,j,k),b1(i,j,k))
     else
      cf1 = cs(i,j,k)
      cf2 = cs(i,j,k)
      cf3 = cs(i,j,k)
     end if

     dti(i,j,k) = dvol(i,j,k) / &
                ( ( cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j,k)) &
                + ( cf2 + abs(v2(i,j,k)) )*off(js,je) * sum(sa2(i,j-1:j,k)) &
                + ( cf3 + abs(v3(i,j,k)) )*off(ks,ke) * sum(sa3(i,j,k-1:k)) )

     cfmax = max(cfmax,cf1,cf2,cf3)
    end do
   end do
  end do
!$omp end parallel do

! Temporary (for workaround mesh)
!!$  if(sphrn>0.and.crdnt==2)then
!!$   do i = is, is+2
!!$    dtdist(i,js:je,ks:ke,2) = sum( dtdist(i,js:je,ks:ke,2) )
!!$   end do
!!$   do i = is+3, sphrn
!!$    do j = js, je,40
!!$     dtdist(i,j:j+39,ks:ke,2) = sum( dtdist(i,j:j+39,ks:ke,2) )
!!$    end do
!!$   end do
!!$   do i = sphrn+1, sphrn+trnsn1
!!$    do j = js, je,8
!!$     dtdist(i,j:j+7,ks:ke,2) = sum( dtdist(i,j:j+7,ks:ke,2) )    
!!$    end do
!!$   end do
!!$   do i = sphrn+trnsn1+1, sphrn+trnsn1+trnsn2
!!$    do j = js, je,4
!!$     dtdist(i,j:j+3,ks:ke,2) = sum( dtdist(i,j:j+3,ks:ke,2) )
!!$    end do
!!$   end do
!!$   do i = sphrn+trnsn1+trnsn2+1, sphrn+trnsn1+trnsn2+trnsn3
!!$    do j = js, je,2
!!$     dtdist(i,j:j+1,ks:ke,2) = sum( dtdist(i,j:j+1,ks:ke,2) )
!!$    end do
!!$   end do
!!$  end if

  if(sphrn>0.and.crdnt==2.and.ke==ks)then
   do i = is, is+sphrn-1
    j=js;k=ks
    select case(eostype)
    case(0:1)
     call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                   p(i,j,k), cs(i,j,k), ierr=ierr )
    case(2)
     call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                   p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
    end select

    if(mag_on)then
     cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
    else
     cf1 = cs(i,j,k)
    end if
    dti(i,js:je,ks:ke) = sum(dvol(i,js:je,ks:ke)) / &
                ( ( cf1+abs(v1(i,j,k)) )*off(is,ie)*sum(sa1(i-1:i,js:je,k)))
   end do
   do i = is+sphrn, is+sphrn+trnsn16-1
    do j = js, je,16
     k=ks
     select case(eostype)
     case(0:1)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), ierr=ierr )
     case(2)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
     end select

     if(mag_on)then
      cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
      cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))
     else
      cf1 = cs(i,j,k)
      cf2 = cs(i,j,k)
     end if
     dti(i,j:j+15,ks:ke) = sum(dvol(i,j:j+15,ks:ke)) / &
      ( (cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j:j+15,k)) &
      + (cf2 + abs(v2(i,j,k)) )*off(js,je) * (sa2(i,j-1,k)+sa2(i,j+15,k))  )
    end do
   end do
   do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
    do j = js, je,8
     k=ks
     select case(eostype)
     case(0:1)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), ierr=ierr )
     case(2)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
     end select

     if(mag_on)then
      cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
      cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))
     else
      cf1 = cs(i,j,k)
      cf2 = cs(i,j,k)
     end if
     dti(i,j:j+7,ks:ke) = sum(dvol(i,j:j+7,ks:ke)) / &
      ( (cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j:j+7,k)) &
      + (cf2 + abs(v2(i,j,k)) )*off(js,je) * (sa2(i,j-1,k)+sa2(i,j+7,k))  )
    end do
   end do
   do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
    do j = js, je,4
     k=ks
     select case(eostype)
     case(0:1)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), ierr=ierr )
     case(2)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
     end select

     if(mag_on)then
      cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
      cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))
     else
      cf1 = cs(i,j,k)
      cf2 = cs(i,j,k)
     end if
     dti(i,j:j+3,ks:ke) = sum(dvol(i,j:j+3,ks:ke)) / &
      ( (cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j:j+3,k)) &
      + (cf2 + abs(v2(i,j,k)) )*off(js,je) * (sa2(i,j-1,k)+sa2(i,j+3,k))  )
    end do
   end do
   do i = is+sphrn+trnsn16+trnsn8+trnsn4, is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
    do j = js, je,2
     k=ks
     select case(eostype)
     case(0:1)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), ierr=ierr )
     case(2)
      call eos_p_cs(d(i,j,k), eint(i,j,k), T(i,j,k), imu(i,j,k), &
                    p(i,j,k), cs(i,j,k), spc(1,i,j,k), spc(2,i,j,k), ierr=ierr )
     end select

     if(mag_on)then
      cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))
      cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))
     else
      cf1 = cs(i,j,k)
      cf2 = cs(i,j,k)
     end if
     dti(i,j:j+1,ks:ke) = sum(dvol(i,j:j+1,ks:ke)) / &
      ( (cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j:j+1,k)) &
      + (cf2 + abs(v2(i,j,k)) )*off(js,je) * (sa2(i,j-1,k)+sa2(i,j+1,k))  )
    end do
   end do
  end if

!  dt = minval( dtdist(is:ie,js:je,ks:ke,1:3) )
  dt = minval( dti(is:ie,js:je,ks:ke) )

  cfmax = maxval(abs(cs))

  dt = courant * dt

  if(outstyle==1) dt = min(dt,t_out-time)

  ch = cfmax

 return
 end subroutine timestep

 pure function off(is,ie)
  integer,intent(in)::is,ie
  real(8):: off
  if(ie==is)then
   off = 0d0
  else
   off = 1d0
  end if
 end function off

end module timestep_mod
