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

  use settings,only:courant,outstyle,wtime,itim
  use grid
  use physval
  use omp_lib

  real(8),allocatable:: dti(:,:,:)
  real(8):: cfmax0,cfmax
  integer:: ierr, jb,kb
  
!-------------------------------------------------------------------------

  wtime(itim) = wtime(itim) - omp_get_wtime()

  allocate( dti(is:ie,js:je,ks:ke) )

  cfmax = 0d0
!$omp parallel do private(i,j,k,cfmax0) reduction(max:cfmax) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     call dti_cell(i,j,k,dti,cfmax=cfmax0)
     cfmax = max(cfmax,cfmax0)
    end do
   end do
  end do
!$omp end parallel do

  if(sphrn>0.and.crdnt==2)then
!$omp parallel
! Spherical symmetry for innermost few cells
   jb=je;kb=ke
!$omp do private(i,j,k)
   do i = is, is+sphrn-1
    j=js;k=ks
    call dti_cell(i,j,k,dti,jb=jb,kb=kb)
   end do
!$omp end do

! Smear out over jb x kb cells for inner cells
   jb = 16 ; kb = 16
!$omp do private(i,j,k) collapse(3)
   do i = is+sphrn, is+sphrn+trnsn16-1
    do k = ks, ke, kb
     do j = js, je, jb
      call dti_cell(i,j,k,dti,jb=jb,kb=kb)
     end do
    end do
   end do
!$omp end do
   jb = 8 ; kb = 8
!$omp do private(i,j,k) collapse(3)
   do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
    do k = ks, ke, kb
     do j = js, je, jb
      call dti_cell(i,j,k,dti,jb=jb,kb=kb)
     end do
    end do
   end do
!$omp end do
   jb = 4 ; kb = 4
!$omp do private(i,j,k) collapse(3)
   do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
    do k = ks, ke, kb
     do j = js, je, jb
      call dti_cell(i,j,k,dti,jb=jb,kb=kb)
     end do
    end do
   end do
!$omp end do
   jb = 2 ; kb = 2
!$omp do private(i,j,k) collapse(3)
   do i = is+sphrn+trnsn16+trnsn8+trnsn4, is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
    do k = ks, ke, kb
     do j = js, je, jb
      call dti_cell(i,j,k,dti,jb=jb,kb=kb)
     end do
    end do
   end do
!$omp end do
!$omp end parallel
  end if

!  dt = minval( dtdist(is:ie,js:je,ks:ke,1:3) )
  dt = minval( dti(is:ie,js:je,ks:ke) )

!  cfmax = maxval(abs(cs))

  dt = courant * dt

  if(outstyle==1) dt = min(dt,t_out-time)

  ch = cfmax

  wtime(itim) = wtime(itim) + omp_get_wtime()
  
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

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE DTI_CELL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute time step for a cell smeared over j-k coordinates

subroutine dti_cell(i,j,k,dti,jb,kb,cfmax)

 use settings,only:mag_on,eostype
 use grid,only:is,ie,js,je,ks,ke,sa1,sa2,sa3,dvol
 use physval,only:d,eint,T,imu,p,cs,v1,v2,v3,b1,b2,b3,spc
 use pressure_mod,only:eos_p_cs,get_cf

 integer,intent(in):: i,j,k
 real(8),allocatable,intent(inout)::dti(:,:,:)
 integer,intent(in),optional:: jb,kb
 real(8),intent(out),optional:: cfmax
 integer:: jn,kn,ierr
 real(8):: cf1,cf2,cf3
!-----------------------------------------------------------------------------

 jn=0;kn=0
 if(present(jb)) jn = min(jb-1,je-js)
 if(present(kb)) kn = min(kb-1,ke-ks)

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
  cf3 = get_cf(d(i,j,k),cs(i,j,k),b3(i,j,k),b1(i,j,k),b2(i,j,k))
 else
  cf1 = cs(i,j,k)
  cf2 = cs(i,j,k)
  cf3 = cs(i,j,k)
 end if

 dti(i,j:j+jn,k:k+kn) = sum(dvol(i,j:j+jn,k:k+kn)) / &
  ( (cf1 + abs(v1(i,j,k)) )*off(is,ie) * sum(sa1(i-1:i,j:j+jn,k:k+kn)) &
  + (cf2 + abs(v2(i,j,k)) )*off(js,je) * sum(sa2(i,j-1:j+jn:jn+1,k:k+kn)) &
  + (cf3 + abs(v3(i,j,k)) )*off(ks,ke) * sum(sa3(i,j:j+jn,k-1:k+kn:kn+1)) )

 if(present(cfmax))cfmax = max(cf1,cf2,cf3)

return
end subroutine dti_cell
 
end module timestep_mod
