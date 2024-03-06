


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

  use settings,only:courant,outstyle,HGfac,hgcfl
  use grid
  use physval
  use gravmod,only:gravswitch,dtgrav,cgrav2
  use profiler_mod

  real(8),allocatable,dimension(:,:,:):: dti,dtg
  real(8):: cfmax0,cfmax
  integer:: i,j,k,n,jb,kb
  
!-------------------------------------------------------------------------

  call start_clock(wttim)

  allocate( dti(is:ie,js:je,ks:ke) )
  if(gravswitch==3)allocate(dtg,mold=dti)

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

! Compute dtgrav if using hyperbolic self-gravity
  if(gravswitch==3)then
   cgrav2 = HGfac*cfmax
!$omp parallel do private(i,j,k) collapse(3)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      call dtgrav_cell(i,j,k,dtg,cgrav2)
     end do
    end do
   end do
!$omp end parallel do
  end if

! Use longer time step if using nested grids
  if(fmr_max>0.and.crdnt==2)then
!$omp parallel
   do n = 1, fmr_max
    if(fmr_lvl(n)==0)cycle
    if(n==1)then
     jb=je;kb=ke
    else
     jb=min(2**(fmr_max-n+1),je) ; kb=min(2**(fmr_max-n+1),ke)
    end if
!$omp do private(i,j,k) collapse(3)
    do k = ks, ke, kb
     do j = js, je, jb
      do i = is+sum(fmr_lvl(0:n-1)), is+sum(fmr_lvl(0:n))-1
       call dti_cell(i,j,k,dti,jb=jb,kb=kb)
       if(gravswitch==3)call dtgrav_cell(i,j,k,dtg,cgrav2,jb=jb,kb=kb)
      end do
     end do
    end do
!$omp end do
   end do
!$omp end parallel
  end if

  dt = courant * minval(dti)

  if(outstyle==1) dt = min(dt,t_out-time)

  ch = cfmax

  if(gravswitch==3)then
   dtgrav = hgcfl * minval(dtg)
   if(outstyle==1) dtgrav = min(dtgrav,t_out-time)
   cgrav2 = cgrav2**2
  end if

  call stop_clock(wttim)
  
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
  cf1 = get_cf(d(i,j,k),cs(i,j,k),b1(i,j,k),b2(i,j,k),b3(i,j,k))+abs(v1(i,j,k)) 
  cf2 = get_cf(d(i,j,k),cs(i,j,k),b2(i,j,k),b1(i,j,k),b3(i,j,k))+abs(v2(i,j,k)) 
  cf3 = get_cf(d(i,j,k),cs(i,j,k),b3(i,j,k),b1(i,j,k),b2(i,j,k))+abs(v3(i,j,k)) 
 else
  cf1 = cs(i,j,k)+abs(v1(i,j,k))
  cf2 = cs(i,j,k)+abs(v2(i,j,k))
  cf3 = cs(i,j,k)+abs(v3(i,j,k))
 end if

 dti(i,j:j+jn,k:k+kn) = sum(dvol(i,j:j+jn,k:k+kn)) / &
                      ( cf1*off(is,ie) * sum(sa1(i-1:i   ,j:j+jn,k:k+kn)) &
                      + cf2*off(js,je) * sum(sa2(i,j-1:j+jn:jn+1,k:k+kn)) &
                      + cf3*off(ks,ke) * sum(sa3(i,j:j+jn,k-1:k+kn:kn+1)) )

 if(present(cfmax))cfmax = max(cf1,cf2,cf3)

return
end subroutine dti_cell

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE DTGRAV_CELL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gravity time step for a cell smeared over j-k coordinates

subroutine dtgrav_cell(i,j,k,dtg,cgrav,jb,kb)

 use grid,only:is,ie,js,je,ks,ke,sa1,sa2,sa3,dvol

 integer,intent(in):: i,j,k
 real(8),allocatable,intent(inout)::dtg(:,:,:)
 real(8),intent(in):: cgrav
 integer,intent(in),optional:: jb,kb
 integer:: jn,kn

!-----------------------------------------------------------------------------

 jn=0;kn=0
 if(present(jb)) jn = min(jb-1,je-js)
 if(present(kb)) kn = min(kb-1,ke-ks)

 dtg(i,j:j+jn,k:k+kn) = sum(dvol(i,j:j+jn,k:k+kn)) / cgrav / &
                      ( off(is,ie) * sum(sa1(i-1:i   ,j:j+jn,k:k+kn)) &
                      + off(js,je) * sum(sa2(i,j-1:j+jn:jn+1,k:k+kn)) &
                      + off(ks,ke) * sum(sa3(i,j:j+jn,k-1:k+kn:kn+1)) )

return
end subroutine dtgrav_cell
 
end module timestep_mod
