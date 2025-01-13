module timestep_mod
 implicit none
 private:: off, dti_cell
 public:: timestep, dtgrav_cell
contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE TIMESTEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set dt.

 subroutine timestep

  use mpi_utils,only:allreduce_mpi
  use mpi_domain,only:partially_my_domain
  use constants,only:tiny
  use settings,only:courant,outstyle,HGfac,hgcfl,maxtngrv
  use grid
  use physval
  use gravmod,only:gravswitch,dtgrav,cgrav,cgrav2,dtg_unit
  use smear_mod
  use profiler_mod

  integer:: i,j,k,l,n,jb,kb
  real(8),allocatable,dimension(:,:,:):: dti
  real(8):: cfmax0,cfmax

!-------------------------------------------------------------------------

  call start_clock(wttim)

  allocate( dti(is:ie,js:je,ks:ke) )

  cfmax = 0d0
  cgrav = tiny
!$omp parallel do private(i,j,k,cfmax0) reduction(max:cfmax,cgrav) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     call dti_cell(i,j,k,dti,cfmax=cfmax0)
     cfmax = max(cfmax,cfmax0)
     cgrav = max(cgrav,cfmax0)
    end do
   end do
  end do
!$omp end parallel do

  call allreduce_mpi('max',cfmax)
!  call allreduce_mpi('max',cgrav)

! Use longer time step if using nested grids
  if(fmr_max>0.and.crdnt==2)then
!$omp parallel do private(i,j,k,l,jb,kb) schedule(dynamic)
   do n = 1, nsmear
    l = lijk_from_id(0,n)
    if(fmr_lvl(l)==0)cycle
    i = lijk_from_id(1,n)
    j = lijk_from_id(2,n)
    k = lijk_from_id(3,n)
    jb = block_j(l)
    kb = block_k(l)
    if(partially_my_domain(i,j,k,1,jb,kb))call dti_cell(i,j,k,dti,jb=jb,kb=kb)
   end do
!$omp end parallel do
  end if

  dt = courant * minval(dti)

  if(outstyle==1)then
   if(t_out-time-dt>0d0.and.t_out-time-2d0*dt<0d0)then
    dt = 0.5d0*(t_out-time)
   else
    dt = min(dt,t_out-time)
   end if
  end if

  call allreduce_mpi('min',dt)

  ch = cfmax

! Compute dtgrav if using hyperbolic self-gravity
  if(gravswitch==3)then
   cgrav = cfmax
   cgrav = HGfac*cgrav
! Limit number of substeps if HG becomes prohibitively slow (use with care!!)
   if(maxtngrv>0) cgrav = min(cgrav,hgcfl*dtg_unit*dble(maxtngrv)/dt)
   dtgrav = min(dt,hgcfl*dtg_unit/cgrav)
   cgrav2 = cgrav**2
  end if

  call stop_clock(wttim)

 return
 end subroutine timestep

 pure function off(solve)
  logical,intent(in)::solve
  real(8):: off
  if(solve)then
   off = 1d0
  else
   off = 0d0
  end if
 end function off

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE DTI_CELL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute time step for a cell smeared over j-k coordinates

subroutine dti_cell(i,j_,k_,dti,jb,kb,cfmax)

 use settings,only:mag_on,eostype,solve_i,solve_j,solve_k
 use grid,only:js,je,ks,ke,js_global,je_global,ks_global,ke_global,&
               sa1,sa2,sa3,dvol
 use physval,only:d,eint,T,imu,p,cs,v1,v2,v3,b1,b2,b3,spc
 use pressure_mod,only:eos_p_cs,get_cf

 integer,intent(in):: i,j_,k_
 real(8),allocatable,intent(inout)::dti(:,:,:)
 integer,intent(in),optional:: jb,kb
 real(8),intent(out),optional:: cfmax
 integer:: j,k,jn,kn,ierr,jl,jr,kl,kr
 real(8):: cf1,cf2,cf3

!-----------------------------------------------------------------------------

 jn=0;kn=0
 if(present(jb)) jn = min(jb-1,je_global-js_global)
 if(present(kb)) kn = min(kb-1,ke_global-ks_global)

 jl = max(j_,js); jr = min(j_+jn,je)
 kl = max(k_,ks); kr = min(k_+kn,ke)
 j = jl ; k = kl

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

 j=j_;k=k_
 dti(i,jl:jr,kl:kr) = sum(dvol(i,j:j+jn,k:k+kn)) &
                    / ( cf1*off(solve_i) * sum(sa1(i-1:i   ,j:j+jn,k:k+kn)) &
                      + cf2*off(solve_j) * sum(sa2(i,j-1:j+jn:jn+1,k:k+kn)) &
                      + cf3*off(solve_k) * sum(sa3(i,j:j+jn,k-1:k+kn:kn+1)) )

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

 use settings,only:solve_i,solve_j,solve_k
 use grid,only:js,je,ks,ke,js_global,je_global,ks_global,ke_global,&
               sa1,sa2,sa3,dvol

 integer,intent(in):: i,j,k
 real(8),allocatable,intent(inout)::dtg(:,:,:)
 real(8),intent(in):: cgrav
 integer,intent(in),optional:: jb,kb
 integer:: jn,kn,jl,jr,kl,kr

!-----------------------------------------------------------------------------

 jn=0;kn=0
 if(present(jb)) jn = min(jb-1,je_global-js_global)
 if(present(kb)) kn = min(kb-1,ke_global-ks_global)

 jl = max(j,js); jr = min(j+jn,je)
 kl = max(k,ks); kr = min(k+kn,ke)

 dtg(i,jl:jr,kl:kr) = sum(dvol(i,j:j+jn,k:k+kn)) &
                      / ( cgrav * &
                         ( off(solve_i) * sum(sa1(i-1:i   ,j:j+jn,k:k+kn)) &
                         + off(solve_j) * sum(sa2(i,j-1:j+jn:jn+1,k:k+kn)) &
                         + off(solve_k) * sum(sa3(i,j:j+jn,k-1:k+kn:kn+1)) ) )

return
end subroutine dtgrav_cell

end module timestep_mod
