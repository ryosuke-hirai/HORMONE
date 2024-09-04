module smear_mod
 implicit none

 public:: smear_setup,smear
 private:: angular_smear,angular_smear_global

 integer,public:: nsmear
 integer,allocatable,public:: block_j(:),block_k(:),lijk_from_id(:,:)
 real(8),allocatable,private:: dvol_block(:)

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SMEAR_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up effective grid for nested-grid-like feature

 subroutine smear_setup

  use grid,only:is_global,js_global,je_global,ks_global,ke_global,&
                fmr_max,fmr_lvl,crdnt,dvol

  integer:: n,i,j,k,jb,kb

!-----------------------------------------------------------------------------

  fmr_lvl(0) = 0

  if(crdnt==2)then

! Determine block structure for each level
   allocate(block_j(1:fmr_max))
   block_j = 0
   allocate(block_k,source=block_j)

   do n = 1, fmr_max
    if(fmr_lvl(n)==0)cycle
    if(n==1)then
     block_j(n) = je_global-js_global+1
     block_k(n) = ke_global-ks_global+1
    else
     block_j(n) = min(2**(fmr_max-n+1),je_global-js_global+1)
     block_k(n) = min(2**(fmr_max-n+1),ke_global-ks_global+1)
    end if
   end do

! Count number of smeared effective cells
   nsmear = fmr_lvl(1)
   do n = 2, fmr_max
    nsmear = nsmear + fmr_lvl(n) &
                      * (je_global-js_global+1)/block_j(n) &
                      * (ke_global-ks_global+1)/block_k(n)
   end do

! Record effective cell properties
   allocate(dvol_block(1:nsmear),lijk_from_id(0:3,1:nsmear))
   do n = 1, fmr_max
    if(fmr_lvl(n)==0)cycle
    jb = block_j(n)
    kb = block_k(n)
    do k = ks_global, ke_global, kb
     do j = js_global, je_global, jb
      do i = is_global+sum(fmr_lvl(0:n-1)), is_global+sum(fmr_lvl(0:n))-1
       dvol_block(get_id(i,j,k)) = sum(dvol(i,j:j+jb-1,k:k+kb-1))
       lijk_from_id(0,get_id(i,j,k)) = n
       lijk_from_id(1,get_id(i,j,k)) = i
       lijk_from_id(2,get_id(i,j,k)) = j
       lijk_from_id(3,get_id(i,j,k)) = k
      end do
     end do
    end do
   end do

  end if

  return
 end subroutine smear_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To smear out quantities for nested-grid-like feature

 subroutine smear(which)

  use grid,only:is_global,js_global,je_global,ks_global,ke_global,&
                fmr_max,fmr_lvl,dim,crdnt
  use profiler_mod
  use mpi_domain,only:fully_my_domain
  use mpi_utils,only:allreduce_mpi

  character(len=*),intent(in):: which
  integer:: i,j,k,l,n,jb,kb,wtind,wtin2
  integer,allocatable:: smeared(:)

!-----------------------------------------------------------------------------

  if(fmr_max==0.or.dim==1)return

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.

  select case(which)
  case('hydro')
   wtind = wtsmr ; wtin2 = wtsm2
  case('grav')
   wtind = wtgsm ; wtin2 = wtgs2
  end select

  call start_clock(wtind)

  if(crdnt==2)then

   allocate( smeared(1:nsmear) )
   smeared = 0

! First sweep
!!$!$omp parallel
!!$   do n = 1, fmr_max
!!$    if(fmr_lvl(n)==0)cycle
!!$    jb = block_j(n)
!!$    kb = block_k(n)
!!$!$omp do private(i,j,k) collapse(3) schedule(dynamic)
!!$    do k = ks_global, ke_global, kb
!!$     do j = js_global, je_global, jb
!!$      do i = is_global+sum(fmr_lvl(0:n-1)), is_global+sum(fmr_lvl(0:n))-1
!!$       if(fully_my_domain(i,j,k,0,jb,kb))then
!!$        select case(which)
!!$        case('hydro')
!!$         call angular_smear(i,j,j+jb-1,k,k+kb-1)
!!$        case('grav')
!!$         call angular_smear_grav(i,j,j+jb-1,k,k+kb-1)
!!$        end select
!!$        smeared(get_id(i,j,k)) = 1
!!$       end if
!!$      end do
!!$     end do
!!$    end do
!!$!$omp end do
!!$   end do
!!$!$omp end parallel

!$omp parallel do private(i,j,k,l,jb,kb) schedule(dynamic)
   do n = 1, nsmear
    l = lijk_from_id(0,n)
    if(fmr_lvl(l)==0)cycle
    i = lijk_from_id(1,n)
    j = lijk_from_id(2,n)
    k = lijk_from_id(3,n)
    jb = block_j(l)
    kb = block_k(l)
    if(fully_my_domain(i,j,k,0,jb,kb))then
     select case(which)
     case('hydro')
      call angular_smear(i,j,j+jb-1,k,k+kb-1)
     case('grav')
      call angular_smear_grav(i,j,j+jb-1,k,k+kb-1)
     end select
     smeared(n) = 1
    end if
   end do
!$omp end parallel do

   call allreduce_mpi('sum',smeared)

! Second sweep (for effective cells that span over multiple MPI ranks)
!!$   if(minval(smeared)==0)then
!!$    call start_clock(wtin2)
!!$    do n = 1, fmr_max
!!$     if(fmr_lvl(n)==0)cycle
!!$     jb = block_j(n)
!!$     kb = block_k(n)
!!$     do k = ks_global, ke_global, kb
!!$      do j = js_global, je_global, jb
!!$       do i = is_global+sum(fmr_lvl(0:n-1)), is_global+sum(fmr_lvl(0:n))-1
!!$        if(smeared(get_id(i,j,k))==0)then
!!$         select case(which)
!!$         case('hydro')
!!$          call angular_smear_global(i,j,j+jb-1,k,k+kb-1)
!!$         case('grav')
!!$          call angular_smear_grav_global(i,j,j+jb-1,k,k+kb-1)
!!$         end select
!!$        end if
!!$       end do
!!$      end do
!!$     end do
!!$    end do
!!$    call stop_clock(wtin2)
!!$   end if

   if(minval(smeared)==0)then
    call start_clock(wtin2)
    do n = 1, nsmear
     if(smeared(n)==0)then
      l = lijk_from_id(0,n)
      if(fmr_lvl(l)==0)cycle
      i = lijk_from_id(1,n)
      j = lijk_from_id(2,n)
      k = lijk_from_id(3,n)
      jb = block_j(l)
      kb = block_k(l)
      select case(which)
      case('hydro')
       call angular_smear_global(i,j,j+jb-1,k,k+kb-1)
      case('grav')
       call angular_smear_grav_global(i,j,j+jb-1,k,k+kb-1)
      end select
     end if
    end do
    call stop_clock(wtin2)
   end if

   deallocate(smeared)

  end if

  call stop_clock(wtind)

  return
 end subroutine smear

! Get smeared cell id from i,j,k \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 function get_id(i,j,k) result(id)

  use grid,only:is_global,js_global,je_global,ks_global,ke_global,&
                fmr_max,fmr_lvl

  integer,intent(in):: i,j,k
  integer:: n, id, jn,kn

  do n = 1, fmr_max
   if(i<=is_global+sum(fmr_lvl(0:n))-1)exit

   if(n==1)then
    id = fmr_lvl(1)
   else
    jn = (je_global-js_global+1)/block_j(n)
    kn = (ke_global-ks_global+1)/block_k(n)
    id = id + fmr_lvl(n)*jn*kn
   end if
  end do

  if(n>fmr_max)error stop 1

  if(n==1)then
   id = i
  else
   jn = (je_global-js_global+1)/block_j(n)
   kn = (ke_global-ks_global+1)/block_k(n)
   id = id + (i-sum(fmr_lvl(0:n-1))-1) * jn * kn &
           + int(j/block_j(n)) + int(k/block_k(n))*jn + 1
  end if

  return
 end function get_id

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE ANGULAR_SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out hydro quantities over several cells in the jk direction

 subroutine angular_smear(i,js_,je_,ks_,ke_)

  use settings,only:spn,compswitch
  use grid,only:x3,dvol,car_x
  use physval,only:u,spc,v1,v2,v3,icnt,iene,imo1,imo2,imo3
  use utils
  use gravmod,only:totphi,gravswitch

  implicit none

  integer,intent(in)::i,js_,je_,ks_,ke_
  integer:: n,j,k
  real(8):: mtot, etot, spctot, vol, dave
  real(8),dimension(1:3):: momtot, compen, tempsum, element, vcar, vave

!-----------------------------------------------------------------------------

  vol = dvol_block(get_id(i,js_,ks_))
  mtot = sum(u(i,js_:je_,ks_:ke_,icnt)*dvol(i,js_:je_,ks_:ke_))

  if(compswitch>=2)then
   do n = 1, spn
    spctot = sum(spc(n,i,js_:je_,ks_:ke_)*u(i,js_:je_,ks_:ke_,icnt)&
                                         *dvol(i,js_:je_,ks_:ke_))
    spc(n,i,js_:je_,ks_:ke_) = spctot / mtot
   end do
  end if

  momtot=0d0;etot=0d0;compen=0d0
  do k = ks_, ke_
   do j = js_, je_
    call get_vcar(car_x(:,i,j,k),x3(k),&
                  u(i,j,k,imo1),u(i,j,k,imo2),u(i,j,k,imo3),vcar)
! Use Kahan summation algorithm to minimize roundoff error
    element = vcar*dvol(i,j,k)
    tempsum = momtot + element
    compen = get_compensation(compen,element,tempsum,momtot)
    momtot = tempsum ! add up momenta using the Kahan algorithm

    etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
    if(gravswitch>0)& ! and gravitational energy
     etot = etot + u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
   end do
  end do

  dave = mtot/vol    ! get average density
  vave = momtot/mtot ! get average cartesian velocity

  do k = ks_, ke_
   do j = js_, je_
    u(i,j,k,icnt) = dave ! density
    call get_vpol(car_x(:,i,j,k),x3(k),vave,&
                  v1(i,j,k),v2(i,j,k),v3(i,j,k))
    u(i,j,k,imo1) = v1(i,j,k)*dave
    u(i,j,k,imo2) = v2(i,j,k)*dave
    u(i,j,k,imo3) = v3(i,j,k)*dave
    if(gravswitch>0)&
     etot = etot - u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
   end do
  end do

  u(i,js_:je_,ks_:ke_,iene) = etot / vol

 end subroutine angular_smear

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE ANGULAR_SMEAR_GRAV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out gravity quantities over several cells in the jk direction

 subroutine angular_smear_grav(i,js_,je_,ks_,ke_)

  use grid,only:dvol
  use gravmod,only:grvphi,grvpsi

  implicit none

  integer,intent(in)::i,js_,je_,ks_,ke_
  real(8):: vol

!-----------------------------------------------------------------------------

  vol = dvol_block(get_id(i,js_,ks_))

  grvphi(i,js_:je_,ks_:ke_) = sum( grvphi(i,js_:je_,ks_:ke_) &
                                 * dvol  (i,js_:je_,ks_:ke_) ) / vol
  grvpsi(i,js_:je_,ks_:ke_) = sum( grvpsi(i,js_:je_,ks_:ke_) &
                                 * dvol  (i,js_:je_,ks_:ke_) ) / vol

 end subroutine angular_smear_grav

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                     SUBROUTINE ANGULAR_SMEAR_GLOBAL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out quantities over several cells in the angular direction
!          Particularly for smearing over multiple MPI ranks

 subroutine angular_smear_global(i,js_,je_,ks_,ke_)

  use settings,only:spn,compswitch
  use grid,only:x3,dvol,js,je,ks,ke,car_x
  use physval,only:u,spc,v1,v2,v3,icnt,iene,imo1,imo2,imo3
  use utils
  use gravmod,only:totphi,gravswitch
  use mpi_utils,only:allreduce_mpi
  use mpi_domain,only:sum_global_array,partially_my_domain

  implicit none

  integer,intent(in)::i,js_,je_,ks_,ke_
  integer:: n,j,k,jl,jr,kl,kr
  real(8):: mtot, etot, vol, dave, arr_sum
  real(8),allocatable:: spctot(:), comm_chunk(:)
  real(8),dimension(1:3):: momtot, compen, tempsum, element, vcar, vave
  logical:: overlap

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)
  overlap = partially_my_domain(i,js_,ks_,0,je_-js_+1,ke_-ks_+1)

  vol = dvol_block(get_id(i,js_,ks_))
  mtot = sum_global_array(u,i,i,js_,je_,ks_,ke_,icnt,weight=dvol)

! First smear chemical elements
  if(compswitch>=2)then
   allocate(spctot(1:spn), comm_chunk(1:spn+3))
!$omp parallel
   do n = 1, spn
    arr_sum = 0d0
    if(overlap)then
!$omp do private(i,j,k) collapse(2) reduction(+:arr_sum)
     do k = kl, kr
      do j = jl, jr
       arr_sum = arr_sum + u(i,j,k,icnt)*dvol(i,j,k)*spc(n,i,j,k)
      end do
     end do
!$omp end do
    end if
    spctot(n) = arr_sum
!!$    spctot = sum_global_array(u,i,i,js_,je_,ks_,ke_,icnt, &
!!$                              l_weight2=n, weight=dvol, weight2=spc )
!!$    if (overlap) spc(n,i,jl:jr,kl:kr) = spctot / mtot
   end do
!$omp end parallel
!   call allreduce_mpi('sum',spctot)
   comm_chunk(4:spn+3) = spctot(1:spn)
  else
   allocate(comm_chunk(1:3))
  end if

  momtot=0d0;etot=0d0;compen=0d0
  if (overlap) then
!$omp parallel do private(j,k,element,vcar) collapse(2) reduction(+:etot)
   do k = kl, kr
    do j = jl, jr
     call get_vcar(car_x(:,i,j,k),x3(k),&
                   u(i,j,k,imo1),u(i,j,k,imo2),u(i,j,k,imo3),vcar)
! Use Kahan-Babuska-Neumaier summation algorithm to minimize roundoff error
     element = vcar*dvol(i,j,k)
!$omp critical
     tempsum = momtot + element
     compen = get_compensation(compen,element,tempsum,momtot)
     momtot = tempsum ! add up momenta using the Kahan algorithm
!$omp end critical
     etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
     if(gravswitch>0)& ! and gravitational energy
      etot = etot + u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
    end do
   end do
!$omp end parallel do
   momtot = momtot + compen
  end if

  dave = mtot/vol
  comm_chunk(1:3) = momtot(1:3) 
!  call allreduce_mpi('sum',momtot)
  call allreduce_mpi('sum',comm_chunk)
  momtot = comm_chunk(1:3)
  vave = momtot/mtot ! get average cartesian velocity

  if(compswitch>=2.and.overlap)then
   spctot(1:spn) = comm_chunk(4:spn+3)
!$omp parallel do private(n)
   do n = 1, spn
    spc(n,i,jl:jr,kl:kr) = spctot(n) / mtot    
   end do
!$omp end parallel do
  end if

  if (overlap) then
!$omp parallel do private(j,k) collapse(2)
   do k = kl, kr
    do j = jl, jr
     u(i,j,k,icnt) = dave ! density
     call get_vpol(car_x(:,i,j,k),x3(k),vave,&
                   v1(i,j,k),v2(i,j,k),v3(i,j,k))
     u(i,j,k,imo1) = v1(i,j,k)*dave
     u(i,j,k,imo2) = v2(i,j,k)*dave
     u(i,j,k,imo3) = v3(i,j,k)*dave
     if(gravswitch>0)&
      etot = etot - u(i,j,k,icnt)*totphi(i,j,k)*dvol(i,j,k)
    end do
   end do
!$omp end parallel do
  end if

  call allreduce_mpi('sum',etot)
  if (overlap) u(i,jl:jr,kl:kr,iene) = etot / vol


!!$  Jtot = 0d0; Itot=0d0; momtot=0d0; vol = sum(dvol(i,js:je,ks:ke))
!!$  do k = ks, ke
!!$   do j = js, je
!!$    Jtot = Jtot + u(i,j,k,imo3)*x1(i)*sinc(j)*dvol(i,j,k)! add up angular momenta
!!$    Itot = Itot + (x1(i)*sinc(j))**2*dvol(i,j,k)! add up moment of inertia
!!$    momtot(1) = momtot(1) + u(i,j,k,imo1)*dvol(i,j,k)! add up radial momenta
!!$    momtot(2) = momtot(2) + u(i,j,k,imo2)*dvol(i,j,k)! add up polar momenta
!!$    etot = etot + u(i,j,k,iene)*dvol(i,j,k)! add up energy
!!$    eave = eave + eint(i,j,k)*dvol(i,j,k)
!!$    if(gravswitch>0)then
!!$     etot = etot + u(i,j,k,icnt)*grvphi(i,j,k)*dvol(i,j,k)! and gravitational ene
!!$    end if
!!$    if(include_extgrv)then
!!$     etot = etot + u(i,j,k,icnt)*extgrv(i,j,k)*dvol(i,j,k)! and external gravity
!!$    end if
!!$   end do
!!$  end do
!!$  mtot = sum( u(i,js:je,ks:ke,icnt)*dvol(i,js:je,ks:ke) )
!!$  u(i,js:je,ks:ke,icnt) = mtot/vol ! density
!!$  Itot = Itot * u(i,js,ks,icnt)
!!$  vave(1:2) = momtot(1:2)/mtot
!!$  vave(3) = Jtot/Itot ! get average angular velocity
!!$  eave = eave / vol
!!$
!!$  do k = ks, ke
!!$   do j = js, je
!!$    v3(i,j,k) = vave(3)*x1(i)*sinc(j)
!!$    u(i,j,k,imo3) = u(i,j,k,icnt) * v3(i,j,k)
!!$    if(gravswitch>0)then
!!$     etot = etot - u(i,j,k,icnt)*grvphi(i,j,k)*dvol(i,j,k)
!!$    end if
!!$    if(include_extgrv)then
!!$     etot = etot - u(i,j,k,icnt)*extgrv(i,j,k)*dvol(i,j,k)
!!$    end if
!!$   end do
!!$  end do
!!$
!!$  vave(1:2) = vave(1:2) / norm2(vave(1:2))
!!$  etot = etot - eave*vol - 0.5d0*Itot*vave(3)**2
!!$  if(etot<0d0)then
!!$   if(abs(etot/(eave*vol))>1d-2)then
!!$    print*,'Error in smear second',i,js,je,ks,ke,sum(u(i,js:je,ks:ke,iene)*dvol(i,js:je,ks:ke)),etot,eave*vol,0.5d0*Itot*vave(3)**2
!!$    stop
!!$   else
!!$    etot = 0d0
!!$   end if
!!$  end if
!!$  vave(1:2) = vave(1:2) * sqrt(2d0/u(i,js,ks,icnt)*etot/vol)
!!$  do k = ks, ke
!!$   do j = js, je
!!$    v1(i,j,k) = vave(1)
!!$    v2(i,j,k) = vave(2)
!!$    u(i,j,k,imo1) = v1(i,j,k)*u(i,j,k,icnt)
!!$    u(i,j,k,imo2) = v2(i,j,k)*u(i,j,k,icnt)
!!$   end do
!!$  end do
!!$
!!$  do k = ks, ke
!!$   do j = js, je
!!$    u(i,j,k,iene) = eave + 0.5d0*(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k))*u(i,j,k,icnt)
!!$   end do
!!$  end do

  return
 end subroutine angular_smear_global

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                   SUBROUTINE ANGULAR_SMEAR_GRAV_GLOBAL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out gravity quantities over several cells in the jk direction

 subroutine angular_smear_grav_global(i,js_,je_,ks_,ke_)

  use grid,only:js,je,ks,ke,dvol
  use gravmod,only:grvphi,grvpsi
  use mpi_domain,only:sum_global_array,partially_my_domain

  implicit none

  integer,intent(in)::i,js_,je_,ks_,ke_
  integer:: jl,jr,kl,kr
  real(8):: vol, phiave, psiave

!-----------------------------------------------------------------------------

  vol = dvol_block(get_id(i,js_,ks_))

  phiave = sum_global_array(grvphi,i,i,js_,je_,ks_,ke_,weight=dvol) / vol
  psiave = sum_global_array(grvpsi,i,i,js_,je_,ks_,ke_,weight=dvol) / vol

  if(partially_my_domain(i,js_,ks_,0,je_-js_+1,ke_-ks_+1))then
   jl = max(js_,js); jr = min(je_,je)
   kl = max(ks_,ks); kr = min(ke_,ke)
   grvphi(i,jl:jr,kl:kr) = phiave
   grvpsi(i,jl:jr,kl:kr) = psiave
  end if

 end subroutine angular_smear_grav_global

 pure elemental function get_compensation(c0,y,t,s) result(c)
! Get compensation term for the Kahan-Babuska-Neumaier method
  real(8),intent(in):: c0,y,t,s
  real(8):: c

  if(abs(s)>=abs(y))then
   c = c0+(s-t)+y
  else
   c = c0+(y-t)+s
  end if

 end function get_compensation

end module smear_mod
