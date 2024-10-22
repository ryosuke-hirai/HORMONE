module smear_mod
 implicit none

 public:: smear_setup,smear
 private:: angular_sums_hydro,angular_sums_grav,&
           angular_smear_hydro,angular_smear_grav,angular_smear_e,&
           get_compensation,get_id

 integer,public:: nsmear
 integer,private:: nsmear_mydom, nsmear_split
 integer,allocatable,public:: block_j(:),block_k(:),lijk_from_id(:,:)
 integer,allocatable,private:: list_split(:),list_mydom(:,:)
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
  use mpi_domain,only:fully_my_domain,partially_my_domain
  use mpi_utils,only:allreduce_mpi

  integer:: n,nn,i,j,k,jb,kb
  integer,allocatable:: temp1(:),temp2(:)

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
   allocate(dvol_block(1:nsmear),lijk_from_id(0:3,1:nsmear),temp1(1:nsmear))
   allocate(temp2,mold=temp1)
   nsmear_mydom = 0 ; temp1 = 0 ; temp2 = 0
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
       if(partially_my_domain(i,j,k,1,jb,kb))then
        nsmear_mydom = nsmear_mydom+1
        temp1(nsmear_mydom) = get_id(i,j,k)
       end if
       if(fully_my_domain(i,j,k,1,jb,kb))then
        temp2(get_id(i,j,k)) = 1
       end if
      end do
     end do
    end do
   end do

! Record list of effective cells that are covered by my domain
   allocate(list_mydom(1:2,1:nsmear_mydom))
   list_mydom(1,1:nsmear_mydom) = temp1(1:nsmear_mydom)
   list_mydom(2,1:nsmear_mydom) = 0

! Record list of effective cells that are split over multiple MPI ranks
   call allreduce_mpi('sum',temp2)
   nsmear_split = nsmear - sum(temp2)
   allocate(list_split(1:nsmear_split))
   nsmear_split = 0
   list_split = 0
   do n = 1, nsmear
    if(temp2(n)==0)then
     nsmear_split = nsmear_split + 1
     do nn = 1, nsmear_mydom
      if(list_mydom(1,nn)==n)then
       list_mydom(2,nn) = nsmear_split
       list_split(nsmear_split) = nn
       exit
      end if
     end do
    end if
   end do

  end if

  return
 end subroutine smear_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                             SUBROUTINE SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To smear out quantities for nested-grid-like feature

 subroutine smear(which)

  use settings,only:spn,compswitch
  use grid,only:fmr_max,dim,crdnt
  use profiler_mod
  use mpi_utils,only:allreduce_mpi

  character(len=*),intent(in):: which
  integer:: i,j,k,l,n,nn,m,jb,kb,wtind,wtin2
  real(8),allocatable,dimension(:,:):: exchange,momtot,spctot
  real(8),allocatable,dimension(:):: spctot1,mtot,etot,phitot,psitot
  real(8):: mtot1,momtot1(1:3),etot1,phitot1,psitot1

!-----------------------------------------------------------------------------

  if(fmr_max==0.or.dim==1)return

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.

  select case(which)
  case('hydro')
   wtind = wtsmr ; wtin2 = wtsm2
   allocate(exchange(4+spn,max(nsmear_split,1)),spctot1(max(spn,1)),&
            mtot(nsmear_mydom),etot(nsmear_mydom),momtot(3,nsmear_mydom),&
            spctot(max(spn,1),nsmear_mydom))
  case('grav')
   wtind = wtgsm ; wtin2 = wtgs2
   allocate(exchange(2,max(nsmear_split,1)),&
            phitot(nsmear_mydom),psitot(nsmear_mydom))
  end select
  exchange = 0d0

  call start_clock(wtind)

  if(crdnt==2)then

! Integrate quantities over effective cells ++++++++++++++++++++++++++++++++++
   if(nsmear_mydom>0)then
!$omp parallel do private(i,j,k,l,n,nn,m,jb,kb,mtot1,momtot1,etot1,spctot1,&
!$omp phitot1,psitot1) schedule(dynamic,1)
    do nn = 1, nsmear_mydom
     n = list_mydom(1,nn)
     l = lijk_from_id(0,n)
     i = lijk_from_id(1,n)
     j = lijk_from_id(2,n)
     k = lijk_from_id(3,n)
     jb = block_j(l)
     kb = block_k(l)

     select case(which)
     case('hydro')
      call angular_sums_hydro(i,j,j+jb-1,k,k+kb-1,mtot1,momtot1,etot1,spctot1)
      mtot(nn) = mtot1
      momtot(1:3,nn) = momtot1(1:3)
      etot(nn) = etot1
      if(compswitch>=2) spctot(1:spn,nn) = spctot1(1:spn)

      if(list_mydom(2,nn)>0)then
       m = list_mydom(2,nn)
       exchange(1  ,m) = mtot(nn)
       exchange(2:4,m) = momtot(1:3,nn)
       if(compswitch>=2) exchange(5:4+spn,m) = spctot(1:spn,nn)
      end if

     case('grav')
      call angular_sums_grav(i,j,j+jb-1,k,k+kb-1,phitot1,psitot1)
      phitot(nn) = phitot1
      psitot(nn) = psitot1
      if(list_mydom(2,nn)>0)then
       m = list_mydom(2,nn)
       exchange(1,m) = phitot(nn)
       exchange(2,m) = psitot(nn)
      end if

     end select
    end do
!$omp end parallel do
   end if

! Reduce summed variables over all MPI ranks for split cells +++++++++++++++++
   if(nsmear_split>0)then
    call start_clock(wtin2)
    call allreduce_mpi('sum',exchange)
!$omp parallel do private(nn,n) schedule(dynamic,1)
    do nn = 1, nsmear_split
     n = list_split(nn)
     if(n>0)then
      select case(which)
      case('hydro')
       mtot(n) = exchange(1,nn)
       momtot(1:3,n) = exchange(2:4,nn)
       if(compswitch>=2) spctot(1:spn,n) = exchange(5:4+spn,nn)
      case('grav')
       phitot(n) = exchange(1,nn)
       psitot(n) = exchange(2,nn)
      end select
     end if
    end do
!$omp end parallel do
    if(which=='hydro')then
     deallocate(exchange)
     allocate(exchange(1,nsmear_split))
!$omp parallel workshare
     exchange = 0d0
!$omp end parallel workshare
    end if
    call stop_clock(wtin2)
   end if

! Average out quantities over effective cells ++++++++++++++++++++++++++++++++
   if(nsmear_mydom>0)then
!$omp parallel do private(n,nn,i,j,k,l,m,jb,kb) schedule(dynamic,1)
    do nn = 1, nsmear_mydom
     n = list_mydom(1,nn)
     l = lijk_from_id(0,n)
     i = lijk_from_id(1,n)
     j = lijk_from_id(2,n)
     k = lijk_from_id(3,n)
     jb = block_j(l)
     kb = block_k(l)
     select case(which)
     case('hydro')
      call angular_smear_hydro(i,j,j+jb-1,k,k+kb-1,&
                               mtot(nn),momtot(1:3,nn),etot(nn),&
                               spctot(1:max(spn,1),nn))
      if(list_mydom(2,nn)>0)then
       m = list_mydom(2,nn)
       exchange(1,m) = etot(nn)
      end if

     case('grav')
      call angular_smear_grav(i,j,j+jb-1,k,k+kb-1,phitot(nn),psitot(nn))
     end select
    end do
!$omp end parallel do
   end if

! Adjust energy over effective cells +++++++++++++++++++++++++++++++++++++++++
   if(which=='hydro')then
    if(nsmear_split>0)then
     call allreduce_mpi('sum',exchange)
     do nn = 1, nsmear_split
      n = list_split(nn)
      if(n>0)etot(n) = exchange(1,nn)
     end do
    end if
    if(nsmear_mydom>0)then
!$omp parallel do private(nn,n,i,j,k,l,jb,kb) schedule(dynamic,1)
     do nn = 1, nsmear_mydom
      n = list_mydom(1,nn)
      l = lijk_from_id(0,n)
      i = lijk_from_id(1,n)
      j = lijk_from_id(2,n)
      k = lijk_from_id(3,n)
      jb = block_j(l)
      kb = block_k(l)
      call angular_smear_e(i,j,j+jb-1,k,k+kb-1,etot(nn))
     end do
!$omp end parallel do
    end if
   end if

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
!                       SUBROUTINE ANGULAR_SMEAR_E
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out energy over several cells in the jk direction

 subroutine angular_smear_e(i,js_,je_,ks_,ke_,etot)

  use grid,only:js,je,ks,ke
  use physval,only:u,iene

  integer,intent(in):: i,js_,je_,ks_,ke_
  real(8),intent(in)::etot
  real(8):: vol
  integer:: jl,jr,kl,kr

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)

  vol = dvol_block(get_id(i,js_,ks_))

  u(i,jl:jr,kl:kr,iene) = etot / vol  
  
 end subroutine angular_smear_e
 
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                     SUBROUTINE ANGULAR_SMEAR_HYDRO
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out hydro quantities over several cells in the jk direction

 subroutine angular_smear_hydro(i,js_,je_,ks_,ke_,mtot,momtot,etot,spctot)

  use settings,only:spn,compswitch
  use grid,only:x3,dvol,car_x,js,je,ks,ke
  use physval,only:u,spc,v1,v2,v3,icnt,imo1,imo2,imo3
  use utils,only:get_vpol
  use gravmod,only:totphi,gravswitch

  integer,intent(in)::i,js_,je_,ks_,ke_
  real(8),intent(in):: mtot,momtot(1:3),spctot(1:max(spn,1))
  real(8),intent(inout):: etot
  integer:: n,j,k,jl,jr,kl,kr
  real(8):: vol, dave, vave(1:3)

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)

  vol = dvol_block(get_id(i,js_,ks_))
  dave = mtot/vol    ! get average density
  vave = momtot/mtot ! get average cartesian velocity

  if(compswitch>=2)then
   do n = 1, spn
    spc(n,i,jl:jr,kl:kr) = spctot(n) / mtot
   end do
  end if
  
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

 end subroutine angular_smear_hydro
 
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE ANGULAR_SUMS_HYDRO
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate total mass/momentum/energy/species over effective cells

 subroutine angular_sums_hydro(i,js_,je_,ks_,ke_,mtot,momtot,etot,spctot)

  use settings,only:spn,compswitch
  use grid,only:x3,dvol,car_x,js,je,ks,ke
  use physval,only:u,spc,icnt,iene,imo1,imo2,imo3
  use utils
  use gravmod,only:totphi,gravswitch

  integer,intent(in)::i,js_,je_,ks_,ke_
  real(8),intent(out):: mtot,momtot(1:3),etot
  real(8),allocatable,intent(inout):: spctot(:)
  integer:: n,j,k,jl,jr,kl,kr
  real(8),dimension(1:3):: compen, tempsum, element, vcar

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)

  mtot = sum(u(i,jl:jr,kl:kr,icnt)*dvol(i,jl:jr,kl:kr))

  if(compswitch>=2)then
   do n = 1, spn
    spctot(n) = sum(spc(n,i,jl:jr,kl:kr)*u(i,jl:jr,kl:kr,icnt)&
                                     *dvol(i,jl:jr,kl:kr))
   end do
  end if

  momtot=0d0;etot=0d0;compen=0d0
  do k = kl, kr
   do j = jl, jr
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
  momtot = momtot + compen

 end subroutine angular_sums_hydro

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE ANGULAR_SUMS_GRAV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Integrate the gravity quantities over angular effective cells

 subroutine angular_sums_grav(i,js_,je_,ks_,ke_,phitot,psitot)

  use grid,only:dvol,js,je,ks,ke
  use gravmod,only:grvphi,grvpsi

  integer,intent(in)::i,js_,je_,ks_,ke_
  real(8),intent(out):: phitot,psitot
  integer:: jl,jr,kl,kr

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)

  phitot = sum( grvphi(i,jl:jr,kl:kr)*dvol(i,jl:jr,kl:kr) )
  psitot = sum( grvpsi(i,jl:jr,kl:kr)*dvol(i,jl:jr,kl:kr) )

 end subroutine angular_sums_grav

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE ANGULAR_SMEAR_GRAV
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Smear the gravity quantities over angular effective cells

 subroutine angular_smear_grav(i,js_,je_,ks_,ke_,phitot,psitot)

  use grid,only:js,je,ks,ke
  use gravmod,only:grvphi,grvpsi

  integer,intent(in)::i,js_,je_,ks_,ke_
  real(8),intent(in):: phitot,psitot
  integer:: jl,jr,kl,kr
  real(8):: vol

!-----------------------------------------------------------------------------

  jl = max(js_,js); jr = min(je_,je)
  kl = max(ks_,ks); kr = min(ke_,ke)

  vol = dvol_block(get_id(i,js_,ks_))

  grvphi(i,jl:jr,kl:kr) = phitot / vol
  grvpsi(i,jl:jr,kl:kr) = psitot / vol

 end subroutine angular_smear_grav


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
