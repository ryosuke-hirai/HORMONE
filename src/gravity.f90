module gravity_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE GRAVITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate gravitational forces.

subroutine gravity

 use settings,only:eq_sym,wtime,igrv,courant
 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use utils,only:masscoordinate
 use omp_lib

 integer:: l, flgcg, gin, tngrav, grungen, jb, kb
 real(8):: rr, rrold, pAp, alpha, beta, phih, cgrav2, dtgrav, mind, h
 real(8),dimension(1:lmax):: gsrc, pp, absrob
 real(8),allocatable,dimension(:,:,:):: lapphi,newphi
 real(8),allocatable,dimension(:):: intphi
 real(8),allocatable,dimension(:):: x,y,z,r,aw
 real(8):: faco, facn, fact

!-----------------------------------------------------------------------------

 if(gravswitch==0.or.gravswitch==1)return
 
 wtime(igrv) = wtime(igrv) - omp_get_wtime()
 
 gin = gie - gis + 1

 if(gravswitch==2.or.(gravswitch==3.and.tn==0))then
  allocate( x(1:lmax), y(1:lmax), z(1:lmax), r(1:lmax), aw(1:lmax) )
 
  if(grav_init_other.and.gravswitch==3)return
! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  call gravbound

! cylindrical (equatorial+axial symmetry) ####################################
  if(je==1.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
   mind = minval(d(is:ie,js:je,ks:ke))
   if(gravswitch==3)then
    mind = sum(d(is:ie,js:je,ks:ke)*dvol(is:ie,js:je,ks:ke))&
         / (pi*xi1e**2*(xi3e-xi3s)) * 1d-2
    hgsrc=mind
   end if
   do l = 1, lmax
    i = modlimax(l) +gis-1
    k = (l-modlimax(l))/gin + gks
    x(l) = grvphi(i,js,k)
    gsrc(l) = 4d0*pi*G*mind*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*0.5d0)
    if(i>=is)then;if(i<=ie)then;if(k>=ks)then;if(k<=ke)then
     gsrc(l) = 4d0*pi*G*d(i,js,k)*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*0.5d0)
    end if;end if;end if;end if
    if(k==gks) gsrc(l) = gsrc(l) - x1 (i)*dxi1(i)*idx3(k  )*phi3i(i,k-1)
    if(i==gie) gsrc(l) = gsrc(l) - xi1(i)*dxi3(k)*idx1(i+1)*phi1o(i+1,k)
    if(k==gke) gsrc(l) = gsrc(l) - x1 (i)*dxi1(i)*idx3(k+1)*phi3o(i,k+1)
   end do

! spherical (axial symmetry) #################################################
  elseif(ke==1.and.crdnt==2.and.dim==2)then

! calculating b for Ax=b
   mind = minval(d(is:ie,js:je,ks:ke))
   do l=1,lmax
    i = modlimax(l) +gis-1
    j = (l-modlimax(l))/gin + gjs
    k = ks
    x(l) = grvphi(i,j,k)
    gsrc(l) = 4d0*pi*G*mind*x1(i)*x1(i)*sinc(j)*dxi1(i)*dxi2(j)
    if(i>=is)then;if(i<=ie)then;if(j>=js)then;if(j<=je)then
     gsrc(l) = 4d0*pi*G*d(i,j,k)*x1(i)*x1(i)*sinc(j)*dxi1(i)*dxi2(j)
    end if;end if;end if;end if
    if(i==gie) gsrc(l)= gsrc(l) - xi1(i)*xi1(i)*sinc(j)*dxi2(j)*idx1(i+1)&
                                  *phiio(i+1,j)
!   if(i==gis) gsrc(l)= gsrc(l) - xi1(i-1)*xi1(i-1)*sinc(j)*dxi2(j)*idx1(i)&
!        *phiii(i-1,j)
   end do

  end if


! set initial r0=p0
  call Avec(x)

  r = gsrc - aw

  call cctr(r)

  pp = z
  rrold = dot_product(r,z)

! start iteration --------------------------------------------------------

  do n = 1, lmax

   call Avec(pp)

   pAp = dot_product(pp,aw) ! (pk,Apk)

   alpha = rrold / pAp

!$omp parallel do private(l)
   do l = 1, lmax
    x(l) = x(l) + alpha*pp(l) !x_k+1
    r(l) = r(l) - alpha*aw(l) !r_k+1
   end do
!$omp end parallel do

   call cctr(r)

   rr = dot_product(z,r) ! ((CC^T)^{-1}rk+1,rk+1)

   beta = rr / rrold

   rrold = rr !((CC^T)^{-1}rk,rk)

!++++++ criterion for convergence ++++++!
   flgcg = 0                            !
   do l = 1,lmax                        !
    absrob(l) = abs(r(l)/gsrc(l))       !
    if(absrob(l)>cgerr)then             !
     flgcg = 1                          !
     exit                               !
    end if                              !
   end do                               !
   if(flgcg==0)exit                     !
!+++++++++++++++++++++++++++++++++++++++!

!$omp parallel do private(l)
   do l = 1, lmax
    pp(l) = z(l) + beta * pp(l)
   end do
!$omp end parallel do

  end do

 if(n>=lmax)then
  print *,'Error from miccg.f',tn
 end if
!-------------------------------------------------------------------------
 if(je==1.and.crdnt==1.and.dim==2)then ! for cylindrical coordinates
! convert x to phi
  do k = gks, gke
   do i = gis, gie
    grvphi(i,js,k) = x(i-is+1+(k-gks)*gie)
   end do
  end do

  do k = gks,gke
   do j = js,je
    grvphi(is-2,j,k) = grvphi(is+1,j,k)
    grvphi(is-1,j,k) = grvphi(is  ,j,k)
   end do
  end do
  grvphi(gis:gie,js:je,gks-2) = grvphi(gis:gie,js:je,gks)
  grvphi(gis:gie,js:je,gks-1) = grvphi(gis:gie,js:je,gks)

  if(gravswitch==3)grvphiold = grvphi

 elseif(ke==1.and.crdnt==2.and.dim==2)then ! for spherical coordinates
  do j = js, je
   do i = is, gie
    grvphi(i,j,ks) = x(i-is+1+(j-js)*gie)
   end do
  end do

  grvphi(is-1:gie+1,js-2,ks) = grvphi(is-1:gie+1,js,ks)
  grvphi(is-1:gie+1,js-1,ks) = grvphi(is-1:gie+1,js,ks)
  grvphi(is-1:gie+1,je+1,ks) = grvphi(is-1:gie+1,je,ks)
  grvphi(is-1:gie+1,je+2,ks) = grvphi(is-1:gie+1,je,ks)

  grvphi(is-1,:,:) = grvphi(is,:,:)
  grvphi(is-2,:,:) = grvphi(is+1,:,:)
  grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                   ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)

  if(gravswitch==3)grvphiold = grvphi

 end if
endif


if(gravswitch==3.and.tn/=0)then
! Hyperbolic Self-Gravity $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(crdnt==1.and.je==1)then
! Cartoon mesh method for axially symmetric cylindrical coordinates %%%%%%%%%%
  allocate(newphi,mold=hgsrc)
  allocate( intphi(1:4) )

  cgrav2 = HGfac*max(maxval(cs(is:ie,js,ks:ke)+abs(v1(is:ie,js,ks:ke))), &
                     maxval(cs(is:ie,js,ks:ke)+abs(v3(is:ie,js,ks:ke))) )
  dtgrav = hgcfl*hg_dx/cgrav2
  cgrav2 = cgrav2**2

  hgsrc(is:ie,js:je,ks:ke) = d(is:ie,js:je,ks:ke)

  if(gbtype==0)call gravbound

  do while (grvtime<time+dt)

   if(grvtime+dtgrav>time+dt)dtgrav=time+dt-grvtime

   grvphi(gis-1,js,gks:gke) = grvphi(gis,js,gks:gke)
   if(eq_sym)grvphi(gis:gie,js,gks-1) = grvphi(gis:gie,js,gks)

   if(gbtype==1)then
!$omp parallel workshare
    grvphi(gie+1,js,gks:gke) = orgdis(gie,js,gks:gke)/orgdis(gie+1,js,gks:gke) &
                             * grvphi(gie,js,gks:gke)
    grvphi(gis:gie,js,gke+1) = orgdis(gis:gie,js,gke)/orgdis(gis:gie,js,gke+1) &
                             * grvphi(gis:gie,js,gke)
    grvphi(gis:gie,js,gks-1) = orgdis(gis:gie,js,gks)/orgdis(gis:gie,js,gks-1) &
                             * grvphi(gis:gie,js,gks)
!$omp end parallel workshare
   end if

!$omp parallel
!$omp do private(i,j,k,h,phih,intphi)
   do k = gks, gke
    do j = js, je
     do i = gis, gie
      h = min(dx1(i),dx1(i+1),dx3(k),dx3(k+1))
      intphi(1) = sum( lag11(-1:1,i,j,k)*grvphi(i-1:i+1,j,k) ) ! phi(i-1,j,k)
      intphi(2) = sum( lag12(-1:1,i,j,k)*grvphi(i-1:i+1,j,k) ) ! phi(i+1,j,k)
      intphi(3) = sum( lag31(-1:1,i,j,k)*grvphi(i,j,k-1:k+1) ) ! phi(i,j,k-1)
      intphi(4) = sum( lag32(-1:1,i,j,k)*grvphi(i,j,k-1:k+1) ) ! phi(i,j,k+1)
      phih = sum( lag21(-1:1,i,j,k)*grvphi(i-1:i+1,j,k) ) ! phi(i,j-1/j+1,k)
      newphi(i,j,k) = 0.5d0*cgrav2*dtgrav*(dtgrav+dt_old)* &
                    ( (sum(intphi(1:4))+2d0*phih-6d0*grvphi(i,j,k))/h**2 &
                      - 4d0*pi*G*hgsrc(i,j,k) ) & ! source term
                    - dtgrav/dt_old*grvphiold(i,j,k) &
                    + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
     end do
    end do
   end do
!$omp end do
!$omp workshare
   grvphiold = grvphi
   grvphi(gis:gie,js:je,gks:gke) = newphi(gis:gie,js:je,gks:gke)
   grvphi(gis-1,js:je,gks:gke) = grvphi(gis,js:je,gks:gke)
!$omp end workshare
!$omp end parallel

   if(eq_sym)then ! for plane symmetry
    grvphi(gis:gie,js:je,ks-1) = newphi(gis:gie,js:je,ks)
    grvphi(gis:gie,js:je,ks-2) = newphi(gis:gie,js:je,ks+1)
   end if

   dt_old  = dtgrav
   grvtime = grvtime + dtgrav

   if(maxval(grvphi(gis:gie,js:je,gks:gke))>=0d0)then
    print*,'Error in gravity: Positive gravitational potential'
    print*,'e.g. (i,j,k)=',maxloc(grvphi(gis:gie,js:je,gks:gke))
    stop
   end if

  end do

  deallocate(newphi,intphi)

 elseif(crdnt==2.and.ke==ks)then
! Axisymmetric spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  allocate(newphi,mold=hgsrc)
  
  k = gks

  cgrav2 = HGfac*max(maxval(cs(is:ie,js:je,k)+abs(v1(is:ie,js:je,k)) &
                                             +abs(v3(is:ie,js:je,k))) , &
                     maxval(cs(is:ie,js:je,k)+abs(v2(is:ie,js:je,k))) )
  dtgrav = hgcfl*hg_dx/cgrav2
  tngrav = ceiling((time+dt-grvtime)/dtgrav)
  dtgrav = (time+dt-grvtime)/dble(tngrav)
  cgrav2 = cgrav2**2

  if(gbtype==0)call gravbound

!$omp parallel private(i,j,l)
!$omp workshare
  hgsrc(is:ie,js:je,k) = 4d0*pi*G*d(is:ie,js:je,k)
!$omp end workshare

  do l = 1, tngrav
!$omp do collapse(2) schedule(static)
   do j = gjs, gje
    do i = gis, gie
     newphi(i,j,k) = 0.5d0*cgrav2*dtgrav*(dtgrav+dt_old)* &
                   ( hg11(i)*grvphi(i+1,j,k) + hg12(i)*grvphi(i-1,j,k) + &
                    (hg21(j)*grvphi(i,j+1,k) + hg22(j)*grvphi(i,j-1,k)) &
                     / x1(i)**2 + &
                     hg123(i,j,k)*grvphi(i,j,k) &
                   - hgsrc(i,j,k) ) & ! source term
                   - dtgrav/dt_old*grvphiold(i,j,k) &
                   + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
    end do
   end do
!$omp end do

! This has to be written this way for OpenMP optimization on NUMA processors
!$omp do collapse(2) schedule(static)
   do j = gjs, gje
    do i = gis, gie
     grvphiold(i,j,k) = grvphi(i,j,k)
     grvphi(i,j,k) = newphi(i,j,k)
    end do
   end do
!$omp end do 

!$omp single
   dt_old  = dtgrav
   grvtime = grvtime + dtgrav
!$omp end single nowait

  end do

! Boundary conditions
!$omp do
  do j = gjs, gje
   grvphi(gis-1,j,gks) = grvphi(gis,j,gks)   
  end do
!$omp end do nowait
!$omp do
  do i = gis, gie
   grvphi(i,gjs-1,gks) = grvphi(i,gjs,gks)
   grvphi(i,gje+1,gks) = grvphi(i,gje,gks)
  end do
!$omp end do nowait

  if(gbtype==1)then
!$omp workshare
   grvphi(gie+1,gjs:gje,k) = x1(gie)/x1(gie+1) * grvphi(gie,gjs:gje,k)
!$omp end workshare nowait
  end if
!$omp end parallel

  deallocate(newphi)

 elseif(crdnt==2.and.dim==3)then
! 3D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  allocate(newphi,mold=hgsrc)

!!$  cgrav2 = HGfac*max(maxval(cs(is:ie,js:je,ks:ke)+abs(v1(is:ie,js:je,ks:ke))), &
!!$                     maxval(cs(is:ie,js:je,ks:ke)+abs(v2(is:ie,js:je,ks:ke))), &
!!$                     maxval(cs(is:ie,js:je,ks:ke)+abs(v3(is:ie,js:je,ks:ke))) )
  cgrav2 = HGfac*ch
!  dtgrav = hgcfl*hg_dx/cgrav2
  dtgrav = hgcfl/courant/HGfac*dt
  tngrav = ceiling((time+dt-grvtime)/dtgrav)
  dtgrav = (time+dt-grvtime)/dble(tngrav)
  cgrav2 = cgrav2**2
 
  if(gbtype==0)call gravbound
  call masscoordinate
!$omp parallel
  do l = 1, tngrav
!$omp workshare
   grvphi(is-1,js:je,ks:ke) = grvphi(is,js:je,ks:ke)
   grvphi(is:ie,js-1,ks:ke) = grvphi(is:ie,js,ks:ke)
   grvphi(is:ie,je+1,ks:ke) = grvphi(is:ie,je,ks:ke)   
   grvphi(is:ie,js:je,ks-1) = grvphi(is:ie,js:je,ke)
   grvphi(is:ie,js:je,ke+1) = grvphi(is:ie,js:je,ks)
!$omp end workshare
   do grungen = 1, 3
! First set flux and source term
!$omp do private (i,j,k) collapse(3)
   do k = ks-1, ke
    do j = js-1, je
     do i = is-1, ie
      grv1(i,j,k) = -cgrav2*(grvphi(i+1,j,k)-grvphi(i,j,k))*idx1(i+1)
      grv2(i,j,k) = -cgrav2*(grvphi(i,j+1,k)-grvphi(i,j,k))*idx2(j+1)/x1(i)
      grv3(i,j,k) = -cgrav2*(grvphi(i,j,k+1)-grvphi(i,j,k))*idx3(k+1)/x1(i)/sin(x2(j))

      if(i>is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1)then
      elseif(i<=is+sphrn-1)then
       grv2(i,j,k) = 0d0;grv3(i,j,k) = 0d0
      elseif(i<=is+sphrn+trnsn16-1)then
       grv2(i,j,k) = grv2(i,j,k)/16d0 ; grv3(i,j,k) = grv3(i,j,k)/16d0
      elseif(i<=is+sphrn+trnsn16+trnsn8-1)then
       grv2(i,j,k) = grv2(i,j,k)/8d0 ; grv3(i,j,k) = grv3(i,j,k)/8d0
      elseif(i<=is+sphrn+trnsn16+trnsn8+trnsn4-1)then
       grv2(i,j,k) = grv2(i,j,k)/4d0 ; grv3(i,j,k) = grv3(i,j,k)/4d0
      elseif(i<=is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1)then
       grv2(i,j,k) = grv2(i,j,k)/2d0 ; grv3(i,j,k) = grv3(i,j,k)/2d0
      end if
      if(i==is-1.or.j==js-1.or.k==ks-1)cycle
      hgsrc(i,j,k) = -cgrav2*4d0*pi*G*d(i,j,k)
     end do
    end do
   end do
!$omp end do

   grk3_number: select case (grungen)
   case(1) grk3_number
    faco = 1d0 ; fact = 1d0 ; facn = 0d0
   case(2) grk3_number
    faco = 0.75d0 ; fact = 0.25d0 ; facn = fact
   case(3) grk3_number
    faco = 1d0/3d0 ; fact = 2d0/3d0 ; facn = fact
   end select grk3_number

!$omp do private (i,j,k) collapse(3)
    do k = ks,ke
     do j = js,je
      do i = is,ie
       if(grungen==1)then
        grvphiorg(i,j,k,1) = grvphi(i,j,k)
        grvphiorg(i,j,k,2) = grvphidot(i,j,k)
       end if
       grvphi(i,j,k) = faco*grvphiorg(i,j,k,1) + facn*grvphi(i,j,k) &
                     + fact*dtgrav*grvphidot(i,j,k)
       grvphidot(i,j,k) = faco*grvphiorg(i,j,k,2) + facn*grvphidot(i,j,k) &
             + fact*dtgrav * &
             ( idetg1(i) * &
               (detg1(i-1  )*grv1(i-1,j,k)-detg1(i  )*grv1(i,j,k)) &
             + idetg2(i,j) * &
               (detg2(i,j-1)*grv2(i,j-1,k)-detg2(i,j)*grv2(i,j,k)) &
             + idetg3(i,j,k) * (grv3(i,j,k-1)-grv3(i,j,k)) &
             + hgsrc(i,j,k) )
      end do
     end do
    end do
!$omp end do

!$omp single
    do i = sphrn, is, -1
     grvphi(i,js:je,ks:ke) = sum(grvphi(i+1,js:je,ks:ke)*dvol(i+1,js:je,ks:ke))&
      / sum(dvol(i+1,js:je,ks:ke)) &
      - G*mc(i)/xi1(i)**2*dx1(i+1)
    end do
!$omp end single
    
!!$!$omp do private(i)
!!$   do i = is, is+sphrn-1
!!$    grvphi(i,js:je,ks:ke) = sum(grvphi(i,js:je,ks:ke)*dvol(i,js:je,ks:ke)) &
!!$                          / sum(dvol(i,js:je,ks:ke))
!!$    grvphidot(i,js:je,ks:ke) = sum(grvphidot(i,js:je,ks:ke)*dvol(i,js:je,ks:ke))&
!!$                             / sum(dvol(i,js:je,ks:ke))
!!$   end do
!!$!$omp end do
!!$   jb=15;kb=15
!!$   if(ke==1)kb=0
!!$!$omp do private(i,j,k) collapse(3)
!!$   do i = is+sphrn, is+sphrn+trnsn16-1
!!$    do k = ks, ke, 16
!!$     do j = js, je, 16
!!$      grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$      grvphidot(i,j:j+jb,k:k+kb) = sum(grvphidot(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do
!!$   jb=7;kb=7
!!$   if(ke==1)kb=0
!!$!$omp do private(i,j,k) collapse(3)
!!$   do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
!!$    do k = ks, ke, 8
!!$     do j = js, je, 8
!!$      grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$      grvphidot(i,j:j+jb,k:k+kb) = sum(grvphidot(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do
!!$   jb=3;kb=3
!!$   if(ke==1)kb=0
!!$!$omp do private(i,j,k) collapse(3)
!!$   do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
!!$    do k = ks, ke, 4
!!$     do j = js, je, 4
!!$      grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$      grvphidot(i,j:j+jb,k:k+kb) = sum(grvphidot(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do
!!$   jb=1;kb=1
!!$   if(ke==1)kb=0
!!$!$omp do private(i,j,k) collapse(3)
!!$   do i = is+sphrn+trnsn16+trnsn8+trnsn4, is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
!!$    do k = ks, ke, 2
!!$     do j = js, je, 2
!!$      grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$      grvphidot(i,j:j+jb,k:k+kb) = sum(grvphidot(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/sum(dvol(i,j:j+jb,k:k+kb))
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do
  end do
!$omp single
  grvtime = grvtime + dtgrav
!$omp end single
 end do
!$omp end parallel


!!$!$omp parallel
!!$!$omp workshare
!!$  hgsrc(is:ie,js:je,ks:ke) = 4d0*pi*G*d(is:ie,js:je,ks:ke)
!!$!$omp end workshare
!!$
!!$  do l = 1, tngrav
!!$
!!$!$omp do private(i,j) collapse(2)
!!$   do j = gjs, gje
!!$    do i = gis, gie
!!$     grvphi(i,j,gks-1) = grvphi(i,j,gke)
!!$     grvphi(i,j,gke+1) = grvphi(i,j,gks)
!!$    end do
!!$   end do
!!$!$omp end do
!!$   
!!$!$omp do private(i,j,k) collapse(3) schedule(static)
!!$   do k = gks, gke
!!$    do j = gjs, gje
!!$     do i = gis, gie
!!$      newphi(i,j,k) = 0.5d0*cgrav2*dtgrav*(dtgrav+dt_old)* &
!!$                    ( hg11(i)*grvphi(i+1,j,k) + hg12(i)*grvphi(i-1,j,k)  &
!!$                    +(hg21(j)*grvphi(i,j+1,k) + hg22(j)*grvphi(i,j-1,k)) &
!!$                      / x1(i)**2 &
!!$                    +(hg31(k)*grvphi(i,j,k+1) + hg32(k)*grvphi(i,j,k-1)) &
!!$                      / (x1(i)*sinc(j))**2 &
!!$                    + hg123(i,j,k)*grvphi(i,j,k) &
!!$                    - hgsrc(i,j,k) ) & ! source term
!!$                    - dtgrav/dt_old*grvphiold(i,j,k) &
!!$                    + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do
!!$
!!$! This has to be written this way for OpenMP optimization on NUMA processors
!!$!$omp do private(i,j,k) schedule(static) collapse(3)
!!$   do k = gks, gke
!!$    do j = gjs, gje
!!$     do i = gis, gie
!!$      grvphiold(i,j,k) = grvphi(i,j,k)
!!$      grvphi(i,j,k)    = newphi(i,j,k)
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do nowait
!!$
!!$!$omp single
!!$   dt_old  = dtgrav
!!$   grvtime = grvtime + dtgrav
!!$!$omp end single nowait
!!$  end do
!!$
!!$!$omp workshare
!!$  grvphi(gis-1,gjs:gje,gks:gke) = grvphi(gis,gjs:gje,gks:gke)
!!$  grvphi(gis:gie,gjs-1,gks:gke) = grvphi(gis:gie,gjs,gks:gke)
!!$  grvphi(gis:gie,gje+1,gks:gke) = grvphi(gis:gie,gje,gks:gke)
!!$  grvphi(gis:gie,gjs:gje,gks-1) = grvphi(gis:gie,gjs:gje,gke)
!!$  grvphi(gis:gie,gjs:gje,gke+1) = grvphi(gis:gie,gjs:gje,gks)
!!$!$omp end workshare
!!$!$omp end parallel
!!$  
!!$  deallocate(newphi)
  
 end if


end if

wtime(igrv) = wtime(igrv) + omp_get_wtime()

contains

! \\\\\\\\\\\\\\\\\\\\\\\\\\\\ SUBROUTINE Avec \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
   subroutine Avec(w)

     use gravmod

     implicit none

  ! Subroutine to calculate A*vector

     real(8),intent(in),dimension(1:lmax)::w

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    aw(1) = a1(1)*w(1) + a2(1)*w(2) + a3(1)*w(1+gin)

   do l = 2,gin
    aw(l) = a2(l-1)*w(l-1) + a1(l)*w(l) + a2(l)*w(l+1) + a3(l)*w(l+gin)
   end do

   do l = gin+1,lmax-gin
    aw(l) = a3(l-gin)*w(l-gin) + a2(l-1)*w(l-1) + a1(l)*w(l) &
                                 + a2(l)*w(l+1) + a3(l)*w(l+gin)
   end do

   do l = lmax-gin+1,lmax-1
    aw(l) = a3(l-gin)*w(l-gin) + a2(l-1)*w(l-1) + a1(l)*w(l) + a2(l)*w(l+1)
   end do

    aw(lmax) = a3(lmax-gin)*w(lmax-gin) + a2(lmax-1)*w(lmax-1) &
               + a1(lmax)*w(lmax)

  end subroutine Avec



!\\\\\\\\\\\\\\\\\\\\\\\ SUBROUTINE (CC^T)^{-1}r \\\\\\\\\\\\\\\\\\\\\\\\\

  subroutine cctr(rrr)

    use gravmod

! subroutine to calculate (CC^T)^{-1}r for preconditioned CG method

    real(8),intent(in),dimension(1:lmax):: rrr

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! solving U^{T}DUr~=r by LU decomposition
! first solve U^{T}Dy=r by forward substitution
    y(1) = rrr(1)

    do l = 2, gin-1
     y(l) = rrr(l) - precb(l-1)/preca(l-1)*y(l-1)
    end do

    y(gin) = rrr(gin) - prece(1)/preca(1)*y(1) &
                      - precb(gin-1)/preca(gin-1)*y(gin-1)

    do l = gin+1, lmax
     y(l) = rrr(l) - precc(l-gin)/preca(l-gin)*y(l-gin) &
                   - prece(l-gin+1)/preca(l-gin+1)*y(l-gin+1) &
                   - precb(l-1)/preca(l-1)*y(l-1) 
    end do

! next solve Uz=y by backward substitution
    z(lmax) = y(lmax)/preca(lmax)

    do l = lmax-1, lmax-gin+2, -1
     z(l) = (y(l)-precb(l)*z(l+1))/preca(l)
    end do

    z(lmax-gin+1) = (y(lmax-gin+1) - precb(lmax-gin+1)*z(lmax-gin+2) &
                                     - prece(lmax-gin+1)*z(lmax))  &
                      / preca(lmax-gin+1)

    do l = lmax-gin, 1, -1
     z(l) = (y(l) - precb(l)*z(l+1) - prece(l)*z(l+gin-1) &
                  - precc(l)*z(l+gin)) / preca(l)
    end do

 end subroutine cctr
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
end subroutine gravity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRAVSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the matrix for Poisson equation

subroutine gravsetup
  
 use settings,only:eq_sym,courant
 use constants,only:huge
 use grid
 use gravmod
  
 integer:: l,gin
 real(8):: h

!-----------------------------------------------------------------------------
 gin = gie-gis+1

 if(gravswitch==2.or.gravswitch==3)then
! Extracting i and j from l
  do l = 1, lmax
   modlimax(l) = mod(l,gin)
   if(modlimax(l)==0) modlimax(l) = gie
  end do

! Constructing matrix A
!  2D cylindrical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(je==1.and.crdnt==1.and.dim==2)then
  
   do l = 1,lmax-gin ! calculating diagonal elements
    i = modlimax(l) +gis-1
    k = (l-modlimax(l))/gin + gks
    a1(l) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
              + 2d0*x1(i)*dxi1(i) / (dx3(k)*dx3(k+1)) ) &
            * ( (dx3(k)+dx3(k+1)) / 2d0 )
    a2(l) = xi1(i)*(dx3(k)+dx3(k+1))/2d0/dx1(i+1)
    a3(l) = x1(i)*dxi1(i)/dx3(k+1)
    if(mod(l,gin)==0)a2(l)=0d0
   end do

   do l=lmax-gin+1,lmax-1 ! k+1 line disappears
    i = modlimax(l) +gis-1
    k = (l-modlimax(l))/gin + gks
    a1(l) = -( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
             + 2d0*x1(i)*dxi1(i)/(dx3(k)*dx3(k+1)) ) &
          * ( (dx3(k)+dx3(k+1))/2.d0 )
    a2(l) = xi1(i)*(dx3(k)+dx3(k+1))/2d0/dx1(i+1)
    a3(l) = 0d0
    if(mod(l,gin)==0)a2(l)=0.0d0
   end do

   i = modlimax(lmax) +gis-1 ! i+1 line disappears
   k = (lmax-modlimax(lmax))/gin + gks
   a1(lmax) = -(xi1(i)/dx1(i+1)+xi1(i-1)/dx1(i) &
               +2d0*x1(i)*dxi1(i)/(dx3(k)*dx3(k+1))) &
               *((dx3(k)+dx3(k+1))/2d0)
   a2(lmax) = 0d0
   a3(lmax) = 0d0

   if(eq_sym)then ! for Neumann boundary at bc3i (equatorial symmetry)
    do l = 1, lmax
     i = modlimax(l) +gis-1
     k = (l-modlimax(l))/gin + gks
     if(k==ks)then
      a1(l) = a1(l) + x1(i)*dxi1(i)/dx3(k+1)
     end if
    end do
   end if

   call mic

! 2D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(ke==1.and.crdnt==2.and.dim==2)then
  
   do l = 1,lmax-gin ! calculating diagonal elements
    i = modlimax(l) +gis-1
    j = (l-modlimax(l))/gin + gjs
    a1(l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
            * sinc(j)*dxi2(j) &
            +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    a2(l) = xi1(i)**2 *sinc(j)*dxi2(j)/dx1(i+1)
    a3(l) = sini(j)*dxi1(i)/dx2(j+1)
    if(i==gis.and.xi1s>0d0)a1(l)=a1(l)+xi1(i-1)**2*sinc(j)*dxi2(j)/dx1(i)
    if(i==gie)a2(l)=0d0
   end do

   do l=lmax-gin+1,lmax-1 ! j+1 line disappears
    i = modlimax(l) +gis-1
    j = (l-modlimax(l))/gin + gjs
    a1(l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
            * sinc(j)*dxi2(j) &
            +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    a2(l) = xi1(i)**2 *sinc(j)*dxi2(j)/dx1(i+1)
    a3(l) = 0d0
    if(i==gis.and.xi1s>0d0)a1(l)=a1(l)+xi1(i-1)**2*sinc(j)*dxi2(j)/dx1(i)
    if(i==gie)a2(l)=0d0
   end do

   i = modlimax(lmax) +gis-1 ! i+1 line disappears
   j = (lmax-modlimax(lmax))/gin + gjs
   a1(lmax) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
               * sinc(j)*dxi2(j) &
               +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
   a2(lmax) = 0d0
   a3(lmax) = 0d0

   if(eq_sym)then ! for Neumann boundary at bc2o (equatorial symmetry)
    do l = 1, lmax
     i = modlimax(l) +gis-1
     j = (l-modlimax(l))/gin + gjs
     if(j==je)then
      a1(l) = a1(l) + sini(j)*dxi1(i)/dx2(j+1)
     end if
    end do
   end if


   call mic

  end if

! set initial x0

  if(tn==0.and.maxval(grvphi(gis:gie,gjs:gje,gks:gke))==0d0) grvphi = -1d3

  if(tn==0) dt_old = dt

 end if

! For Hyperbolic gravity solver ----------------------------------------------
 if(gravswitch==3)then
! for axisymmetric cylindrical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(je==1.and.dim==2.and.crdnt==1)then
! Cartoon mesh method
   allocate( hg123(gis-1:gie+1,1:1,gks-1:gke+1), &
             lag(-1:1,gis-1:gie+1),&
             orgdis(gis-1:gie+1,1:1,gks-1:gke+1) )
   allocate( hg11,hg12,hg21,hg22, mold=x1 )
   allocate( hg31,hg32, mold=x3 )

   do i = gis-1, gie+1
    hg11(i) = idx1(i+1)/sum(dx1(i:i+1))
    hg12(i) = idx1(i  )/sum(dx1(i:i+1))
   end do
   do i = gis, gie+1
    hg21(i) = idx1(i+1)*idx1(i+1) ! 1/h**2
    hg22(i) = sqrt(x1(i)*x1(i)+dx1(i+1)*dx1(i+1)) ! x(h)
! for Lagrange interpolation (second order)
    lag(-1,i) = (hg22(i)-x1(i))*(hg22(i)-x1(i+1))*idx1(i)/sum(dx1(i:i+1))
    lag( 0,i) = -(hg22(i)-x1(i-1))*(hg22(i)-x1(i+1))*idx1(i)*idx1(i+1)
    lag( 1,i) = (hg22(i)-x1(i-1))*(hg22(i)-x1(i))/sum(dx1(i:i+1))*idx1(i+1)
   end do
   do k = gks-1, gke+1
    hg31(k) = idx3(k+1)/sum(dx3(k:k+1))
    hg32(k) = idx3(k  )/sum(dx3(k:k+1))
   end do
   do k = gks-1, gke+1
    do i = gis-1, gie+1
     hg123(i,js,k) = idx1(i+1)*( idx1(i)+idx1(i+1) ) + idx3(k)*idx3(k+1)
     orgdis(i,js,k)= sqrt( x1(i)*x1(i) + x3(k)*x3(k) )
    end do
   end do

!Experimental for mapping non-uniform mesh to uniform mesh Laplacian
   allocate( lag11(-1:1,gis:gie,js:je,gks:gke) )
   allocate( lag12,lag21,lag31,lag32,mold=lag11 )
   do k = gks, gke
    do j = js, je
     do i = gis, gie
      h = min( dx1(i),dx1(i+1),dx3(k),dx3(k+1) )

      lag11(-1,i,j,k) = h*(dx1(i+1)+h)*idx1(i  )/sum(dx1(i:i+1))
      lag11( 0,i,j,k) = (dx1(i)-h)*(dx1(i+1)+h)*idx1(i)*idx1(i+1)
      lag11( 1,i,j,k) =-h*(dx1(i  )-h)*idx1(i+1)/sum(dx1(i:i+1))
      lag12(-1,i,j,k) =-h*(dx1(i+1)-h)*idx1(i  )/sum(dx1(i:i+1))
      lag12( 0,i,j,k) = (dx1(i)+h)*(dx1(i+1)-h)*idx1(i)*idx1(i+1)
      lag12( 1,i,j,k) = h*(dx1(i  )+h)*idx1(i+1)/sum(dx1(i:i+1))

      lag31(-1,i,j,k) = h*(dx3(k+1)+h)*idx3(k  )/sum(dx3(k:k+1))
      lag31( 0,i,j,k) = (dx3(k)-h)*(dx3(k+1)+h)*idx3(k)*idx3(k+1)
      lag31( 1,i,j,k) =-h*(dx3(k  )-h)*idx3(k+1)/sum(dx3(k:k+1))
      lag32(-1,i,j,k) =-h*(dx3(k+1)-h)*idx3(k  )/sum(dx3(k:k+1))
      lag32( 0,i,j,k) = (dx3(k)+h)*(dx3(k+1)-h)*idx3(k)*idx3(k+1)
      lag32( 1,i,j,k) = h*(dx3(k  )+h)*idx3(k+1)/sum(dx3(k:k+1))

      h = sqrt( x1(i)*x1(i)+h*h )! x(h)

      lag21(-1,i,j,k) = (h-x1(i))*(h-x1(i+1))*idx1(i  )/sum(dx1(i:i+1))
      lag21( 0,i,j,k) = (h-x1(i-1))*(x1(i+1)-h)*idx1(i)*idx1(i+1)
      lag21( 1,i,j,k) = (h-x1(i-1))*(h-x1(i))*idx1(i+1)/sum(dx1(i:i+1))

     end do
    end do
   end do

   hg_dx = huge
   do k = ks, ke
    do i = is, ie
     hg_dx = min(hg_dx,dxi1(i)*dxi3(k)/(dxi1(i)+dxi3(k)))
    end do
   end do

! for axisymmetrical spherical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.ke==1.and.dim==2)then
! Normal discretization
   allocate( hg123(gis:gie,gjs:gje,1:1) )
   allocate( hg11(gis:gie),hg12(gis:gie) )
   allocate( hg21(gjs:gje),hg22(gjs:gje) )

!$omp parallel
!$omp do private(i) schedule(static)
   do i = gis, gie
    hg11(i) = 2d0*(x1(i)+dx1(i  ))/x1(i)*idx1(i+1)/sum(dx1(i:i+1))
    hg12(i) = 2d0*(x1(i)-dx1(i+1))/x1(i)*idx1(i  )/sum(dx1(i:i+1))
   end do
!$omp end do
!$omp do private(j) schedule(static)
   do j = gjs, gje
    hg21(j) = ( dx2(j  )/tan(x2(j))+2d0)*idx2(j+1)/sum(dx2(j:j+1))
    hg22(j) = (-dx2(j+1)/tan(x2(j))+2d0)*idx2(j  )/sum(dx2(j:j+1))
   end do
!$omp end do
   k = ks
!$omp do private(i,j) schedule(static) collapse(2)
   do j = gjs, gje
    do i = gis, gie
     hg123(i,j,k) = ( 2d0*(dx1(i+1)-dx1(i)-x1(i))*idx1(i)*idx1(i+1) &
                    + ((dx2(j+1)-dx2(j))/tan(x2(j))-2d0) &
                      *idx2(j)*idx2(j+1)/x1(i) ) &
                  / x1(i)
     if(i==gis) hg123(i,j,k) = hg123(i,j,k) + hg12(i)
     if(j==gjs) hg123(i,j,k) = hg123(i,j,k) + hg22(j)/x1(i)**2
     if(j==gje) hg123(i,j,k) = hg123(i,j,k) + hg21(j)/x1(i)**2
     if(gbtype==1.and.i==gie)hg123(i,j,k) = hg123(i,j,k) + hg11(i)*x1(i)/x1(i+1)
    end do
   end do
!$omp end do
!$omp end parallel
   hg12(gis) = 0d0
   if(gbtype==1)hg11(gie) = 0d0
   hg21(gje) = 0d0
   hg22(gjs) = 0d0
   
   hg_dx = huge
   do j = js, je
    do i = is, ie
     hg_dx = min(hg_dx,dxi1(i)*x1(i)*dxi2(j)/(dxi1(i)+x1(i)*dxi2(j)))
    end do
   end do

! for 3D spherical %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.dim==3)then
! Normal discretization
   allocate( hg11(gis:gie),hg12(gis:gie) )
   allocate( hg21(gjs:gje),hg22(gjs:gje) )
   allocate( hg31(gks:gke),hg32(gks:gke) )
   allocate( hg123, mold=hgsrc )

!$omp parallel
!$omp do private(i)
   do i = gis, gie
    hg11(i) = 2d0*(x1(i)+dx1(i  ))/x1(i)*idx1(i+1)/sum(dx1(i:i+1))
    hg12(i) = 2d0*(x1(i)-dx1(i+1))/x1(i)*idx1(i  )/sum(dx1(i:i+1))
   end do
!$omp end do
!$omp do private(j)
   do j = gjs, gje
    hg21(j) = ( dx2(j  )/tan(x2(j))+2d0)*idx2(j+1)/sum(dx2(j:j+1))
    hg22(j) = (-dx2(j+1)/tan(x2(j))+2d0)*idx2(j  )/sum(dx2(j:j+1))
   end do
!$omp end do
!$omp do private(k)
   do k = gks, gke
    hg31(k) = 2d0*idx3(k+1)/sum(dx3(k:k+1))
    hg32(k) = 2d0*idx3(k  )/sum(dx3(k:k+1))
   end do
!$omp end do

!$omp do private(i,j,k) collapse(3)
   do k = gks, gke
    do j = gjs, gje
     do i = gis, gie
      hg123(i,j,k) = ( 2d0*(dx1(i+1)-dx1(i)-x1(i))*idx1(i)*idx1(i+1) &
                     + ((dx2(j+1)-dx2(j))/tan(x2(j))-2d0) &
                       *idx2(j)*idx2(j+1)/x1(i)  &
                     - 2d0*idx3(k)*idx3(k+1)/x1(i)/sinc(j)**2 ) &
                   / x1(i)
     if(i==gis) hg123(i,j,k) = hg123(i,j,k) + hg12(i)
     if(j==gjs) hg123(i,j,k) = hg123(i,j,k) + hg22(j)/x1(i)**2
     if(j==gje) hg123(i,j,k) = hg123(i,j,k) + hg21(j)/x1(i)**2
     if(gbtype==1.and.i==gie)hg123(i,j,k) = hg123(i,j,k) + hg11(i)*x1(i)/x1(i+1)
     end do
    end do
   end do
!$omp end do
!$omp end parallel
   hg12(gis) = 0d0
   if(gbtype==1)hg11(gie) = 0d0
   hg21(gje) = 0d0
   hg22(gjs) = 0d0

   hg_dx = huge
   do k = ks, ke
    do j = js, je
     do i = is, ie
      hg_dx = min(hg_dx,dxi1(i)*x1(i)*dxi2(j)*x1(i)*dxi3(k) &
                        / ( dxi1(i)*x1(i)*dxi2(j) &
                          + x1(i)*dxi2(j)*x1(i)*dxi3(k)  &
                          + x1(i)*dxi3(k)*dxi1(i) ) &
                 )
     end do
    end do
   end do

  else
   print *,'coefficients for Hyperbolic gravity under construction (gravsetup)'
   stop
  end if

! call gravbound
  if(tn==0)dt_old = dt / (courant*HGfac) * hgcfl
  hgsrc = 0d0

 end if

 return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                             SUBROUTINE MIC
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Modified Incomplete Cholesky Decomposition (1,2)

 subroutine mic

  use grid,only:ie
  use gravmod

  implicit none

  real(8):: micpara
  integer:: gin

!----------------------------------------------------------------------------

  gin = gie-gis+1

  micpara = 0.95d0 ! Should be set <1

! l=1
  preca(1) = a1(1)
  precb(1) = a2(1)
  precc(1) = a3(1)
  prece(1) = 0d0
  precd(1) = 1d0/preca(1)

  do l=2,gie-gis!ie-1
   preca(l) = a1(l) - precb(l-1)**2 *precd(l-1) &
            - micpara*precb(l-1)*prece(l-1)*precd(l-1)
   precb(l) = a2(l)
   precc(l) = a3(l)
   prece(l) =-a3(l-1)*precb(l-1)*precd(l-1)
   precd(l) = 1d0/preca(l) 
  end do

! l=gin
   preca(gin) = a1(gin) - precb(gin-1)**2 *precd(gin-1) &
                        - prece(    1)**2 *precd(    1) &
              - micpara*precb(gin-1)*prece(gin-1)*precd(gin-1)
   precb(gin) = a2(gin) - a3(1)*prece(1)*precd(1)
   precc(gin) = a3(gin)
   prece(gin) =-a3(gin-1)*precb(gin-1)*precd(gin-1)
   precd(gin) = 1d0/preca(gin) 

  do l=gin+1,lmax
   preca(l) = a1(l) - precb(l-1)**2 *precd(l-1) - a3(l-gin)**2 *precd(l-gin) &
                    - prece(l-gin+1)**2 *precd(l-gin+1)&
                    - micpara*(precb(l-1)     *prece(l-1)     *precd(l-1)&
                           + precb(l-gin+1)*prece(l-gin+1)*precd(l-gin+1))
   precb(l) = a2(l) - a3(l-gin+1)*prece(l-gin+1)*precd(l-gin+1)
   precc(l) = a3(l)
   prece(l) =-a3(l-1)*precb(l-1)*precd(l-1)
   precd(l) = 1d0/preca(l)
  end do

return
end subroutine mic

end subroutine gravsetup


end module gravity_mod
