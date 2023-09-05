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
 use miccg_mod,only:cg=>cg_grv,miccg,ijk_from_l,l_from_ijk
 use omp_lib

 integer:: l, flgcg, gin, gjn, gkn, tngrav, grungen, jb, kb
 real(8):: phih, cgrav2, dtgrav, mind, h
 real(8),allocatable,dimension(:,:,:):: lapphi,newphi
 real(8),allocatable,dimension(:):: intphi
 real(8),allocatable,dimension(:):: x,y,z,r,gsrc
 real(8):: faco, facn, fact, vol

!-----------------------------------------------------------------------------

 if(gravswitch==0.or.gravswitch==1)return
 
 wtime(igrv) = wtime(igrv) - omp_get_wtime()
 
 gin = gie - gis + 1
 gjn = gje - gjs + 1
 gkn = gke - gks + 1

 if(gravswitch==2.or.(gravswitch==3.and.tn==0.and.dim==2))then
  allocate( x(1:cg%lmax), gsrc(1:cg%lmax) )
 
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
!$omp parallel do private(i,j,k,l)
   do l = 1, cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
    x(l) = grvphi(i,j,k)
    gsrc(l) = 4d0*pi*G*mind*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*0.5d0)
    if(i>=is)then;if(i<=ie)then;if(k>=ks)then;if(k<=ke)then
     gsrc(l) = 4d0*pi*G*d(i,j,k)*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*0.5d0)
    end if;end if;end if;end if
    if(k==gks) gsrc(l) = gsrc(l) - x1 (i)*dxi1(i)*idx3(k  )*phi3i(i,k-1)
    if(i==gie) gsrc(l) = gsrc(l) - xi1(i)*dxi3(k)*idx1(i+1)*phi1o(i+1,k)
    if(k==gke) gsrc(l) = gsrc(l) - x1 (i)*dxi1(i)*idx3(k+1)*phi3o(i,k+1)
   end do
!$omp end parallel do

! spherical (axial symmetry) #################################################
  elseif(crdnt==2.and.dim==2)then
print*,'kita'
! calculating b for Ax=b
   mind = minval(d(is:ie,js:je,ks:ke))
!$omp parallel do private(i,j,k,l)
   do l=1,cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
    x(l) = grvphi(i,j,k)
    gsrc(l) = 4d0*pi*G*mind*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)
    if(i>=is)then;if(i<=ie)then;if(j>=js)then;if(j<=je)then
     gsrc(l) = 4d0*pi*G*d(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)
    end if;end if;end if;end if
    if(i==gie) gsrc(l)= gsrc(l) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)&
                                  *phiio(i+1,j)
!   if(i==gis) gsrc(l)= gsrc(l) - xi1(i-1)*xi1(i-1)*sinc(j)*dxi2(j)*idx1(i)&
!        *phiii(i-1,j)
   end do
!$omp end parallel do

! spherical (3D) #############################################################
  elseif(crdnt==2.and.dim==3)then

! calculating b for Ax=b
   mind = minval(d(is:ie,js:je,ks:ke))
!$omp parallel do private(i,j,k,l)
   do l=1,cg%lmax
    call ijk_from_l(l,cg%is,cg%js,cg%ks,cg%in,cg%jn,cg%kn,i,j,k)
    x(l) = grvphi(i,j,k)
    gsrc(l) = 4d0*pi*G*mind*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)*dxi3(k)
    if(i>=is)then;if(i<=ie)then;if(j>=js)then;if(j<=je)then
     gsrc(l) = 4d0*pi*G*d(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)*dxi3(k)
    end if;end if;end if;end if
    if(i==gie) phiio(i+1,j) = -G*1.989d33/x1(gie+1)
    if(i==gie) gsrc(l)= gsrc(l) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)*dxi3(k)&
                                  *phiio(i+1,j)
!   if(i==gis) gsrc(l)= gsrc(l) - xi1(i-1)*xi1(i-1)*sinc(j)*dxi2(j)*idx1(i)&
!        *phiii(i-1,j)
   end do
!$omp end parallel do
  end if

  call miccg(cg,gsrc,x)

!-------------------------------------------------------------------------

! convert x to phi
!$omp parallel do private(i,j,k,l) collapse(3)
  do k = gks, gke
   do j = gjs, gje
    do i = gis, gie
     l = l_from_ijk(i,j,k,gis,gjs,gks,gie-gis+1,gje-gjs+1,gke-gks+1)
     grvphi(i,j,k) = x(l)
    end do
   end do
  end do
!$omp end parallel do


  if(je==1.and.crdnt==1.and.dim==2)then ! for cylindrical coordinates

   do k = gks,gke
    do j = js,je
     grvphi(is-2,j,k) = grvphi(is+1,j,k)
     grvphi(is-1,j,k) = grvphi(is  ,j,k)
    end do
   end do
   grvphi(gis:gie,js:je,gks-2) = grvphi(gis:gie,js:je,gks)
   grvphi(gis:gie,js:je,gks-1) = grvphi(gis:gie,js:je,gks)

  elseif(ke==1.and.crdnt==2.and.dim==2)then ! for spherical coordinates (2D)

   grvphi(is-1:gie+1,js-2,ks) = grvphi(is-1:gie+1,js+1,ks)
   grvphi(is-1:gie+1,js-1,ks) = grvphi(is-1:gie+1,js,ks)
   grvphi(is-1:gie+1,je+1,ks) = grvphi(is-1:gie+1,je,ks)
   grvphi(is-1:gie+1,je+2,ks) = grvphi(is-1:gie+1,je-1,ks)

   grvphi(is-1,:,:) = grvphi(is,:,:)
   grvphi(is-2,:,:) = grvphi(is+1,:,:)
   grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                   ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)

  elseif(crdnt==2.and.dim==3)then ! for spherical coordinates (3D)

   grvphi(is-1:gie+1,js-2,ks:ke) = grvphi(is-1:gie+1,js,ks:ke)
   grvphi(is-1:gie+1,js-1,ks:ke) = grvphi(is-1:gie+1,js,ks:ke)
   grvphi(is-1:gie+1,je+1,ks:ke) = grvphi(is-1:gie+1,je,ks:ke)
   grvphi(is-1:gie+1,je+2,ks:ke) = grvphi(is-1:gie+1,je,ks:ke)

   grvphi(is-1,:,:) = grvphi(is,:,:)
   grvphi(is-2,:,:) = grvphi(is+1,:,:)
   grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                   ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)

  end if

  if(gravswitch==3)grvphiold = grvphi
  
 endif


if(gravswitch==3.and.tn/=0)then
! Hyperbolic Self-Gravity $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(crdnt==1.and.je==js)then
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
! Boundary conditions
!$omp do private(j,k) collapse(2)
    do k = ks, ke
     do j = js, je
      grvphi(is-1,j,k) = grvphi(is,j,k)
      grvphi(ie+1,j,k) = grvphi(ie,j,k)*x1(ie)/x1(ie+1)
     end do
    end do
!$omp end do nowait
!$omp do private(i,k) collapse(2)
    do k = ks, ke
     do i = is, ie
      grvphi(i,js-1,k) = grvphi(i,js,k)
      grvphi(i,je+1,k) = grvphi(i,je,k)
     end do
    end do
!$omp end do nowait
!$omp do private(i,j) collapse(2)
    do j = js, je
     do i = is, ie
      grvphi(i,j,ks-1) = grvphi(i,j,ke)
      grvphi(i,j,ke+1) = grvphi(i,j,ks)
     end do
    end do
!$omp end do

  do l = 1, tngrav
   do grungen = 1, 3
! First set flux and source term
!$omp do private (i,j,k) collapse(3)
    do k = ks-1, ke
     do j = js-1, je
      do i = is-1, ie
       grv1(i,j,k) = -cgrav2*(grvphi(i+1,j,k)-grvphi(i,j,k))*idx1(i+1)
       grv2(i,j,k) = -cgrav2*(grvphi(i,j+1,k)-grvphi(i,j,k))*idx2(j+1)/x1(i)
       grv3(i,j,k) = -cgrav2*(grvphi(i,j,k+1)-grvphi(i,j,k))*idx3(k+1)/x1(i)/sin(x2(j))

       if(i<=is+sum(fmr_lvl(1:fmr_max))-1)then
        if(i<=is+fmr_lvl(1)-1)then
         grv2(i,j,k) = 0d0;grv3(i,j,k) = 0d0
        else
         fmr_loop: do n = 2, fmr_max
          if(i<=is+sum(fmr_lvl(1:n))-1)then
           grv2(i,j,k) = grv2(i,j,k)/dble(2**(fmr_max-n+1))
           exit fmr_loop
          end if
         end do fmr_loop
        end if
       end if
       if(i==is-1.or.j==js-1.or.k==ks-1)cycle
       hgsrc(i,j,k) = -cgrav2*4d0*pi*G*d(i,j,k)
      end do
     end do
    end do
!$omp end do

!$omp single
    grk3_number: select case (grungen)
    case(1) grk3_number
     faco = 1d0 ; fact = 1d0 ; facn = 0d0
    case(2) grk3_number
     faco = 0.75d0 ; fact = 0.25d0 ; facn = fact
    case(3) grk3_number
     faco = 1d0/3d0 ; fact = 2d0/3d0 ; facn = fact
    end select grk3_number
!$omp end single

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

    do n = 2, fmr_max
     if(fmr_lvl(n)==0)cycle
     jb=min(2**(fmr_max-n+1),je)-1 ; kb=min(2**(fmr_max-n+1),ke)-1
!$omp do private(i,j,k,vol) collapse(3)
     do k = ks, ke, kb+1
      do j = js, je, jb+1
       do i = is+sum(fmr_lvl(0:n-1)), is+sum(fmr_lvl(0:n))-1
        vol = sum(dvol(i,j:j+jb,k:k+kb))
        grvphi(i,j:j+jb,k:k+kb) = sum(grvphi(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/vol
        grvphidot(i,j:j+jb,k:k+kb) = sum(grvphidot(i,j:j+jb,k:k+kb)*dvol(i,j:j+jb,k:k+kb))/vol
       end do
      end do
     end do
!$omp end do
    end do
!$omp do private(i,j,k) collapse(2)
    do k = ks, ke
     do j = js, je
      do i = fmr_lvl(1)-1, is, -1
       grvphi(i,j,k) = grvphi(i+1,j,k) &
                     - G*(mc(i)-mc(is-1))/xi1(i)**2*dx1(i+1)
      end do
     end do
    end do
!$omp end do
   
! Boundary conditions
!$omp do private(j,k) collapse(2)
    do k = ks, ke
     do j = js, je
      grvphi(is-1,j,k) = grvphi(is,j,k)
      grvphi(ie+1,j,k) = grvphi(ie,j,k)*x1(ie)/x1(ie+1)
     end do
    end do
!$omp end do nowait
!$omp do private(i,k) collapse(2)
    do k = ks, ke
     do i = is, ie
      grvphi(i,js-1,k) = grvphi(i,js,k)
      grvphi(i,je+1,k) = grvphi(i,je,k)
     end do
    end do
!$omp end do nowait
!$omp do private(i,j) collapse(2)
    do j = js, je
     do i = is, ie
      grvphi(i,j,ks-1) = grvphi(i,j,ke)
      grvphi(i,j,ke+1) = grvphi(i,j,ks)
     end do
    end do
!$omp end do nowait
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

end subroutine gravity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SETUP_GRVCG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the A matrix elements for Laplacian

subroutine setup_grvcg(is,ie,js,je,ks,ke,cg)

 use settings,only:crdnt,eq_sym
 use grid,only:xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,idx2,dx3,dxi3,idx3,sini,sinc
 use miccg_mod,only:cg_set,ijk_from_l,get_preconditioner

 integer,intent(in)::is,ie,js,je,ks,ke
 type(cg_set),intent(out):: cg
 integer:: in,jn,kn,lmax,dim,i,j,k,l

!-----------------------------------------------------------------------------

 cg%is=is;cg%ie=ie; cg%js=js;cg%je=je; cg%ks=ks;cg%ke=ke
 in = ie-is+1; jn = je-js+1; kn = ke-ks+1
 cg%in=in; cg%jn=jn; cg%kn=kn
 lmax = in*jn*kn; cg%lmax=lmax
 if(ie>is.and.je>js.and.ke>ks)then
  dim=3
 elseif(ie>is.and.(je>js.or.ke>ks))then
  dim=2
 elseif(ie>is.and.je==js.and.ke==ks.and.crdnt==2)then
  dim=1
 else
  print*,'Error in setup_grvcg, dimension is not supported; dim=',dim
 end if

 select case(dim)
 case(1) ! 1D
! 1D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==1)then
   cg%Adiags = 2
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
    cg%A(1,l) = -( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) )
    cg%A(2,l) = xi1(i)**2/dx1(i+1)
    if(i==is.and.xi1s>0d0)cg%A(1,l)=cg%A(1,l)+xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(i==ie)cg%A(2,l) = 0d0

   end do

! Pre-conditioner matrix with MICCG(1,1) method
   cg%cdiags = 2
   allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
   cg%ic(1) = 0
   cg%ic(2) = 1
   cg%alpha = 0.99d0

  end if

 case(2) ! 2D
! 2D cylindrical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==1.and.je==js.and.ke>ks)then

   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

! Equation matrix
   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
    cg%A(1,l) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
                      + 2d0*x1(i)*dxi1(i) / (dx3(k)*dx3(k+1)) ) &
                    * 0.5d0*sum(dx3(k:k+1))
    cg%A(2,l) = 0.5d0*xi1(i)*sum(dx3(k:k+1))*idx1(i+1)
    cg%A(3,l) = x1(i)*dxi1(i)/dx3(k+1)
    if(i==ie)cg%A(2,l) = 0d0
    if(k==ke)cg%A(3,l) = 0d0
    
    if(eq_sym.and.k==ks)then! for Neumann boundary at bc3i (equatorial symmetry)
     cg%A(1,l) = cg%A(1,l) + x1(i)*dxi1(i)/dx3(k+1)
    end if
   
   end do

! 2D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(crdnt==2.and.ie>is.and.je>js.and.ke==ks)then

   cg%Adiags = 3
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in

   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
! Equation matrix
    cg%A(1,l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                  * sinc(j)*dxi2(j) &
               + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    cg%A(2,l) = xi1(i)**2 *sinc(j)*dxi2(j)/dx1(i+1)
    cg%A(3,l) = sini(j)*dxi1(i)/dx2(j+1)
    if(i==is.and.xi1s>0d0)cg%A(1,l)=cg%A(1,l)+xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(i==ie)cg%A(2,l) = 0d0
    if(j==je)then ! Reflection boundary at axis or equatorial plane (if eq_sym)
     cg%A(1,l) = cg%A(1,l) + sini(j)*dxi1(i)*idx2(j+1)
     cg%A(3,l) = 0d0
    end if

   end do

  end if

! Pre-conditioner matrix with MICCG(1,2) method
  cg%cdiags = 4
  allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
  cg%ic(1) = 0
  cg%ic(2) = 1
  cg%ic(3) = in-1
  cg%ic(4) = in
  cg%alpha = 0.99d0
  

 case(3) ! 3D 
! 3D spherical coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(crdnt==2)then

   cg%Adiags = 5
   allocate(cg%ia(1:cg%Adiags),cg%A(1:cg%Adiags,1:lmax))
   cg%ia(1) = 0
   cg%ia(2) = 1
   cg%ia(3) = in
   cg%ia(4) = in*jn
   cg%ia(5) = in*jn*(kn-1)

   do l = 1, lmax
    call ijk_from_l(l,is,js,ks,in,jn,kn,i,j,k)
    cg%A(1,l) = -( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) &
                  * sinc(j)*dxi2(j)*dxi3(k) &
                  + ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i)*dxi3(k) &
                  + ( idx3(k)+idx3(k-1) )*dxi1(i)*dxi2(j) )
    cg%A(2,l) = xi1(i)**2 *sinc(j)*dxi2(j)*dxi3(k)/dx1(i+1)
    cg%A(3,l) = sini(j)*dxi1(i)*dxi3(k)/dx2(j+1)
    cg%A(4,l) = dxi1(i)*dxi2(j)/dx3(k+1)
    cg%A(5,l) = dxi1(i)*dxi2(j)/dx3(k  )
    if(eq_sym.and.j==je)then ! for Neumann boundary at bc2o (equatorial symmetry)
     cg%A(1,l) = cg%A(1,l) + cg%A(3,l)
    end if
    if(i==is.and.xi1s>0d0)cg%A(1,l) = cg%A(1,l) &
                                     + xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)
    if(i==ie)cg%A(2,l) = 0d0
    if(j==je)cg%A(3,l) = 0d0
    if(k==ke)cg%A(4,l) = 0d0
    if(k/=ks)cg%A(5,l) = 0d0

   end do

! Pre-conditioner matrix with ICCG(1,1) method
   cg%cdiags = 5
   allocate(cg%ic(1:cg%cdiags),cg%c(1:cg%cdiags,1:lmax))
   cg%ic(1) = 0
   cg%ic(2) = 1
   cg%ic(3) = in
   cg%ic(4) = in*jn
   cg%ic(5) = in*jn*(kn-1)
   cg%alpha = 0.999d0 ! No modification
  end if

 end select

 call get_preconditioner(cg)

return
end subroutine setup_grvcg

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
 use miccg_mod,only:cg_grv

 integer:: l,gin
 real(8):: h

!-----------------------------------------------------------------------------

! set initial x0

 if(tn==0.and.maxval(grvphi(gis:gie,gjs:gje,gks:gke))==0d0) grvphi = -1d3

 if(tn==0) dt_old = dt

 if(gravswitch==2.or.(gravswitch==3.and.tn==0))then
  call setup_grvcg(gis,gie,gjs,gje,gks,gke,cg_grv)
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

end subroutine gravsetup


end module gravity_mod
