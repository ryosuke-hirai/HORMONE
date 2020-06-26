!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE GRAVITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate gravitational forces.

subroutine gravity

  use settings,only:eq_sym,courant
  use grid
  use constants
  use physval
  use gravmod

  implicit none

  integer flgcg
  real*8 rr, rrold, pAp, alpha, beta, phih, cgrav2, dtgrav, mind
  real*8,dimension(1:lmax):: gsrc, pp, absrob
  real*8,allocatable,dimension(:,:,:):: newphi!, fric
  real*8,allocatable,dimension(:):: intphi

!-----------------------------------------------------------------------------
gin = gie - gis + 1

if(gravswitch==0.or.gravswitch==1)then
 grvphi = 0.d0
elseif(gravswitch==2.or.(gravswitch==3.and.tn==0))then
!elseif(gravswitch==2)then
! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 call gravbound

! cylindrical (equatorial+axial symmetry) ####################################
 if(je==1.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
  mind = minval(d(is:ie,js:je,ks:ke))
  do l = 1, lmax
   i = modlimax(l) +gis-1
   k = (l-modlimax(l))/gin + gks
   x(l) = grvphi(i,js,k)
   gsrc(l) = 4d0*pi*G*mind*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*5d-1)
   if(i>=is)then;if(i<=ie)then;if(k>=ks)then;if(k<=ke)then
    gsrc(l) = 4d0*pi*G*d(i,js,k)*x1(i)*dxi1(i)*((dx3(k)+dx3(k+1))*5d-1)
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

   do l = 1, lmax
    x(l) = x(l) + alpha*pp(l) !x_k+1
    r(l) = r(l) - alpha*aw(l) !r_k+1
   end do

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

   do l = 1, lmax
    pp(l) = z(l) + beta * pp(l)
   end do

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

!  gphidt = 0d0
!  phicg     = grvphi
  grvphiold = grvphi

 elseif(je==1.and.crdnt==1.and.dim==2)then ! for cylindrical coordinates
! convert x to phi
  do k = ks, ke
   do i = is, ie
    phicg(i,js,k) = x(i-is+1+(k-ks)*ie)
   end do
  end do

  do k = ks,ke
   do j = js,je
    phicg(is-2,j,k) = phicg(is+1,j,k)
    phicg(is-1,j,k) = phicg(is  ,j,k)
   end do
  end do
  phicg(is:ie,js:je,ks-2) = phicg(is:ie,js:je,ks)
  phicg(is:ie,js:je,ks-1) = phicg(is:ie,js:je,ks)

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

  grvphiold = grvphi

 end if
endif
if(gravswitch==3.and.tn/=0)then
! Hyperbolic Self-Gravity $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(crdnt==1.and.je==1)then
! Cartoon mesh method for axially symmetric cylindrical coordinates
  allocate( newphi(gis:gie,js:je,gks:gke), intphi(1:4) )

!  fric = 0.5d0*d(is:ie,js:je,ks:ke)! temp

!  dtgrav = dt / (courant*HGfac) * hgcfl

  cgrav2 = HGfac*max(maxval(cf(is:ie,js,ks:ke)+abs(v1(is:ie,js,ks:ke))), &
                     0d0*maxval(cf(is:ie,js,ks:ke)+abs(v2(is:ie,js,ks:ke))), &
                     maxval(cf(is:ie,js,ks:ke)+abs(v3(is:ie,js,ks:ke))) )
  dtgrav = min(minval(dxi1(gis:gie)),minval(dxi3(gks:gke)))/cgrav2
  dtgrav = dtgrav*hgcfl
  cgrav2 = cgrav2*cgrav2

  hgsrc(is:ie,js:je,ks:ke) = d(is:ie,js:je,ks:ke)
!  do n = 1, int(HGfac)
  do while (grvtime<time+dt)

   grvphi(gis-1,js,gks:gke) = grvphi(gis,js,gks:gke)
   if(eq_sym)grvphi(gis:gie,js,gks-1) = grvphi(gis:gie,js,gks)

!$omp parallel
!$omp do private(k)
   do k = gks, gke
    grvphi(gie+1,js,k) = orgdis(gie,js,k)/orgdis(gie+1,js,k) * grvphi(gie,js,k)
   end do
!$omp end do
!$omp do private(i)
   do i = gis, gie
    grvphi(i,js,gke+1) = orgdis(i,js,gke)/orgdis(i,js,gke+1) * grvphi(i,js,gke)
    grvphi(i,js,gks-1) = orgdis(i,js,gks)/orgdis(i,js,gks-1) * grvphi(i,js,gks)
   end do
!$omp end do
!!$!$omp do private(i,j,k,phih)
!!$   do k = gks, gke
!!$    do j = js, je
!!$     do i = gis, gie
!!$      phih = sum( lag(-1:1,i)*grvphi(i-1:i+1,j,k) )
!!$      newphi(i,j,k) = cgrav2*dtgrav*(dtgrav+dt_old)* &
!!$                    ( hg11(i)*grvphi(i+1,j,k) + hg12(i)*grvphi(i-1,j,k) + &
!!$                      hg31(k)*grvphi(i,j,k+1) + hg32(k)*grvphi(i,j,k-1) + &
!!$                      hg21(i)*phih - hg123(i,j,k)*grvphi(i,j,k) &
!!$                    - 2d0*pi*G*hgsrc(i,j,k) ) & ! source term
!!$                    - dtgrav/dt_old*grvphiold(i,j,k) &
!!$                    + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
!!$     end do
!!$    end do
!!$   end do
!!$!$omp end do

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
                    ( (sum(intphi(1:4))+2d0*phih-6d0*grvphi(i,j,k))/(h*h) &
                      - 4d0*pi*G*hgsrc(i,j,k) ) & ! source term
                    - dtgrav/dt_old*grvphiold(i,j,k) &
                    + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
     end do
    end do
   end do
!$omp end do
!$omp single
   grvphiold = grvphi
   grvphi(gis:gie,js:je,gks:gke) = newphi(gis:gie,js:je,gks:gke)
!$omp end single
!$omp do private(k)
   do k = gks, gke
    grvphi(gie+1,js,k) = orgdis(gie,js,k)/orgdis(gie+1,js,k) * newphi(gie,js,k)
   end do
!$omp end do
!$omp do private(i)
   do i = gis, gie
    grvphi(i,js,gke+1) = orgdis(i,js,gke)/orgdis(i,js,gke+1) * newphi(i,js,gke)
    grvphi(i,js,gks-1) = orgdis(i,js,gks)/orgdis(i,js,gks-1) * newphi(i,js,gks)
   end do
!$omp end do
!$omp end parallel

   grvphi(gis-1,js:je,gks:gke) = grvphi(gis,js:je,gks:gke)

   if(eq_sym)then ! for plane symmetry
    grvphi(is:gie,js:je,ks-1) = newphi(is:gie,js:je,ks)
    grvphi(is:gie,js:je,ks-2) = newphi(is:gie,js:je,ks+1)
   end if
   dt_old = dtgrav
   grvtime = grvtime + dtgrav

   if(grvtime>=time+dt)then
    exit
   end if


  end do

  deallocate(newphi,intphi)

 elseif(crdnt==2.and.ke==1)then
! Axisymmetric spherical coordinates
  allocate( newphi(gis:gie,gjs:gje,gks:gke) )

  cgrav2 = HGfac*max(maxval(cf(is:ie,js:je,ks)+abs(v1(is:ie,js:je,ks))), &
                     maxval(cf(is:ie,js:je,ks)+abs(v2(is:ie,js:je,ks))) )

  dtgrav = min(minval(dxi1(gis:gie)),x1(is)*minval(dxi2(gjs:gje)))/cgrav2
  dtgrav = dtgrav*hgcfl
  cgrav2 = cgrav2*cgrav2

  hgsrc(is:ie,js:je,ks:ke) = d(is:ie,js:je,ks:ke)

!  do n = 1, int(HGfac)
  do while (grvtime<time+dt)
!$omp parallel
   k = gks
!$omp workshare
   grvphi(gis-1,gjs:gje,k) = grvphi(gis,gjs:gje,k)
   grvphi(gis:gie,gjs-1,k) = grvphi(gis:gie,gjs,k)
   grvphi(gis:gie,gje+1,k) = grvphi(gis:gie,gje,k)
!$omp end workshare

!$omp do private(j)
   do j = gjs, gje
    grvphi(gie+1,j,k) = x1(gie)/x1(gie+1) * grvphi(gie,j,k)
   end do
!$omp end do

!$omp do private(i,j)
   do j = gjs, gje
    do i = gis, gie
     newphi(i,j,k) = 0.5d0*cgrav2*dtgrav*(dtgrav+dt_old)* &
                   ( hg11(i)*grvphi(i+1,j,k) + hg12(i)*grvphi(i-1,j,k) + &
                    (hg21(j)*grvphi(i,j+1,k) + hg22(j)*grvphi(i,j-1,k)) &
                     / (x1(i)*x1(i)) + &
                     hg123(i,j,k)*grvphi(i,j,k) &
                   - 4d0*pi*G*hgsrc(i,j,k) ) & ! source term
                   - dtgrav/dt_old*grvphiold(i,j,k) &
                   + (1d0+dtgrav/dt_old)*grvphi(i,j,k)
    end do
   end do
!$omp end do
!$omp workshare
   grvphiold(gis:gie,js:je,gks:gke) = grvphi(gis:gie,js:je,gks:gke)
   grvphi(gis:gie,js:je,gks:gke) = newphi(gis:gie,js:je,gks:gke)
!$omp end workshare
!$omp do private(j)
   do j = gjs, gje
    grvphi(gie+1,j,gks) = x1(gie)/x1(gie+1) * grvphi(gie,j,gks)
   end do
!$omp end do
!$omp end parallel

   dt_old = dtgrav
   grvtime = grvtime + dtgrav

!!$   if(grvtime>=time+dt)then
!!$    exit
!!$   end if

  end do

  grvphi(gis-1,gjs:gje,gks) = grvphi(gis,gjs:gje,gks)
  grvphi(gis:gie,gjs-1,gks) = grvphi(gis:gje,gjs,gks)
  grvphi(gis:gie,gje+1,gks) = grvphi(gis:gje,gje,gks)

  deallocate(newphi)

 end if

! dt_old = dt


end if


contains

! \\\\\\\\\\\\\\\\\\\\\\\\\\\\ SUBROUTINE Avec \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
   subroutine Avec(w)

     use gravmod

     implicit none

  ! Subroutine to calculate A*vector

     real*8,intent(in),dimension(1:lmax)::w

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

    implicit none

 ! subroutine to calculate (CC^T)^{-1}r for preconditioned CG method

    real*8,intent(in),dimension(1:lmax):: rrr

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

