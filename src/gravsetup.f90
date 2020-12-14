!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GRAVSETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up the matrix for Poisson equation

subroutine gravsetup
  
  use settings,only:eq_sym,courant
  use grid
  use gravmod
  
  implicit none

!-----------------------------------------------------------------------------
 gin = gie-gis+1

 if(gravswitch==2.or.gravswitch==3)then
! Extracting i and j from l
 do l = 1, lmax
!  modlimax(l) = mod(l,ie)
  modlimax(l) = mod(l,gin)
!  if(modlimax(l)==0)modlimax(l)=ie
  if(modlimax(l)==0) modlimax(l) = gie
 end do

! Constructing matrix A
 if(je==1.and.crdnt==1.and.dim==2)then
  
   do l = 1,lmax-gin ! calculating diagonal elements
    i = modlimax(l) +gis-1
    k = (l-modlimax(l))/gin + gks
    a1(l) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
              + 2.d0*x1(i)*dxi1(i) / (dx3(k)*dx3(k+1)) ) &
            * ( (dx3(k)+dx3(k+1)) / 2.d0 )
    a2(l) = xi1(i)*(dx3(k)+dx3(k+1))/2.d0/dx1(i+1)
    a3(l) = x1(i)*dxi1(i)/dx3(k+1)
!    if(mod(l,ie)==0)a2(l)=0.0d0
    if(mod(l,gin)==0)a2(l)=0.0d0
   end do

   do l=lmax-gin+1,lmax-1 ! k+1 line disappears
    i = modlimax(l) +gis-1
    k = (l-modlimax(l))/gin + gks
    a1(l) = -( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) &
             + 2.d0*x1(i)*dxi1(i)/(dx3(k)*dx3(k+1)) ) &
          * ( (dx3(k)+dx3(k+1))/2.d0 )
    a2(l) = xi1(i)*(dx3(k)+dx3(k+1))/2.d0/dx1(i+1)
    a3(l) = 0.0d0
    if(mod(l,gin)==0)a2(l)=0.0d0
   end do

    i = modlimax(lmax) +gis-1 ! i+1 line disappears
    k = (lmax-modlimax(lmax))/gin + gks
    a1(lmax) = -(xi1(i)/dx1(i+1)+xi1(i-1)/dx1(i) &
                +2.d0*x1(i)*dxi1(i)/(dx3(k)*dx3(k+1))) &
                *((dx3(k)+dx3(k+1))/2.d0)
    a2(lmax) = 0.0d0
    a3(lmax) = 0.0d0

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

 elseif(ke==1.and.crdnt==2.and.dim==2)then
  
   do l = 1,lmax-gin ! calculating diagonal elements
    i = modlimax(l) +gis-1
    j = (l-modlimax(l))/gin + gjs
    a1(l) = -( ( xi1(i)**2d0/dx1(i+1) + xi1(i-1)**2d0/dx1(i) ) &
            * sinc(j)*dxi2(j) &
            +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    a2(l) = xi1(i)**2d0 *sinc(j)*dxi2(j)/dx1(i+1)
    a3(l) = sini(j)*dxi1(i)/dx2(j+1)
    if(i==gis.and.xi1s>0d0)a1(l)=a1(l)+xi1(i-1)**2d0*sinc(j)*dxi2(j)/dx1(i)
    if(i==gie)a2(l)=0d0
   end do

   do l=lmax-gin+1,lmax-1 ! j+1 line disappears
    i = modlimax(l) +gis-1
    j = (l-modlimax(l))/gin + gjs
    a1(l) = -( ( xi1(i)**2d0/dx1(i+1) + xi1(i-1)**2d0/dx1(i) ) &
            * sinc(j)*dxi2(j) &
            +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    a2(l) = xi1(i)**2.d0 *sinc(j)*dxi2(j)/dx1(i+1)
    a3(l) = 0d0
    if(i==gis.and.xi1s>0d0)a1(l)=a1(l)+xi1(i-1)**2d0*sinc(j)*dxi2(j)/dx1(i)
    if(i==gie)a2(l)=0d0
   end do

    i = modlimax(lmax) +gis-1 ! i+1 line disappears
    j = (lmax-modlimax(lmax))/gin + gjs
    a1(lmax) = -( ( xi1(i)**2d0/dx1(i+1) + xi1(i-1)**2d0/dx1(i) ) &
               * sinc(j)*dxi2(j) &
               +( sini(j)/dx2(j+1) + sini(j-1)/dx2(j)) * dxi1(i) )
    a2(lmax) = 0d0
    a3(lmax) = 0d0

   call mic

 end if

! set initial x0

 if(tn==0.and.maxval(grvphi(gis:gie,gjs:gje,gks:gke))==0d0) grvphi = -1d3

 if(tn==0) dt_old = dt

end if

! For Hyperbolic gravity solver ----------------------------------------------
if(gravswitch==3)then
! for axisymmetric cylindrical
 if(je==1.and.dim==2.and.crdnt==1)then
! Cartoon mesh method
  gin = gie + 2
  allocate( hg11(gis-1:gin),hg12(gis-1:gin),hg21(gis-1:gin),hg22(gis-1:gin),&
            hg31(gks-1:gkn),hg32(gks-1:gkn),hg123(gis-1:gin,1:1,gks-1:gkn), &
            lag(-1:1,gis-1:gin),&
            orgdis(gis-1:gin,1:1,gks-1:gkn) )
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
  allocate( lag11(-1:1,gis:gie,js:je,gks:gke), &
            lag12(-1:1,gis:gie,js:je,gks:gke), &
            lag21(-1:1,gis:gie,js:je,gks:gke), &
            lag31(-1:1,gis:gie,js:je,gks:gke), &
            lag32(-1:1,gis:gie,js:je,gks:gke) )
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

 elseif(crdnt==2.and.ke==1.and.dim==2)then
! Normal discretization
  gin = gie + 2; gjn = gje + 2
  allocate( hg11(gis-1:gin),hg12(gis-1:gin),hg21(gjs-1:gjn),hg22(gjs-1:gjn),&
            hg123(gis-1:gin,gjs-1:gjn,1:1) )

  do i = gis-1, gie+1
   hg11(i) = 2d0*(x1(i)+dx1(i  ))/x1(i)*idx1(i+1)/sum(dx1(i:i+1))
   hg12(i) = 2d0*(x1(i)-dx1(i+1))/x1(i)*idx1(i  )/sum(dx1(i:i+1))
  end do
  do j = gjs-1, gje+1
   hg21(j) = ( dx2(j  )/tan(x2(j))+2d0)*idx2(j+1)/sum(dx2(j:j+1))
   hg22(j) = (-dx2(j+1)/tan(x2(j))+2d0)*idx2(j  )/sum(dx2(j:j+1))
  end do
  k = ks
  do j = gjs-1, gje+1
   do i = gis-1, gie+1
    hg123(i,j,k) = ( 2d0*(dx1(i+1)-dx1(i)-x1(i))*idx1(i)*idx1(i+1) &
                   + ((dx2(j+1)-dx2(j))/tan(x2(j))-2d0) &
                     *idx2(j)*idx2(j+1)/x1(i) ) &
                 / x1(i)
   end do
  end do

 else
  print *,'coefficients for Hyperbolic gravity under construction (gravsetup)'
  stop
 end if

! call gravbound
 if(tn==0)dt_old = dt / (courant*HGfac) * hgcfl

end if

hgsrc = 0d0

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

  real*8 micpara

!----------------------------------------------------------------------------

  gin = gie-gis+1

  micpara = 0.95d0 ! Should be set <1

! l=1
  preca(1) = a1(1)
  precb(1) = a2(1)
  precc(1) = a3(1)
  prece(1) = 0.0d0
  precd(1) = 1.0d0/preca(1)

  do l=2,gie-gis!ie-1
   preca(l) = a1(l) - precb(l-1)**2 *precd(l-1) &
            - micpara*precb(l-1)*prece(l-1)*precd(l-1)
   precb(l) = a2(l)
   precc(l) = a3(l)
   prece(l) =-a3(l-1)*precb(l-1)*precd(l-1)
   precd(l) = 1.0d0/preca(l) 
  end do

! l=gin
   preca(gin) = a1(gin) - precb(gin-1)**2 *precd(gin-1) &
                        - prece(    1)**2 *precd(    1) &
              - micpara*precb(gin-1)*prece(gin-1)*precd(gin-1)
   precb(gin) = a2(gin) - a3(1)*prece(1)*precd(1)
   precc(gin) = a3(gin)
   prece(gin) =-a3(gin-1)*precb(gin-1)*precd(gin-1)
   precd(gin) = 1.0d0/preca(gin) 

  do l=gin+1,lmax
   preca(l) = a1(l) - precb(l-1)**2 *precd(l-1) - a3(l-gin)**2 *precd(l-gin) &
                    - prece(l-gin+1)**2 *precd(l-gin+1)&
                    - micpara*(precb(l-1)     *prece(l-1)     *precd(l-1)&
                           + precb(l-gin+1)*prece(l-gin+1)*precd(l-gin+1))
   precb(l) = a2(l) - a3(l-gin+1)*prece(l-gin+1)*precd(l-gin+1)
   precc(l) = a3(l)
   prece(l) =-a3(l-1)*precb(l-1)*precd(l-1)
   precd(l) = 1.0d0/preca(l)
  end do

return
end subroutine mic

end subroutine gravsetup
