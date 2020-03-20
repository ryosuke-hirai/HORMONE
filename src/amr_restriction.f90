!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE AMR_RESTRICTION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To fix parent blocks so that they maintain conservation

subroutine amr_restriction(level)

 use amr_module
 use pressure_mod

 implicit none

 integer,intent(in):: level
 type(leaf_contents),dimension(1:cib):: clf

!-----------------------------------------------------------------------------

 do bkn = 1, totbloks
  if(bk(bkn)%level==level-1)then!.and.bk(bkn)%child(1)/=0)then
! Fix physical values using child blocks
   if(bk(bkn)%child(1)/=0)then
    do lid = 1, cib
     clf(lid) = lf(bk(bkn)%child(lid))
    end do
!print *,'before',bkn,sum(lf(bkn)%d(1:ib,1,1)*lf(bkn)%dvol(1:ib,1,1))
    call amr_restrict(lf(bkn),clf)
!print *,'after ',bkn,sum(lf(bkn)%d(1:ib,1,1)*lf(bkn)%dvol(1:ib,1,1))
! Fix coarse-fine interface cells to maintain conservation
   else
    do lfn = 1, face
     if(bk(bkn)%ldif(lfn)==1)then
      do lid = 1, cib
       clf(lid) = lf(bk(bk(bkn)%neigh(lfn))%child(lid))
      end do
!print *,bkn,lfn,clf(1)%ptot(:,1,1)
!      call amr_update(lfn,lf(bkn),clf)
!print *,'after',lf
     end if
    end do

   end if
  end if
 end do

return

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE AMR_RESTRICT
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To fix parent block using child blocks

subroutine amr_restrict(lf,clf)

 use settings,only:crdnt
 use grid,only:is,js,ks
 use amr_templates
 use amr_module,only:cib,ib,jb,kb,cuts,bkn

 implicit none

 integer i,j,k,i1,j1,k1,i2,j2,k2,lid, ii,jj,kk
 type(leaf_contents),intent(in),dimension(1:cib):: clf
 type(leaf_contents),intent(inout):: lf
 real*8 sintp, costp, sinph, cosph, sintc, costc
 real*8,allocatable,dimension(:,:,:,:):: vec

!-----------------------------------------------------------------------------

! Average conserved values ++++++++++++++++++++++++++++++++++++++++++++++++++
if(crdnt==0)then ! for Cartesian coordinates
 do k = ks, kb
  do j = js, jb
   do i = is, ib
    lid = (i-1)*cuts/ib + (j-1)*cuts/jb*cuts + (k-1)*cuts/kb*cuts*cuts + 1
    i1 = mod(i,max(1,ib/cuts))
    i1 = cuts*i1 + (1-sign(1,i1-1))/2*ib - cuts + 1
    i1 = max(  1,(1+sign(1,i1))/2 *  i1         + (1-sign(1,i1))/2 )
    i2 = min( ib,(1+sign(1,i1))/2 * (i1+cuts-1) + (1-sign(1,i1))/2 )
    j1 = mod(j,max(1,jb/cuts))
    j1 = cuts*j1 + (1-sign(1,j1-1))/2*jb - cuts + 1
    j1 = max(  1,(1+sign(1,j1))/2 *  j1         + (1-sign(1,j1))/2 )
    j2 = min( jb,(1+sign(1,j1))/2 * (j1+cuts-1) + (1-sign(1,j1))/2 )
    k1 = mod(k,max(1,kb/cuts))
    k1 = cuts*k1 + (1-sign(1,k1-1))/2*kb - cuts + 1
    k1 = max(  1,(1+sign(1,k1))/2 *  k1         + (1-sign(1,k1))/2 )
    k2 = min( kb,(1+sign(1,k1))/2 * (k1+cuts-1) + (1-sign(1,k1))/2 )

    lf% d(i,j,k) = sum( clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf%v1(i,j,k) = sum( clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k) * lf%d(i,j,k) )

    lf%v2(i,j,k) = sum( clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k) * lf%d(i,j,k) )

    lf%v3(i,j,k) = sum( clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k) * lf%d(i,j,k) )

    lf%b1(i,j,k) = sum( clf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf%b2(i,j,k) = sum( clf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf%b3(i,j,k) = sum( clf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf% e(i,j,k) = sum( clf(lid)%   e(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf%phi(i,j,k)= sum( clf(lid)% phi(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)
   end do
  end do
 end do

elseif(crdnt==1)then ! for cylindrical coordinates
  print *,"vecave for cylindrical is under construction"
  stop

elseif(crdnt==2)then ! for spherical coordinates
 do k = ks, kb
  do j = js, jb
   do i = is, ib
    lid = (i-1)*cuts/ib + (j-1)*cuts/jb*cuts + (k-1)*cuts/kb*cuts*cuts + 1
    i1 = mod(i,max(1,ib/cuts))
    i1 = cuts*i1 + (1-sign(1,i1-1))/2*ib - cuts + 1
    i1 = max(  1,(1+sign(1,i1))/2 *  i1         + (1-sign(1,i1))/2 )
    i2 = min( ib,(1+sign(1,i1))/2 * (i1+cuts-1) + (1-sign(1,i1))/2 )
    j1 = mod(j,max(1,jb/cuts))
    j1 = cuts*j1 + (1-sign(1,j1-1))/2*jb - cuts + 1
    j1 = max(  1,(1+sign(1,j1))/2 *  j1         + (1-sign(1,j1))/2 )
    j2 = min( jb,(1+sign(1,j1))/2 * (j1+cuts-1) + (1-sign(1,j1))/2 )
    k1 = mod(k,max(1,kb/cuts))
    k1 = cuts*k1 + (1-sign(1,k1-1))/2*kb - cuts + 1
    k1 = max(  1,(1+sign(1,k1))/2 *  k1         + (1-sign(1,k1))/2 )
    k2 = min( kb,(1+sign(1,k1))/2 * (k1+cuts-1) + (1-sign(1,k1))/2 )

    lf% d(i,j,k) = sum( clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf% e(i,j,k) = sum( clf(lid)%   e(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

    lf%phi(i,j,k)= sum( clf(lid)% phi(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / lf%dvol(i,j,k)

! transform bases for vector values
! tp : theta_parent, tc : theta_child, ph : phi_child - phi_parent
    sintp = sin(lf%x2(j)) ; costp = cos(lf%x2(j))
    allocate( vec(1:3,i1:i2,j1:j2,k1:k2) )
    do kk = k1, k2
     do jj = j1, j2
      do ii = i1, i2
       sintc = sin(clf(lid)%x2(jj)) ; costc = cos(clf(lid)%x2(jj))
       sinph = sin(clf(lid)%x3(kk)-lf%x3(k))
       cosph = cos(clf(lid)%x3(kk)-lf%x3(k))

       vec(1,ii,jj,kk) = clf(lid)%v1(ii,jj,kk) &
                       * ( costp*costc + sintp*sintc*cosph ) &
                       + clf(lid)%v2(ii,jj,kk) &
                       * (-costp*sintc + sintp*costc*cosph ) &
                       + clf(lid)%v3(ii,jj,kk) * (-sintp*sinph)

       vec(2,ii,jj,kk) = clf(lid)%v1(ii,jj,kk) &
                       * (-sintp*costc + costp*sintc*cosph ) &
                       + clf(lid)%v2(ii,jj,kk) &
                       * ( sintp*sintc + costp*costc*cosph ) &
                       + clf(lid)%v3(ii,jj,kk) * (-costp*sinph)

       vec(3,ii,jj,kk) = clf(lid)%v1(ii,jj,kk) * sintc*sinph &
                       + clf(lid)%v2(ii,jj,kk) * costc*sinph &
                       + clf(lid)%v3(ii,jj,kk) * cosph
      end do
     end do
    end do

    lf%v1(i,j,k) = sum( vec(1,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )

    lf%v2(i,j,k) = sum( vec(2,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )

    lf%v3(i,j,k) = sum( vec(3,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )

    do kk = k1, k2
     do jj = j1, j2
      do ii = i1, i2
       sintc = sin(clf(lid)%x2(jj)) ; costc = cos(clf(lid)%x2(jj))
       sinph = sin(clf(lid)%x3(kk)-lf%x3(k))
       cosph = cos(clf(lid)%x3(kk)-lf%x3(k))

       vec(1,ii,jj,kk) = clf(lid)%b1(ii,jj,kk) &
                       * ( costp*costc + sintp*sintc*cosph ) &
                       + clf(lid)%b2(ii,jj,kk) &
                       * (-costp*sintc + sintp*costc*cosph ) &
                       + clf(lid)%b3(ii,jj,kk) * (-sintp*sinph)

       vec(2,ii,jj,kk) = clf(lid)%b1(ii,jj,kk) &
                       * (-sintp*costc + costp*sintc*cosph ) &
                       + clf(lid)%b2(ii,jj,kk) &
                       * ( sintp*sintc + costp*costc*cosph ) &
                       + clf(lid)%b3(ii,jj,kk) * (-costp*sinph)

       vec(3,ii,jj,kk) = clf(lid)%b1(ii,jj,kk) * sintc*sinph &
                       + clf(lid)%b2(ii,jj,kk) * costc*sinph &
                       + clf(lid)%b3(ii,jj,kk) * cosph
      end do
     end do
    end do

    lf%b1(i,j,k) = sum( vec(1,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )

    lf%b2(i,j,k) = sum( vec(2,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )

    lf%b3(i,j,k) = sum( vec(3,i1:i2,j1:j2,k1:k2) &
                      * clf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / ( lf%dvol(i,j,k)*lf%d(i,j,k) )
    deallocate( vec )

   end do
  end do
 end do

end if

 call pressure_block(lf)

return
end subroutine amr_restrict


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE AMR_UPDATE
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To update coarse/fine interface cells

subroutine amr_update(fn,lf,clf)

 use settings,only:crdnt
 use grid,only:is,js,ks
 use physval,only:gamma
 use amr_templates
 use amr_module,only:cib,ib,jb,kb,cuts,schedule,lvl

 implicit none

 integer,intent(in):: fn
 integer i,j,k,i1,j1,k1,i2,j2,k2,lid, ii,jj,kk, l,m
 type(leaf_contents),intent(inout),dimension(1:cib):: clf
 type(leaf_contents),intent(inout):: lf
 real*8,dimension(1:9):: newflux
 real*8 bsq, vsq, sintc,costc,sintp,costp,sinph,cosph
 real*8,allocatable,dimension(:,:,:):: vec

!-----------------------------------------------------------------------------

! for tflx numbering
 if(lvl==1)then
  l = 1 ; m = 2
 elseif(schedule(lvl-1)==schedule(lvl-2))then
  l = 3 ; m = 4
 else
  l = 1 ; m = 2
 end if

 if(fn == 1)then ! x1 inner boundary of coarse cell
  do k = ks, kb
   do j = js, jb
    lid = cuts*int((j-1)*cuts/jb+1) + cuts**2*int((k-1)*cuts/kb)
    j2  = mod(cuts*j,jb) ; j2 = j2 + (1-sign(1,j2-1))/2*jb
    j1  = max(1,j2 - cuts + 1)
    k2  = mod(cuts*k,kb) ; k2 = k2 + (1-sign(1,k2-1))/2*kb
    k1  = max(1,k2 - cuts + 1)

    newflux = 0.d0 ; if(crdnt/=0)allocate( vec(1:3,j1:j2,k1:k2) )
    do kk = k1, k2 ; do jj = j1, j2
     ! transform bases for vector fluxes -------------------------------------
     if(crdnt==1)then ! cylindrical coordinates
      print *,'flux update for cylindrical coordinates is under construction'
      stop
     elseif(crdnt==2)then ! spherical coordinates
      sintp = sin(lf%x2(j)) ; costp = cos(lf%x2(j))
      sintc = sin(clf(lid)%x2(jj)) ; costc = cos(clf(lid)%x2(jj))
      sinph = sin(clf(lid)%x3(kk)-lf%x3(k))
      cosph = cos(clf(lid)%x3(kk)-lf%x3(k))
      ! momentum fluxes
      vec(1,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,2) + clf(lid)%tflx1(4,jj,kk,2) ) &
                   * ( costp*costc + sintp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,3) + clf(lid)%tflx1(4,jj,kk,3) ) &
                   * (-costp*sintc + sintp*costc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,4) + clf(lid)%tflx1(4,jj,kk,4) ) &
                   * (-sintp*sinph)

      vec(2,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,2) + clf(lid)%tflx1(4,jj,kk,2) ) &
                   * (-sintp*costc + costp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,3) + clf(lid)%tflx1(4,jj,kk,3) ) &
                   * ( sintp*sintc + costp*costc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,4) + clf(lid)%tflx1(4,jj,kk,4) ) &
                   * (-costp*sinph)

      vec(3,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,2) + clf(lid)%tflx1(4,jj,kk,2) ) &
                   * sintc*sinph &
                   + ( clf(lid)%tflx1(2,jj,kk,3) + clf(lid)%tflx1(4,jj,kk,3) ) &
                   * costc*sinph &
                   + ( clf(lid)%tflx1(2,jj,kk,4) + clf(lid)%tflx1(4,jj,kk,4) ) &
                   * cosph

!      clf(lid)%tflx1(m,jj,kk,2:4) = vec(1:3,jj,kk)
      newflux(2:4) = newflux(2:4) &
                   + clf(lid)%idetg1(ib)*clf(lid)%detg1(ib) * vec(1:3,jj,kk)

      ! magnetic field fluxes
      vec(1,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,5) + clf(lid)%tflx1(4,jj,kk,5) ) &
                   * ( costp*costc + sintp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,6) + clf(lid)%tflx1(4,jj,kk,6) ) &
                   * (-costp*sintc + sintp*costc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,7) + clf(lid)%tflx1(4,jj,kk,7) ) &
                   * (-sintp*sinph)

      vec(2,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,5) + clf(lid)%tflx1(4,jj,kk,5) ) &
                   * (-sintp*costc + costp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,6) + clf(lid)%tflx1(4,jj,kk,6) ) &
                   * ( sintp*sintc + costp*costc*cosph ) &
                   + ( clf(lid)%tflx1(2,jj,kk,7) + clf(lid)%tflx1(4,jj,kk,7) ) &
                   * (-costp*sinph)

      vec(3,jj,kk) = ( clf(lid)%tflx1(2,jj,kk,5) + clf(lid)%tflx1(4,jj,kk,5) ) &
                   * sintc*sinph &
                   + ( clf(lid)%tflx1(2,jj,kk,6) + clf(lid)%tflx1(4,jj,kk,6) ) &
                   * costc*sinph &
                   + ( clf(lid)%tflx1(2,jj,kk,7) + clf(lid)%tflx1(4,jj,kk,7) ) &
                   * cosph

!      clf(lid)%tflx1(m,jj,kk,5:7) = vec(1:3,jj,kk)
      newflux(5:7) = newflux(5:7) &
                   + clf(lid)%idetg1(ib)*clf(lid)%detg1(ib) * vec(1:3,jj,kk)
     end if
     !------------------------------------------------------------------------
     newflux(1) = newflux(1) &
           + clf(lid)%dvol(ib,jj,kk) &
           * clf(lid)%idetg1(ib)*clf(lid)%detg1(ib) &
           * ( clf(lid)%tflx1(2,jj,kk,1) + clf(lid)%tflx1(4,jj,kk,1) )

     newflux(8:9) = newflux(8:9) &
           + clf(lid)%dvol(ib,jj,kk) &
           * clf(lid)%idetg1(ib)*clf(lid)%detg1(ib) &
           * ( clf(lid)%tflx1(2,jj,kk,8:9) + clf(lid)%tflx1(4,jj,kk,8:9) )

    end do; end do
    deallocate(vec)
    newflux = newflux * dt_local*5.d-1 / lf%dvol(is,j,k)

    lf%u(is,j,k,:) = lf%u(is,j,k,:) &
                 - dt_local*lf%idetg1(is)*lf%detg1(is-1)*lf%tflx1(l,j,k,:) &
                 + newflux(1:9)!*dt_local

    !temporary
!    lf%u(is,j,k,2) = lf%u(is,j,k,2) - lf%ptot(is,j,k)*lf%sx1(is)

    lf% d(is,j,k) = lf%u(is,j,k,1)
    lf%v1(is,j,k) = lf%u(is,j,k,2) / lf%u(is,j,k,1)
    lf%v2(is,j,k) = lf%u(is,j,k,3) / lf%u(is,j,k,1)
    lf%v3(is,j,k) = lf%u(is,j,k,4) / lf%u(is,j,k,1)
    lf%b1(is,j,k) = lf%u(is,j,k,5)
    lf%b2(is,j,k) = lf%u(is,j,k,6)
    lf%b3(is,j,k) = lf%u(is,j,k,7)
    lf% e(is,j,k) = lf%u(is,j,k,8)
    lf%phi(is,j,k)= lf%u(is,j,k,9)
! pressure and total pressure
    vsq = lf%v1(is,j,k)*lf%v1(is,j,k) + lf%v2(is,j,k)*lf%v2(is,j,k) &
        + lf%v3(is,j,k)*lf%v3(is,j,k)
    bsq = lf%b1(is,j,k)*lf%b1(is,j,k) + lf%b2(is,j,k)*lf%b2(is,j,k) &
        + lf%b3(is,j,k)*lf%b3(is,j,k)
    lf%p(is,j,k) = (gamma-1.d0) &
                 * (lf%e(is,j,k) - 0.5d0*lf%d(is,j,k)*vsq - 0.5d0 * bsq )
    lf%ptot(is,j,k) = lf%p(is,j,k) + 0.5d0 * bsq
    !temporary
    lf%u(is,j,k,2) = lf%u(is,j,k,2) + lf%ptot(is,j,k)*lf%sx1(is)
   end do
  end do
 elseif(fn == 2)then ! x1 outer boundary of coarse cell
  do k = ks, kb
   do j = js, jb
    lid = cuts*int((j-1)*cuts/jb+1) + cuts**2*int((k-1)*cuts/kb) - 1
    j2  = mod(cuts*j,jb) ; j2 = j2 + (1-sign(1,j2-1))/2*jb
    j1  = max(1,j2 - cuts + 1)
    k2  = mod(cuts*k,kb) ; k2 = k2 + (1-sign(1,k2-1))/2*kb
    k1  = max(1,k2 - cuts + 1)

    newflux = 0.d0 ; if(crdnt/=0)allocate( vec(1:3,j1:j2,k1:k2) )
    do kk = k1, k2 ; do jj = j1, j2
     ! transform bases for vector fluxes -------------------------------------
     if(crdnt==1)then ! cylindrical coordinates
      print *,'flux update for cylindrical coordinates is under construction'
      stop
     elseif(crdnt==2)then ! spherical coordinates
      sintp = sin(lf%x2(j)) ; costp = cos(lf%x2(j))
      sintc = sin(clf(lid)%x2(jj)) ; costc = cos(clf(lid)%x2(jj))
      sinph = sin(clf(lid)%x3(kk)-lf%x3(k))
      cosph = cos(clf(lid)%x3(kk)-lf%x3(k))
      ! momentum fluxes
      vec(1,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,2) + clf(lid)%tflx1(3,jj,kk,2) ) &
                   * ( costp*costc + sintp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,3) + clf(lid)%tflx1(3,jj,kk,3) ) &
                   * (-costp*sintc + sintp*costc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,4) + clf(lid)%tflx1(3,jj,kk,4) ) &
                   * (-sintp*sinph)

      vec(2,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,2) + clf(lid)%tflx1(3,jj,kk,2) ) &
                   * (-sintp*costc + costp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,3) + clf(lid)%tflx1(3,jj,kk,3) ) &
                   * ( sintp*sintc + costp*costc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,4) + clf(lid)%tflx1(3,jj,kk,4) ) &
                   * (-costp*sinph)

      vec(3,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,2) + clf(lid)%tflx1(3,jj,kk,2) ) &
                   * sintc*sinph &
                   + ( clf(lid)%tflx1(1,jj,kk,3) + clf(lid)%tflx1(3,jj,kk,3) ) &
                   * costc*sinph &
                   + ( clf(lid)%tflx1(1,jj,kk,4) + clf(lid)%tflx1(3,jj,kk,4) ) &
                   * cosph

!      clf(lid)%tflx1(1,jj,kk,2:4) = vec(1:3,jj,kk)
      newflux(2:4) = newflux(2:4) &
                   + clf(lid)%dvol(is,jj,kk) &
                   * clf(lid)%idetg1(is)*clf(lid)%detg1(is-1) * vec(1:3,jj,kk)

      ! magnetic field fluxes
      vec(1,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,5) + clf(lid)%tflx1(3,jj,kk,5) ) &
                   * ( costp*costc + sintp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,6) + clf(lid)%tflx1(3,jj,kk,6) ) &
                   * (-costp*sintc + sintp*costc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,7) + clf(lid)%tflx1(3,jj,kk,7) ) &
                   * (-sintp*sinph)

      vec(2,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,5) + clf(lid)%tflx1(3,jj,kk,5) ) &
                   * (-sintp*costc + costp*sintc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,6) + clf(lid)%tflx1(3,jj,kk,6) ) &
                   * ( sintp*sintc + costp*costc*cosph ) &
                   + ( clf(lid)%tflx1(1,jj,kk,7) + clf(lid)%tflx1(3,jj,kk,7) ) &
                   * (-costp*sinph)

      vec(3,jj,kk) = ( clf(lid)%tflx1(1,jj,kk,5) + clf(lid)%tflx1(3,jj,kk,5) ) &
                   * sintc*sinph &
                   + ( clf(lid)%tflx1(1,jj,kk,6) + clf(lid)%tflx1(3,jj,kk,6) ) &
                   * costc*sinph &
                   + ( clf(lid)%tflx1(1,jj,kk,7) + clf(lid)%tflx1(3,jj,kk,7) ) &
                   * cosph

!      clf(lid)%tflx1(1,jj,kk,5:7) = vec(1:3,jj,kk)
      newflux(5:7) = newflux(5:7) &
                   + clf(lid)%dvol(is,jj,kk) &
                   * clf(lid)%idetg1(is)*clf(lid)%detg1(is-1) * vec(1:3,jj,kk)
     end if
     !-----------------------------------------------------------------------
     newflux(1) = newflux(1) &
          + clf(lid)%dvol(is,jj,kk) &
          * clf(lid)%idetg1(is)*clf(lid)%detg1(is-1)&
          * ( clf(lid)%tflx1(1,jj,kk,1) + clf(lid)%tflx1(3,jj,kk,1) )

     newflux(8:9) = newflux(8:9) &
          + clf(lid)%dvol(is,jj,kk) &
          * clf(lid)%idetg1(is)*clf(lid)%detg1(is-1)&
          * ( clf(lid)%tflx1(1,jj,kk,8:9) + clf(lid)%tflx1(3,jj,kk,8:9) )
    end do ; end do

    deallocate(vec)
    newflux = newflux * dt_local*5.d-1 / lf%dvol(ib,j,k)

    lf%u(ib,j,k,:) = lf%u(ib,j,k,:) &
                 + dt_local*lf%idetg1(ib)*lf%detg1(ib)*lf%tflx1(m,j,k,:) &
                 - newflux! * dt_local

    lf% d(ib,j,k) = lf%u(ib,j,k,1)
    lf%v1(ib,j,k) = lf%u(ib,j,k,2) / lf%u(ib,j,k,1)
    lf%v2(ib,j,k) = lf%u(ib,j,k,3) / lf%u(ib,j,k,1)
    lf%v3(ib,j,k) = lf%u(ib,j,k,4) / lf%u(ib,j,k,1)
    lf%b1(ib,j,k) = lf%u(ib,j,k,5)
    lf%b2(ib,j,k) = lf%u(ib,j,k,6)
    lf%b3(ib,j,k) = lf%u(ib,j,k,7)
    lf% e(ib,j,k) = lf%u(ib,j,k,8)
    lf%phi(ib,j,k)= lf%u(ib,j,k,9)
! pressure and total pressure
    vsq = lf%v1(ib,j,k)*lf%v1(ib,j,k) + lf%v2(ib,j,k)*lf%v2(ib,j,k) &
        + lf%v3(ib,j,k)*lf%v3(ib,j,k)
    bsq = lf%b1(ib,j,k)*lf%b1(ib,j,k) + lf%b2(ib,j,k)*lf%b2(ib,j,k) &
        + lf%b3(ib,j,k)*lf%b3(ib,j,k)
    lf%p(ib,j,k) = (gamma-1.d0) &
                 * (lf%e(ib,j,k) - 0.5d0*lf%d(ib,j,k)*vsq - 0.5d0 * bsq )
    lf%ptot(ib,j,k) = lf%p(ib,j,k) + 0.5d0 * bsq
   end do
  end do
 elseif(fn == 3)then ! x2 inner boundary of coarse cell
  do k = ks, kb
   do i = is, ib
    lid = cuts + 1 + int((i-1)*cuts/ib) + cuts*cuts*int((k-1)*cuts/kb)
    i2  = mod(cuts*i,ib) ; i2 = i2 + (1-sign(1,i2-1))/2*ib
    i1  = max(1,i2 - cuts + 1)
    k2  = mod(cuts*k,kb) ; k2 = k2 + (1-sign(1,k2-1))/2*kb
    k1  = max(1,k2 - cuts + 1)

    newflux = 0.d0
    do kk = k1, k2 ; do ii = i1, i2
     newflux = newflux &
          + clf(lid)%dvol(ii,jb,kk) &
          * clf(lid)%idetg2(ii,jb)*clf(lid)%detg2(ii,jb)*clf(lid)%tflx2(ii,2,kk,:)
    end do ; end do
    newflux = newflux / lf%dvol(i,js,k) * dt_local*5.d-1

    lf%u(i,js,k,:) = lf%u(i,js,k,:) &
             - dt_local*lf%idetg2(i,js)*lf%detg2(i,js-1)*lf%tflx2(i,1,k,:) &
             + newflux

    lf% d(i,js,k) = lf%u(i,js,k,1)
    lf%v1(i,js,k) = lf%u(i,js,k,2) / lf%u(i,js,k,1)
    lf%v2(i,js,k) = lf%u(i,js,k,3) / lf%u(i,js,k,1)
    lf%v3(i,js,k) = lf%u(i,js,k,4) / lf%u(i,js,k,1)
    lf%b1(i,js,k) = lf%u(i,js,k,5)
    lf%b2(i,js,k) = lf%u(i,js,k,6)
    lf%b3(i,js,k) = lf%u(i,js,k,7)
    lf% e(i,js,k) = lf%u(i,js,k,8)
    lf%phi(i,js,k)= lf%u(i,js,k,9)
! pressure and total pressure
    vsq = lf%v1(i,js,k)*lf%v1(i,js,k) + lf%v2(i,js,k)*lf%v2(i,js,k) &
        + lf%v3(i,js,k)*lf%v3(i,js,k)
    bsq = lf%b1(i,js,k)*lf%b1(i,js,k) + lf%b2(i,js,k)*lf%b2(i,js,k) &
        + lf%b3(i,js,k)*lf%b3(i,js,k)
    lf%p(i,js,k) = (gamma-1.d0) &
                 * (lf%e(i,js,k) - 0.5d0*lf%d(i,js,k)*vsq - 0.5d0 * bsq )
    lf%ptot(i,js,k) = lf%p(i,js,k) + 0.5d0 * bsq
   end do
  end do
 elseif(fn == 4)then ! x2 outer boundary of coarse cell
  do k = ks, kb
   do i = is, ib
    lid = 1 + int((i-1)*cuts/ib) + cuts*cuts*int((k-1)*cuts/kb)
    i2  = mod(cuts*i,ib) ; i2 = i2 + (1-sign(1,i2-1))/2*ib
    i1  = max(1,i2 - cuts + 1)
    k2  = mod(cuts*k,kb) ; k2 = k2 + (1-sign(1,k2-1))/2*kb
    k1  = max(1,k2 - cuts + 1)

    newflux = 0.d0
    do kk = k1, k2 ; do ii = i1, i2
     newflux = newflux &
      + clf(lid)%dvol(ii,js,kk) &
      * clf(lid)%idetg2(ii,js)*clf(lid)%detg2(ii,js-1)*clf(lid)%tflx2(ii,1,kk,:)
    end do ; end do
    newflux = newflux / lf%dvol(i,jb,k) * dt_local*5.d-1

    lf%u(i,jb,k,:) = lf%u(i,jb,k,:) &
             + dt_local*lf%idetg2(i,jb)*lf%detg2(i,jb)*lf%tflx2(i,2,k,:) &
             - newflux

    lf% d(i,jb,k) = lf%u(i,jb,k,1)
    lf%v1(i,jb,k) = lf%u(i,jb,k,2) / lf%u(i,jb,k,1)
    lf%v2(i,jb,k) = lf%u(i,jb,k,3) / lf%u(i,jb,k,1)
    lf%v3(i,jb,k) = lf%u(i,jb,k,4) / lf%u(i,jb,k,1)
    lf%b1(i,jb,k) = lf%u(i,jb,k,5)
    lf%b2(i,jb,k) = lf%u(i,jb,k,6)
    lf%b3(i,jb,k) = lf%u(i,jb,k,7)
    lf% e(i,jb,k) = lf%u(i,jb,k,8)
    lf%phi(i,jb,k)= lf%u(i,jb,k,9)
! pressure and total pressure
    vsq = lf%v1(i,jb,k)*lf%v1(i,jb,k) + lf%v2(i,jb,k)*lf%v2(i,jb,k) &
        + lf%v3(i,jb,k)*lf%v3(i,jb,k)
    bsq = lf%b1(i,jb,k)*lf%b1(i,jb,k) + lf%b2(i,jb,k)*lf%b2(i,jb,k) &
        + lf%b3(i,jb,k)*lf%b3(i,jb,k)
    lf%p(i,jb,k) = (gamma-1.d0) &
                 * (lf%e(i,jb,k) - 0.5d0*lf%d(i,jb,k)*vsq - 0.5d0 * bsq )
    lf%ptot(i,jb,k) = lf%p(i,jb,k) + 0.5d0 * bsq
   end do
  end do
 elseif(fn == 5)then ! x3 inner boundary of coarse cell
  do j = js, jb
   do i = is, ib
    lid = cuts*cuts + 1 + int((i-1)*cuts/ib) + cuts*int((j-1)*cuts/jb)
    i2  = mod(cuts*i,ib) ; i2 = i2 + (1-sign(1,i2-1))/2*ib
    i1  = max(1,i2 - cuts + 1)
    j2  = mod(cuts*j,jb) ; j2 = j2 + (1-sign(1,j2-1))/2*jb
    j1  = max(1,j2 - cuts + 1)

    newflux = 0.d0
    do jj = j1, j2 ; do ii = i1, i2
     newflux = newflux &
             + clf(lid)%dvol(ii,jj,kb) &
             * clf(lid)%idetg3(ii,jj,kb)*clf(lid)%tflx3(ii,jj,2,:)
    end do ; end do
    newflux = newflux / lf%dvol(i,j,ks) * dt_local*5.d-1

    lf%u(i,j,ks,:) = lf%u(i,j,ks,:) &
         - dt_local*lf%idetg3(i,j,ks)*lf%tflx3(i,j,1,:) &
         + newflux

    lf% d(i,j,ks) = lf%u(i,j,ks,1)
    lf%v1(i,j,ks) = lf%u(i,j,ks,2) / lf%u(i,j,ks,1)
    lf%v2(i,j,ks) = lf%u(i,j,ks,3) / lf%u(i,j,ks,1)
    lf%v3(i,j,ks) = lf%u(i,j,ks,4) / lf%u(i,j,ks,1)
    lf%b1(i,j,ks) = lf%u(i,j,ks,5)
    lf%b2(i,j,ks) = lf%u(i,j,ks,6)
    lf%b3(i,j,ks) = lf%u(i,j,ks,7)
    lf% e(i,j,ks) = lf%u(i,j,ks,8)
    lf%phi(i,j,ks)= lf%u(i,j,ks,9)
! pressure and total pressure
    vsq = lf%v1(i,j,ks)*lf%v1(i,j,ks) + lf%v2(i,j,ks)*lf%v2(i,j,ks) &
        + lf%v3(i,j,ks)*lf%v3(i,j,ks)
    bsq = lf%b1(i,j,ks)*lf%b1(i,j,ks) + lf%b2(i,j,ks)*lf%b2(i,j,ks) &
        + lf%b3(i,j,ks)*lf%b3(i,j,ks)
    lf%p(i,j,ks) = (gamma-1.d0) &
                 * (lf%e(i,j,ks) - 0.5d0*lf%d(i,j,ks)*vsq - 0.5d0 * bsq )
    lf%ptot(i,j,ks) = lf%p(i,j,ks) + 0.5d0 * bsq
   end do
  end do
 elseif(fn == 6)then ! x3 outer boundary of coarse cell
  do j = js, jb
   do i = is, ib
  lid = 1 + int((i-1)*cuts/ib) + cuts*int((j-1)*cuts/jb)
    i2  = mod(cuts*i,ib) ; i2 = i2 + (1-sign(1,i2-1))/2*ib
    i1  = max(1,i2 - cuts + 1)
    j2  = mod(cuts*j,jb) ; j2 = j2 + (1-sign(1,j2-1))/2*jb
    j1  = max(1,j2 - cuts + 1)

    newflux = 0.d0
    do jj = j1, j2 ; do ii = i1, i2
     newflux = newflux &
             + clf(lid)%dvol(ii,jj,ks) &
             * clf(lid)%idetg3(ii,jj,ks)*clf(lid)%tflx3(ii,jj,1,:)
    end do ; end do
    newflux = newflux / lf%dvol(i,j,kb) * dt_local*5.d-1

    lf%u(i,j,kb,:) = lf%u(i,j,kb,:) &
         + dt_local*lf%idetg3(i,j,kb)*lf%tflx3(i,j,2,:) &
         - newflux

    lf% d(i,j,kb) = lf%u(i,j,kb,1)
    lf%v1(i,j,kb) = lf%u(i,j,kb,2) / lf%u(i,j,kb,1)
    lf%v2(i,j,kb) = lf%u(i,j,kb,3) / lf%u(i,j,kb,1)
    lf%v3(i,j,kb) = lf%u(i,j,kb,4) / lf%u(i,j,kb,1)
    lf%b1(i,j,kb) = lf%u(i,j,kb,5)
    lf%b2(i,j,kb) = lf%u(i,j,kb,6)
    lf%b3(i,j,kb) = lf%u(i,j,kb,7)
    lf% e(i,j,kb) = lf%u(i,j,kb,8)
    lf%phi(i,j,kb)= lf%u(i,j,kb,9)
! pressure and total pressure
    vsq = lf%v1(i,j,kb)*lf%v1(i,j,kb) + lf%v2(i,j,kb)*lf%v2(i,j,kb) &
        + lf%v3(i,j,kb)*lf%v3(i,j,kb)
    bsq = lf%b1(i,j,kb)*lf%b1(i,j,kb) + lf%b2(i,j,kb)*lf%b2(i,j,kb) &
        + lf%b3(i,j,kb)*lf%b3(i,j,kb)
    lf%p(i,j,kb) = (gamma-1.d0) &
                 * (lf%e(i,j,kb) - 0.5d0*lf%d(i,j,kb)*vsq - 0.5d0 * bsq )
    lf%ptot(i,j,kb) = lf%p(i,j,kb) + 0.5d0 * bsq
   end do
  end do
 else
  print *,"Error from face number",fn,"@amr_restriction.f"
  stop
 end if

return
end subroutine amr_update

end subroutine amr_restriction
