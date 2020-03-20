module amr_derefinement_mod

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE AMR_DEREFINEMENT
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To derefine leaves into a larger leaf with restrictions.

subroutine amr_derefinement(lf,olf)

 use grid,only:is,js,ks
 use amr_templates
 use amr_module,only:cib,ib,jb,kb,cuts
 use pressure_mod
 use amr_interpolation_mod

 implicit none

 integer i,j,k,i1,j1,k1,i2,j2,k2,lid
 type(leaf_contents),intent(in),dimension(1:cib):: olf
 type(leaf_contents),intent(out):: lf

!-----------------------------------------------------------------------------

! Set coordinates ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 lf%xi1(is-1) = olf(1)%xi1(is-1) ; lf%xi1(ib) = olf(cib)%xi1(ib)
 lf%xi2(js-1) = olf(1)%xi2(js-1) ; lf%xi2(jb) = olf(cib)%xi2(jb)
 lf%xi3(ks-1) = olf(1)%xi3(ks-1) ; lf%xi3(kb) = olf(cib)%xi3(kb)

 lf%dx1 = ( lf%xi1(ib)-lf%xi1(is-1) ) / dble(ib) ; lf%idx1 = 1.d0 / lf%dx1
 lf%dx2 = ( lf%xi2(jb)-lf%xi2(js-1) ) / dble(jb) ; lf%idx2 = 1.d0 / lf%dx2
 lf%dx3 = ( lf%xi3(kb)-lf%xi3(ks-1) ) / dble(kb) ; lf%idx3 = 1.d0 / lf%dx3

 lf%xi1(is-2) = lf%xi1(is-1) - lf%dx1(is-1)
 lf%x1(is-1)  = lf%xi1(is-1) - 5.d-1 * lf%dx1(is-1)
 lf%x1(is-2)  = lf%x1(is-1)  - lf%dx1(is-2)

 lf%xi2(js-2) = lf%xi2(js-1) - lf%dx2(js-1)
 lf%x2(js-1)  = lf%xi2(js-1) - 5.d-1 * lf%dx2(js-1)
 lf%x2(js-2)  = lf%x2(js-1)  - lf%dx2(js-2)

 lf%xi3(ks-2) = lf%xi3(ks-1) - lf%dx3(ks-1)
 lf%x3(ks-1)  = lf%xi3(ks-1) - 5.d-1 * lf%dx3(ks-1)
 lf%x3(ks-2)  = lf%x3(ks-1)  - lf%dx3(ks-2)

 do i = is, ib+2
  lf%x1(i) = lf%x1(i-1) + lf%dx1(i) ; lf%xi1(i) = lf%xi1(i-1) + lf%dx1(i)
 end do

 do j = js, jb+2
  lf%x2(j) = lf%x2(j-1) + lf%dx2(j) ; lf%xi2(j) = lf%xi2(j-1) + lf%dx2(j)
 end do

 do k = ks, kb+2
  lf%x3(k) = lf%x3(k-1) + lf%dx3(k) ; lf%xi3(k) = lf%xi3(k-1) + lf%dx3(k)
 end do

! Set metrics +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 call amr_block_metric(lf)

! Average conserved values ++++++++++++++++++++++++++++++++++++++++++++++++++

 do k = ks, kb
  do j = js, jb
   do i = is, ib
    lid = (i-1)*cuts/ib + (j-1)*cuts/jb*cuts + (k-1)*cuts/kb*cuts*cuts + 1
    i1 = mod(i,max(1,ib/cuts))
    i1 = cuts*i1 + (1-sign(1,i1-1))/2*ib - cuts + 1
    i1 = (1+sign(1,i1))/2 *  i1         + (1-sign(1,i1))/2
    i2 = (1+sign(1,i1))/2 * (i1+cuts-1) + (1-sign(1,i1))/2
    j1 = mod(j,max(1,jb/cuts))
    j1 = cuts*j1 + (1-sign(1,j1-1))/2*jb - cuts + 1
    j1 = (1+sign(1,j1))/2 *  j1         + (1-sign(1,j1))/2
    j2 = (1+sign(1,j1))/2 * (j1+cuts-1) + (1-sign(1,j1))/2
    k1 = mod(k,max(1,kb/cuts))
    k1 = cuts*k1 + (1-sign(1,k1-1))/2*kb - cuts + 1
    k1 = (1+sign(1,k1))/2 *  k1         + (1-sign(1,k1))/2
    k2 = (1+sign(1,k1))/2 * (k1+cuts-1) + (1-sign(1,k1))/2

    lf% d(i,j,k) = sum( olf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )

    lf%v1(i,j,k) = sum( olf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%  v1(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) / lf%d(i,j,k)

    lf%v2(i,j,k) = sum( olf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%  v2(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) / lf%d(i,j,k)

    lf%v3(i,j,k) = sum( olf(lid)%   d(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%  v3(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) / lf%d(i,j,k)

    lf%b1(i,j,k) = sum( olf(lid)%  b1(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )

    lf%b2(i,j,k) = sum( olf(lid)%  b2(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )

    lf%b3(i,j,k) = sum( olf(lid)%  b3(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )

    lf% e(i,j,k) = sum( olf(lid)%   e(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )

    lf%phi(i,j,k)= sum( olf(lid)% phi(i1:i2,j1:j2,k1:k2)   &
                      * olf(lid)%dvol(i1:i2,j1:j2,k1:k2) ) &
                 / sum( olf(lid)%dvol(i1:i2,j1:j2,k1:k2) )
   end do
  end do
 end do

 call pressure_block(lf)

return
end subroutine amr_derefinement


end module amr_derefinement_mod
