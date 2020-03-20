!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE EXTERNALFIELD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To add external gravitational field to the gravitational potential.

subroutine externalfield

 use grid
 use physval,only:d
 use gravmod

 implicit none

!-----------------------------------------------------------------------------

 do k = ks, ke
  do j = js, je
   do i = is, ie
    grv1(i,j,k) = grv1(i,j,k) &
                 -( (dx1(i  )*idx1(i+1)*extgrv(i+1,j,k)    &
                   - dx1(i+1)*idx1(i  )*extgrv(i-1,j,k) )  &
                    /sum(dx1(i:i+1)) &
                  + (dx1(i+1)-dx1(i))*idx1(i)*idx1(i+1)*extgrv(i,j,k) ) &
               * d(i,j,k)
    grv2(i,j,k) = grv2(i,j,k) &
                 -( (dx2(j  )*idx2(j+1)*extgrv(i,j+1,k)   &
                   - dx2(j+1)*idx2(j  )*extgrv(i,j-1,k) ) &
                    /sum(dx2(j:j+1)) &
                  + (dx2(j+1)-dx2(j))*idx2(j)*idx2(j+1)*extgrv(i,j,k) ) &
               * d(i,j,k)
    grv3(i,j,k) = grv3(i,j,k) &
                 -( (dx3(k  )*idx3(k+1)*extgrv(i,j,k+1) &
                   - dx3(k+1)*idx3(k  )*extgrv(i,j,k-1) ) &
                    /sum(dx3(k:k+1)) &
                  + (dx3(k+1)-dx3(k))*idx3(k)*idx3(k+1)*extgrv(i,j,k) ) &
               * d(i,j,k)
   end do
  end do
 end do

return
end subroutine externalfield
