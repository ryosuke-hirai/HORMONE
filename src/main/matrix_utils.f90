module matrix_utils
  implicit none

contains

pure function l_from_ijk(i,j,k,is,js,ks,in,jn) result(l)
! PURPOSE: Compute a single index l from i,j,k
! is,js,ks: Starting index for i,j,k directions
! in,jn,kn: Number of grid points in i,j,k directions
 integer,intent(in)::i,j,k,is,js,ks,in,jn
 integer::l
 l = i-is+1 + in*(j-js) + in*jn*(k-ks)
end function l_from_ijk

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE IJK_FROM_L
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get i,j,k from l
! is,js,ks: Starting index for i,j,k directions
! in,jn,kn: Number of grid points in i,j,k directions

pure subroutine ijk_from_l(l,is,js,ks,in,jn,i,j,k)

 integer,intent(in)::l,is,js,ks,in,jn
 integer,intent(out)::i,j,k

!-----------------------------------------------------------------------------

 i = mod(l,in)
 if(i==0)i=in
 j = mod(l-i,in*jn)/in
 k = (l-i-j*in)/(in*jn)

 i=is+i-1
 j=js+j
 k=ks+k

return
end subroutine ijk_from_l

end module matrix_utils