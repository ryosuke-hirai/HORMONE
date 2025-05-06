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

pure function get_raddim(in, jn, kn) result(dim)
! PURPOSE: Get the dimension of the radiation problem
  integer,intent(in)::in,jn,kn
  integer::dim

  if (in>1.and.jn>1.and.kn==1)then
    dim=3
  elseif(in>1.and.jn>1.and.kn==1)then
    dim=21
  elseif(in>1.and.jn==1.and.kn>1)then
    dim=22
  elseif(in==1.and.jn>1.and.kn>1)then
    dim=23
  elseif(in>1.and.jn==1.and.kn==1)then
    dim=11
  elseif(in==1.and.jn>1.and.kn==1)then
    dim=12
  elseif(in==1.and.jn==1.and.kn>1)then
    dim=13
  endif

end function get_raddim


end module matrix_utils