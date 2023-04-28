module smear_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To smear out quantities for nested-grid-like feature

subroutine smear

 use settings,only:compswitch,spn
 use grid
 use physval
 use composition_mod

 implicit none

 integer:: jb,kb

!-----------------------------------------------------------------------------

 if(sphrn==0.or.dim==1)return

! Average out central cells in spherical coordinates
! -> This avoids severe Courant conditions at the centre.


 if(crdnt==2)then
! First average mass fractions
  if(compswitch>=2)then

   do i = is, is+sphrn-1
    do n = 1, spn
     spc(n,i,js:je,ks:ke) = sum( u(i,js:je,ks:ke,1)*spc(n,i,js:je,ks:ke) &
                                *dvol(i,js:je,ks:ke) ) &
                          / sum( u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke) )
    end do
   end do
   jb=15;kb=15
   if(ke==1)kb=0
   do i = is+sphrn, is+sphrn+trnsn16-1
    do n = 1, spn
     do k = ks, ke,16
      do j = js, je,16
       spc(n,i,j:j+jb,k:k+kb)=sum( u(i,j:j+jb,k:k+kb,1)*spc(n,i,j:j+jb,k:k+kb) &
                                  *dvol(i,j:j+jb,k:k+kb) ) &
                             / sum( u(i,j:j+jb,k:k+kb,1)*dvol(i,j:j+jb,k:k+kb) )
      end do
     end do
    end do
   end do
   jb=7;kb=7
   if(ke==1)kb=0
   do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
    do n = 1, spn
     do k = ks, ke,8
      do j = js, je,8
       spc(n,i,j:j+jb,k:k+kb) = sum( u(i,j:j+jb,k:k+kb,1)*spc(n,i,j:j+jb,k:k+kb) &
                                  *dvol(i,j:j+jb,k:k+kb) ) &
                            / sum( u(i,j:j+jb,k:k+kb,1)*dvol(i,j:j+jb,k:k+kb) )
      end do
     end do
    end do
   end do
   jb=3;kb=3
   if(ke==1)kb=0
   do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
    do n = 1, spn
     do k = ks, ke,4
      do j = js, je,4
       spc(n,i,j:j+jb,k:k+kb) = sum( u(i,j:j+jb,k:k+kb,1)*spc(n,i,j:j+jb,k:k+kb) &
                                  *dvol(i,j:j+jb,k:k+kb) ) &
                            / sum( u(i,j:j+jb,k:k+kb,1)*dvol(i,j:j+jb,k:k+kb) )
      end do
     end do
    end do
   end do
   jb=1;kb=1
   if(ke==1)kb=0
   do i = is+sphrn+trnsn16+trnsn8+trnsn4,is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
    do n = 1, spn
     do k = ks, ke,2
      do j = js, je,2
       spc(n,i,j:j+jb,k:k+kb) = sum( u(i,j:j+jb,k:k+kb,1)*spc(n,i,j:j+jb,k:k+kb) &
                                  *dvol(i,j:j+jb,k:k+kb) ) &
                            / sum( u(i,j:j+jb,k:k+kb,1)*dvol(i,j:j+jb,k:k+kb) )
      end do
     end do
    end do
   end do
  end if

! Then average other conserved quantities
  do i = is, is+sphrn-1
   call angular_smear(i,js,je,ks,ke)
  end do
  jb=15;kb=15
  if(ke==1)kb=0
  do i = is+sphrn, is+sphrn+trnsn16-1
   do k = ks, ke, 16
    do j = js, je, 16
     call angular_smear(i,j,j+jb,k,k+kb)
    end do
   end do
  end do
  jb=7;kb=7
  if(ke==1)kb=0
  do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
   do k = ks, ke, 8
    do j = js, je, 8
     call angular_smear(i,j,j+jb,k,k+kb)
    end do
   end do
  end do
  jb=3;kb=3
  if(ke==1)kb=0
  do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
   do k = ks, ke, 4
    do j = js, je, 4
     call angular_smear(i,j,j+jb,k,k+kb)
    end do
   end do
  end do
  jb=1;kb=1
  if(ke==1)kb=0
  do i = is+sphrn+trnsn16+trnsn8+trnsn4, is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
   do k = ks, ke, 2
    do j = js, je, 2
     call angular_smear(i,j,j+jb,k,k+kb)
    end do
   end do
  end do
  
!!$  do ufn = 1, 9
!!$   do i = is, is+sphrn-1
!!$    u(i,js:je,ks:ke,ufn) = sum( u(i,js:je,ks:ke,ufn)*dvol(i,js:je,ks:ke) ) &
!!$                         / sum( dvol(i,js:je,ks:ke) )
!!$   end do
!!$   jb=15;kb=15
!!$   if(ke==1)kb=0
!!$   do i = is+sphrn, is+sphrn+trnsn16-1
!!$    do k = ks, ke, 16
!!$     do j = js, je, 16
!!$      u(i,j:j+jb,k:k+kb,ufn) = sum( u(i,j:j+jb,k:k+kb,ufn)*dvol(i,j:j+jb,k:k+kb) ) & 
!!$                            / sum( dvol(i,j:j+jb,k:k+kb) )
!!$     end do
!!$    end do
!!$   end do
!!$   jb=7;kb=7
!!$   if(ke==1)kb=0
!!$   do i = is+sphrn+trnsn16, is+sphrn+trnsn16+trnsn8-1
!!$    do k = ks, ke,8
!!$     do j = js, je,8
!!$      u(i,j:j+jb,k:k+kb,ufn) = sum( u(i,j:j+jb,k:k+kb,ufn)*dvol(i,j:j+jb,k:k+kb) ) &
!!$                           / sum( dvol(i,j:j+jb,k:k+kb) )
!!$     end do
!!$    end do
!!$   end do
!!$   jb=3;kb=3
!!$   if(ke==1)kb=0
!!$   do i = is+sphrn+trnsn16+trnsn8, is+sphrn+trnsn16+trnsn8+trnsn4-1
!!$    do k = ks, ke,4
!!$     do j = js, je,4
!!$      u(i,j:j+jb,k:k+kb,ufn) = sum( u(i,j:j+jb,k:k+kb,ufn)*dvol(i,j:j+jb,k:k+kb) ) &
!!$                           / sum( dvol(i,j:j+jb,k:k+kb) )
!!$     end do
!!$    end do
!!$   end do
!!$   jb=1;kb=1
!!$   if(ke==1)kb=0
!!$   do i = is+sphrn+trnsn16+trnsn8+trnsn4,is+sphrn+trnsn16+trnsn8+trnsn4+trnsn2-1
!!$    do k = ks, ke,2
!!$     do j = js, je,2
!!$      u(i,j:j+jb,k:k+kb,ufn) = sum( u(i,j:j+jb,k:k+kb,ufn)*dvol(i,j:j+jb,k:k+kb) ) &
!!$                           / sum( dvol(i,j:j+jb,k:k+kb) )
!!$     end do
!!$    end do
!!$   end do
!!$  end do
  
 end if
 
 return
end subroutine smear


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE ANGULAR_SMEAR
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Average out quantities over several cells in the angular direction

 subroutine angular_smear(i,js,je,ks,ke)

  use grid,only:x1,x2,x3,dvol
  use physval,only:u,v1,v2,v3
  use utils
  use gravmod,only:grvphi,extgrv,gravswitch,include_extgrv
  

  implicit none

  integer,intent(in)::i,js,je,ks,ke
  integer:: j,k
  real(8):: mtot, etot
  real(8),dimension(1:3):: vcar, xcar, momtot, vave
  
!-----------------------------------------------------------------------------

  momtot=0d0;etot=0d0
  do j = js, je
   do k = ks, ke
    xcar = polcar([x1(i),x2(j),x3(k)])
    call get_vcar(xcar,x3(k),u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),vcar)
    momtot = momtot + vcar*dvol(i,j,k)! add up momenta
    etot = etot + u(i,j,k,8)*dvol(i,j,k)! add up energy
    if(gravswitch>0)then
     etot = etot + u(i,j,k,1)*grvphi(i,j,k)*dvol(i,j,k)! and gravitational ene
    end if
    if(include_extgrv)then
     etot = etot + u(i,j,k,1)*extgrv(i,j,k)*dvol(i,j,k)! and external gravity
    end if
   end do
  end do
  mtot = sum( u(i,js:je,ks:ke,1)*dvol(i,js:je,ks:ke) )
  vave = momtot/mtot ! get average cartesian velocity
  u(i,js:je,ks:ke,1) = mtot/sum( dvol(i,js:je,ks:ke) ) ! density

  do j = js, je
   do k = ks, ke
    xcar = polcar([x1(i),x2(j),x3(k)])
    call get_vpol(xcar,x3(k),vave,v1(i,j,k),v2(i,j,k),v3(i,j,k))
    u(i,j,k,2) = v1(i,j,k)*u(i,j,k,1)
    u(i,j,k,3) = v2(i,j,k)*u(i,j,k,1)
    u(i,j,k,4) = v3(i,j,k)*u(i,j,k,1)
    if(gravswitch>0)then
     etot = etot - u(i,j,k,1)*grvphi(i,j,k)*dvol(i,j,k)
    end if
    if(include_extgrv)then
     etot = etot - u(i,j,k,1)*extgrv(i,j,k)*dvol(i,j,k)
    end if
   end do
  end do

  u(i,js:je,ks:ke,8) = etot / sum( dvol(i,js:je,ks:ke) )
!  u(i,js:je,ks:ke,8) = sum( u(i,js:je,ks:ke,8)*dvol(i,js:je,ks:ke) )&
!                      / sum( dvol(i,js:je,ks:ke) )
  
 return
end subroutine angular_smear


end module smear_mod
