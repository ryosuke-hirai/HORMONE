!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE HEADONCOLLISION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for head-on collision of stars.

subroutine headoncollision

 use funcs
 use grid
 use settings,only:dt_out
 use physval
 use constants
 use ejectamod
 use gravmod!,only:extgrv,grvtime,include_extgrv,rdis,grvphi

 implicit none

 real*8 mass, radius, coredis, mnow, mold, rnow, dr, dbg, shell, shelld
 real*8,allocatable,dimension(:):: dpol,ppol,rpol,mpol
 integer N, nn, nol
 real*8,allocatable,dimension(:,:):: dat

!-----------------------------------------------------------------------------

 dbg = 1d-10
 d(is:ie,js:je,ks:ke) = 0d0

 open(unit=1230,file='../stardata/s20-9Rsun.dat',status='old')
 read(1230,'()'); read(1230,'()')
 read(1230,*)i,nol
 allocate(dat(1:nol,1:23), dpol(1:nol), ppol(1:nol), rpol(1:nol), mpol(1:nol))
 read(1230,'()'); read(1230,'()'); read(1230,'()')
 do i = 1, nol
  read(1230,*) dat(i,1:23)
 end do
 mass = dat(1,2)*msun
 radius = dat(1,10)

 dpol(1:nol) = dat(nol:1:-1,12)
 ppol(1:nol) = dat(nol:1:-1,14)
 rpol(1:nol) = dat(nol:1:-1,10)
 mpol(1:nol) = dat(nol:1:-1, 2)*msun

 dr = dx1(is)*1.d0
 rnow = dx1(is) ; mold = 0d0
 do
  shell = 0d0
  ! Choose cells in current mass shell
  do k = ks, ke
   do i = is, ie
    if(rdis(i,k)<rnow.and.d(i,js,k)<=1d-99)then
     shell = shell + dvol(i,js,k)
    end if
   end do
  end do
  ! Evaluate the mass coordinate
  do j = 1, nol
   if(rnow>=rpol(j).and.rnow<rpol(j+1))then
    mnow = ( mpol(j+1)*(rnow-rpol(j))+mpol(j)*(rpol(j+1)-rnow) )&
         / ( rpol(j+1)-rpol(j) )
    shelld = (mnow-mold)/shell
    exit
   end if
  end do
  ! Set density for chosen cells
  do k = ks, ke
   do i = is, ie
    if(rdis(i,k)<=rnow.and.d(i,js,k)<=1d-99)then
     d(i,js,k) = shelld
     do j = 1, nol
      if(rdis(i,k)>=rpol(j).and.rdis(i,k)<rpol(j+1))then
       p(i,js,k) = ( ppol(j+1)*(rdis(i,k)-rpol(j))+ppol(j)*(rpol(j+1)-rdis(i,k)) )&
                 / (rpol(j+1)-rpol(j))
       exit
      end if
     end do
!!$!     p(i,js,k) = G*mnow*(mnow-mold)/shell/rdis(i,k)
    end if
   end do
  end do
  mold = mnow ; rnow = rnow + dr
  if(rnow>radius)then
   shell = 0d0
   do k = ks, ke
    do i = is, ie
     if(rdis(i,k)<=rnow.and.d(i,js,k)<=1d-99)then
      shell = shell + dvol(i,js,k)
     end if
    end do
   end do
   do k = ks, ke
    do i = is, ie
     if(rdis(i,k)<=rnow.and.d(i,js,k)<=1d-99)then
      d(i,js,k) = (mass-mold)/shell
      p(i,js,k) = G*mass*(mass-mold)/shell/rdis(i,k)
     elseif(rdis(i,k)>rnow)then
      d(i,js,k) = dbg
     end if
    end do
   end do
   exit
  end if
 end do

!!$
!!$ call gravsetup
!!$ call gravity
!!$
!!$ do k = ks, ke
!!$  do i = is, ie
!!$   p(i,js,k) = -grvphi(i,js,k)*d(i,js,k)
!!$  end do
!!$ end do

!!$ do k = ks, ke
!!$  do i = is, ie
!!$   coredis = sqrt(x1(i)*x1(i)+x3(k)*x3(k))
!!$   if(coredis<radius)then
!!$    do j = 1, nol
!!$     if(coredis>=rpol(j).and.coredis<rpol(j+1))then
!!$      d(i,js,k) = ( dpol(j+1)*(coredis-rpol(j))+dpol(j)*(rpol(j+1)-coredis) )&
!!$                / (rpol(j+1)-rpol(j))
!!$      p(i,js,k) = ( ppol(j+1)*(coredis-rpol(j))+ppol(j)*(rpol(j+1)-coredis) )&
!!$                / (rpol(j+1)-rpol(j))
!!$      if(d(i,js,k)<1d-10)d(i,js,k) = 1d-10
!!$      exit
!!$     end if
!!$    end do
!!$   end if
!!$  end do
!!$ end do

 

 do k = ks, ke
  do i = is, ie
   if(rdis(i,k)>=radius)then
    p(i,js,k) = G*mass*d(i,js,k)/rdis(i,k)!-1d5+ G*mass*d(i,js,k)/coredis
   end if
  end do
 end do


! set ejecta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 time = tstart / ejectadistance * (sep-radius*1.01d0)
 t_out = (int(time / dt_out) + 1 )* dt_out

 j = js

 do k = ks,ke
  do i = is,ie

   if(nsdis(i,j,k)<=(sep-radius*1.01d0))then
    do nn = 1, count-1
     if(time*nsdfr(i,j,k)>=t_ej(nn).and.time*nsdfr(i,j,k)<t_ej(nn+1))then
      d(i,j,k)  = (d_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  d_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * nsdfr(i,j,k)**3.d0
      p(i,j,k)  = (p_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  p_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * nsdfr(i,j,k)**(3.d0)
      v1(i,j,k) = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * nssin(i,j,k)! * nsdfr(i,j,k)
      v3(i,j,k) = (v_ej(nn)   * (t_ej(nn+1)-time*nsdfr(i,j,k))  &
                -  v_ej(nn+1) * (t_ej(nn)  -time*nsdfr(i,j,k))) &
                / (t_ej(nn+1)-t_ej(nn)) &
                * (-nscos(i,j,k))! * nsdfr(i,j,k)
      p(i,j,k) = 0.5d0*d(i,j,k)*(pw(2,v1(i,j,k))+pw(2,v3(i,j,k)))*1d-1
      if(p_ej(nn)>pmax*1.d-3)then
       p(i,j,k) = pmax !for shock condition
!       e(i,j,k) = p(i,j,k) / (gamma-1.d0)
      end if
     end if

    end do
   end if

  end do
 end do

! set external gravitational field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(include_extgrv)then
 do k = ks-2, ke+2
  do i = is-2, ie+2
   extgrv(i,j,k) = -1.4d0*msun*G/nsdis(i,j,k)
  end do
 end do
 extgrv(is-1,js:je,ks:ke) = extgrv(is+1,js:je,ks:ke)
 extgrv(is-2,js:je,ks:ke) = extgrv(is,js:je,ks:ke)
end if

 grvtime = time

!call output
 
return
end subroutine headoncollision

