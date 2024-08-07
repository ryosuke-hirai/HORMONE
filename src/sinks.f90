module sink_mod
 implicit none

 type sink_prop
  sequence
  integer:: i,j,k
  real(8):: mass, softfac, lsoft, locres, dt
  real(8),dimension(1:3):: x,v,a,xpol
 end type sink_prop
 integer,public:: nsink
 type(sink_prop),allocatable,public:: sink(:)
 real(8),allocatable,public:: snkphi(:,:,:)

 public:: sink_motion,sinkfield,get_sink_acc,get_sink_loc,get_sinkgas_acc,&
          get_sinksink_acc

contains
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE SINK_MOTION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To update position and velocity of sink particles

subroutine sink_motion

 use grid,only:dt,frame_acc
 use profiler_mod

 integer:: n

!-----------------------------------------------------------------------------

 call start_clock(wtsnk)

 do n = 1, nsink
  sink(n)%a = sink(n)%a + frame_acc
  sink(n)%v = sink(n)%v + sink(n)%a*dt
  sink(n)%x = sink(n)%x + sink(n)%v*dt
 end do

 call stop_clock(wtsnk)

 return
end subroutine sink_motion

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SINKFIELD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the field from sinks

subroutine sinkfield

 use constants,only:G
 use grid
 use physval
 use utils,only:polcar,softened_pot
 use gravmod,only:totphi

 integer:: n, i,j,k
 real(8):: dis, xcar(1:3)

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,n,dis,xcar) collapse(3)
 do k = gks-1, gke+1
  do j = gjs-1, gje+1
   do i = gis-1, gie+1
    xcar = polcar([x1(i),x2(j),x3(k)])
    snkphi(i,j,k) = 0d0
    do n = 1, nsink
     dis = norm2(xcar-sink(n)%x)
     snkphi(i,j,k) = snkphi(i,j,k) &
                   + G*sink(n)%mass &
                     *softened_pot(dis,max(sink(n)%lsoft,sink(n)%locres))
    end do
    totphi(i,j,k) = totphi(i,j,k) + snkphi(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine sinkfield

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE GET_SINK_ACC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate teh total acceleration on sinks

subroutine get_sink_acc(sink)

 use constants,only:huge
 use grid,only:dt
 use profiler_mod

 type(sink_prop),allocatable,intent(inout):: sink(:)
 integer:: n
 real(8):: dtsink

!-----------------------------------------------------------------------------

 call start_clock(wtsnk)

 dtsink = huge
 do n = 1, nsink
  call get_sink_loc(sink(n))
  call get_sinkgas_acc(sink(n))
  dtsink = min(dtsink,sink(n)%dt)
 end do

 call get_sinksink_acc(sink)

 dt = min(dtsink,dt) ! update dt

 call stop_clock(wtsnk)

 return
end subroutine get_sink_acc

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GET_SINKGAS_ACC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate the acceleration from sink-gas interactions

subroutine get_sinkgas_acc(sink)

 use constants,only:tiny
 use settings,only:crdnt,eq_sym,courant
 use utils,only:carpol
 use grid,only:idx1,idx3,x1,g22,dxi1,dxi2,dxi3,x1,x3
 use gravmod,only:gphi=>grvphi

 type(sink_prop),intent(inout):: sink
 real(8):: xpol(1:3),acc(1:3)
 integer:: i,j,k

!-----------------------------------------------------------------------------

 i = sink%i
 j = sink%j
 k = sink%k
 xpol = sink%xpol

 select case(crdnt)
 case(2) ! for spherical coordinates

  if(eq_sym)then ! for equatorial symmetry

   acc(1) = (-(gphi(i+1,j,k  )-gphi(i,j,k  ))*idx1(i+1)*(x3(k+1)-xpol(3))&
             -(gphi(i+1,j,k+1)-gphi(i,j,k+1))*idx1(i+1)*(xpol(3)-x3(k)))&
          *idx3(k+1)
   acc(2) = 0d0
   acc(3) = (-(gphi(i  ,j,k+1)-gphi(i  ,j,k))*idx3(k)/x1(i)*(x1(i+1)-xpol(1))&
             -(gphi(i+1,j,k+1)-gphi(i+1,j,k))*idx3(k)/x1(i+1)*(xpol(1)-x1(i)))&
          *idx1(i+1)

   sink%a(1) = acc(1)*cos(xpol(3)) - acc(3)*sin(xpol(3))
   sink%a(2) = acc(1)*sin(xpol(3)) + acc(3)*cos(xpol(3))
   sink%a(3) = 0d0

  else
   print*,'Error in sinks.f90: get_sinkgas_acc'
   stop 'Sink particles currently only implemented for eq_sym=.true.'
  end if

 case default
  print*,'Error in sinks.f90: get_sinkgas_acc'
  stop 'Sink particles currently only implemented for spherical coordinates'

 end select

 sink%dt = min(dxi1(i),g22(i)*dxi2(j),g22(i)*dxi3(k)) &
         / max(norm2(sink%v),tiny)*courant

return
end subroutine get_sinkgas_acc

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GET_SINKSINK_ACC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get acceleration from sink-sink interactions

subroutine get_sinksink_acc(sink)

 use constants,only:G,tiny
 use utils,only:softened_acc

 type(sink_prop),allocatable,intent(inout):: sink(:)
 integer:: n,m
 real(8):: dis, softacc, dir(1:3)

!-----------------------------------------------------------------------------

 if(nsink==1)return ! No sink-sink interactions if there is only one sink

 do n = 2, nsink
  do m = 1, n-1
   dir = sink(n)%x-sink(m)%x
   dis = max(norm2(dir),tiny)
   dir = dir/dis
   softacc = softened_acc(dis,max(sink(n)%lsoft,sink(m)%lsoft))
   sink(n)%a = sink(n)%a - G*sink(m)%mass*softacc*dir
   sink(m)%a = sink(m)%a + G*sink(n)%mass*softacc*dir
  end do
 end do

! UPDATE TIME STEP FOR SINK-SINK INTERACTION

return
end subroutine get_sinksink_acc


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE GET_SINK_LOC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get sink location id

subroutine get_sink_loc(sink)

 use settings,only:crdnt
 use grid,only:is,ie,js,je,ks,ke,x1,x2,x3,dxi1,dxi2,dxi3,g22
 use utils,only:carpol

 type(sink_prop),intent(inout):: sink
 integer:: i,j,k

!-----------------------------------------------------------------------------

 select case(crdnt)
 case(0) ! Cartesian coordinates

  sink%xpol = sink%x
  if(ie/=is)then
   do i = is-1, ie
    if(sink%x(1)>=x1(i).and.sink%x(1)<x1(i+1))exit
   end do
   sink%i = i
  end if
  if(je/=js)then
   do j = js-1, je
    if(sink%x(2)>=x2(j).and.sink%x(2)<x2(j+1))exit
   end do
   sink%j = j
  end if
  if(ke/=ks)then
   do k = ks-1, ke
    if(sink%x(3)>=x3(k).and.sink%x(3)<x3(k+1))exit
   end do
   sink%k = k
  end if

 case(2) ! Spherical coordinates

  sink%xpol = carpol(sink%x)
  do i = is-1, ie
   if(sink%xpol(1)>=x1(i).and.sink%xpol(1)<x1(i+1))exit
  end do
  sink%i = i

  if(je/=js)then
   do j = js-1, je
    if(sink%xpol(2)>=x2(j).and.sink%xpol(2)<x2(j+1))exit
   end do
   sink%j = j
  else
   sink%j = js
  end if

  if(ke/=ks)then
   do k = ks-1, ke
    if(sink%xpol(3)>=x3(k).and.sink%xpol(3)<x3(k+1))exit
   end do
   sink%k = k
  else
   sink%k = ks
  end if

 case default
  stop 'Error in get_sink_loc: wrong crdnt'
 end select

 sink%locres = max(dxi1(sink%i),&
                   g22(sink%i)*dxi2(sink%j),&
                   g22(sink%i)*dxi3(sink%k) )&
             * sink%softfac


return
end subroutine get_sink_loc

!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$!                         SUBROUTINE GET_TOTMOM
!!$!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$
!!$! PURPOSE: To calculate total momentum
!!$
!!$ subroutine get_totmom(totmom)
!!$
!!$  use grid
!!$  use physval,only:d,v1,v2,v3
!!$  use gravmod
!!$  use constants,only:G
!!$  use utils
!!$
!!$  implicit none
!!$
!!$  real*8,intent(out)::totmom(1:3)
!!$  real*8,dimension(1:3):: xcar,vcar
!!$  real*8:: nspot
!!$
!!$!-----------------------------------------------------------------------------
!!$
!!$!$omp parallel do private(i,j,k,xcar,vcar,nspot) &
!!$!$omp reduction (+:totmom)
!!$  do k = ks, ke
!!$   do j = js, je
!!$    do i = is, ie
!!$!     call polcar(x1(i),x2(j),x3(k),xcar)
!!$     xcar = polcar((/x1(i),x2(j),x3(k)/))
!!$     call get_vcar(xcar,x3(k),v1(i,j,k),v2(i,j,k),v3(i,j,k),vcar)
!!$     nspot = G*nsmass*softened_pot(norm2(xcar-nspos),nssoft)
!!$
!!$     totmom = totmom + d(i,j,k)*dvol(i,j,k)*vcar
!!$    end do
!!$   end do
!!$  end do
!!$!$omp end parallel do
!!$
!!$  totmom = totmom*2d0
!!$  totmom(3) = 0d0
!!$
!!$  return
!!$ end subroutine get_totmom

end module sink_mod
