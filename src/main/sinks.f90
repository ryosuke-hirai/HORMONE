module sink_mod

 use settings,only:nsink
 use derived_types,only:sink_prop

 implicit none

 type(sink_prop),allocatable,public:: sink(:)
 real(8),allocatable,public:: snkphi(:,:,:)

 public:: sink_motion,sinkfield,get_sink_acc,sink_accretion
 private:: get_sinksink_acc,get_sink_loc,get_sinkgas_acc

contains
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SINK_MOTION
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
!                         SUBROUTINE SINK_ACCRETION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To apply mass accretion onto a sink particle

subroutine sink_accretion

 use settings,only:eostype,eq_sym,mag_on
 use constants,only:G,pi
 use utils,only:softened_pot,get_vcar,cross,get_vpol
 use grid,only:is,ie,js,je,ks,ke,dt,dvol,car_x,x3
 use physval
 use pressure_mod,only:eos_p,eos_e
 use mpi_utils,only:allreduce_mpi
 use profiler_mod

 integer:: n,i,j,k
 real(8):: dis, philoc, itauacc, newd, newe, newp, accmass, dm, spc1, spc2
 real(8):: r_acc, r_sof, phi0, ambmass, kappa, exchange(8), accmass_in
 real(8),dimension(1:3):: accmom, accang, totmom, vcell, vrel
 real(8):: Ejet, mjet, Ijet, ang_to_ns, perp, vjet, r_circ, omegajet(3), jetrotdir(3), exchange2(5)

!-----------------------------------------------------------------------------

 call start_clock(wtacc)

 kappa = 0.d0 ! Temporary: ignore radiative heating

 do n = 1, nsink
  if(sink(n)%mass<=0d0)cycle
  accmass = 0d0
  accmom  = 0d0
  accang  = 0d0
  ambmass = 0d0
  Ijet = 0d0
  mjet = 0d0
  Ejet = 0d0

  r_acc = max(sink(n)%laccr,sink(n)%locres)
  r_sof = max(sink(n)%lsoft,sink(n)%locres)
  phi0 = G*sink(n)%mass*softened_pot(r_sof,r_sof)
!$omp parallel do private(i,j,k,dis) reduction(+:ambmass) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     dis = norm2(car_x(:,i,j,k)-sink(n)%x)
     if(dis<=r_acc) ambmass = ambmass + d(i,j,k)*dvol(i,j,k)
    end do
   end do
  end do
!$omp end parallel do

  call allreduce_mpi('sum',ambmass)

! Roughly estimate heat rate by optical depth
  sink(n)%facc = 1d0 - exp(-kappa*ambmass/(4d0*pi/3d0*r_acc**2))

!$omp parallel do private(i,j,k,dis,philoc,itauacc,newd,vcell,vrel,&
!$omp dm,ang_to_ns,perp) &
!$omp reduction(+:accmass,accmom,accang,Ejet) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     dis = norm2(car_x(:,i,j,k)-sink(n)%x)

     if(dis>r_acc)cycle

! Compute Cartesian velocity in cell
     call get_vcar(car_x(:,i,j,k),x3(k),v1(i,j,k),v2(i,j,k),v3(i,j,k),vcell)
     vrel = vcell-sink(n)%v

! Suck out mass in the accretion radius on the dynamical timescale
     philoc = G*sink(n)%mass * softened_pot(dis,r_sof)
!     if(0.5d0*dot_product(vrel,vrel)+eint(i,j,k)/d(i,j,k)+philoc>0d0)cycle
     itauacc = sqrt(-philoc+phi0)/r_acc
     newd = d(i,j,k)*exp(-dt*itauacc)

! Add up mass/momentum/AM lost from each cell
     dm = (d(i,j,k) - newd)*dvol(i,j,k)
     accmass = accmass + dm
     accmom = accmom + vcell*dm
     accang = accang + cross(car_x(:,i,j,k)-sink(n)%x,vcell-sink(n)%v)*dm

! Calculate jet quantities
     perp = sqrt(dis**2 - (car_x(3,i,j,k)-sink(n)%x(3))**2)
     ang_to_ns = asin(perp/dis)/pi*180d0

     if(ang_to_ns<=sink(n)%jet_ang)then
      Ejet = Ejet + 0.5d0*newd*dvol(i,j,k)*dot_product(vcell,vcell)
     end if

    end do
   end do
  end do
!$omp end parallel do

  exchange = [accmass,accmom(1),accmom(2),accmom(3),accang(1),accang(2),accang(3),Ejet]
  call allreduce_mpi('sum',exchange)
  accmass = exchange(1)
  accmom = exchange(2:4)
  accang = exchange(5:7)
  Ejet = exchange(8)

  if(eq_sym)then
   accmass = 2d0*accmass
   accmom(1:2) = 2d0*accmom(1:2)
   accmom(3) = 0d0
   accang(1:2) = 0d0
   accang(3) = 2d0*accang(3)
   Ejet = 2d0*Ejet
  end if

! Compute circularization radius of accreting material
  r_circ = dot_product(accang,accang)/accmass**2/(G*sink(n)%mass)
  r_circ = max(r_circ,sink(n)%racc)
! Mass accreted onto NS
  accmass_in = accmass * (sink(n)%racc/r_circ)**p_acc
! Jet energy is determined by accretion energy
  Ejet = Ejet + G*sink(n)%mass*accmass_in/sink(n)%racc

  accang = 0d0
!$omp parallel do private(i,j,k,dis,philoc,itauacc,newd,vcell,vrel,spc1,spc2,&
!$omp dm,ang_to_ns,perp,newp,newe) &
!$omp reduction(+:mjet,Ijet,accang) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     dis = norm2(car_x(:,i,j,k)-sink(n)%x)

     if(dis>r_acc)cycle

! Compute Cartesian velocity in cell
     call get_vcar(car_x(:,i,j,k),x3(k),v1(i,j,k),v2(i,j,k),v3(i,j,k),vcell)
     vrel = vcell-sink(n)%v

! Suck out mass in the accretion radius on the dynamical timescale
     philoc = G*sink(n)%mass * softened_pot(dis,r_sof)
!     if(0.5d0*dot_product(vrel,vrel)+eint(i,j,k)/d(i,j,k)+philoc>0d0)cycle
     itauacc = sqrt(-philoc+phi0)/r_acc
     dm = d(i,j,k)*(1d0-exp(-dt*itauacc))*accmass_in/accmass
     newd = d(i,j,k)-dm
     newp = newd/d(i,j,k)*p(i,j,k)
     select case(eostype)
     case(0:1)
      spc1=0d0;spc2=0d0
     case(2)
      spc1=spc(1,i,j,k); spc2=spc(2,i,j,k)
     end select
     newe = eos_e(newd,newp,T(i,j,k),imu(i,j,k),spc1,spc2) &
          + G*sink(n)%mass*(d(i,j,k)-newd)/sink(n)%racc*sink(n)%facc
     newp = eos_p(newd,newe,T(i,j,k),imu(i,j,k),spc1,spc2)

! Add up mass/momentum/AM lost from each cell
     accang = accang + cross(car_x(:,i,j,k)-sink(n)%x,vcell-sink(n)%v)*dm

! Update primitive variables
     d(i,j,k) = newd
     p(i,j,k) = newp
     eint(i,j,k) = newe
     e(i,j,k) = eint(i,j,k) &
              + 0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2)
     if(mag_on) &
      e(i,j,k) = e(i,j,k) + 0.5d0*(b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2)
! Update conservative variables
     u(i,j,k,icnt) = newd
     u(i,j,k,imo1) = newd*v1(i,j,k)
     u(i,j,k,imo2) = newd*v2(i,j,k)
     u(i,j,k,imo3) = newd*v3(i,j,k)
     u(i,j,k,iene) = e(i,j,k)

! Calculate jet quantities
     perp = sqrt(dis**2 - (car_x(3,i,j,k)-sink(n)%x(3))**2)
     ang_to_ns = asin(perp/dis)/pi*180d0

     if(ang_to_ns<=sink(n)%jet_ang)then
      mjet = mjet + d(i,j,k)*dvol(i,j,k)
      Ijet = Ijet + d(i,j,k)*perp**2*dvol(i,j,k)
     end if

    end do
   end do
  end do
!$omp end parallel do

  exchange2 = [mjet,Ijet,accang(1),accang(2),accang(3)]
  call allreduce_mpi('sum',exchange2)
  mjet = exchange2(1)
  Ijet = exchange2(2)
  accang = exchange2(3:5)

  if(eq_sym)then
   mjet = 2d0*mjet
   Ijet = 2d0*Ijet
   accang(1:2) = 0d0
   accang(3) = 2d0*accang(3)
  end if

! The rest goes into outflow
  omegajet = accang/Ijet
  vjet = sqrt(2d0*(Ejet-Ijet*dot_product(omegajet,omegajet))/mjet)

!$omp parallel do private(i,j,k,dis,vcell,ang_to_ns,perp,jetrotdir) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie

     dis = norm2(car_x(:,i,j,k)-sink(n)%x)

     if(dis>r_acc)cycle

     ! Calculate jet quantities
     perp = sqrt(dis**2 - (car_x(3,i,j,k)-sink(n)%x(3))**2)
     ang_to_ns = asin(perp/dis)/pi*180d0

     if(ang_to_ns>sink(n)%jet_ang)cycle
     jetrotdir = cross(accang,car_x(:,i,j,k)-sink(n)%x)
     jetrotdir = jetrotdir/max(norm2(jetrotdir),epsilon(car_x(:,i,j,k)))
     vcell = vjet*(car_x(:,i,j,k)-sink(n)%x)/dis + jetrotdir*perp*omegajet &
           + sink(n)%v
     call get_vpol(car_x(:,i,j,k),x3(k),vcell,v1(i,j,k),v2(i,j,k),v3(i,j,k))

! Update conservative variables
     u(i,j,k,imo1) = d(i,j,k)*v1(i,j,k)
     u(i,j,k,imo2) = d(i,j,k)*v2(i,j,k)
     u(i,j,k,imo3) = d(i,j,k)*v3(i,j,k)
     e(i,j,k) = eint(i,j,k) &
              + 0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2+v3(i,j,k)**2)
     if(mag_on) &
      e(i,j,k) = e(i,j,k) + 0.5d0*(b1(i,j,k)**2+b2(i,j,k)**2+b3(i,j,k)**2)
     u(i,j,k,iene) = e(i,j,k)

    end do
   end do
  end do
!$omp end parallel do

! Update sink properties in response to accretion
  totmom = sink(n)%mass*sink(n)%v + accmom
  sink(n)%mdot = accmass_in/dt
  sink(n)%mass = sink(n)%mass + accmass_in
  sink(n)%acclum = G*sink(n)%mass*accmass_in/sink(n)%racc/dt
  sink(n)%v = totmom / sink(n)%mass
  sink(n)%jdot = accmass_in*sqrt(G*sink(n)%mass*sink(n)%racc)/dt * sign(1d0,accang)
  sink(n)%Jspin = sink(n)%Jspin + sink(n)%jdot*dt

 end do

 call stop_clock(wtacc)

return
end subroutine sink_accretion


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
 use utils,only:softened_pot
 use gravmod,only:totphi

 integer:: n, i,j,k
 real(8):: dis

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,n,dis) collapse(3)
 do k = ks-1, ke+1
  do j = js-1, je+1
   do i = is-1, ie+1
    snkphi(i,j,k) = 0d0
    do n = 1, nsink
     dis = norm2(car_x(:,i,j,k)-sink(n)%x)
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

 use mpi_utils,only:allreduce_mpi
 use mpi_domain,only:is_my_domain
 use constants,only:tiny
 use settings,only:crdnt,eq_sym,courant
 use grid,only:idx1,idx3,x1,g22,dxi1,dxi2,dxi3,x1,x3,ie,ke
 use gravmod,only:gphi=>grvphi

 type(sink_prop),intent(inout):: sink
 real(8),dimension(1:3):: xpol,acc,acar
 integer:: i,j,k

!-----------------------------------------------------------------------------

 i = sink%i
 j = sink%j
 k = sink%k
 xpol = sink%xpol


 acar(1:3) = 0d0

! Only calculate for MPI thread that contains the sink
 if(is_my_domain(i,j,k))then
  select case(crdnt)
  case(2) ! for spherical coordinates

   if(i==0)then
    ! Do nothing if sink is interior to the innermost cell centre

   elseif(eq_sym)then ! for equatorial symmetry

    ! temporary workaround to avoid corner boundary problems
    if(k==ke.and.i==ie)gphi(i+1,j,k+1)=gphi(i,j,k+1)+gphi(i+1,j,k)-gphi(i,j,k)

! acc is in polar coordinates
    acc(1) = (-(gphi(i+1,j,k  )-gphi(i,j,k  ))*(x3(k+1)-xpol(3)) &
              -(gphi(i+1,j,k+1)-gphi(i,j,k+1))*(xpol(3)-x3(k  )))&
           *idx3(k+1)*idx1(i+1)
    acc(2) = 0d0
    acc(3) = (-(gphi(i  ,j,k+1)-gphi(i  ,j,k))/abs(x1(i  ))*(x1(i+1)-xpol(1)) &
              -(gphi(i+1,j,k+1)-gphi(i+1,j,k))/abs(x1(i+1))*(xpol(1)-x1(i  )))&
           *idx1(i+1)*idx3(k+1)

! sink%a is in cartesian coordinates
    acar(1) = acc(1)*cos(xpol(3)) - acc(3)*sin(xpol(3))
    acar(2) = acc(1)*sin(xpol(3)) + acc(3)*cos(xpol(3))
    acar(3) = 0d0

   else
    print*,'Error in sinks.f90: get_sinkgas_acc'
    stop 'Sink particles currently only implemented for eq_sym=.true.'
   end if
  
  case default
   print*,'Error in sinks.f90: get_sinkgas_acc'
   stop 'Sink particles currently only implemented for spherical coordinates'
  
  end select

 end if

 call allreduce_mpi('sum',acar)

 sink%a = acar
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
 use grid
 use utils,only:carpol

 type(sink_prop),intent(inout):: sink
 integer:: i,j,k

!-----------------------------------------------------------------------------

 select case(crdnt)
 case(0) ! Cartesian coordinates

  sink%xpol = sink%x
  if(ie_global/=is_global)then
   do i = is_global-1, ie_global
    if(sink%x(1)>=x1(i).and.sink%x(1)<x1(i+1))exit
   end do
   sink%i = i
  else
   sink%i = is_global
  end if
  if(je_global/=js_global)then
   do j = js_global-1, je_global
    if(sink%x(2)>=x2(j).and.sink%x(2)<x2(j+1))exit
   end do
   sink%j = j
  else
   sink%j = js_global
  end if
  if(ke_global/=ks_global)then
   do k = ks_global-1, ke_global
    if(sink%x(3)>=x3(k).and.sink%x(3)<x3(k+1))exit
   end do
   sink%k = k
  else
   sink%k = ks_global
  end if

 case(2) ! Spherical coordinates

  sink%xpol = carpol(sink%x)
  do i = is_global-1, ie_global
   if(sink%xpol(1)>=x1(i).and.sink%xpol(1)<x1(i+1))exit
  end do
  sink%i = i

  if(je_global/=js_global)then
   do j = js_global-1, je_global
    if(sink%xpol(2)>=x2(j).and.sink%xpol(2)<x2(j+1))exit
   end do
   sink%j = j
  else
   sink%j = js_global
  end if

  if(ke_global/=ks_global)then
   do k = ks_global-1, ke_global
    if(sink%xpol(3)>=x3(k).and.sink%xpol(3)<x3(k+1))exit
   end do
   sink%k = k
  else
   sink%k = ks_global
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
