module radiation_mod

 use matrix_vars,only:cg_set
 use opacity_mod
 use radiation_utils

 implicit none

 public :: radiation,radiation_setup,radiative_force,rad_heat_cool,get_gradE
 private:: get_diffusion_coeff,get_radb
 real(8),allocatable,private:: rsrc(:),gradE(:,:,:,:)

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE RADIATION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute radiation transport with flux limited diffusion

subroutine radiation

 use settings,only:radswitch
 use constants,only:Cv,Rgas
 use grid
 use physval
 use miccg_mod,only:miccg
 use profiler_mod
 use pressure_mod,only:Trad,get_etot_from_eint
 use matrix_solver_mod,only:write_A_rad,solve_system_rad
 use matrix_utils,only:ijk_from_l,l_from_ijk
 use mpi_domain,only:exchange_scalar

 integer:: l,i,j,k,ll
 integer:: in,jn,kn
 integer:: in_global,jn_global,kn_global
 real(8),allocatable:: x(:)
 integer,allocatable:: map(:)

!-----------------------------------------------------------------------------

 call start_clock(wtrad)

 in = ie-is+1; jn = je-js+1; kn = ke-ks+1
 in_global = ie_global-is_global+1
 jn_global = je_global-js_global+1
 kn_global = ke_global-ks_global+1

 allocate(map(1:in*jn*kn))
 ll = 0
 do k = ks, ke
   do j = js, je
     do i = is, ie
       ll = ll + 1
       l = l_from_ijk(i, j, k, is_global, js_global,ks_global, in_global, jn_global)
       map(ll) = l
     end do
   end do
 end do

! Advection and radiative acceleration terms are updated in hydro step

! Update heating/cooling term first if following Moens+2022
 if(radswitch==2)call rad_heat_cool

! Then update the diffusion term
 call exchange_scalar(erad)
 call rad_boundary
 call get_gradE
 call get_diffusion_coeff ! use erad^n for diffusion coefficients

 call write_A_rad
 call get_radb(map) ! Sets up rsrc

 allocate( x(1:in*jn*kn) )
!$omp parallel do private(l,i,j,k,ll)
 do ll = 1, size(x)
  l = map(ll)
  ! Get the i,j,k indices from the local index
  call ijk_from_l(l,is_global,js_global,ks_global,in_global,jn_global,i,j,k)
  x(ll) = erad(i,j,k)
 end do
!$omp end parallel do

 call solve_system_rad(rsrc, x) ! returns erad^{n+1}

! update erad and u
!$omp parallel do private(l,i,j,k)
 do ll = 1, size(x)
  l = map(ll)
  call ijk_from_l(l,is_global,js_global,ks_global,in_global,jn_global,i,j,k)
  erad(i,j,k) = x(ll)
  T   (i,j,k) = update_Tgas(d(i,j,k),erad(i,j,k),T(i,j,k),dt)
  eint(i,j,k) = Cv  *d(i,j,k)*T(i,j,k)*imu(i,j,k)
  p   (i,j,k) = Rgas*d(i,j,k)*T(i,j,k)*imu(i,j,k)
  e   (i,j,k) = get_etot_from_eint(i,j,k)
  u(i,j,k,iene) = e   (i,j,k)
  u(i,j,k,irad) = erad(i,j,k)
 end do
!$omp end parallel do

 call stop_clock(wtrad)

return
end subroutine radiation

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE GET_GRADE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gradE on each cell

subroutine get_gradE

 use utils,only:get_grad
 use physval,only:erad
 use grid,only:is,ie,js,je,ks,ke

 integer:: i,j,k

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k) collapse(3)
 do k = ks-1, ke+1
  do j = js-1, je+1
   do i = is-1, ie+1
    call get_grad(erad,i,j,k,gradE(1:3,i,j,k))
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_gradE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GET_DIFFUSION_COEFF
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute diffusion coefficients on each cell

subroutine get_diffusion_coeff

 use constants,only:clight
 use grid,only:is,ie,js,je,ks,ke
 use physval,only:erad,d,T,erad,radK

 integer:: i,j,k
 real(8):: RR,ll,kappar

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,RR,ll,kappar) collapse(3)
 do k = ks-1, ke+1
  do j = js-1, je+1
   do i = is-1, ie+1
    ! Skip if this is a corner ghost cell, which is uninitialised and unused)
    if ((i == is-1 .or. i == ie+1) .and. &
      (j == js-1 .or. j == je+1)) cycle
    if ((i == is-1 .or. i == ie+1) .and. &
      (k == ks-1 .or. k == ke+1)) cycle
    if ((j == js-1 .or. j == je+1) .and. &
      (k == ks-1 .or. k == ke+1)) cycle

    kappar = kappa_r(d(i,j,k),T(i,j,k))
    RR = norm2(gradE(1:3,i,j,k)) / (d(i,j,k)*kappar*erad(i,j,k))
    ll = lambda(RR)
    radK(i,j,k) = clight*ll/(d(i,j,k)*kappar)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine get_diffusion_coeff


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE GET_RADB
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To get source term (b) for the radiation diffusion equation

subroutine get_radb(map)

 use settings,only:radswitch
 use constants,only:clight,arad
 use grid,only:dvol,dt,is_global,js_global,ks_global,ie_global,je_global,ke_global
 use physval
 use matrix_utils,only:ijk_from_l

 integer,intent(in):: map(:)

 integer:: i,j,k,l,ll
 real(8):: kappap
 integer :: in_global,jn_global,kn_global
!-----------------------------------------------------------------------------

 in_global = ie_global-is_global+1
 jn_global = je_global-js_global+1
 kn_global = ke_global-ks_global+1

!$omp parallel do private(l,i,j,k,kappap)
 do ll = 1, size(rsrc)
   l = map(ll)
   call ijk_from_l(l,is_global,js_global,ks_global,in_global,jn_global,i,j,k)
! Note: Dirichlet boundary conditions should be included here.
   kappap = kappa_p(d(i,j,k),T(i,j,k))
   rsrc(ll) = erad(i,j,k)*dvol(i,j,k)/dt
   if(radswitch==1) &
    rsrc(ll) = rsrc(ll) + clight*kappap*d(i,j,k)*dvol(i,j,k)*arad &
     * (4d0*T(i,j,k)**3*update_Tgas(d(i,j,k),0d0,T(i,j,k),dt)-3d0*T(i,j,k)**4)
  end do
!$omp end parallel do

return
end subroutine get_radb


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE RADIATION_SETUP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set up radiation parameters

subroutine radiation_setup

 use settings,only:radswitch
 use grid,only:is,ie,js,je,ks,ke
 use matrix_vars,only:irad
 use matrix_solver_mod,only:setup_matrix
 use miccg_mod,only:setup_cg
 use physval,only:radK

!-----------------------------------------------------------------------------

 if(radswitch==1.or.radswitch==2)then
  call setup_matrix(irad)
  call get_geo

  allocate(radK(is-1:ie+1,js-1:je+1,ks-1:ke+1),gradE(1:3,is-1:ie+1,js-1:je+1,ks-1:ke+1),&
           rsrc(1:(ie-is+1)*(je-js+1)*(ke-ks+1)) )
 end if

return
end subroutine radiation_setup

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE RADIATIVE_FORCE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute radiative acceleration and associated source terms

subroutine radiative_force

 use utils,only:get_grad
 use grid
 use physval

 integer:: i,j,k,l,m
 real(8):: RR,ll,ff,vdotfrad,Pedd(1:3,1:3),radwork,kappar
 real(8),dimension(1:3):: frad,gradv1,gradv2,gradv3,nn

!-----------------------------------------------------------------------------

 call get_gradE

!$omp parallel do private(i,j,k,RR,ll,ff,frad,vdotfrad,gradv1,gradv2,gradv3,&
!$omp nn,l,m,radwork,Pedd,kappar) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    kappar = kappa_r(d(i,j,k),T(i,j,k))
    if(kappar<=tiny(kappar))cycle
    RR = norm2(gradE(1:3,i,j,k)) / (d(i,j,k)*kappar*erad(i,j,k))
    ll = lambda(RR)
    ff = ll + (ll*RR)**2

    frad(1:3) = -ll*gradE(1:3,i,j,k)
    vdotfrad = frad(1)*v1(i,j,k) + frad(2)*v2(i,j,k) + frad(3)*v3(i,j,k)

    call get_grad(v1,i,j,k,gradv1)
    call get_grad(v2,i,j,k,gradv2)
    call get_grad(v3,i,j,k,gradv3)

    nn(1:3) = gradE(1:3,i,j,k)/max(norm2(gradE(1:3,i,j,k)),epsilon(erad(i,j,k)))
    do m = 1, 3
     do l = 1, 3
      if(l>m)cycle
      Pedd(l,m) = (3d0*ff-1d0)*nn(l)*nn(m)
      if(l==m)Pedd(l,m) = Pedd(l,m) + (1d0-ff)
      Pedd(l,m) = 0.5d0*erad(i,j,k)*Pedd(l,m)
     end do
    end do

    radwork = Pedd(1,1)*gradv1(1) + Pedd(1,2)*gradv1(2) + Pedd(1,3)*gradv1(3) &
            + Pedd(1,2)*gradv2(1) + Pedd(2,2)*gradv2(2) + Pedd(2,3)*gradv2(3) &
            + Pedd(1,3)*gradv3(1) + Pedd(2,3)*gradv3(2) + Pedd(3,3)*gradv3(3)

    src(i,j,k,imo1) = src(i,j,k,imo1) + frad(1)
    src(i,j,k,imo2) = src(i,j,k,imo2) + frad(2)
    src(i,j,k,imo3) = src(i,j,k,imo3) + frad(3)
    src(i,j,k,iene) = src(i,j,k,iene) + vdotfrad
    src(i,j,k,irad) = -radwork

   end do
  end do
 end do
!$omp end parallel do

return
end subroutine radiative_force

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE RAD_HEAT_COOL
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute the radiative heating and cooling
!          Follows the method of Turner & Stone 2001

subroutine rad_heat_cool

 use constants,only:clight,sigma,Cv
 use grid
 use physval
 use pressure_mod,only:get_etot_from_eint,getT_from_de,Trad

 integer:: i,j,k
 real(8):: a1,a2,c1,c2,kappap,eint1

!-----------------------------------------------------------------------------

!$omp parallel do private(i,j,k,a1,a2,c1,c2,kappap,eint1) collapse(3)
 do k = ks, ke
  do j = js, je
   do i = is, ie
    kappap = kappa_p(d(i,j,k),T(i,j,k))
    if(kappap<=tiny(kappap))cycle

    a1 = 4d0*kappap*sigma/(Cv*imu(i,j,k))**4/d(i,j,k)**3*dt
    a2 = clight*kappap*d(i,j,k)*dt
    c1 = (1d0+a2)/a1
    c2 = -c1*eint(i,j,k)-a2/a1*erad(i,j,k)

    eint1 = max(eint(i,j,k),erad(i,j,k))
    call solve_quartic(c1,c2,eint1)
    erad(i,j,k) = (a1*eint1**4+erad(i,j,k))/(1d0+a2)
    eint(i,j,k) = eint1
    e(i,j,k) = get_etot_from_eint(i,j,k)
    call getT_from_de(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))

    u(i,j,k,iene) = e   (i,j,k)
    u(i,j,k,irad) = erad(i,j,k)
   end do
  end do
 end do
!$omp end parallel do

return
end subroutine rad_heat_cool

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE RAD_BOUNDARY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set boundary conditions for erad

subroutine rad_boundary

 use grid
 use physval

 integer:: i,j,k

!-----------------------------------------------------------------------------

! Only zero-flux boundary for now

!$omp parallel do private(j,k) collapse(2)
 do k = ks, ke
  do j = js, je
   if (is == is_global) erad(is-2,j,k) = erad(is+1,j,k)
   if (is == is_global) erad(is-1,j,k) = erad(is  ,j,k)
   if (ie == ie_global) erad(ie+1,j,k) = erad(ie  ,j,k)
   if (ie == ie_global) erad(ie+2,j,k) = erad(ie-1,j,k)
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,k) collapse(2)
 do k = ks, ke
  do i = is, ie
   if (js == js_global) erad(i,js-2,k) = erad(i,js+1,k)
   if (js == js_global) erad(i,js-1,k) = erad(i,js  ,k)
   if (je == je_global) erad(i,je+1,k) = erad(i,je  ,k)
   if (je == je_global) erad(i,je+2,k) = erad(i,je-1,k)
  end do
 end do
!$omp end parallel do
!$omp parallel do private(i,j) collapse(2)
 do j = js, je
  do i = is, ie
   if (ks == ks_global) erad(i,j,ks-2) = erad(i,j,ks+1)
   if (ks == ks_global) erad(i,j,ks-1) = erad(i,j,ks  )
   if (ke == ke_global) erad(i,j,ke+1) = erad(i,j,ke  )
   if (ke == ke_global) erad(i,j,ke+2) = erad(i,j,ke-1)
  end do
 end do
!$omp end parallel do

return
end subroutine rad_boundary

end module radiation_mod
