module hydro_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE ADVANCE_STEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To select integration scheme for advancing one time step

subroutine advance_step

 use settings,only:radswitch

!-----------------------------------------------------------------------------

 select case(radswitch)
 case(2) ! Use IMEX scheme for Moens et al. 2022 style radiation hydro
  call radhydro_imex_step
 case default ! Use TVD Runge-Kutta as default
  call hydro_step
 end select

 return
end subroutine advance_step

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE HYDRO_STEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To solve hydrodynamic equations for one time step via the
!          TVD (Time Variation Diminishing) Runge-Kutta scheme

subroutine hydro_step

 use settings,only:dirichlet_on,rktype,radswitch,solve_hydro
 use grid,only:rungen
 use boundary_mod,only:boundarycondition
 use numflux_mod,only:numflux
 use source_mod,only:source,phidamp
 use rungekutta_mod,only:rungekutta
 use dirichlet_mod,only:dirichletbound
 use shockfind_mod,only:shockfind
 use mpi_domain,only:exchange_mpi
 use profiler_mod,only:start_clock,stop_clock,wthyd
 use radiation_mod,only:radiative_diffusion
use physval
!-----------------------------------------------------------------------------

 if(solve_hydro)then

  call start_clock(wthyd)

  if(dirichlet_on) call dirichletbound
  call shockfind

  do rungen = 1, rktype
   call exchange_mpi
   call boundarycondition
   call numflux
   call source
   call rungekutta
  end do

  call phidamp

  call stop_clock(wthyd)

 end if

 if(radswitch==1)call radiative_diffusion

 return
end subroutine hydro_step

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE RADHYDRO_IMEX_STEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To solve radiation hydrodynamic equations for one time step via the
!          IMEX scheme (Moens et al. 2022)
!          Note that radiation is integrated into the hydro updates here

subroutine radhydro_imex_step

 use settings,only:dirichlet_on,spn,compswitch
 use grid,only:dt,is,ie,js,je,ks,ke
 use physval,only:u,uorg,d,spc,spcorg,irad,ufnmax,T
 use boundary_mod,only:boundarycondition
 use numflux_mod,only:numflux
 use source_mod,only:source
 use rungekutta_mod,only:primitive
 use dirichlet_mod,only:dirichletbound
 use shockfind_mod,only:shockfind
 use mpi_domain,only:exchange_mpi
 use radiation_mod,only:radiative_diffusion

 integer:: ufn,i,j,k,n
 real(8):: dt_global,G_ex_half
 real(8),allocatable:: u_half(:,:,:,:),G_im_half(:,:,:),spc_half(:,:,:,:),u_half_plus(:,:,:,:)

!-----------------------------------------------------------------------------

 allocate(u_half,u_half_plus,mold=u)
 allocate(G_im_half,mold=d)
 if(compswitch>=2)allocate(spc_half,mold=spc)

 dt_global = dt

 if(dirichlet_on) call dirichletbound
 call shockfind

!$omp parallel
!$omp do private (ufn,i,j,k) collapse(4)
 do ufn = 1,ufnmax

  do k = ks,ke
   do j = js,je
    do i = is,ie
     uorg(i,j,k,ufn) = u(i,j,k,ufn)
    end do
   end do
  end do

 end do
!$omp end do
 if(compswitch>=2)then
!$omp do private (i,j,k,n) collapse(4)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     do n = 1, spn
      spcorg(n,i,j,k) = spc(n,i,j,k)*d(i,j,k)
      spc(n,i,j,k) = spcorg(n,i,j,k)
     end do
    end do
   end do
  end do
!$omp end do
 end if
!$omp end parallel

 ! First do half an explicit step
 dt = dt_global
 call imex_explicit_terms

!$omp parallel
!$omp do private (ufn,i,j,k) collapse(4)
 do ufn = 1,ufnmax
  do k = ks,ke
   do j = js,je
    do i = is,ie
     !u_half_plus(i,j,k,ufn) = u(i,j,k,ufn)
     u(i,j,k,ufn) = 0.5d0*(u(i,j,k,ufn)+uorg(i,j,k,ufn))
     u_half_plus(i,j,k,ufn) = u(i,j,k,ufn)
    end do
   end do
  end do
 end do
!$omp end do
 if(compswitch>=2)then
!$omp do private(i,j,k,n) collapse(4)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     do n = 1, spn
      spc_half(n,i,j,k) = spc(n,i,j,k)
     end do
    end do
   end do
  end do
!$omp end do
 end if
!$omp end parallel

 ! Then do implicit diffusion for half a step
 call exchange_mpi
 dt = dt_global
 call radiative_diffusion

 ! Compute the implicit derivative
!$omp parallel do private (ufn,i,j,k) collapse(4)
 do ufn = 1,ufnmax

  do k = ks,ke
   do j = js,je
    do i = is,ie
     if(ufn==irad)G_im_half(i,j,k) = u(i,j,k,ufn)-u_half_plus(i,j,k,ufn)
     u_half(i,j,k,ufn) = 0.5d0*(u(i,j,k,ufn)+u_half_plus(i,j,k,ufn))
     u(i,j,k,ufn) = u_half(i,j,k,ufn)
    end do
   end do
  end do

 end do
!$omp end parallel do

 ! Finally do an explicit step from n+1/2
 dt = dt_global
 call imex_explicit_terms

 ! Compute the explicit derivative
!$omp parallel
!$omp do private (ufn,i,j,k,G_ex_half) collapse(4)
 do ufn = 1,ufnmax
  do k = ks,ke
   do j = js,je
    do i = is,ie
     G_ex_half = u(i,j,k,ufn) - u_half(i,j,k,ufn)
     u(i,j,k,ufn) = uorg(i,j,k,ufn) + G_ex_half
     if(ufn==irad)u(i,j,k,ufn) = u(i,j,k,ufn) + G_im_half(i,j,k)
    end do
   end do
  end do
 end do
!$omp end do
 if(compswitch>=2)then
!$omp do private(n,i,j,k) collapse(3)
  do k = ks, ke
   do j = js, je
    do i = is, ie
     do n = 1, spn
      print*,n,i,j,k
      spc(n,i,j,k) = spcorg(n,i,j,k) + spc(n,i,j,k)-spc_half(n,i,j,k)
     end do
     spc(:,i,j,k) = spc(:,i,j,k) / sum(spc(:,i,j,k))
    end do
   end do
  end do
!$omp end do
 end if
!$omp end parallel

 deallocate(u_half,u_half_plus,G_im_half)
 if(compswitch>=2)deallocate(spc_half)

 call primitive

 return
end subroutine radhydro_imex_step

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE IMEX_EXPLICIT_TERMS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To update explicit terms in the IMEX scheme
!          Corresponds to G_ex in Moens et al. 2022

subroutine imex_explicit_terms

 use settings,only:compswitch,spn,solve_hydro
 use grid,only:dt,is,ie,js,je,ks,ke
 use physval,only:u,spc,src,ufnmax,irad,iene
 use source_mod,only:source,phidamp
 use radiation_mod,only:rad_heat_cool
 use rungekutta_mod,only:primitive,flux_sum,spcflx_sum
 use numflux_mod,only:numflux
 use eos_mod,only:pressure
 use mpi_domain,only:exchange_mpi
 use boundary_mod,only:boundarycondition
 use profiler_mod

 integer:: ufn,i,j,k,n

!-----------------------------------------------------------------------------

 if(solve_hydro)then
  call start_clock(wthyd)

  call boundarycondition
  call source

!$omp parallel do private (ufn,i,j,k) collapse(4)
  do ufn = 1,ufnmax
   do k = ks,ke
    do j = js,je
     do i = is,ie
      u(i,j,k,ufn) = u(i,j,k,ufn) + dt * src(i,j,k,ufn)
     end do
    end do
   end do
  end do
!$omp end parallel do

  call primitive

  call stop_clock(wthyd)
 end if

 call rad_heat_cool

 if(solve_hydro)then
  call start_clock(wthyd)
  call exchange_mpi
  call boundarycondition
  call numflux

!$omp parallel
!$omp do private (ufn,i,j,k) collapse(4)
  do ufn = 1,ufnmax
   do k = ks,ke
    do j = js,je
     do i = is,ie
      u(i,j,k,ufn) = u(i,j,k,ufn) + dt * flux_sum(i,j,k,ufn)
     end do
    end do
   end do
  end do
!$omp end do
  if(compswitch>=2)then
!$omp do private (i,j,k,n) collapse(4)
   do k = ks, ke
    do j = js, je
     do i = is, ie
      do n = 1, spn
       spc(n,i,j,k) = spc(n,i,j,k) + dt * spcflx_sum(n,i,j,k)
      end do
     end do
    end do
   end do
!$omp end do
  end if
!$omp end parallel

  call phidamp

  call primitive

  call stop_clock(wthyd)
 end if

 return
end subroutine imex_explicit_terms

end module hydro_mod
