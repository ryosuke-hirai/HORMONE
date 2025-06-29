module gravity_elliptic_mod

 implicit none
 private

 public :: gravity_elliptic

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE GRAVITY_ELLIPTIC
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To compute gravitational forces using the sparse matrix method.

subroutine gravity_elliptic

 use grid
 use constants
 use physval
 use gravmod
 use gravbound_mod
 use miccg_mod,only:miccg
 use timestep_mod,only:timestep
 use profiler_mod
 use matrix_utils,only:l_from_ijk,ijk_from_l,contiguous_map
 use matrix_solver_mod,only:solve_system_grv
 use matrix_vars,only:lmax_grv
 use mpi_domain,only:exchange_gravity_mpi
 integer:: i,j,k,l,ll
 real(8),allocatable,dimension(:):: x, cgsrc
 integer,allocatable:: map(:)

!-------------------------------------------------------------------------

 call start_clock(wtelg)

 allocate( x(lmax_grv), cgsrc(lmax_grv) )
 call contiguous_map(gis, gie, gjs, gje, gks, gke, &
                     gis_global, gie_global, gjs_global, gje_global, gks_global, map)

! MICCG method to solve Poisson equation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 if(gbtype==0)call gravbound

! cylindrical (equatorial+axial symmetry) ####################################
 if(gje_global==gjs_global.and.crdnt==1.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l,ll)
  do ll=1,lmax_grv
   l = map(ll)
   call ijk_from_l(l,gis_global,gjs_global,gks_global,gie_global-gis_global+1,gje_global-gjs_global+1,i,j,k)
   x(ll) = grvphi(i,j,k)
   cgsrc(ll) = 0d0
   if(i>=is)then;if(i<=ie)then;if(k>=ks)then;if(k<=ke)then
    cgsrc(ll) = 4d0*pi*G*gsrc(i,j,k)*x1(i)*dxi1(i)*sum(dx3(k:k+1))*0.5d0
   end if;end if;end if;end if
   if(gbtype==0)then
    if(k==gks_global) cgsrc(ll) = cgsrc(ll) - x1 (i)*dxi1(i)*idx3(k  )*phi3i(i,k-1)
    if(i==gie_global) cgsrc(ll) = cgsrc(ll) - xi1(i)*dxi3(k)*idx1(i+1)*phi1o(i+1,k)
    if(k==gke_global) cgsrc(ll) = cgsrc(ll) - x1 (i)*dxi1(i)*idx3(k+1)*phi3o(i,k+1)
   end if
  end do
!$omp end parallel do

! spherical (axial symmetry) #################################################
 elseif(crdnt==2.and.dim==2)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l,ll)
  do ll=1,lmax_grv
   l = map(ll)
   call ijk_from_l(l,gis_global,gjs_global,gks_global,gie_global-gis_global+1,gje_global-gjs_global+1,i,j,k)
   x(ll) = grvphi(i,j,k)
   cgsrc(ll) = 0d0
   if(i>=is)then;if(i<=ie)then;if(j>=js)then;if(j<=je)then
    cgsrc(ll) = 4d0*pi*G*gsrc(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)
   end if;end if;end if;end if
   if(gbtype==0)then
    if(i==gie_global) cgsrc(ll) = cgsrc(ll) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)&
                                    *phiio(i+1,j)
    if(i==gis_global) cgsrc(ll) = cgsrc(ll) - xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)&
                                    *phiii(i-1,j)
   end if
  end do
!$omp end parallel do

! spherical (3D) #############################################################
 elseif(crdnt==2.and.dim==3)then

! calculating b for Ax=b
!$omp parallel do private(i,j,k,l,ll)
  do ll=1,lmax_grv
   l = map(ll)
   call ijk_from_l(l,gis_global,gjs_global,gks_global,gie_global-gis_global+1,gje_global-gjs_global+1,i,j,k)
   x(ll) = grvphi(i,j,k)
   cgsrc(ll) = 0d0
   if(i>=is)then;if(i<=ie)then
    cgsrc(ll) = 4d0*pi*G*gsrc(i,j,k)*x1(i)**2*sinc(j)*dxi1(i)*dxi2(j)*dxi3(k)
   end if;end if

   if(gbtype==0)then
    if(i==gie_global)cgsrc(ll) = cgsrc(ll) - xi1(i)**2*sinc(j)*dxi2(j)*idx1(i+1)*dxi3(k)&
                                  *phiio(i+1,j)
    if(i==gis_global)cgsrc(ll) = cgsrc(ll) - xi1(i-1)**2*sinc(j)*dxi2(j)*idx1(i)*dxi3(k)&
                                  *phiii(i-1,j)
   end if
  end do
!$omp end parallel do
 end if

! ############################################################################

! Solve Poisson equation with the conjugate gradient method
  call solve_system_grv(cgsrc, x)

!-------------------------------------------------------------------------

! convert x to phi
!$omp parallel do private(i,j,k,l,ll)
 do ll = 1, lmax_grv
  l = map(ll)
  call ijk_from_l(l,gis_global,gjs_global,gks_global,gie_global-gis_global+1,gje_global-gjs_global+1,i,j,k)
  grvphi(i,j,k) = x(ll)
 end do
!$omp end parallel do

 call exchange_gravity_mpi


 if(gje_global==gjs_global.and.crdnt==1.and.dim==2)then ! for cylindrical coordinates

  if (gis==gis_global) then
    grvphi(gis-2,gjs:gje,gks:gke) = grvphi(gis+1,gjs:gje,gks:gke)
    grvphi(gis-1,gjs:gje,gks:gke) = grvphi(gis  ,gjs:gje,gks:gke)
  end if
  if (gks==gks_global) then
    grvphi(gis:gie,gjs:gje,gks-2) = grvphi(gis:gie,gjs:gje,gks)
    grvphi(gis:gie,gjs:gje,gks-1) = grvphi(gis:gie,gjs:gje,gks)
  end if

 elseif(gke_global==gks_global.and.crdnt==2.and.dim==2)then ! for spherical coordinates (2D)
  if (gks==gks_global) then
    if (gjs==gjs_global) then
      grvphi(gis-1:gie+1,gjs-2,gks) = grvphi(gis-1:gie+1,gjs+1,gks)
      grvphi(gis-1:gie+1,gjs-1,gks) = grvphi(gis-1:gie+1,gjs,gks)
    endif
    if (gje==gje_global) then
      grvphi(gis-1:gie+1,gje+1,gks) = grvphi(gis-1:gie+1,gje,gks)
      grvphi(gis-1:gie+1,gje+2,gks) = grvphi(gis-1:gie+1,gje-1,gks)
    endif
  end if
  if (gis==gis_global) then
    grvphi(gis-1,:,:) = grvphi(gis,:,:)
    grvphi(gis-2,:,:) = grvphi(gis+1,:,:)
  end if
  if (gbtype==1 .and. gie==gie_global) then
    grvphi(gie+1,:,:)= grvphi(gie,:,:) + &
                    ( grvphi(gie,:,:) - grvphi(gie-1,:,:) ) * dx1(gie)/dx1(gie-1)
    grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                    ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)
  end if

 elseif(crdnt==2.and.dim==3)then ! for spherical coordinates (3D)
  if (gjs==gjs_global) then
    grvphi(gis-1:gie+1,gjs-2,gks:gke) = grvphi(gis-1:gie+1,gjs,gks:gke)
    grvphi(gis-1:gie+1,gjs-1,gks:gke) = grvphi(gis-1:gie+1,gjs,gks:gke)
  endif
  if (gje==gje_global) then
    grvphi(gis-1:gie+1,gje+1,gks:gke) = grvphi(gis-1:gie+1,gje,gks:gke)
    grvphi(gis-1:gie+1,gje+2,gks:gke) = grvphi(gis-1:gie+1,gje-1,gks:gke)
  endif
  if (gis==gis_global) then
    grvphi(gis-1,:,:) = grvphi(gis,:,:)
    grvphi(gis-2,:,:) = grvphi(gis+1,:,:)
  endif
  if (gbtype==1 .and. gie==gie_global) then
    grvphi(gie+1,:,:)= grvphi(gie,:,:) + &
                    ( grvphi(gie,:,:) - grvphi(gie-1,:,:) ) * dx1(gie)/dx1(gie-1)
    grvphi(gie+2,:,:)= grvphi(gie+1,:,:) + &
                    ( grvphi(gie+1,:,:) - grvphi(gie,:,:) ) * dx1(gie+1)/dx1(gie)
  end if
 endif

 if(gravswitch==3)then
  call timestep
  cgrav_old = cgrav
 end if

 call stop_clock(wtelg)
end subroutine gravity_elliptic

end module gravity_elliptic_mod
