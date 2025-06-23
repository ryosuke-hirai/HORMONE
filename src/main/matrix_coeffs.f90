module matrix_coeffs

  implicit none
  private

  public :: compute_coeffs, get_matrix_offsets

  contains

  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                     SUBROUTINE COMPUTE_COEFFS
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: Computes the stencil coefficients for a given grid point.

  subroutine compute_coeffs(system, dim, i, j, k, coeffs)
    use matrix_vars, only: igrv, irad
    use grid, only: is_global, ie_global, js_global, je_global, ks_global, ke_global
    integer, intent(in) :: system, dim, i, j, k
    real(8), intent(out) :: coeffs(5)

    integer :: raddim, in, jn, kn

    if (system == igrv) then
      ! Compute coefficients for gravity system
      call compute_coeffs_gravity(dim, i, j, k, coeffs)
    else if (system == irad) then
      ! raddim specfies the dimension and plane of the grid
      ! TODO: move this up one level so that it's not computed for every point
      in = ie_global - is_global + 1
      jn = je_global - js_global + 1
      kn = ke_global - ks_global + 1
      if(in>1.and.jn>1.and.kn>1)then
        raddim=3
      elseif(in>1.and.jn>1.and.kn==1)then
        raddim=21
      elseif(in>1.and.jn==1.and.kn>1)then
        raddim=22
      elseif(in==1.and.jn>1.and.kn>1)then
        raddim=23
      elseif(in>1.and.jn==1.and.kn==1)then
        raddim=11
      elseif(in==1.and.jn>1.and.kn==1)then
        raddim=12
      elseif(in==1.and.jn==1.and.kn>1)then
        raddim=13
      else
        print*,'Error in compute_coeffs, dimension is not supported; raddim=',raddim
      end if

      ! Compute coefficients for radiation system
      call compute_coeffs_radiation(raddim, i, j, k, coeffs)
    end if
  end subroutine compute_coeffs


  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                     SUBROUTINE COMPUTE_COEFFS_GRAVITY
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: Computes the stencil coefficients for gravity

  subroutine compute_coeffs_gravity(dim, i, j, k, coeffs)
    use settings, only: crdnt, eq_sym, gbtype
    use grid,only:gis_global,gie_global,gje_global,gks_global,gke_global,&
    xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,dx3,dxi3,idx3,&
    sini,sinc,rdis

    integer, intent(in) :: dim, i, j, k
    real(8), intent(out), dimension(:) :: coeffs
    real(8) :: sum_dx3

    select case(dim)
    case(1)  ! 1D spherical coordinates (crdnt == 2 assumed)
      ! Expect coeffs(1:2)
      coeffs(1) = - ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) )
      coeffs(2) = xi1(i)**2/dx1(i+1)
      if ( i == gis_global .and. xi1s > 0d0 ) then
        coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j) * dxi2(j) * idx1(i)
      end if
      if ( gbtype == 1 .and. i == gie_global ) then
        coeffs(1) = coeffs(1) + coeffs(2) * x1(i)/x1(i+1)
      end if
      if ( i == gie_global ) then
        coeffs(2) = 0d0
      end if

    case(2)  ! 2D
      if ( crdnt == 1 ) then
        ! 2D cylindrical coordinates
        ! Expect coeffs(1:3)
        ! Use sum_dx3 = dx3(k) + dx3(k+1)
        sum_dx3 = dx3(k) + dx3(k+1)
        coeffs(1) = - ( xi1(i)/dx1(i+1) + xi1(i-1)/dx1(i) + &
        2d0*x1(i)*dxi1(i)/( dx3(k)*dx3(k+1) ) ) * 0.5d0 * sum_dx3
        coeffs(2) = 0.5d0 * xi1(i) * sum_dx3 / dx1(i+1)
        coeffs(3) = x1(i) * dxi1(i) / dx3(k+1)
        if ( gbtype == 1 ) then
          if ( i == gie_global ) then
            coeffs(1) = coeffs(1) + coeffs(2) * rdis(i, k)/rdis(i+1, k)
          end if
          if ( k == gke_global ) then
            coeffs(1) = coeffs(1) + coeffs(2) * rdis(i, k)/rdis(i, k+1)
          end if
        end if
        if ( i == gie_global ) then
          coeffs(2) = 0d0
        end if
        if ( k == gke_global ) then
          coeffs(3) = 0d0
        end if
        if ( eq_sym .and. k == gks_global ) then
          coeffs(1) = coeffs(1) + x1(i)* dxi1(i) / dx3(k+1)
        end if
      else if ( crdnt == 2 ) then
        ! 2D spherical coordinates
        ! Expect coeffs(1:3)
        coeffs(1) = - ( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) * sinc(j)* dxi2(j) + &
        ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j) ) * dxi1(i) )
        coeffs(2) = xi1(i)**2 * sinc(j)* dxi2(j) / dx1(i+1)
        coeffs(3) = sini(j) * dxi1(i) / dx2(j+1)
        if ( i == gis_global .and. xi1s > 0d0 ) then
          coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j)* dxi2(j)* idx1(i)
        end if
        if ( gbtype == 1 .and. i == gie_global ) then
          coeffs(1) = coeffs(1) + coeffs(2)* x1(i)/x1(i+1)
        end if
        if ( i == gie_global ) then
          coeffs(2) = 0d0
        end if
        if ( j == gje_global ) then
          coeffs(1) = coeffs(1) + sini(j) * dxi1(i) / dx2(j+1)
          coeffs(3) = 0d0
        end if
      end if

    case(3)  ! 3D spherical coordinates (crdnt == 2 assumed)
      ! Expect coeffs(1:5)
      coeffs(1) = - ( ( xi1(i)**2/dx1(i+1) + xi1(i-1)**2/dx1(i) ) * sinc(j)* dxi2(j)* dxi3(k) + &
      ( sini(j)/dx2(j+1) + sini(j-1)/dx2(j) ) * dxi1(i)* dxi3(k) + &
      ( idx3(k) + idx3(k-1) ) * dxi1(i)* dxi2(j) )
      coeffs(2) = xi1(i)**2 * sinc(j)* dxi2(j)* dxi3(k) / dx1(i+1)
      coeffs(3) = sini(j) * dxi1(i)* dxi3(k) / dx2(j+1)
      coeffs(4) = dxi1(i)* dxi2(j) / dx3(k+1)
      coeffs(5) = dxi1(i)* dxi2(j) / dx3(k)
      if ( eq_sym .and. j == gje_global ) then
        coeffs(1) = coeffs(1) + coeffs(3)
      end if
      if ( i == gis_global .and. xi1s > 0d0 ) then
        coeffs(1) = coeffs(1) + xi1(i-1)**2 * sinc(j)* dxi2(j)* idx1(i)
      end if
      if ( gbtype == 1 .and. i == gie_global ) then
        coeffs(1) = coeffs(1) + coeffs(2)* x1(i)/x1(i+1)
      end if
      if ( i == gie_global ) then
        coeffs(2) = 0d0
      end if
      if ( j == gje_global ) then
        coeffs(3) = 0d0
      end if
      if ( k == gke_global ) then
        coeffs(4) = 0d0
      end if
      if ( k /= gks_global ) then
        coeffs(5) = 0d0
      end if

    end select

  end subroutine compute_coeffs_gravity

  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                     SUBROUTINE COMPUTE_COEFFS_RADIATION
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: Computes the stencil coefficients for radiation

  subroutine compute_coeffs_radiation(raddim, i, j, k, coeffs)
    use settings, only: crdnt, radswitch
    use constants,only: arad, clight, Cv
    use utils, only: har_mean
    use radiation_utils, only: geo
    use opacity_mod, only: kappa_p
    use grid, only: dvol, dt, is_global, ie_global, js_global, je_global, ks_global, ke_global
    use physval, only: T, imu, radK, d

    integer, intent(in) :: raddim, i, j, k
    real(8), intent(out) :: coeffs(:)

    real(8) :: kappap

    select case(raddim)
    case (11)  ! 1D x-direction
      coeffs(1) = dvol(i,j,k)/dt
      if(i > is_global) then
         coeffs(1) = coeffs(1) + geo(1, i-1, j, k)*har_mean(radK(i-1:i, j, k))
      end if
      if(i < ie_global) then
         coeffs(1) = coeffs(1) + geo(1, i  , j, k)*har_mean(radK(i:i+1, j, k))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      if(i == ie_global) coeffs(2) = 0d0

    case (12)  ! 1D y-direction
      coeffs(1) = dvol(i,j,k)/dt
      if(j > js_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j-1, k)*har_mean(radK(i, j-1:j, k))
      end if
      if(j < je_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j  , k)*har_mean(radK(i, j:j+1, k))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      if(j == je_global) coeffs(2) = 0d0

    case (13)  ! 1D z-direction
      coeffs(1) = dvol(i,j,k)/dt
      if(k > ks_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k-1)*har_mean(radK(i, j, k-1:k))
      end if
      if(k < ke_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k  )*har_mean(radK(i, j, k:k+1))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      if(k == ke_global) coeffs(2) = 0d0

    case (21)  ! 2D xy-plane
      coeffs(1) = dvol(i,j,k)/dt
      if(i > is_global) then
         coeffs(1) = coeffs(1) + geo(1, i-1, j, k)*har_mean(radK(i-1:i, j, k))
      end if
      if(i < ie_global) then
         coeffs(1) = coeffs(1) + geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      end if
      if(j > js_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j-1, k)*har_mean(radK(i, j-1:j, k))
      end if
      if(j < je_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      if(i == ie_global) coeffs(2) = 0d0
      coeffs(3) = - geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      if(j == je_global) coeffs(3) = 0d0

    case (22)  ! 2D xz-plane
      coeffs(1) = dvol(i,j,k)/dt
      if(i > is_global) then
         coeffs(1) = coeffs(1) + geo(1, i-1, j, k)*har_mean(radK(i-1:i, j, k))
      end if
      if(i < ie_global) then
         coeffs(1) = coeffs(1) + geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      end if
      if(k > ks_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k-1)*har_mean(radK(i, j, k-1:k))
      end if
      if(k < ke_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      if(i == ie_global) coeffs(2) = 0d0
      coeffs(3) = - geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      if(k == ke_global) coeffs(3) = 0d0

    case (23)  ! 2D yz-plane
      coeffs(1) = dvol(i,j,k)/dt
      if(j > js_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j-1, k)*har_mean(radK(i, j-1:j, k))
      end if
      if(j < je_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      end if
      if(k > ks_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k-1)*har_mean(radK(i, j, k-1:k))
      end if
      if(k < ke_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      if(j == je_global) coeffs(2) = 0d0
      coeffs(3) = - geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      if(k == ke_global) coeffs(3) = 0d0

    case (3)  ! 3D case
      coeffs(1) = dvol(i,j,k)/dt
      if(i > is_global) then
         coeffs(1) = coeffs(1) + geo(1, i-1, j, k)*har_mean(radK(i-1:i, j, k))
      end if
      if(i < ie_global) then
         coeffs(1) = coeffs(1) + geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      end if
      if(j > js_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j-1, k)*har_mean(radK(i, j-1:j, k))
      end if
      if(j < je_global) then
         coeffs(1) = coeffs(1) + geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      end if
      if(k > ks_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k-1)*har_mean(radK(i, j, k-1:k))
      end if
      if(k < ke_global) then
         coeffs(1) = coeffs(1) + geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      end if
      if(radswitch == 1) then
         kappap = kappa_p(d(i,j,k), T(i,j,k))
         coeffs(1) = coeffs(1) - dvol(i,j,k)*clight*kappap*d(i,j,k) * &
              ( 4d0*arad*T(i,j,k)**3*kappap*clight*dt / (Cv*imu(i,j,k) + &
                4d0*arad*kappap*clight*dt*T(i,j,k)**3) - 1d0 )
      end if
      coeffs(2) = - geo(1, i, j, k)*har_mean(radK(i:i+1, j, k))
      if(i == ie_global) coeffs(2) = 0d0
      coeffs(3) = - geo(2, i, j, k)*har_mean(radK(i, j:j+1, k))
      if(j == je_global) coeffs(3) = 0d0
      coeffs(4) = - geo(3, i, j, k)*har_mean(radK(i, j, k:k+1))
      if(k == ke_global) coeffs(4) = 0d0
      if(crdnt == 2) then
         if(k == 1) then
            ! radK is computed for the first layer of ghost cells, so we can use k-1
            ! which is equivalent to ke_global
            coeffs(5) = - geo(3, i, j, k-1)*har_mean( (/ radK(i,j,k), radK(i,j,k-1) /) )
         else
            coeffs(5) = 0d0
         end if
      end if

    end select

  end subroutine compute_coeffs_radiation

  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                     SUBROUTINE get_matrix_offsets
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: Computes the diagonal offsets for the matrix based on the dimension.
  !          Same for gravity and radiation.

  subroutine get_matrix_offsets(dim, offsets)
    use grid, only: is_global, ie_global, js_global, je_global, ks_global, ke_global
    use settings, only: crdnt
    integer, intent(in)  :: dim
    integer, intent(out) :: offsets(:)

    integer :: in, jn, kn

    ! TODO: store globally
    in = ie_global - is_global + 1
    jn = je_global - js_global + 1
    kn = ke_global - ks_global + 1

    ! As the dimension increases, additional diagonals are added, but the
    ! offsets of existing diagonals are not changed.
    ! In 1D, the offsets are at 0, 1.
    ! In 2D, the offsets are at 0, 1, in
    ! In 3D, the offsets are at 0, 1, in, in*jn, in*jn*(kn-1)
    ! These are the same for gravity and radiation
    ! Special case of 2D, which depends on the plane of the grid.

    offsets(1) = 0
    offsets(2) = 1

    if (dim == 2) then
      if (in > 1) then
        offsets(3) = in
      else
        offsets(3) = jn
      end if
    endif

    if (dim == 3) then
      offsets(3) = in
      offsets(4) = in * jn
      ! The 5th diagonal is only needed for periodic boundaries in spherical coordinates
      if (crdnt == 2) then
        offsets(5) = in * jn * (kn - 1)
      endif
    endif

    ! TODO: Future implementation should check boundary conditions individually:
    ! if( bc1is==0 )then
    !  n=n+1
    !  offsets(n) = in-1
    ! end if
    ! if( bc2is==0 )then
    !  n=n+1
    !  offsets(n) = in * (jn-1)
    ! end if
    ! if( bc3is==0 )then
    !  n=n+1
    !  offsets(n) = in * jn * (kn-1)
    ! end if
    ! This would also require setting cg%Adiags in this subroutine.

  end subroutine get_matrix_offsets


end module matrix_coeffs
