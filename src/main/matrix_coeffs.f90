module matrix_coeffs_mod

  implicit none
  private

  public :: compute_coeffs_gravity

  contains


  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !
  !                     SUBROUTINE COMPUTE_COEFFS_GRAVITY
  !
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ! PURPOSE: Computes the stencil coefficients for a given grid point.

  subroutine compute_coeffs_gravity (dim, i, j, k, ie, ke, coeffs)
    use settings, only: crdnt, eq_sym, gbtype
    use grid,only:gis_global,gie_global,gje_global,gks_global,gke_global,&
    xi1s,x1,xi1,dx1,idx1,dxi1,dx2,dxi2,dx3,dxi3,idx3,&
    sini,sinc,rdis

    integer, intent(in) :: dim, i, j, k, ie, ke
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
      if ( i == ie ) then ! TODO: check MPI
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
          if ( i == ie ) then ! TODO: check MPI
            coeffs(1) = coeffs(1) + coeffs(2) * rdis(i, k)/rdis(i+1, k)
          end if
          if ( k == ke ) then ! TODO: check MPI
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

end module matrix_coeffs_mod