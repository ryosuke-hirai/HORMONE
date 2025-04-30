module matrix_vars
  implicit none

#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
    use petsc
#endif

  ! Integer switch for the gravity and radiation solvers
  integer, parameter :: igrv = 1
  integer, parameter :: irad = 2

  ! Matrix sizes
  integer :: lmax_grv
  integer :: lmax_rad

  ! PETSc arrays
#ifdef USE_PETSC
  Mat :: A_petsc_grv
  Mat :: A_petsc_rad
#else
  ! Dummy declarations for non-PETSc builds
  integer :: A_petsc_grv
  integer :: A_petsc_rad
#endif

  ! MICCG arrays
  type cg_set
    integer :: is, ie, js ,je, ks, ke, in, jn, kn, lmax
    integer :: cdiags, Adiags
    integer, allocatable :: ia(:), ic(:)
    real(8), allocatable, dimension(:,:) :: A, c
    real(8) :: alpha
  end type cg_set
  type(cg_set), target :: cg_grv
  type(cg_set), target :: cg_rad

end module matrix_vars