module matrix_vars

#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
    use petsc
#endif

  implicit none

  ! Integer switch for the gravity and radiation solvers
  integer, parameter :: igrv = 1
  integer, parameter :: irad = 2

  ! PETSc arrays
#ifdef USE_PETSC
  type petsc_set
    integer :: is, ie, js ,je, ks, ke, in, jn, kn, dim
    Mat :: A
    Vec :: x
    Vec :: b
    KSP :: ksp
    PetscInt :: lmax
  end type petsc_set
  type(petsc_set), target :: petsc_grv
  type(petsc_set), target :: petsc_rad
#endif

  ! MICCG arrays
  type cg_set
    integer :: is, ie, js ,je, ks, ke, in, jn, kn, lmax, dim
    integer :: cdiags, Adiags
    integer, allocatable :: ia(:), ic(:)
    real(8), allocatable, dimension(:,:) :: A, c
    real(8) :: alpha
  end type cg_set
  type(cg_set), target :: cg_grv
  type(cg_set), target :: cg_rad

end module matrix_vars