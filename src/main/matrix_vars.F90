module matrix_vars

#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
    use petsc
#endif

  implicit none

  ! Integer switch for the gravity and radiation solvers
  integer, parameter :: igrv = 1
  integer, parameter :: irad = 2

  ! Matrix sizes for allocatioon
  ! Should mirror PETSc and CG lmax
  integer :: lmax_grv
  integer :: lmax_rad

  ! PETSc arrays
#ifdef USE_PETSC
  Mat :: A_petsc_grv
  Mat :: A_petsc_rad
  Vec :: x_petsc_grv
  Vec :: x_petsc_rad
  Vec :: b_petsc_grv
  Vec :: b_petsc_rad
  KSP :: ksp_grv
  KSP :: ksp_rad
  PetscInt :: lmax_petsc_grv
  PetscInt :: lmax_petsc_rad
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