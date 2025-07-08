module matrix_vars

#ifdef USE_PETSC
#include <petsc/finclude/petsc.h>
    use petsc
#endif

  implicit none

  ! Integer switch for the gravity and radiation solvers
  integer, parameter :: igrv = 1
  integer, parameter :: irad = 2

  ! A method-independent variable to store the number of rows owned
  integer :: lmax_rad
  integer :: lmax_grv

  ! Mapping arrays for global indexing (shared by all methods)
  integer, allocatable :: map_grv(:)
  integer, allocatable :: map_rad(:)

  ! PETSc arrays
#ifdef USE_PETSC
  type petsc_set
    integer :: dim
    integer :: is, ie, js ,je, ks, ke, in, jn, kn
    integer :: is_global, ie_global, in_global
    integer :: js_global, je_global, jn_global
    integer :: ks_global, ke_global, kn_global
    integer, allocatable :: my_rows(:)
    Mat :: A
    Vec :: x
    Vec :: b
    Vec :: x_gather ! Local vector for gathering solution
    KSP :: ksp
    IS :: is_from   ! Index set for scatter
    VecScatter :: scatter_ctx ! Scatter context for solution vector
    PetscInt :: lmax, lmax_global
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