module matrix_solver
  use miccg_mod, only: ijk_from_l, get_preconditioner
  use matrix_vars, only: cg_set, igrv, irad, lmax_grv, lmax_rad

  implicit none

  contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SETUP_A_GRV/RAD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Top-level assembly routine for the matrix A
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.

! TODO: rename this to write_A_grv/rad
subroutine setup_A_grv
#ifdef USE_PETSC
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(igrv)
#else
  use miccg_mod, only: setup_cg, write_A_cg
  ! Setup is performed here
  call setup_cg(igrv)
  call write_A_cg(igrv)
#endif

end subroutine setup_A_grv

subroutine setup_A_rad
#ifdef USE_PETSC
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(irad)
#else
  use miccg_mod, only: write_A_cg
  ! Setup already performed at the start of the program
  call write_A_cg(irad)
#endif
end subroutine setup_A_rad


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SOLVE_SYSTEM_GRV/RAD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solve the system of equations for gravity or radiation
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.

subroutine solve_system_grv(cgsrc, x)
  use matrix_vars, only: A_petsc_grv, cg_grv
  use petsc_solver_mod, only: solve_system_petsc
  use miccg_mod, only: miccg
  real(8), allocatable, intent(in) :: cgsrc(:)
  real(8), allocatable, intent(inout) :: x(:)

#ifdef USE_PETSC
  ! PETSc solver for gravity system
  call solve_system_petsc(A_petsc_grv, cgsrc, x)
#else
  ! MICCG solver for gravity system
  call miccg(cg_grv, cgsrc, x)
#endif

end subroutine solve_system_grv

subroutine solve_system_rad(cgsrc, x)
  use matrix_vars, only: A_petsc_rad, cg_rad
  use petsc_solver_mod, only: solve_system_petsc
  use miccg_mod, only: miccg
  real(8), allocatable, intent(in) :: cgsrc(:)
  real(8), allocatable, intent(inout) :: x(:)

#ifdef USE_PETSC
  ! PETSc solver for radiation system
  call solve_system_petsc(A_petsc_rad, cgsrc, x)
#else
  ! MICCG solver for radiation system
  call miccg(cg_rad, cgsrc, x)
#endif

end subroutine solve_system_rad

end module matrix_solver