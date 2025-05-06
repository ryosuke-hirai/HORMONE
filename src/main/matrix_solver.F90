module matrix_solver

  implicit none

contains

subroutine setup_matrix(system)

#ifdef USE_PETSC
  use petsc_solver_mod, only: setup_petsc
  use matrix_vars, only: A_petsc_grv, A_petsc_rad, b_petsc_grv, b_petsc_rad, &
  x_petsc_grv, x_petsc_rad, ksp_grv, ksp_rad, lmax_petsc_grv, lmax_petsc_rad
#else
  use miccg_mod, only: setup_cg
  use matrix_vars, only: cg_grv, cg_rad
#endif
  use matrix_vars, only: lmax_grv, lmax_rad, igrv, irad

  integer, intent(in) :: system

#ifdef USE_PETSC

  ! Set up the PETSc matrix and vectors
  if (system == igrv) then
    call setup_petsc(A_petsc_grv, b_petsc_grv, x_petsc_grv, ksp_grv, lmax_petsc_grv)
    lmax_grv = lmax_petsc_grv
  else if (system == irad) then
    call setup_petsc(A_petsc_rad, b_petsc_rad, x_petsc_rad, ksp_rad, lmax_petsc_rad)
    lmax_rad = lmax_petsc_rad
  end if
#else
  if (system == igrv) then
    call setup_cg(cg_grv)
    lmax_grv = cg_grv%lmax
  else if (system == irad) then
    call setup_cg(cg_rad)
    lmax_rad = cg_rad%lmax
  end if

#endif

end subroutine setup_matrix
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE WRITE_A_GRV/RAD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Top-level assembly routine for the matrix A
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.

subroutine write_A_grv
  use matrix_vars, only: igrv
#ifdef USE_PETSC
  use matrix_vars, only: A_petsc_grv, lmax_petsc_grv, ksp_grv
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(A_petsc_grv, ksp_grv, lmax_petsc_grv, igrv)
#else
  use matrix_vars, only: cg_grv
  use miccg_mod, only: write_A_cg
  call write_A_cg(cg_grv, igrv)
#endif

end subroutine write_A_grv

subroutine write_A_rad
  use matrix_vars, only: irad
#ifdef USE_PETSC
  use matrix_vars, only: A_petsc_rad, lmax_petsc_rad, ksp_rad
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(A_petsc_rad, ksp_rad, lmax_petsc_rad, irad)
#else
  use matrix_vars, only: cg_rad
  use miccg_mod, only: write_A_cg
  call write_A_cg(cg_rad, irad)
#endif
end subroutine write_A_rad


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                          SUBROUTINE SOLVE_SYSTEM_GRV/RAD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solve the system of equations for gravity or radiation
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.

subroutine solve_system_grv(x, b)
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: x_petsc_grv, b_petsc_grv, ksp_grv, lmax_petsc_grv
#else
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_grv
#endif

  real(8), allocatable, intent(inout) :: x(:)
  real(8), allocatable, intent(in)    :: b(:)

#ifdef USE_PETSC
  ! PETSc solver for gravity system
  call solve_system_petsc(x_petsc_grv, b_petsc_grv, ksp_grv, lmax_petsc_grv, x, b)
#else
  ! MICCG solver for gravity system
  call miccg(cg_grv, b, x)
#endif

end subroutine solve_system_grv

subroutine solve_system_rad(x, b)
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: x_petsc_rad, b_petsc_rad, ksp_rad, lmax_petsc_rad
#else
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_rad
#endif

real(8), allocatable, intent(inout) :: x(:)
real(8), allocatable, intent(in)    :: b(:)

#ifdef USE_PETSC
  ! PETSc solver for radiation system
call solve_system_petsc(x_petsc_rad, b_petsc_rad, ksp_rad, lmax_petsc_rad, x, b)
#else
  ! MICCG solver for radiation system
  call miccg(cg_rad, b, x)
#endif

end subroutine solve_system_rad

end module matrix_solver