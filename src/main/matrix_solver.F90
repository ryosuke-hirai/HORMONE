module matrix_solver

  implicit none

contains

subroutine setup_matrix(system)
  use settings, only: matrix_solver
  use grid, only: is, ie, js, je, ks, ke, gis, gie, gjs, gje, gks, gke
#ifdef USE_PETSC
  use petsc_solver_mod, only: setup_petsc
  use matrix_vars, only: petsc_grv, petsc_rad
#endif
  use miccg_mod, only: setup_cg
  use matrix_vars, only: cg_grv, cg_rad
  use matrix_vars, only: igrv, irad

  integer, intent(in) :: system

#ifndef USE_PETSC
  if (matrix_solver == 1) then
    print *, "Error: matrix_solver=1 (MICCG) is not available. Please compile with USE_PETSC."
    stop
  end if
#endif

  print*, "Setting up", trim(adjustl(merge("MICCG", "PETSc", matrix_solver == 0))), "solver for", trim(adjustl(merge("gravity", "radiation", system == igrv)))

  if (matrix_solver == 0) then
    if (system == igrv) then
      call setup_cg(cg_grv, gis, gie, gjs, gje, gks, gke)
    else if (system == irad) then
      call setup_cg(cg_rad, is, ie, js, je, ks, ke)
    end if
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    if (system == igrv) then
      call setup_petsc(petsc_grv, gis, gie, gjs, gje, gks, gke)
    else if (system == irad) then
      call setup_petsc(petsc_rad, is, ie, js, je, ks, ke)
    end if
#endif
  end if

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
  use settings, only: matrix_solver
  use matrix_vars, only: igrv, cg_grv
  use miccg_mod, only: write_A_cg
#ifdef USE_PETSC
  use matrix_vars, only: petsc_grv
  use petsc_solver_mod, only: write_A_petsc
#endif

  if (matrix_solver == 0) then
    call write_A_cg(cg_grv, igrv)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call write_A_petsc(petsc_grv, igrv)
#endif
  end if

end subroutine write_A_grv

subroutine write_A_rad
  use settings, only: matrix_solver
  use matrix_vars, only: irad, cg_rad
  use miccg_mod, only: write_A_cg
#ifdef USE_PETSC
  use matrix_vars, only: petsc_rad
  use petsc_solver_mod, only: write_A_petsc
#endif

  if (matrix_solver == 0) then
    call write_A_cg(cg_rad, irad)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call write_A_petsc(petsc_rad, irad)
#endif
  endif

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
  use settings, only: matrix_solver
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_grv
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: petsc_grv
#endif

  real(8), allocatable, intent(inout) :: x(:)
  real(8), allocatable, intent(in)    :: b(:)

  if (matrix_solver == 0) then
    call miccg(cg_grv, b, x)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call solve_system_petsc(petsc_grv, x, b)
#endif
  endif

end subroutine solve_system_grv

subroutine solve_system_rad(x, b)
  use settings, only: matrix_solver
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_rad
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: petsc_rad
#endif

  real(8), allocatable, intent(inout) :: x(:)
  real(8), allocatable, intent(in)    :: b(:)

  if (matrix_solver == 0) then
    call miccg(cg_rad, b, x)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call solve_system_petsc(petsc_rad, x, b)
#endif
  endif

end subroutine solve_system_rad

end module matrix_solver