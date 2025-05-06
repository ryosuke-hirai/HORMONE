module matrix_solver

  implicit none

contains

subroutine setup_matrix(system)
  use grid, only: is, ie, js, je, ks, ke, gis, gie, gjs, gje, gks, gke
#ifdef USE_PETSC
  use petsc_solver_mod, only: setup_petsc
  use matrix_vars, only: petsc_grv, petsc_rad
#else
  use miccg_mod, only: setup_cg
  use matrix_vars, only: cg_grv, cg_rad
#endif
  use matrix_vars, only: igrv, irad

  integer, intent(in) :: system

#ifdef USE_PETSC

  ! Set up the PETSc matrix and vectors
  if (system == igrv) then
    call setup_petsc(petsc_grv, gis, gie, gjs, gje, gks, gke)
  else if (system == irad) then
    call setup_petsc(petsc_rad, is, ie, js, je, ks, ke)
  end if
#else
  if (system == igrv) then
    call setup_cg(cg_grv, gis, gie, gjs, gje, gks, gke)
  else if (system == irad) then
    call setup_cg(cg_rad, is, ie, js, je, ks, ke)
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
  use matrix_vars, only: petsc_grv
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(petsc_grv, igrv)
#else
  use matrix_vars, only: cg_grv
  use miccg_mod, only: write_A_cg
  call write_A_cg(cg_grv, igrv)
#endif

end subroutine write_A_grv

subroutine write_A_rad
  use matrix_vars, only: irad
#ifdef USE_PETSC
  use matrix_vars, only: petsc_rad
  use petsc_solver_mod, only: write_A_petsc
  call write_A_petsc(petsc_rad, irad)
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
  use matrix_vars, only: petsc_grv
#else
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_grv
#endif

  real(8), allocatable, intent(inout) :: x(:)
  real(8), allocatable, intent(in)    :: b(:)

#ifdef USE_PETSC
  ! PETSc solver for gravity system
  call solve_system_petsc(petsc_grv, x, b)
#else
  ! MICCG solver for gravity system
  call miccg(cg_grv, b, x)
#endif

end subroutine solve_system_grv

subroutine solve_system_rad(x, b)
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: petsc_rad
#else
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_rad
#endif

real(8), allocatable, intent(inout) :: x(:)
real(8), allocatable, intent(in)    :: b(:)

#ifdef USE_PETSC
  ! PETSc solver for radiation system
call solve_system_petsc(petsc_rad, x, b)
#else
  ! MICCG solver for radiation system
  call miccg(cg_rad, b, x)
#endif

end subroutine solve_system_rad

end module matrix_solver