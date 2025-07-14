module matrix_solver_mod

  implicit none

  public :: setup_matrix, write_A_grv, write_A_rad
  public :: solve_system_grv, solve_system_rad

contains


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE SETUP_MATRIX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Setup the matrix solver (MICCG or PETSc) for either gravity or radiation

subroutine setup_matrix(system)
  use settings, only: matrix_solver
  use mpi_utils, only: nprocs, stop_mpi, myrank
  use grid, only: is, ie, js, je, ks, ke, gis, gie, gjs, gje, gks, gke, &
                  is_global, ie_global, js_global, je_global, ks_global, &
                  gis_global, gie_global, gjs_global, gje_global, gks_global

#ifdef USE_PETSC
  use grid, only: ke_global, gke_global
  use petsc_solver_mod, only: setup_petsc
  use matrix_vars, only: petsc_grv, petsc_rad
#endif
  use miccg_mod, only: setup_cg
  use matrix_utils, only: contiguous_map
  use matrix_vars, only: cg_grv, cg_rad
  use matrix_vars, only: igrv, irad, lmax_grv, lmax_rad, map_grv, map_rad

  integer, intent(in) :: system

#ifndef USE_PETSC
  if (matrix_solver == 1) then
    if (myrank == 0) then
      print*, "Error: matrix_solver=1 (PETSc) is not available. Please compile with USE_PETSC."
    end if
    call stop_mpi(1)
  end if
#endif

  ! MICCG solver does not work with MPI
  if (matrix_solver == 0 .and. (nprocs >= 2)) then
    if (myrank == 0) then
      print*, "Error: MICCG solver is not compatible with MPI. Please use PETSc (matrix_solver=1)."
    end if
    call stop_mpi(1)
  end if

  if (myrank == 0) then
    print*, "Setting up ", trim(adjustl(merge("MICCG", "PETSc", matrix_solver == 0))), " solver for ", trim(adjustl(merge("gravity  ", "radiation", system == igrv)))
    print*, ""
  end if

  select case (system)
  case (igrv)
    ! Set up solver for gravity
    select case (matrix_solver)
    case (0)
      call setup_cg(gis, gie, gjs, gje, gks, gke, cg_grv)
    case (1)
#ifdef USE_PETSC
      call setup_petsc(gis, gie, gjs, gje, gks, gke, gis_global, &
                       gie_global, gjs_global, gje_global, gks_global, gke_global, &
                       petsc_grv)
#endif
    end select

    ! The local number of rows on the matrix for gravity
    lmax_grv = (gie - gis + 1) * (gje - gjs + 1) * (gke - gks + 1)

    ! Create the mapping from local to global indices for gravity
    if (allocated(map_grv)) deallocate(map_grv)
    call contiguous_map(gis, gie, gjs, gje, gks, gke, gis_global, gie_global, &
                       gjs_global, gje_global, gks_global, map_grv)

  case (irad)
    ! Set up solver for radiation
    select case (matrix_solver)
    case (0)
      call setup_cg(is, ie, js, je, ks, ke, cg_rad)
    case (1)
#ifdef USE_PETSC
      call setup_petsc(is, ie, js, je, ks, ke, &
                       is_global, ie_global, js_global, je_global, ks_global, ke_global, &
                       petsc_rad)
#endif
    end select

    ! The local number of rows on the matrix for radiation
    lmax_rad = (ie - is + 1) * (je - js + 1) * (ke - ks + 1)

    ! Create the mapping from local to global indices for radiation
    if (allocated(map_rad)) deallocate(map_rad)
    call contiguous_map(is, ie, js, je, ks, ke, is_global, ie_global, &
                      js_global, je_global, ks_global, map_rad)
  end select

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
    call write_A_cg(igrv, cg_grv)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call write_A_petsc(igrv, petsc_grv)
#endif
  end if

end subroutine write_A_grv

subroutine write_A_rad
  use settings, only: matrix_solver
  use matrix_vars, only: irad, cg_rad
  use miccg_mod, only: write_A_cg
  use profiler_mod
#ifdef USE_PETSC
  use matrix_vars, only: petsc_rad
  use petsc_solver_mod, only: write_A_petsc
#endif

  if (matrix_solver == 0) then
    call start_clock(wtmra)
    call write_A_cg(irad, cg_rad)
    call stop_clock(wtmra)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call start_clock(wtpra)
    call write_A_petsc(irad, petsc_rad)
    call stop_clock(wtpra)
#endif
  endif

end subroutine write_A_rad


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE SOLVE_SYSTEM_GRV/RAD
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Solve the system of equations for gravity or radiation
!          It calls either the MICCG or PETSc branch depending on the
!          compilation flag USE_PETSC.

subroutine solve_system_grv(b, x)
  use settings, only: matrix_solver
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_grv
  use profiler_mod
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: petsc_grv
#endif

  real(8), allocatable, intent(in)    :: b(:)
  real(8), allocatable, intent(inout) :: x(:)

  if (matrix_solver == 0) then
    call start_clock(wtmig)
    call miccg(cg_grv, b, x)
    call stop_clock(wtmig)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call start_clock(wtpeg)
    call solve_system_petsc(petsc_grv, b, x)
    call stop_clock(wtpeg)
#endif
  endif

end subroutine solve_system_grv

subroutine solve_system_rad(b, x)
  use settings, only: matrix_solver
  use miccg_mod, only: miccg
  use matrix_vars, only: cg_rad
  use profiler_mod
#ifdef USE_PETSC
  use petsc_solver_mod, only: solve_system_petsc
  use matrix_vars, only: petsc_rad
#endif

  real(8), allocatable, intent(in)    :: b(:)
  real(8), allocatable, intent(inout) :: x(:)

  if (matrix_solver == 0) then
    call start_clock(wtmir)
    call miccg(cg_rad, b, x)
    call stop_clock(wtmir)
  else if (matrix_solver == 1) then
#ifdef USE_PETSC
    call start_clock(wtper)
    call solve_system_petsc(petsc_rad, b, x)
    call stop_clock(wtper)
#endif
  endif

end subroutine solve_system_rad

end module matrix_solver_mod