module hydro_mod
 implicit none

contains


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE HYDRO_STEP
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To solve hydrodynamic equations for one time step

subroutine hydro_step

 use settings,only:dirichlet_on,rktype,radswitch
 use grid,only:rungen
 use boundary_mod,only:boundarycondition
 use numflux_mod,only:numflux
 use source_mod,only:source
 use rungekutta_mod,only:rungekutta
 use dirichlet_mod,only:dirichletbound
 use shockfind_mod,only:shockfind
 use mpi_domain,only:exchange_mpi
 use profiler_mod,only:start_clock,stop_clock,wthyd

!-----------------------------------------------------------------------------

 call start_clock(wthyd)

 if(dirichlet_on) call dirichletbound
 call shockfind

 do rungen = 1, rktype
  call exchange_mpi
  call boundarycondition
  call numflux
  call source
  call rungekutta
 end do

 ! if solving radiation, the boudary and ghost cells need to be updated
  if (radswitch > 0) then
    call exchange_mpi
    call boundarycondition
  endif

 call stop_clock(wthyd)

 return
end subroutine hydro_step

end module hydro_mod
