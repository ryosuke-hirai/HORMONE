python -m numpy.f2py -h horpy.pyf -m horpy horpy.f90
cd ../../src/main
python -m numpy.f2py -c ../../work/scripts/horpy.pyf -m horpy --fcompiler=gfortran --f90flags='-fopenmp' -lgomp mpi_utils.F90 modules.f90 profiler.f90 mpi_domain.F90 io.F90 utils.f90 conserve.f90 ionization.f90 eos.f90 matrix_vars.F90 opacity.f90 radiation_utils.f90 matrix_coeffs.f90 matrix_utils.f90 miccg.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 petsc_solver.F90 matrix_solver.F90 radiation.f90 cooling.f90 sinks.f90 externalforce.f90 composition.f90 star.f90 shockfind.f90 output.f90 smear.f90 input.f90 setup.f90 readbin.f90 timestep.f90 tests.f90 checksetup.f90 ../../work/scripts/horpy.f90

mv horpy*so ../../work/horpy.so
