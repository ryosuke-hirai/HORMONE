# Some versions of make default FC to f77
ifeq ($(FC),f77)
    FC = gfortran
endif
FC ?= gfortran
FCFLAG=
FPPFLAG=

# for gfortran
ifeq (${FC},gfortran)
ifeq (${MPI},yes)
override FC := mpifort
FPPFLAG += -DMPI
endif
MODFLAG= -J
FCFLAG:= $(CFLAG) -fopenmp -ffree-line-length-512 #-fmax-stack-var-size=32768
DBFLAG= -g -ffpe-trap=invalid,zero,overflow -fbacktrace -O0 -Wall -Wextra -fcheck=all
PRFLAG= -pg -save-temps
endif

# for Intel fortran compiler
ifeq (${FC},ifx)
ifeq (${MPI},yes)
override FC := mpiifx
FPPFLAG += -DMPI
endif
MODFLAG= -module # trailing space is important
FCFLAG:= $(CFLAG) -qopenmp -fp-model precise #-heap-arrays
DBFLAG= -traceback -fpe0 -CB -warn all -O0 -std #-pg #-check all
PRFLAG= -pg
endif

# debug mode
ifeq (${DEBUG},yes)
FCFLAG+=$(DBFLAG)
endif

# profile mode
ifeq (${PROFILE},yes)
FCFLAG+=$(PRFLAG)
endif

# PETSc settings
ifeq (${USE_PETSC},yes)
# Useful as an escape hatch if libpetsc.so was not compiled properly
# and has missing symbols from external libraries.
ifeq ($(PETSC_PC_STATIC), yes)
STATIC = "--static"
else
STATIC =
endif

PETSC_PC = $(shell pkg-config --path petsc || pkg-config --path PETSc || if [ -n "${PETSC_DIR}" ] && [ -d "${PETSC_DIR}" ]; then find "${PETSC_DIR}" -type f -iname petsc.pc 2>/dev/null | head -n 1; fi)
PETSC_INCLUDE ?= $(shell pkg-config --cflags $(PETSC_PC) 2>/dev/null)
PETSC_LDFLAGS ?= $(shell pkg-config $(STATIC) --libs $(PETSC_PC) 2>/dev/null)

LDFLAGS += $(PETSC_LDFLAGS)
FPPFLAG += -DUSE_PETSC $(PETSC_INCLUDE)
endif

# directories
OBJ_DIR = obj
SRC_DIR = ../src/main
INI_DIR = ../src/init
POS_DIR = ../src/post
BIN_DIR = $(HORMONE_DIR)/bin

# installer settings
SCRIPTS_DIR := ./scripts
SED_INPLACE := $(shell if sed --version >/dev/null 2>&1; then echo "sed -i"; else echo "sed -i ''"; fi)

TARGET = hormone

TARGETA = ascii_output

TARGETL = lightcurve

# compile modules that are being depended on first
MOD = mpi_utils.F90 derived_types.f90 modules.f90 profiler.f90 mpi_domain.F90 io.F90 utils.f90 conserve.f90 ionization.f90 eos.f90 matrix_vars.F90 opacity.f90 radiation_utils.f90 matrix_coeffs.f90 matrix_utils.f90 miccg.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 petsc_solver.F90 matrix_solver.F90 radiation.f90 cooling.f90 sinks.f90 externalforce.f90 composition.f90 star.f90 shockfind.f90 output.f90 smear.f90 input.f90 setup.f90 readbin.f90  restart.f90 fluxboundary.f90 rungekutta.f90 timestep.f90 tests.f90 gridset.f90 gravbound.f90 gravity_hyperbolic.f90 gravity_elliptic.f90 interpolation.f90 numflux.f90 source.f90
MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))

# Initial condition routines should be compiled before initialcondition.f90
INI = $(sort $(wildcard $(INI_DIR)/*.f90))

# Post-processing routines
POS = utils_analysis.f90
POSF = $(patsubst %,$(POS_DIR)/%,$(POS))

TARGETF = $(patsubst %,$(SRC_DIR)/%.f90,$(TARGET))
TARGETAF= $(patsubst %,$(POS_DIR)/%.f90,$(TARGETA))
TARGETLF= $(patsubst %,$(POS_DIR)/%.f90,$(TARGETL))
SRC = $(MODF) $(INI) \
      $(filter-out $(MODF) $(TARGETF), \
      $(sort $(wildcard $(SRC_DIR)/*.f90 $(SRC_DIR)/*.F90)))
SRCPOS = $(POSF) $(filter-out $(POSF), $(sort $(wildcard $(POS_DIR)/*.f90)))
OBJ  = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ := $(patsubst $(INI_DIR)/%.f90,%.o,$(OBJ))
OBJ := $(patsubst $(SRC_DIR)/%.F90,%.o,$(OBJ))
OBJ_BASE = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))
OBJ_POS  = $(filter-out $(TARGETAF) $(TARGETLF),$(SRCPOS))
OBJ_POS := $(patsubst $(POS_DIR)/%.f90,$(OBJ_DIR)/%.o,$(OBJ_POS))

OBJ_F  = $(OBJ_BASE) $(OBJ_DIR)/$(TARGET).o
OBJA_F = $(OBJ_BASE) $(OBJ_DIR)/$(TARGETA).o
OBJL_F = $(OBJ_BASE) $(OBJ_POS) $(OBJ_DIR)/$(TARGETL).o

all: $(TARGET)

debug: FCFLAG+=$(DBFLAG)
debug: $(TARGET)

profile: FCFLAG+=$(PRFLAG)
profile: $(TARGET)

ascii: $(TARGETA)

lc: $(TARGETL)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $@ $(OBJ_F) $(LDFLAGS)

$(TARGETA): $(OBJA_F)
	$(FC) $(FCFLAG) -o $@ $(OBJA_F)

$(TARGETL): $(OBJL_F)
	$(FC) $(FCFLAG) -o $@ $(OBJL_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	$(FC) $(FPPFLAG) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(INI_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

$(OBJ_DIR)/%.o: $(POS_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

.PHONY: clean
clean:
	rm -f *~ $(SRC_DIR)/*~ $(INI_DIR)/*~ $(POS_DIR)/*~ $(TARGET) $(TARGETA) $(TARGETL) $(OBJ_DIR)/*


install:
	@echo "Installing setup_hormone..."

# Detect user shell and rc file
	@read -p "Enter Shell environment [default: $$SHELL]: " SHELL_INPUT; \
	if [ -z "$$SHELL_INPUT" ]; then SHELL_INPUT="$$SHELL"; fi; \
	if [ -n "$$SHELL_INPUT" ]; then \
	  case "$$SHELL_INPUT" in \
	    *zsh) RCFILE="$$HOME/.zshrc"; SHELLTYPE=zsh;; \
	    *bash) RCFILE="$$HOME/.bashrc"; SHELLTYPE=bash;; \
	    *tcsh) RCFILE="$$HOME/.tcshrc"; SHELLTYPE=csh;; \
	    *csh) RCFILE="$$HOME/.cshrc"; SHELLTYPE=csh;; \
	    *) RCFILE="$$HOME/.profile"; SHELLTYPE=sh;; \
	  esac; \
	else \
	  RCFILE="$$HOME/.bashrc"; SHELLTYPE=bash; \
	fi; \
	echo "Using rc file: $$RCFILE ($$SHELLTYPE)"; \
	read -p "Enter Fortran compiler [default: gfortran]: " FC_INPUT; \
	FC=$${FC_INPUT:-gfortran}; \
	if [ -z "$$FC_INPUT" ]; then FC=gfortran; fi; \
	read -p "Enter number of OpenMP threads [default: 8]: " OMP_INPUT; \
	OMP=$${OMP_INPUT:-8}; \
	if [ -z "$$OMP_INPUT" ]; then OMP=8; fi; \
	echo "FC=$$FC, OMP_NUM_THREADS=$$OMP"; \
	$(SED_INPLACE) '/### hormone-setup begin ###/,/### hormone-setup end ###/d' $${RCFILE}; \
	sed -e "s/@FC@/$$FC/" \
	    -e "s/@OMP_NUM_THREADS@/$$OMP/" \
	    $(SCRIPTS_DIR)/hormone_setup | cat >> $${RCFILE}; \
	if [ "$$SHELLTYPE" = "csh" ]; then \
	  $(SED_INPLACE) '/## hormone-setup for bash ##/,/## end hormone-setup for bash ##/d' $${RCFILE}; \
	else \
	  $(SED_INPLACE) '/## hormone-setup for csh ##/,/## end hormone-setup for csh ##/d' $${RCFILE}; \
	fi; \
	if [ "$$SHELLTYPE" != "zsh" ]; then \
	  $(SED_INPLACE) '/## hormone-setup for zsh ##/,/## end hormone-setup for zsh ##/d' $${RCFILE}; \
	fi; \
	echo "Updated $$RCFILE"; \
	cp $(SCRIPTS_DIR)/setup ./setup; \
	chmod +x ./setup; \
	echo "Done. Restart your shell or run 'source $$RCFILE' to enable hormone-setup."

uninstall:
	@echo "Uninstalling setup_hormone..."
	@read -p "Enter Shell environment [default: $$SHELL]: " SHELL_INPUT; \
	if [ -z "$$SHELL_INPUT" ]; then SHELL_INPUT="$$SHELL"; fi; \
	if [ -n "$$SHELL_INPUT" ]; then \
	  case "$$SHELL_INPUT" in \
	    *zsh) RCFILE="$$HOME/.zshrc"; SHELLTYPE=zsh;; \
	    *bash) RCFILE="$$HOME/.bashrc"; SHELLTYPE=bash;; \
	    *tcsh) RCFILE="$$HOME/.tcshrc"; SHELLTYPE=csh;; \
	    *csh) RCFILE="$$HOME/.cshrc"; SHELLTYPE=csh;; \
	    *) RCFILE="$$HOME/.profile"; SHELLTYPE=sh;; \
	  esac; \
	else \
	  RCFILE="$$HOME/.bashrc"; SHELLTYPE=bash; \
	fi; \
	echo "Using rc file: $$RCFILE ($$SHELLTYPE)"; \
	$(SED_INPLACE) '/### hormone-setup begin ###/,/### hormone-setup end ###/d' $${RCFILE}
	@rm -f ./setup
	@echo "Uninstalled."
