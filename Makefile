#FC  = ifort
#MODFLAG = -module 
#CFLAG =  -qopenmp #-heap-arrays -traceback -fpe0 -CB #-pg -check all # for ifort

FC = gfortran
MODFLAG = -J
FCFLAG = -fopenmp -ffree-line-length-512 -ffpe-trap=invalid,zero,overflow -fbacktrace #-fmax-stack-var-size=32768 # for gfortran

OBJ_DIR = obj
SRC_DIR = src

TARGET = hormone

# compile modules first
MOD = modules.f90 recombination.f90 eos.f90 fluxlimiter.f90 hlldflux.f90 eigen.f90
MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))
SRC = $(MODF) $(filter-out $(MODF), $(sort $(wildcard $(SRC_DIR)/*.f90)))
#SRC = $(SRC_DIR)/$(MODF) $(filter-out $(SRC_DIR)/$(MODF), \
      $(sort $(wildcard $(SRC_DIR)/*.f90)))
OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ_F = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))

all: $(TARGET)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $@ $(OBJ_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

clean:
	rm -f obj/*.o *~ src/*~ obj/*.mod hormone
