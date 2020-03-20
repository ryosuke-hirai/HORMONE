FC  = ifort #gfortran
#-J for gfortran, -module for ifort
MODFLAG = -module 
#FCFLAG = -fopenmp -ffree-line-length-512 #-fmax-stack-var-size=32768 # for gfortran
FCFLAG =  -openmp #-heap-arrays #-traceback -fpe0 -CB #-pg -check all # for ifort

OBJ_DIR = obj
SRC_DIR = src

TARGET = hormone

# compile modules first
MODF= modules.f90
SRC = $(SRC_DIR)/$(MODF) $(filter-out $(SRC_DIR)/$(MODF), \
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
