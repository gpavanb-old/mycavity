ARMADILLO_DIR=/Users/gpavanb/cfd/numLibs/armadillo
ARMADILLO_LIB=$(ARMADILLO_DIR)/lib/
ARMADILLO_INC=$(ARMADILLO_DIR)/include/

SUPERLU_DIR=/Users/gpavanb/cfd/numLibs/SuperLU
SUPERLU_INC=$(SUPERLU_DIR)/SRC
SUPERLU_LIB=$(SUPERLU_DIR)/lib

OPT_FLAGS=-O3

.PHONY = clean clean_all

all: mycavity

mycavity: mycavity.cpp
	g++ mycavity.cpp -o mycavity $(OPT_FLAGS) -I$(ARMADILLO_INC) -I$(SUPERLU_INC) -L$(ARMADILLO_LIB) -larmadillo -L$(SUPERLU_LIB) -lsuperlu_mt_PTHREAD

clean:
	rm -f mycavity

clean_all:
	rm -f *.dat
	rm -f *.out
	rm -f mycavity 
