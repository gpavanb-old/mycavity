HOME_DIR=/home/pavang

ARMADILLO_DIR=$(HOME_DIR)/numLibs/armadillo
ARMADILLO_LIB=$(ARMADILLO_DIR)/lib/
ARMADILLO_INC=$(ARMADILLO_DIR)/include/

SUPERLU_DIR=$(HOME_DIR)/numLibs/SuperLU
SUPERLU_INC=$(SUPERLU_DIR)/include
SUPERLU_LIB=$(SUPERLU_DIR)/lib

OPT_FLAGS=-O3

.PHONY = clean clean_all

all: mycavity

mycavity: mycavity.cpp
	g++ mycavity.cpp -o mycavity $(OPT_FLAGS) -I$(ARMADILLO_INC) -I$(SUPERLU_INC) -L$(ARMADILLO_LIB) -larmadillo

clean:
	rm -f mycavity

clean_all:
	rm -f *.dat
	rm -f *.out
	rm -f mycavity 
