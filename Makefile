ARCH=-arch sm_30

COMP=nvcc $(ARCH)

PARALUTION_DIR=~/Documents/paralution-1.1.0
PARALUTION_LIB=$(PARALUTION_DIR)/build/lib/
PARALUTION_INC=$(PARALUTION_DIR)/build/inc/

OPT=-O3

mycavityGPU: mycavityGPU.o
	$(COMP) $(OPT) -L $(PARALUTION_LIB)  -o x.mycavityGPU mycavityGPU.o -lparalution

mycavityGPU.o: mycavityGPU.cpp
	$(COMP) -I $(PARALUTION_INC) -c mycavityGPU.cpp

clean:
	rm -rf mycavityGPU.o x.mycavityGPU mycavityGPU.out mycavityCPU.out mycavity_tecplot_gpu_0.dat mycavity_tecplot_cpu_0.data