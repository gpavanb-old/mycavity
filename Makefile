include $(HOME)/numLibs/elemental/conf/ElVars

CXX=CC
FLAGS=-g

.PHONY = clean clean_all

all: mycavity

mycavity: mycavity.cpp
	$(CXX) $(FLAGS) $(EL_COMPILE_FLAGS) $< -o $@ $(EL_LINK_FLAGS) $(EL_LIBS)

clean:
	rm -f mycavity

clean_all:
	rm -f *.dat
	rm -f *.out
	rm -f mycavity 
