include $(HOME)/numLibs/elemental/conf/ElVars

CXX=CC
FLAGS=-g

all: test

test: mycavity.cpp
	$(CXX) $(FLAGS) $(EL_COMPILE_FLAGS) $< -o $@ $(EL_LINK_FLAGS) $(EL_LIBS)
