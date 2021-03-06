# Mycavity #

This C++ code is a 2D, structured, massively parallel CFD code for a thermal lid-driven cavity problem. This code uses the following packages

* Armadillo[1] for the sparse matrix support routines and SuperLU connection
* SuperLU[2] library for the sequential sparse-direct solver
* Elemental[3] for the multi-CPU parallelization 

For more information, please refer to the [report](gpavanb_eboigne_cfd_report.pdf) or [presentation](gpavanb_me469_presentation.pptx) in the repository.

## How do I get set up? ##

### Sequential ###
* Install [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/), specifically v5.2.0 from the specified site. This is a CMake based installation

* Install [Armadillo](http://arma.sourceforge.net/). Note that v7.800 was tested for this code and it specifically used the above SuperLU version. It is recommended that the same version of Armadillo also be used.

* Change `ARMADILLO_DIR` in the `Makefile`. The code can now be compiled using `make` and run using `./mycavity`.

### Parallel - CPU ###

* Install [Elemental](http://libelemental.org/) and its list of dependencies. Note that v0.87 was used for this code and many sparse matrix routines are on the verge of deprecation.

* Run `make` in the `mycavity` directory. The code `mycavity` can now be run using the batch management system of your cluster.

* Additional options include `--verbose true` and `--timing true`. These can be used to obtain stepwise prompts for the code or measure time respectively.

### Parallel - GPU ###

* Download [Paralution](http://www.paralution.com) and its list of dependencies. Note that v1.1.0 was used for this code.

* Before compilling the source code, the following changes are necessary:
  - Make the private attribute `initial_residual_` of the IterationControl class public. To do so, in `/paralution-1.1.0/src/solvers/iter_ctrl.hpp`, move the line `double initial_residual_;` from the private section to the public one.
  - Make the protected attribute `iter_ctrl_` of the Solver class public. To do so, in `/paralution-1.1.0/src/solvers/solver.hpp`, move the line `IterationControl iter_ctrl_;` from the private section to the public one.
  
* Compile the code as detailled in the Paralution instructions.

* Run `make` in the `mycavity` directory.

## Who do I talk to? ##

* Repo owner or admin : [Pavan Bharadwaj](https://github.com/gpavanb)
* Other community or team contact : The code was developed at Stanford University. Please direct any official queries to [Prof. Gianluca Iaccarino](mailto:jops@stanford.edu)

## References ##

[1] Sanderson, Conrad. Armadillo: An open source C++ linear algebra library for fast prototyping and computationally intensive experiments. Technical report, NICTA, 2010.

[2] Demmel, James W., John R. Gilbert, and Xiaoye S. Li. SuperLU user's guide. Computer Science Division, University of California, 1997.

[3] Poulson, Jack, et al. "Elemental: A new framework for distributed memory dense matrix computations." ACM Transactions on Mathematical Software (TOMS) 39.2 (2013): 13.

[4] PARALUTION Labs, Paralutionv1.1.0, [http://www.paralution.com/], 2016.
