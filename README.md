# dadaC
## Description
Code repository for the thesis project *Density-based clustering application to substructures identification in cosmological simulations*, Francesco Tomba, 2023 @ University of Trieste.

`dadaC` is the porting and optimization of the implementation of ADP (Laio et al. 2021) which is present in the python package [`dadaPy`](https://github.com/sissa-data-science/DADApy).
In particular dadaC implements at the moment:

- k-NN search using a kd-tree
- TWO-NN Intrinsic Dimension estimator
- k*-NN density estimator
- ADP Heuristics

On the same input dadaC achieves a one to one match on the results, obtaining a factor 40 speedup on the whole procedure w.r.t. Python/Cython implementation. 
NOTE: this is an experimental implementation, program may crash if a dataset requires a memory usage larger than the memory available on the machine. This issue will be fixed soon. 

## Usage

dadaC comes with an example driver file `driver.c`. Data is expected to be a matrix of floats of type `FLOAT_TYPE` (defined at compile time) of dimensions `N x data_dims` with `data_dims` being a global variable defined at compile time.

It returns a text file where for each point are reported:

- `k*`: number of neighbors used for computing the density value 
- `cluster_idx`: cluster assignement of the point
- `rho`: value of the density
- `is_center`: flag for identifying cluster centers

Once parameters are setted `dadaC` can be launched using:

`./driver [input file] [output file] [z value] [halo (0 or 1)] [k]`

## Compiling

dadaC comes with a make file which compiles the executable `driver` and the shared library `bin/libclustering.so` from which ADP methods can be linked to.

## TODO

MUCH MORE

- Complete porting of all density estimation methods in `dadaPy` (i.e. pAK)
- Integration of C compiled code into `dadaPy` library
- Improvement of memory imprinting for big datasets 
- Adaptive strategies on k-NN search and density estimation
- Add compilation targets to makefile to have different float type compilation


