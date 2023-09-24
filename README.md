# dadaC
## Description
Code repository for the thesis project *Density-based clustering application to substructures identification in cosmological simulations*, Francesco Tomba, 2023 @ University of Trieste.

`dadaC` is the porting and optimization of the implementation of ADP (Laio et al. 2021) which is present in the python package [`dadaPy`](https://github.com/sissa-data-science/DADApy).
In particular dadaC implements at the moment:

- TWO-NN Intrinsic Dimension estimator
- k-NN search using a kd-tree
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

`./driver [input file] [output file] [z value] [halo (0 or 1)] [k] [sparse impl. flag (s or d)]`

- `input/output file`: file paths to import and export data
- `z`: Z value to use during the clustering procedure 
- `halo [0 or 1]`: compute clustering accounting for halo points
- `k`: number of neighbors to compute for each point
- `sparse impl. flag [s or d]`: use `s` to leverage sparse implementation of cluster borders (suggested for big datasets) `d` for dense implementation (to use for small dataset when cluster topology is an objective)


## Compiling

dadaC comes with a make file which compiles the executable `driver` and the shared library `bin/libclustering.so` from which ADP methods can be linked to.
dadaC supports compilation to use `float` or `double` to store data and kNN distances, and `uint32` or `uint64` to store indexes. 
Add `-DUSE_FLOAT32` or `-DUSE_INT32` to compile with support to 32bit types. By default 64bit ones are used. This feature is important for the application on big datasets, allowing of course some sort of rounding error to happen when processing data. 
Implementation with 64bit types results are binary equal w.r.t. `dadaPy`.


## Python interface

dadaC comes also with a python interface build with `ctypes` which leverages the capabilities of the C-compiled library. In order to use it, build the package and then import `dadaC` module from python

## TODO

MUCH MORE

- Complete porting of all density estimation methods in `dadaPy` (i.e. pAK)
- Integration of C compiled code into `dadaPy` library
- Improvement of memory imprinting for big datasets 
- Adaptive strategies on k-NN search and density estimation
- Add compilation targets to makefile to have different float type compilation


