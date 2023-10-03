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

dadaC comes with an example driver file `driver.c`. Data is expected to be a matrix of floats of type `FLOAT_TYPE` (defined at compile time) of dimensions `N x d`.

It returns a text file where for each point are reported:

- `k*`: number of neighbors used for computing the density value 
- `cluster_idx`: cluster assignement of the point
- `rho`: value of the density
- `is_center`: flag for identifying cluster centers

Once parameters are set `dadaC` can be launched using:
`./driver i=[INPUT_FILE] o=[OUTPUT_FILE] d=[d] t=[t] z=[Z] h=[HALO] k=[k] s=[s] t=[t]`

- `INPUT_FILE `: input file, file path
- `OUTPUT_FILE`: output file, file path
- `d`        : Lenght of the data vectors (number of columns of the data matrix)
- `Z`	     : Z value, float
- `HALO`     : Assign halo, bool 0/1 
- `k`	     : Number of neighbors to use, int (>0)
- `s`	     : Use sparse borders implementation, y [sparse/dense]
- `t`	     : Input binary is in Float32, y/n [float/double]



## Compiling

dadaC comes with a make file which compiles the executable `driver` and the shared library `bin/libdadac.so` from which ADP methods can be linked to.
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


