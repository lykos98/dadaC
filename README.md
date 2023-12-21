# dadaC
## NOTE
This implementation is in progress, currently clenaning up code and working on documentation. 

## Description
Code repository for the thesis project *Density-based clustering application to substructures identification in cosmological simulations*, Francesco Tomba, 2023 @ University of Trieste.

`dadaC` is the porting and optimization of the implementation of ADP (Laio et al. 2021) which is present in the python package [`dadaPy`](https://github.com/sissa-data-science/DADApy).
In particular dadaC implements at the moment:

- TWO-NN Intrinsic Dimension estimator
- k-NN search using a kd-tree (or vp-tree with custom callable metric)
- k*-NN density estimator
- ADP Heuristics

On the same input dadaC achieves a one to one match on the results, obtaining up to a factor 40 speedup on the whole procedure w.r.t. Python/Cython implementation. 

## Benchmarks

Validation of the procedure is done against the implementation of the original package with and without the computation of the halo. Hereafter some of the results, complete ones in the directory benchmarks.

On AMD Ryzen 7 7735HS @ 4.80GHz (my laptop) (8 cores - 16 Threads, 16GB RAM)

_from dadapy examples_ 

**Fig1.dat:  N = 20.0k, D = 2**


| Method   | part             | time   |
|----------|------------------|--------|
| py       | ngbh and density | 1.53s  |
| py       | ADP              | 0.13s  |
| C        | ngbh and density | 0.26s  |
| C        | ADP              | 0.04s  |

_from dadapy examples_ 

**Fig2.dat:  N = 38.4k, D = 2**

| Method   | part             | time   |
|----------|------------------|--------|
| py       | ngbh and density | 2.76s  |
| py       | ADP              | 0.94s  |
| C        | ngbh and density | 0.22s  |
| C        | ADP              | 0.03s  |

**CosmoSim (sub)Set1:   N = 100.0k,   D = 5 (~500 Clusters)**

| Method   | part             | time   |
|----------|------------------|--------|
| py       | ngbh and density | 21.93s |
| py       | ADP              | 4.08s  |
| C        | ngbh and density | 1.37s  |
| C        | ADP              | 0.07s  |

On Intel Xeon Gold 5118 CPU @ 2.30GHz (4 sockets x 12 cores - 48 Threads, 512GB RAM)

**CosmoSim Set1:   N = 1.8M,     D = 5 (~2000 Clusters)**

| Method   | part             | time     |
|----------|------------------|----------|
| py       | ngbh and density | 414.83s  |
| py       | ADP              | 3282.80s |
| C        | ngbh and density | 47.52s   |
| C        | ADP              | 7.46s    |

**MNIST N = 70.0k,  D = 784**

| Method   | part             | time   |
|----------|------------------|--------|
| py       | ngbh and density | 11.62s |
| py       | ADP              | 2.21s  |
| C        | ngbh and density | 12.07s |
| C        | ADP              | 0.18s  |



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
- `HALO`     : Assign halo, y/n [yes/no] 
- `k`	     : Number of neighbors to use, int (>0)
- `s`	     : Use sparse borders implementation, y/n [sparse/dense]
- `t`	     : Input binary is in Float32, y/n [float/double]

Relies on kd-Tree or vantage point tree alogirthms in order to perform neighborhood search. _V2_ versions are optimized ones, use them. KNN search is validated against `scipy.spatial.KDtree` query time is comparable with state of the art libraries except for bruteforce search used when the dimensionality D > 15. 

## Python interface

dadaC comes also with a python interface build with `ctypes` which leverages the capabilities of the C-compiled library. In order to use it, build the package and then import `dadaC` module from python

## Compiling

dadaC comes with a make file which compiles the executable `driver` and the shared library `bin/libdadac.so` from which ADP methods can be linked to.
dadaC supports compilation to use `float` or `double` to store data and kNN distances, and `uint32` or `uint64` to store indexes. 
Add `-DUSE_FLOAT32` or `-DUSE_INT32` to compile with support to 32bit types. By default 64bit ones are used. This feature is important for the application on big datasets, allowing of course some sort of rounding error to happen when processing data. 
Implementation with 64bit types results are binary equal w.r.t. `dadaPy`.



**REQUIRES** Numpy package to properly work

## TODO

MUCH MORE

- Complete porting of all density estimation methods in `dadaPy` (i.e. pAK)
- Adaptive strategies on k-NN search and density estimation


