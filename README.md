# dadaC
## Description
Code repository for the thesis project *Density-based clustering application to substructures identification in cosmological simulations*, Francesco Tomba, 2023 @ University of Trieste.

`dadaC` is the porting and optimization of the implementation of ADP (Laio et al. 2021) which is present in the python package `dadaPy`.
In particular dadaC implements at the moment:

    - k-NN search using a kd-tree
    - TWO-NN Intrinsic Dimension estimator
    - k*-NN density estimator
    - ADP Heuristics

## Usage

dadaC comes with an example driver file `driver.c`. Data is expected to be a matrix of floats of type `FLOAT_TYPE` (defined at compile time) of dimensions `N x data_dims` with `data_dims` being a global variable defined at compile time.

It returns a text file where for each point are reported:

    - `k*`: number of neighbors used for computing the density value 
    - `cluster_idx`: cluster assignement of the point
    - `rho`: value of the density
    - `is_center`: flag for identifying cluster centers

Once parameters are setted `dadaC` can be launched using:

`./driver [input file] [output file] [z value] [halo (0 or 1)] [k]`


