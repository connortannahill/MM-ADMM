This is the source code for the paper

MM-ADMM: Implicit Integration of MMPDEs in Parallel.
Connor Tannahill and Justin Wan.

This code generates an adapted mesh in areas indicated by a user-provided monitor function. The code requires as input some initial mesh (some initialization code for simple geometries is provided). A simple example of the code running on a circular geometry is given below. For more complex examples see the corresponding paper (to be submitted in 2022).

![Simple example of mesh adaptation in a circular domain](./MeshGif.gif)

Note that this code is still experimental and has not been optimized for any use besides our personal experiments. Documentation is sparse but feel free to contact the authors for any assistance.

And yes, we are embarrassed by the commit messages.

## Install

git clone https://github.com/connortannahill/OptMeshes.git  

## Dependencies

Install Eigen with the latest version at the standard library.  

git clone https://gitlab.com/libeigen/eigen.git  

And copy the /lib/Eigen/ subdirectory to the standard path:  

sudo cp -rf Eigen/ /usr/include/  

Also clone [nanoflann](https://github.com/jlblancoc/nanoflann.git) into the /lib/ directory.

## Compilation

Simply run make within the /OptMeshes directory to compile the project. If using the parallel capabilities you must add the -fopenmp flag in the CFLAGS variable and assign the DEFS variable to DEFS = -D THREADS

