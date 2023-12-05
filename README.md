# cufQUMC
CUDA Fortran version of QU-fitting with replica exchange MCMC method.

## Prerequisites
* Python 3
* CUDA Fortran 23.1 or greater.
  * You can also use the NVHPC container: https://catalog.ngc.nvidia.com/orgs/nvidia/containers/nvhpc

## Installation
* Compile the Fortran code
<pre>
  $ nvfortran -cuda -cudalib=curand -gpu=cc80,managed -o libcufqumc libcufqumc.f90
</pre>
Here, "cc" means Compute Capability and "cc80" is for NVIDIA A100. See [this page](https://developer.nvidia.com/cuda-gpus) to check the compute capability of the GPU you use.

## Usage
* run Python code to perform QUMC
<pre>
  $ python cufqumc.py [par_file]
</pre>
The format of data files is [frequency (in Hz), Stokes Q, Stokes U, error in Q, error in U] (five columns).
<pre>
usage: cufqumc.py [-h] [-d DATADIR] [-m FITMODEL] -b NBURNIN -s NSAMPLE [-o OUTPUTDIR] [-n NBETA] [-r] [-t NTHREADS]

CUDA Fortran QUMC.

optional arguments:
  -h, --help    show this help message and exit
  -d DATADIR    Directory containing data files.
  -m FITMODEL   Fitting Model. ["delta", "gauss", "tophat", "osul17"] are available.
                e.g. One "delta" and two "gauss" components: -m="{'delta':1,'gauss':2}".
  -b NBURNIN    Number of chain for Burn-in (step widths adjustment) [Required].
  -s NSAMPLE    Number of chain for sampling [Required].
  -o OUTPUTDIR  Output directory ["./out"].
  -n NBETA      The number of Replica [100]. beta=[1] when Nbeta=1.
  -r            Resume MCMC from the last results in "OutputDir" directory [False]. 
                When activated, only "Nburnin", "Nsample", "OutputDir" and "Nthreads" are needed to specify.
  -t NTHREADS   Number of threads per block [1024].
</pre>

## Fitting Models
* delta
* gauss
* tophat
* osul17

## Notes
This code runs QUMC and outputs parameter chains at beta=1.0 and chisq chains.
