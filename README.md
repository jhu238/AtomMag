# AtomMag
This repo contains the source code of a GPU-parallel Atomistic Spin Dynamics model AtomMag developed by [Prof.Jiamin Hu's group](https://mesomod.weebly.com/software.html) at University of Wisconsin Madison. The simulation results using this code is validated through comparison with [Fidimag](https://fidimag.readthedocs.io/en/latest/ipynb/steepest_descent_atomistic.html) and analytical results. 66x speed up can be achieved on NVIDIA Tesla P100 compared to serial code. 

## How to run the code

### 1. Install or load PGI compiler. For example, on Bridges
```
module load pgi
```

### 2. Compile the code 
```
pgf90 -acc -Minfo=accel -o AtomMag AtomMag.f90
```

### 3.Run the code
```
./AtomMag
```
