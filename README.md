# AtomMag
This repo contains the source code of a GPU-parallel Atomistic Spin Dynamics model AtomMag developed by [Prof.Jiamian Hu's group](https://mesomod.weebly.com/software.html) at University of Wisconsin-Madison. The simulation results using this code is validated through comparison with [Fidimag](https://fidimag.readthedocs.io/en/latest/ipynb/steepest_descent_atomistic.html) and analytical results. 66x speed up can be achieved on NVIDIA Tesla P100 compared to serial code. 
## Features
* Various Lattice type provided including simple cubic, face-centered cubic and 2D hexagonal lattice. Self-defined atomic position can also be input through 'atomposition.in'
* Heisenberg exchange interaction, Dzyaloshinskii–Moriya interaction, uniaxial magnetocrystalline energy, dipole-dipole interaction and zeeman energy are considered
* Spin orbit torque can be simulated
* Fortran code paralleled using OpenACC

## How to run the code

### 1. Install or load PGI compiler. For example, on Bridges
```
module load pgi
```

### 2. Compile
```
pgf90 -acc -Minfo=accel -o AtomMag AtomMag.f90
```

### 3.Run
```
./AtomMag
```
