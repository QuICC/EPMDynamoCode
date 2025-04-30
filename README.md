## Overview
EPMDynamoCode is a simulation framework for convection and dynamos in a full sphere using a fully spetral method. For details of the nuemrical scheme, please refer to the paper
- Marti, P., & Jackson, A. (2016). A fully spectral methodology for magnetohydrodynamic calculations in a whole sphere. Journal of Computational Physics, 305, 403–422. https://doi.org/10.1016/J.JCP.2015.10.056
## Key Features
- Dynamo and precession dynamo simulations

- Parallel computation support

- Spectral methods implementation

- HDF5 state file output

- Visualization tools for ParaView/Visit

## Requirements
- Git

- GCC

- FFTW3

- LAPACK

- HDF5 (with parallel IO for MPI)

- CMake

## Installation
### Clone the repository:

```
git clone git@github.com:QuICC/EPMDynamoCode.git
cd EPMDynamoCode
git submodule init
git submodule update
```
### Create a build directory and configure with CMake:
```
mkdir build && cd build
ccmake ../EPMDynamoCode
```
Configure the build options as needed.
### Compile the desired simulation:
```
make DynamoSimulation -j4  # Example for dynamo simulation with 4 cores
```
## Running Simulations
### Basic Requirements
- parameters.cfg configuration file

- Initial state file (state_initial.hdf5)

- Simulation binary (e.g., DynamoSimulation)

### Example Commands
```
# For dynamo simulation
./DynamoSimulation

# For precession dynamo simulation
./PrecessionDynamoSimulation
```
## Documentation
- Technical documentation can be generated with make CodeDoc

- Configuration examples in Config/ directory

## Note
The code has been upgraded to the QuICC code https://github.com/QuICC/QuICC-Solver. This repository is not supported and updated anymore. 
