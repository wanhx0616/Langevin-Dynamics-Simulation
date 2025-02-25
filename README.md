### Langevin Dynamics Simulation on CPUs for Cell Entry of Soft Active Particles ###

# Overview

This project implements a Langevin Dynamics simulation for cell entry of soft active particles. 
The simulation is designed to study the nonequilibrium interaction and dynamics in endocytosis process. 
The program is written in Fortran and is optimized for CPU execution.

The algorithms used in the program is based on the following publications:
Vutukuri, H. R. et. al, Active particles induce large shape deformations in giant lipid vesicles. 
Nature 2020, 586, 52.
Yuan H. et. al, One-particle-thick, solvent-free, coarse-grained model for biological and biomimetic fluid membranes.
Phys. Rev. E 2010, 82, 011905.
Gnan N. et. al, The microscopic role of deformation in the dynamics of soft colloids.
Nat. Phys. 2019 15, 683.

# System Requirements
## Hardware requirements

This package requires a standard computer with a multi-core CPU to support efficient parallel computation.

## Software requirements
### OS Requirements
This package is supported for *Linux* and *Windows* operating systems. The package has been tested on the following systems:

Linux: Ubuntu >16.04
Windows: Windows 10

### Compiler Requirements
GNU Fortran (gfortran)

# Installation Guide:

This directory contains the Fortran source code for the Langevin Dynamics simulation. 
The main program is in the file bd.f90, while parameters are in the file module.f90.

To compile the code, run:

gfortran bd.f90 module.f90 -o3 -o bd

Once the code is built, you can execute it using:

./bd

# System Description

There are four types of beads in the system:
Type 1: beads that make up the membrane
Type 2: beads that make up the soft particle
Type 3: the bead on the surface of the soft particle according to the given initial orientation
Type 4: the center of mass of the soft particle

# File Description:

'aout.xyz': orientations of membrane beads
'auvw': evolution velocity and acceleration of the orientation of membrane beads
'data': used for storing data
'out.xyz': positions of membrane beads
'outnp.xyz': positions of soft particle beads
'PARA': bonds information of soft particle beads
'uvw': evolution velocity and acceleration of the position of membrane beads
'uvwnp': evolution velocity and acceleration of soft particle beads

# Parameter Adjustment：

Number of cell membrane bead：#define memnum 10501
Number of soft active particle bead：#define nanonum 676
Diameter of soft active particle：#define diameter1 16.0

# DEMO

- To run demo:
  - `gfortran bd.f90 module.f90 -o3 -o bd`
  - `./bd`
  then the output can be found in the current directory.

# License

This project is covered under the **MIT License**.
