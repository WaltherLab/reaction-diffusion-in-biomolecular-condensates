# Readme for MD simulations with oxDNA



## System Requirements

Reproducing the MD simulations with the oxDNA model from the input files provided in the subdirectories requires a local instance of the oxDNA software (v3.6.0) with CUDA support in order to use an NVIDIA GPU.

Alternatively (e.g. if no machine with an NVIDIA GPU is available locally), the service provided via oxDNA.org can be used to run the simulations online, cf. https://oxdna.org/



## Installation guide

cf. https://lorenzo-rovigatti.github.io/oxDNA/install.html

Installation time strongly depends on necessary troubleshooting. It is advised to work under Linux and ensure compatibility of GCC and NVCC versions, cf. https://gist.github.com/ax3l/9489132 and https://github.com/NVlabs/instant-ngp/issues/119#issuecomment-1034701258 for possible problems during installation.



## Demo

The standalone version of oxDNA comes with a variety of example simulations. Because of this, no additional demo is provided in this case.



## Instructions for use

After successful installation of oxDNA, MD simulations -- each one starting from a relaxed configuration -- can be done by executing
```
oxDNA input
```
from within the ssDNA subdirectory, and
```
oxDNA input2
```
from within the dsDNA subdirectory, respectively. Each one will take about three days, depending on your computer system.

The expected output are files containing particle trajectories and energy values. The latter one can be used to ensure the system is at equilibrium. The trajectories can be evaluated using the self-made data analysis tool gyration.py.
