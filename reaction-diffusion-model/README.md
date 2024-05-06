# Reaction-Diffusion-in-Biomolecular-Condensates
Reaction-diffusion numerical model for non-Fickian diffusion and sharp signaling fronts in all-DNA condensates

## System Requirements
The code was tested on\
  **Operating system:** Debian GNU/Linux10 (buster)\
  **SUNDIALS 2.5.0:** SUite of Nonlinear and DIfferential/ALgebraic equation Solvers\
  **Gnuplot 5.4.3**\
No special software or non-standard hardware is required. The above items are free and open source.

## Installation Guide
Download the SUNDIALS software package (https://computing.llnl.gov/projects/sundials/cvode).
Download the Gnuplot graphing program (http://www.gnuplot.info/).
Follow the installation guidelines of the publishers. Typical installation time is not more than some hours.

## Demo and Instructions for Use
  **[1]** Copy the folder **general-pc** to your computer.\
  **[2]** In the **Makefile**: change the file locations of the solver according to your folder structure.\
  **[3]** Use the **Terminal** to run the commands.\
\
  **3.1** Rebuild the executable by typing
  ```
  make
  ```
  This has to be done once, but has to be done again if the **general-pc.c** was updated.\
\
  **3.2** Run the code by typing
  ```
  ./run
  ```
  The **res** file contains the results in the following columns:
  time, space, (time- and space dependent) concentrations of A, B, C, and D, respectively.
  The units are in SI.\
\
\
  **[4]** Simulation parameters\
  Here you find a list of the parameters (with file locations in brackets) that can be modified to run parameter sweeps.\
\
  **4.1** Numerical parameters\
    - Total time: tend (**input_numpar**)\
    - Time stepping: N·delta (**input_numpar**)\
    - Spatial size: NX·dx, where NX (**general-pc.c**) is the mesh size. The number of lines in **input_osfr** has to match the mesh size.\
    - Spatial resolution: dx (**input_numpar**)\
    - Relative and absolute error tolerances: reltol and abstol (**input_numpar**)\
\
  **4.2** Chemical and kinetic parameters\
    - Rate constants: k1, k2 (**input_par**)\
    - Initial concentrations: Ai, Bi (**input_par**) matching the values in **input_osfr**\
    - Diffusion coefficients: D0, D1 (**input_par**)\
    - Threshold parameter: Q (**input_par**)
