# Readme for oxDNA data analysis with gyration.py



## System Requirements

Running the script requires an interpreter for Python 3.x, and having NumPy and matplotlib installed.

It was developed and tested under Ubuntu 22.04 using conda 3.

No non-standard hardware is required.



## Installation guide

If Python is installed (including NumPy and matplotlib), no further installation procedure is required.

The provided Python code will directly be executed by the interpreter.



## Demo

Copy the script to the same directory as demo.top and demo.dat. Open a terminal, enter the directory, and execute
```
python3 gyration.py demo.top demo.dat 1
```

The expected output is:
```
looking for strand with index 1
trajectory file will be evaluated from the beginning
trajectory file will be evaluated until the end

demo.top -> evaluating only nucleotides in strand 1
demo.top -> number of evaluated nucleotides: 1025
demo.top -> index  of reference strand for end-to-end distance: 1
demo.top -> length of reference strand for end-to-end distance [nt]: 1025
demo.top -> base composition of evaluated nucleotides:
	 650 A
	 175 T
	 50 C
	 150 G
	 0 unknown
	 -----------------------------------------------------------
	 total mass of evaluated nucleotides [g/mol]: 320662.0
	 mean  mass of evaluated nucleotides [g/mol]: 312.8409756097561

demo.dat -> found 258 configurations
demo.dat -> evaluating 258 configurations

estimated CPU time to evaluate trajectory [min]: 0.05724310929999976
required  CPU time to evaluate trajectory [min]: 0.05492885955

root mean square radius of gyration for homopolymer [nm]:
	 51.053384049424785
root mean square radius of gyration for heteropolymer [nm]:
	 51.05356899038488
root mean square end-to-end distance [nm]:
	 148.22835959105706
r.m.s. radius of gyration for homopolymer inferred from r.m.s. end-to-end distance [nm]:
	 60.5139744013118
gyr2.png <- trajectory plotted
```
The output file gyr2.png should be found in the same directory and show a diagram of square radii of gyration as calculated for a homopolymer and a heteropolymer over simulation steps as well as a sixth of the square end-to-end distance, the latter quantity showing stronger fluctuations. As can be easily seen from the plot, the molecule is not at equilibrium, but undergoing a contraction process.

Expected run time for this demo is around five seconds on a "normal" desktop computer.



## Instructions for use

The easiest way to run the analysis is to copy the script to the same directory as the raw data and execute
```
	python gyration.py xxx.top yyy.dat X Y Z
```
with
```
	python	an interpreter for Python 3.x (usually python3 or python if using conda),
	xxx.top	the name of the topology file (can be a path as well),
	yyy.dat	the name of the configuration or trajectory file (can be a path as well),
	X	the index of the strand for which radii of gyration are calculated (defaults to 0):
		by setting this parameter to an existing strand index, only nucleotides of this strand will be evaluated
		(this avoids artifacts introduced by strand-wise shifts from oxDNA diffusion correction in long simulations),
		otherwise, all nucleotides will be treated as if they are part of one strand/duplex,
	Y	number of configurations to skip at the beginning,
	Z	number of configurations to skip at the end.
```
Ensemble-averaged results will be printed to standard output.
Additionally, a trajectory over configurations will be plotted as gyr2.png, so you can check for ongoing relaxation processes.
