'''
This script can be used to calculate the radius of gyration of all nucleotides or a specific strand in an oxDNA configuration or trajectory file.
Radii of gyration will be calculated both in a sequence-averaged and a sequence-dependent way,
thereby treating the analyzed nucleotides as monomeric units in a homopolymer or a heteropolymer, respectively.

For a homopolymer, root mean square radius of gyration is calculated as
	R_G = sqrt(< sum_i (r_i - r)^2 / N >)
with
	R_G	radius of gyration,
	N	chain length in nucleotides,
	r_i	3d position of nucleotide i,
	r	3d mean position of all nucleotides:
		r = sum_i r_i / N,
	sum_i	sum over i = 1, ..., N,
	sqrt	square root,
	< >	ensemble average.
In this case, mean square radius of gyration < R_G^2 > and mean square end-to-end distance < R^2 > are related via
	< R_G^2 > = < R^2 > / 6
in the limit of infinite chain length for the unperturbed chain (Cantor and Schimmel: Biophysical Chemistry, 1980, 983),
thus, the root mean square end-to-end distance is calculated, too.
If all nucleotides are set to be considered part of a duplex (as is the default), R will be determined for the longest strand (the reference strand).

For a heteropolymer, root mean square radius of gyration is calculated as
	R_G = sqrt(< sum_i m_i*s_i^2 / M >)
with
	m_i	mass of nucleotide i,
	s_i	distance of nucleotide i from the common center of mass
		s_i = r_i - (sum_j m_j*r_j / M),
	M	total mass of all nucleotides
		M = sum_j m_j,
	and everything else as before.

The easiest way to run the analysis is to copy the script to the same directory as the raw data and execute
	python gyration.py xxx.top yyy.dat X Y Z
with
	python	an interpreter for Python 3.x (usually python3 or python if using conda),
	xxx.top	the name of the topology file (can be a path as well),
	yyy.dat	the name of the configuration or trajectory file (can be a path as well),
	X	the index of the strand for which radii of gyration are calculated (defaults to 0):
		by setting this parameter to an existing strand index, only nucleotides of this strand will be evaluated
		(this avoids artifacts introduced by strand-wise shifts from oxDNA diffusion correction in long simulations),
		otherwise, all nucleotides will be treated as if they are part of one strand/duplex,
	Y	number of configurations to skip at the beginning,
	Z	number of configurations to skip at the end.

Ensemble-averaged results will be printed to standard output.
Additionally, a trajectory over configurations will be plotted as gyr2.png, so you can check for ongoing relaxation processes.
Please note that in the case of Monte Carlo simulations, the ordinate does not correspond to physical time.

This code is cc-by-sa. It is provided in a ready-to-use form without any warranty.
If you use it for your research, please give appropriate credit by citing the paper in which it was published:
	Weixiang Chen, Brigitta Dúzs, Pablo G. Argudo, Sebastian V. Bauer, Wei Liu, Avik Samanta, Sapun H. Parekh, Mischa Bonn, Andreas Walther:
	Ultrasharp and ballistic diffusion fronts in biomolecular condensates (2024)

-- Sebastian V. Bauer, University of Mainz, 18.03.2024
-- sebastian.bauer@uni-mainz.de
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
import time

###############################################
# customized functions to evaluate oxDNA data #
###############################################

# determine center of mass for a single nucleotide from an oxDNA configuration file
def get_oxDNA_cm_site(line):
	line_arr = (np.asarray(line.split())).astype(float)
	return line_arr[0:3]

# determine time of a configuration from an oxDNA trajectory file
def get_oxDNA_step(line):
	line_arr = np.asarray(line.split())
	return int(line_arr[2])

# determine mass of nucleotide (values from Tseytlin: Advanced Mechanical Models of DNA Elasticity, 2016, 32)
def get_mass_for_nucleotide(base):
	if base == "A":
		return 313.21
	elif base == "T":
		return 304.2
	elif base == "C":
		return 289.18
	elif base == "G":
		return 329.21
	else:
		return 325.0

# determine base composition from a sequence array
def get_composition(sequence):
	bases = np.zeros(5, int)
	for base in sequence:
		if base == "A":
			bases[0] += 1
		elif base == "T":
			bases[1] += 1
		elif base == "C":
			bases[2] += 1
		elif base == "G":
			bases[3] += 1
		else:
			bases[4] += 1
	return bases

# determine absolute mass from composition
def get_mass_from_composition(bases):
	m_A = bases[0] * 313.21
	m_T = bases[1] * 304.2
	m_C = bases[2] * 289.18
	m_G = bases[3] * 329.21
	m_X = bases[4] * 325.0
	return m_A + m_T + m_C + m_G + m_X

# convert length from oxDNA simulation units to nanometer
def su2nm(x):
	return 0.8518 * x

# convert time from oxDNA simulation units to picoseconds
def su2ps(t):
	return 3.031 * t



###################
# data evaluation #
###################

# parse positional arguments
if len(sys.argv) < 3:
	print("python gyration.py topology.top trajectory.dat X Y Z")
	print("You can choose only to evaluate the strand with index X (defaults to 0, meaning all nucleotides will be evaluated).")
	print("You can choose whether to skip the first Y and last Z configurations by setting the last positional arguments (defaults to 0 0).")
	sys.exit(2)
## mandatory arguments
topfile_path  = sys.argv[1]
trajfile_path = sys.argv[2]
## optional arguments
try:
	strand = int(sys.argv[3])
except:
	strand = 0
	print("no valid argument for a strand index received")
else:
	print("looking for strand with index", strand)
try:
	skip_first = int(sys.argv[4])
except:
	skip_first = 0
	print("trajectory file will be evaluated from the beginning")
else:
	print("the first", skip_first, "configurations in the trajectory file will be skipped")
try:
	skip_last = int(sys.argv[5])
except:
	skip_last = 0
	print("trajectory file will be evaluated until the end")
else:
	print("the last", skip_last, "configurations in the trajectory file will be skipped")	
print()



# evaluate topology file
## extract data
topfile = open(topfile_path, 'r')
line = topfile.readline()
line_lst = line.split()
number_of_nucleotides = int(line_lst[0])
number_of_strands = int(line_lst[1])
strand_lengths = np.zeros(number_of_strands + 1, int)
strand_lengths[0] = number_of_nucleotides
sequence = []
mass = np.zeros(number_of_nucleotides, float)
for n in range(number_of_nucleotides):
	line = topfile.readline()
	strand_index = int(line.split()[0])
	strand_lengths[strand_index] += 1
	base = line.split()[1]
	sequence.append(base)
	mass[n] = get_mass_for_nucleotide(base)
topfile.close()
## check if the strand specified per optional argument exists
if strand < 0 or strand > number_of_strands:
	print(topfile_path, "-> strand index out of range, defaulting to evaluation of all nucleotides")
	strand = 0
elif strand == 0:
	print(topfile_path, "-> evaluating all nucleotides as part of one molecule")
else:
	print(topfile_path, "-> evaluating only nucleotides in strand", strand)
## determine nucleotide indices for evaluation
eva_nt = strand_lengths[strand]
offset = 0
for i in range(1, strand):
	offset += strand_lengths[i]
## get base sequence and masses
composition = get_composition(sequence[offset:(offset+eva_nt)])
total_mass = get_mass_from_composition(composition)
mean_nucleotide_mass = total_mass / eva_nt
relative_mass = mass / total_mass
## manage indices of nucleotides within reference strand for evaluation of end-to-end distance
e2erefstrand_index = strand
if strand == 0:
	e2erefstrand_index = 1 + np.argmax(strand_lengths[1:])
e2erefstrand_length = strand_lengths[e2erefstrand_index]
e2erefstrand_offset = 0
if strand == 0:
	for i in range(1, e2erefstrand_index):
		e2erefstrand_offset += strand_lengths[i]
## tell user
print(topfile_path, "-> number of evaluated nucleotides:", eva_nt)
print(topfile_path, "-> index  of reference strand for end-to-end distance:", e2erefstrand_index)
print(topfile_path, "-> length of reference strand for end-to-end distance [nt]:", e2erefstrand_length)
print(topfile_path, "-> base composition of evaluated nucleotides:")
print("\t", composition[0], "A")
print("\t", composition[1], "T")
print("\t", composition[2], "C")
print("\t", composition[3], "G")
print("\t", composition[4], "unknown")
print("\t", "-----------------------------------------------------------")
print("\t", "total mass of evaluated nucleotides [g/mol]:", total_mass)
print("\t", "mean  mass of evaluated nucleotides [g/mol]:", mean_nucleotide_mass)
print()



# browse whole trajectory file to get number of configurations
all_configs = 0
trajfile = open(trajfile_path, 'r')
line = trajfile.readline()
while line:
	for n in range (number_of_nucleotides + 2):
		line = trajfile.readline()
	all_configs += 1
	line = trajfile.readline()
trajfile.close()
num_configs = all_configs - (skip_first + skip_last)
print(trajfile_path, "-> found", all_configs, "configurations")
if num_configs < 1:
	print(trajfile_path, "-> no configurations left after skipping, aborting.")
	sys.exit(1)
else:
	print(trajfile_path, "-> evaluating", num_configs, "configurations\n")



# evaluate configurations in trajectory file
pos           = np.zeros((eva_nt, 3), float)
averaged_gyr2 = np.zeros(num_configs, float)
weighted_gyr2 = np.zeros(num_configs, float)
e2e2          = np.zeros(num_configs, float)
steps         = np.zeros(num_configs, int)
time_cpu      = np.zeros(4, float)
trajfile = open(trajfile_path, 'r')
time_cpu[0] = time.process_time()
for k in range(skip_first):
	## skip unwanted configurations
	for l in range(number_of_nucleotides + 3):
		line = trajfile.readline()
for k in range(num_configs):
	## measure CPU time to evaluate first configuration
	if k == 0:
		time_cpu[1] = time.process_time()
	
	## get sim time at which the current configuration was recorded
	line = trajfile.readline()
	steps[k] = get_oxDNA_step(line)

	## skip box size and energy
	for l in range(2):
		line = trajfile.readline()
	
	## skip nucleotides not to be evaluated
	for n in range(offset):
		line = trajfile.readline()
	
	## get nucleotide positions
	for n in range(eva_nt):
		line = trajfile.readline()
		pos[n] = get_oxDNA_cm_site(line)
	
	## calculate common center of mass for the evaluated nucleotides (for homopolymer and heteropolymer, respectively)
	averaged_com = np.mean(pos, axis=0)
	weighted_positions = mass[offset:(offset+eva_nt), None] * pos
	weighted_com = np.sum(weighted_positions, axis=0) / total_mass
	
	## calculate square radii of gyration (for homopolymer and heteropolymer, respectively)
	averaged_difference = pos - averaged_com
	averaged_gyr2[k] = np.mean(np.einsum('ij,ij->i', averaged_difference, averaged_difference))
	weighted_difference = pos - weighted_com
	weighted_gyr2[k] = np.sum(np.einsum('ij,ij->i', relative_mass[offset:(offset+eva_nt), None] * weighted_difference, weighted_difference))
	
	## calculate square end-to-end distance
	e2e = pos[e2erefstrand_offset + e2erefstrand_length - 1] - pos[e2erefstrand_offset]
	e2e2[k] = np.dot(e2e, e2e)
	
	## scroll down to next configuration
	for l in range (number_of_nucleotides - (offset + eva_nt)):
		line = trajfile.readline()
	
	## get estimate for CPU time to evaluate whole trajectory
	if k == 0:
		time_cpu[2] = time.process_time()
		print("estimated CPU time to evaluate trajectory [min]:", num_configs * (time_cpu[2] - time_cpu[1]) / 60)
	
	## don't evaluate remaining configurations
trajfile.close()
time_cpu[3] = time.process_time()
print("required  CPU time to evaluate trajectory [min]:", (time_cpu[3] - time_cpu[0]) / 60)
print()



###############
# get results #
###############

# calculate ensemble means
averaged_gyr = np.sqrt(np.mean(averaged_gyr2))
weighted_gyr = np.sqrt(np.mean(weighted_gyr2))
e2e = np.sqrt(np.mean(e2e2))
gyr_from_e2e = np.sqrt(np.mean(e2e2)/6)

# output results
## numerical output
print("root mean square radius of gyration for homopolymer [nm]:\n\t", su2nm(averaged_gyr))
print("root mean square radius of gyration for heteropolymer [nm]:\n\t", su2nm(weighted_gyr))
print("root mean square end-to-end distance [nm]:\n\t", su2nm(e2e))
print("r.m.s. radius of gyration for homopolymer inferred from r.m.s. end-to-end distance [nm]:\n\t", su2nm(gyr_from_e2e))
## plot
fig1, ax1 = plt.subplots()
ax1.set_xlabel("simulation steps")
ax1.set_ylabel("square radius of gyration [nm²]")
ax1.plot(steps, su2nm(averaged_gyr2), color="blue", alpha = 0.75, label="homopolymer (averaged)")
ax1.plot(steps, su2nm(weighted_gyr2), color="red",  alpha = 0.75, label="heteropolymer (weighted)")
ax1.plot(steps, su2nm(e2e2/6),        color="grey", alpha = 0.50, label="R² / 6")
ax1.legend()
plt.grid()
fig1.savefig(f"gyr2.png")
plt.cla()
plt.clf()
print("gyr2.png <- trajectory plotted")
print()



