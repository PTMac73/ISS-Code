#!/bin/bash
# Sorts through the gen_run files and writes the desired quantities in the PTMonitors.C TSelector
# code to files
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# LAST EDITED: 26/10/18
# =============================================================================================== #
# GLOBAL VARIABLES
DIR=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/bufferSizeCheck/root_data

# SWITCHES
GEN_SWITCH=1
FIN_SWITCH=1

# PROCESS GEN_RUN FILES
if [ $GEN_SWITCH -eq 1 ]
then
	# Loop over files
	for f in $DIR/gen_run25-*.root
	do
		# Run PTMonitors code for the root file
		root -q -b -l "sort_runs.C(\"${f}\", 0)"

		# Rename the fin.root file
		mv fin.root "${DIR}/fin${f##*/gen}"
	done
fi

# PROCESS FIN_RUN FILES
if [ $FIN_SWITCH -eq 1 ]
then
	# Print initial table line
	echo -e "\nFILE NAME\tEVZ\tEXE"
	# Now get the entries in each histogram
	for g in $DIR/fin_run25-*.root
	do
		# Run the root script to get entries in the histograms
		root -q -b -l "sort_runs.C(\"${g}\", 1)"
	done
fi
