#!/bin/bash
# Sorts through the gen_run files and writes the desired quantities in the PTMonitors.C TSelector
# code to files
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# LAST EDITED: 24/10/18
# =============================================================================================== #
# GLOBAL VARIABLES
DIR=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/bufferSizeCheck/root_data


# Loop over files
for f in $DIR/gen_run*450.root
do
	# Run PTMonitors code for the root file
	root -q -b "sort_runs.C(${f##/*/})"

	# Rename the fin.root file
	mv fin.root "fin${f##*/gen}"
done
