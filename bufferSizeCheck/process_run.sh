#!/bin/bash
# Generates root files for varying buffer sizes (use on run 51)
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# LAST EDITED: 25/10/18
# =============================================================================================== #
exp=iss000
dir=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis

GEBDIR=$dir/GEBSort
MERGDIR=$dir/bufferSizeCheck/merged_data
ROOTDIR=$dir/bufferSizeCheck/root_data
DATADIR=$dir/data
#MERGECHAT=$dir/bufferSizeCheck/GEBMerge.chat
SORTCHAT=$dir/bufferSizeCheck/GEBSort.chat

RUN=51

# Loop over the .chat files
for f in $dir/bufferSizeCheck/chat_files/*.chat
do
	# Define the merge chat file
	MERGECHAT=$f
	
	# Extract the bigbufsize number
	BBNUM=${MERGECHAT##*GEBMerge-}
	BBNUM=${BBNUM%.chat}

	# Run GEBMerge
	echo "RUN $RUN: GEBMerge started at `date`"
	$GEBDIR/GEBMerge $MERGECHAT  $MERGDIR/GEBMerged_run$RUN-$BBNUM.gtd `ls $DATADIR/${exp}_run_$RUN.gtd*` > $MERGDIR/GEBMerge_run$RUN-$BBNUM.log
	echo "RUN $RUN: GEBMerge DONE at `date`"

	echo "GEBSort started sorting run $RUN at `date`"
	$GEBDIR/GEBSort_nogeb -input disk $MERGDIR/GEBMerged_run$RUN-$BBNUM.gtd_000 -rootfile $ROOTDIR/run$RUN-$BBNUM.root RECREATE -chat $SORTCHAT 
	echo "GEBSort DONE at `date`"
	echo "saved root file -->  "  $ROOTDIR/run$RUN-$BBNUM.root 

	echo "=================="
	root -q -b "process_run.C(${RUN},${BBNUM})"

	cp gen.root ./root_data/gen_run$RUN-$BBNUM.root

	echo "============================================="
	echo " done, $ROOTDIR/gen_run${RUN}-${BBNUM}.root"
	echo "============================================="

done
