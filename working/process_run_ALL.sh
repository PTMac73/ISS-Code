#!/bin/bash
# Sorts all of the data files (takes a long time!)
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# LAST EDITED: 01/11/18
# =============================================================================================== #
SCRIPT_DIR=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/working
START=10
STOP=125
LOG_FILE=$SCRIPT_DIR/process_run_ALL.log
TOTAL_TIME=0

if [ -e $LOG_FILE ]
then
	echo "Log file exists. Clearing contents."
	> $LOG_FILE
else
	echo "No log file detected. Making a new file."
	touch $LOG_FILE
fi

# Write start time to the log file
echo "`date +"%Y-%m-%d %T"`    SORTING RUNS BEGINS" >> $LOG_FILE

# Loop over all the runs
for (( i=$START; i<=$STOP; i++ ))
do
	# Start the timer and reset time
	SECONDS=0
	DURATION=0
	
	# Print start message
	if [ $i -lt 100 ]
	then
		echo "`date +"%Y-%m-%d %T"`    Run 0${i} start" >> $LOG_FILE
	else
		echo "`date +"%Y-%m-%d %T"`    Run ${i} start" >> $LOG_FILE
	fi
	
	# Run the actual script
	$SCRIPT_DIR/process_run_simple.sh $i
	
	# Stop the timer
	DURATION=$SECONDS

	# Add the duration to the total time
	TOTAL_TIME=$(($TOTAL_TIME + $DURATION))

	# Print friendly message to console
	echo "Run $i finished after ${DURATION}s"
	echo "Total time elapsed: ${TOTAL_TIME}s"
	
	# Print finish message to log
	if [ $i -lt 100 ]
	then
		echo "`date +"%Y-%m-%d %T"`    Run 0${i} finish (${DURATION}s)" >> $LOG_FILE
	else
		echo "`date +"%Y-%m-%d %T"`    Run ${i} finish (${DURATION}s)" >> $LOG_FILE
	fi
done

# Write finishing time
echo "`date +"%Y-%m-%d %T"`    SORTING RUNS FINISHES (${TOTAL_TIME}s)" >> $LOG_FILE
