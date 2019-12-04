#!/bin/bash
# doublet-fit-all-OMs fits all of the optical models using the doublet fitting code in root
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# =============================================================================================== #

# FUNCTIONS
CreateFileIfItDoesntExist(){
	if [ -e $1 ]
	then
		rm $1
	fi
	touch $1
}

CreateFileIfItDoesntExistNoDelete(){
	if [ ! -e $1 ]
	then
		touch $1
	fi
}



# GLOBAL VARIABLES
PT_DIR="/home/ptmac/Documents/07-CERN-ISS-Mg/Mg-Analysis/PtolemyOUT"
FIT_CODE_DIR="/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/analysis-codes/Doublet"
TEMPLATE_LOC="${FIT_CODE_DIR}/doublet_data_template.dat"
DATA_DIR="${FIT_CODE_DIR}/model-data"

SWITCH_MAKE_ARRAYS=1
SWITCH_CREATE_INPUT_FILES=1
SWITCH_RUN_ROOT=1
SWITCH_MAKE_TABLE=1

line_array=()
angle_array=()

# MAKE ARRAYS
if [[ $SWITCH_MAKE_ARRAYS -eq 1 ]]
then
	# Extract the lines from the template file and store them
	while read -r line
	do
		# Find the angle for a particular line
		ANGLE=${line%%	*}

		# Store the lines and angles in an array
		line_array+=("${line}")
		angle_array+=($ANGLE)
	done < $TEMPLATE_LOC
fi

# CREATE ALL OF THE INPUT DATA FILES
# Loop over all the directories
for i in "${PT_DIR}"/*/*-0.000.out-clean
do
	# Extract the model name
	FILE_NAME=${i##*/}
	L_STRIP=${FILE_NAME#*-}
	MOD_NAME=${L_STRIP%%-*}

	# Create the data file name
	DATA_FILE_NAME="${DATA_DIR}/doublet_data_${MOD_NAME}.dat"

	if [[ $SWITCH_CREATE_INPUT_FILES -eq 1 ]]
	then
		# Create the data file
		CreateFileIfItDoesntExist ${DATA_FILE_NAME}

		# Get the components of the L values
		index=0
		for j in "${angle_array[@]}"
		do
			# Add the first lines
			printf "%s" "${line_array[$index]}" >> ${DATA_FILE_NAME}
		
			# Round the angle
			ROUND_ANGLE=`printf "%.1f\n" $j`
			
			# Now search for that angle using regular expression (need to replace . with \. so that the
			# regexp works). Then select cut out the second field with " " as the delimiter
			PT_VALUES=($(grep -e "${ROUND_ANGLE//./\.}0 " $i | cut -f2- -d " "))

			# Put the angles in the data file as well [N.B. L = 0,1,2,2,3,4]
			for k in "${PT_VALUES[@]}"
			do
				printf "\t%s" "${k}" >> ${DATA_FILE_NAME}
			done

			# Add a new line character
			printf "%s\n" "" >> ${DATA_FILE_NAME}
			index=$((index+1))
		done
	fi

	
	# RUN THE DOUBLET DATA CODE
	# Create the output file name
	OUTPUT_FILE_NAME="${DATA_DIR}/doublet_output_${MOD_NAME}.dat"
	OUTPUT_FILE_NAME_TEMP="${OUTPUT_FILE_NAME//doub//temp_doub}"	

	if [[ $SWITCH_RUN_ROOT -eq 1 ]]
	then
		# Create the data file
		CreateFileIfItDoesntExist ${OUTPUT_FILE_NAME}
		CreateFileIfItDoesntExist ${OUTPUT_FILE_NAME_TEMP}

		root -l -q "${FIT_CODE_DIR}/DoubletFitter.C++(\"${DATA_FILE_NAME}\")" >> $OUTPUT_FILE_NAME_TEMP

		# Fix the output file by removing the first two lines (1,2d means delete lines 1 and 2)
		sed '1,2d' "${OUTPUT_FILE_NAME_TEMP}" > "${OUTPUT_FILE_NAME}"

		# Print output message to the console
		echo -e "${MOD_NAME} parameters calculated.\n"

	fi
done

# MAKE A MASSIVE TABLE OF RESULTS
# Delete all the CSV files
rm "${DATA_DIR}/"*.csv

if [[ $SWITCH_MAKE_TABLE -eq 1 ]]
then
	# Create CSV template name
	CSV_LOC_TEMP="${DATA_DIR}/doublet_output"

	# Organise model names
	for i in ${DATA_DIR}/doublet_o*.dat
	do
		# Get the name of the start of the model
		MOD_NAME=$( echo "${i##*/doublet_output_}" | cut -d "." -f 1 )
		MOD_NAME_START=$(echo "${MOD_NAME}" | cut -d "_" -f 1 )

		# Create CSV file if it doesn't exist
		CSV_FILE_NAME="${CSV_LOC_TEMP}-${MOD_NAME_START}.csv"
		CreateFileIfItDoesntExistNoDelete "${CSV_FILE_NAME}"

		# Append the data
		printf "%s,\u21132,\u03b71,\u03b72,\u03c7^2,\n" "${MOD_NAME}" >> "${CSV_FILE_NAME}"
		cat "${i}" | sed -r 's/\t/,/g' | sed -r 's/$/,/g' >> "${CSV_FILE_NAME}"
		printf ",,,,,\n" >> "${CSV_FILE_NAME}"
	done

	# Make a super-massive spreadsheet
	MEGA_CSV_LOC="${DATA_DIR}/mega_doublet_output.csv"
	CreateFileIfItDoesntExist "${MEGA_CSV_LOC}"
	FILE_LIST=""
	for i in ${DATA_DIR}/doublet_output-*.csv
	do
		FILE_LIST+="${i} "
	done

	paste -d ',' ${FILE_LIST} > ${MEGA_CSV_LOC}

fi













































