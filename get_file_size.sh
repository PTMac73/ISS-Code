# Prints out the file sizes for the ISS raw data and calculates the largest set of files
# =============================================================================================== #
# Patrick MacGregor
# Nuclear Physics Research Group
# School of Physics and Astronomy
# The University of Manchester
# LAST EDITED: 24/10/18
# =============================================================================================== #
# GLOBAL VARIABLES
PREFIX="/home/ptmac/Documents/07-CERN-ISS-Mg/analysis/data/iss000_run_"
SUFFIX=".gtd*"

MAX_SIZE=0
MAX_INDEX=0

ARRAY=()
I_ARRAY=()
ARRAY_TOT=()

# Loop over the desired run numbers
for i in {20..33}
do
	FILECHECK="${PREFIX}${i}${SUFFIX}"

	# Extract the total
	TEST=$(du -c ${FILECHECK} &>> /dev/null )
	if [ $? -eq 0 ]
	then
		TAIL=$(du -c ${FILECHECK} | tail -1)


		# Strip the gubbins
		SIZE="${TAIL%total}"

		# Store in an array
		echo -e "$i \t $SIZE"

		# Find max size
		if [ $SIZE -gt $MAX_SIZE ]
		then
			MAX_SIZE=$SIZE
			MAX_INDEX=$i
		fi
		ARRAY+=($SIZE)
		I_ARRAY+=($i)
	fi
done

echo -e "\nMAX: $MAX_INDEX \t $MAX_SIZE"
