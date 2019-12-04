#!/bin/bash
STEP=1000
FILE="testfile.dat"
if [ -e $FILE ]
then
	rm $FILE
fi

touch $FILE

EX=4.31613
Z_START=-29.43
STEP_SIZE=0.00001

for (( i=0; i<$STEP; i++ ))
do
	printf "${EX}\t%14.12f\n" "$( echo "${Z_START} + ${STEP_SIZE}*$i" | bc )"
done >> $FILE
