#!/bin/sh

if [ $# -eq 0 ]
then
    read -p 'Please enter the run number you would like to process: ' RUN
fi

RUN=$1

exp=iss631
expDir=/home/helios/experiments/${exp}/analysis

#ATTEMPTS=0

#while [ "${ATTEMPTS}" -lt "10000" ];
#do

#rsync -rtuh --delete --progress rsync://helios@anldaqrouter:12000/digiosdata /Users/calemhoffman/Research/isolde/data/.
# 
#not needed when sync daemon running
#/Users/heliosdigios/Applications/get_digios_data.sh $RUN $WHERE
#rsync -rtuavh --delete --progress helios@anldaqrouter:/media/DIGIOSDATA3/data/${exp}/*run_${RUN}* /Users/heliosdigios/experiments/${exp}/data/*

#/Users/heliosdigios/Applications/get_digios_data.sh $RUN 3
${expDir}/working/gebmerge_local.sh $RUN
${expDir}/working/gebsortmerged_local.sh $RUN

echo Just created root file run${RUN}.root in ${expDir}/root_data/
ls -ltrh ${expDir}/root_data

root -q -b "process_run.C(${RUN},0)"
cp gen.root ${expDir}/root_data/gen_run${RUN}.root
echo copied gen.root to gen_run${RUN}.root

echo ----Done with Processing Run Number ${RUN}----
#    echo "process event attempt number = ${ATTEMPTS}"

#    ATTEMPTS=$(expr ${ATTEMPTS} + 1)
#    sleep 5
    
#done

