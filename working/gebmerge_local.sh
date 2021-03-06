#!/bin/sh

if [ $# -ne 1 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

RUN=$1

exp=iss000
dir=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis

echo "RUN $RUN: GEBMerge started at `date`"
 
$dir/GEBSort/GEBMerge $dir/working/GEBMerge.chat  $dir/merged_data/GEBMerged_run$RUN.gtd `ls $dir/data/${exp}_run_$RUN.gtd*`
 
echo "RUN $RUN: GEBMerge DONE at `date`"
