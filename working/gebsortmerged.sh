#!/bin/sh

if [ $# -ne 1 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

RUN=$1

dir=/Users/heliosdigios/experiments
exp=iss631

#cd $dir/GEBSort_working
GEBDIR=$dir/${exp}/analysis/GEBSort
MERGDIR=$dir/${exp}/anlysis/merged_data
ROOTDIR=$dir/${exp}/anlysis/root_data
CHATDIR=$dir/${exp}/analysis/working

echo "GEBSort started sorting run $RUN at `date`"
$GEBDIR/GEBSort_nogeb -input disk $MERGDIR/GEBMerged_run$RUN.gtd_000 -rootfile $ROOTDIR/run$RUN.root RECREATE -chat $CHATDIR/GEBSort.chat 
echo "GEBSort DONE at `date`"

#exit

