#!/bin/sh

if [ $# -ne 1 ] 
  then
   echo "Specify only the run number to sort"
  exit 1
fi

exp=iss000
dir=/home/ptmac/Documents/07-CERN-ISS-Mg/analysis

GEBDIR=$dir/GEBSort
MERGDIR=$dir/merged_data
ROOTDIR=$dir/root_data
DATADIR=$dir/data
MERGECHAT=$dir/working/GEBMerge.chat
SORTCHAT=$dir/working/GEBSort.chat


RUN=$1

#for RUN in {223..237}
#do
echo "============================================="
echo "============ RUN $RUN ======================="
echo "============================================="

du -hSc $DATADIR/${exp}_run_$RUN.gtd*

echo "RUN $RUN: GEBMerge started at `date`"
$GEBDIR/GEBMerge $MERGECHAT  $MERGDIR/GEBMerged_run$RUN.gtd `ls $DATADIR/${exp}_run_$RUN.gtd*` > $MERGDIR/GEBMerge_run$RUN.log
echo "RUN $RUN: GEBMerge DONE at `date`"

echo "GEBSort started sorting run $RUN at `date`"
$GEBDIR/GEBSort_nogeb -input disk $MERGDIR/GEBMerged_run$RUN.gtd_000 -rootfile $ROOTDIR/run$RUN.root RECREATE -chat $SORTCHAT 
echo "GEBSort DONE at `date`"

echo "saved root file -->  "  $ROOTDIR/run$RUN.root 

echo "=================="
root -q -b "process_run.C(${RUN},0)"

cp gen.root ../root_data/gen_run$RUN.root

echo "============================================="
echo " done, $ROOTDIR/gen_run${RUN}.root"
echo "============================================="
#done

#exit
