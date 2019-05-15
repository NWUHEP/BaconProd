#!/bin/bash 

set -e

targName=$1
inputdir=$2
outputdir=$3

echo "Copying \"${targName}\" ntuples from input directory ${inputdir}"
echo "to output directory ${outputdir}"

nFiles=$(eos root://cmseos.fnal.gov ls -1 ${inputdir}/${targName}_*.root | wc -l)
echo "Directory contains ${nFiles} root files"

for filename in `eos root://cmseos.fnal.gov ls -1 ${inputdir}/${targName}_*.root`
do
    echo "Copying file ""$filename"
    xrdcp -f root://cmseos.fnal.gov/"$inputdir""/""$filename" root://cmseos.fnal.gov/"$outputdir""/""$filename"
done
