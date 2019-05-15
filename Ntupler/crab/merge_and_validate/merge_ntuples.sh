#!/bin/bash 

#set -e

nTarget=$1
srcName=$2
targName=$3
srcDir=$4
targDir=$5

files=($(eval xrdfs root://cmseos.fnal.gov ls -u ${srcDir} \| grep '${srcName}'))
nSource=${#files[@]}
#nSource=$(eos root://cmseos.fnal.gov ls -1 ${srcDir}/${srcName}_*.root | wc -l)
nHadd=$((nSource / nTarget))
remainder=$((nSource - (nTarget * nHadd)))

echo "Source directory contains ${nSource} \"${srcName}\" files"
if [ ${remainder} -ne 0 ]
then
    nHadd1=$nHadd
    nHadd2=$((nHadd + 1))
    echo "${nTarget} target \"${targName}\" files => ${nHadd1} or ${nHadd2} source \"${srcName}\" files each"
else
    nHadd1=$nHadd
    echo "${nTarget} target \"${targName}\" files => ${nHadd1} source \"${srcName}\" files each"
fi


#files=($(eos root://cmseos.fnal.gov ls -1 ${srcDir}/${srcName}_*.root))
startIndex=0



for ((i=0; i<nTarget; i+=1))
do
    echo ""

    # Pick nHadd to ensure most even file splitting
    if [ "$remainder" -ne 0 ]
    then
        nHadd=$nHadd2
        remainder=$((remainder - 1))
    else
        nHadd=$nHadd1
    fi


#   for filename in ${files[@]:startIndex:nHadd}
#   do
#       echo "Making local copy of ${filename}"
#       xrdcp root://cmseos.fnal.gov/${srcDir}/${filename} .
#   done
    hadd ${targName}_${i}.root ${files[@]:startIndex:nHadd}
    echo "Copying ${targName}_${i}.root to ${targDir}"
    if xrdcp ${targName}_${i}.root root://cmseos.fnal.gov/${targDir}/${targName}_${i}.root
    then
        echo "Deleting local copy of ${targName}_${i}.root"
        rm ${targName}_${i}.root

#       for filename in ${files[@]:startIndex:nHadd}
#       do
#           echo "Deleting local copy of ${filename}"
#           rm ${filename}
#       done
    fi


    startIndex=$((startIndex + nHadd))
done


if [ $startIndex == $nSource ]
then
    echo "Done!  Don't forget to check the output :)"
else
    echo "There was a problem :("
fi
