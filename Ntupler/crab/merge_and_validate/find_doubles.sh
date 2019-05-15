#!/bin/bash 

set -e
# You should probably do cmsenv in the correct CMSSW area before trying to execute this

filePref=$1     # Prefix of ROOT files NOT including underscore
                #   (e.g. "Output" for files called "Output_*.root")

startIdx=$2     # Number index of first file; files must be consecutive!!!
                #   (e.g. "1" if we start with "Output_1.root") 

inputdir=$3     # EOS directory where files are located; NO trailing slash
                #   (e.g. "/store/user/.../MYDIR")


# Get the number of "${filePref}_*.root" files in the directory
nFiles=$(eos root://cmseos.fnal.gov ls -1 ${inputdir}/${filePref}_*.root | wc -l)


# Call ROOT macro to loop over all of the files
#   Writes out Info->evtNum and the file's number index to "${filePref}_evtNum.txt"
#   Also prints number of events in each file and total over all files
root.exe -q -b "printEvtNum.cc(\"${filePref}\", ${startIdx}, ${nFiles}, \"${inputdir}\")"


# Call a numpy script to search for duplicates across all files
#   Reads in "${filePref}_evtNum.txt"
#   Writes out only the lines for duplicated events to "${filePref}_dupes.txt"
#   Also prints total number of duplicates found
python getDupes.py ${filePref}
