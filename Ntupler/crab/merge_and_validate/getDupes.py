from __future__ import print_function
import numpy as np
import sys
import time


inFilePref = sys.argv[1]

inFileName = inFilePref + "_EvtNum.txt"
outFileName = inFilePref + "_dupes.txt"


# Load array of (evtNum, file index)
print("Loading", inFileName, "into numpy array...", end='') 
sys.stdout.flush()
stTime = time.clock()
#allEvts = np.loadtxt(inFileName, dtype=[('entryIdx', np.uint64), ('fileIdx', np.uint16), ('evtNum', np.uint64)])
allEvts = np.loadtxt(inFileName, dtype=[('entryIdx', np.uint64), ('fileIdx', np.uint16), ('runNum', np.uint64), ('lumiSec', np.uint64), ('evtNum', np.uint64)])
elTime = time.clock() - stTime
print(elTime, "seconds")


# Get array of duplicates 
print("Processing numpy.unique...", end='')
sys.stdout.flush()
stTime = time.clock()
#_, uniqIdx = np.unique(allEvts['evtNum'], return_index=True)
#print(allEvts[['runNum', 'lumiSec', 'evtNum']][:5])
_, uniqIdx = np.unique(allEvts[['runNum', 'lumiSec', 'evtNum']], return_index=True)
elTime = time.clock() - stTime
print(elTime, "seconds")
dupEvts = np.delete(allEvts, uniqIdx)


# Print 
print("Found", len(dupEvts), "duplicated events") 
np.savetxt(outFileName, dupEvts, fmt='%i', delimiter='\t')
print("Wrote duplicates to", outFileName)
