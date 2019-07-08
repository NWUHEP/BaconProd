BaconProd
=========

Branch for producting bacon ntuples for 2016-2018 data and MC. 

 * runs on CMS LPC
 * uses CMSSW\_10\_2\_13
 * depends on NWUHEP/BaconAna jbueghly\_prod branch

All objects are declared in BaconAna and filled in BaconProd/Ntupler, see e.g.:

[TJet for Jets](https://github.com/ksung25/BaconAna/blob/master/DataFormats/interface/TJet.hh)
[FillerJet for Jets](https://github.com/ksung25/BaconProd/blob/master/Ntupler/src/FillerJet.cc)

Setup
----------

```Shell
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git cms-merge-topic cmantill:baconprod-10213-v15
git clone -b jbueghly_prod git@github.com:NWUHEP/BaconProd
git clone -b jbueghly_prod git@github.com:NWUHEP/BaconAna
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210
git cms-addpkg EgammaAnalysis/ElectronTools
rm -rf EgammaAnalysis/ElectronTools/data
git clone git@github.com:cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd - 
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
scram b -j 12
```
Production
----------

To run, navigate to BaconProd/Ntupler/config. 

```
cmsRun makingBacon_25ns_MINIAOD.py isData=<isData> doHLTFilter=<doHLTFilter> era=<era>
```
 * isData: True for data, False for MC
 * doHLTFilter: if True, filter events that don't pass Bacon triggers
    + usually True for Data and False for MC
 * era: 2016, 2017, or 2018

Running with CRAB:
----------

Source crab:
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

Navigate to the BaconProd/Ntupler/crab directory. CRAB submission settings are controlled by the multicrab script. 

Specify the relevant splitting, configuration arguments, and datasets. Then run with:
```
./multicrab -c submit
``` 

and check status with 
```
./multicrab -c status -w <path to crab workArea>
```

Merging Output and Checking for Duplicates:
----------

Output ntuples can be merged by navigating to the BaconProd/Ntupler/crab/merge\_and\_validate directory and running:

```
./merge_ntuples.sh $nTarget $srcName $targName $srcDir $targDir
```

 * nTarget: final number of merged files
 * srcName: prefix in ntuple filenames: for example 'Output' or 'ntuple'
 * targName: prefix in output merged files
 * srcDir: directory containing files to be merged
 * targDir: directory for output merged ntuples

After merging, one may wish to validate that no events are duplicated in the output (from an errant merge). 
This can be done by running:

```
./find_doubles.sh $filePref $startIdx $inputdir
```

 * filePref: prefix for files to check: for example 'Output' or 'ntuple'
 * startIdx: initial number for filecount: for example 0 for Output\_0, Output\_1, etc.
 * inputdir: eos directory where files are located; NO trailing slash; e.g. '/store/user/.../MYDIR'

