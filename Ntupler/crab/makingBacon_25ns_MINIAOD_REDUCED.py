import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as vp
import os

process = cms.Process('MakingBacon')

options = vp.VarParsing('analysis')

# register command line options
options.register('isData',
                 True, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('doHLTFilter',
                 True, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('doAlpaca',
                 False, # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )

options.register('era',
                 '2016', # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.string, 
                 )

options.register('isLocal',
                 'False', # default value
                 vp.VarParsing.multiplicity.singleton,
                 vp.VarParsing.varType.bool, 
                 )
# user input values
options.parseArguments()

hlt_filename  = "BaconAna/DataFormats/data/HLTFile_25ns"   # list of relevant triggers

cmssw_base = os.environ['CMSSW_BASE']
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if options.era == '2016':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('94X_dataRun2_v10')
      JECTag='Summer16_23Sep2016AllV4_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('94X_mcRun2_asymptotic_v3')
      JECTag='Summer16_23Sep2016V4_MC'

elif options.era == '2017':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('94X_dataRun2_v11')
      JECTag='Fall17_17Nov2017_V32_102X_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v17')
      JECTag='Fall17_17Nov2017_V32_102X_MC'

elif options.era == '2018':
    if options.isData:
      process.GlobalTag.globaltag = cms.string('102X_dataRun2_Sep2018ABC_v2')
      JECTag='Fall17_17Nov2017_V32_102X_DATA'
    else:
      process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v18')
      JECTag='Fall17_17Nov2017_V32_102X_MC'

from BaconProd.Ntupler.myJecFromDB_cff    import setupJEC
setupJEC(process,options.isData,JECTag)
if options.isData:
	if options.isLocal:
		process.jec.connect = cms.string('sqlite_file:/uscms/home/corderom/nobackup/2016/CMSSW_10_2_13/src/BaconProd/Utils/data/'+JECTag+'.db')
	else:
		process.jec.connect = cms.string('sqlite:///src/BaconProd/Utils/data/'+JECTag+'.db')
else:
  	if options.isLocal:
		process.jec.connect = cms.string('sqlite_file:/uscms/home/corderom/nobackup/2016/CMSSW_10_2_13/src/BaconProd/Utils/data/'+JECTag+'.db')
	else:
  		process.jec.connect = cms.string('sqlite:///src/BaconProd/Utils/data/'+JECTag+'.db')
#--------------------------------------------------------------------------------
# Import of standard configurations
#================================================================================
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

if options.era == '2017' or options.era == '2018':
    process.load('RecoBTag.SoftLepton.SoftLeptonByMVAComputers_cff')

process.pfNoPileUpJME = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
process.load('BaconProd/Ntupler/myPUPPICorrections_cff')
process.load('BaconProd/Ntupler/myCHSCorrections_cff')
process.load('BaconProd/Ntupler/myCorrections_cff')

#--------------------------------------------------------------------------------
# Import custom configurations
#================================================================================
# custom jet stuff (incl. GenJets, b-tagging, grooming, njettiness)
process.load('BaconProd/Ntupler/myGenJets_cff')
process.load('BaconProd/Ntupler/myJetExtrasAK4CHS_cff')

process.load('BaconProd/Ntupler/myJetExtrasAK4Puppi_cff')

process.load("RecoBTag.ImpactParameter.impactParameter_cff")
process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
process.load("RecoBTag.SoftLepton.softLepton_cff")
process.load("RecoBTag.Combined.combinedMVA_cff")
process.load("RecoBTag.CTagging.cTagging_cff")
process.load("RecoBTag.Combined.deepFlavour_cff")
if options.era == '2016':
    process.pfDeepCSVJetTags.NNConfig = cms.FileInPath('BaconProd/Utils/data/DeepFlavourNoSL.json')

from BaconProd.Ntupler.myBtagging_cff           import addBTagging,addBTaggingAK4CHS
from BaconProd.Ntupler.myGenJets_cff            import setMiniAODGenJets
from BaconProd.Ntupler.myJetExtrasAK4CHS_cff    import setMiniAODAK4CHS

from BaconProd.Ntupler.myJetExtrasAK4Puppi_cff  import setMiniAODAK4Puppi

process.btagging = cms.Sequence()
addBTagging(process,'AK4PFJetsPuppi' ,0.4,'AK4' ,'Puppi')
addBTaggingAK4CHS(process,'updatedPatJets'    ,0.4,'AK4' ,'CHS',True,True)

setMiniAODGenJets(process)
setMiniAODAK4CHS(process)

setMiniAODAK4Puppi (process)

#METFilters
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#CHS
process.chs = cms.EDFilter("CandPtrSelector",
                           src = cms.InputTag('packedPFCandidates'),
                           cut = cms.string('fromPV')
                           )
if options.isData:
  process.AK4QGTaggerCHS.jec  = cms.InputTag("ak4chsL1FastL2L3ResidualCorrector")
  process.AK4QGTaggerSubJetsCHS.jec  = cms.InputTag("ak4chsL1FastL2L3ResidualCorrector")

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if options.era == '2016':
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2016-Legacy', #era is new to select between 2016 / 2017,  it defaults to 2017
                           phoIDModules=[]) #bug with default modules for photon VID; off for now
elif options.era == '2017':
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2017-Nov17ReReco', #era is new to select between 2016 / 2017,  it defaults to 2017
                           phoIDModules=[]) #bug with default modules for photon VID; off for now

elif options.era == '2018':
    setupEgammaPostRecoSeq(process,
                           runVID=True,
                           era='2018-Prompt', #era is new to select between 2016 / 2017,  it defaults to 2017
                           phoIDModules=[]) #bug with default modules for photon VID; off for now

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
if options.era == '2016':
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
        DataEra = cms.string("2016BtoH"), 
        UseJetEMPt = cms.bool(False),
        PrefiringRateSystematicUncty = cms.double(0.2),
        SkipWarnings = False)
elif options.era == '2017' or options.era=='2018':
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
        DataEra = cms.string("2017BtoF"), 
        UseJetEMPt = cms.bool(False),
        PrefiringRateSystematicUncty = cms.double(0.2),
        SkipWarnings = False)

# PF MET corrections
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=options.isData,
                           manualJetConfig=True,
                           jetCorLabelL3="ak4chsL1FastL2L3Corrector",
                           jetCorLabelRes="ak4chsL1FastL2L3ResidualCorrector",
                           reclusterJets=True,
                           recoMetFromPFCs=True,
                           postfix="V2"
                           )

# PUPPI Woof Woof
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD (process, True )
process.puppi.useExistingWeights = True
runMetCorAndUncFromMiniAOD(process,
                           isData=options.isData,
                           manualJetConfig=True,
                           metType="Puppi",
                           pfCandColl=cms.InputTag("puppiForMET"),
                           recoMetFromPFCs=True,
                           jetFlavor="AK4PFPuppi",
                           jetCorLabelL3="ak4PuppiL1FastL2L3Corrector",
                           jetCorLabelRes="ak4PuppiL1FastL2L3ResidualCorrector",
                           reclusterJets=True,
                           postfix="Puppi"
                           )

if options.isData:
  process.AK4QGTaggerPuppi.jec           = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")
  process.AK4QGTaggerSubJetsPuppi.jec    = cms.InputTag("ak4PuppiL1FastL2L3ResidualCorrector")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
  btagDiscriminators = ['pfDeepCMVAJetTags:probb','pfDeepCMVAJetTags:probc','pfDeepCMVAJetTags:probudsg','pfDeepCMVAJetTags:probbb',
                        'pfDeepCSVJetTags:probb' ,'pfDeepCSVJetTags:probc' ,'pfDeepCSVJetTags:probudsg' ,'pfDeepCSVJetTags:probbb']
  )

# ALPACA
alpacaMet = ''
alpacaPuppiMet = ''
if options.doAlpaca: 
  alpacaMet      = ('pfMetAlpacaData'        if options.isData else 'pfMetAlpacaMC' )
  alpacaPuppiMet = ('pfMetPuppiAlpacaData'   if options.isData else 'pfMetPuppiAlpacaMC' ) 

#--------------------------------------------------------------------------------
# input settings
#================================================================================
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

if options.era == '2016':
    if options.isData:
        #test_file = cms.untracked.vstring('/store/data/Run2016C/DoubleEG/MINIAOD/03Feb2017-v1/810000/D8A12591-5DED-E611-BAF8-02163E019C24.root')
        test_file = cms.untracked.vstring('/store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver1-v1/20000/F6F9090E-C395-E811-8593-0242AC1C0503.root')
    else:
        test_file = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/GluGluHToZG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/702BC677-2AC8-E611-A5B3-02163E011463.root')
elif options.era == '2017':
    if options.isData:
        #test_file = cms.untracked.vstring('/store/data/Run2017F/DoubleEG/MINIAOD/31Mar2018-v1/90001/F61EA338-8E37-E811-A203-0025905C4262.root')
        test_file = cms.untracked.vstring('/store/data/Run2017B/DoubleMuon/MINIAOD/31Mar2018-v1/00000/06893E33-AB37-E811-A32E-0019B9CAB9CB.root')
    else:
        #test_file = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/46FB5EDE-F708-E811-A50F-0025905C53A4.root')
        test_file = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/005DC030-D3F4-E711-889A-02163E01A62D.root')
elif options.era == '2018':
    if options.isData:
        test_file = cms.untracked.vstring('/store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/00000/0DFD591F-DB0C-F447-9FBF-BAC1BF96D9B1.root')
    else:
        test_file = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/90000/F6F754E3-9026-CC48-9018-FFBB087DADA5.root')

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/GluGluHToZG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/702BC677-2AC8-E611-A5B3-02163E011463.root')
                            fileNames = test_file
                            )
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

#--------------------------------------------------------------------------------
# Reporting
#================================================================================
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE'),
)

#--------------------------------------------------------------------------------
# Bacon making settings
#================================================================================
process.ntupler = cms.EDAnalyzer('NtuplerMod_REDUCED',
  skipOnHLTFail     = cms.untracked.bool(options.doHLTFilter),
  useTrigger        = cms.untracked.bool(True),
  #useTriggerObject  = cms.untracked.bool(False),
  useTriggerObject  = cms.untracked.bool(True),

  TriggerObject     = cms.untracked.string("selectedPatTrigger"),
  TriggerFile       = cms.untracked.string(hlt_filename),
  useAOD            = cms.untracked.bool(False),
  outputName        = cms.untracked.string('Output.root'),
  edmPVName         = cms.untracked.string('offlineSlimmedPrimaryVertices'),
  edmGenRunInfoName = cms.untracked.string('generator'),
  
  Info = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    edmPFCandName        = cms.untracked.string('packedPFCandidates'),
    edmPileupInfoName    = cms.untracked.string('slimmedAddPileupInfo'),
    edmBeamspotName      = cms.untracked.string('offlineBeamSpot'),
    edmMETName           = cms.untracked.string('slimmedMETs'),
    edmPFMETName         = cms.untracked.InputTag('slimmedMETsV2'),
    edmMVAMETName        = cms.untracked.string(''),
    edmPuppETName        = cms.untracked.InputTag('slimmedMETsPuppi'),
    edmAlpacaMETName     = cms.untracked.string(alpacaMet),
    edmPupAlpacaMETName  = cms.untracked.string(alpacaPuppiMet),
    edmRhoForIsoName     = cms.untracked.string('fixedGridRhoFastjetAll'),
    edmRhoForJetEnergy   = cms.untracked.string('fixedGridRhoFastjetAll'),
    doFillMETFilters     = cms.untracked.bool(True),
    doFillMET            = cms.untracked.bool(True)
  ),
  
  GenInfo = cms.untracked.PSet(
    isActive                = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    edmGenEventInfoName     = cms.untracked.string('generator'),
    edmGenParticlesName     = cms.untracked.string('prunedGenParticles'),
    edmGenPackParticlesName = cms.untracked.string('packedGenParticles'),
    fillAllGen              = cms.untracked.bool(True),
    fillLHEWeights          = cms.untracked.bool(True)
  ),

  GenJet  = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False)),
    isActiveFatJet      = ( cms.untracked.bool(False)),
    #isActive            = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    #isActiveFatJet      = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    edmGenParticlesName = cms.untracked.string('prunedGenParticles'),
    genJetName          = cms.untracked.string('AK4GenJetsCHS'),
    fillAllGen          = cms.untracked.bool(False)
  ),
                                   
  PV = cms.untracked.PSet(
    isActive      = cms.untracked.bool(True),   
    edmName       = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    minNTracksFit = cms.untracked.uint32(0),
    minNdof       = cms.untracked.double(4),
    maxAbsZ       = cms.untracked.double(24),
    maxRho        = cms.untracked.double(2)
  ),
  
  Electron = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(7),
    edmName                   = cms.untracked.string('slimmedElectrons'),
    edmSCName                 = cms.untracked.InputTag('reducedEgamma','reducedSuperClusters'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppiNoLep'),
    usePuppi                  = cms.untracked.bool(True),
    #useTriggerObject          = cms.untracked.bool(False),
    useTriggerObject          = cms.untracked.bool(True),

    edmEcalPFClusterIsoMapTag = cms.untracked.InputTag('electronEcalPFClusterIsolationProducer'),
    edmHcalPFClusterIsoMapTag = cms.untracked.InputTag('electronHcalPFClusterIsolationProducer'),

    edmEleMVASpring16         = cms.untracked.string('ElectronMVAEstimatorRun2Spring16GeneralPurposeV1'),
    edmEleMVAFall17V1Iso      = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV1'),
    edmEleMVAFall17V1NoIso    = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV1'),
    edmEleMVAFall17V2Iso      = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV2'),
    edmEleMVAFall17V2NoIso    = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV2'),
    edmEleMVASpring16HZZ      = cms.untracked.string('ElectronMVAEstimatorRun2Spring16HZZV1'),

    edmEleMediumMVA           = cms.untracked.string('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    edmEleTightMVA            = cms.untracked.string('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),

    #edmEleMVA                 = cms.untracked.string('ElectronMVAEstimatorRun2Spring16GeneralPurposeV1'),
    #edmEleMediumMVAIso        = cms.untracked.string(''), 
    #edmEleTightMVAIso         = cms.untracked.string(''),
    #edmEleMVAIso              = cms.untracked.string(''),
    #storeSecondMVA            = cms.untracked.bool(False),
    #storeHZZMVA               = cms.untracked.bool(True),
  ),
  
  Muon = cms.untracked.PSet(
    isActive                  = cms.untracked.bool(True),
    minPt                     = cms.untracked.double(3),
    edmName                   = cms.untracked.string('slimmedMuons'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppiNoLep'),
    usePuppi                  = cms.untracked.bool(True),
    #useTriggerObject          = cms.untracked.bool(False),    
    useTriggerObject          = cms.untracked.bool(True),
  ),
  
  Photon = cms.untracked.PSet(
    isActive              = cms.untracked.bool(True),
    minPt                 = cms.untracked.double(10),
    edmName               = cms.untracked.string('slimmedPhotons'),
    edmSCName             = cms.untracked.InputTag('reducedEgamma','reducedSuperClusters'),
    edmPhoMVASpring16     = cms.untracked.string('PhotonMVAEstimatorRun2Spring16NonTrigV1'),
    edmPhoMVAFall17V1     = cms.untracked.string('PhotonMVAEstimatorRunIIFall17v1'),
    edmPhoMVAFall17V2     = cms.untracked.string('PhotonMVAEstimatorRunIIFall17v2'),

    eeReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits"), 
    ebReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits"), 
    esReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedESRecHits"),

    useTriggerObject      = cms.untracked.bool(False),
  ),
  
  
  AK4CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(False),
    useTriggerObject     = cms.untracked.bool(False),
    minPt                = cms.untracked.double(15),
    maxEta               = cms.untracked.double(3),
    coneSize             = cms.untracked.double(0.4),
    addPFCand            = cms.untracked.bool(False),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak4chsL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak4chsL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK4chs')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),

    # names of various jet-related collections
    jetName              = cms.untracked.string('selectedUpdatedPatJets'),
    genJetName           = cms.untracked.string('slimmedGenJets'),
    csvBTagName          = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    mvaBTagName          = cms.untracked.string('pfCombinedMVAV2BJetTags'),
    cvlcTagName          = cms.untracked.string('pfCombinedCvsLJetTags'),
    cvbcTagName          = cms.untracked.string('pfCombinedCvsBJetTags'),
    qgLikelihood         = cms.untracked.string('QGTagger'),
    deepCSVBTagName      = cms.untracked.string('pfDeepCSVJetTags'),
    deepCMVABTagName     = cms.untracked.string('pfDeepCMVAJetTags'),
    ),

  AK4Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    useTriggerObject     = cms.untracked.bool(False),
    applyJEC             = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    maxEta               = cms.untracked.double(3),
    coneSize             = cms.untracked.double(0.4),
    addPFCand            = cms.untracked.bool(True),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doComputeSVInfo      = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if options.isData else cms.untracked.bool(True) ),
    showerDecoConf       = cms.untracked.string(''),
    
    edmPVName   = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    jecName     = (cms.untracked.string('ak4PuppiL1FastL2L3ResidualCorrector') if options.isData else cms.untracked.string('ak4PuppiL1FastL2L3Corrector') ),
    jecUncName  = (cms.untracked.string('AK4Puppi')),
    edmRhoName  = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('', 'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    BRegNNFileName          = cms.untracked.string('BaconProd/Utils/data/breg_training_2016_updated.pb'),
    BRegNNMean              = cms.untracked.double(1.0454729795455933),
    BRegNNStd               = cms.untracked.double(0.31628304719924927),
    # names of various jet-related collections
    jetName            = cms.untracked.string('AK4PFJetsPuppi'),
    genJetName         = cms.untracked.string('AK4GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK4FlavorPuppi'),
    prunedJetName      = cms.untracked.string('AK4caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('AK4caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    mvaBTagName        = cms.untracked.string('AK4PFCombinedMVAV2BJetTagsPuppi'),
    cvlcTagName        = cms.untracked.string('AK4PFCombinedCvsLJetTagsPuppi'),
    cvbcTagName        = cms.untracked.string('AK4PFCombinedCvsBJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    csvDoubleBTagName  = cms.untracked.string('AK4PFBoostedDoubleSecondaryVertexBJetTagsPuppi'),
    deepCSVBTagName    = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi'),
    deepCSVBTagNameb   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probb'),
    deepCSVBTagNamec   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probc'),
    deepCSVBTagNamel   = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probudsg'),
    deepCSVBTagNamebb  = cms.untracked.string('AK4PFDeepCSVJetTagsPuppi:probbb'),
    deepCMVABTagName   = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi'),
    deepCMVABTagNameb  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probb'),
    deepCMVABTagNamec  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probc'),
    deepCMVABTagNamel  = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probudsg'),
    deepCMVABTagNamebb = cms.untracked.string('AK4PFDeepCMVAJetTagsPuppi:probbb'),
    deepDoubleBvLTagName = cms.untracked.string('AK4PFBoostedDeepDoubleBvLJetTagsPuppi:probHbb'),
    deepDoubleCvLTagName = cms.untracked.string('AK4PFBoostedDeepDoubleCvLJetTagsPuppi:probHcc'),
    deepDoubleCvBTagName = cms.untracked.string('AK4PFBoostedDeepDoubleCvBJetTagsPuppi:probHcc'),
    deepDoubleBvLNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleBvLNoMassSculptPenJetTagsPuppi:probHbb'),
    deepDoubleCvLNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleCvLNoMassSculptPenJetTagsPuppi:probHcc'),
    deepDoubleCvBNoMassSculptPenTagName   = cms.untracked.string('AK4PFBoostedDeepDoubleCvBNoMassSculptPenJetTagsPuppi:probHcc'),
    boostedDoubleSVTagInfoName = cms.untracked.string('AK4PFBoostedDoubleSVTagInfosPuppi'), 
    secVertices        = cms.untracked.string('slimmedSecondaryVertices'),
    edmMuonName        = cms.untracked.string('slimmedMuons'),
    edmElectronName    = cms.untracked.string('slimmedElectrons'),
    softPFMuonTagInfoName     = cms.untracked.string('AK4PFSoftPFMuonsTagInfosPuppi'),
    softPFElectronTagInfoName = cms.untracked.string('AK4PFSoftPFElectronsTagInfosPuppi'),
    jettiness          = cms.untracked.string('AK4NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('AK4QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('AK4QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('')
  ),

                                
  
  PFCand = cms.untracked.PSet(
    isActive       = cms.untracked.bool(False),
    edmName        = cms.untracked.string('packedPFCandidates'),
    edmPVName      = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    doAddDepthTime = cms.untracked.bool(False)
  )
)

# overwrite parameters for different eras 
if options.era == '2017' or options.era == '2018':
    #process.ntupler.TriggerObject = cms.untracked.string("slimmedPatTrigger")

    #electron
    process.ntupler.Electron.edmEleMediumMVA = cms.untracked.string('mvaEleID-Fall17-noIso-V1-wp90')
    process.ntupler.Electron.edmEleTightMVA = cms.untracked.string('mvaEleID-Fall17-noIso-V1-wp80')
    #process.ntupler.Electron.edmEleMVA = cms.untracked.string('ElectronMVAEstimatorRun2Fall17NoIsoV1')
    #process.ntupler.Electron.edmEleMediumMVAIso = cms.untracked.string('mvaEleID-Fall17-iso-V1-wp90')
    #process.ntupler.Electron.edmEleTightMVAIso = cms.untracked.string('mvaEleID-Fall17-iso-V1-wp80')
    #process.ntupler.Electron.edmEleMVAIso = cms.untracked.string('ElectronMVAEstimatorRun2Fall17IsoV1')
    #process.ntupler.Electron.edmEleMVAHZZ = cms.untracked.string('')
    #process.ntupler.Electron.storeSecondMVA = cms.untracked.bool(True)
    #process.ntupler.Electron.storeHZZMVA = cms.untracked.bool(False)

    #AK4CHS
    process.ntupler.AK4CHS.BRegNNFileName     = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK4CHS.BRegNNMean         = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK4CHS.BRegNNStd          = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK4CHS.deepCMVABTagNameb  = cms.untracked.string('ak4pfDeepCMVAJetTags:probb')
    process.ntupler.AK4CHS.deepCMVABTagNamec  = cms.untracked.string('ak4pfDeepCMVAJetTags:probc')
    process.ntupler.AK4CHS.deepCMVABTagNamel  = cms.untracked.string('ak4pfDeepCMVAJetTags:probudsg')
    process.ntupler.AK4CHS.deepCMVABTagNamebb = cms.untracked.string('ak4pfDeepCMVAJetTags:probbb')

    #AK4Puppi
    process.ntupler.AK4Puppi.BRegNNFileName = cms.untracked.string('BaconProd/Utils/data/breg_training_2017.pb')
    process.ntupler.AK4Puppi.BRegNNMean     = cms.untracked.double(1.0610932111740112)
    process.ntupler.AK4Puppi.BRegNNStd      = cms.untracked.double(0.39077115058898926)
    process.ntupler.AK4Puppi.genJetName     = cms.untracked.string('slimmedGenJets')



if options.isData:
    process.ntupler.GenInfo.fillLHEWeights = cms.untracked.bool(False)
    process.ntupler.AK4CHS.deepCMVABTagName = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNameb = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamec = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamel = cms.untracked.string('')
    process.ntupler.AK4CHS.deepCMVABTagNamebb = cms.untracked.string('')
    

if options.isData:
    process.baconSequence = cms.Sequence(
                                         process.BadPFMuonFilter                     *
                                         process.BadChargedCandidateFilter           *
                                         process.ak4chsL1FastL2L3ResidualChain       *
                                         process.ak4PuppiL1FastL2L3ResidualChain     *
                                         process.ak8PuppiL1FastL2L3ResidualChain     *
                                         process.ak4chsL1FastL2L3Corrector           *
                                         process.ak4PuppiL1FastL2L3Corrector         *
                                         process.pfNoPileUpJME                       *
                                         process.egammaPostRecoSeq                   *
                                         process.prefiringweight                     *
                                         process.puppiMETSequence                    *
                                         process.AK4jetsequencePuppiData             *
                                         process.patJetCorrFactors                   *
                                         process.updatedPatJets                      *
                                         process.btagging                            *
                                         process.fullPatMetSequenceV2                *
                                         process.fullPatMetSequencePuppi             *
                                         process.patJetCorrFactorsTransientCorrected *
                                         process.updatedPatJetsTransientCorrected    *
                                         process.selectedUpdatedPatJets              *
                                         process.QGTagger                            *
                                         process.ntupler
                                         )
else:
    if options.era == '2016':
        process.baconSequence = cms.Sequence(
                                             process.BadPFMuonFilter          *
                                             process.BadChargedCandidateFilter*
                                             process.ak4chsL1FastL2L3Chain    *
                                             process.ak4PuppiL1FastL2L3Chain  *
                                             process.ak8PuppiL1FastL2L3Chain  *
                                             process.pfNoPileUpJME            *
                                             process.egammaPostRecoSeq        *
                                             process.prefiringweight          *
                                             process.puppiMETSequence         *
                                             process.genjetsequence           *
                                             process.AK4genjetsequenceCHS     *
                                             process.AK4jetsequencePuppi      *
                                             process.patJetCorrFactors*
                                             process.updatedPatJets*
                                             process.btagging                 *
                                             process.fullPatMetSequenceV2     *
                                             process.fullPatMetSequencePuppi  *
                                             process.patJetCorrFactorsTransientCorrected*
                                             process.updatedPatJetsTransientCorrected*
                                             process.selectedUpdatedPatJets*
                                             process.QGTagger              *
                                             process.ntupler)
    elif options.era == '2017' or options.era == '2018':
        process.baconSequence = cms.Sequence(
                                             process.BadPFMuonFilter          *
                                             process.BadChargedCandidateFilter*
                                             process.ak4chsL1FastL2L3Chain    *
                                             process.ak4PuppiL1FastL2L3Chain  *
                                             process.ak8PuppiL1FastL2L3Chain  *
                                             process.pfNoPileUpJME            *
                                             process.egammaPostRecoSeq        *
                                             process.prefiringweight          *
                                             process.puppiMETSequence         *
                                             process.genjetsequence           *
                                             process.AK4genjetsequenceCHS     *
                                             process.AK4jetsequencePuppi      *
                                             process.patJetCorrFactors*
                                             process.updatedPatJets*
                                             process.btagging                 *
                                             process.fullPatMetSequenceV2     *
                                             #process.fullPatMetSequencePuppi  *
                                             process.patJetCorrFactorsTransientCorrected*
                                             process.updatedPatJetsTransientCorrected*
                                             process.selectedUpdatedPatJets*
                                             process.QGTagger              *
                                             process.ntupler)

#--------------------------------------------------------------------------------
# apply trigger filter, if necessary
#================================================================================
if options.doHLTFilter:
  process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
  process.hltHighLevel.throw = cms.bool(False)
  process.hltHighLevel.HLTPaths = cms.vstring()
  hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
  for line in hlt_file.readlines():
    line = line.strip()              # strip preceding and trailing whitespaces
    if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
      hlt_path = line.split()[0]
      process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))
  process.p = cms.EndPath(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.EndPath(process.baconSequence)

#--------------------------------------------------------------------------------
# simple checks to catch some mistakes...
#================================================================================
if options.isData:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
  assert process.ntupler.AK4CHS.doGenJet  == cms.untracked.bool(False)

#process.out = cms.OutputModule("PoolOutputModule",                                                                                                                                                                                                                
#                               outputCommands = cms.untracked.vstring('keep *'),                                                                                                                                                                                
#                               fileName       = cms.untracked.string ("test.root")                                                                                                                                                                              
#                               )                                                                                                                                                                                                                                
#process.endpath = cms.EndPath(process.out)    

