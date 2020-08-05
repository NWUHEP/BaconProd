#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h" 
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h" 
//#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
//#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
//#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include <TClonesArray.h>
#include <TLorentzVector.h>
//#include <TVector3.h>
#include <TMath.h>
#include <map>
#include <vector>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerPhoton::FillerPhoton(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt             (iConfig.getUntrackedParameter<double>("minPt",10)),
  fPhotonName        (iConfig.getUntrackedParameter<std::string>("edmName","gedPhotons")),
  fPFCandName        (iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fBSName            (iConfig.getUntrackedParameter<std::string>("edmBeamspotName","offlineBeamSpot")),
  fEleName           (iConfig.getUntrackedParameter<std::string>("edmElectronName","gedGsfElectrons")),
  fConvName          (iConfig.getUntrackedParameter<std::string>("edmConversionName","allConversions")),
  fSCName            (iConfig.getUntrackedParameter<edm::InputTag>("edmSCName")),//,"particleFlowEGamma")),
  fMVASpring16       (iConfig.getUntrackedParameter<std::string>("edmPhoMVASpring16")),
  fMVAFall17V1       (iConfig.getUntrackedParameter<std::string>("edmPhoMVAFall17V1")),
  fMVAFall17V2       (iConfig.getUntrackedParameter<std::string>("edmPhoMVAFall17V2")),

  //feeReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection")),
  //febReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection")),
  //fesReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection")),
 
  feeReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection")),
  fesReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection")),
  febReducedRecHitCollection          (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection")),

  fphotonChargedIsolationCollection         (iConfig.getParameter<edm::InputTag>("phoChargedIsolationCollection")),
  fphotonPhotonIsolationCollection          (iConfig.getParameter<edm::InputTag>("phoPhotonIsolationCollection")),
  fphotonNeutralHadronIsolationCollection   (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolationCollection")),
  fphotonWorstChargeIsolationCollection     (iConfig.getParameter<edm::InputTag>("phoWorstChargeIsolationCollection")),


  fUseTO             (iConfig.getUntrackedParameter<bool>("useTriggerObject",false)),
  fUseAOD            (useAOD)
{
//  fPhotonMVA = new PhotonMVACalculator();
  if(fUseAOD)  fTokPhotonName    =  iC.consumes<reco::PhotonCollection>     (fPhotonName);
  if(!fUseAOD) fTokPatPhotonName =  iC.consumes<pat::PhotonCollection>      (fPhotonName);
  fTokPFCandName =  iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  fTokBSName     =  iC.consumes<reco::BeamSpot>             (fBSName);
  fTokEleName    =  iC.consumes<reco::GsfElectronCollection>(fEleName);
  fTokConvName   =  iC.consumes<reco::ConversionCollection>(fConvName);
  fTokSCName     =  iC.consumes<reco::SuperClusterCollection>(fSCName);

  fTokeeReducedRecHitCollection  =  iC.consumes<EcalRecHitCollection>(feeReducedRecHitCollection);
  fTokesReducedRecHitCollection  =  iC.consumes<EcalRecHitCollection>(fesReducedRecHitCollection);
  fTokebReducedRecHitCollection  =  iC.consumes<EcalRecHitCollection>(febReducedRecHitCollection);

  fTokphotonChargedIsolationCollection       =  iC.consumes<edm::ValueMap<float>>(fphotonChargedIsolationCollection);
  fTokphotonPhotonIsolationCollection        =  iC.consumes<edm::ValueMap<float>>(fphotonPhotonIsolationCollection);
  fTokphotonNeutralHadronIsolationCollection =  iC.consumes<edm::ValueMap<float>>(fphotonNeutralHadronIsolationCollection);
  fTokphotonWorstChargeIsolationCollection   =  iC.consumes<edm::ValueMap<float>>(fphotonWorstChargeIsolationCollection);
}

//--------------------------------------------------------------------------------------------------
FillerPhoton::~FillerPhoton(){}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerPhoton::fill(TClonesArray *array, 
                        const edm::Event &iEvent, const edm::EventSetup &iSetup,
		        const std::vector<TriggerRecord> &triggerRecords,
		        const trigger::TriggerEvent &triggerEvent)
{
  assert(array);
//  if(!fPhotonReg->IsInitialized()) { 
//      std::string cmssw_base_utils = getenv("CMSSW_BASE"); 
//      cmssw_base_utils+="/src/BaconProd/Utils/data/";
//      fPhotonReg->Initialize(iSetup,cmssw_base_utils+"gbrv3ph_52x.root");
//      fPhotonMVA->initialize(cmssw_base_utils+"2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15.weights.xml",cmssw_base_utils+"/2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15.weights.xml");
//  }

  
  // Get photon collection
  edm::Handle<reco::PhotonCollection> hPhotonProduct;
  iEvent.getByToken(fTokPhotonName,hPhotonProduct);
  assert(hPhotonProduct.isValid());
  const reco::PhotonCollection *photonCol = hPhotonProduct.product();
  /*
  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  //iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  */
  // Get beamspot
  edm::Handle<reco::BeamSpot> hBeamSpotProduct;
  iEvent.getByToken(fTokBSName,hBeamSpotProduct);
  assert(hBeamSpotProduct.isValid());
  const reco::BeamSpot *bs = hBeamSpotProduct.product();

  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> hEleProduct;
  iEvent.getByToken(fTokEleName,hEleProduct);
  assert(hEleProduct.isValid());
  
  // Get conversions collection
  edm::Handle<reco::ConversionCollection> hConvProduct;
  iEvent.getByToken(fTokConvName,hConvProduct);
  assert(hConvProduct.isValid());

  // Get SuperCluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fTokSCName,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();

  
  //edm::Handle<edm::View<pat::Photon> > photonHandle;
  //iEvent.getByToken(fPhotonName, photonHandle);
  edm::Handle<edm::View<pat::PhotonCollection>> handlePhoton;
  iEvent.getByToken(fTokPatPhotonName,handlePhoton);
  assert(handlePhoton.isValid());

  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;

  iEvent.getByToken(fTokphotonChargedIsolationCollection      , phoChargedIsolationMap      );
  iEvent.getByToken(fTokphotonPhotonIsolationCollection       , phoNeutralHadronIsolationMap);
  iEvent.getByToken(fTokphotonNeutralHadronIsolationCollection, phoPhotonIsolationMap       );
  iEvent.getByToken(fTokphotonWorstChargeIsolationCollection  , phoWorstChargedIsolationMap );
  //-------------------------------------------------------------
  //--------------------Added by Joseph Cordero -----------------
  // Get ebEcalHit Collection
  //edm::Handle<EcalRecHitCollection> hebCalProduct;
  //iEvent.getByToken(fTokebReducedRecHitCollection, hebCalProduct);
  //assert(hebCalProduct.isValid());
  //const EcalRecHitCollection *ebCalCol = hebCalProduct.product();

  //// Get eeEcalHit Collection
  //edm::Handle<EcalRecHitCollection> heeCalProduct;
  //iEvent.getByToken(fTokeeReducedRecHitCollection, heeCalProduct);
  //assert(heeCalProduct.isValid());
  //const EcalRecHitCollection *eeCalCol = heeCalProduct.product();

  //// Get esEcalHit Collection
  //edm::Handle<EcalRecHitCollection> hesCalProduct;
  //iEvent.getByToken(fTokesReducedRecHitCollection, hesCalProduct);
  //assert(hesCalProduct.isValid());
  //const EcalRecHitCollection *esCalCol = hesCalProduct.product();
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  
 
  //EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  //EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebCalCol, eeCalCol, esCalCol);
  //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebCalCol, eeCalCol, esCalCol);

  EcalClusterLazyTools       lazyTool    (iEvent, iSetup, fTokebReducedRecHitCollection, fTokeeReducedRecHitCollection, fTokesReducedRecHitCollection);
  noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, fTokebReducedRecHitCollection, fTokeeReducedRecHitCollection, fTokesReducedRecHitCollection);

  unsigned int nPhoton;
  for(reco::PhotonCollection::const_iterator itPho = photonCol->begin(); itPho!=photonCol->end(); ++itPho) {
    nPhoton = itPho - photonCol->begin();

    // Photon cuts
    if(itPho->pt() < fMinPt) continue;
    
    // construct object and place in array
    TClonesArray &rPhotonArr = *array;
    assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
    const int index = rPhotonArr.GetEntries();
    new(rPhotonArr[index]) baconhep::TPhoton();
    baconhep::TPhoton *pPhoton = (baconhep::TPhoton*)rPhotonArr[index];

    const reco::SuperClusterRef sc = itPho->superCluster();

    // ref to access value maps
    edm::RefToBase<reco::Photon> phoBaseRef( edm::Ref<reco::PhotonCollection>(hPhotonProduct, itPho - photonCol->begin()) );


    //
    // Kinematics
    //==============================
    pPhoton->pt  = itPho->pt();
    pPhoton->eta = itPho->eta();
    pPhoton->phi = itPho->phi();

    pPhoton->scEt  = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pPhoton->scEta = sc->eta();
    pPhoton->scPhi = sc->phi();


    //
    // Isolation
    //==============================
    pPhoton->trkIso  = itPho->trkSumPtHollowConeDR04();
    pPhoton->ecalIso = itPho->ecalRecHitSumEtConeDR04();
    pPhoton->hcalIso = itPho->hcalTowerSumEtConeDR04();
 
    //Isolation for Photon MVA
//    pPhoton->chHadIso03SelVtx  = -1;
//    pPhoton->chHadIso03WstVtx  = -1;
//    computeVtxIso(*itPho,*pfCandCol,*vtxCol,
//		  pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx);

    //
    // Identification
    //==============================
    pPhoton->hovere     = itPho->hadronicOverEm();
    pPhoton->sthovere   = itPho->hadTowOverEm();
    pPhoton->sieie      = itPho->full5x5_sigmaIetaIeta();
    //pPhoton->sipip      = 0;  // (!) todo (lazy tools)
    pPhoton->r9         = itPho->r9();  // (!) change to full5x5 after 7_2_0 MC?
    pPhoton->r9_full5x5 = itPho->full5x5_r9();

    //---------------Added by Joseph Cordero Mercado ---------------
    //----------Date: 2019/07/30 ----------------------------------
    //----------email: jcordero@u.northwestern.edu------------------
     
    //pPhoton->e2x2       = lazyTool->e2x2(sc->seed());
    pPhoton->e2x2       = lazyTool.e2x2(*(sc->seed()));
    pPhoton->e5x5       = itPho->e5x5();
    pPhoton->scEtaWidth = sc->etaWidth();
    pPhoton->scPhiWidth = sc->phiWidth();
    pPhoton->scBrem     = sc->phiWidth()/sc->etaWidth();
    pPhoton->scRawE     = sc->rawEnergy();
    pPhoton->scESEn     = sc->preshowerEnergy();

    pPhoton->phoChIso      =  (phoChargedIsolationMap->begin()      )[nPhoton]; 
    pPhoton->phoNeuHadIso  =  (phoNeutralHadronIsolationMap->begin())[nPhoton]; 
    pPhoton->phoPhIso      =  (phoPhotonIsolationMap->begin()       )[nPhoton]; 
    pPhoton->phoWorstChIso =  (phoWorstChargedIsolationMap->begin() )[nPhoton];

    std::vector<float> vCov = lazyToolnoZS.localCovariances(*(sc->seed()));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    pPhoton->srr        = lazyToolnoZS.eseffsirir(*(sc));
    pPhoton->sipip      = spp;
    pPhoton->sieip      = sep;
    //--------------------------------------------------------------

    pPhoton->fiducialBits=0;
    if(itPho->isEB())        pPhoton->fiducialBits |= kIsEB;
    if(itPho->isEE())        pPhoton->fiducialBits |= kIsEE;
    if(itPho->isEBEEGap())   pPhoton->fiducialBits |= kIsEBEEGap;
    if(itPho->isEBEtaGap())  pPhoton->fiducialBits |= kIsEBEtaGap;
    if(itPho->isEBPhiGap())  pPhoton->fiducialBits |= kIsEBPhiGap;
    if(itPho->isEEDeeGap())  pPhoton->fiducialBits |= kIsEEDeeGap;
    if(itPho->isEERingGap()) pPhoton->fiducialBits |= kIsEERingGap;

    pPhoton->typeBits=0;
    if(itPho->isStandardPhoton()) pPhoton->typeBits |= baconhep::kEGamma;
    if(itPho->isPFlowPhoton())    pPhoton->typeBits |= baconhep::kPFPhoton;

    // Obtain a supercluster ID, unique per event. The SC ID is the index in the SC collection
    pPhoton->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator itSC = scCol->begin(); itSC!=scCol->end(); ++itSC) {
      scIndex++;
      if(itPho->superCluster().get() == &(*itSC)) {
        pPhoton->scID = scIndex;
	break;
      }
    }
    
    pPhoton->hasPixelSeed     = itPho->hasPixelSeed();
    pPhoton->isConv           = ConversionTools::hasMatchedPromptElectron(itPho->superCluster(), hEleProduct, hConvProduct, bs->position(), 2.0, 1e-6, 0);
    pPhoton->passElectronVeto = !(pPhoton->isConv); // here for backwards compatibility
 
    if(fUseTO) pPhoton->hltMatchBits = TriggerTools::matchHLT(pPhoton->eta, pPhoton->phi, triggerRecords, triggerEvent);


  }
}

// === filler for MINIAOD ===
void FillerPhoton::fill(TClonesArray *array,
                        const edm::Event &iEvent, const edm::EventSetup &iSetup,
                        const std::vector<TriggerRecord> &triggerRecords,
                        const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);

  // Get photon collection
  edm::Handle<pat::PhotonCollection> hPhotonProduct;
  iEvent.getByToken(fTokPatPhotonName,hPhotonProduct);
  assert(hPhotonProduct.isValid());
  const pat::PhotonCollection *photonCol = hPhotonProduct.product();

  // Get supercluster collection
  edm::Handle<reco::SuperClusterCollection> hSCProduct;
  iEvent.getByToken(fTokSCName,hSCProduct);
  assert(hSCProduct.isValid());
  const reco::SuperClusterCollection *scCol = hSCProduct.product();



  edm::Handle<edm::View<pat::PhotonCollection>> handlePhoton;
  iEvent.getByToken(fTokPatPhotonName,handlePhoton);
  assert(handlePhoton.isValid());

  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;

  iEvent.getByToken(fTokphotonChargedIsolationCollection      , phoChargedIsolationMap      );
  iEvent.getByToken(fTokphotonPhotonIsolationCollection       , phoNeutralHadronIsolationMap);
  iEvent.getByToken(fTokphotonNeutralHadronIsolationCollection, phoPhotonIsolationMap       );
  iEvent.getByToken(fTokphotonWorstChargeIsolationCollection  , phoWorstChargedIsolationMap );

 
  //EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  //EcalClusterLazyTools       lazyTool    (iEvent, iSetup, ebCalCol, eeCalCol, esCalCol);
  //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebCalCol, eeCalCol, esCalCol);

  EcalClusterLazyTools       lazyTool    (iEvent, iSetup, fTokebReducedRecHitCollection, fTokeeReducedRecHitCollection, fTokesReducedRecHitCollection);
  noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, fTokebReducedRecHitCollection, fTokeeReducedRecHitCollection, fTokesReducedRecHitCollection)
;
  unsigned int nPhoton;
  //for(unsigned int nPhoton = 0; nPhoton < photonCol->size(); ++nPhoton) {
    //pat::PhotonCollection& itPho = (* photonCol)[i]

  for(pat::PhotonCollection::const_iterator itPho = photonCol->begin(); itPho!=photonCol->end(); ++itPho) {
    nPhoton = itPho-photonCol->begin();

    // Photon cuts
    if(itPho->pt() < fMinPt) continue;

    // construct object and place in array
    TClonesArray &rPhotonArr = *array;
    assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
    const int index = rPhotonArr.GetEntries();
    new(rPhotonArr[index]) baconhep::TPhoton();
    baconhep::TPhoton *pPhoton = (baconhep::TPhoton*)rPhotonArr[index];

    const reco::SuperClusterRef sc = itPho->superCluster();

    // ref to access value maps
    edm::RefToBase<pat::Photon> phoBaseRef( edm::Ref<pat::PhotonCollection>(hPhotonProduct, itPho - photonCol->begin()) );


    //
    // Kinematics
    //==============================
    pPhoton->pt  = itPho->pt();
    pPhoton->eta = itPho->eta();
    pPhoton->phi = itPho->phi();

    pPhoton->scEt  = (sc->energy())*(sc->position().Rho())/(sc->position().R());
    pPhoton->scEta = sc->eta();
    pPhoton->scPhi = sc->phi();
    
    auto corrP4  = itPho->p4() * itPho->userFloat("ecalEnergyPostCorr") / itPho->energy();
    pPhoton->calibPt = corrP4.Pt();
    pPhoton->calibE = corrP4.E();
    pPhoton->ecalEnergyErrPostCorr = itPho->userFloat("ecalEnergyErrPostCorr");
    pPhoton->energyScaleStatUp = itPho->userFloat("energyScaleStatUp");
    pPhoton->energyScaleStatDown = itPho->userFloat("energyScaleStatDown");
    pPhoton->energyScaleSystUp = itPho->userFloat("energyScaleSystUp");
    pPhoton->energyScaleSystDown = itPho->userFloat("energyScaleSystDown");
    pPhoton->energyScaleGainUp = itPho->userFloat("energyScaleGainUp");
    pPhoton->energyScaleGainDown = itPho->userFloat("energyScaleGainDown");
    pPhoton->energySigmaRhoUp = itPho->userFloat("energySigmaRhoUp");
    pPhoton->energySigmaRhoDown = itPho->userFloat("energySigmaRhoDown");
    pPhoton->energySigmaPhiUp = itPho->userFloat("energySigmaPhiUp");
    pPhoton->energySigmaPhiDown = itPho->userFloat("energySigmaPhiDown");
    
    pPhoton->eRes = itPho->getCorrectedEnergyError(itPho->getCandidateP4type())/itPho->getCorrectedEnergy(itPho->getCandidateP4type());

    //
    // Isolation
    //==============================
    pPhoton->trkIso  = itPho->trkSumPtHollowConeDR04();
    pPhoton->ecalIso = itPho->ecalRecHitSumEtConeDR04();
    pPhoton->hcalIso = itPho->hcalTowerSumEtConeDR04();

    pPhoton->chHadIso  = itPho->chargedHadronIso();
    pPhoton->gammaIso  = itPho->photonIso();       //(*hGammaIsoMap)[phoBaseRef];
    pPhoton->neuHadIso = itPho->neutralHadronIso();//(*hNeuHadIsoMap)[phoBaseRef];

    //Isolation for Photon MVA
//    pPhoton->chHadIso03SelVtx  = -1;
//    pPhoton->chHadIso03WstVtx  = -1;
//    computeVtxIso(*itPho,*pfCandCol,*vtxCol,
//                pPhoton->chHadIso03SelVtx,pPhoton->chHadIso03WstVtx);


    //
    // Identification
    //==============================
    pPhoton->hovere     = itPho->hadronicOverEm();
    pPhoton->sthovere   = itPho->hadTowOverEm();
    pPhoton->sieie      = itPho->full5x5_sigmaIetaIeta();
    //pPhoton->sipip      = 0;  // (!) todo (lazy tools)
    pPhoton->r9         = itPho->r9();  // (!) change to full5x5 after 7_2_0 MC?
    pPhoton->r9_full5x5 = itPho->full5x5_r9();

    //---------------Added by Joseph Cordero Mercado ---------------
    //----------Date: 2019/07/30 ----------------------------------
    //----------email: jcordero@u.northwestern.edu------------------
     
    pPhoton->e2x2       = lazyTool.e2x2(*(sc->seed()));
    pPhoton->e5x5       = itPho->e5x5();
    pPhoton->scEtaWidth = sc->etaWidth();
    pPhoton->scPhiWidth = sc->phiWidth();
    pPhoton->scBrem     = sc->phiWidth()/sc->etaWidth();
    pPhoton->scRawE     = sc->rawEnergy();
    pPhoton->scESEn     = sc->preshowerEnergy();

    pPhoton->phoChIso      =  (phoChargedIsolationMap->begin()      )[nPhoton]; 
    pPhoton->phoNeuHadIso  =  (phoNeutralHadronIsolationMap->begin())[nPhoton]; 
    pPhoton->phoPhIso      =  (phoPhotonIsolationMap->begin()       )[nPhoton]; 
    pPhoton->phoWorstChIso =  (phoWorstChargedIsolationMap->begin() )[nPhoton]; 

    std::vector<float> vCov = lazyToolnoZS.localCovariances(*(sc->seed()));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    pPhoton->srr        = lazyToolnoZS.eseffsirir(*(sc));
    pPhoton->sipip      = spp;
    pPhoton->sieip      = sep;
    //--------------------------------------------------------------
    
    pPhoton->fiducialBits=0;
    if(itPho->isEB())        pPhoton->fiducialBits |= kIsEB;
    if(itPho->isEE())        pPhoton->fiducialBits |= kIsEE;
    if(itPho->isEBEEGap())   pPhoton->fiducialBits |= kIsEBEEGap;
    if(itPho->isEBEtaGap())  pPhoton->fiducialBits |= kIsEBEtaGap;
    if(itPho->isEBPhiGap())  pPhoton->fiducialBits |= kIsEBPhiGap;
    if(itPho->isEEDeeGap())  pPhoton->fiducialBits |= kIsEEDeeGap;
    if(itPho->isEERingGap()) pPhoton->fiducialBits |= kIsEERingGap;

    pPhoton->typeBits=0;
    if(itPho->isStandardPhoton()) pPhoton->typeBits |= baconhep::kEGamma;
    if(itPho->isPFlowPhoton())    pPhoton->typeBits |= baconhep::kPFPhoton;

    // Obtain a supercluster ID, unique per event. The SC ID is the index in the SC collection
    pPhoton->scID = -1;
    int scIndex = -1;
    for(reco::SuperClusterCollection::const_iterator itSC = scCol->begin(); itSC!=scCol->end(); ++itSC) {
      scIndex++;
      if(itPho->superCluster().get() == &(*itSC)) {
        pPhoton->scID = scIndex;
        break;
      }
    }

    pPhoton->hasPixelSeed     = itPho->hasPixelSeed();
    pPhoton->isConv           = !itPho->passElectronVeto();
    pPhoton->passElectronVeto = itPho->passElectronVeto(); // here for backwards compatibility

    // Photon MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2
    //pPhoton->mva = itPho->userFloat(fMVA+"Values"); 
    //pPhoton->mvaCat = itPho->userInt(fMVA+"Categories"); 
    //
    //// v2 photon MVA
    //pPhoton->mvaV2 = itPho->userFloat(fMVAV2+"Values");
    //pPhoton->mvaV2Cat = itPho->userInt(fMVAV2+"Categories");
    
    pPhoton->mvaSpring16    = itPho->userFloat(fMVASpring16+"Values");
    pPhoton->mvaSpring16Cat = itPho->userInt(fMVASpring16+"Categories"); 
    pPhoton->mvaFall17V1    = itPho->userFloat(fMVAFall17V1+"Values");
    pPhoton->mvaFall17V1Cat = itPho->userInt(fMVAFall17V1+"Categories");
    pPhoton->mvaFall17V2    = itPho->userFloat(fMVAFall17V2+"Values");
    pPhoton->mvaFall17V2Cat = itPho->userInt(fMVAFall17V2+"Categories");

    if(fUseTO) pPhoton->hltMatchBits = TriggerTools::matchHLT(pPhoton->eta, pPhoton->phi, triggerRecords, triggerObjects);

  }
}
/*
//--------------------------------------------------------------------------------------------------
void FillerPhoton::computeVtxIso(const reco::Photon &photon,
				 const std::vector<reco::PFCandidate>        &pf,
				 const std::vector<reco::Vertex>             &iVertex,
				 float &out_chHadIsoWvtx,float &out_chHadIsoFirstVtx) const 
{
  double extRadius = 0.3;
  double intRadius = 0.02;
  for(unsigned int iVtx=0; iVtx<iVertex.size(); iVtx++) { 
    const reco::Vertex vertex = iVertex.at(iVtx);
    double pIso = 0; 
    for(unsigned int ipf=0; ipf<pf.size(); ipf++) {      
      const reco::PFCandidate pfcand = pf.at(ipf);
      if(pfcand.particleId() != reco::PFCandidate::h) continue;
      // Add p_T to running sum if PFCandidate is close enough
      TVector3 pVec; 
      pVec.SetXYZ(photon.superCluster()->position().x()-vertex.position().x(),
		  photon.superCluster()->position().y()-vertex.position().y(),
		  photon.superCluster()->position().z()-vertex.position().z());
      
      double dr = reco::deltaR(pfcand.eta(), pfcand.phi(), pVec.Eta(), pVec.Phi());
      if(dr >= extRadius || dr < intRadius) continue;
      double dZ = pfcand.trackRef()   ->dz (vertex.position());
      double d0 = pfcand.trackRef()   ->dxy(vertex.position());
      if(fabs(dZ) > 0.2) continue;
      if(fabs(d0) > 0.1) continue;
      pIso  += pfcand.pt(); 
    }
    if(iVtx == 0)               out_chHadIsoFirstVtx = pIso;
    if(out_chHadIsoWvtx < pIso) out_chHadIsoWvtx     = pIso;
  }
}
*/
