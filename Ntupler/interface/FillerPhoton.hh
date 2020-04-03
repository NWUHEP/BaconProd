#ifndef BACONPROD_NTUPLER_FILLERPHOTON_HH
#define BACONPROD_NTUPLER_FILLERPHOTON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
//#include "BaconProd/Utils/interface/PhotonMVACalculator.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"



class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerPhoton
  {
    public:
      FillerPhoton(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC);
      ~FillerPhoton();

      // === filler for AOD ===
      void fill(TClonesArray			 *array,	                // output array to be filled
                const edm::Event		 &iEvent,	                // event info
	        const edm::EventSetup		 &iSetup,	                // event setup info
	        const std::vector<TriggerRecord> &triggerRecords,               // list of trigger names and objects
	        const trigger::TriggerEvent	 &triggerEvent);                // event trigger objects

      // === filler for MINIAOD ===
      void fill(TClonesArray                                 *array,            // output array to be filled
                const edm::Event                             &iEvent,           // event info
                const edm::EventSetup                        &iSetup,           // event setup info
                const std::vector<TriggerRecord>             &triggerRecords,   // list of trigger names and objects
                const pat::TriggerObjectStandAloneCollection &triggerObjects);  // event trigger objects

    protected:
/*            
      void computeVtxIso    (const reco::Photon &photon,
			     const std::vector<reco::PFCandidate>   &pf,
			     const std::vector<reco::Vertex>        &iVetex,
			     float &out_chHadIsoWvtx,float &out_chHadIsoFirstVtx) const;
      std::vector<float>  getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row);
      std::vector<float>  getESShape(std::vector<float> ESHits0);      
*/

      // Photon cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fPhotonName;
      std::string fPFCandName;
      std::string fBSName;
      std::string fEleName;
      std::string fConvName;
      edm::InputTag fSCName;

      edm::EDGetTokenT<reco::PhotonCollection>       fTokPhotonName;
      edm::EDGetTokenT<pat::PhotonCollection>        fTokPatPhotonName;
      edm::EDGetTokenT<reco::PFCandidateCollection>  fTokPFCandName;
      edm::EDGetTokenT<reco::BeamSpot>               fTokBSName;
      edm::EDGetTokenT<reco::GsfElectronCollection>  fTokEleName;
      edm::EDGetTokenT<reco::ConversionCollection>   fTokConvName;
      edm::EDGetTokenT<reco::SuperClusterCollection> fTokSCName;

      edm::EDGetTokenT<EcalRecHitCollection> fTokebReducedRecHitCollection;
      edm::EDGetTokenT<EcalRecHitCollection> fTokeeReducedRecHitCollection;
      edm::EDGetTokenT<EcalRecHitCollection> fTokesReducedRecHitCollection;

      edm::EDGetTokenT<edm::ValueMap<float>> fTokphotonChargedIsolationCollection;
      edm::EDGetTokenT<edm::ValueMap<float>> fTokphotonPhotonIsolationCollection; 
      edm::EDGetTokenT<edm::ValueMap<float>> fTokphotonNeutralHadronIsolationCollection;
      edm::EDGetTokenT<edm::ValueMap<float>> fTokphotonWorstChargeIsolationCollection;

      std::string fMVASpring16;
      std::string fMVAFall17V1;
      std::string fMVAFall17V2;

      edm::InputTag feeReducedRecHitCollection;
      edm::InputTag fesReducedRecHitCollection;
      edm::InputTag febReducedRecHitCollection;

      edm::InputTag fphotonChargedIsolationCollection;
      edm::InputTag fphotonPhotonIsolationCollection;
      edm::InputTag fphotonNeutralHadronIsolationCollection;
      edm::InputTag fphotonWorstChargeIsolationCollection;

      bool fUseTO;
      bool fUseAOD;

  };
}
#endif
