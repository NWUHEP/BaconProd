#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"      // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"     // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"       // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include <string>                                        // string class

// forward class declarations
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


class TFile;
class TH1D;
class TTree;
class TClonesArray;
namespace edm {
  class TriggerResults;
  class TriggerNames;
}
namespace baconhep {
  class TEventInfo;
  class TGenEventInfo;
  class TTrigger;
  class FillerEventInfo;
  class FillerGenInfo;
  class FillerGenJets;
  class FillerVertex;
  class FillerElectron;
  class FillerMuon;
  class FillerPhoton;
  class FillerCaloJet;
  class FillerJet_REDUCED;
  class FillerGenJet;
  class FillerPF;
  class FillerRH;
}


class NtuplerMod_REDUCED : public edm::EDAnalyzer {
  public:
    explicit NtuplerMod_REDUCED(const edm::ParameterSet &iConfig);
    ~NtuplerMod_REDUCED();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void endRun  (const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void endLuminosityBlock  (const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void respondToOpenInputFile(edm::FileBlock const& fb);
    // specify trigger paths of interest
    void setTriggers(bool iUseTriggers);
    
    // initialization from HLT menu; needs to be called on every change in HLT menu
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);
    

    //--------------------------------------------------------------------------------------------------
    //  data members
    //==================================================================================================   
    bool fSkipOnHLTFail;
    bool fUseAOD;
    
    // variables to handle triggers
    edm::ParameterSetID fTriggerNamesID;
    std::string         fHLTEnd;
    edm::InputTag	fHLTTag;
    edm::InputTag       fHLTObjTag;
    std::string         fHLTFile;

    // Collection names
    std::string fPVName;
    std::string fGenRunInfoName;
    bool        fUseRunInfo;
    edm::EDGetTokenT<GenRunInfoProduct>      fTokGenRunInfo;
    edm::EDGetTokenT<edm::TriggerResults>    fTokTrgRes       ;
    edm::EDGetTokenT<trigger::TriggerEvent>    fTokTrgEvt       ;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>    fTokTrgObj       ;

    bool fFillLHEWgt;
    bool fUseAODJet;  
    bool fUseAODPuppiJet; 

    bool fComputeFullJetInfo; 
    bool fComputeFullPuppiJetInfo; 

    bool fComputeFullSVInfo; 
    bool fComputeFullPuppiSVInfo;
 
    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerGenJets   *fFillerGenJet;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerElectron  *fFillerEle;
    baconhep::FillerMuon      *fFillerMuon;
    baconhep::FillerPhoton    *fFillerPhoton;
    baconhep::FillerCaloJet   *fFillerCaloJet;
    baconhep::FillerJet_REDUCED       *fFillerJet_REDUCED;
    baconhep::FillerJet_REDUCED       *fFillerPuppiJet;
    baconhep::FillerPF        *fFillerPF;    
    baconhep::FillerRH        *fFillerRH;    
    
    baconhep::TTrigger        *fTrigger;
    
    bool fIsActiveEvtInfo;
    bool fIsActiveGenInfo;
    bool fIsActiveGenJet;
    bool fIsActivePV;
    bool fIsActiveEle;
    bool fIsActiveMuon;
    bool fIsActivePhoton;
    bool fIsActiveCaloJet;
    bool fIsActiveJet; 
    bool fIsActivePuppiJet; 
    bool fIsActivePF;
    bool fIsActiveRH;
    bool fUseTrigger;
    bool fUseTriggerObject;
  
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TH1D                    *fTotalEvents;
    float                    fXS;
    TTree                   *fEventTree;
    baconhep::TEventInfo    *fEvtInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fLHEWgtArr;
    TClonesArray            *fGenParArr;
    TClonesArray            *fGenJetArr;
    TClonesArray	    *fEleArr;
    TClonesArray	    *fMuonArr;
    TClonesArray	    *fCaloJetArr;
    TClonesArray	    *fJetArr; 
    TClonesArray	    *fPuppiJetArr;
    TClonesArray	    *fPhotonArr;
    TClonesArray	    *fPVArr;
    TClonesArray	    *fAddJetArr;
    TClonesArray	    *fAddPuppiJetArr;
    TClonesArray	    *fPFParArr;
    TClonesArray	    *fRHParArr;
    TClonesArray            *fSVArr;
};
