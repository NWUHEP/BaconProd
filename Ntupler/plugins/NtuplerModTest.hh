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
  class FillerCaloJet;
  class FillerGenJet;
  class FillerJet;
  class FillerPF;
  class FillerRH;
}


class NtuplerModTest : public edm::EDAnalyzer {
  public:
    explicit NtuplerModTest(const edm::ParameterSet &iConfig);
    ~NtuplerModTest();

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
    bool fUseAODJet,      fUseAODFatJet,      fUseAODFatterJet;
    bool fUseAODPuppiJet, fUseAODFatPuppiJet, fUseAODFatterPuppiJet;

    bool fComputeFullJetInfo,      fComputeFullFatJetInfo,      fComputeFullFatterJetInfo;
    bool fComputeFullPuppiJetInfo, fComputeFullFatPuppiJetInfo, fComputeFullFatterPuppiJetInfo;

    bool fComputeFullSVInfo,      fComputeFullFatSVInfo,      fComputeFullFatterSVInfo;
    bool fComputeFullPuppiSVInfo, fComputeFullFatPuppiSVInfo, fComputeFullFatterPuppiSVInfo;
 
    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerGenJets   *fFillerGenJet,*fFillerGenFatJet;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerCaloJet   *fFillerCaloJet;
    baconhep::FillerJet       *fFillerFatJet, *fFillerFatterJet;
    baconhep::FillerJet       *fFillerFatPuppiJet, *fFillerFatterPuppiJet;
    baconhep::FillerPF        *fFillerPF;    
    baconhep::FillerRH        *fFillerRH;    
    
    baconhep::TTrigger        *fTrigger;
    
    bool fIsActiveEvtInfo;
    bool fIsActiveGenInfo;
    bool fIsActiveGenJet;
    bool fIsActiveGenFatJet;
    bool fIsActivePV;
    bool fIsActiveCaloJet;
    bool fIsActiveFatJet, fIsActiveFatterJet;
    bool fIsActiveFatPuppiJet, fIsActiveFatterPuppiJet;
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
    TClonesArray            *fGenFatJetArr;
    TClonesArray	    *fCaloJetArr;
    TClonesArray	    *fFatJetArr, *fFatterJetArr;
    TClonesArray	    *fFatPuppiJetArr, *fFatterPuppiJetArr;
    TClonesArray	    *fPVArr;
    TClonesArray	    *fAddFatJetArr, *fAddFatterJetArr;
    TClonesArray	    *fAddFatPuppiJetArr, *fAddFatterPuppiJetArr;
    TClonesArray	    *fPFParArr;
    TClonesArray	    *fRHParArr;
    TClonesArray            *fSVArr;
};
