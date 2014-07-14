#ifndef BACONPROD_UTILS_JETTOOLS_HH
#define BACONPROD_UTILS_JETTOOLS_HH

// forward class declarations

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "fastjet/PseudoJet.hh"
#include "TLorentzVector.h"
#include <TVector2.h>


namespace baconhep {

  class JetTools
  {
    public:
      
      // fraction of pT contributed by particles associated with the primary vertex
      static double beta(const reco::PFJet &jet, const reco::Vertex &pv, const double dzCut=0.2);
      
      // fraction of pT contributed by particles associated with a pile-up vertex
      static double betaStar(const reco::PFJet &jet, const reco::Vertex &pv, const reco::VertexCollection *pvCol, const double dzCut=0.2);
      
      // mean dR of constituents from jet axis
      static double dRMean(const reco::PFJet &jet, const int pfType=-1);
      
      // mean dR^2 of constituents from jet axis
      static double dR2Mean(const reco::PFJet &jet, const int pfType=-1);
      
      // fraction of pT contribution in (dRMax-0.1) < dR < dRMax from constituents of a specific PF-type (-1 => any type)
      static double frac(const reco::PFJet &jet, const double dRMax, const int pfType=-1);
      
      // dz of the leading charged constituent in the jet 
      static double jetDz(const reco::PFJet &jet, const reco::Vertex &pv);
      
      // d0 of the leading charged constituent in the jet
      static double jetD0(const reco::PFJet &jet, const reco::Vertex &pv);
      
      // jet width variables
      static double jetWidth(const reco::PFJet  &jet, const int varType=0, const int pfType=-1);
      static double jetWidth(fastjet::PseudoJet &jet, const int varType=0); 
      
      
      // Check if PF jet passes loose ID cuts
      static bool passPFLooseID(const reco::PFJet &jet);
    
      //Jet Chrge pT weighted 
      static double jetCharge(const reco::PFJet &jet, const bool iSquare=false); 
      static double jetCharge(const reco::PFJet &jet, const double kappa); 
   
      //Sub Jet Quark Gluon
      static double* subJetQG(const reco::PFJet &jet,edm::Handle<reco::PFJetCollection> &subJets,const edm::ValueMap<float> iQGLikelihood,double iConeSize);     
  
      //Sub Jet BTag
      static double* subJetBTag(const reco::PFJet &jet,reco::JetTagCollection &subJetVal,double iConeSize );
 
      //Jet Pull Vector and measure of color flow
      // type == 0: jet pull
      // type == 1: charged component of jet pull
      // type == 2: neutral component of jet pull
      static TVector2 jetPull(const reco::PFJet &jet, const int type=0);

      //Jet Pull Angle an event better measure of color flow
      static double jetPullAngle(const reco::PFJet &jet ,edm::Handle<reco::PFJetCollection> &subJets,double iConeSize);

      //Q-jet volatility http://arxiv.org/pdf/1201.1914.pdf
      static double qJetVolatility(std::vector<fastjet::PseudoJet> &iConstits, int iQJetsN, int iSeed);
      static float findRMS ( std::vector<float> &iQJetMass);
      static float findMean( std::vector<float> &iQJetMass);
  };
}
#endif
