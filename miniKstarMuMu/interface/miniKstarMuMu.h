#ifndef MINIKSTARMUMU_H
#define MINIKSTARMUMU_H

// system include files
#include <memory>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "Utils.h"
#include "B0KstMuMuTreeContent.h"


//
// class declaration
//

class miniKstarMuMu : public edm::EDAnalyzer {
   public:
      explicit miniKstarMuMu (const edm::ParameterSet&);
      ~miniKstarMuMu();

      std::string getMuCat (reco::Muon const& muon);

   private:
      virtual void beginJob ();
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob ();

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection>         vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection>            muonToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> trackToken_;
      edm::EDGetTokenT<reco::BeamSpot>                 beamSpotToken_;

	  // ####################
	  // # HLT-trigger cuts #
	  // ####################
	  double CLMUMUVTX;
	  double LSMUMUBS;
	  double DCAMUMU;
	  double DCAMUBS;
	  double COSALPHAMUMUBS;
	  double MUMINPT;
	  double MUMAXETA;
	  double MINMUMUPT;
	  double MINMUMUINVMASS;
	  double MAXMUMUINVMASS;

	  // ######################
	  // # Pre-selection cuts #
	  // ######################
	  double B0MASSLOWLIMIT;
	  double B0MASSUPLIMIT;
	  double CLB0VTX;
	  double KSTMASSWINDOW;
	  double HADDCASBS;
	  double MINHADPT;
	  double MAXB0PREMASS;

      bool printMsg;

      TTree* theTree;
      B0KstMuMuTreeContent* NTuple;
      Utils* Utility;

};

#endif
