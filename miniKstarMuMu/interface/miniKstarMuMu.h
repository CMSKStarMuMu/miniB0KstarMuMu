#ifndef MINIKSTARMUMU_H
#define MINIKSTARMUMU_H

// system include files
#include <memory>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TMath.h"
#include "TLorentzVector.h"

#include "Utils.h"
#include "B0KstMuMuTreeContent.h"
#include "miniB0KstarMuMu/miniKstarMuMu/interface/miniHLTObj.h"


//
// class declaration
//

class miniKstarMuMu : public edm::EDAnalyzer {
   public:
      explicit miniKstarMuMu (const edm::ParameterSet&);
      ~miniKstarMuMu();

      std::string getMuCat (pat::Muon const& muon);

   private:
      virtual void beginJob ();
      virtual void beginRun(const edm::Run &, const edm::EventSetup &);
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob ();
      void         MonteCarloStudies(  const edm::Event& );
      bool         isAncestor(const reco::Candidate*, const reco::Candidate * );
      bool         skipOscillations (const reco::GenParticle &, edm::Handle<reco::GenParticleCollection>);

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection>         vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection>            muonToken_;
      edm::EDGetTokenT<reco::TrackCollection>          trackToken_;
      edm::EDGetTokenT<reco::BeamSpot>                 beamSpotToken_;

      edm::EDGetTokenT<edm::TriggerResults>                    triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
      std::vector<std::string> trigTable_;
      std::vector<std::string> l1Table_;

      edm::EDGetTokenT<reco::GenParticleCollection>         prunedGenToken_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection>    packedGenToken_;

      edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;

      HLTPrescaleProvider hltPrescaleProvider_;
      edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_;
      edm::EDGetTokenT<GlobalExtBlkBxCollection> l1ext_;
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
      l1t::L1TGlobalUtil *fGtUtil;

};

#endif
