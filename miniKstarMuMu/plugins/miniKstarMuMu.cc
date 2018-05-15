// system include files
#include <memory>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "../interface/miniKstarMuMu.h"
#include "../src/B0KstMuMuTreeContent.cc"
#include "../src/Utils.cc"

#define TRKMAXR 110.0 // [cm]
#define TRKMAXZ 280.0 // [cm]

#define MUVARTOLE 0.01  // [From 0 to 1]
#define HADVARTOLE 0.10 // [From 0 to 1]



miniKstarMuMu::miniKstarMuMu(const edm::ParameterSet& iConfig):
    vtxToken_       (consumes<reco::VertexCollection>          (iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_      (consumes<pat::MuonCollection>             (iConfig.getParameter<edm::InputTag>("muons"))),
    trackToken_     (consumes<pat::PackedCandidateCollection>  (iConfig.getParameter<edm::InputTag>("tracks"))),
    beamSpotToken_  (consumes<reco::BeamSpot>                  (iConfig.getParameter<edm::InputTag>("beamSpot"))),

    // # Load HLT-trigger cuts #
    CLMUMUVTX          ( iConfig.getUntrackedParameter<double>("MuMuVtxCL")      ),
    LSMUMUBS           ( iConfig.getUntrackedParameter<double>("MuMuLsBS")       ),
    DCAMUMU            ( iConfig.getUntrackedParameter<double>("DCAMuMu")        ),
    DCAMUBS            ( iConfig.getUntrackedParameter<double>("DCAMuBS")        ),
    COSALPHAMUMUBS     ( iConfig.getUntrackedParameter<double>("cosAlphaMuMuBS") ),
    MUMINPT            ( iConfig.getUntrackedParameter<double>("MinMupT")        ),
    MUMAXETA           ( iConfig.getUntrackedParameter<double>("MuEta")          ),
    MINMUMUPT          ( iConfig.getUntrackedParameter<double>("MuMupT")         ),
    MINMUMUINVMASS     ( iConfig.getUntrackedParameter<double>("MinMuMuMass")    ),
    MAXMUMUINVMASS     ( iConfig.getUntrackedParameter<double>("MaxMuMuMass")    ),
   
    // # Load pre-selection cuts #
    B0MASSLOWLIMIT     ( iConfig.getUntrackedParameter<double>("MinB0Mass")      ),
    B0MASSUPLIMIT      ( iConfig.getUntrackedParameter<double>("MaxB0Mass")      ),
    CLB0VTX            ( iConfig.getUntrackedParameter<double>("B0VtxCL")        ),
    KSTMASSWINDOW      ( iConfig.getUntrackedParameter<double>("KstMass")        ),
    HADDCASBS          ( iConfig.getUntrackedParameter<double>("HadDCASBS")      ),
    MINHADPT           ( iConfig.getUntrackedParameter<double>("HadpT")          ),
    MAXB0PREMASS       ( iConfig.getUntrackedParameter<double>("MaxB0RoughMass") ),
  
    printMsg        ( iConfig.getUntrackedParameter<bool>("printMsg",false) ),    
    theTree(0)
{
  NTuple = new B0KstMuMuTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();
  
  Utility = new Utils();

}

void
miniKstarMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// Get magnetic field
	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
    // Get BeamSpot
    edm::Handle<reco::BeamSpot> beamSpotH;
    iEvent.getByToken(beamSpotToken_, beamSpotH);
    reco::BeamSpot beamSpot = *beamSpotH;

    reco::TrackRef muTrackm;
    reco::TrackRef muTrackp;
    const reco::Track *Trackm;
    const reco::Track *Trackp;
    double pT;
    double eta;
    double chi;
    double ndf;

    double LSBS;
    double LSBSErr;
    double cosAlphaBS;
    double cosAlphaBSErr;
    
    TrajectoryStateClosestToPoint theDCAXBS;
    ClosestApproachInRPhi ClosestApp;
    GlobalPoint XingPoint;

    float muonMassErr           = Utility->muonMassErr;
    float pionMassErr           = Utility->pionMassErr;
    float kaonMassErr           = Utility->kaonMassErr;
    const ParticleMass muonMass = Utility->muonMass;
    const ParticleMass pionMass = Utility->pionMass;
    const ParticleMass kaonMass = Utility->kaonMass;


    TLorentzVector  mum_lv;
    TLorentzVector  mup_lv;
    TLorentzVector jpsi_lv;
  
    TLorentzVector tkm_lv;
    TLorentzVector tkp_lv;
    TLorentzVector kst_lv;
    TLorentzVector kstbar_lv;

    KinematicParticleFactoryFromTransientTrack partFactory;
    AdaptiveVertexFitter theVtxFitter;                              // Vertex fitter in nominal reconstruction
    KinematicParticleVertexFitter PartVtxFitter;                    // Vertex fit with vtx constraint
  
    std::vector<float> mum_isovec, mup_isovec, trkm_isovec, trkp_isovec; 

    
//     std::cout <<iEvent.id() << std::endl;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    // Get PAT Tracks
    edm::Handle<pat::PackedCandidateCollection> tracks;
    iEvent.getByToken(trackToken_, tracks);

    for (const pat::Muon &mum : *muons) {
        muTrackm = mum.innerTrack();
        if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1)) continue;

        pT  = muTrackm -> pt();
        eta = muTrackm -> eta();

        if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
        {
          std::cout << __LINE__ << " : break --> too low pT of mu- : " << pT  << " or " << mum.pt() << " or too high eta : " << eta << std::endl;
          break;
        }

        const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle));
        theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
        if (theDCAXBS.isValid() == false)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
          continue;
        }
        double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
        double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
        if (fabs(DCAmumBS) > DCAMUBS)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
          continue;
        }

        for (const pat::Muon &mup : *muons) {

            muTrackp = mup.innerTrack();
            if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1)) continue;

            pT  = muTrackp->pt();
            eta = muTrackp->eta();
            if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
            {
              if (printMsg) std::cout << __LINE__ << " : break --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
              break;
            }

            const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));

            // ###############################
            // # Compute mu+ DCA to BeamSpot #
            // ###############################
            theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
            if (theDCAXBS.isValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
              continue;
            }
            double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
            double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
            if (fabs(DCAmupBS) > DCAMUBS)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS << std::endl;
              continue;
            }
          
            // ############################################
            // # Check goodness of muons closest approach #
            // ############################################
            ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
            if (ClosestApp.status() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
              continue;
            }
            XingPoint = ClosestApp.crossingPoint();
            if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
              continue;
            }
      
            double mumuDCA = ClosestApp.distance();
            if (mumuDCA > DCAMUMU)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA << std::endl;
              continue;
            }

            // ############################################
            // # Cut on the dimuon inviariant mass and pT #
            // ############################################
            jpsi_lv.SetPxPyPzE( muTrackmTT.track().px() + muTrackpTT.track().px(), 
                                muTrackmTT.track().py() + muTrackpTT.track().py(),
                                muTrackmTT.track().pz() + muTrackpTT.track().pz(),
                                sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ) + sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) )
                              );
            if ((jpsi_lv.Pt() < (MINMUMUPT*(1.0-MUVARTOLE))) || (jpsi_lv.M() < (MINMUMUINVMASS*(1.0-MUVARTOLE))) || (jpsi_lv.M() > (MAXMUMUINVMASS*(1.0+MUVARTOLE))))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << jpsi_lv.Pt() << "\tinv. mass: " << jpsi_lv.M() << std::endl;
              continue;
            }

            chi = 0.;
            ndf = 0.;
            // ####################################################
            // # Try to vertex the two muons to get dimuon vertex #
            // ####################################################
            std::vector<RefCountedKinematicParticle> muonParticles;
            muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi,ndf, muonMassErr));
            muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass,chi,ndf, muonMassErr));
        
            RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
            if (mumuVertexFitTree->isValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
              continue; 
            }
        
            mumuVertexFitTree->movePointerToTheTop();
            RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
            if (mumu_KV->vertexIsValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
              continue;
            }
            if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) < CLMUMUVTX)
              {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: " << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) << std::endl;
              continue;
            }
  
  
            // #########################################################
            // # Extract the re-fitted tracks after the dimuon vtx fit #
            // #########################################################
            RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
            mumuVertexFitTree->movePointerToTheTop();
  
            mumuVertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle refitMum  = mumuVertexFitTree->currentParticle();
            const reco::TransientTrack refitMumTT = refitMum->refittedTransientTrack();
  
            mumuVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle refitMup  = mumuVertexFitTree->currentParticle();
            const reco::TransientTrack refitMupTT = refitMup->refittedTransientTrack();
  
  
            // ########################
            // # Muon pT and eta cuts #
            // ########################
            pT  = refitMupTT.track().pt();
            eta = refitMupTT.track().eta();
            if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
              continue;
            }
  
            pT  = refitMumTT.track().pt();
            eta = refitMumTT.track().eta();
            if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> too low pT of mu- : " << pT << " or too high eta : " << eta << std::endl;
//               skip = true;
              continue;
            }

			// ############################################
			// # Cut on the dimuon invariant mass and pT #
			// ############################################
			pT = sqrt((refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) * (refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) +
					  (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()) * (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()));
			double MuMuInvMass = mumu_KP->currentState().mass();
			if ((pT < MINMUMUPT) || (MuMuInvMass < MINMUMUINVMASS) || (MuMuInvMass > MAXMUMUINVMASS))
			{
			  if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << pT << "\tinv. mass: " << MuMuInvMass << std::endl;
			  continue;
			}

            // ######################################################
            // # Compute the distance between mumu vtx and BeamSpot #
            // ######################################################
            double MuMuLSBS;
            double MuMuLSBSErr;
            Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
                        beamSpot.position().x(),beamSpot.position().y(),0.0,
                        mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
                        mumu_KV->error().matrix()(0,1),0.0,0.0,
                        beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
                        beamSpot.covariance()(0,1),0.0,0.0,
                        &MuMuLSBS,&MuMuLSBSErr);
            if (MuMuLSBS/MuMuLSBSErr < LSMUMUBS)     
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-" << MuMuLSBSErr << std::endl;
              continue;
            }
            // ###################################################################
            // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
            // ###################################################################
            double MuMuCosAlphaBS;
            double MuMuCosAlphaBSErr;
            Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),mumu_KP->currentState().globalMomentum().y(),0.0,
                                      mumu_KV->position().x() - beamSpot.position().x(),mumu_KV->position().y() - beamSpot.position().y(),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),mumu_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                      mumu_KV->error().cxx() + beamSpot.covariance()(0,0),mumu_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                      mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
            if (MuMuCosAlphaBS < COSALPHAMUMUBS)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS << "+/-" << MuMuCosAlphaBSErr << std::endl;
              continue;
            }
  
  
            // ########################### convert KinFit vertex to reco Vertex
            reco::Vertex::Point mumu_GP  = reco::Vertex::Point(mumu_KV->position().x(), mumu_KV->position().y(), mumu_KV->position().z());
            const reco::Vertex::Error mumu_error = mumu_KV->vertexState().error().matrix();
            float mumu_chi2      = mumu_KV -> chiSquared();
            float mumu_ndof      = mumu_KV -> degreesOfFreedom();
            reco::Vertex mumu_rv =  reco::Vertex( mumu_GP, mumu_error, mumu_chi2, mumu_ndof, 2 );

            for (const pat::PackedCandidate &tkm : *tracks) {
                Trackm = tkm.bestTrack();
                if (Trackm->charge() != -1) continue;
//                 if ((Trackm->isNull() 	== true) || (Trackm->charge() != -1)) continue;

                const reco::TransientTrack TrackmTT((*Trackm), &(*bFieldHandle));

                // ##########################
                // # Track- pT and eta cuts #
                // ##########################
                pT = TrackmTT.track().pt();
                if (pT < (MINHADPT*(1.0-HADVARTOLE)))
                {
                  if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
                  break;
                }
                
                // ######################################
                // # Compute K*0 track- DCA to BeamSpot #
                // ######################################
                theDCAXBS = TrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
                if (theDCAXBS.isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track-" << std::endl;
                  continue;
                }
                double DCAKstTrkmBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
                double DCAKstTrkmBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
                if (fabs(DCAKstTrkmBS/DCAKstTrkmBSErr) < HADDCASBS)
                {
//                   if (printMsg) std::cout << __LINE__ << " : continue --> track- DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkmBS << "+/-" << DCAKstTrkmBSErr << std::endl;
                  continue;
                }

				for (const pat::PackedCandidate &tkp : *tracks) {
					Trackp = tkp.bestTrack();
                    if (Trackp->charge() != +1) continue;

                    const reco::TransientTrack TrackpTT((*Trackp), &(*bFieldHandle));

                    pT = TrackpTT.track().pt();
                    if (pT < (MINHADPT*(1.0-HADVARTOLE)))
                    {
                      if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
                      break;
                    }
    
                    // ######################################
                    // # Compute K*0 track+ DCA to BeamSpot #
                    // ######################################
                    theDCAXBS = TrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
                    if (theDCAXBS.isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track+" << std::endl;
                      continue;
                    }
                    double DCAKstTrkpBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
                    double DCAKstTrkpBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
                    if (fabs(DCAKstTrkpBS/DCAKstTrkpBSErr) < HADDCASBS)
                    {
//                       if (printMsg) std::cout << __LINE__ << " : continue --> track+ DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkpBS << "+/-" << DCAKstTrkpBSErr << std::endl;
                      continue;
                    }
     
     
                    // ##############################################
                    // # Check goodness of hadrons closest approach #
                    // ##############################################
                    ClosestApp.calculate(TrackpTT.initialFreeState(),TrackmTT.initialFreeState());
                    if (ClosestApp.status() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
                      continue;
                    }
                    XingPoint = ClosestApp.crossingPoint();
                    if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
                      continue;
                    }



                    // ######################################################
                    // # Check if K*0 (OR K*0bar) mass is within acceptance #
                    // ######################################################
                    tkm_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                                      TrackmTT.track().momentum().y(),
                                      TrackmTT.track().momentum().z(),
                                      sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->pionMass,2) ));
                    tkp_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                                      TrackpTT.track().momentum().y(),
                                      TrackpTT.track().momentum().z(),
                                      sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->kaonMass,2) ) );
                    double kstInvMass = (tkm_lv + tkp_lv).M();
                    
                    tkm_lv.SetE(sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->kaonMass,2) ));
                    tkp_lv.SetE(sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->pionMass,2) ));
                    double kstBarInvMass = (tkm_lv + tkp_lv).M();
                    
                    if ((fabs(kstInvMass - Utility->kstMass)    > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))) &&
                        (fabs(kstBarInvMass - Utility->kstMass) > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))))
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
                      continue;
                    }
    
    
    
                    // ####################################################
                    // # @@@ Make K* and implement pre-selection cuts @@@ #
                    // ####################################################
                    if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ I have 2 good oppositely-charged tracks. I'm trying to vertex them @@@" << std::endl;
    
                    chi = 0.;
                    ndf = 0.;
                    // ##############################################################################
                    // # Try to vertex the two Tracks to get K*0 vertex: pion = track- | k = track+ #
                    // ##############################################################################
                    std::vector<RefCountedKinematicParticle> kstParticles;
                    kstParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
                    kstParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));
                    
                    RefCountedKinematicTree kstVertexFitTree = PartVtxFitter.fit(kstParticles);
                    if (kstVertexFitTree->isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
                      continue;
                    }
                    
                    kstVertexFitTree->movePointerToTheTop();
                    RefCountedKinematicParticle kst_KP = kstVertexFitTree->currentParticle();
                    RefCountedKinematicVertex kst_KV   = kstVertexFitTree->currentDecayVertex();
                    if (kst_KV->vertexIsValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
                      continue;
                    }
    
    
                    chi = 0.;
                    ndf = 0.;
                    // #################################################################################
                    // # Try to vertex the two Tracks to get K*0bar vertex: pion = track+ | k = track- #
                    // #################################################################################
                    std::vector<RefCountedKinematicParticle> kstBarParticles;
                    kstBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
                    kstBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
                    
                    RefCountedKinematicTree kstBarVertexFitTree = PartVtxFitter.fit(kstBarParticles);
                    if (kstBarVertexFitTree->isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
                      continue;
                    }
                    
                    kstBarVertexFitTree->movePointerToTheTop();
                    RefCountedKinematicParticle kstBar_KP = kstBarVertexFitTree->currentParticle();
                    RefCountedKinematicVertex kstBar_KV   = kstBarVertexFitTree->currentDecayVertex();
                    if (kstBar_KV->vertexIsValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
                      continue;
                    }
      
                    // ######################################################
                    // # Check if K*0 (OR K*0bar) mass is within acceptance #
                    // ######################################################
                    kstInvMass    = kst_KP->currentState().mass();
                    kstBarInvMass = kstBar_KP->currentState().mass();
                    if ((fabs(kstInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma) && (fabs(kstBarInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma))
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
                      continue;
                    }
      
      
                    // ###########################################################
                    // # Extract the re-fitted tracks after the dihadron vtx fit #
                    // ###########################################################
                    kstVertexFitTree->movePointerToTheTop();
      
                    kstVertexFitTree->movePointerToTheFirstChild();
                    RefCountedKinematicParticle refitTrkm  = kstVertexFitTree->currentParticle();
                    const reco::TransientTrack refitTrkmTT = refitTrkm->refittedTransientTrack();
      
                    kstVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle refitTrkp  = kstVertexFitTree->currentParticle();
                    const reco::TransientTrack refitTrkpTT = refitTrkp->refittedTransientTrack();
    
    
                    // ##########################
                    // # Hadron pT and eta cuts #
                    // ##########################
                    pT = refitTrkpTT.track().pt();
                    if (pT < MINHADPT)
                    {
                      if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
                      continue;
                    }
      
                    pT = refitTrkmTT.track().pt();
                    if (pT < MINHADPT)
                    {
                      if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
//                       skip = true;
                      continue;
                    }
    
                    // ####################################################
                    // # @@@ Make B0 and implement pre-selection cuts @@@ #
                    // ####################################################
                    if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ I have 4 good charged tracks. I'm trying to vertex them @@@" << std::endl;
      
                    TLorentzVector a_lv, b_lv, c_lv, d_lv, tot_lv;
                    a_lv.SetPxPyPzE(muTrackmTT.track().momentum().x(),
                                    muTrackmTT.track().momentum().y(),
                                    muTrackmTT.track().momentum().z(),
                                    sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ));
                    b_lv.SetPxPyPzE(muTrackpTT.track().momentum().x(),
                                    muTrackpTT.track().momentum().y(),
                                    muTrackpTT.track().momentum().z(),
                                    sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) ));
                    c_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                                    TrackmTT.track().momentum().y(),
                                    TrackmTT.track().momentum().z(),
                                    sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->pionMass,2) ));
                    d_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                                    TrackpTT.track().momentum().y(),
                                    TrackpTT.track().momentum().z(),
                                    sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->kaonMass,2) ));
                                
                    tot_lv = a_lv + b_lv + c_lv + d_lv;
                    if  (tot_lv.M() > MAXB0PREMASS) {
                      if (printMsg) std::cout << __LINE__ << " : continue --> b0 mass before fit is > max value" << std::endl;
                        continue;    
                    }         
  
  
  
                    // #################################################
                    // # Check if the hadron Track- is actually a muon #
                    // #################################################
//                     MuMCat.clear();
//                     MuMCat = "NotMatched";
//                     for (const pat::Muon &imutmp : *muons) {
//                     {
//                       muTrackTmp = imutmp.innerTrack();
//                       if ((muTrackTmp.isNull() == true) || (muTrackTmp->charge() != -1)) continue;
//                       if (Trackm == muTrackTmp)
//                       {
//                         MuMCat.clear();
//                         MuMCat.append(getMuCat(*iMuon));
//                         if (printMsg) std::cout << __LINE__ << " : negative charged hadron is actually a muon (momentum: " << Trackm->p() << ") whose category is: " << MuMCat.c_str() << std::endl;
//                         break;
//                       }
//                     }
                    // same for track+
                    
 
 
                    chi = 0.;
                    ndf = 0.;
                    // #####################################
                    // # B0 vertex fit with vtx constraint #
                    // #####################################
                    std::vector<RefCountedKinematicParticle> bParticles;
                    bParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
                    bParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
                    bParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
                    bParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));
      
                    RefCountedKinematicTree bVertexFitTree = PartVtxFitter.fit(bParticles);
                    if (bVertexFitTree->isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
                      continue;
                    }
      
                    bVertexFitTree->movePointerToTheTop();
                    RefCountedKinematicParticle b_KP = bVertexFitTree->currentParticle();
                    RefCountedKinematicVertex b_KV   = bVertexFitTree->currentDecayVertex();
                    if (b_KV->vertexIsValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
                      continue;
                    }
      
      
                    chi = 0.;
                    ndf = 0.;
                    // ########################################
                    // # B0bar vertex fit with vtx constraint #
                    // ########################################
                    std::vector<RefCountedKinematicParticle> bBarParticles;
                    bBarParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
                    bBarParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
                    bBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
                    bBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
                        
                    RefCountedKinematicTree bBarVertexFitTree = PartVtxFitter.fit(bBarParticles);
                    if (bBarVertexFitTree->isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
                      continue;
                    }
      
                    bBarVertexFitTree->movePointerToTheTop();
                    RefCountedKinematicParticle bBar_KP = bBarVertexFitTree->currentParticle();
                    RefCountedKinematicVertex bBar_KV   = bBarVertexFitTree->currentDecayVertex();
                    if (bBar_KV->vertexIsValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
                      continue;
                    }
      
                    // ########################
                    // # Cuts on B0 AND B0bar #
                    // ########################
                    if (((b_KP->currentState().mass()    < B0MASSLOWLIMIT) || (b_KP->currentState().mass()    > B0MASSUPLIMIT)) &&
                        ((bBar_KP->currentState().mass() < B0MASSLOWLIMIT) || (bBar_KP->currentState().mass() > B0MASSUPLIMIT)))
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> bad B0 mass: " << b_KP->currentState().mass() << " AND B0bar mass: " << bBar_KP->currentState().mass() << std::endl;
                      continue;
                    }
                    if ((TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))) < CLB0VTX) &&
                        (TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) < CLB0VTX))
                    {
                      if (printMsg)
                      {
                        std::cout << __LINE__ << " : continue --> bad vtx CL from B0 fit: " << TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom())));
                        std::cout << " AND bad vtx CL from B0bar fit: " << TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) << std::endl;
                      }
                      continue;
                    }
    
    
                    // ###############################
                    // # Cuts on B0 L/sigma BeamSpot #
                    // ###############################
                    Utility->computeLS (b_KV->position().x(),b_KV->position().y(),0.0,
                                        beamSpot.position().x(),beamSpot.position().y(),0.0,
                                        b_KV->error().cxx(),b_KV->error().cyy(),0.0,
                                        b_KV->error().matrix()(0,1),0.0,0.0,
                                        beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
                                        beamSpot.covariance()(0,1),0.0,0.0,
                                        &LSBS,&LSBSErr);
      
      
                    // ##############################
                    // # Compute B0 DCA to BeamSpot #
                    // ##############################
                    theDCAXBS = b_KP->refittedTransientTrack().trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
                    if (theDCAXBS.isValid() == false)
                    {
                      if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for B0" << std::endl;
                      continue;
                    }      
                    double DCAB0BS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
                    double DCAB0BSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
      
              
                    // #####################################
                    // # Compute B0 cos(alpha) to BeamSpot #
                    // #####################################
                    Utility->computeCosAlpha (b_KP->currentState().globalMomentum().x(),b_KP->currentState().globalMomentum().y(),0.0,
                                              b_KV->position().x() - beamSpot.position().x(),b_KV->position().y() - beamSpot.position().y(),0.0,
                                              b_KP->currentState().kinematicParametersError().matrix()(3,3),b_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                              b_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                              b_KV->error().cxx() + beamSpot.covariance()(0,0),b_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                              b_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                              &cosAlphaBS,&cosAlphaBSErr);






                    // #######################################
                    // # @@@ Fill B0-candidate variables @@@ #
                    // #######################################
                    if (printMsg) std::cout << __LINE__ << " : @@@ Filling B0 candidate variables @@@\n\n" << std::endl;
      
      
                    // ############
                    // # Save: B0 #
                    // ############
                    NTuple->nB++;
      
                    NTuple->bMass->push_back(b_KP->currentState().mass());
                    NTuple->bMassE->push_back(sqrt(b_KP->currentState().kinematicParametersError().matrix()(6,6)));
                    NTuple->bBarMass->push_back(bBar_KP->currentState().mass());
                    NTuple->bBarMassE->push_back(sqrt(bBar_KP->currentState().kinematicParametersError().matrix()(6,6)));
      
                    NTuple->bPx->push_back(b_KP->currentState().globalMomentum().x());
                    NTuple->bPy->push_back(b_KP->currentState().globalMomentum().y());
                    NTuple->bPz->push_back(b_KP->currentState().globalMomentum().z());          
      
                    NTuple->bVtxCL->push_back(TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))));
                    NTuple->bVtxX->push_back(b_KV->position().x());
                    NTuple->bVtxY->push_back(b_KV->position().y());
                    NTuple->bVtxZ->push_back(b_KV->position().z());
          
                    NTuple->bCosAlphaBS->push_back(cosAlphaBS);
                    NTuple->bCosAlphaBSE->push_back(cosAlphaBSErr);
      
                    NTuple->bLBS->push_back(LSBS);
                    NTuple->bLBSE->push_back(LSBSErr);
      
                    NTuple->bDCABS->push_back(DCAB0BS);
                    NTuple->bDCABSE->push_back(DCAB0BSErr);
             
      
                    // #############
                    // # Save: K*0 #
                    // #############
                    NTuple->kstMass->push_back(kstInvMass);
                    NTuple->kstMassE->push_back(sqrt(kst_KP->currentState().kinematicParametersError().matrix()(6,6)));
                    NTuple->kstBarMass->push_back(kstBarInvMass);
                    NTuple->kstBarMassE->push_back(sqrt(kstBar_KP->currentState().kinematicParametersError().matrix()(6,6)));
      
                    NTuple->kstPx->push_back(kst_KP->currentState().globalMomentum().x());
                    NTuple->kstPy->push_back(kst_KP->currentState().globalMomentum().y());
                    NTuple->kstPz->push_back(kst_KP->currentState().globalMomentum().z());
      
                    NTuple->kstVtxCL->push_back(TMath::Prob(static_cast<double>(kst_KV->chiSquared()), static_cast<int>(rint(kst_KV->degreesOfFreedom()))));
                    NTuple->kstVtxX->push_back(kst_KV->position().x());
                    NTuple->kstVtxY->push_back(kst_KV->position().y());
                    NTuple->kstVtxZ->push_back(kst_KV->position().z());
      
      
                    // #################
                    // # Save: mu+ mu- #
                    // #################
                    NTuple->mumuMass->push_back(MuMuInvMass);
                    NTuple->mumuMassE->push_back(sqrt(mumu_KP->currentState().kinematicParametersError().matrix()(6,6)));
      
                    NTuple->mumuPx->push_back(mumu_KP->currentState().globalMomentum().x());
                    NTuple->mumuPy->push_back(mumu_KP->currentState().globalMomentum().y());
                    NTuple->mumuPz->push_back(mumu_KP->currentState().globalMomentum().z());
      
                    NTuple->mumuVtxCL->push_back(TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
                    NTuple->mumuVtxX->push_back(mumu_KV->position().x());
                    NTuple->mumuVtxY->push_back(mumu_KV->position().y());
                    NTuple->mumuVtxZ->push_back(mumu_KV->position().z());
      
                    NTuple->mumuCosAlphaBS->push_back(MuMuCosAlphaBS);
                    NTuple->mumuCosAlphaBSE->push_back(MuMuCosAlphaBSErr);
                    NTuple->mumuLBS->push_back(MuMuLSBS);
                    NTuple->mumuLBSE->push_back(MuMuLSBSErr);
                    NTuple->mumuDCA->push_back(mumuDCA);
      
      
                    // #############
                    // # Save: mu- #
                    // #############
                    NTuple->mumHighPurity->push_back( (int)muTrackm->quality(reco::Track::highPurity));
                    NTuple->mumCL->push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
                    NTuple->mumNormChi2->push_back(muTrackm->normalizedChi2());
                    NTuple->mumPx->push_back(refitMumTT.track().momentum().x());
                    NTuple->mumPy->push_back(refitMumTT.track().momentum().y());
                    NTuple->mumPz->push_back(refitMumTT.track().momentum().z());
      
//                     NTuple->mumDCAVtx->push_back(DCAmumVtx);
//                     NTuple->mumDCAVtxE->push_back(DCAmumVtxErr);
                    NTuple->mumDCABS->push_back(DCAmumBS);
                    NTuple->mumDCABSE->push_back(DCAmumBSErr);
      
//                     NTuple->mumKinkChi2->push_back(iMuonM->combinedQuality().trkKink);
                    NTuple->mumFracHits->push_back(static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
                                                                       muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                       muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                                                       muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
//                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
//                     NTuple->mumdxyVtx->push_back(theDCAXVtx.second.value());
//                     NTuple->mumdzVtx->push_back(muTrackmTT.track().dz(bestVtxReFit.position()));
      
                    NTuple->mumCat->push_back(getMuCat(mum));
      
                    NTuple->mumNPixHits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
                    NTuple->mumNPixLayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());  
                    NTuple->mumNTrkHits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
                    NTuple->mumNTrkLayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
//                     if (mum->isGlobalMuon() == true) NTuple->mumNMuonHits->push_back(iMuonM->globalTrack()->hitPattern().numberOfValidMuonHits());
//                     else NTuple->mumNMuonHits->push_back(0);
                    NTuple->mumNMatchStation->push_back(mum.numberOfMatchedStations());
      
      
      
                    // #############
                    // # Save: mu+ #
                    // #############
                    NTuple->mupHighPurity->push_back( (int) muTrackp->quality(reco::Track::highPurity));
                    NTuple->mupCL->push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
                    NTuple->mupNormChi2->push_back(muTrackp->normalizedChi2());
                    NTuple->mupPx->push_back(refitMupTT.track().momentum().x());
                    NTuple->mupPy->push_back(refitMupTT.track().momentum().y());
                    NTuple->mupPz->push_back(refitMupTT.track().momentum().z());
      
                    NTuple->mupDCABS->push_back(DCAmupBS);
                    NTuple->mupDCABSE->push_back(DCAmupBSErr);
                    
//                     NTuple->mupKinkChi2->push_back(iMuonP->combinedQuality().trkKink);
                    NTuple->mupFracHits->push_back(static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackp->hitPattern().numberOfValidHits() +
                                                                       muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                       muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                                                       muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
//                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
//                     NTuple->mupdxyVtx->push_back(theDCAXVtx.second.value());
//                     NTuple->mupdzVtx->push_back(muTrackpTT.track().dz(bestVtxReFit.position()));
      
                    NTuple->mupCat->push_back(getMuCat(mup));
      
                    NTuple->mupNPixHits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
                    NTuple->mupNPixLayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());  
                    NTuple->mupNTrkHits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
                    NTuple->mupNTrkLayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
                    if (mup.isGlobalMuon() == true) NTuple->mupNMuonHits->push_back(mup.globalTrack()->hitPattern().numberOfValidMuonHits());
                    else NTuple->mupNMuonHits->push_back(0);
                    NTuple->mupNMatchStation->push_back(mup.numberOfMatchedStations());
      
      
                    // ################
                    // # Save: Track- #
                    // ################
                    NTuple->kstTrkmHighPurity->push_back( (int)Trackm->quality(reco::Track::highPurity));
                    NTuple->kstTrkmCL->push_back(TMath::Prob(TrackmTT.chi2(), static_cast<int>(rint(TrackmTT.ndof()))));
                    NTuple->kstTrkmNormChi2->push_back(Trackm->normalizedChi2());
                    NTuple->kstTrkmPx->push_back(refitTrkmTT.track().momentum().x());
                    NTuple->kstTrkmPy->push_back(refitTrkmTT.track().momentum().y());
                    NTuple->kstTrkmPz->push_back(refitTrkmTT.track().momentum().z());
      
//                     NTuple->kstTrkmDCAVtx->push_back(DCAKstTrkmVtx);
//                     NTuple->kstTrkmDCAVtxE->push_back(DCAKstTrkmVtxErr);
                    NTuple->kstTrkmDCABS->push_back(DCAKstTrkmBS);
                    NTuple->kstTrkmDCABSE->push_back(DCAKstTrkmBSErr);
      
                    NTuple->kstTrkmFracHits->push_back(static_cast<double>(Trackm->hitPattern().numberOfValidHits()) / static_cast<double>(Trackm->hitPattern().numberOfValidHits() +
                                                                           Trackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                           Trackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
//                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackmTT, bestVtxReFit);
//                     NTuple->kstTrkmdxyVtx->push_back(theDCAXVtx.second.value());
//                     NTuple->kstTrkmdzVtx->push_back(TrackmTT.track().dz(bestVtxReFit.position()));
      
//                     NTuple->kstTrkmMuMatch->push_back(MuMCat);
      
                    // I do NOT include the number of missing outer hits because the hadron might interact
                    NTuple->kstTrkmNPixHits->push_back(Trackm->hitPattern().numberOfValidPixelHits());
                    NTuple->kstTrkmNPixLayers->push_back(Trackm->hitPattern().pixelLayersWithMeasurement());  
                    NTuple->kstTrkmNTrkHits->push_back(Trackm->hitPattern().numberOfValidTrackerHits());
                    NTuple->kstTrkmNTrkLayers->push_back(Trackm->hitPattern().trackerLayersWithMeasurement());
      
      
                    // ################
                    // # Save: Track+ #
                    // ################
                    NTuple->kstTrkpHighPurity->push_back((int)Trackp->quality(reco::Track::highPurity));
                    NTuple->kstTrkpCL->push_back(TMath::Prob(TrackpTT.chi2(), static_cast<int>(rint(TrackpTT.ndof()))));
                    NTuple->kstTrkpNormChi2->push_back(Trackp->normalizedChi2());
                    NTuple->kstTrkpPx->push_back(refitTrkpTT.track().momentum().x());
                    NTuple->kstTrkpPy->push_back(refitTrkpTT.track().momentum().y());
                    NTuple->kstTrkpPz->push_back(refitTrkpTT.track().momentum().z());
      
//                     NTuple->kstTrkpDCAVtx->push_back(DCAKstTrkpVtx);
//                     NTuple->kstTrkpDCAVtxE->push_back(DCAKstTrkpVtxErr);
                    NTuple->kstTrkpDCABS->push_back(DCAKstTrkpBS);
                    NTuple->kstTrkpDCABSE->push_back(DCAKstTrkpBSErr);
      
                    NTuple->kstTrkpFracHits->push_back(static_cast<double>(Trackp->hitPattern().numberOfValidHits()) / static_cast<double>(Trackp->hitPattern().numberOfValidHits() +
                                                                           Trackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                           Trackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
//                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackpTT, bestVtxReFit);
//                     NTuple->kstTrkpdxyVtx->push_back(theDCAXVtx.second.value());
//                     NTuple->kstTrkpdzVtx->push_back(TrackpTT.track().dz(bestVtxReFit.position()));
      
//                     NTuple->kstTrkpMuMatch->push_back(MuPCat);
      
                    // I do NOT include the number of missing outer hits because the hadron might interact
                    NTuple->kstTrkpNPixHits->push_back(Trackp->hitPattern().numberOfValidPixelHits());
                    NTuple->kstTrkpNPixLayers->push_back(Trackp->hitPattern().pixelLayersWithMeasurement());  
                    NTuple->kstTrkpNTrkHits->push_back(Trackp->hitPattern().numberOfValidTrackerHits());
                    NTuple->kstTrkpNTrkLayers->push_back(Trackp->hitPattern().trackerLayersWithMeasurement());
      
/*      
                    // Save trigger matching for the 4 tracks
                    tmpString1.clear(); tmpString2.clear(); tmpString3.clear(); tmpString4.clear();
                    const pat::Muon* patMuonM = &(*iMuonM);
                    const pat::Muon* patMuonP = &(*iMuonP);
                    const pat::GenericParticle* patTrkm = &(*iTrackM);
                    const pat::GenericParticle* patTrkp = &(*iTrackP);
                    for (unsigned int i = 0; i < TrigTable_.size(); i++)
                      {
                        myString.clear(); myString.str(""); myString << TrigTable_[i].c_str() << "*";
                        if (patMuonM->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString1.append(TrigTable_[i]+" ");
                        if (patMuonP->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString2.append(TrigTable_[i]+" ");
                        if (patTrkm ->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString3.append(TrigTable_[i]+" ");
                        if (patTrkp ->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString4.append(TrigTable_[i]+" ");
                      }
                    if (tmpString1.size() == 0) tmpString1.append("NotInTable");
                    if (tmpString2.size() == 0) tmpString2.append("NotInTable");
                    if (tmpString3.size() == 0) tmpString3.append("NotInTable");
                    if (tmpString4.size() == 0) tmpString4.append("NotInTable");
                    NTuple->mumTrig    ->push_back(tmpString1);
                    NTuple->mupTrig    ->push_back(tmpString2);
                    NTuple->kstTrkmTrig->push_back(tmpString3);
                    NTuple->kstTrkpTrig->push_back(tmpString4);
      
                    // save minimum IP from any PV (2D and 3D)
                    mumMind0  = mupMind0 = 100;
                    TrkmMind0 = 100;
                    TrkpMind0 = 100;
                    mumMinIP  = 100;
                    mupMinIP  = 100;
                    TrkmMinIP = 100;
                    TrkpMinIP = 100;
                    std::pair<double,double>  IPPair;
    
                    for (std::vector<reco::Vertex>::const_iterator ipv = recVtx->begin(); ipv != recVtx->end(); ipv++) { 
                        vert = GlobalPoint(ipv->x(), ipv->y(), ipv->z());
    
                        traj = muTrackmTT.trajectoryStateClosestToPoint(vert );
                        if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mumMind0){
                          mumMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                          mumMind0E = traj.perigeeError().transverseImpactParameterError();
                        }  
    
                        traj = muTrackpTT.trajectoryStateClosestToPoint(vert );
                        if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mupMind0){
                          mupMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                          mupMind0E = traj.perigeeError().transverseImpactParameterError();
                        }  
    
                        traj = TrackmTT.trajectoryStateClosestToPoint(vert );
                        if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkmMind0){
                          TrkmMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                          TrkmMind0E = traj.perigeeError().transverseImpactParameterError();
                        }  
    
                        traj = TrackpTT.trajectoryStateClosestToPoint(vert );
                        if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkpMind0){
                          TrkpMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                          TrkpMind0E = traj.perigeeError().transverseImpactParameterError();
                        }  
                        
                        IPPair = pionImpactParameter(muTrackmTT,*ipv);
                        if (IPPair.first < mumMinIP){
                          mumMinIP  = IPPair.first;
                          mumMinIPE = IPPair.second;
                        }
                        IPPair = pionImpactParameter(muTrackpTT,*ipv);
                        if (IPPair.first < mupMinIP){
                          mupMinIP  = IPPair.first;
                          mupMinIPE = IPPair.second;
                        }
                        IPPair = pionImpactParameter(TrackmTT,*ipv);
                        if (IPPair.first < TrkmMinIP){
                          TrkmMinIP  = IPPair.first;
                          TrkmMinIPE = IPPair.second;
                        }
                        IPPair = pionImpactParameter(TrackpTT,*ipv);
                        if (IPPair.first < TrkpMinIP){
                          TrkpMinIP  = IPPair.first;
                          TrkpMinIPE = IPPair.second;
                        }
                    }
    
                    // isolation
                    double iso;
                    chi = 0; ndf = 0;
                    mum_isovec.clear(); mup_isovec.clear(); trkm_isovec.clear(); trkp_isovec.clear();
                    float bestVtx_x = bestVtx.x();
                    float bestVtx_y = bestVtx.y();
                    float bestVtx_z = bestVtx.z();
    
                     for (std::vector<pat::GenericParticle>::const_iterator iTrackIso = thePATTrackHandle->begin(); iTrackIso != thePATTrackHandle->end(); iTrackIso++)
                     {
                       if (iTrackIso == iTrackM || iTrackIso == iTrackP) continue;
                       if (muTrackm == iTrackIso->track() || muTrackp == iTrackIso->track()) continue;
                      const reco::TransientTrack TrackIsoTT(iTrackIso->track(), &(*bFieldHandle));
                      
                      iso = calculateIsolation(ClosestApp, muTrackmTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                      if (iso > 0 ) mum_isovec.push_back(iso);
    
                      iso = calculateIsolation(ClosestApp, muTrackpTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                      if (iso > 0 ) mup_isovec.push_back(iso);
    
                      iso = calculateIsolation(ClosestApp, TrackmTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                      if (iso > 0 ) trkm_isovec.push_back(iso);
    
                      iso = calculateIsolation(ClosestApp, TrackpTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                      if (iso > 0 ) trkp_isovec.push_back(iso);
                     } 
    
                    NTuple->mumIso      -> push_back(mum_isovec);
                    NTuple->mupIso      -> push_back(mup_isovec);
                    NTuple->kstTrkmIso  -> push_back(trkm_isovec);
                    NTuple->kstTrkpIso  -> push_back(trkp_isovec);
    
    
                    NTuple->mumMinIP2D      -> push_back(mumMind0);
                    NTuple->mumMinIP2DE     -> push_back(mumMind0E);
                    NTuple->mupMinIP2D      -> push_back(mupMind0);
                    NTuple->mupMinIP2DE     -> push_back(mupMind0E);
                    NTuple->kstTrkmMinIP2D  -> push_back(TrkmMind0);
                    NTuple->kstTrkmMinIP2DE -> push_back(TrkmMind0E);
                    NTuple->kstTrkpMinIP2D  -> push_back(TrkpMind0);
                    NTuple->kstTrkpMinIP2DE -> push_back(TrkpMind0E);
    
                    NTuple->mumMinIP      -> push_back(mumMinIP);
                    NTuple->mumMinIPE     -> push_back(mumMinIPE);
                    NTuple->mupMinIP      -> push_back(mupMinIP);
                    NTuple->mupMinIPE     -> push_back(mupMinIPE);
                    NTuple->kstTrkmMinIP  -> push_back(TrkmMinIP);
                    NTuple->kstTrkmMinIPE -> push_back(TrkmMinIPE);
                    NTuple->kstTrkpMinIP  -> push_back(TrkpMinIP);
                    NTuple->kstTrkpMinIPE -> push_back(TrkpMinIPE);
     
*/      
                    // #####################
                    // # Clear all vectors #
                    // #####################
//                     vertexTracks.clear();
                    bParticles.clear();
                    bBarParticles.clear();
                    kstParticles.clear();
                    kstBarParticles.clear();


                } // end for trackp
            }    
        }
    }


//     edm::Handle jets;
//     iEvent.getByToken(jetToken_, jets);
//     int ijet = 0;
//     for (const pat::Jet &j : *jets) {
//         if (j.pt() < 20) continue;
//         printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
//             j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
//         if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
//             std::vector daus(j.daughterPtrVector());
//             std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
//             for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
//                 const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
//                 printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
//             }
//         }
//     }
// 
//     edm::Handle fatjets;
//     iEvent.getByToken(fatjetToken_, fatjets);
//     for (const pat::Jet &j : *fatjets) {
//         printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f softdrop, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
//             j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSSoftDropMass"), j.userFloat("ak8PFJetsCHSPrunedMass"), j.userFloat("ak8PFJetsCHSTrimmedMass"), j.userFloat("ak8PFJetsCHSFilteredMass"), j.userFloat("cmsTopTagPFJetsCHSMassAK8"));
// 
//         // To get the constituents of the AK8 jets, you have to loop over all of the
//         // daughters recursively. To save space, the first two constituents are actually
//         // the Soft Drop SUBJETS, which will then point to their daughters.
//         // The remaining constituents are those constituents removed by soft drop but
//         // still in the AK8 jet.
//    std::vector constituents;
//         for ( unsigned ida = 0; ida < j.numberOfDaughters(); ++ida ) {
//      reco::Candidate const * cand = j.daughter(ida);
//      if ( cand->numberOfDaughters() == 0 )
//        constituents.push_back( cand ) ;
//      else {
//        for ( unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda ) {
//          reco::Candidate const * cand2 = cand->daughter(jda);
//          constituents.push_back( cand2 );
//        }
//      }
//    }
//    std::sort( constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda){return ida->pt() > jda->pt();} );
// 
//    for ( unsigned int ida = 0; ida < constituents.size(); ++ida ) {
//      const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ida]);
//      printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", ida,cand.pt(),cand.dz(PV.position()),cand.pdgId());
//    }
// 
//     }
// 

    theTree->Fill();
    NTuple->ClearNTuple();

}


std::string miniKstarMuMu::getMuCat (reco::Muon const& muon)
{
  std::stringstream muCat;
  muCat.str("");

  if (muon.isGlobalMuon() == true)
  {
    muCat << " GlobalMuon";
    if (muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) == true) muCat << " GlobalMuonPromptTight";
  }
  if (muon.isTrackerMuon() == true)
  {
    muCat << " TrackerMuon";
    if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) == true) muCat << " TrackerMuonArbitrated";
    if (muon::isGoodMuon(muon, muon::TMLastStationTight)     == true) muCat << " TMLastStationTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationLoose)     == true) muCat << " TMLastStationLoose";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) == true) muCat << " TM2DCompatibilityTight";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) == true) muCat << " TM2DCompatibilityLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationTight)      == true) muCat << " TMOneStationTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationLoose)      == true) muCat << " TMOneStationLoose";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)  == true) muCat << " TMLastStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngLoose)  == true) muCat << " TMLastStationAngLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngTight)   == true) muCat << " TMOneStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngLoose)   == true) muCat << " TMOneStationAngLoose";
  }
  if (muon.isStandAloneMuon() == true) muCat << " StandAloneMuon";
  if (muon.isCaloMuon()       == true) muCat << " CaloMuon";
  if ((muon.isGlobalMuon() == false) && (muon.isTrackerMuon() == false) && (muon.isStandAloneMuon() == false) && (muon.isCaloMuon() == false)) muCat << " NotInTable";

  return muCat.str();
}


void miniKstarMuMu::beginJob ()
{
    edm::Service<TFileService> outfile_;
    theTree = outfile_->make<TTree>("B0KstMuMuNTuple","B0KstMuMuNTuple");
    NTuple->MakeTreeBranches(theTree);
}
void miniKstarMuMu::endJob ()
{
    theTree->GetDirectory()->cd();
    theTree->Write();
}

miniKstarMuMu::~miniKstarMuMu()
{
  delete NTuple;
  delete Utility;

}



//define this as a plug-in
DEFINE_FWK_MODULE(miniKstarMuMu);