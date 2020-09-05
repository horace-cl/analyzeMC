#ifndef _Bsphi_PAT_h
#define _Bsphi_PAT_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"


//
// class decleration
//

class Bsphi_PAT : public edm::EDAnalyzer {
public:
  explicit Bsphi_PAT(const edm::ParameterSet&);
  ~Bsphi_PAT();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);

  //int const getMuCat(reco::Muon const& muon) const;
  //bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  //void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  double GetLifetime(TLorentzVector, TVector3, TVector3);
  bool hasFirstLayerPixelHits(const reco::TransientTrack& track);

  // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<std::vector<pat::GenericParticle>> trakCollection_label;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
 
  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;

  TTree*      tree_;
  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;

  std::vector<int>         *tri_JpsiTkTk, *tri_JpsiTk; 

  int                      muAcc, muTrig, weight;
  // *************************************
  unsigned int             nB;
  unsigned int             nMu;

  // this is the information for tracks after Bs vertex fit
  std::vector<float>       *J_px1, *J_py1, *J_pz1;
  std::vector<float>       *J_px2, *J_py2, *J_pz2;
  std::vector<float>       *B_k1_px, *B_k1_py, *B_k1_pz;
  std::vector<float>       *B_k2_px, *B_k2_py, *B_k2_pz;

  std::vector<int>         *pi1_trackerhits, *pi2_trackerhits;
  std::vector<int>         *pi1_pixelhits, *pi2_pixelhits;
  std::vector<float>       *pi1dz, *pi1dzE,*pi2dz, *pi2dzE;
  std::vector<float>       *pi1dxy, *pi1dxyE,*pi2dxy, *pi2dxyE;
  std::vector<float>       *pi1d0;

  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz;
  std::vector<float>       *B_phi_mass;
  std::vector<int>         *B_k_charge1, *B_k_charge2; 
  std::vector<float>       *B_k_px_track, *B_k_py_track, *B_k_pz_track;
  std::vector<float>       *B_k2_px_track, *B_k2_py_track, *B_k2_pz_track;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;

  std::vector<float>       *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  // vertice primario CON mayor Pt
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;

   // vertice primario CON mejor pointin-angle
   std::vector<float>          *pVtxIPX,  *pVtxIPY, *pVtxIPZ, *pVtxIPXE, *pVtxIPYE, *pVtxIPZE, *pVtxIPCL;
   std::vector<float>          *pVtxIPXYE,  *pVtxIPXZE, *pVtxIPYZE;

  // refitting the primary without the tracks in the B reco candidate
  std::vector<float>       *priRfVtxX, *priRfVtxY, *priRfVtxZ, *priRfVtxXE, *priRfVtxYE, *priRfVtxZE, *priRfVtxCL;
  std::vector<float>       *priRfVtxXYE, *priRfVtxXZE, *priRfVtxYZE;
  std::vector<int>         *priRfNTrkDif;
  
  // ********************************** ************************************************************************

  std::vector<float>       *B_Prob, *B_J_Prob;

  std::vector<float>       *B_DecayVtxX,  *B_DecayVtxY,  *B_DecayVtxZ;
  std::vector<double>      *B_DecayVtxXE, *B_DecayVtxYE, *B_DecayVtxZE;
  std::vector<double>      *B_DecayVtxXYE, *B_DecayVtxXZE, *B_DecayVtxYZE;
 
  std::vector<float>       *B_J_DecayVtxX,   *B_J_DecayVtxY,   *B_J_DecayVtxZ;
  std::vector<float>       *B_J_DecayVtxXE,  *B_J_DecayVtxYE,  *B_J_DecayVtxZE;
  std::vector<float>       *B_J_DecayVtxXYE, *B_J_DecayVtxXZE, *B_J_DecayVtxYZE;

  UInt_t trigger;

  int  run, event;
  int   lumiblock;

  TLorentzVector gen_b_p4,gen_phi_p4,gen_pion1_p4,gen_pion2_p4,gen_jpsi_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_b_vtx,gen_jpsi_vtx;
  float          gen_b_ct;


};
#endif