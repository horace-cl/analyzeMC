#include <memory>
#include "Bsphi_PAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

typedef math::Error<3>::type CovarianceMatrix;

Bsphi_PAT::Bsphi_PAT(const edm::ParameterSet& iConfig):
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),

  tree_(0), 

  tri_JpsiTkTk(0), tri_JpsiTk(0),
 
  nB(0), nMu(0),

  // this is the information for tracks after Lb mass contrain
  J_px1(0), J_py1(0), J_pz1(0),
  J_px2(0), J_py2(0), J_pz2(0),
  B_k1_px(0), B_k1_py(0), B_k1_pz(0),
  B_k2_px(0), B_k2_py(0), B_k2_pz(0),
 
  pi1_trackerhits(0), pi2_trackerhits(0),
  pi1_pixelhits(0), pi2_pixelhits(0),
  pi1dz(0), pi1dzE(0), pi2dz(0), pi2dzE(0),
  pi1dxy(0), pi1dxyE(0), pi2dxy(0), pi2dxyE(0),
  pi1d0(0),

  B_mass(0), B_px(0), B_py(0), B_pz(0),
  B_phi_mass(0),
  B_k_charge1(0), B_k_charge2(0),
  B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  B_k2_px_track(0), B_k2_py_track(0), B_k2_pz_track(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),

  B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),
  
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  pVtxIPX(0),   pVtxIPY(0),   pVtxIPZ(0), pVtxIPXE(0),   pVtxIPYE(0),   pVtxIPZE(0), pVtxIPCL(0),
  pVtxIPXYE(0),   pVtxIPXZE(0),   pVtxIPYZE(0),

  priRfVtxX(0), priRfVtxY(0), priRfVtxZ(0), priRfVtxXE(0), priRfVtxYE(0), priRfVtxZE(0), priRfVtxCL(0),
  priRfVtxXYE(0), priRfVtxXZE(0), priRfVtxYZE(0),
  priRfNTrkDif(0),
 
  B_Prob(0), B_J_Prob(0), 
 
  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  B_J_DecayVtxX(0),   B_J_DecayVtxY(0),   B_J_DecayVtxZ(0),
  B_J_DecayVtxXE(0),  B_J_DecayVtxYE(0),  B_J_DecayVtxZE(0),
  B_J_DecayVtxXYE(0), B_J_DecayVtxXZE(0), B_J_DecayVtxYZE(0),

  trigger(0),
 
  run(0), event(0),
  lumiblock(0)

{}


Bsphi_PAT::~Bsphi_PAT(){}

void Bsphi_PAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
 using std::vector;
 using namespace edm;
 using namespace reco;
 using namespace std;
 
 // Get event content information
 
 edm::ESHandle<TransientTrackBuilder> theB; 
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
 
 edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle;
 iEvent.getByToken(trakCollection_label,thePATTrackHandle);
 
 edm::Handle< View<pat::Muon> > thePATMuonHandle;
 iEvent.getByToken(dimuon_Label,thePATMuonHandle);
 
 edm::Handle<edm::TriggerResults> triggerResults_handle;
 iEvent.getByToken(triggerResults_Label, triggerResults_handle);
 
 edm::Handle<reco::GenParticleCollection> pruned;
 iEvent.getByToken(genCands_, pruned);
 
 edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

 lumiblock = iEvent.id().luminosityBlock();
 run = iEvent.id().run();
 event = iEvent.id().event();

 gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.); 
 gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
 gen_b_vtx.SetXYZ(0.,0.,0.);
 gen_jpsi_vtx.SetXYZ(0.,0.,0.);
 gen_b_ct = -9999.;

 if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 531) ) { //&& (dau->status() == 2) ) {
            foundit++;
            gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
            gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
            //int npion=0;
            for (size_t k=0; k<dau->numberOfDaughters(); k++) {
              const reco::Candidate *gdau = dau->daughter(k);
              if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
                foundit++;
                gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
                gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
                int nm=0;
                for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
                  const reco::Candidate *mm = gdau->daughter(l);
                  if (mm->pdgId()==13) { foundit++;
                     if (mm->status()!=1) {
                        for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                           const reco::Candidate *mu = mm->daughter(m);
                           if (mu->pdgId()==13 ) { //&& mu->status()==1) {
                              nm++;
                              gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                              break;
                           }
                        }
                     } else {
                       gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                       nm++;
                     }
                  }
                  if (mm->pdgId()==-13) { foundit++;
                     if (mm->status()!=1) {
                        for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                           const reco::Candidate *mu = mm->daughter(m);
                           if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
                              nm++;
                              gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                              break;
                           }
                        }
                     } else {
                       gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                       nm++;
                     }
                  }
                }
                if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                else foundit-=nm;
              }
	      if (gdau->pdgId()==333 ) {// pdgi for phi(1020)=333
		foundit++;
		gen_phi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
		for (size_t lk=0; lk<gdau->numberOfDaughters(); lk++) {
                  const reco::Candidate *kk = gdau->daughter(lk);
		  if (kk->pdgId()==321) { foundit++;
		    gen_pion1_p4.SetPtEtaPhiM(kk->pt(),kk->eta(),kk->phi(),kk->mass());
		  }
		  if (kk->pdgId()==-321) { foundit++;
		    gen_pion2_p4.SetPtEtaPhiM(kk->pt(),kk->eta(),kk->phi(),kk->mass());
		  }		  
		}//for lk
	      }
            } // for (size_t k
      }   // if (abs(dau->pdgId())==531 )
      if (foundit>=7) break;
    } // for i
    if (foundit!=7) {
       gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
       gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
       gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
       gen_b_vtx.SetXYZ(0.,0.,0.);
       gen_jpsi_vtx.SetXYZ(0.,0.,0.);
       gen_b_ct = -9999.;
       std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
 }
 
 nB = 0; nMu = 0;
 trigger = 0;

  if ( OnlyGen_ ) { 
    tree_->Fill();
    return;
 }

 
 iEvent.getByToken(triggerObjects_, triggerObjects);
 //std::stringstream ss0;
 std::string ss0 = "HLT_DoubleMu4_JpsiTrk_Displaced";
 std::string cc0 = "hltJpsiTkAllConeTracksIter";
 std::vector<float> obj_eta, obj_phi;

 if ( triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      unsigned int NTRIGGERS = 4;
       // para el 2017-2018
      std::string TriggersToTest[NTRIGGERS] = {
	"HLT_Dimuon25_Jpsi","HLT_Dimuon20_Jpsi_Barrel_Seagulls",
	"HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced"};
      /*
         // para el 2016
      std::string TriggersToTest[NTRIGGERS] = {
	 "HLT_Dimuon16_Jpsi","HLT_Dimuon20_Jpsi",
	 "HLT_DoubleMu4_JpsiTrk_Displaced","HLT_Dimuon10_Jpsi_Barrel"};
      */
	 
      for (unsigned int i = 0; i < NTRIGGERS; i++) {
	for (int version = 1; version < 19; version++) {
	  std::stringstream ss;
	  ss << TriggersToTest[i] << "_v" << version;
	  unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
	  if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
	    trigger += (1<<i);
	    break;
	  }
	}
      }

      if (triggerObjects.isValid() && trigger != 0 ) {
	//std::cout << "will try to match trigger object with track " << triggerObjects->size() << " with " << triggerResults_handle->size() << " trigger " << trigger << endl;
	for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	  std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	  std::string cc1 = obj.collection();
	  for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	    if ( pathNamesAll[h].find(ss0) != std::string::npos && cc1.find(cc0) != std::string::npos ) {
	      obj_eta.push_back(obj.eta());
	      obj_phi.push_back(obj.phi());
	      /*
		std::cout << pathNamesAll[h] << std::endl;
		std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
		std::cout << "\t   Collection: " << obj.collection() << std::endl;
		std::cout << "\t   Type IDs:   ";
		for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
		std::cout << std::endl;*/
	    }
	  }
	}
      } else std::cout << "*** NO triggerObjects found " << iEvent.id().run() << "," << iEvent.id().event() << " for trigger bits = " << trigger << std::endl; 
      
 } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(BSLabel_, beamSpotHandle);
  //iEvent.getByLabel(BSLabel_, beamSpotHandle);
  if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
  else std::cout << "No beam spot available from EventSetup" << endl;
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  //edm::Handle<std::vector<reco::Vertex> > primaryVertices_handle;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  //priVtxXE = bestVtx.xError();
  //priVtxYE = bestVtx.yError();
  //priVtxZE = bestVtx.zError();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);

  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  nVtx = primaryVertices_handle->size(); 

  //*****************************************
  //Let's begin by looking for J/psi

  unsigned int nMu_tmp = thePATMuonHandle->size();
  nMu = nMu_tmp;

  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;

	  // "softmuons". If commented is because the condition is the configuration file (.py)
	  //if( !(iMuon1->isSoftMuon(bestVtx)) ) continue;
	  //if( !(iMuon2->isSoftMuon(bestVtx)) ) continue;	 
	  
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************
	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;
	  
	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  VirtualKinematicParticleFactory vFactory;
		  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	  muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	 	  
	  KinematicParticleVertexFitter fitter;   	  
	  RefCountedKinematicTree psiVertexFitTree;
	  try{
	  psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }
	  //psiVertexFitTree = fitter.fit(muonParticles); 	  
	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }
	  
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  	  
	  //if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	  double Omb_J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(Omb_J_Prob_tmp<0.01)continue;
	  	  
	  //Now that we have a J/psi candidate, we look for K^+ P^- candidates

	  for(std::vector<pat::GenericParticle>::const_iterator iTrack1 = thePATTrackHandle->begin();
	       iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) 
	     {
	       // ************* offline selection for pions ************
	       //pat::GenericParticle patTrack1 = *iTrack1;			       
	       if(iTrack1->track()->charge()==0) continue;
	       if(iTrack1->track()->pt()<1.0) continue;// Muy fuerte? Es que si lo relajamos tarda mucho. Y para ser honesto no creo que este por debajo de 1MeV los cortes finales que usemos
	       if(iTrack1->track()->hitPattern().numberOfValidPixelHits()<1)continue;
	       if(iTrack1->track()->numberOfValidHits()<5)continue;
	       if(!(iTrack1->track()->quality(reco::TrackBase::highPurity))) continue;
	       //Now let's checks if our muons do not use the same tracks as we are using now
	       if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
	       
	       for(std::vector<pat::GenericParticle>::const_iterator iTrack2 = iTrack1+1;
		   iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) {
		 
		 if(iTrack2->track()==iTrack1->track()) continue;
		 if(iTrack2->track()->charge()==0) continue;
		 if(iTrack2->track()->pt()<1.0) continue;
		 if(iTrack2->track()->hitPattern().numberOfValidPixelHits()<1)continue;
		 if(iTrack2->track()->numberOfValidHits()<5)continue;
		 if(!(iTrack2->track()->quality(reco::TrackBase::highPurity))) continue;
		 //std::cout << "Pt pion2:  "<< iTrack2->pt() << std::endl;
		 if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;

		 //opposite charge?
		 //if(iTrack1->track()->charge()==iTrack2->track()->charge()) continue;	
		 		   		   
		 reco::TransientTrack pion1TT((*theB).build(iTrack1->track()));
		 reco::TransientTrack pion2TT((*theB).build(iTrack2->track()));

		 ParticleMass kaon_mass = 0.493677;
		 float kaon_sigma = kaon_mass*1.e-6;		 
		   
		 float chi = 0.;
		 float ndf = 0.;

		 // *************************************************
		 // pipi invariant mass (before kinematic vertex fit)
		 // *************************************************
		 TLorentzVector pion14V,pion24V,Jpsi4V; 
		 pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass);
		 pion24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),kaon_mass);
		 
		 if((pion14V + pion24V).M()<0.970 || (pion14V + pion24V).M()>1.070) continue;
		 //if((pion14V + pion24V).M()<0.990 || (pion14V + pion24V).M()>1.050) continue;

		 // ***********************************************
		 // Bs invariant mass (before kinematic vertex fit)
		 // ***********************************************
		 		 
		 Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());

		 // Podemos discutir este corte si desean. Es solo para acelerar el proceso en la GRID
		 // ( (pion14V + pion24V + Jpsi4V).M()>7.0 ) continue;
		 if ( (pion14V + pion24V + Jpsi4V).M()<4.6 || (pion14V + pion24V + Jpsi4V).M()>6.4 ) continue;
		 
		 //Now we are ready to combine!
		 // JPsi mass constraint is applied in the final Lambdab(Lb) fit,
		 
		 vector<RefCountedKinematicParticle> vFitMCParticles;
		 vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		 vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		 vFitMCParticles.push_back(pFactory.particle(pion1TT,kaon_mass ,chi,ndf,kaon_sigma));
		 vFitMCParticles.push_back(pFactory.particle(pion2TT,kaon_mass ,chi,ndf,kaon_sigma));
		 
		 MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
		 KinematicConstrainedVertexFitter kcvFitter;
		 RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
		 if (!vertexFitTree->isValid()) {
		   //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		   continue;
		 }
		 vertexFitTree->movePointerToTheTop();
		 
		 RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		 RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		 if (!bDecayVertexMC->vertexIsValid()){
		   // cout << "B MC fit vertex is not valid" << endl;
		   continue;
		 }

		 // esta ventana de masa esta bien? son del orden de 700 MeV de cada lado comparado a la masa del PDG
		 if ( (bCandMC->currentState().mass() < 5.0) || (bCandMC->currentState().mass() > 6.0) ) {
		   // cout << "continue from bmass > 6.4 or < 4.8 = " << bCandMC->currentState().mass() << endl;
		   continue;
		 }
		 
		 double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		 if(B_Prob_tmp<0.01)continue;

		 // get children from final Lambdab fit with mass constrain		 
		 vertexFitTree->movePointerToTheFirstChild();
		 RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		 
		 vertexFitTree->movePointerToTheNextChild();
		 RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
		 
		 vertexFitTree->movePointerToTheNextChild();
		 RefCountedKinematicParticle T1CandMC = vertexFitTree->currentParticle();
		 
		 vertexFitTree->movePointerToTheNextChild();
		 RefCountedKinematicParticle T2CandMC = vertexFitTree->currentParticle();
		 
		 KinematicParameters Mu1KP = mu1CandMC->currentState().kinematicParameters();
		 KinematicParameters Mu2KP = mu2CandMC->currentState().kinematicParameters();		   
		 KinematicParameters Pi1KP = T1CandMC->currentState().kinematicParameters();// this is the kaon1 
		 KinematicParameters Pi2KP = T2CandMC->currentState().kinematicParameters();// this is the Kaon2		 
		   
		 // ********************* loop over all the primary vertices and we choose the one with the best pointing angle **************** 
		 reco::Vertex bestVtxIP;
		 
		 Double_t pVtxIPX_temp = -10000.0;
		 Double_t pVtxIPY_temp = -10000.0;
		 Double_t pVtxIPZ_temp = -10000.0;
		 Double_t pVtxIPXE_temp = -10000.0;
		 Double_t pVtxIPYE_temp = -10000.0;
		 Double_t pVtxIPZE_temp = -10000.0;
		 Double_t pVtxIPXYE_temp = -10000.0;
		 Double_t pVtxIPXZE_temp = -10000.0;
		 Double_t pVtxIPYZE_temp = -10000.0;
		 Double_t pVtxIPCL_temp = -10000.0;	
		 Double_t lip1 = -1000000.0;
		 for(size_t i = 0; i < primaryVertices_handle->size(); ++i) {
		   const Vertex &vtx = (*primaryVertices_handle)[i];
		   
		   Double_t dx1 = (*bDecayVertexMC).position().x() - vtx.x(); 
		   Double_t dy1 = (*bDecayVertexMC).position().y() - vtx.y();
		   Double_t dz1 = (*bDecayVertexMC).position().z() - vtx.z();
		   float cosAlphaXYb1 = ( bCandMC->currentState().globalMomentum().x() * dx1 + bCandMC->currentState().globalMomentum().y()*dy1 + bCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* bCandMC->currentState().globalMomentum().mag() );
		   
		   if(cosAlphaXYb1>lip1)
		     {
		       lip1 = cosAlphaXYb1 ;
		       pVtxIPX_temp = vtx.x();
		       pVtxIPY_temp = vtx.y();
		       pVtxIPZ_temp = vtx.z();
		       pVtxIPXE_temp = vtx.covariance(0, 0);
		       pVtxIPYE_temp = vtx.covariance(1, 1);
		       pVtxIPZE_temp = vtx.covariance(2, 2);
		       pVtxIPXYE_temp = vtx.covariance(0, 1);
		       pVtxIPXZE_temp = vtx.covariance(0, 2);
		       pVtxIPYZE_temp = vtx.covariance(1, 2);
		       pVtxIPCL_temp = (TMath::Prob(vtx.chi2(),(int)vtx.ndof()) );
		       
		       bestVtxIP = vtx;
		       
		     }		   
		 }
		 
		 // try refitting the primary without the tracks in the B reco candidate		     		     
		 const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(iMuon1->originalObject());
		 const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(iMuon2->originalObject());
		 reco::TrackRef patTrack1_1 = iTrack1->track();
		 reco::TrackRef patTrack2_2 = iTrack2->track();
		 
		 // first get tracks from the original primary
		 vector<reco::TransientTrack> vertexTracks;
		 
		 for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtxIP.tracks_begin();
		       iTrack != bestVtxIP.tracks_end(); ++iTrack) {
		   // compare primary tracks to check for matches with B cand
		   reco::TrackRef trackRef = iTrack->castTo<TrackRef>();
		   
		   // the 4 tracks in the Lb candidate are  patTrack1_1, patTrack2_2, rmu1 and rmu2 
		   if (  !( (patTrack1_1.key()==trackRef.key()) || (patTrack2_2.key()==trackRef.key()) ||
			    (rmu1->track().key()==trackRef.key()) || (rmu2->track().key()==trackRef.key()) ) ) {
		     
		     //TransientTrack tt(trackRef, &(*bFieldHandle) );
		     reco::TransientTrack tt((*theB).build(trackRef));
		     vertexTracks.push_back(tt);
		   }//else { std::cout << "found track match with primary" << endl;}
		 }
		 
		 // *** if no tracks in primary or no reco track included in primary then don't do anything ***
		 reco::Vertex bestVtxRf = bestVtxIP;
		 GlobalPoint PVRfP = GlobalPoint( bestVtxIP.x(), bestVtxIP.y(), bestVtxIP.z() );
		 
		 if (  vertexTracks.size()>0 && (bestVtxIP.tracksSize()!=vertexTracks.size()) ) {
		   AdaptiveVertexFitter theFitter;
		   TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
		   if ( v.isValid() ) {		    
		     //set bestVtxRf as new best vertex to fill variables for refitting PV
		     bestVtxRf = reco::Vertex(v);
		     
		   }
		 }
		 
		 // ************ fill candidate variables now
		 
		 // You can get the momentum components from the final Bs childrens		 
		 J_px1->push_back( Mu1KP.momentum().x() );
		 J_py1->push_back( Mu1KP.momentum().y() );
		 J_pz1->push_back( Mu1KP.momentum().z() );
		 
		 J_px2->push_back(Mu2KP.momentum().x());
		 J_py2->push_back(Mu2KP.momentum().y());
		 J_pz2->push_back(Mu2KP.momentum().z());
		 
		 B_k1_px->push_back( Pi1KP.momentum().x() );
		 B_k1_py->push_back( Pi1KP.momentum().y() );
		 B_k1_pz->push_back( Pi1KP.momentum().z() );

		 B_k2_px->push_back( Pi2KP.momentum().x() );
		 B_k2_py->push_back( Pi2KP.momentum().y() );
		 B_k2_pz->push_back( Pi2KP.momentum().z() );		 
		 
		 // now we gill fill our "nominal" variables
		 
		 pi1_trackerhits->push_back(iTrack1->track()->numberOfValidHits() );
		 pi1_pixelhits->push_back(iTrack1->track()->hitPattern().numberOfValidPixelHits() );
		 pi2_trackerhits->push_back(iTrack2->track()->numberOfValidHits() );
		 pi2_pixelhits->push_back(iTrack2->track()->hitPattern().numberOfValidPixelHits() );
		 
		 pi1dz->push_back(iTrack1->track()->dz(bestVtxIP.position()) );
		 pi1dzE->push_back(iTrack1->track()->dzError() );
		 pi1dxy->push_back(iTrack1->track()->dxy(bestVtxIP.position()) );
		 pi1dxyE->push_back(iTrack1->track()->dxyError() );
		 pi1d0->push_back(iTrack1->track()->dxy(beamSpot.position()) );

		 pi2dz->push_back(iTrack2->track()->dz(bestVtxIP.position()) );
		 pi2dzE->push_back(iTrack2->track()->dzError() );
		 pi2dxy->push_back(iTrack2->track()->dxy(bestVtxIP.position()) );
		 pi2dxyE->push_back(iTrack2->track()->dxyError() );
		 		 
		 B_mass->push_back(bCandMC->currentState().mass());
		 B_px->push_back(bCandMC->currentState().globalMomentum().x());
		 B_py->push_back(bCandMC->currentState().globalMomentum().y());
		 B_pz->push_back(bCandMC->currentState().globalMomentum().z());

		 B_phi_mass->push_back( (pion14V + pion24V).M() );
		 
		 B_k_px_track->push_back(iTrack1->px() );
		 B_k_py_track->push_back(iTrack1->py() );
		 B_k_pz_track->push_back(iTrack1->pz() );
		 B_k_charge1->push_back(iTrack1->charge() );
		 
		 B_k2_px_track->push_back(iTrack2->px() );
		 B_k2_py_track->push_back(iTrack2->py() );
		 B_k2_pz_track->push_back(iTrack2->pz() );
		 B_k_charge2->push_back(iTrack2->charge() );
		 
		 B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		 B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		 B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		 B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
		 
		 B_J_px1->push_back(iMuon1->track()->px());
		 B_J_py1->push_back(iMuon1->track()->py());
		 B_J_pz1->push_back(iMuon1->track()->pz());		    
		 B_J_charge1->push_back(iMuon1->charge());
		 
		 B_J_px2->push_back(iMuon2->track()->px());
		 B_J_py2->push_back(iMuon2->track()->py());
		 B_J_pz2->push_back(iMuon2->track()->pz());
		 B_J_charge2->push_back(iMuon2->charge());	   
		 
		 B_Prob    ->push_back(B_Prob_tmp);
		 B_J_Prob  ->push_back(Omb_J_Prob_tmp);			    
		 
		 B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
		 B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
		 B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
		 
		 B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
		 B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
		 B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
		 B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
		 B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
		 B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());		  
		 
		 B_J_DecayVtxX ->push_back( psi_vFit_vertex_noMC->position().x() );
		 B_J_DecayVtxY ->push_back( psi_vFit_vertex_noMC->position().y() );
		 B_J_DecayVtxZ ->push_back( psi_vFit_vertex_noMC->position().z() );
		 
		 B_J_DecayVtxXE ->push_back( psi_vFit_vertex_noMC->error().cxx() );
		 B_J_DecayVtxYE ->push_back( psi_vFit_vertex_noMC->error().cyy() );
		 B_J_DecayVtxZE ->push_back( psi_vFit_vertex_noMC->error().czz() );
		 B_J_DecayVtxXYE ->push_back( psi_vFit_vertex_noMC->error().cyx() );
		 B_J_DecayVtxXZE ->push_back( psi_vFit_vertex_noMC->error().czx() );
		 B_J_DecayVtxYZE ->push_back( psi_vFit_vertex_noMC->error().czy() );
		 
		 pVtxIPX->push_back( pVtxIPX_temp);
		 pVtxIPY->push_back(  pVtxIPY_temp);	    
		 pVtxIPZ->push_back(  pVtxIPZ_temp);
		 pVtxIPXE->push_back( pVtxIPXE_temp);
		 pVtxIPYE->push_back( pVtxIPYE_temp);	    
		 pVtxIPZE->push_back( pVtxIPZE_temp);
		 pVtxIPXYE->push_back( pVtxIPXYE_temp);
		 pVtxIPXZE->push_back( pVtxIPXZE_temp);	    
		 pVtxIPYZE->push_back( pVtxIPYZE_temp);
		 pVtxIPCL->push_back(  pVtxIPCL_temp);
		 
		 priRfVtxX->push_back( bestVtxRf.x() );
		 priRfVtxY->push_back( bestVtxRf.y() );
		 priRfVtxZ->push_back( bestVtxRf.z() );
		 priRfVtxXE->push_back( bestVtxRf.covariance(0, 0) );
		 priRfVtxYE->push_back( bestVtxRf.covariance(1, 1) );
		 priRfVtxZE->push_back( bestVtxRf.covariance(2, 2) );
		 priRfVtxXYE->push_back( bestVtxRf.covariance(0, 1) );
		 priRfVtxXZE->push_back( bestVtxRf.covariance(0, 2) );
		 priRfVtxYZE->push_back( bestVtxRf.covariance(1, 2) );		  
		 priRfVtxCL->push_back( ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );
		 priRfNTrkDif->push_back( bestVtxIP.tracksSize() - vertexTracks.size() );
		 
		 
		 // Pion Trigger mach (just for JpsiTk trigger). Probably we can do the same for the other trigger but it could be the next interaction
		 float dr0 = 99999.;
		 for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		   float dp = iTrack1->phi() - obj_phi[ii];
		   float de = iTrack1->eta() - obj_eta[ii];
		   if (dp>float(M_PI)) dp-=float(2*M_PI);  
		   float dr = std::sqrt(de*de + dp*dp);
		   if (dr < dr0) dr0 = dr;
		   //std::cout << "\tTrigger object: eta " << obj_eta[ii] << " phi " << obj_phi[ii] << " -> " << dr0 << std::endl;
		 }
		 
		 // here we will check for muon-trigger-machint  (iMuon1 and iMuon2)			    
		 const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
		 const pat::TriggerObjectStandAloneCollection muHLTMatches2_t2 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
		 
		 const pat::TriggerObjectStandAloneCollection muHLTMatches1_t4 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");
		 const pat::TriggerObjectStandAloneCollection muHLTMatches2_t4 = iMuon2->triggerObjectMatchesByFilter("hltJpsiTkTkVertexFilterPhiKstar");
		 
		 int tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		 
		 //if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1;
		 if (muHLTMatches1_t2.size() > 0 && muHLTMatches2_t2.size() > 0) tri_JpsiTk_tmp = 1+int(1000*dr0);
		 if (muHLTMatches1_t4.size() > 0 && muHLTMatches2_t4.size() > 0) tri_JpsiTkTk_tmp = 1;
		 
		 tri_JpsiTk->push_back( tri_JpsiTk_tmp );
		 tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );
		 
		 nB++;	       
		 muonParticles.clear();
		 vFitMCParticles.clear();
		 
	       }
	     }
	}
    }
    
   if (nB > 0 ) tree_->Fill();

   nB = 0; nMu = 0;
   trigger = 0;

   J_px1->clear();  J_py1->clear();  J_pz1->clear(); 
   J_px2->clear();  J_py2->clear();  J_pz2->clear();
   B_k1_px->clear(); B_k1_py->clear(); B_k1_pz->clear();
   B_k2_px->clear(); B_k2_py->clear(); B_k2_pz->clear();

   pi1_trackerhits->clear(); pi2_trackerhits->clear(); 
   pi1_pixelhits->clear(); pi2_pixelhits->clear();
   pi1dz->clear(); pi1dzE->clear();
   pi1dxy->clear(); pi1dxyE->clear();
   pi1d0->clear();
   pi2dz->clear(); pi2dzE->clear();
   pi2dxy->clear(); pi2dxyE->clear();
   
   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   B_phi_mass->clear(); 
   B_k_charge1->clear(); B_k_px_track->clear(); B_k_py_track->clear(); B_k_pz_track->clear();
   B_k_charge2->clear(); B_k2_px_track->clear(); B_k2_py_track->clear(); B_k2_pz_track->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(); B_J_charge1->clear();
   B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(); B_J_charge2->clear();

   B_Prob->clear(); B_J_Prob->clear();

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();

   B_J_DecayVtxX->clear();   B_J_DecayVtxY->clear();   B_J_DecayVtxZ->clear();
   B_J_DecayVtxXE->clear();  B_J_DecayVtxYE->clear();  B_J_DecayVtxZE->clear();
   B_J_DecayVtxXYE->clear(); B_J_DecayVtxXZE->clear(); B_J_DecayVtxYZE->clear();

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear();

   priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   priRfVtxZE->clear(); priRfVtxXYE->clear(); priRfVtxXZE->clear(); priRfVtxYZE->clear(); priRfVtxCL->clear(); 
   priRfNTrkDif->clear(); 
   
   tri_JpsiTkTk->clear(); tri_JpsiTk->clear();

}

bool Bsphi_PAT::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method to check the presence of pixel hits  ------------
bool Bsphi_PAT::hasFirstLayerPixelHits(const reco::TransientTrack& track)
{
  using namespace reco;
  const HitPattern& p = track.hitPattern();      
  for (int i=0; i<p.numberOfAllHits(HitPattern::TRACK_HITS); i++) {
    uint32_t pattern = p.getHitPattern(HitPattern::TRACK_HITS, i);   
    if (p.pixelBarrelHitFilter(pattern) || p.pixelEndcapHitFilter(pattern) ) {
      if (p.getLayer(pattern) == 1) {
    if (p.validHitFilter(pattern)) {
      return true;
    }
      }
    }
  }
  return false;
} 

// ------------ method called once each job just before starting event loop  ------------

void
Bsphi_PAT::beginJob()
{
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Lb ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  
  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  
  tree_->Branch("B_k1_px", &B_k1_px);
  tree_->Branch("B_k1_py", &B_k1_py);
  tree_->Branch("B_k1_pz", &B_k1_pz);
  
  tree_->Branch("B_k2_px", &B_k2_px);
  tree_->Branch("B_k2_py", &B_k2_py);
  tree_->Branch("B_k2_pz", &B_k2_pz);

  tree_->Branch("pi1_trackerhits", &pi1_trackerhits);
  tree_->Branch("pi2_trackerhits", &pi2_trackerhits);
  tree_->Branch("pi1_pixelhits", &pi1_pixelhits);
  tree_->Branch("pi2_pixelhits", &pi2_pixelhits);

  tree_->Branch("pi1dz", &pi1dz);
  tree_->Branch("pi1dzE", &pi1dzE);
  tree_->Branch("pi1dxy", &pi1dxy);
  tree_->Branch("pi1dxyE", &pi1dxyE);
  tree_->Branch("pi1d0", &pi1d0);

  tree_->Branch("pi2dz", &pi2dz);
  tree_->Branch("pi2dzE", &pi2dzE);
  tree_->Branch("pi2dxy", &pi2dxy);
  tree_->Branch("pi2dxyE", &pi2dxyE);

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_phi_mass", &B_phi_mass);

  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_k_charge2", &B_k_charge2);
  tree_->Branch("B_k2_px_track", &B_k2_px_track);
  tree_->Branch("B_k2_py_track", &B_k2_py_track);
  tree_->Branch("B_k2_pz_track", &B_k2_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",    &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",    &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",    &B_DecayVtxYZE);
 
  tree_->Branch("B_J_DecayVtxX",   &B_J_DecayVtxX);
  tree_->Branch("B_J_DecayVtxY",   &B_J_DecayVtxY);
  tree_->Branch("B_J_DecayVtxZ",   &B_J_DecayVtxZ);
  tree_->Branch("B_J_DecayVtxXE",  &B_J_DecayVtxXE);
  tree_->Branch("B_J_DecayVtxYE",  &B_J_DecayVtxYE);
  tree_->Branch("B_J_DecayVtxZE",  &B_J_DecayVtxZE);
  tree_->Branch("B_J_DecayVtxXYE",  &B_J_DecayVtxXYE);
  tree_->Branch("B_J_DecayVtxXZE",  &B_J_DecayVtxXZE);
  tree_->Branch("B_J_DecayVtxYZE",  &B_J_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

  tree_->Branch("priRfVtxX",&priRfVtxX);
  tree_->Branch("priRfVtxY",&priRfVtxY);
  tree_->Branch("priRfVtxZ",&priRfVtxZ);
  tree_->Branch("priRfVtxXE",&priRfVtxXE);
  tree_->Branch("priRfVtxYE",&priRfVtxYE);
  tree_->Branch("priRfVtxZE",&priRfVtxZE);
  tree_->Branch("priRfVtxXYE",&priRfVtxXYE);
  tree_->Branch("priRfVtxXZE",&priRfVtxXZE);
  tree_->Branch("priRfVtxYZE",&priRfVtxYZE);
  tree_->Branch("priRfVtxCL",&priRfVtxCL);
  tree_->Branch("priRfNTrkDif",&priRfNTrkDif);
 
  tree_->Branch("nVtx", &nVtx);
  tree_->Branch("run", &run,       "run/I");
  tree_->Branch("event", &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
      
  // *************************

  tree_->Branch("trigger", &trigger, "trigger/i");
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);

  // gen
  if (isMC_) {
     tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
     tree_->Branch("gen_phi_p4",   "TLorentzVector",  &gen_phi_p4);
     tree_->Branch("gen_pion1_p4",  "TLorentzVector",  &gen_pion1_p4);
     tree_->Branch("gen_pion2_p4",  "TLorentzVector",  &gen_pion2_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
  }

}

double Bsphi_PAT::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

// ------------ method called once each job just after ending the event loop  ------------
void Bsphi_PAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bsphi_PAT);