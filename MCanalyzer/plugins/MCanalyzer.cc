// -*- C++ -*-
//
// Package:    Analyze/MCanalyzer
// Class:      MCanalyzer
//
/**\class MCanalyzer MCanalyzer.cc Analyze/MCanalyzer/plugins/MCanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Horacio Crotte Ledesma
//         Created:  Wed, 02 Sep 2020 19:00:43 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

/*
HCL
 I will try to follow all information given here:
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
*/
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//WE NEED THESE ONES FOR MAKING THE NTUPLES
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <utility>
#include <string>
#include "Math/GenVector/Boost.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include <Math/VectorUtil.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//
// FROM JHOVANNYS CODE
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

using reco::TrackCollection;

//class MCanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class MCanalyzer : public edm::EDAnalyzer {
   public:
      explicit MCanalyzer(const edm::ParameterSet&);
      ~MCanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::HepMCProduct> hepmcproduct_;
      //edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genCands_;
      
      TLorentzVector gen_b_p4,gen_phi_p4,gen_kaon_p4,gen_muon1_p4,gen_muon2_p4, gen_gamma1_p4, gen_gamma2_p4;
      TLorentzVector gen_b_p4J,gen_phi_p4J,gen_kaon_p4J,gen_muon1_p4J,gen_muon2_p4J, gen_gamma1_p4J, gen_gamma2_p4J;
      TVector3       gen_b_vtx;
      TTree*         tree_;
      std::vector<std::vector<int>>    daughter_id;
     // std::vector<int> number_daughters;
      int number_daughters;
      float costhetaL, costhetaKL;
};


//
//genCands_(consumes<rec::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MCanalyzer::MCanalyzer(const edm::ParameterSet& iConfig)
 :number_daughters(0), costhetaL(0.0), costhetaKL(0.0)
{
  std::cout << "INITIALIZER?" << std::endl;
  genCands_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  hepmcproduct_ = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));
  
}


MCanalyzer::~MCanalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}



//
// member functions
//

// ------------ method called for each event  ------------
void
MCanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool JHOVANNYS=false;
  bool debug = false;

  gen_b_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_kaon_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon1_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon2_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_gamma1_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma2_p4.SetPxPyPzE(0.,0.,0.,0.);

  gen_b_p4J.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaon_p4J.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4J.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4J.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_gamma1_p4J.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_gamma2_p4J.SetPtEtaPhiM(0.,0.,0.,0.);


  if (debug) std::cout << "HELLO FROM ANALYZER! " << std::endl;
 
  //edm::Handle<reco::GenParticleCollection> pruned;
  edm::Handle<std::vector<reco::GenParticle>> pruned; 
  iEvent.getByToken(genCands_, pruned);

  edm::Handle<edm:: HepMCProduct > genEvtHandle;
  iEvent.getByToken(hepmcproduct_, genEvtHandle);



  //JHOVANNYS
  //JHOVANNYS
  //JHOVANNYS
  //JHOVANNYS
  if (debug) std::cout << "PRUNED? \n";
  std::cout << "SIZE = " << pruned->size() << std::endl;
  if ( pruned.isValid() ) {
    if (debug) std::cout << "VALID SIZE = " << pruned->size() << std::endl;
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      //GETTING DAUGHTERS!
      const reco::Candidate *dau = &(*pruned)[i];
      //ONLY LOOKING FOR B+-
      if ( (abs(dau->pdgId()) == 521) ) { //&& (dau->status() == 2) ) {
            //foundit++;
            gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
            gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
            //int npion=0;
            std::cout << "NUMBER OF GARND?DAUGHTERS : "<< dau->numberOfDaughters() << std::endl;
            if (dau->numberOfDaughters()>5) continue;

            for (size_t k=0; k<dau->numberOfDaughters(); k++) {
              //GETTING GRANDAUGHTERS
              const reco::Candidate *gdau = dau->daughter(k);
              //LOOK FOR GRANDAUGHTERS TO BE K+-  443
              if ( abs(gdau->pdgId())==443 ) { //&& gdau->status()==2) { HERE JHOVANNY WAS LOOKING FOR THE JPSI(443)
                //foundit++;
                gen_kaon_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                //int nm=0;
              }
              //LOOK FOR GRANDAUGHTERS TO BE MU+- 13
              else if( abs(gdau->pdgId())==13){
                if (dau->pdgId()*gdau->pdgId()<0){
                  gen_muon1_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }
                else {
                  gen_muon2_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }
              }
              //LOOK FOR ANY DAMN PHOTON
              else if(dau->pdgId()==22){
                if (debug) std::cout << "foundit : "<< foundit<< std::endl;
                if (foundit==0){
                  gen_gamma1_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }
                else{
                  gen_gamma2_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }

                foundit++;
              }
            }
      }
    }
  }
 


  // from 
  // Calibration/HcalCalibAlgos/plugins/SimAnalyzerMinbias.cc
  if (!genEvtHandle.isValid()) 
  {
      if (debug) std::cout << " ------------->  no HepMCProduct found" << std::endl;    
  } 
  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(genEvtHandle->GetEvent()));
  if (debug) std::cout << "Event with : \n"; 
  std::cout << myGenEvent->particles_size() << " particles \n";
  if (debug) std::cout << myGenEvent->vertices_size() << " vertices\n";



  // from
  // Alignment/OfflineValidation/plugins/ValidationMisalignedTracker.cc:
  // Iterate over all particles
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) { 
    if (JHOVANNYS) break;
    if (abs((*p)->pdg_id())!=521){
      //if (debug) std::cout << "\t\tNot B+\n";
      continue;
    }

    if (debug) std::cout << "\tPDG ID : " << (*p)->pdg_id() << std::endl;
    if (debug) std::cout << "\tSTATUS : " << (*p)->status() << std::endl;

    if (debug) std::cout << "PX " << (*p)->momentum().px() << std::endl;
    if (debug) std::cout << "PY " << (*p)->momentum().py() << std::endl;
    if (debug) std::cout << "PZ " << (*p)->momentum().pz() << std::endl;
    if (debug) std::cout << "ENERGY " << (*p)->momentum().e() << std::endl;
    if (debug) std::cout << "MASS " << (*p)->momentum().m() << std::endl; 
    std::vector<int> ids;

    gen_b_p4.SetPxPyPzE((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());

    //if (debug) std::cout << "\tDaugthers : " << (*p)->numberOfDaughters() << std::endl;
    std::cout << "\tDaugthers : " << std::endl;
    int photons=0;

    
    //Ierate over its daughters
    for(
      HepMC::GenVertex::particle_iterator aDaughter=(*p)->end_vertex()->particles_begin(HepMC::descendants); 
      aDaughter !=(*p)->end_vertex()->particles_end(HepMC::descendants);
      aDaughter++)
    {
      //Just a vector to have control over all daugther particles
      ids.push_back((*aDaughter)->pdg_id());
      if (debug) std::cout << "\t\tPDG ID : " << (*aDaughter)->pdg_id() << std::endl;
      if (debug) std::cout << "\t\tSTATUS : " << (*aDaughter)->status() << std::endl;

      //KAON MUST HAVE THE SAME SIGN AS THE B CHARGED MESON
      if (abs((*aDaughter)->pdg_id())==321){
          gen_kaon_p4.SetPxPyPzE((*aDaughter)->momentum().px(),(*aDaughter)->momentum().py(),(*aDaughter)->momentum().pz(),(*aDaughter)->momentum().e());}
      // muons - INDEX 1 WILL BE ASSIGNED TO THE MUON WITH OPPOSITE SIGN WRT KAON
      // THETA L IS THE AGNLE BEWTEEN THESE TWO (KAON - MUON1)
	    else if (abs((*aDaughter)->pdg_id())==13){
        //CHECK IS THE CURRENT MUON HAS THE SAME SIGN AS THE B MESON
        if ((*p)->pdg_id()*(*aDaughter)->pdg_id() > 0){
          gen_muon2_p4.SetPxPyPzE((*aDaughter)->momentum().px(),(*aDaughter)->momentum().py(),(*aDaughter)->momentum().pz(),(*aDaughter)->momentum().e());
        }
        else {
          gen_muon1_p4.SetPxPyPzE((*aDaughter)->momentum().px(),(*aDaughter)->momentum().py(),(*aDaughter)->momentum().pz(),(*aDaughter)->momentum().e());
        }
      }
      // IN CASE THERE ARE ANY PHOTONS
      else if ((*aDaughter)->pdg_id()==22){
        if (photons==0){
          gen_gamma1_p4.SetPxPyPzE((*aDaughter)->momentum().px(),(*aDaughter)->momentum().py(),(*aDaughter)->momentum().pz(),(*aDaughter)->momentum().e());
        	photons=1;
        }
        else{
          gen_gamma2_p4.SetPxPyPzE((*aDaughter)->momentum().px(),(*aDaughter)->momentum().py(),(*aDaughter)->momentum().pz(),(*aDaughter)->momentum().e());
        }
      //if (debug) std::cout << "\t\tGrandDaughters : " << (*aDaughter)->numberOfDaughters() << std::endl;
      }
    }

    // NOW CREATE THE BOOST TO DILEPTON CM FRAME
    math::XYZTLorentzVector muon1(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(), gen_muon1_p4.E());
    math::XYZTLorentzVector muon2(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(), gen_muon2_p4.E());
    math::XYZTLorentzVector kaon(gen_kaon_p4.Px(), gen_kaon_p4.Py(), gen_kaon_p4.Pz(), gen_kaon_p4.E());
   
    math::XYZTLorentzVector dilep = muon1+muon2;
    ROOT::Math::Boost cmboost(dilep.BoostToCM());

    math::XYZTLorentzVector kaonCM(  cmboost( kaon )  );
    math::XYZTLorentzVector muonCM1(  cmboost( muon1 )  );
    math::XYZTLorentzVector muonCM2(  cmboost( muon2 )  );


    costhetaL = ( muonCM1.x()*muonCM2.x() 
                       + muonCM1.y()*muonCM2.y() 
                       + muonCM1.z()*muonCM2.z() ) / (muonCM1.P()*muonCM2.P() );

    costhetaKL = ( muonCM1.x()*kaonCM.x()
                       + muonCM1.y()*kaonCM.y()
                       + muonCM1.z()*kaonCM.z() ) / (muonCM1.P()*kaonCM.P() );

  
  if (debug) std::cout << "Number of Daugthers : "<< ids.size() <<std::endl;
  daughter_id.push_back(ids);
  number_daughters= ids.size();

  tree_->Fill();

  daughter_id.clear();
  // for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
  // p != myGenEvent->particles_end(); ++p ) 
  //   {
  //     // phiParticle = (*p)->momentum().phi();
  //     // etaParticle = (*p)->momentum().eta();
  //     // double pt  = (*p)->momentum().perp();
  //     // mom_MC = (*p)->momentum().rho();
  //     // if(pt > maxPt) { npart++; maxPt = pt; /*phi_MC = phiParticle; eta_MC = etaParticle;*/ }
  //     // GlobalVector mom ((*p)->momentum().x(),(*p)->momentum().y(),(*p)->momentum().z());
  //   }

  



  // VERTEX ITERATOR
  //FROM TWIKI 
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
  // int i=0;
  // int j=0;
  // if (debug) std::cout << "HANDLE OBTAINED" << std::endl;
  // for ( HepMC::GenEvent::vertex_const_iterator
  //           itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx )
  //   {
  //         i++;
  //         j=0;
  //         //
  //         // this is an example loop over particles coming out of each vertex in the loop
  //         //

  //         std::cout << "VERTEX ITERATOR " << i << std::endl;
  //         for ( HepMC::GenVertex::particles_out_const_iterator
  //                 itPartOut=(*itVtx)->particles_out_const_begin();
  //                 itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut )
  //           {
  //             j+=1;
  //             std::cout << "PARTICLES " << j << std::endl;
  //              // and more of your code...
  //           }
    }


}


// ------------ method called once each job just before starting event loop  ------------
void
MCanalyzer::beginJob()
{
  std::cout << "Beginning analyzer job" << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","B+->K+ mu mu ntuple");

  tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
  tree_->Branch("gen_kaon_p4",  "TLorentzVector",  &gen_kaon_p4);
  tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
  tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
  tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
  tree_->Branch("gen_gamma1_p4",  "TLorentzVector",  &gen_muon1_p4);
  tree_->Branch("gen_gamma2_p4",  "TLorentzVector",  &gen_muon1_p4);

  tree_->Branch("gen_b_p4J",     "TLorentzVector",  &gen_b_p4);
  tree_->Branch("gen_kaon_p4J",  "TLorentzVector",  &gen_kaon_p4);
  tree_->Branch("gen_muon1_p4J",  "TLorentzVector",  &gen_muon1_p4);
  tree_->Branch("gen_muon2_p4J",  "TLorentzVector",  &gen_muon2_p4);
  tree_->Branch("gen_gamma1_p4J",  "TLorentzVector",  &gen_muon1_p4);
  tree_->Branch("gen_gamma2_p4J",  "TLorentzVector",  &gen_muon1_p4);
  
  tree_->Branch("daughter_id",   "vector", &daughter_id);
  tree_->Branch("number_daughters",  &number_daughters);
  tree_->Branch("costhetaL",  &costhetaL);
  tree_->Branch("costhetaKL",  &costhetaKL);
}

// ------------ method called once each job just after ending the event loop  ------------
void
MCanalyzer::endJob()
{
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MCanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCanalyzer);
