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
      TLorentzVector gen_b_p4CM,gen_phi_p4CM,gen_kaon_p4CM,gen_muon1_p4CM,gen_muon2_p4CM, gen_gamma1_p4CM, gen_gamma2_p4CM;
      TLorentzVector gen_b_p4CMJ,gen_phi_p4CMJ,gen_kaon_p4CMJ,gen_muon1_p4CMJ,gen_muon2_p4CMJ, gen_gamma1_p4CMJ, gen_gamma2_p4CMJ;
      TVector3       gen_b_vtx;
      TTree*         tree_;
      std::vector<std::vector<int>>    daughter_id;
     // std::vector<int> number_daughters;
      int number_daughters, number_daughtersJ, bplus;
      float costhetaL, costhetaKL, costhetaLJ, costhetaKLJ;
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
 :number_daughters(0), number_daughtersJ(0), bplus(0), costhetaL(0.0), costhetaKL(0.0), costhetaLJ(0.0), costhetaKLJ(0.0)
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

  gen_b_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_kaon_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon1_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon2_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma1_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma2_p4CM.SetPxPyPzE(0.,0.,0.,0.);

  gen_b_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);
  gen_kaon_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon1_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon2_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma1_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma2_p4CMJ.SetPxPyPzE(0.,0.,0.,0.);


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
    int bplus_ = 0;
    
    for (size_t i=0; i<pruned->size(); i++) {
      //GETTING DAUGHTERS!
      int kaon_D = 0;
      int muon_D = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      //ONLY LOOKING FOR B+-
      if ( (abs(dau->pdgId()) == 521) && (dau->status() == 2) ) {
            //foundit++;
            bplus_++;
            gen_b_p4J.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
            gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
            //int npion=0;
            std::cout << "NUMBER OF GARND?DAUGHTERS : "<< dau->numberOfDaughters() << std::endl;
            number_daughtersJ= dau->numberOfDaughters();
            if (dau->numberOfDaughters()!=3) continue;

            for (size_t k=0; k<dau->numberOfDaughters(); k++) {
              //GETTING GRANDAUGHTERS
              const reco::Candidate *gdau = dau->daughter(k);
              //LOOK FOR GRANDAUGHTERS TO BE K+-  321
              if ( (abs(gdau->pdgId())==321)  && (gdau->status() == 1) ) { //&& gdau->status()==2) { HERE JHOVANNY WAS LOOKING FOR THE JPSI(443)
                kaon_D++;
                gen_kaon_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                //int nm=0;
              }
              //LOOK FOR GRANDAUGHTERS TO BE MU+- 13
              else if( (abs(gdau->pdgId())==13) && (gdau->status() == 1) ){
                muon_D++;
                if (dau->pdgId()*gdau->pdgId()<0){
                  gen_muon1_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }
                else {
                  gen_muon2_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                }
              }
              //LOOK FOR ANY DAMN PHOTON
              //VAMOS A COMNETARLO POR EL MOMENTO
              // else if((dau->pdgId()==22) && false){
              //   if (debug) std::cout << "foundit : "<< foundit<< std::endl;
              //   if (foundit==0){
              //     gen_gamma1_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
              //   }
              //   else{
              //     gen_gamma2_p4J.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
              //   }
              //   foundit++;
              // }
            }

            if ((kaon_D==1) && (muon_D==2)){
              math::XYZTLorentzVector muon1J(gen_muon1_p4J.Px(), gen_muon1_p4J.Py(), gen_muon1_p4J.Pz(), gen_muon1_p4J.E());
              math::XYZTLorentzVector muon2J(gen_muon2_p4J.Px(), gen_muon2_p4J.Py(), gen_muon2_p4J.Pz(), gen_muon2_p4J.E());
              math::XYZTLorentzVector kaonJ(gen_kaon_p4J.Px(), gen_kaon_p4J.Py(), gen_kaon_p4J.Pz(), gen_kaon_p4J.E());
              math::XYZTLorentzVector bmesonJ(gen_b_p4J.Px(), gen_b_p4J.Py(), gen_b_p4J.Pz(), gen_b_p4J.E());
              math::XYZTLorentzVector gamma1J(gen_gamma1_p4J.Px(), gen_gamma1_p4J.Py(), gen_gamma1_p4J.Pz(), gen_gamma1_p4J.E());
              math::XYZTLorentzVector gamma2J(gen_gamma2_p4J.Px(), gen_gamma2_p4J.Py(), gen_gamma2_p4J.Pz(), gen_gamma2_p4J.E());


              math::XYZTLorentzVector dilepJ = muon1J+muon2J;
              ROOT::Math::Boost dileptonCMBoost(dilepJ.BoostToCM());

              math::XYZTLorentzVector kaonCMJ(  dileptonCMBoost( kaonJ )  );
              math::XYZTLorentzVector muonCM1J(  dileptonCMBoost( muon1J )  );
              math::XYZTLorentzVector muonCM2J(  dileptonCMBoost( muon2J )  );
              math::XYZTLorentzVector bmesonCMJ(  dileptonCMBoost( bmesonJ )  );
              math::XYZTLorentzVector gamma1CMJ(  dileptonCMBoost( gamma1J )  );
              math::XYZTLorentzVector gamma2CMJ(  dileptonCMBoost( gamma2J )  );

              gen_b_p4CMJ.SetPxPyPzE(bmesonCMJ.x(), bmesonCMJ.y(), bmesonCMJ.z(), bmesonCMJ.t() ) ;
              gen_kaon_p4CMJ.SetPxPyPzE(kaonCMJ.x(), kaonCMJ.y(), kaonCMJ.z(), kaonCMJ.t() ) ;
              gen_muon1_p4CMJ.SetPxPyPzE(muonCM1J.x(), muonCM1J.y(), muonCM1J.z(), muonCM1J.t() ) ;
              gen_muon2_p4CMJ.SetPxPyPzE(muonCM2J.x(), muonCM2J.y(), muonCM2J.z(), muonCM2J.t() ) ;
              gen_gamma1_p4CMJ.SetPxPyPzE(gamma1CMJ.x(), gamma1CMJ.y(), gamma1CMJ.z(), gamma1CMJ.t() ) ;
              gen_gamma2_p4CMJ.SetPxPyPzE(gamma2CMJ.x(), gamma2CMJ.y(), gamma2CMJ.z(), gamma2CMJ.t() ) ;


              costhetaLJ = ( muonCM1J.x()*muonCM2J.x() 
                                 + muonCM1J.y()*muonCM2J.y() 
                                 + muonCM1J.z()*muonCM2J.z() ) / (muonCM1J.P()*muonCM2J.P() );

              costhetaKLJ = ( muonCM1J.x()*kaonCMJ.x()
                                 + muonCM1J.y()*kaonCMJ.y()
                                 + muonCM1J.z()*kaonCMJ.z() ) / (muonCM1J.P()*kaonCMJ.P() );
          }
      }
    }
   
  bplus=bplus_;
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
    math::XYZTLorentzVector bmeson(gen_b_p4.Px(), gen_b_p4.Py(), gen_b_p4.Pz(), gen_b_p4.E());
    math::XYZTLorentzVector gamma1(gen_gamma1_p4.Px(), gen_gamma1_p4.Py(), gen_gamma1_p4.Pz(), gen_gamma1_p4.E());
    math::XYZTLorentzVector gamma2(gen_gamma2_p4.Px(), gen_gamma2_p4.Py(), gen_gamma2_p4.Pz(), gen_gamma2_p4.E()); 
    

    math::XYZTLorentzVector dilep = muon1+muon2;
    ROOT::Math::Boost dileptonCMBoost(dilep.BoostToCM());



    math::XYZTLorentzVector kaonCM(  dileptonCMBoost( kaon )  );
    math::XYZTLorentzVector muonCM1(  dileptonCMBoost( muon1 )  );
    math::XYZTLorentzVector muonCM2(  dileptonCMBoost( muon2 )  );
    math::XYZTLorentzVector bmesonCM(  dileptonCMBoost( bmeson )  );
    math::XYZTLorentzVector gamma1CM(  dileptonCMBoost( gamma1 )  );
    math::XYZTLorentzVector gamma2CM(  dileptonCMBoost( gamma2 )  );


    gen_b_p4CM.SetPxPyPzE(bmesonCM.x(), bmesonCM.y(), bmesonCM.z(), bmesonCM.t() ) ;
    gen_kaon_p4CM.SetPxPyPzE(kaonCM.x(), kaonCM.y(), kaonCM.z(), kaonCM.t() ) ;
    gen_muon1_p4CM.SetPxPyPzE(muonCM1.x(), muonCM1.y(), muonCM1.z(), muonCM1.t() ) ;
    gen_muon2_p4CM.SetPxPyPzE(muonCM2.x(), muonCM2.y(), muonCM2.z(), muonCM2.t() ) ;
    gen_gamma1_p4CM.SetPxPyPzE(gamma1CM.x(), gamma1CM.y(), gamma1CM.z(), gamma1CM.t() ) ;
    gen_gamma2_p4CM.SetPxPyPzE(gamma2CM.x(), gamma2CM.y(), gamma2CM.z(), gamma2CM.t() ) ;


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

  tree_->Branch("gen_b_p4J",     "TLorentzVector",  &gen_b_p4J);
  tree_->Branch("gen_kaon_p4J",  "TLorentzVector",  &gen_kaon_p4J);
  tree_->Branch("gen_muon1_p4J",  "TLorentzVector",  &gen_muon1_p4J);
  tree_->Branch("gen_muon2_p4J",  "TLorentzVector",  &gen_muon2_p4J);
  tree_->Branch("gen_gamma1_p4J",  "TLorentzVector",  &gen_muon1_p4J);
  tree_->Branch("gen_gamma2_p4J",  "TLorentzVector",  &gen_muon1_p4J);

  tree_->Branch("gen_b_p4CM",     "TLorentzVector",  &gen_b_p4CM);
  tree_->Branch("gen_kaon_p4CM",  "TLorentzVector",  &gen_kaon_p4CM);
  tree_->Branch("gen_muon1_p4CM",  "TLorentzVector",  &gen_muon1_p4CM);
  tree_->Branch("gen_muon2_p4CM",  "TLorentzVector",  &gen_muon2_p4CM);
  tree_->Branch("gen_gamma1_p4CM",  "TLorentzVector",  &gen_muon1_p4CM);
  tree_->Branch("gen_gamma2_p4CM",  "TLorentzVector",  &gen_muon1_p4CM);

  tree_->Branch("gen_b_p4CMJ",     "TLorentzVector",  &gen_b_p4CMJ);
  tree_->Branch("gen_kaon_p4CMJ",  "TLorentzVector",  &gen_kaon_p4CMJ);
  tree_->Branch("gen_muon1_p4CMJ",  "TLorentzVector",  &gen_muon1_p4CMJ);
  tree_->Branch("gen_muon2_p4CMJ",  "TLorentzVector",  &gen_muon2_p4CMJ);
  tree_->Branch("gen_gamma1_p4CMJ",  "TLorentzVector",  &gen_muon1_p4CMJ);
  tree_->Branch("gen_gamma2_p4CMJ",  "TLorentzVector",  &gen_muon1_p4CMJ);
  
  tree_->Branch("daughter_id",   "vector", &daughter_id);
  tree_->Branch("number_daughters",  &number_daughters);
  tree_->Branch("costhetaL",  &costhetaL);
  tree_->Branch("costhetaKL",  &costhetaKL);

  tree_->Branch("number_daughtersJ",  &number_daughtersJ);
  tree_->Branch("costhetaLJ",  &costhetaLJ);
  tree_->Branch("costhetaKLJ",  &costhetaKLJ);

  tree_->Branch("Nbplus", &bplus);
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
