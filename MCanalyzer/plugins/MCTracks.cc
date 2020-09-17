// -*- C++ -*-
//
// Package:    Analyze/MCTracks
// Class:      MCTracks
//
/**\class MCTracks MCTracks.cc Analyze/MCTracks/plugins/MCTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Horacio Crotte Ledesma
//         Created:  Wed, 17 Sep 2020 19:00:43 GMT
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

//class MCTracks : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class MCTracks : public edm::EDAnalyzer {
   public:
      explicit MCTracks(const edm::ParameterSet&);
      ~MCTracks();

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
      
      TTree*         tree_;

			std::vector<TLorentzVector>      tracksP4;
			std::vector<int>                 pdgID;
			std::vector<int>                 motherID;
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
MCTracks::MCTracks(const edm::ParameterSet& iConfig)
 //:number_daughters(0), number_daughtersJ(0), bplus(0), costhetaL(0.0), costhetaKL(0.0), costhetaLJ(0.0), costhetaKLJ(0.0)
{
  std::cout << "INITIALIZER?" << std::endl;
  genCands_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  hepmcproduct_ = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));
  std::cout << "INITIALIZED\n";
}


MCTracks::~MCTracks()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}



//
// member functions
//

// ------------ method called for each event  ------------
void
MCTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool debug = true;

  if (debug) std::cout << "HELLO FROM ANALYZER! " << std::endl;
 
  //edm::Handle<reco::GenParticleCollection> pruned;
  edm::Handle<std::vector<reco::GenParticle>> pruned; 
  iEvent.getByToken(genCands_, pruned);

  edm::Handle<edm:: HepMCProduct > genEvtHandle;
  iEvent.getByToken(hepmcproduct_, genEvtHandle);

  if (debug) std::cout << "PRUNED? \n";
  //std::cout << "SIZE = " << pruned->size() << std::endl;
  if ( pruned.isValid() ) {
    if (debug) std::cout << "VALID SIZE = " << pruned->size() << std::endl;
    
    for (size_t i=0; i<pruned->size(); i++) {
      //GETTING DAUGHTERS!
      const reco::Candidate *dau = &(*pruned)[i];
      //ONLY LOOKING FOR B+-
     // if ((abs(dau->pdgId()) == 521) && (dau->status() == 2) ) {
		  if (dau->status() == 2){
            
            //if (dau->numberOfDaughters()3) continue;

            for (size_t k=0; k<dau->numberOfDaughters(); k++) {
              //GETTING GRANDAUGHTERS
              const reco::Candidate *gdau = dau->daughter(k);
              if ((gdau->status() == 1) && gdau->charge()!=0){
                motherID.push_back(dau->pdgId());
                pdgID.push_back(gdau->pdgId());
                TLorentzVector p4;
                p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                tracksP4.push_back(p4);

              }
            }
              
      }
    }
  }
 
  tree_->Fill();

 
  tracksP4.clear();
  pdgID.clear();
  motherID.clear();
}


// ------------ method called once each job just before starting event loop  ------------
void
MCTracks::beginJob()
{
  std::cout << "Beginning analyzer job" << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","B+->K+ mu mu ntuple");

  tree_->Branch("tracksP4",   "vector", &tracksP4);
  tree_->Branch("pdgID",   "vector", &pdgID);
  tree_->Branch("motherID",   "vector", &motherID);
}

// ------------ method called once each job just after ending the event loop  ------------
void
MCTracks::endJob()
{
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MCTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(MCTracks);
