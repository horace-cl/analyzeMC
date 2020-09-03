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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class MCanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MCanalyzer::MCanalyzer(const edm::ParameterSet& iConfig)
 // :
 //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{

  //std::cout << "INITIALIZER?" << std::endl;
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

  std::cout << "HELLO FROM ANALYZER! " << std::endl;
  

  edm::Handle<edm::HepMCProduct> genEvtHandle;
  iEvent.getByToken(hepmcproduct_, genEvtHandle);

  const HepMC::GenEvent* Evt = genEvtHandle->GetEvent() ;
  //
  // this is an example loop over the hierarchy of vertices
  //

  // from 
  // Calibration/HcalCalibAlgos/plugins/SimAnalyzerMinbias.cc
  if (!genEvtHandle.isValid()) 
  {
      std::cout << " ------------->  no HepMCProduct found" << std::endl;    
  } 
  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(genEvtHandle->GetEvent()));
  std::cout << "Event with : \n"; 
  std::cout << myGenEvent->particles_size() << " particles \n";
  std::cout << myGenEvent->vertices_size() << " vertices\n";



  //from
  // Alignment/OfflineValidation/plugins/ValidationMisalignedTracker.cc:
  // Iterate over all particles
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) { 

    std::cout << "\tPDG ID : " << (*p)->pdg_id() << std::endl;
    std::cout << "\tSTATUS : " << (*p)->status() << std::endl;

    std::cout << "\tDaugthers : \n"<< std::endl;
    //Ierate over its daughters
    for( HepMC::GenVertex::particle_iterator aDaughter=(*p)->end_vertex()->particles_begin(HepMC::descendants); aDaughter !=(*p)->end_vertex()->particles_end(HepMC::descendants);aDaughter++){
      std::cout << "\t\tPDG ID : " << (*aDaughter)->pdg_id() << std::endl;
      std::cout << "\t\tSTATUS : " << (*aDaughter)->status() << std::endl;
    }



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

  
  HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();



  // VERTEX ITERATOR
  //FROM TWIKI 
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideDataFormatGeneratorInterface
  // int i=0;
  // int j=0;
  // std::cout << "HANDLE OBTAINED" << std::endl;
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
}

// ------------ method called once each job just after ending the event loop  ------------
void
MCanalyzer::endJob()
{
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
