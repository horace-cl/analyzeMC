#include "OnlyOneParticleFilter.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <iostream>
#include <vector>

using namespace edm;
using namespace std;
using namespace Pythia8;


OnlyOneParticleFilter::OnlyOneParticleFilter(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter("verbose",0)),
  token_(consumes<edm::HepMCProduct>(edm::InputTag(iConfig.getUntrackedParameter("moduleLabel",std::string("generator")),"unsmeared"))),
  particleID(iConfig.getUntrackedParameter("ParticleID", 0)),
  chargeconju(iConfig.getUntrackedParameter("ChargeConjugation", true))
{
  edm::LogInfo("OnlyOneParticleFilter") << "----------------------------------------------------------------------" << endl;
  edm::LogInfo("OnlyOneParticleFilter") << "--- OnlyOneParticleFilter" << endl;
  edm::LogInfo("OnlyOneParticleFilter") << "Creating pythia8 instance for particle properties" << endl;
  if(!fLookupGen.get()) fLookupGen.reset(new Pythia());
}


OnlyOneParticleFilter::~OnlyOneParticleFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool OnlyOneParticleFilter::filter(edm::StreamID,edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  using namespace edm;
  bool accepted = false;
  Handle<HepMCProduct> evt;
  iEvent.getByToken(token_, evt);

  int nParticles(0); 

  
  HepMC::GenEvent *myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
  
  
  if (fVerbose > 5) {
    edm::LogInfo("OnlyOneParticleFilter") << "looking for " << particleID << endl;
  }
  
  
  for (HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p) {
    
    if ( abs((*p)->pdg_id()) == particleID) {
      nParticles++;
      accepted = true;
    }

    if (nParticles>1){
      accepted = false;
      break;
    }

  }



  delete myGenEvent; 
  return accepted;
}
 //    // -- Check for mother of this particle
 //    if (0 != motherID) {
 //      OK = 0; 
 //      for (HepMC::GenVertex::particles_in_const_iterator des = (*p)->production_vertex()->particles_in_const_begin(); 
	//    des != (*p)->production_vertex()->particles_in_const_end();
	//    ++des) {
	// if (fVerbose > 10) {
	// 	edm::LogInfo("OnlyOneParticleFilter") << "mother: " << (*des)->pdg_id() << " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta() << endl;
	// }
	// if (abs(motherID) == abs((*des)->pdg_id())) {
	//   OK = 1; 
	//   break;
	// }
 //      }
 //    }
 //    if (0 == OK) continue; 

 //    // -- check for daugthers
 //    int ndauac = 0;
 //    int ndau = 0;     
 //    if (fVerbose > 5) {
	//     edm::LogInfo("OnlyOneParticleFilter") << "found ID: " << (*p)->pdg_id() << " pT: " << (*p)->momentum().perp() << " eta: " << (*p)->momentum().eta() << endl;
 //    }
 //    vparticles.push_back((*p)->pdg_id()); 
 //    if ((*p)->end_vertex()) {	
 //      for (HepMC::GenVertex::particle_iterator des=(*p)->end_vertex()->particles_begin(HepMC::children);
	//    des != (*p)->end_vertex()->particles_end(HepMC::children);
	//    ++des) {
	// ++ndau;       
	// if (fVerbose > 5) {
	// 	edm::LogInfo("OnlyOneParticleFilter") << "ID: " << (*des)->pdg_id() << " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta() << endl;
	// }
	// for (unsigned int i=0; i<dauIDs.size(); ++i) {
	//   if ((*des)->pdg_id() != dauIDs[i] ) continue ;
	//   if (fVerbose > 5) {
	// 	  edm::LogInfo("OnlyOneParticleFilter") << "i = " << i << " pT = " << (*des)->momentum().perp() << " eta = " << (*des)->momentum().eta() << endl;
	//   }
	//   if ((*des)->momentum().perp() >  minptcut[i]  &&
	//       (*des)->momentum().perp() <  maxptcut  &&
	//       (*des)->momentum().eta()  >  minetacut[i] && 
	//       (*des)->momentum().eta()  <  maxetacut[i] ) {
	//     ++ndauac;
	//     vparticles.push_back((*des)->pdg_id()); 
	//     if (fVerbose > 2) {
	// 	    edm::LogInfo("OnlyOneParticleFilter") << "  accepted this particle " <<  (*des)->pdg_id()
	// 	   << " pT = " << (*des)->momentum().perp() << " eta = " << (*des)->momentum().eta() << endl;
	//     }
	//     break;
	//   } 
	// }	       		     
 //      }
 //    }  

 //    // -- allow photons
 //    if (ndau >=  ndaughters && ndauac == ndaughters) {
 //      accepted = true;
 //      if (fVerbose > 0) {
 //         edm::LogInfo("OnlyOneParticleFilter") << "  accepted this decay: ";
	// for (unsigned int iv = 0; iv < vparticles.size(); ++iv) edm::LogInfo("OnlyOneParticleFilter") << vparticles[iv] << " "; 
	// edm::LogInfo("OnlyOneParticleFilter") << " from mother = " << motherID << endl;
 //      }
 //      break;
 //    }    
    
 //  }
  
  
 //  if (!accepted && chargeconju ) {
 //    OK = 1; 

 //    for (HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	//  p != myGenEvent->particles_end(); ++p) {
      
 //      if ((*p)->pdg_id() != -particleID) continue ;
      
 //      // -- Check for mother of this particle
 //      if (0 != motherID) {
	// OK = 0; 
	// for (HepMC::GenVertex::particles_in_const_iterator des = (*p)->production_vertex()->particles_in_const_begin(); 
	//      des != (*p)->production_vertex()->particles_in_const_end();
	//      ++des) {
	//   if (fVerbose > 10) {
	// 	  edm::LogInfo("OnlyOneParticleFilter") << "mother: " << (*des)->pdg_id() << " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta() << endl;
	//   }
	//   if (abs(motherID) == abs((*des)->pdg_id())) {
	//     OK = 1; 
	//     break;
	//   }
	// }
 //      }
 //      if (0 == OK) continue; 
      
 //      if (fVerbose > 5) {
	//       edm::LogInfo("OnlyOneParticleFilter") << "found ID: " << (*p)->pdg_id() << " pT: " << (*p)->momentum().perp() << " eta: " << (*p)->momentum().eta() << endl;
 //      }
 //      vparticles.push_back((*p)->pdg_id()); 
 //      int ndauac = 0;
 //      int ndau = 0;     
 //      if ((*p)->end_vertex()) {
	// for (HepMC::GenVertex::particle_iterator des=(*p)->end_vertex()->particles_begin(HepMC::children);
	//      des != (*p)->end_vertex()->particles_end(HepMC::children);
	//      ++des) {
	//   ++ndau;
	//   if (fVerbose > 5) {
	// 	  edm::LogInfo("OnlyOneParticleFilter") << "ID: " << (*des)->pdg_id() << " pT: " << (*des)->momentum().perp() << " eta: " << (*des)->momentum().eta() << endl;
	//   }
	//   for (unsigned int i=0; i<dauIDs.size(); ++i) {
	//     int IDanti = -dauIDs[i];
 //            if ( !(fLookupGen->particleData.isParticle( IDanti )) ) IDanti = dauIDs[i]; 
	//     if ((*des)->pdg_id() != IDanti) continue ;
	//     if (fVerbose > 5) {
	// 	    edm::LogInfo("OnlyOneParticleFilter") << "i = " << i << " pT = " << (*des)->momentum().perp() << " eta = " << (*des)->momentum().eta() << endl;
	//     }
	//     if ((*des)->momentum().perp() >  minptcut[i]  &&
	// 	(*des)->momentum().perp() <  maxptcut  &&
	// 	(*des)->momentum().eta()  >  minetacut[i] && 
	// 	(*des)->momentum().eta()  <  maxetacut[i] ) {
	//       ++ndauac;
	//       vparticles.push_back((*des)->pdg_id()); 
	//       if (fVerbose > 2) {
	// 	      edm::LogInfo("OnlyOneParticleFilter") << "  accepted this particle " <<  (*des)->pdg_id()
	// 	     << " pT = " << (*des)->momentum().perp() << " eta = " << (*des)->momentum().eta() << endl;
	//       }
	//       break;
	//     } 
	//   }	       		     
	// }
 //      }
 //      if (ndau >=  ndaughters && ndauac == ndaughters ) {
	// accepted = true;
	// if (fVerbose > 0) {
	//   edm::LogInfo("OnlyOneParticleFilter") << "  accepted this anti-decay: ";
	//   for (unsigned int iv = 0; iv < vparticles.size(); ++iv) edm::LogInfo("OnlyOneParticleFilter") << vparticles[iv] << " "; 
	//   edm::LogInfo("OnlyOneParticleFilter") << " from mother = " << motherID << endl;
	// }
	// break;
 //      }    
 //    }
    
 //  }    
  
 //  delete myGenEvent; 
 //  return accepted;


//define this as a plug-in
DEFINE_FWK_MODULE(OnlyOneParticleFilter);
