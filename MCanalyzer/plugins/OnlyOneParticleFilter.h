#ifndef ONLYONEPARTICLEFILTER_h
#define ONLYONEPARTICLEFILTER_h
// -*- C++ -*-
//
// Package:    OnlyOneParticle
// Class:      OnlyOneParticle
// 
/**\class OnlyOneParticle OnlyOneParticle.cc 

 Description: Filter events with only one particle

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Horacio Crotte Ledesma
//         Created:  Sep 7 2020
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Pythia8/Pythia.h"

//
// class decleration
//
namespace edm {
  class HepMCProduct;
}

class OnlyOneParticle : public edm::global::EDFilter<> {
 public:
  explicit OnlyOneParticle(const edm::ParameterSet&);
  ~OnlyOneParticle() override;
  
  
  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
 private:
  const int fVerbose;  
  const edm::EDGetTokenT<edm::HepMCProduct> token_;
  //std::vector<int> dauIDs;
  const int particleID;
  const bool chargeconju; 
  //const int ndaughters;
  //std::vector<double> minptcut;
  //const double maxptcut;
  // std::vector<double> minetacut;
  // std::vector<double> maxetacut;
  std::unique_ptr<Pythia8::Pythia> fLookupGen; // this instance is for accessing particleData information
};
#endif