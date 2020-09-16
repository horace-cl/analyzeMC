 // Original Author:  Andrea RIZZI
 //         Created:  Mon, 07 Jul 2014 07:56:38 GMT

 // system include files
 #include <memory>

 // user include files
 #include "FWCore/Framework/interface/Frameworkfwd.h"
 #include "FWCore/Framework/interface/EDAnalyzer.h"

 #include "FWCore/Framework/interface/Event.h"
 #include "FWCore/Framework/interface/MakerMacros.h"

 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
 #include "DataFormats/Candidate/interface/Candidate.h"
 #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
 //
 // class declaration
 //

class MiniAODGenPartAnalyzer : public edm::EDAnalyzer {
public:
explicit MiniAODGenPartAnalyzer(const edm::ParameterSet&);
~MiniAODGenPartAnalyzer();
bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle, int &calls);



private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;

edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
};

MiniAODGenPartAnalyzer::MiniAODGenPartAnalyzer(const edm::ParameterSet& iConfig):
prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
}


MiniAODGenPartAnalyzer::~MiniAODGenPartAnalyzer()
{
}

//Check recursively if any ancestor of particle is the given one
bool MiniAODGenPartAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle, int & calls)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
  // we also increase the counter by one
  calls+=1;
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(isAncestor(ancestor,particle->mother(i), calls)) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
    return false;
}

void
MiniAODGenPartAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace pat;

  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);

  //let's try to find all status1 originating directly from a B meson decay 

  for(size_t i=0; i<pruned->size();i++){
    if(abs((*pruned)[i].pdgId()) == 521){
    //if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600){
      const Candidate * bMeson = &(*pruned)[i];
      std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
      std::cout << "  found daugthers: " << std::endl;
      for(size_t j=0; j<packed->size();j++){
        int Ncalls=0;
        //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
        const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
        if(motherInPrunedCollection != nullptr && isAncestor( bMeson , motherInPrunedCollection, Ncalls)){
          std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << std::endl;
          std::cout << "           calls: " << Ncalls << std::endl;
        }
      }
    }

  }


}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODGenPartAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODGenPartAnalyzer::endJob() 
{
}

DEFINE_FWK_MODULE(MiniAODGenPartAnalyzer);
