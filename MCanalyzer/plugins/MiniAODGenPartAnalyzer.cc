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

 //#include "TLorentzVector.h"
 //#include "TTree.h"
 //#include "Math/GenVector/Boost.h"
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

class MiniAODGenPartAnalyzer : public edm::EDAnalyzer {
  public:
    explicit MiniAODGenPartAnalyzer(const edm::ParameterSet&);
    ~MiniAODGenPartAnalyzer();
    bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle, int &calls);



  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    TLorentzVector gen_b_p4, gen_phi_p4, gen_kaon_p4, gen_muon1_p4, gen_muon2_p4, gen_gamma1_p4, gen_gamma2_p4;
    TLorentzVector gen_b_p4CM, gen_phi_p4CM, gen_kaon_p4CM, gen_muon1_p4CM, gen_muon2_p4CM, gen_gamma1_p4CM, gen_gamma2_p4CM;
    TTree*         tree_;

    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    float costhetaL, costhetaKL;
};

MiniAODGenPartAnalyzer::MiniAODGenPartAnalyzer(const edm::ParameterSet& iConfig):
prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
costhetaL(-2.0), costhetaKL(-2.0)
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

  bool debug = true;
  bool iskaon, isMuon1, isMuon2;


  gen_b_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_kaon_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon1_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon2_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma1_p4.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma2_p4.SetPxPyPzE(0.,0.,0.,0.);

  gen_b_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_kaon_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon1_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_muon2_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma1_p4CM.SetPxPyPzE(0.,0.,0.,0.);
  gen_gamma2_p4CM.SetPxPyPzE(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);

  //let's try to find all status1 originating directly from a B meson decay 

  for(size_t i=0; i<pruned->size();i++){
    iskaon=false;
    isMuon2=false;
    isMuon1=false;
    if(abs((*pruned)[i].pdgId()) == 521){
    //if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600){
      const Candidate * bMeson = &(*pruned)[i];
      gen_b_p4.SetPtEtaPhiM(bMeson->pt(),bMeson->eta(),bMeson->phi(),bMeson->mass());

      if (debug){
        std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
        std::cout << "  found daugthers: " << std::endl;
      }
      
      for(size_t j=0; j<packed->size();j++){
        int Ncalls=0;
        //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
        const Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
        //Lets check if the daughter particles are K+(321) and muons(13)
        if(motherInPrunedCollection != nullptr && isAncestor( bMeson , motherInPrunedCollection, Ncalls)){
          if (abs((*packed)[j].pdgId()) == 321){
            iskaon=true;
            gen_kaon_p4.SetPtEtaPhiM((*packed)[j].pt(),(*packed)[j].eta(),(*packed)[j].phi(),(*packed)[j].mass());
          }
          else if(abs((*packed)[j].pdgId()) == 13){
            if((*packed)[j].pdgId()*bMeson->pdgId()<0){
              isMuon1=true;
              gen_muon1_p4.SetPtEtaPhiM((*packed)[j].pt(),(*packed)[j].eta(),(*packed)[j].phi(),(*packed)[j].mass());        
            }
            else {
              isMuon2=true;
              gen_muon2_p4.SetPtEtaPhiM((*packed)[j].pt(),(*packed)[j].eta(),(*packed)[j].phi(),(*packed)[j].mass()); 
            }
          }
          std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << std::endl;
          std::cout << "           calls: " << Ncalls << std::endl;
        }
      }

      if (iskaon && isMuon1 && isMuon2){
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
      }
    }

  }
 tree_->Fill();
 
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODGenPartAnalyzer::beginJob()
{
  std::cout << "Beginning analyzer job" << std::endl;

  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("ntuple","B+->K+ mu mu ntuple");

  tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
  tree_->Branch("gen_kaon_p4",  "TLorentzVector",  &gen_kaon_p4);
  tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
  tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
  tree_->Branch("gen_gamma1_p4",  "TLorentzVector",  &gen_gamma1_p4);
  tree_->Branch("gen_gamma2_p4",  "TLorentzVector",  &gen_gamma2_p4);

  tree_->Branch("gen_b_p4CM",     "TLorentzVector",  &gen_b_p4CM);
  tree_->Branch("gen_kaon_p4CM",  "TLorentzVector",  &gen_kaon_p4CM);
  tree_->Branch("gen_muon1_p4CM",  "TLorentzVector",  &gen_muon1_p4CM);
  tree_->Branch("gen_muon2_p4CM",  "TLorentzVector",  &gen_muon2_p4CM);
  tree_->Branch("gen_gamma1_p4CM",  "TLorentzVector",  &gen_gamma2_p4CM);
  tree_->Branch("gen_gamma2_p4CM",  "TLorentzVector",  &gen_gamma2_p4CM);
  

  //tree_->Branch("daughter_id",   "vector", &daughter_id);

  //tree_->Branch("number_daughters",  &number_daughters);
  tree_->Branch("costhetaL",  &costhetaL);
  tree_->Branch("costhetaKL",  &costhetaKL);

  //tree_->Branch("Nbplus", &bplus);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODGenPartAnalyzer::endJob() 
{
}

DEFINE_FWK_MODULE(MiniAODGenPartAnalyzer);
