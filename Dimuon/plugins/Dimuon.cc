// -*- C++ -*-
//
// Package:    GenStudy/Dimuon
// Class:      Dimuon
// 
/**\class Dimuon Dimuon.cc GenStudy/Dimuon/plugins/Dimuon.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Shawn Gregory Zaleski
//         Created:  Tue, 30 Jun 2015 14:01:36 GMT
// Edited by Ming Huang to look at lepton jets originating from dark photons
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_4/src/FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <TTree.h>
#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include <vector>
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class Dimuon : public edm::EDAnalyzer {
public:
  explicit Dimuon(const edm::ParameterSet&);
  ~Dimuon();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const reco::Candidate* getDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getLastDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getBoson( const reco::GenParticleCollection& genParts);
  const reco::Candidate* getMother(const reco::Candidate* part, int pid);
  //const reco::Candidate* getDYBoson(const reco::Candidate* part int pid)
  bool isBoson(int pid);
  bool isMuon(int pid);
  bool checkBosonStatus(const reco::GenParticleCollection& genParts);


  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  void gammavsSort(std::vector<const reco::Candidate*> fullGammavs, std::vector<const reco::Candidate*> &emptyGammavs,
			   std::vector<const reco::Candidate*> &eleSet, std::vector<const reco::Candidate*> &notGammavs,
			   std::vector<const reco::Candidate*> &fullGammavsRef);
  void notGammavsSort(std::vector<const reco::Candidate*> fullNotGammavs, std::vector<const reco::Candidate*> &emptyNotGammavs,
			      std::vector<const reco::Candidate*> &eleSet, std::vector<const reco::Candidate*> &gammavs,
			      std::vector<const reco::Candidate*> &fullNotGammavsRef);

  int leptonJetReco(std::vector<const reco::Candidate*> eles);
  const reco::Candidate* findBiggestPT(std::vector<const reco::Candidate*> leptons, std::vector<const reco::Candidate*> &smallerPTLeptons);

  struct P4Struct {
    float energy,et,eta,phi,pt,mass,theta;
    void fill(const math::XYZTLorentzVector& p4){
      if(p4.Pt()!=0 && p4.Et()!=0){
	energy = p4.E();
	et = p4.Et();
	eta = p4.Eta();
	phi = p4.Phi();
	pt =p4.Pt();
	mass = p4.mag();
	theta = p4.Theta();
      }else clear();
    }
    void clear(){energy=-999;et=-999;eta=-999;phi=-0;pt=-999;mass=-999;theta=-999;}
    static std::string contents(){return "energy/F:et:eta:phi:pt:mass:theta";}
  };

  TTree* tree_;


  TH1F *h_darkPhotonNum, *h_darkPhotonPT;
  TH1F *h_eleFromGammavNum, *h_eleFromGammavPT;
  TH1F *h_elePhi, *h_eleEta, *h_eleR, *h_eleInvariantMass;
  TH1F *h_eleDeltaPhi, *h_eleDeltaEta;
  TH1F *h_eleBiggestPT, *h_eleBiggestPTEta, *h_eleBiggestPTPhi;
  TH1F *h_neutralinoNum;
  TH1F *h_eleSetPT, *h_eleSetE, *h_eleSetSigmaPhi, *h_eleSetSigmaEta;
  TH1F *h_recoLeptonJetNum;
  //  TH1F *h_nuEleFromGammavNum, *h_nuEleFromGammavPT;

  //  TH1F *h_gammavMissingDaughters;
  //  TH1F *h_gammavExtraDaughters; // make sure the dark photon only has two daughters

  // TH1F * h_Zmass, *h_Zpt,*h_Zeta,*h_Zphi,*h_Zcharge;
  // TH1F *h_muMinusmass,*h_muMinuspt,*h_muMinuseta,*h_muMinusphi,*h_muMinuscharge;
  // TH1F *h_muPlusmass,*h_muPluspt,*h_muPluseta,*h_muPlusphi,*h_muPluscharge;
  // TH1F *h_dphi,*h_dtheta, *h_dr, *h_thetaMuMinus,*h_thetaMuPlus;
  // TH1F *h_massInvar, *h_dimuonPt, *h_dimuonEta, *h_dimuonPhi;
  // TH1F *h_cosTheta, *h_tanPhi, *h_csTheta, *h_csPhi;
  // TH1F *h_cosThetaPlusInvariantMass, *h_cosThetaMinusInvariantMass;

  // TH2F *h2_pt1_vs_pt2,*h2_eta1_vs_eta2,*h2_phi1_vs_phi2;

  P4Struct bosonP4_; // as a sanity check we have the right event...
  P4Struct muMinusP4_;
  P4Struct muPlusP4_;
  int muMinusPID_;
  int muPlusPID_;
  int bosonId_;
  double crossSec, cosTheta, tanPhi, csTheta, csPhi;
  double mCosThetaPlus, mCosThetaMinus;

  int debug_;
  edm::InputTag genPartsTag_;
  int decayParticlePID_;
  edm::InputTag genInfoProduct_;
  edm::EDGetTokenT<GenRunInfoProduct> genInfoProductToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;

  // ----------member data ---------------------------
  
};


void Dimuon::beginJob()
{
  edm::Service<TFileService> fs;
  h_darkPhotonNum = fs->make<TH1F>("darkPhotonNum", "Amount of dark photons per event", 15, -.5, 14.5);
  h_darkPhotonPT = fs->make<TH1F>("darkPhotonPT", "Dark Photon PT", 100, 0, 400);

  h_eleFromGammavNum = fs->make<TH1F>("eleFromGammavNum", "Number of electrons from dark photons in each event", 14, 0, 14);
  h_eleFromGammavPT = fs->make<TH1F>("eleFromGammavPT", "PT of electrons from dark photons", 100, 0, 400);
  h_elePhi = fs->make<TH1F>("elePhi", "Phi of electrons", 100, -3.2, 3.2);
  h_eleEta = fs->make<TH1F>("eleEta", "Eta of electrons", 100, -10, 10);
  h_eleInvariantMass = fs->make<TH1F>("eleInvariantMass", "Invariant mass of electrons", 100, -0.5, 1.5);

  h_eleDeltaPhi = fs->make<TH1F>("eleDeltaPhi", "difference in phi between two electrons that came from the same dark photon", 100, -10, 10);
  h_eleDeltaEta = fs->make<TH1F>("eleDeltaEta", "delta eta of electrons from the same dark photon", 100, -10, 10);
  h_eleR = fs->make<TH1F>("eleR", "distance", 100, -0.5, 10);

  h_eleBiggestPT = fs->make<TH1F>("eleBiggestPT", "PT of electron with biggest PT", 100, 0, 400);
  h_eleBiggestPTEta = fs->make<TH1F>("eleBiggestPTEta", "Eta of electron with biggest PT", 100, -10, 10);
  h_eleBiggestPTPhi = fs->make<TH1F>("eleBigestPTPhi", "Phi of electron with biggest PT", 100, -3.2, 3.2);

  h_neutralinoNum = fs->make<TH1F>("neutralinoNum", "Number of neutralinos per event", 14, 0, 14);

  h_eleSetPT = fs->make<TH1F>("eleSetPT", "The total PT of a set of all the electrons that came from the same neutralino", 100, 0, 400);
  h_eleSetE = fs->make<TH1F>("eleSetE", "The total energy of a set of all the electrons that came from the same neutralino", 100, 0, 400);
  h_eleSetSigmaPhi = fs->make<TH1F>("eleSetSigmaPhi", "The standard deviation of the phis in a set of all the electrons that came from the same neutralino", 100, -0.1, 3.2);
  h_eleSetSigmaEta = fs->make<TH1F>("eleSetSigmaEta", "The standard deviation of the etas in a set of all the electrons that came from the same neutralino", 150, -.5, 2.5);

  h_recoLeptonJetNum = fs->make<TH1F>("recoLeptonJetNum", "The number of lepton jets that there appear to be based off of only observing electrons and their angle", 15, -0.5, 14.5);

  //  h_nuEleFromGammavNum = fs->make<TH1F>("nuEleFromGammavNum", "Number of electron neutralinos from dark photons", 14, 0, 14);
  //  h_nuEleFromGammavPT = fs->make<TH1F>("nuEleFromGammavPT", "PT of electron neutralinos from dark photons", 100, 0, 400);

  //  h_gammavMissingDaughters = fs->make<TH1F>("gammavMissingDaughters", "Number of dark photon with less than two daughters", 14, 0, 14);
  //  h_gammavExtraDaughters = fs->make<TH1F>("gammavExtraDaughters", "Numbers of dark photons with more than two daughters", 14, 0, 14);

  // h_Zmass = fs->make<TH1F>("Zmass" , "m", 1000, 0., 600);
  // h_Zpt  = fs->make<TH1F>( "Zpt"  , "p_{t}", 500,  0., 2500. );
  // h_Zeta = fs->make<TH1F>( "Zeta" , "#eta" , 100, -10., 10.    );
  // h_Zphi = fs->make<TH1F>( "Zphi" , "#phi" , 100,  -3.20, 3.20   );
  // h_Zcharge = fs->make<TH1F>( "Zcharge" , "Q" ,3,  -1.5, 1.5    );
  // h_muMinusmass = fs->make<TH1F>("muMinusmass" , "m", 1000, 0., 500);
  // h_muMinuspt  = fs->make<TH1F>( "muMinuspt"  , "p_{t}", 500,  0., 2500. );
  // h_muMinuseta = fs->make<TH1F>( "muMinuseta" , "#eta" , 100, -5., 5.    );
  // h_muMinusphi = fs->make<TH1F>( "muMinusphi" , "#phi" , 100,  -3.15, 3.15   );
  // h_muMinuscharge = fs->make<TH1F>( "muMinuscharge" , "Q" ,3,  -1.5, 1.5    );

  // h_muPlusmass = fs->make<TH1F>("muPlusmass" , "m", 1000, 0., 500);
  // h_muPluspt  = fs->make<TH1F>( "muPluspt"  , "p_{t}", 500,  0., 2500. );
  // h_muPluseta = fs->make<TH1F>( "muPluseta" , "#eta" , 100, -5., 5.    );
  // h_muPlusphi = fs->make<TH1F>( "muPlusphi" , "#phi" , 100,  -3.15, 3.15   );
  // h_muPluscharge = fs->make<TH1F>( "muPluscharge" , "Q" ,3,  -1.5, 1.5    );

  // h_dphi = fs->make<TH1F>("delta phi", "#delta #phi", 100, -3.15, 3.15 );       
  // h_dtheta = fs->make<TH1F>("delta theta", "#delta #theta", 100, -3.15, 3.15); 
  // h_dr = fs->make<TH1F>("delta r", "#delta r", 100, 0, 10);
  // h_thetaMuMinus = fs->make<TH1F>("theta muMinus", "#theta", 100, -3.15, 3.15);      
  // h_thetaMuPlus = fs->make<TH1F>("theta muPlus", "#theta", 100, -3.15, 3.15); 
  // h_massInvar = fs->make<TH1F>("Invariant mass", "Invariant mass", 350, 0., 3500.);
  // h_dimuonPt = fs->make<TH1F>("Dimuon Pt", "Dimuon Pt", 500, 0, 2500);
  // h_dimuonEta = fs->make<TH1F>("Dimuon eta", "Dimuon #eta", 100, -5, 5);
  // h_dimuonPhi = fs->make<TH1F>("Dimuon Phi", "Dimuon #phi", 100, -3.15, 3.15);

  // h_cosTheta = fs->make<TH1F>("cosTheta", "cos #theta", 100, -1.01, 1.01);
  // h_tanPhi = fs->make<TH1F>("tanPhi", "tan #phi", 100, -1000.0, 1000.0);
  // h_csTheta = fs->make<TH1F>("csTheta", "#theta_{CS}", 100, -3.15, 3.15);
  // h_csPhi = fs->make<TH1F>("csPhi", "#phi_{CS}", 100, -3.15, 3.15);
  // h_cosThetaMinusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaMinus", "InvariantMass_cosThetaMinus", 350, 0., 3500.);
  // h_cosThetaPlusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaPlus", "InvariantMass_cosThetaPlus", 350, 0., 3500.);

  // h2_pt1_vs_pt2   = fs->make<TH2F>( "pt1_vs_pt2"   , "p_{t,1} vs. p_{t,2}"   , 500,  0., 2500., 500,  0., 2500.);
  // h2_eta1_vs_eta2 = fs->make<TH2F>( "eta1_vs_eta2" , "#eta_{1} vs. #eta_{2}" , 100, -5., 5.   , 100, -5., 5.   );
  // h2_phi1_vs_phi2 = fs->make<TH2F>( "phi1_vs_phi2" , "#phi_{1} vs. #phi_{2}" , 100,  -3.15, 3.15  , 100,  -3.15, 3.15  );

  tree_= fs->make<TTree>("pdfTree","PDF Tree");
  // tree_->Branch("evtId",&evtId_,EventId::contents().c_str());
  tree_->Branch("bosonP4",&bosonP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1P4",&muMinusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay2P4",&muPlusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1PID",&muMinusPID_,"decay1PID/I");
  tree_->Branch("decay2PID",&muPlusPID_,"decay2PID/I");
  tree_->Branch("bosonPID",&bosonId_,"bosonPID/I");
  tree_->Branch("crossSec", &crossSec, "crossSec/D");
  tree_->Branch("cosTheta", &cosTheta, "cosTheta/D");
  tree_->Branch("tanPhi", &tanPhi, "tanPhi/D");
  tree_->Branch("csTheta", &csTheta, "csTheta/D");
  tree_->Branch("csPhi", &csPhi, "csPhi/D");
  tree_->Branch("mCosThetaPlus", &mCosThetaPlus, "mCosThetaPlus/D");
  tree_->Branch("mCosThetaMinus", &mCosThetaMinus, "mCosThetaMinus/D");
  // tree_->Branch("pdfInfo",&pdfInfo_,PDFInfo::contents().c_str());
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
Dimuon::Dimuon(const edm::ParameterSet& iConfig)

{


  debug_=iConfig.getParameter<int>("debug");
  genPartsTag_=iConfig.getParameter<edm::InputTag>("genPartsTag");
  decayParticlePID_ = iConfig.getParameter<int>("decayParticlePID");
  genInfoProduct_ = iConfig.getParameter<edm::InputTag>("genInfoProduct");
  
  //now do what ever initialization is needed

  genInfoProductToken_ = consumes<GenRunInfoProduct,edm::InRun>(genInfoProduct_);
  genPartsToken_ = consumes<reco::GenParticleCollection>(genPartsTag_);
}


Dimuon::~Dimuon()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Dimuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genPartsHandle;
  iEvent.getByToken(genPartsToken_,genPartsHandle);
  const reco::GenParticleCollection& genParts = *genPartsHandle;

  bosonId_=0;
  bosonP4_.clear();
  muMinusP4_.clear();
  muPlusP4_.clear();
  muMinusPID_=0;
  muPlusPID_=0;

  //  const reco::Candidate* boson;
  //  const reco::Candidate* mother1;
  //  const reco::Candidate* mother2;
  //  const reco::Candidate* muMinus;
  //  const reco::Candidate* muPlus;
  math::XYZTLorentzVectorD dimuon;
  //  double dimuonPx, dimuonPy, dimuonPz, dimuonPt, pseudorapidity, Phi, mu1Energy, mu2Energy, dimuonQ;
  //  double thetaCos, thetaCS, phiTan, phiCS;
  //  double muPlusKPlus, muPlusKMinus, muMinusKPlus, muMinusKMinus, invariantK;

  double numDarkPhotons = 0;
  double numNeutralinos = 0;

  std::vector<const reco::Candidate*> eles;
  std::vector<const reco::Candidate*> recoEles;

  for(auto &part : genParts){
    if(abs(part.pdgId()) == 11 && part.numberOfDaughters() == 0){
      const reco::Candidate* mother = part.mother();
      if(mother->daughter(0)->pdgId() == part.pdgId()){
	recoEles.push_back(mother->daughter(0));
      }
      else if(mother->daughter(1)->pdgId() == part.pdgId()){
	recoEles.push_back(mother->daughter(1));
      }
    }
    if(part.pdgId() == 4900022){ // if particle is a dark photon
      numDarkPhotons++; // add one to the dark photon amount counter
      h_darkPhotonPT->Fill(part.pt()); // note the PT in histogram
      
      if(part.numberOfDaughters() == 2){ // if particle has two daughters
	
	// if both daughters are either electron or positron and they have no daughters thus they are final state
	if(abs(part.daughter(0)->pdgId()) == 11 && abs(part.daughter(1)->pdgId()) == 11 && part.daughter(0)->numberOfDaughters() == 0 && part.daughter(1)->numberOfDaughters() == 0){
	  h_eleFromGammavPT->Fill(part.daughter(0)->pt()); // enter the pt of the electrons
	  h_eleFromGammavPT->Fill(part.daughter(1)->pt());

	  eles.push_back(part.daughter(0)); // stick the electrons in a storage system for later use
	  eles.push_back(part.daughter(1));

	  h_elePhi->Fill(part.daughter(0)->phi()); // put the phi and eta of the two electrons into their histograms
	  h_eleEta->Fill(part.daughter(0)->eta());

	  h_elePhi->Fill(part.daughter(1)->phi());
	  h_eleEta->Fill(part.daughter(1)->eta());

	  const reco::Candidate* daughter1 = part.daughter(0); // invariant mass, delta phi, delta eta, delta r calculations
	  const reco::Candidate* daughter2 = part.daughter(1);

	  double invariantMass = sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi()))));
	  h_eleInvariantMass->Fill(invariantMass);

	  double deltaEta = daughter2->eta()-daughter1->eta();
	  double deltaPhi = daughter2->phi()-daughter1->phi();

	  h_eleDeltaPhi->Fill(deltaPhi);
	  h_eleDeltaEta->Fill(deltaEta);

	  double deltaR = std::sqrt((deltaEta)*(deltaEta)-(deltaPhi)*(deltaPhi));
	  h_eleR->Fill(deltaR);
	}
      }
    }

    // lepton groups analysis leading to reconstruction
    // more specifically sorting neutralino daughters as electron or not electron and then finding the standard deviation of the phi and eta for the electrons
    std::vector<const reco::Candidate*> gammavs;
    std::vector<const reco::Candidate*> gammavs2;
    std::vector<const reco::Candidate*> notGammavs;
    std::vector<const reco::Candidate*> notGammavs2;
    std::vector<const reco::Candidate*> eleSet;
    if(abs(part.pdgId() == 1000022)){ // if particle is a neutralino
      numNeutralinos++; // count neutralino number

      if(part.numberOfDaughters() == 2 && part.daughter(0)->pdgId() != 1000022 && part.daughter(1)->pdgId() != 1000022){ // if neutralino has two daughters and none of those are other neutralinos
	// sorting neutralino daughters as gammav or not gammav for further analysis
	if(part.daughter(0)->pdgId() == 4900022){ // if first daughter is a dark photon or not
	  gammavs.push_back(part.daughter(0));
	}
	else{
	  notGammavs.push_back(part.daughter(0));
	}
	if(part.daughter(1)->pdgId() == 4900022){ // if second daughter is a dark photon or not
	  gammavs.push_back(part.daughter(1));
	}
	else{
	  notGammavs.push_back(part.daughter(1));
	}
      }
    }
    while(!notGammavs.empty() || !gammavs.empty()){ // as long as there are non electrons that aren't sorted
      gammavsSort(gammavs, gammavs2, eleSet, notGammavs, gammavs); 
      gammavsSort(gammavs2, gammavs, eleSet, notGammavs, gammavs2);
      
      notGammavsSort(notGammavs, notGammavs2, eleSet, gammavs, notGammavs);
      notGammavsSort(notGammavs2, notGammavs, eleSet, gammavs, notGammavs2);
    }
    // total pt, total energy, of the set of electrons
    double eleSetPT = 0;
    double eleSetE = 0;
    for(auto &ele : eleSet){
      eleSetPT += ele->pt();
      eleSetE += ele->energy();
    }
    if(eleSetPT != 0){
      h_eleSetPT->Fill(eleSetPT);
      h_eleSetE->Fill(eleSetE);
    
      double eleSetSumPhi = 0;
      double eleSetSumEta = 0;
      // standard deviation of set of electrons
      for(auto &ele : eleSet){ // sum the phi/eta
	eleSetSumPhi += ele->phi();
	eleSetSumEta += ele->eta();
      }
      double eleSetMeanPhi = eleSetSumPhi/eleSet.size(); // get the mean of phi/eta
      double eleSetMeanEta = eleSetSumEta/eleSet.size();
      
      double phiMinusMeanSquaresSum = 0;
      double etaMinusMeanSquaresSum = 0;
      for(auto &ele:eleSet){ // get the phi/eta minus their mean and then square that and then add it to a running total (sum)
	phiMinusMeanSquaresSum += (ele->phi() - eleSetMeanPhi) * (ele->phi() - eleSetMeanPhi);
	etaMinusMeanSquaresSum += (ele->eta() - eleSetMeanEta) * (ele->eta() - eleSetMeanEta);
      }
      double phiMinusMeanSquaresMean = phiMinusMeanSquaresSum/eleSet.size(); // get the mean of the previous step
      double etaMinusMeanSquaresMean = etaMinusMeanSquaresSum/eleSet.size();
      
      h_eleSetSigmaPhi->Fill(std::sqrt(phiMinusMeanSquaresMean)); // get the standard deviation 
      h_eleSetSigmaEta->Fill(std::sqrt(etaMinusMeanSquaresMean));
    }
    h_darkPhotonNum->Fill(numDarkPhotons); // fill the number of dark photons for this event
    h_eleFromGammavNum->Fill(eles.size()); // fill the number of electrons for this event
    h_neutralinoNum->Fill(numNeutralinos); // fill the number of neutralinos for this event
  }

  // reconstructing electrons
  if(!recoEles.empty()){
    h_recoLeptonJetNum->Fill(leptonJetReco(recoEles));
  }
  else{
    std::cout << "No Electrons for this event" << std::endl;
  }

  const reco::Candidate* bigPTEle;
  double biggestElePt = 0;

  for(auto &ele : eles){ // for all electrons that came from dark photons in this event
    if(ele->pt() > biggestElePt){ // get the biggest pt electron, and enter all electron pts into histogram
      bigPTEle = ele;
      biggestElePt = ele->pt();
      
    }
  }
  if(biggestElePt != 0){ // get the pt, eta, phi for the electron with the biggest pt into histograms
    h_eleBiggestPT->Fill(bigPTEle->pt());
    h_eleBiggestPTEta->Fill(bigPTEle->eta());
    h_eleBiggestPTPhi->Fill(bigPTEle->phi());
  }


	// if((part.pdgId() == 1 || part.pdgId() == 2 || part.pdgId() == 3 || part.pdgId() == 4 || part.pdgId() == 5 || part.pdgId() == 6) && 
  // 	   (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){
  // 	  if(debug_ > 0){ std::cout << "\nFound the quark! " << "\nQuark is: " << part.pdgId() << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
  // 	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
  // 	    //	    mother1 = getMother(part.mother(0), 2212);
  // 	    //std::cout << "\nQuark mother is:" << mother1->pdgId() << std::endl;
  // 	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
  // 	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
  // 	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
  // 	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
  // 	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
  // 	  }
  // 	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
  // 	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
  // 	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
  // 		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
  // 	  boson = nullptr;
  // 	  if(!daughter1 || !daughter2){
  // 	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
  // 	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
  // 	  }
  // 	}

  // 	else if(part.pdgId() == 23 && (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){
  // 	  if(debug_ > 0){std::cout << "\nFound the Z boson! " << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
  // 	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
  // 	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
  // 	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
  // 	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
  // 	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
  // 	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
  // 	  }
  // 	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
  // 	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
  // 	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
  // 		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
  // 	  mother1 = &part;
  // 	  boson = mother1;
  // 	  if(!boson || !daughter1 || !daughter2){
  // 	    std::cout<<"boson::0x"<<std::hex<<boson<<std::dec<<std::endl;
  // 	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
  // 	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
  // 	  }

  // 	}
  //   }
      
  


  //   if(debug_ > 2){
  //     std::cout << "Eta of daughter1 is: " << daughter1->eta() << "\n";
  //     std::cout << "Eta of daughter2 is: " << daughter2->eta() << "\n";
  //   }

  //     if(boson){
  // 	bosonId_=boson->pdgId();
  // 	bosonP4_.fill(boson->p4());
	
  // 	h_Zmass->Fill(boson->mass());
  // 	h_Zpt->Fill(boson->pt());
  // 	h_Zeta ->Fill(boson->eta());
  // 	h_Zphi ->Fill(boson->phi());
  // 	h_Zcharge->Fill(boson->charge());
  //     }

  //   if(daughter1->charge() > 0 && daughter2->charge() < 0){
  //     muMinus = daughter2;
  //     muPlus = daughter1;
  //   }
  //   else if(daughter1->charge() < 0 && daughter2->charge() > 0){
  //     muMinus = daughter1;
  //     muPlus = daughter2;
  //   }

  //   else return;

  //   if(debug_ > 0){  
  //     std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
  // 	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
  //     std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
  // 	muPlus->phi() << "\tq = " << muPlus->charge();
  //   }
    
  // muMinusP4_.fill(muMinus->p4());
  //   muMinusPID_=muMinus->pdgId();
  //   if(debug_ > 0){  
  //     std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
  // 	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
  //   std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
  //     muPlus->phi() << "\tq = " << muPlus->charge();
  //   }


  //   h_muMinusmass->Fill(muMinus->mass());
  //   h_muMinuspt->Fill(muMinus->pt());
  //   h_muMinuseta->Fill(muMinus->eta());
  //   h_muMinusphi->Fill(muMinus->phi());
  //   h_muMinuscharge->Fill(muMinus->charge());
  //   h_thetaMuMinus->Fill(muMinus->theta());  
  
  //   muPlusP4_.fill(muPlus->p4());
  //   muPlusPID_=muPlus->pdgId();

  //   h_muPlusmass->Fill(muPlus->mass());
  //   h_muPluspt->Fill(muPlus->pt());
  //   h_muPluseta->Fill(muPlus->eta());
  //   h_muPlusphi->Fill(muPlus->phi());
  //   h_muPluscharge->Fill(muPlus->charge());
  //   h_thetaMuPlus->Fill(muPlus->theta());    

  //   muPlusKPlus = (1/sqrt(2))*(muPlus->energy() + muPlus->pz());
  //   muPlusKMinus = (1/sqrt(2))*(muPlus->energy() - muPlus->pz());
  //   muMinusKPlus = (1/sqrt(2))*(muMinus->energy() + muMinus->pz());
  //   muMinusKMinus = (1/sqrt(2))*(muMinus->energy() - muMinus->pz());

  //   invariantK = (muPlusKPlus*muMinusKMinus - muMinusKPlus*muPlusKMinus);
  //   std::cout << "\n\nInvariantK is: " << invariantK << std::endl;

  //   dimuon = muMinus->p4() + muPlus->p4();

  //   dimuonPt =dimuon.pt();
  //   dimuonPz = dimuon.pz();
  //   pseudorapidity = asinh(dimuonPz/dimuonPt);
  //   dimuonPx = dimuon.px();
  //   dimuonPy = dimuon.py();
  //   Phi = acos(dimuonPx/dimuonPt);
  //   dimuonQ = sqrt(pow(dimuon.energy(),2) - pow(dimuon.pt(),2) - pow(dimuon.pz(),2));
  //   std::cout << "\n\nDimuon Energy is: " << dimuon.energy() << std::endl << std::endl;
    
  //   double denominatorTheta, denominatorPhi1, denominatorPhi2, numeratorPhi1, numeratorPhi2;
  //   double denominatorPhi, numeratorPhi;
  //   double deltaX, deltaY;
  //   double invariantMass;

  //   denominatorTheta = dimuonQ*sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
  //   thetaCos = (dimuon.pz()/fabs(dimuon.pz()))*(2/denominatorTheta)*invariantK;
  //   thetaCS = acos(thetaCos);

  //   denominatorPhi1 = dimuonQ*dimuon.pt();
  //   numeratorPhi1 = sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
  //   deltaX = muPlus->px() - muMinus->px();
  //   deltaY = muPlus->py() - muMinus->py();
  //   denominatorPhi2 = ((deltaX*dimuon.px()) + (deltaY*dimuon.py()));
  //   numeratorPhi2 = ((deltaX*dimuon.py()) + (deltaY*dimuon.px()));
  //   numeratorPhi = numeratorPhi1*numeratorPhi2;
  //   denominatorPhi = denominatorPhi1 * denominatorPhi2;

  //   phiTan = numeratorPhi/denominatorPhi;
  //   phiCS = atan(phiTan);


  //   mu1Energy = muPlus->energy();
  //   mu2Energy = muMinus->energy();
  //   std::cout << "\n\nmuon Energies are: " << mu1Energy << "__" << mu2Energy << std::endl << std::endl;
  //   std::cout << "\ndimuon px_py_pz are: "<< dimuonPx << "_" << dimuonPy << "_" << dimuonPz << std::endl;
 
  //   cosTheta=thetaCos;
  //   tanPhi=phiTan;
  //   csTheta=thetaCS;
  //   csPhi=phiCS;


  //   h_cosTheta->Fill(thetaCos);
  //   h_csTheta->Fill(thetaCS);
  //   h_tanPhi->Fill(phiTan);
  //   h_csPhi->Fill(phiCS);

    

  //   std::cout << "\n\n\ncos(Theta_CS) = " << thetaCos << "\tThetaCS = " << thetaCS << std::endl;
  //   std::cout << "\n\n\nTan(phi_CS) = " << phiTan << "\tPhiCS = " << phiCS << std::endl;

  //   invariantMass = sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi()))));


  //   if(thetaCos < 0.0){
  //     h_cosThetaMinusInvariantMass->Fill(invariantMass);
  //     mCosThetaMinus = invariantMass;
  //   }
  //   else{
  //     h_cosThetaPlusInvariantMass->Fill(invariantMass);
  //     mCosThetaPlus = invariantMass;
  //   }


  //   h_dphi->Fill(TVector2::Phi_mpi_pi(muMinus->phi()- muPlus->phi()));
  //   h_dtheta->Fill(TVector2::Phi_mpi_pi(muMinus->theta()- muPlus->theta()));
  //   h_dr->Fill(reco::deltaR(muMinus->p4(),muPlus->p4()));
  //   h_massInvar->Fill(sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi())))));
  //   h_dimuonPt->Fill(dimuonPt);
  //   h_dimuonEta->Fill(pseudorapidity);
  //   h_dimuonPhi->Fill(Phi);
  //   h2_phi1_vs_phi2->Fill(muMinus->phi(),muPlus->phi());  
  //   h2_eta1_vs_eta2->Fill(muMinus->eta(),muPlus->eta());
  //   h2_pt1_vs_pt2->Fill(muMinus->pt(),muPlus->pt());


    //  }

  std::cout << "\n\n===========================================================================================================" << std::endl;
  tree_->Fill();  
}

void Dimuon::gammavsSort(std::vector<const reco::Candidate*> fullGammavs, std::vector<const reco::Candidate*> &emptyGammavs,
			 std::vector<const reco::Candidate*> &eleSet, std::vector<const reco::Candidate*> &notGammavs,
			 std::vector<const reco::Candidate*> &fullGammavsRef){
  for(auto &gammav : fullGammavs){
    if(abs(gammav->daughter(0)->pdgId()) == 11){ // first daughter sorting for gammavs
      eleSet.push_back(gammav->daughter(0));
    }
    else if(abs(gammav->daughter(0)->pdgId()) == 4900022){
      emptyGammavs.push_back(gammav->daughter(0));
    }
    else{
      notGammavs.push_back(gammav->daughter(0));
    }
    
    if(abs(gammav->daughter(1)->pdgId()) == 11){ // second daughter sorting for gammavs
      eleSet.push_back(gammav->daughter(1));
    }
    else if(abs(gammav->daughter(1)->pdgId()) == 4900022){
      emptyGammavs.push_back(gammav->daughter(1));
    }
    else{
      notGammavs.push_back(gammav->daughter(1));
    }  
  }
  fullGammavsRef.clear();
}

void Dimuon::notGammavsSort(std::vector<const reco::Candidate*> fullNotGammavs, std::vector<const reco::Candidate*> &emptyNotGammavs,
			    std::vector<const reco::Candidate*> &eleSet, std::vector<const reco::Candidate*> &gammavs,
			    std::vector<const reco::Candidate*> &fullNotGammavsRef){
  for(auto &notGammav : fullNotGammavs){
    if(notGammav->pdgId() == 4900002 || notGammav->pdgId() == 4900004){
      continue;
    }
    else{
      if(notGammav->numberOfDaughters() == 2){
	if(abs(notGammav->daughter(0)->pdgId()) == 4900022){
	  gammavs.push_back(notGammav->daughter(0));
	}
	else{
	  emptyNotGammavs.push_back(notGammav->daughter(0));
	}
	if(abs(notGammav->daughter(1)->pdgId()) == 4900022){
	  gammavs.push_back(notGammav->daughter(1));
	}
	else{
	  emptyNotGammavs.push_back(notGammav->daughter(1));
	}
      }
      else{
	std::cout << "Particle " << notGammav->pdgId() << " doesn't have 2 daughter" << std::endl;
      }
    }
  }
  fullNotGammavsRef.clear();
}


int Dimuon::leptonJetReco(std::vector<const reco::Candidate*> leptons){
  double leptonJets = 0;
  double deltaRCutoff = 0.5;

  // find the lepton with the highest pt
  std::vector<const reco::Candidate*> smallerPTLeptons;
  std::vector<const reco::Candidate*> smallerPTLeptons2;
  const reco::Candidate* biggestPTLepton = findBiggestPT(leptons, smallerPTLeptons); 
  

  // find the other leptons that are in the same jet
  while(!smallerPTLeptons.empty()){
    bool isJet = false;
    for(auto &lepton : smallerPTLeptons){
      double deltaEta = biggestPTLepton->eta() - lepton->eta(); // find deltar between the biggest pt lepton and the other leptons
      double deltaPhi = biggestPTLepton->phi() - lepton->phi();
      
      double deltaR = std::sqrt((deltaEta)*(deltaEta)-(deltaPhi)*(deltaPhi));

      if(deltaR > deltaRCutoff){
	smallerPTLeptons2.push_back(lepton);
      }
      else{
	isJet = true;
      }
    }
    if(isJet){
      leptonJets++;
    }
    smallerPTLeptons.clear();

    biggestPTLepton = findBiggestPT(smallerPTLeptons2, smallerPTLeptons);
    smallerPTLeptons2.clear();
  }
    
  return leptonJets;
}

const reco::Candidate* Dimuon::findBiggestPT(std::vector<const reco::Candidate*> leptons, std::vector<const reco::Candidate*> &smallerPTLeptons){
  double biggestPT = 0;
  const reco::Candidate* biggestPTLepton;

  for(auto &lepton : leptons){
    if(lepton->pt() > biggestPT){
      if(biggestPT != 0){
	smallerPTLeptons.push_back(biggestPTLepton);
      }
      biggestPTLepton = lepton;
      biggestPT = lepton->pt();
    }
    else{
      smallerPTLeptons.push_back(lepton);
    }
  }
  return biggestPTLepton;
}

bool Dimuon::isBoson(int pid)
{
  if(pid==23 || abs(pid)==22 || pid==32){
    if(debug_ > 0) std::cout << "\n\nFound Boson\n";
    return true;
  }
  else return false;
}

bool Dimuon::isMuon(int pid){
  if(abs(pid)==11 || abs(pid) ==13){
    if(debug_ > 0) std::cout << "\n\nFound A Muon!\n";
    return true;
  }
  else return false;
}

bool Dimuon::checkBosonStatus( const reco::GenParticleCollection& genParts){
  const reco::Candidate* boson = getBoson(genParts);
  if(boson == nullptr){
    if(debug_ > 0) std::cout << "\nBoson is: "  << boson;
    return false;
  }

  else if( boson->status() != 22){
    if(debug_ > 0)  std::cout <<"\nBoson Status is: "<< boson->status();
    return false;
  }
 
    return true;
 }

const reco::Candidate* Dimuon::getBoson( const reco::GenParticleCollection& genParts)
{
  for(auto &part : genParts){
    if(isBoson(part.pdgId())){
      if(debug_ > 1){
      std::cout << "\npId is: " << part.pdgId();
      std::cout << "\nStatus is: " << part.status();
      }
      return getLastDaughter(&part,part.pdgId());
    }
  }
  return nullptr;
}


const reco::Candidate* Dimuon::getLastDaughter(const reco::Candidate* part,int pid)
{
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return getLastDaughter(part->daughter(partNr),pid);
  }
  return part;
}
       
const reco::Candidate* Dimuon::getDaughter(const reco::Candidate* part,int pid)
{  
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return part->daughter(partNr);
  }
  return nullptr;
}

 const reco::Candidate* Dimuon::getMother(const reco::Candidate* part, int pid)
{
  for(size_t partNr = 0; part && partNr < part->numberOfMothers(); partNr++){
    if(part->mother(partNr)->pdgId() == pid) return getMother(part->mother(partNr),pid);
  
    else if(abs(part->mother(partNr)->pdgId()) == 1 || abs(part->mother(partNr)->pdgId()) == 2 ||
	    abs(part->mother(partNr)->pdgId()) == 3 || abs(part->mother(partNr)->pdgId()) == 4 ||
		 abs(part->mother(partNr)->pdgId()) == 5 || abs(part->mother(partNr)->pdgId()) == 6 ||
	    abs(part->mother(partNr)->pdgId()) == 7 || abs(part->mother(partNr)->pdgId()) == 8 ||
	    abs(part->mother(partNr)->pdgId()) == 23 || abs(part->mother(partNr)->pdgId()) == 32  || 
	    abs(part->mother(partNr)->pdgId()) == 22) return part->mother(partNr);
  }  
  return nullptr;
  
}


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
Dimuon::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  
/*void 
Dimuon::beginRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
{
   edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByToken(genInfoProductToken_, genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();
  //  tree_->Fill();
  std::cout<< "Cross Section is: "  << crossSec << std::endl;  
 

}
*/

// ------------ method called when ending the processing of a run  ------------

void 
  Dimuon::endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
  {
   edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByToken(genInfoProductToken_, genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();
  std::cout<< "Cross Section is: "  << crossSec << std::endl;  
 
  }

  
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  Dimuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  Dimuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Dimuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Dimuon);
