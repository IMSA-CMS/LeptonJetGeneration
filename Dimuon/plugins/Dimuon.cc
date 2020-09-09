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
#include "DataFormats/Math/interface/normalizedPhi.h"
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

  const void printGenParticleCollection(const reco::GenParticleCollection& genParts);
  const int getIndexOf(const reco::Candidate* part, const reco::GenParticleCollection& genParts);

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

  // Neutralino and dark photon histograms
  TH1F *h_neutralinoNum, *h_chiR;
  TH1F *h_darkPhotonNum, *h_darkPhotonPT;

  // Electron analysis histograms
  TH1F *h_eleFromGammavNum, *h_eleFromGammavPT;
  TH1F *h_elePhi, *h_eleEta, *h_eleR, *h_eleInvariantMass;
  TH1F *h_eleDeltaPhi, *h_eleDeltaEta;

  // High pT electron histograms
  TH1F *h_30GevOrMoreLeptons;
  TH1F *h_eleBiggestPT, *h_eleBiggestPTEta, *h_eleBiggestPTPhi;

  // Electron set analysis histograms
  TH1F *h_eleSetPT, *h_eleSetE, *h_eleSetSigmaPhi, *h_eleSetSigmaEta;

  // Reconstruction histograms
  TH1F *h_recoLeptonJetNum, *h_actualJetNum, *h_jetNumDiff;
  TH1F *h_recoLeptonNumPerEvent, *h_recoLeptonNumInOneJetEvents;
  TH1F *h_notInJetLepton;

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
  // Neutralino histograms
  h_neutralinoNum = fs->make<TH1F>("neutralinoNum", "Number of neutralinos per event", 14, 0, 14);
  h_chiR = fs->make<TH1F>("chiR", "delta r of the two neutralinos in the event", 100, 0, 10);

  // Dark photon histograms
  h_darkPhotonNum = fs->make<TH1F>("darkPhotonNum", "Amount of dark photons per event", 15, -.5, 14.5);
  h_darkPhotonPT = fs->make<TH1F>("darkPhotonPT", "Dark Photon PT", 100, 0, 400);

  // Electrons that came from a dark photon histograms
  h_eleFromGammavNum = fs->make<TH1F>("eleFromGammavNum", "Number of electrons from dark photons in each event", 14, 0, 14);
  h_eleFromGammavPT = fs->make<TH1F>("eleFromGammavPT", "PT of electrons from dark photons", 100, 0, 400);

  // Electron histograms
  h_elePhi = fs->make<TH1F>("elePhi", "Phi of electrons", 100, -3.2, 3.2);
  h_eleEta = fs->make<TH1F>("eleEta", "Eta of electrons", 100, -10, 10);
  h_eleInvariantMass = fs->make<TH1F>("eleInvariantMass", "Invariant mass of electrons", 100, -0.5, 1.5);
  h_eleDeltaPhi = fs->make<TH1F>("eleDeltaPhi", "difference in phi between two electrons that came from the same dark photon", 100, -10, 10);
  h_eleDeltaEta = fs->make<TH1F>("eleDeltaEta", "delta eta of electrons from the same dark photon", 100, -10, 10);
  h_eleR = fs->make<TH1F>("eleR", "Distance between a pair of two electrons that originated from the same dark photon in phi-eta space", 300, -0.1, 2.9);

  // High pT electron histograms
  h_30GevOrMoreLeptons = fs->make<TH1F>("30GevOrMoreLeptons", "The number of events that have a lepton 30 gev or greater", 2, 0, 2);
  h_eleBiggestPT = fs->make<TH1F>("eleBiggestPT", "PT of the electron with the biggest PT in its event", 100, 0, 400);
  h_eleBiggestPTEta = fs->make<TH1F>("eleBiggestPTEta", "Eta of electron with biggest PT", 100, -10, 10);
  h_eleBiggestPTPhi = fs->make<TH1F>("eleBigestPTPhi", "Phi of electron with biggest PT", 100, -3.2, 3.2);

  // Set of electrons that came from the same neutralino histograms
  h_eleSetPT = fs->make<TH1F>("eleSetPT", "The total PT of a set of all the electrons that came from the same neutralino", 100, 0, 400);
  h_eleSetE = fs->make<TH1F>("eleSetE", "The total energy of a set of all the electrons that came from the same neutralino", 100, 0, 400);
  h_eleSetSigmaPhi = fs->make<TH1F>("eleSetSigmaPhi", "The standard deviation of the phis in a set of all the electrons that came from the same neutralino", 100, -0.1, 3.2);
  h_eleSetSigmaEta = fs->make<TH1F>("eleSetSigmaEta", "The standard deviation of the etas in a set of all the electrons that came from the same neutralino", 150, -.5, 2.5);

  // Reconstruction and reconstruction checking histograms
  h_recoLeptonJetNum = fs->make<TH1F>("recoLeptonJetNum", "The number of lepton jets that there appear to be based off of only observing electrons and their angle", 15, -0.5, 14.5);
  h_actualJetNum = fs->make<TH1F>("actualJetNum", "The actual number of lepton jets in the event, based on whether neutrinos produced electrons", 15, -0.5, 14.5);
  h_jetNumDiff = fs->make<TH1F>("jetNumDiff", "The difference between the number of lepton jets reconstructed, and how many lepton jets there actually were", 5, 0, 5);
  h_recoLeptonNumPerEvent = fs->make<TH1F>("recoLeptonNumPerEvent", "The number of leptons being reconstructed in an event", 15, 0, 15);
  h_recoLeptonNumInOneJetEvents = fs->make<TH1F>("recoLeptonNumInOneJetEvents", "The number of leptons in events where one jet is reconstructed", 15, 0, 15);
  h_notInJetLepton = fs->make<TH1F>("notInJetLepton", "The number of leptons that didn't get put into a jet in each event", 15, 0, 15);
 
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

  math::XYZTLorentzVectorD dimuon;

  double numDarkPhotons = 0;
  double numNeutralinos = 0;

  std::vector<const reco::Candidate*> eles;   // Analysis electron container
  std::vector<const reco::Candidate*> recoEles;  // Reconstruction electron container

  bool thirtyGevOrMoreLeptonExists = false;

  const reco::Candidate* chi1;   // Variable to place first neutralino while looking for the second
  bool chi1Initialized = false;

  int actualJetNum = 0;

  for(auto &part : genParts) // for every particle (part) in this event (genParts)
  {
    // Neutralino delta r calculations
    if(abs(part.pdgId()) == 1000022 && part.numberOfDaughters() > 1) // if particle (part) is a neutralino and particle (part) has more than one daughter 
    {
      if(chi1Initialized)      // if the first neutralino has been found and stored
      {
        ROOT::Math::LorentzVector d1 = chi1->p4();
	const reco::Candidate* chi2 = &part;
        ROOT::Math::LorentzVector d2 = chi2->p4();
	h_chiR->Fill(reco::deltaR(d1, d2));
      }
      else // else save the first neutralino
      {
	chi1 = &part;
	chi1Initialized = true;
      }
    }

    // Counting events with electrons with a Gev >=30
    if(part.numberOfDaughters() == 0 && abs(part.pdgId()) == 11 && part.pt() >= 30)
    {
      thirtyGevOrMoreLeptonExists = true;
    }
    
    // Filling reconstruction electrons container
    if(abs(part.pdgId()) == 11 && part.numberOfDaughters() == 0)
    {
      recoEles.push_back(&part);
    }

    // Electron and dark photon analysis
    if(part.pdgId() == 4900022) // if particle (part) is a dark photon
    {      
      if(part.numberOfDaughters() == 2) // if particle (part) has two daughters
      {
	numDarkPhotons++;
	h_darkPhotonPT->Fill(part.pt()); // Enter pT into histogram
	
	const reco::Candidate* daughter1; // Container for final state electrons
	const reco::Candidate* daughter2;

	if(abs(part.daughter(0)->pdgId()) == 11 && abs(part.daughter(1)->pdgId()) == 11) // if both daughters of the dark photon (part) are electrons
	{
	  bool stillLooking = true;
	  const reco::Candidate* checkNext; // Container for the next particle to be examined
	  
	  checkNext = part.daughter(0); // Set the first particle to be checked as first daughter of dark photon
	  while(stillLooking) // while the final state electron coming from the first daughter hasn't been found
	  {
	    if(checkNext->numberOfDaughters() == 0) // if the particle being checked has no daughters
	    {
	      daughter1 = checkNext; // store this particle as the final state of the first daughter
	      stillLooking = false;
	    }
	    else // else check the daughter of the currently being checked particle, next
	    {
	      checkNext = checkNext->daughter(0);
	    }
	  }

	  // Finding final state particle of the second daughter of the dark photon
	  stillLooking = true;
	  checkNext = part.daughter(1);
	  while(stillLooking)
	  {
	    if(checkNext->numberOfDaughters() == 0)
	    {
	      daughter2 = checkNext;
	      stillLooking = false;
	    }
	    else
	    {
	      checkNext = checkNext->daughter(0);
	    }
	  }

	  // Calculations involving the two final state daughters of the dark photon
	  h_eleFromGammavPT->Fill(daughter1->pt()); // Filling pT of electrons
	  h_eleFromGammavPT->Fill(daughter2->pt());
	  
	  eles.push_back(daughter1); // Save the electrons for further analysis later
	  eles.push_back(daughter2);
	  
	  h_elePhi->Fill(daughter1->phi()); // Filling phi and eta for both electrons
	  h_eleEta->Fill(daughter1->eta());
	  
	  h_elePhi->Fill(daughter2->phi());
	  h_eleEta->Fill(daughter2->eta());
	  
	  double invariantMass = sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi()))));
	  h_eleInvariantMass->Fill(invariantMass);

	  double deltaEta = daughter2->eta()-daughter1->eta();
	  double deltaPhi = reco::deltaPhi(daughter2->phi(), daughter1->phi());
	  
	  h_eleDeltaPhi->Fill(deltaPhi);
	  h_eleDeltaEta->Fill(deltaEta);

	  ROOT::Math::LorentzVector d1 = daughter1->p4();
	  ROOT::Math::LorentzVector d2 = daughter2->p4();
	  
	  //double deltaR = std::sqrt((deltaEta)*(deltaEta)-(deltaPhi)*(deltaPhi));
	  h_eleR->Fill(reco::deltaR(d1, d2));
	}
	else
	{
	  std::cout << "Dark photon daughters are not both electrons" << std::endl;
	}
      }
      else
      {
	std::cout << "Dark photon doesn't have 2 daughters" << std::endl;
      }
    }

    // sorting neutralino daughters as electron or not electron and then finding the standard deviation of the phi and eta for the electrons
    std::vector<const reco::Candidate*> gammavs;
    std::vector<const reco::Candidate*> gammavs2;
    std::vector<const reco::Candidate*> notGammavs;
    std::vector<const reco::Candidate*> notGammavs2;
    std::vector<const reco::Candidate*> eleSet;
    if(abs(part.pdgId()) == 1000022) // if particle is a neutralino
    {

      if(part.numberOfDaughters() == 2 && part.daughter(0)->pdgId() != 1000022 && part.daughter(1)->pdgId() != 1000022) // if neutralino has two daughters and none of those are other neutralinos
      {
	//numNeutralinos++; // count neutralino number
	// sorting neutralino daughters as gammav or not gammav for further analysis
	if(part.daughter(0)->pdgId() == 4900022) // if first daughter is a dark photon or not
	{	  
	  gammavs.push_back(part.daughter(0));
	}
	else
	{
	  notGammavs.push_back(part.daughter(0));
	}
	if(part.daughter(1)->pdgId() == 4900022) // if second daughter is a dark photon or not
	{
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

    // starting from the neutralino, find out how many lepton jets were produced
    // also store all the electrons from this jet into a collection for analysis
    if(abs(part.pdgId()) == 1000022) // if particle is a neutralino
    {
      if(part.numberOfDaughters() == 2 && part.daughter(0)->pdgId() != 1000022 && part.daughter(1)->pdgId() != 1000022) // if neutralino has two daughters and none of those are other neutralinos
      {
	numNeutralinos++;
	// std::cout << "First daughter is " << part.daughter(0)->pdgId() << std::endl;
	// std::cout << "Second daughter is " << part.daughter(1)->pdgId() << std::endl;
	
	if(part.daughter(0)->pdgId() != 4900002 && part.daughter(1)->pdgId() != 4900002)
	{
	  actualJetNum++;
	}
	else if(part.daughter(0)->pdgId() ==4900002 || part.daughter(1)->pdgId() == 4900002)
	{
	  std::cout << "First daughter is " << part.daughter(0)->pdgId() << std::endl;
	  std::cout << "Second daughter is " << part.daughter(1)->pdgId() << std::endl;
	}
      }
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
  }
  
  h_actualJetNum->Fill(actualJetNum);
  if(thirtyGevOrMoreLeptonExists)
  {
    h_30GevOrMoreLeptons->Fill(1);
  }
  if(numDarkPhotons != 0)
  {
    h_darkPhotonNum->Fill(numDarkPhotons); // fill the number of dark photons for this event
  }
  else
  {
    std::cout << "No dark photons in this event" << std::endl;
  }
  if(eles.size() != 0)
  {
    h_eleFromGammavNum->Fill(eles.size()); // fill the number of electrons for this event
  }
  else
  {
    std::cout << "No electrons in this event" << std::endl;
  }
  if(numNeutralinos != 0)
  {
    h_neutralinoNum->Fill(numNeutralinos); // fill the number of neutralinos for this event
  }
  else
  {
    std::cout << "No neutralinos in this event" << std::endl;
  }

  // reconstructing electrons
  if(!recoEles.empty()){
    int recoJetNum = leptonJetReco(recoEles);
    h_recoLeptonJetNum->Fill(recoJetNum);
    h_jetNumDiff->Fill(abs(recoJetNum-actualJetNum));
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

  //printGenParticleCollection(genParts);
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
  double deltaRCutoff = 0.2;

  int notInJetLeptonNum = 0;
  int originalLeptonNum = leptons.size();

  h_recoLeptonNumPerEvent->Fill(leptons.size());

  // find the lepton with the highest pt
  std::vector<const reco::Candidate*> smallerPTLeptons;
  std::vector<const reco::Candidate*> notInJetLeptons;
  const reco::Candidate* biggestPTLepton = findBiggestPT(leptons, smallerPTLeptons); 
  

  // find the other leptons that are in the same jet
  while(!smallerPTLeptons.empty()){
    std::cout << "size of smallerPTLeptons " << smallerPTLeptons.size() << std::endl;

    bool isJet = false;
    for(auto &lepton : smallerPTLeptons){
      // double deltaEta = biggestPTLepton->eta() - lepton->eta(); // find deltar between the biggest pt lepton and the other leptons
      // double deltaPhi = biggestPTLepton->phi() - lepton->phi();
      
      // double deltaR = std::sqrt((deltaEta)*(deltaEta)-(deltaPhi)*(deltaPhi));
      // checking github update
      ROOT::Math::LorentzVector e1 = biggestPTLepton->p4();
      ROOT::Math::LorentzVector e2 = lepton->p4();
      double deltaR = reco::deltaR(e1, e2);

      if(deltaR > deltaRCutoff){
	notInJetLeptons.push_back(lepton);
      }
      else{
	isJet = true;
      }
    }
    if(isJet){
      leptonJets++;
    }
    else{
      notInJetLeptonNum++;
      std::cout << "not a jet with " << smallerPTLeptons.size() << " electrons observed" << std::endl;
    }
    smallerPTLeptons.clear();

    biggestPTLepton = findBiggestPT(notInJetLeptons, smallerPTLeptons);
    notInJetLeptons.clear();
  }
  std::cout << "num of lepton jets " << leptonJets << std::endl;
  if(leptonJets == 1)
  {
    h_recoLeptonNumInOneJetEvents->Fill(originalLeptonNum);
  }
  h_notInJetLepton->Fill(notInJetLeptonNum);

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

const void Dimuon::printGenParticleCollection(const reco::GenParticleCollection& genParts)
{
  const reco::Candidate* daughter1 = nullptr;
  const reco::Candidate* daughter2 = nullptr;

  const reco::Candidate* mother1 = nullptr;
  const reco::Candidate* mother2 = nullptr;

  int eventIndex = 0;

  //Format
  std::cout << std::left << std::setw(10) << "index, " << std::setw(10) << "pdfId, " 
	    << std::setw(10) << "status; " << std::setw(10) << "mother1; "
	    << std::setw(10) << "mother2; " << std::setw(10) << "daughter1; " 
	    << std::setw(10) << "daughter2, " << std::setw(10) << "px, " 
	    << std::setw(10) << "py, " << std::setw(10) << "pz, " 
	    << std::setw(10) << "E, " << std::setw(5) << "mass\n";

  //Prints out all of the particles
  for(auto &part : genParts){

    //Distinguishing Mothers
    if(part.numberOfMothers() == 2)
    {
      mother1 = part.mother(0);
      mother2 = part.mother(1);
    }
    else if(part.numberOfMothers() == 1)
    {
      mother1 = part.mother(0);
      mother2 = part.mother(0);
    }
    
    //Distinguishing Daughters (There can be more than 2 daughters, max of three)
    if(part.numberOfDaughters() == 2)
    {
      daughter1 = part.daughter(0);
      daughter2 = part.daughter(1);
    }
    else if(part.numberOfDaughters() == 1)
    {
      daughter1 = part.daughter(0);
      daughter2 = part.daughter(0);
    }

    std::cout << std::setw(10) << eventIndex << ", " << std::setw(10) << part.pdgId() << ", " << std::setw(10) << part.status() << "; ";

    if(part.numberOfMothers() != 0)
    {
      std::cout << std::setw(10) << getIndexOf(mother1, genParts) << "; " << std::setw(10) << getIndexOf(mother2, genParts) << "; ";
    }
    else
    {
      std::cout << std::setw(10) << 00 << "; " << std::setw(10) << 00 << "; ";
    }

    if(part.numberOfDaughters() != 0 && part.numberOfDaughters() != 3)
    {
      std::cout << std::setw(10) << getIndexOf(daughter1, genParts) << "; " << std::setw(10) << getIndexOf(daughter2, genParts) << "; ";
    }
    else
    {
      std::cout << std::setw(10) << 00 << "; " << std::setw(10) << 00 << "; ";
    }

    std::cout << std::setw(10) << part.px() << ", " << std::setw(10) << part.py() << ", " << std::setw(10) << part.pz() << ", " << std::setw(10) << part.energy() << ", " << std::setw(10) << part.mass() << "\n";
    
    eventIndex++;
  }
}

const int Dimuon::getIndexOf(const reco::Candidate* part, const reco::GenParticleCollection& genParts)
{
  int indexOf = 0;
  for(auto &possiblePart : genParts)
  {
    if(&possiblePart == part)
      return indexOf;
    indexOf++;
  }
  return -1;
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
