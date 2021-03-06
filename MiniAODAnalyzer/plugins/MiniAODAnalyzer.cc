// -*- C++ -*-
//
// Package:    FastAnalyzer/MiniAODAnalyzer
// Class:      MiniAODAnalyzer
// 
/**\class MiniAODAnalyzer MiniAODAnalyzer.cc FastAnalyzer/MiniAODAnalyzer/plugins/MiniAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Carlos Felipe Gonzalez Hernandez
//         Created:  Wed, 13 Jul 2016 12:02:35 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/normalizedPhi.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TTree.h"
#include "TBranch.h"


#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"


using namespace std;
using namespace pat;
using namespace reco;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MiniAODAnalyzer(const edm::ParameterSet&);
      ~MiniAODAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  
  bool isGoodVertex(const reco::Vertex& vtx);	

  edm::Service<TFileService> fs;

  //Objects Token
  edm::EDGetTokenT<std::vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT<std::vector<pat::Tau> > theTauToken;
  edm::EDGetTokenT<edm::View<pat::Electron> > theElectronToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > theJetToken;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken;  

  //Electron Token
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken;

  //MET
  edm::EDGetTokenT<pat::METCollection> metToken_;

  //Trigger Token
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1max_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1min_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigger_objects;

  //Trigger Paths
  string HLTPath1;/* HLTFilter1a, HLTFilter1b;
  string HLTPath2, HLTFilter2a,HLTFilter2b,HLTFilter2c,HLTFilter2d;   
  string HLTPath3, HLTFilter3a, HLTFilter3b;
  string HLTPath4 , HLTFilter4a, HLTFilter4b;  	
		  */
  // Muon histos;
  TH1F *muptdis,*muetadis,*muphidis,*muenergydis,*muchdis;
  //Jet histos
  TH1F *jetptdis,*jetetadis,*jetphidis,*jetenergydis;
  //Tau histos
  TH1F *tauptdis,*tauetadis,*tauphidis,*tauenergydis,*tauchdis;
  //Electron histos
  TH1F *electronptdis,*electronetadis,*electronphidis,*electronenergydis,*electronchdis;
  //MET
  TH1F *Met_type1PF_pt,*Met_type1PF_px, *Met_type1PF_py,*Met_type1PF_pz,*Met_type1PF_phi,*Met_type1PF_sumEt; 
   
  TTree *myTree;
  double totalweight;
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
MiniAODAnalyzer::MiniAODAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //Vertex
  vertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));  

  //Object
  theMuonToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  theTauToken = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"));
  theJetToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  theElectronToken = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  //MET
  metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));


  //Electron Vetos
  electronVetoIdMapToken =consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"));
  electronLooseIdMapToken = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"));
  electronMediumIdMapToken = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"));
  electronTightIdMapToken = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"));
  eleHEEPIdMapToken = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"));

  //Trigger
  triggerResultsTag  = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"));
  triggerPrescales   = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("trigger_prescale"));
  trigger_objects = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"));
  triggerPrescalesL1max_ = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1max"));
  triggerPrescalesL1min_ = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1min"));

  //Trigger Paths
  //  HLTPath1 = iConfig.getParameter<string>("HLTPath1");

  /*  HLTFilter1a = iConfig.getParameter<string>("HLTFilter1a");
  HLTFilter1b = iConfig.getParameter<string>("HLTFilter1b");
  HLTPath2 = iConfig.getParameter<string>("HLTPath2");
  HLTFilter2a = iConfig.getParameter<string>("HLTFilter2a");
  HLTFilter2b = iConfig.getParameter<string>("HLTFilter2b");
  */


   //now do what ever initialization is needed

}


MiniAODAnalyzer::~MiniAODAnalyzer()
{
 

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Vertex Collection
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken, vertices);

   //Object Collections
   edm::Handle<pat::MuonCollection> muons;
   edm::Handle<pat::TauCollection> taus;
   edm::Handle<edm::View<pat::Electron> > electrons;
   edm::Handle<edm::View<pat::Jet> > jets;
   iEvent.getByToken(theJetToken, jets);
   iEvent.getByToken(theMuonToken, muons); 
   iEvent.getByToken(theTauToken, taus);  
   iEvent.getByToken(theElectronToken, electrons);
   //MET
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();

   //Electron Vetos
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
   iEvent.getByToken(electronVetoIdMapToken,veto_id_decisions);
   iEvent.getByToken(electronLooseIdMapToken,loose_id_decisions);
   iEvent.getByToken(electronMediumIdMapToken,medium_id_decisions);
   iEvent.getByToken(electronTightIdMapToken,tight_id_decisions);  
   iEvent.getByToken(eleHEEPIdMapToken, heep_id_decisions);

   //Trigger
   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::PackedTriggerPrescales> trigger_prescale;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerResultsTag, triggerBits);
   iEvent.getByToken(triggerPrescales, trigger_prescale);
   iEvent.getByToken(triggerPrescalesL1min_, triggerPrescalesL1min);
   iEvent.getByToken(triggerPrescalesL1max_, triggerPrescalesL1max);
   iEvent.getByToken(trigger_objects, triggerObjects);   

   // Vertex
   reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
   for (reco::VertexCollection::const_iterator it = vertices->begin(); it != firstGoodVertex; it++)
     {
       isGoodVertex(*it);
       firstGoodVertex = it;
       break;
     }
   // require a good vertex 
   if (firstGoodVertex == vertices->end()) return;
  
   for(edm::View<pat::Electron>::const_iterator electron = electrons->begin(); electron != electrons->end(); electron++) {
     cout << "Electron pt " << electron->pt() << endl;
     electronptdis->Fill(electron->pt());
     electronetadis->Fill(electron->eta());
     electronphidis->Fill(electron->phi());
     electronenergydis->Fill(electron->energy());
     electronchdis->Fill(electron->charge());
   }
   /*
   //Trigger Information
   bool trigfired = false;
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
     if((names.triggerName(i) == HLTPath1) && (triggerBits->accept(i) == 1)) {trigfired = true;} 
   }
   if(!(trigfired)) return;
   */

   // Object Information (Jets, Tau, Electron, Muons)
   for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
     cout << "Jet pt " << jet->pt() << endl;
     jetptdis->Fill(jet->pt());
     jetetadis->Fill(jet->eta());
     jetphidis->Fill(jet->phi());
     jetenergydis->Fill(jet->energy());
   }

   if (!jets.isValid()) {
     edm::LogWarning("MiniAODAnalyzer") << "no pat::Jets (AK4) in event";
     return;
   }

   for (pat::TauCollection::const_iterator tau = taus->begin();  tau != taus->end(); ++tau){                   
     cout << "Tau pt " << tau->pt() << endl;

     tauptdis->Fill(tau->pt());
     tauetadis->Fill(tau->eta());
     tauphidis->Fill(tau->phi());
     tauenergydis->Fill(tau->energy());
     tauchdis->Fill(tau->charge());

   }                                                                                      

   for (pat::MuonCollection::const_iterator muon = muons->begin();  muon != muons->end(); ++muon){                   

     cout << "Muon pt " << muon->pt() << endl;
     muptdis->Fill(muon->pt());
     muetadis->Fill(muon->eta());
     muphidis->Fill(muon->phi());
     muenergydis->Fill(muon->energy());
     muchdis->Fill(muon->charge());
                                                               
   }                                                                                      

   //MET
   Met_type1PF_pt->Fill(met.pt());
   Met_type1PF_px->Fill(met.px());
   Met_type1PF_py->Fill(met.py());
   Met_type1PF_pz->Fill(met.pz());
   Met_type1PF_phi->Fill(met.phi());
   Met_type1PF_sumEt->Fill(met.sumEt());  
   

   /*

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODAnalyzer::beginJob()
{


  myTree = fs->make<TTree>("Tree","Tree");
  //muon histos
  muptdis  = fs->make<TH1F>("muptdis","muptdis",500,0,2000);
  muetadis  = fs->make<TH1F>("muetadis","muetadis",20,-4,4);
  muphidis  = fs->make<TH1F>("muphidis","muphidis",20,-3.5,3.5);
  muenergydis  = fs->make<TH1F>("muenergydis","muenergydis",500,0,2000);
  muchdis  = fs->make<TH1F>("muchdis","muchdis",1,-2,2);
  //jet histos
  jetptdis  = fs->make<TH1F>("jetptdis","jetptdis",500,0,2000);
  jetetadis  = fs->make<TH1F>("jetetadis","jetetadis",20,-4,4);
  jetphidis  = fs->make<TH1F>("jetphidis","jetphidis",20,-3.5,3.5);
  jetenergydis  = fs->make<TH1F>("jetenergydis","jetenergydis",500,0,2000);
  //tau histos
  tauptdis  = fs->make<TH1F>("tauptdis","tauptdis",500,0,2000);
  tauetadis  = fs->make<TH1F>("tauetadis","tauetadis",20,-4,4);
  tauphidis  = fs->make<TH1F>("tauphidis","tauphidis",20,-3.5,3.5);
  tauenergydis  = fs->make<TH1F>("tauenergydis","tauenergydis",500,0,2000);
  tauchdis  = fs->make<TH1F>("tauchdis","tauchdis",4,-2,2);
  //electron histos
  electronptdis  = fs->make<TH1F>("electronptdis","electronptdis",500,0,2000);
  electronetadis  = fs->make<TH1F>("electronetadis","electronetadis",20,-4,4);
  electronphidis  = fs->make<TH1F>("electronphidis","electronphidis",20,-3.5,3.5);
  electronenergydis  = fs->make<TH1F>("electronenergydis","electronenergydis",500,0,2000);
  electronchdis  = fs->make<TH1F>("electronchdis","electronchdis",4,-2,2);
  //MET histos
   Met_type1PF_pt = fs->make<TH1F>("METptdis","METptdis",500,0,2000);
   Met_type1PF_px = fs->make<TH1F>("METpt_x_dis","METpt_x_dis",500,0,2000);
   Met_type1PF_py = fs->make<TH1F>("METpt_y_dis","METpt_y_dis",500,0,2000);
   Met_type1PF_pz = fs->make<TH1F>("METpt_z_dis","METpt_z_dis",500,0,2000);
   Met_type1PF_phi = fs->make<TH1F>("METphidis","METphidis",20,-3.5,3.5);
   Met_type1PF_sumEt = fs->make<TH1F>("METenergydis","METenergydis",500,0,2000);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODAnalyzer);

bool 
MiniAODAnalyzer::isGoodVertex(const reco::Vertex& vtx)
{      
	if (vtx.isFake()) return false;
	if (vtx.ndof() < 4.) return false;
	if (vtx.position().Rho() > 2.) return false;
	if (fabs(vtx.position().Z()) > 24) return false;
	return true;
}   
