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
  //Muon Token
  edm::EDGetTokenT<std::vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken;  
  edm::EDGetTokenT<edm::View<pat::Jet> > jetsToken;



  // Muon histos;
  TH1F *muptdis,*muetadis,*muphidis,*muenergydis,*muchdis;


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

  // Muon token definition
  theMuonToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  /*
  vertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex_inputtag"));  
  jetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jet_inputtag"));
  */

  //  muonTag_ = iConfig.getParameter<edm::InputTag>("muons");
  /*
  muptdis = iConfig.getParameter<double>("muptdis");
  muetadis = iConfig.getParameter<double>("muetadis");
  muphidis = iConfig.getParameter<double>("muphidis");
  muenergydis = iConfig.getParameter<double>("muenergydis");
  muchdis = iConfig.getParameter<double>("muchdis");
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
   /*
   edm::Handle<pat::TauCollection> taus;
   edm::Handle<pat::METCollection> met;
   
   edm::Handle<edm::View<pat::Electron> > electron_pat;
   */
   //Vertex Collection

   //Muon Collection
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByLabel("slimmedMuons",muons);
   /*
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByLabel("slimmedJets", jets);

   edm::Handle<reco::VertexCollection> vtx_h;
   iEvent.getByLabel("offlinePrimaryVertices", vtx_h);
   */
   /*
   
   reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
   for (reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++)
     {
       isGoodVertex(*it);
       firstGoodVertex = it;
       break;
     }
   */
   // require a good vertex 
   // if (firstGoodVertex == vtx_h->end()) return;
  


   for (pat::MuonCollection::const_iterator muon = muons->begin();  muon != muons->end(); ++muon){                   

     muptdis->Fill(muon->pt());
     muetadis->Fill(muon->eta());
     muphidis->Fill(muon->phi());
     muenergydis->Fill(muon->energy());
     muchdis->Fill(muon->charge());
                                                               
   }                                                                                      

                                                                                                  
     // pT spectra of muons                                                                                                                                                              


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
  muptdis  = fs->make<TH1F>("muptdis","muptdis",500,0,2000);
  muetadis  = fs->make<TH1F>("muetadis","muetadis",500,-4,4);
  muphidis  = fs->make<TH1F>("muphidis","muphidis",500,-3.5,3.5);
  muenergydis  = fs->make<TH1F>("muenergydis","muenergydis",500,0,2000);
  muchdis  = fs->make<TH1F>("muchdis","muchdis",5,-2,2);
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
