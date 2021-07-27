// -*- C++ -*-
//
// Package:    diphoton-analysis/ExoEfficiencyAnalyzer
// Class:      ExoEfficiencyAnalyzer
//
/**\class ExoEfficiencyAnalyzer ExoEfficiencyAnalyzer.cc diphoton-analysis/ExoEfficiencyAnalyzer/plugins/ExoEfficiencyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bhim Bam
//         Created:  Fri, 05 Mar 2021 23:29:50 GMT
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


 // from our CommonClasses
 #include "diphoton-analysis/CommonClasses/interface/EventInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/BeamSpotInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/VertexInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/TriggerInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/JetInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/PhotonID.h"
 #include "diphoton-analysis/CommonClasses/interface/PhotonInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/GenParticleInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/DiPhotonInfo.h"
 #include "diphoton-analysis/CommonClasses/interface/PileupInfo.h"

 // for TFileService, trees
 #include "CommonTools/UtilAlgos/interface/TFileService.h"
 #include "FWCore/ServiceRegistry/interface/Service.h"
 #include "TTree.h"

 // for ECAL topology
 #include "Geometry/CaloTopology/interface/CaloTopology.h"
 #include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
 #include "FWCore/Framework/interface/ESHandle.h"
 #include "FWCore/Framework/interface/EventSetup.h"
 #include "DataFormats/EcalDetId/interface/EBDetId.h"
 #include "DataFormats/EcalDetId/interface/EEDetId.h"

 // for EGM ID
 #include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

 // for photons
 #include "DataFormats/PatCandidates/interface/Photon.h"
 #include "DataFormats/PatCandidates/interface/Electron.h"

 // for jets
 #include "DataFormats/PatCandidates/interface/Jet.h"

 // for genParticles
 #include "DataFormats/HepMCCandidate/interface/GenParticle.h"

 // for genEventInfo
 #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

 // for deltaR
 #include "DataFormats/Math/interface/deltaR.h"

 // for trigger and filter decisions
 #include "DataFormats/Common/interface/TriggerResults.h"
 #include "FWCore/Common/interface/TriggerNames.h"
 #include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
////-------------------------------------------------------------------------------------------------------------------

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ExoEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
   public:
      explicit ExoEfficiencyAnalyzer(const edm::ParameterSet&);
      ~ExoEfficiencyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------


  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;

  edm::EDGetTokenT<GenEventInfoProduct>           genInfoToken_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupToken_;

  edm::EDGetToken beamHaloSummaryToken_;

  edm::EDGetToken photonsMiniAODToken_;

  edm::EDGetTokenT<double> rhoToken_;
  double rho_;

  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;


  bool isSat = false; // Currently set to false but need to checkout global tag stuffs to do it correctly


  edm::Service<TFileService> fs;

  //edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_;
  // edm::EDGetTokenT<reco::ConversionCollection> hConversionsToken_;

  edm::EDGetToken electronToken_;
  edm::EDGetToken hConversionsToken_;




  TTree *fTree;
  ExoDiPhotons::eventInfo_t fEventInfo;
  ExoDiPhotons::genParticleInfo_t fGenPhoton1; // leading
  ExoDiPhotons::genParticleInfo_t fGenPhoton2; // subleading
  ExoDiPhotons::photonInfo_t fPhoton1Info;
  ExoDiPhotons::photonInfo_t fPhoton2Info;
  ExoDiPhotons::vertexInfo_t fVertex0Info;
  ExoDiPhotons::vertexInfo_t fPrimaryVertexInfo;
  ExoDiPhotons::beamSpotInfo_t fBeamSpotInfo;

  int nPV_;







};
////----------------------------------------------------------------------------------------------------------------------------------
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ExoEfficiencyAnalyzer::ExoEfficiencyAnalyzer(const edm::ParameterSet& iConfig)
 //:
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
    :genParticlesToken_ (consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparticles"))),
    genInfoToken_ (consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"))),
    pileupToken_ (consumes<std::vector<PileupSummaryInfo> >( edm::InputTag("slimmedAddPileupInfo"))),
    beamHaloSummaryToken_ (consumes<reco::BeamHaloSummary>( edm::InputTag("BeamHaloSummary"))),
    rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
    verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot")))
    // electronToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSource"))),
    // hConversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversionSource")))
    // electronToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronSource")))
    // hConversionsToken_(consumes<reco::Conversion>(iConfig.getParameter<edm::InputTag>("conversionSource")))

{
    // genParticlesToken_ = consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparticles"));
    // genInfoToken_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"));
    // pileupToken_ = consumes<std::vector<PileupSummaryInfo> >( edm::InputTag("slimmedAddPileupInfo") );
    // beamHaloSummaryToken_ = consumes<reco::BeamHaloSummary>( edm::InputTag("BeamHaloSummary") );
    photonsMiniAODToken_ = consumes<edm::View<pat::Photon>> (iConfig.getParameter<edm::InputTag>("photonsMiniAOD"));
    electronToken_ = consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("electronSource"));
    hConversionsToken_ = consumes<edm::View<reco::Conversion>> (iConfig.getParameter<edm::InputTag>("conversionSource"));
    // recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
    // recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
    // recHitsEBToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEBTag_);
    // recHitsEEToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEETag_);
    //fMinPt = iConfig.getParameter<double>("minPhotonPt");
    fTree = fs->make<TTree>("fTree", "GENDiphotonTree");
    fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventBranchDefString.c_str());
    fTree->Branch("GenPhoton1",&fGenPhoton1,ExoDiPhotons::genParticleBranchDefString.c_str());
    fTree->Branch("GenPhoton2",&fGenPhoton2,ExoDiPhotons::genParticleBranchDefString.c_str());
    fTree->Branch("Photon1",&fPhoton1Info,ExoDiPhotons::photonBranchDefString.c_str());
    fTree->Branch("Photon2",&fPhoton2Info,ExoDiPhotons::photonBranchDefString.c_str());
    fTree->Branch("Vertex0",&fVertex0Info,ExoDiPhotons::vertexBranchDefString.c_str());
    fTree->Branch("PrimaryVertex",&fPrimaryVertexInfo,ExoDiPhotons::vertexBranchDefString.c_str());
    fTree->Branch("nPV", &nPV_);
    fTree->Branch("BeamSpot",&fBeamSpotInfo,ExoDiPhotons::beamSpotBranchDefString.c_str());








}
////--------------------------------------------------------------------------------------------------------------------------------
ExoEfficiencyAnalyzer::~ExoEfficiencyAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
////---------------------------------------------------------------------------------------------------------------------------------

//
// member functions
//

// ------------ method called for each event  ------------
void
ExoEfficiencyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 using namespace edm;

  // Handle<TrackCollection> tracks;
  // iEvent.getByToken(tracksToken_, tracks);
  // for(TrackCollection::const_iterator itTrack = tracks->begin();
  //     itTrack != tracks->end();
  //     ++itTrack) {
  //   // do something with track parameters, e.g, plot the charge.
  //   // int charge = itTrack->charge();
  // }

  nPV_ = 0;

  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  edm::Handle<GenEventInfoProduct>           genInfo;
  iEvent.getByToken(genInfoToken_,         genInfo);

  edm::Handle<std::vector< PileupSummaryInfo > > puSummary;
  iEvent.getByToken(pileupToken_, puSummary);

  edm::Handle< reco::BeamHaloSummary > bhsHandle;
  iEvent.getByToken(beamHaloSummaryToken_, bhsHandle);
  const reco::BeamHaloSummary* bhs = &(*bhsHandle);
  // get photon collection
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByToken(photonsMiniAODToken_,photons);

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(verticesToken_,vertices);



  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);


  // edm::Handle<reco::GsfElectronCollection> hElectrons;
  // iEvent.getByToken(electronToken_, hElectrons);
  //
  // edm::Handle<reco::ConversionCollection> hConversions;
  // iEvent.getByToken(hConversionsToken_, hConversions);

  edm::Handle<edm::View<pat::Electron>> hElectrons;
  iEvent.getByToken(electronToken_, hElectrons);

  edm::Handle<edm::View<reco::Conversion>> hConversions;
  iEvent.getByToken(hConversionsToken_, hConversions);







  ExoDiPhotons::InitEventInfo(fEventInfo);
  ExoDiPhotons::InitGenParticleInfo(fGenPhoton1);
  ExoDiPhotons::InitGenParticleInfo(fGenPhoton2);
  ExoDiPhotons::InitPhotonInfo(fPhoton1Info);
  ExoDiPhotons::InitPhotonInfo(fPhoton2Info);
  ExoDiPhotons::InitVertexInfo(fVertex0Info);
  ExoDiPhotons::InitVertexInfo(fPrimaryVertexInfo);
  ExoDiPhotons::InitBeamSpotInfo(fBeamSpotInfo);




  ExoDiPhotons::FillBasicEventInfo(fEventInfo, iEvent);
  ExoDiPhotons::FillGenEventInfo(fEventInfo, &(*genInfo));
  ExoDiPhotons::FillPileupInfo(fEventInfo, &(*puSummary));
  ExoDiPhotons::FillBeamHaloEventInfo(fEventInfo, bhs);
  ExoDiPhotons::FillVertexInfo(fVertex0Info,&(vertices->at(0)));
  ExoDiPhotons::FillBeamSpotInfo(fBeamSpotInfo,&(*beamSpotHandle));


  std::vector< edm::Ptr<const reco::GenParticle> > genPhotons;
  std::cout<<"Numer of gen Particles: "<<genParticles->size()<<std::endl;
// for Photons------------------------------------------------
  for (size_t i = 0; i < genParticles->size(); ++i)
  {
    edm::Ptr<const reco::GenParticle> gen = genParticles->ptrAt(i);
    if (gen->status()==1 && gen->pdgId() == 22)
    {
        genPhotons.push_back(gen);// status still need to be make sure == 1 or == 3
        std::cout<<"genPhoton  pt  :  " <<gen->pt()<<" ; eta  :  "<< gen->eta() <<" ; phi   :  "<< gen->phi()<< ";  status   :  "<<gen->status() <<std::endl;
    }
  }
  std::cout<<"Numer of gen photons:   "<<genPhotons.size()<<std::endl;

//for electrons------------------------------------------------
  // for (size_t i = 0; i < genParticles->size(); ++i)
  // {
  //   edm::Ptr<const reco::GenParticle> gen = genParticles->ptrAt(i);
  //   if (gen->status()==1 && (gen->pdgId() == 11 || gen->pdgId() == -11))
  //   {
  //     genPhotons.push_back(gen);// status still need to be make sure
  //     std::cout<<"genElectron  pt  :  " <<gen->pt()<<" ; eta  :  "<< gen->eta() <<" ; phi   :  "<< gen->phi()<< " ; pdgId  :  "<< gen->pdgId()<<";  status   :  "<<gen->status() <<std::endl;
  //   }
  // }
  // std::cout<<"Numer of gen Electrons:   "<<genPhotons.size()<<std::endl;
//------------------------------------------------------------------------


  float mindeltaR = 0.5;
  std::cout<<"Numer of pat Photons: "<<photons->size()<<std::endl;

// Match with GenPhoton1

  const pat::Photon *matchpatPhoton = NULL;
  if (genPhotons.size() >= 1)
  {
    const reco::GenParticle *genPhoton1 = &(*genPhotons.at(0));
    ExoDiPhotons::FillGenParticleInfo(fGenPhoton1, genPhoton1);
    for (size_t i = 0; i < photons->size(); i++)
    {
      std::cout<<"pat Photon  pt:  " <<photons->ptrAt(i)->pt()<<" ; eta:  "<< photons->ptrAt(i)->eta() <<" ; phi:  "<< photons->ptrAt(i)->phi()<< std::endl;
      float deltaR = reco::deltaR(genPhoton1->eta(), genPhoton1->phi(), photons->ptrAt(i)->eta(), photons->ptrAt(i)->phi());
      if (deltaR <= mindeltaR)
      {
        mindeltaR = deltaR;
        matchpatPhoton = &(*photons->ptrAt(i));
      }
    }
    if (matchpatPhoton)
    {
    std::cout<<"Photon1  pt:  "<<matchpatPhoton->pt()<<" ; eta:  "<<matchpatPhoton->eta()<<" ; phi:  "<<matchpatPhoton->phi()<<std::endl;
    ExoDiPhotons::FillBasicPhotonInfo(fPhoton1Info, matchpatPhoton);
    ExoDiPhotons::FillPhotonIDInfo(fPhoton1Info, matchpatPhoton, rho_, isSat);
    ////----------------------------------------------CSEV Check------------------------------------------------------------------------------------------------
    bool passmyElectronVeto = !ConversionTools::hasMatchedPromptElectron(matchpatPhoton->superCluster(), hElectrons, hConversions, beamSpotHandle->position());
    std::cout<< "Def Veto Result  :" << matchpatPhoton->passElectronVeto() << "          my electron Veto Result  :" << passmyElectronVeto << std::endl;

    ////--------------------------------------------------------------------------------------------------------------------------------------------------------------
    }
  }

  //Match with GenPhoton2
  mindeltaR = 0.5;
  matchpatPhoton = NULL;
  if (genPhotons.size() >= 2)
  {
    const reco::GenParticle *genPhoton2 = &(*genPhotons.at(1));
    ExoDiPhotons::FillGenParticleInfo(fGenPhoton2, genPhoton2);
    for (size_t i = 0; i < photons->size(); i++)
    {
      float deltaR = reco::deltaR(genPhoton2->eta(), genPhoton2->phi(), photons->ptrAt(i)->eta(), photons->ptrAt(i)->phi());
      if (deltaR <= mindeltaR)
      {
        mindeltaR = deltaR;
        matchpatPhoton = &(*photons->ptrAt(i));

      }
    }
    if (matchpatPhoton)
    {
    std::cout<<"Photon2:   pt  "<<matchpatPhoton->pt()<<"  eta  "<<matchpatPhoton->eta()<<"  phi  "<<matchpatPhoton->phi()<<std::endl;
    ExoDiPhotons::FillBasicPhotonInfo(fPhoton2Info, matchpatPhoton);
    ExoDiPhotons::FillPhotonIDInfo(fPhoton2Info, matchpatPhoton, rho_, isSat);
    }
  }



// for the vertex

  const reco::Vertex *myPV = &(vertices->front());
  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++)
  {
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())
    {
      if (!foundPV)
      {
	myPV = &(vertices->at(i));
	foundPV = true;
      }
      nPV_++;
    }
  }
  ExoDiPhotons::FillVertexInfo(fPrimaryVertexInfo,&(*myPV));


  fTree->Fill();







  #ifdef THIS_IS_AN_EVENT_EXAMPLE
     Handle<ExampleData> pIn;
     iEvent.getByLabel("example",pIn);
  #endif

  #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
     ESHandle<SetupData> pSetup;
     iSetup.get<SetupRecord>().get(pSetup);
  #endif
}


// ------------ method called once each job just before starting event loop  ------------
////-------------------------------------------------------------------------------------------------------------------------------------------
void
ExoEfficiencyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ExoEfficiencyAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExoEfficiencyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
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
DEFINE_FWK_MODULE(ExoEfficiencyAnalyzer);
