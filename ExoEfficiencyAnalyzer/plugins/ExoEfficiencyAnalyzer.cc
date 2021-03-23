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



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ExoEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ExoEfficiencyAnalyzer(const edm::ParameterSet&);
      ~ExoEfficiencyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
//edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file


edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
edm::EDGetTokenT<GenEventInfoProduct>           genInfoToken_;
edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupToken_;
edm::EDGetToken beamHaloSummaryToken_;
edm::EDGetToken photonsMiniAODToken_;
//double fMinPt;
edm::Service<TFileService> fs;
// output file name
 //TString outputFile_;
 // number of events in sample
 //uint32_t nEventsSample_;

TTree *fTree;
ExoDiPhotons::eventInfo_t fEventInfo;
ExoDiPhotons::genParticleInfo_t fGenPhoton1; // leading
ExoDiPhotons::genParticleInfo_t fGenPhoton2; // subleading
ExoDiPhotons::photonInfo_t fPhoton1Info;
ExoDiPhotons::photonInfo_t fPhoton2Info;

// struct eventInfo_t {
//   Long64_t run;
//   Long64_t LS;
//   Long64_t evnum;
// };
//eventInfo_t fEventInfo;
// struct photonInfo_t{
//   double pt;
//   double eta;
//   double phi;
//
// };
// photonInfo_t fPhoton1Info;
// photonInfo_t fPhoton2Info;



// edm::Service<TFileService> fs;
//edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
// edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
// edm::EDGetTokenT<GenEventInfoProduct>           genInfoToken_;
//
// edm::InputTag genParticles_;
// TTree *fTree;
//
//
//  ExoDiPhotons::eventInfo_t         fEventInfo;
//  ExoDiPhotons::genParticleInfo_t   fGenPhotonInfo; // Each entry is an in an individual photons not whole event
//  int fGenPhotonNumber;


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
ExoEfficiencyAnalyzer::ExoEfficiencyAnalyzer(const edm::ParameterSet& iConfig)
 //:
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
    :genParticlesToken_ (consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparticles"))),
    genInfoToken_ (consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"))),
    pileupToken_ (consumes<std::vector<PileupSummaryInfo> >( edm::InputTag("slimmedAddPileupInfo"))),
    beamHaloSummaryToken_ (consumes<reco::BeamHaloSummary>( edm::InputTag("BeamHaloSummary")))

{
  // genParticlesToken_ = consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparticles"));
  // genInfoToken_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"));
  // pileupToken_ = consumes<std::vector<PileupSummaryInfo> >( edm::InputTag("slimmedAddPileupInfo") );
  // beamHaloSummaryToken_ = consumes<reco::BeamHaloSummary>( edm::InputTag("BeamHaloSummary") );
  photonsMiniAODToken_ = consumes<edm::View<pat::Photon> > (iConfig.getParameter<edm::InputTag>("photonsMiniAOD"));
  //fMinPt = iConfig.getParameter<double>("minPhotonPt");
  fTree = fs->make<TTree>("fTree", "GENDiphotonTree");
  fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventBranchDefString.c_str());
  fTree->Branch("GenPhoton1",&fGenPhoton1,ExoDiPhotons::genParticleBranchDefString.c_str());
  fTree->Branch("GenPhoton2",&fGenPhoton2,ExoDiPhotons::genParticleBranchDefString.c_str());
  fTree->Branch("Photon1",&fPhoton1Info,ExoDiPhotons::photonBranchDefString.c_str());
  fTree->Branch("Photon2",&fPhoton2Info,ExoDiPhotons::photonBranchDefString.c_str());
  // fTree->Branch("Photon1", &fPhoton1Info, "pt/D:eta:phi");
  // fTree->Branch("Photon2", &fPhoton2Info, "pt/D:eta:phi");






   //now do what ever initialization is needed
}

ExoEfficiencyAnalyzer::~ExoEfficiencyAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


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


    ExoDiPhotons::InitEventInfo(fEventInfo);
    ExoDiPhotons::InitGenParticleInfo(fGenPhoton1);
    ExoDiPhotons::InitGenParticleInfo(fGenPhoton2);
    ExoDiPhotons::InitPhotonInfo(fPhoton1Info);
    ExoDiPhotons::InitPhotonInfo(fPhoton2Info);


    // fEventInfo.run = -99999.99;
    // fEventInfo.LS = -99999.99;
    // fEventInfo.evnum = -99999.99;

    // fPhoton1Info.pt = -99999.99;
    // fPhoton1Info.eta =-99999.99;
    // fPhoton1Info.phi =-99999.99;
    // fPhoton2Info.pt =-99999.99;
    // fPhoton2Info.eta =-99999.99;
    // fPhoton2Info.phi =-99999.99;

    ExoDiPhotons::FillBasicEventInfo(fEventInfo, iEvent);
    ExoDiPhotons::FillGenEventInfo(fEventInfo, &(*genInfo));
    ExoDiPhotons::FillPileupInfo(fEventInfo, &(*puSummary));
    ExoDiPhotons::FillBeamHaloEventInfo(fEventInfo, bhs);
    //ExoDiPhotons::FillEventWeights(fEventInfo, outputFile_, nEventsSample_);

    // for(size_t i=0; i<genParticles->size(); i++) {
    // const auto gen = genParticles->ptrAt(i);
    // if(gen->pdgId()==22){
    // //genphotoncount++;
    //   double pt = gen->pt();
    //   double eta = gen->eta();
    //   double phi = gen->phi();
    //   if (pt> fGenPhoton1.pt){
    //       fGenPhoton2.pt = fGenPhoton1.pt;
    //       fGenPhoton2.eta = fGenPhoton1.eta;
    //       fGenPhoton2.phi = fGenPhoton1.phi;
    //       fGenPhoton1.pt = pt;
    //       fGenPhoton1.eta = eta;
    //       fGenPhoton1.phi = phi;
    //     }
    //   if ((pt< fGenPhoton1.pt) && (pt> fGenPhoton2.pt)){
    //       fGenPhoton2.pt = pt;
    //       fGenPhoton2.eta = eta;
    //       fGenPhoton2.phi = phi;
    //     }}}

    std::vector< edm::Ptr<const reco::GenParticle> > genPhotons;
    std::cout<<"Numer of gen Particle  ="<<genParticles->size()<<std::endl;
    for (size_t i = 0; i < genParticles->size(); ++i)
    {  // begain for loop
      edm::Ptr<const reco::GenParticle> gen = genParticles->ptrAt(i);
      if (gen->status()==1 && gen->pdgId() == 22)
      { genPhotons.push_back(gen);// status still need to be make sure == 1 or == 3?
        // std::cout<<"Numer of gen photons   ="<<genPhotons.size()<<std::endl;
        std::cout<<"genPhoton pt  :" <<gen->pt()<<" ; eta  :"<< gen->eta() <<"phi   :"<< gen->phi()<<std::endl;//"; pdgId:"<<gen->pdgId() <<
        //  " status:" <<gen->status()<<std::endl;//" ; motherStatus:"<< gen->motherStatus() <<" ; motherPdgId:"<< gen->motherPdgId()<<"; grandmotherStatus:"<<gen->grandmotherStatus()<< std::endl;
      }
    } // end of genParticle loop.
    std::cout<<"Numer of gen photons   ="<<genPhotons.size()<<std::endl;

       //std::cout<<"  pt   :"<<genPhotons[0]->pt()<<std::endl;
    std::vector<edm::Ptr<const pat::Photon>> patPhotons;
    std::cout<<"Numer of pat Photons  ="<<photons->size()<<std::endl;
    for (size_t i = 0; i < photons->size(); ++i)
    {   // begain for loop
      edm::Ptr<const pat::Photon> pho = photons->ptrAt(i);
      patPhotons.push_back(pho);
      std::cout<<"pat Photon pt:" <<pho->pt()<<" ; eta:"<< pho->eta() <<" ; phi:"<< pho->phi()<< std::endl;
    }//end for loop

// if there is only one GenPatricle that passes pdgId and ststus.
    if(genPhotons.size() == 1)
    {
      const reco::GenParticle *genPhoton1 = &(*genPhotons.at(0));
      ExoDiPhotons::FillGenParticleInfo(fGenPhoton1, genPhoton1);
      std::vector<float> deltaR0;
      for (size_t i = 0; i < patPhotons.size(); i++)
      {
      float deltar0 = reco::deltaR(genPhotons[0]->eta(), genPhotons[0]->phi(), patPhotons[i]->eta(), patPhotons[i]->phi());

      deltaR0.push_back(deltar0);

      }
      float deltaR0min = *std::min_element(deltaR0.begin(),deltaR0.end());


      for(size_t i =0; i < deltaR0.size(); i++)
      {
        if(deltaR0min == deltaR0[i])
        {
          const pat::Photon *patPhoton1 = &(*patPhotons[i]);
          ExoDiPhotons::FillBasicPhotonInfo(fPhoton1Info, patPhoton1);
        }
      }

    }

// if there are more than 1 GenParticles that passes pdgId and status
    if(genPhotons.size() > 1)
    {
      sort(genPhotons.begin(), genPhotons.end(), ExoDiPhotons::comparePhotonsByPt);
      const reco::GenParticle *genPhoton1 = &(*genPhotons.at(0));
      const reco::GenParticle *genPhoton2 = &(*genPhotons.at(1));
      ExoDiPhotons::FillGenParticleInfo(fGenPhoton1, genPhoton1);
      ExoDiPhotons::FillGenParticleInfo(fGenPhoton2, genPhoton2);
        if(patPhotons.size()==1)
        {
          sort(patPhotons.begin(), patPhotons.end(), ExoDiPhotons::comparePhotonsByPt);
          const pat::Photon *patPhoton1 = &(*patPhotons.at(0));
          double deltaR0  = reco::deltaR(genPhoton1->eta(), genPhoton1->phi(), patPhoton1->eta(), patPhoton1->phi());
          double deltaR1  = reco::deltaR(genPhoton2->eta(), genPhoton2->phi(), patPhoton1->eta(), patPhoton1->phi());
            if (deltaR0 <= deltaR1)
            { ExoDiPhotons::FillBasicPhotonInfo(fPhoton1Info, patPhoton1);
            }
            else
            { ExoDiPhotons::FillBasicPhotonInfo(fPhoton2Info, patPhoton1);
            }
        }

       if(patPhotons.size() > 1)
       {
          std::vector<float> deltaR0, deltaR1;
          for (size_t i = 0; i < patPhotons.size(); i++)
          {
          float deltar0 = reco::deltaR(genPhotons[0]->eta(), genPhotons[0]->phi(), patPhotons[i]->eta(), patPhotons[i]->phi());
          float deltar1 = reco::deltaR(genPhotons[1]->eta(), genPhotons[1]->phi(), patPhotons[i]->eta(), patPhotons[i]->phi());
          deltaR0.push_back(deltar0);
          deltaR1.push_back(deltar1);
          }
          float deltaR0min = *std::min_element(deltaR0.begin(),deltaR0.end());
          float deltaR1min = *std::min_element(deltaR1.begin(),deltaR1.end());

          for(size_t i =0; i < deltaR0.size(); i++)
          {
            if(deltaR0min == deltaR0[i])
            {
              const pat::Photon *patPhoton1 = &(*patPhotons[i]);
              ExoDiPhotons::FillBasicPhotonInfo(fPhoton1Info, patPhoton1);
              std::cout<<"Photon1  " <<patPhotons[i]->pt()<<std::endl;
            }
          }
          for(size_t i =0; i < deltaR1.size(); i++)
          {
            if(deltaR1min == deltaR1[i])
            {
              const pat::Photon *patPhoton2 = &(*patPhotons[i]);
              ExoDiPhotons::FillBasicPhotonInfo(fPhoton2Info, patPhoton2);
              std::cout<<"Photon2  " <<patPhotons[i]->pt()<<std::endl;
            }
          }
        }
      }

            //    std::cout << "pt1:" << pt1 <<std::endl;
            //    bool isEB        = eta1 < 1.4442;
            //      bool isEE        = 1.566 < eta1 && eta2 < 2.5;
            //       bool isEBEtaGap  = patPho1->isEBEtaGap();
            //       bool isEBPhiGap  = patPho1->isEBPhiGap();
            //       bool isEERingGap = patPho1->isEERingGap();
            //       bool isEEDeeGap  = patPho1->isEEDeeGap();
            //       bool isEBEEGap   = patPho1->isEBEEGap();
            //     std::cout << isEB << " " << isEBEtaGap << " " << isEBPhiGap << " " << isEBEEGap << std::endl;
        //std::cout<<"genPhoton pt:" <<pt<<"; eta:"<< eta <<"; phi:"<< phi <<std::endl;
        //FillBasicPhotonInfo(fPhoton2Info, patPho1);
      //}
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
ExoEfficiencyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
