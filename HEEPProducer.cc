//Stores useful HEEP ID and Isolation Variables
//        [Notes on implementation]
//        */
//        //
//        // Original Author:  Muzamil Ahmad Bhat
//        //         Created:  Thu, 17 Sep 2016 06:55:01 GMT
//                   
//        //
//        //
//
//
//        // system includes files

#include <memory>
// user include files
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
//#include "FWCore/Utilities/interface/InputTag.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//#include "TreeMaker/Utils/interface/genParticlesProducer.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/GeometryVector/interface/Vector3DBase.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//#include "LeptoQuarkTreeMaker/Utils/interface/CutNrs.h"
//#include "LeptoQuarkTreeMaker/Utils/interface/VIDCutCodes.h"
#include "HEEP/VID/interface/CutNrs.h"
#include "HEEP/VID/interface/VIDCutCodes.h"
#include "TMath.h"
#include "TVector2.h"
#include "TLorentzVector.h"

class HEEPProducer : public edm::EDProducer {
   public:
      explicit HEEPProducer(const edm::ParameterSet&);
  //    ~TrackAndPointsProducer()
         ~HEEPProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

    edm::InputTag eletag_;
    edm::InputTag gentag_;
    edm::InputTag vtxtag_, RhoTag_;
    //edm::InputTag srcJet_;

 edm::EDGetTokenT<edm::View<pat::Electron>> ElecTok_;
 edm::EDGetTokenT<reco::VertexCollection> PrimVtxTok_;
 edm::EDGetTokenT<double> RhoTok_;


 const edm::EDGetTokenT<reco::GenParticleCollection> genPartInputToken_;

edm::EDGetTokenT<edm::ValueMap<float> > heep70trkIsolMapToken_;
edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleVetoIdCutFlowResultMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleLooseIdCutFlowResultMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdCutFlowResultMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdCutFlowResultMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleHEEPIdCutFlowResultMapToken_;
 edm::EDGetTokenT<edm::ValueMap<unsigned int> > vidBitmapToken_;

typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;




};


HEEPProducer::HEEPProducer(const edm::ParameterSet& iConfig):
 genPartInputToken_ (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gentag"))),
 heep70trkIsolMapToken_            (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("heep70trkIsolMap"))),
 electronVetoIdMapToken_      (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronVetoIdMap"))),
  electronLooseIdMapToken_     (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronLooseIdMap"))),
  electronMediumIdMapToken_    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronMediumIdMap"))),
  electronTightIdMapToken_     (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronTightIdMap"))),
 electronHEEPIdMapToken_      (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronHEEPIdMap"))),
  eleVetoIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdCutFlowResultMap"))),
  eleLooseIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdCutFlowResultMap"))),
  eleMediumIdCutFlowResultMapToken_ (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdCutFlowResultMap"))),
  eleTightIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdCutFlowResultMap"))),
  eleHEEPIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdCutFlowResultMap"))),
 vidBitmapToken_   (consumes<edm::ValueMap<unsigned int> >(iConfig.getParameter<edm::InputTag>("vidBitmap")))


{ 
    eletag_  = iConfig.getParameter<edm::InputTag>( "eletag" );
   
    vtxtag_  = iConfig.getParameter<edm::InputTag>( "vtxtag" );

    RhoTag_ = edm::InputTag("fixedGridRhoFastjetCentralNeutral");
    //gentag_  = iConfig.getParameter<edm::InputTag>( "gentag");

  //genPartInputToken_ = consumes<reco::GenParticleCollection>(gentag_);

  ElecTok_ = consumes<edm::View<pat::Electron>>(eletag_);
  PrimVtxTok_ = consumes<reco::VertexCollection>(vtxtag_);
  RhoTok_ = consumes<double>(RhoTag_);
    produces<std::vector<double>>( "Eta" );
    //produces<std::vector<double>>( "Phi" );
    produces<std::vector<double>>( "Et" );
    //produces<std::vector<double>>( "Eta" );
    produces<std::vector<bool>>( "ecalDriven" );
    produces<std::vector<double>>( "DeltaEtain" );
    produces<std::vector<double>>( "DeltaPhiin" );
    produces<std::vector<double>>( "HbE" );
    produces<std::vector<double>>( "SiEtaiEta" );
    produces<std::vector<double>>( "Ecaliso" );
    produces<std::vector<double>>( "HD1iso" );
    produces<std::vector<double>>( "HD2iso" );
    produces<std::vector<double>>( "trackiso" );
    produces<std::vector<double>>( "e25max");
    produces<std::vector<double>>( "e55");
    produces<std::vector<double>>( "e25bye55");
    produces<std::vector<double>>( "DeltaEtaSeed");
    produces<std::vector<double>>( "rho");
    produces<std::vector<int>>( "Charge");
    produces<std::vector<double>>( "ePt");
    produces<std::vector<double>>( "e15");
    produces<std::vector<double>>( "ecalEnergy");
    produces<std::vector<double>>( "full55SiEtaiEta");
    produces<std::vector<double>>( "sce25max");
    produces<std::vector<double>>( "sce55");
    produces<std::vector<double>>( "sce25bye55");
    produces<std::vector<double>>( "e15bye55");
    produces<std::vector<double>>( "Fullsce25bye55");
    produces<std::vector<double>>( "Fulle15bye55");

    produces<std::vector<double>>( "DeltaEtaSeedscandTrack");
    produces<std::vector<double>>( "Phi");
    produces<std::vector<double>>( "eEnergy");
    produces<std::vector<double>>( "dxy");
    produces<std::vector<int>>( "losthits");
    produces<std::vector<double>>( "ePz");
    produces<std::vector<double>>( "eTheta");
    produces<std::vector<double>>( "ePx");
    produces<std::vector<double>>( "ePy");
    produces<std::vector<double>>( "normalizedChi2");
   //produces< std::vector< TLorentzVector > >("genParticle");
    produces< std::vector< int > >("PDGID");
    produces< std::vector< int > >("gencharge");
    produces< std::vector< double > >("genPt");
    produces< std::vector< double > >("genEta");
    produces< std::vector< double > >("genPhi");
    produces< std::vector< double > >("genEnergy");
    produces<std::vector<double>>( "PtHEEP");
    produces<std::vector<double>>( "scEtaa" );
    produces<std::vector<double>>( "scEta" );
    produces<std::vector<double>>( "scEnergy" );


    /*
    produces< std::vector< int > >("refpdgid");
    produces< std::vector< double > >("refe");
    produces< std::vector< double > >("refpt");
    produces< std::vector< double > >("refeta");
    produces< std::vector< double > >("refphi");
   */
   produces <std::vector<double> > ("heep70TrkIso" );
 
    produces< std::vector< int > >("motherPDGID");
    produces< std::vector< int > >("elstatus");
    /* produces< std::vector< double > >("refm");
    produces< std::vector< double > >("refy");
    produces< std::vector< double > >(" refarea");
    produces< std::vector< double > >("refdrjt");
    produces<std::vector<bool>>( "isMatched");

*/   
    produces<std::vector<bool>>( "passShowerShape");
    produces<std::vector<bool>>( "passDeltaEta");
    produces<std::vector<bool>>( "passDeltaPhi");
    produces<std::vector<bool>>( "passEMHD1iso");
    produces<std::vector<bool>>( "passHoverE");
    produces<std::vector<bool>>( "passDXY");
    produces<std::vector<bool>>( "passMissingHits");
    produces<std::vector<bool>>( "passEcaldriven");
    produces<std::vector<bool>>( "passN1TrkIso");
//    produces <float>            ("WorZSystemPt");
produces <std::vector<double> > ("worzsystempt");

}


HEEPProducer::~HEEPProducer()
{}



void HEEPProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm; 
    using namespace reco; 
    using namespace std;
edm::Handle<edm::ValueMap<float> > heep70trkIsolMapHandle;
iEvent.getByToken(heep70trkIsolMapToken_,heep70trkIsolMapHandle);


  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_cutflow_data;
  edm::Handle<edm::ValueMap<unsigned int> > vidBitmap;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  iEvent.getByToken(electronHEEPIdMapToken_,heep_id_decisions);
  iEvent.getByToken(eleVetoIdCutFlowResultMapToken_,veto_id_cutflow_data);
  iEvent.getByToken(eleLooseIdCutFlowResultMapToken_,loose_id_cutflow_data);
  iEvent.getByToken(eleMediumIdCutFlowResultMapToken_,medium_id_cutflow_data);
  iEvent.getByToken(eleTightIdCutFlowResultMapToken_,tight_id_cutflow_data);
  iEvent.getByToken(eleHEEPIdCutFlowResultMapToken_,heep_id_cutflow_data);
  iEvent.getByToken(vidBitmapToken_,vidBitmap);


  edm::Handle<edm::View<pat::Electron> > electron;
  iEvent.getByToken(ElecTok_, electron);

  edm::Handle< double > rho_;
  iEvent.getByToken(RhoTok_, rho_); // Central rho recommended for SUSY
  double rh = *rho_;

  edm::Handle<reco::VertexCollection> vtxh;
  iEvent.getByToken(PrimVtxTok_, vtxh);    
   const reco::Vertex vtx  = vtxh->at(0);




    std::auto_ptr<std::vector<double> >Eta(new std::vector<double>());
    std::auto_ptr<std::vector<double> > Et(new std::vector<double>());
    std::auto_ptr<std::vector<bool> > ecalDriven(new std::vector<bool>());
    std::auto_ptr<std::vector<double> > DeltaEtain(new std::vector<double>());
    std::auto_ptr<std::vector<double> > DeltaPhiin(new std::vector<double>());
    std::auto_ptr<std::vector<double> > HbE(new std::vector<double>());
    std::auto_ptr<std::vector<double> > SiEtaiEta(new std::vector<double>());
    std::auto_ptr<std::vector<double> >Ecaliso(new std::vector<double>());
    std::auto_ptr<std::vector<double> > HD1iso(new std::vector<double>());
    std::auto_ptr<std::vector<double> > HD2iso(new std::vector<double>());
    std::auto_ptr<std::vector<double> > trackiso(new std::vector<double>());
    std::auto_ptr<std::vector<double> > e25max(new std::vector<double>());
    std::auto_ptr<std::vector<double> > e55(new std::vector<double>());
    std::auto_ptr<std::vector<double> > e25bye55(new std::vector<double>());
    std::auto_ptr<std::vector<double> > DeltaEtaSeed(new std::vector<double>());
    std::auto_ptr<std::vector<double> > Fullsce25bye55(new std::vector<double>());

    std::auto_ptr<std::vector<double > > rho(new std::vector<double>());
    std::auto_ptr<std::vector<int > > Charge(new std::vector<int>());
    std::auto_ptr<std::vector<double > > ePt(new std::vector<double>());
    std::auto_ptr<std::vector<double > >e15 (new std::vector<double>());
    std::auto_ptr<std::vector<double > > ecalEnergy(new std::vector<double>());
    std::auto_ptr<std::vector<double > > full55SiEtaiEta(new std::vector<double>());
    std::auto_ptr<std::vector<double > > sce25max(new std::vector<double>());
    std::auto_ptr<std::vector<double > > sce55(new std::vector<double>());
    std::auto_ptr<std::vector<double > > sce25bye55(new std::vector<double>());
    std::auto_ptr<std::vector<double > >e15bye55 (new std::vector<double>());
    std::auto_ptr<std::vector<double > >Fulle15bye55 (new std::vector<double>());


    std::auto_ptr<std::vector<double > >DeltaEtaSeedscandTrack (new std::vector<double>());
    std::auto_ptr<std::vector<double > >Phi (new std::vector<double>());
    std::auto_ptr<std::vector<double > >eEnergy (new std::vector<double>());
    std::auto_ptr<std::vector<double > >dxy (new std::vector<double>());
    std::auto_ptr<std::vector<int > >losthits (new std::vector<int>());
    std::auto_ptr<std::vector<double > >ePz (new std::vector<double>());
    std::auto_ptr<std::vector<double > >eTheta (new std::vector<double>());
    std::auto_ptr<std::vector<double > >ePx (new std::vector<double>());
    std::auto_ptr<std::vector<double > >ePy (new std::vector<double>());
    std::auto_ptr<std::vector<double > >normalizedChi2 (new std::vector<double>());
    std::auto_ptr<std::vector<double > > PtHEEP(new std::vector<double>());
    std::auto_ptr<std::vector<double> > scEtaa(new std::vector<double>());
    std::auto_ptr<std::vector<double> > scEta(new std::vector<double>());
    std::auto_ptr<std::vector<double> > scEnergy(new std::vector<double>());

     
    std::auto_ptr< std::vector< int > > PDGID( new std::vector< int > () );
    std::auto_ptr< std::vector< int > >gencharge( new std::vector< int > () );
    std::auto_ptr<std::vector<double > >genPt (new std::vector<double>());
    std::auto_ptr<std::vector<double > >genEta (new std::vector<double>());
    std::auto_ptr<std::vector<double > >genPhi (new std::vector<double>());
    std::auto_ptr<std::vector<double > >genEnergy (new std::vector<double>());
    std::auto_ptr< std::vector< int > > motherPDGID( new std::vector< int > () );
    std::auto_ptr< std::vector< int > > elstatus( new std::vector< int > () );
  std::auto_ptr<std::vector<double> > heep70TrkIso ( new std::vector<double>() );
   





std::auto_ptr<std::vector<bool> > passShowerShape ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passDeltaEta ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passDeltaPhi ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passEMHD1iso ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passHoverE ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passDXY ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passMissingHits ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passEcaldriven ( new std::vector<bool>() );
std::auto_ptr<std::vector<bool> > passN1TrkIso ( new std::vector<bool>() );


  std::auto_ptr<std::vector<double> > worzsystempt ( new std::vector<double>() );



    for(edm::View<pat::Electron>::const_iterator elect=electron->begin(); elect!=electron->end(); ++elect){

//if(elect->pt() > 0.0){//pt cut

     Eta->push_back(elect->eta());
     Et->push_back(elect->caloEnergy()*sin(elect->p4().theta()));
     ecalDriven->push_back(elect->ecalDrivenSeed());
     DeltaEtain->push_back(elect->deltaEtaSuperClusterTrackAtVtx());
     DeltaPhiin->push_back(elect->deltaPhiSuperClusterTrackAtVtx());
     HbE->push_back(elect->hadronicOverEm());
     SiEtaiEta->push_back(elect->scSigmaIEtaIEta());

     Ecaliso->push_back(elect->dr03EcalRecHitSumEt());
     HD1iso->push_back(elect->dr03HcalDepth1TowerSumEt());
     HD2iso->push_back(elect->dr03HcalDepth2TowerSumEt());
     trackiso->push_back(elect->dr03TkSumPt());
     e25max->push_back(elect->e2x5Max());
     e55->push_back(elect->e5x5());
     e25bye55->push_back((elect->e5x5() > 0) ? ((elect->e2x5Max())/(elect->e5x5())) : 0);
     DeltaEtaSeed->push_back(elect->deltaEtaSuperClusterTrackAtVtx()-elect->superCluster()->eta() + elect->superCluster()->seed()->eta());  
     Charge->push_back(elect->charge());
     ePt->push_back(elect->pt());
     e15->push_back(elect->e1x5());
     ecalEnergy->push_back(elect->ecalEnergy());
     full55SiEtaiEta->push_back(elect->full5x5_sigmaIetaIeta());
     sce25max->push_back(elect->scE2x5Max());
     sce55->push_back(elect->scE5x5());
     sce25bye55->push_back((elect->scE5x5() > 0) ? (elect->scE2x5Max()/elect->scE5x5()) : 0);
     e15bye55->push_back((elect->e5x5() > 0) ? ((elect->e1x5())/(elect->e5x5())) : 0 );
     Fullsce25bye55->push_back((elect->full5x5_e5x5() > 0) ? (elect->full5x5_e2x5Max()/elect->full5x5_e5x5()) : 0);
     Fulle15bye55->push_back((elect->full5x5_e5x5() > 0) ? ((elect->full5x5_e1x5())/(elect->full5x5_e5x5())) : 0);
     scEnergy->push_back(elect->superCluster()->energy() );
     DeltaEtaSeedscandTrack->push_back(elect->deltaEtaSeedClusterTrackAtVtx());
     Phi->push_back(elect->phi());
     eEnergy->push_back(elect->energy());
     dxy->push_back(elect->gsfTrack()->dxy(vtx.position()));
     
     constexpr reco::HitPattern::HitCategory missingHitType = reco::HitPattern::MISSING_INNER_HITS;
     losthits->push_back(elect->gsfTrack()->hitPattern().numberOfHits(missingHitType));

     ePz->push_back(elect->pz());
     eTheta->push_back(elect->theta());
     ePx->push_back(elect->px());
     ePy->push_back(elect->py());
     normalizedChi2->push_back(elect->gsfTrack()->normalizedChi2());
     PtHEEP->push_back(elect->caloEnergy()*sin(elect->p4().theta()) );
     scEtaa->push_back( elect->superCluster()->energy()/cosh(elect->superCluster()->eta()) );
     scEta->push_back( elect->superCluster()->eta() );



const edm::Ptr<pat::Electron> elPtr(electron, elect - electron->begin() );
 unsigned int heepV70Bitmap = (*vidBitmap)[elPtr];

using HEEPV70 = VIDCutCodes<cutnrs::HEEPV70>; 

const bool passShowerShape1 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::SIGMAIETAIETA,HEEPV70::E2X5OVER5X5});
passShowerShape->push_back(passShowerShape1);
const bool passDeltaEta1 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DETAINSEED});
passDeltaEta->push_back(passDeltaEta1);
const bool passDeltaPhi1 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DPHIIN});
passDeltaPhi->push_back(passDeltaPhi1);
const bool passEMHD1 = HEEPV70::pass(heepV70Bitmap,{HEEPV70::EMHADD1ISO});
passEMHD1iso->push_back(passEMHD1);
const bool passHE = HEEPV70::pass(heepV70Bitmap,{HEEPV70::HADEM});
passHoverE->push_back(passHE);
const bool passDX = HEEPV70::pass(heepV70Bitmap,{HEEPV70::DXY});
passDXY->push_back(passDX);
const bool passMHits = HEEPV70::pass(heepV70Bitmap,{HEEPV70::MISSHITS});
passMissingHits->push_back(passMHits);
const bool passEdriven = HEEPV70::pass(heepV70Bitmap,{HEEPV70::ECALDRIVEN});
passEcaldriven->push_back(passEdriven);
const bool passN1TrkIso1 = HEEPV70::pass(heepV70Bitmap,HEEPV70::TRKISO,HEEPV70::IGNORE);
passN1TrkIso->push_back(passN1TrkIso1);


     if(!heep70trkIsolMapHandle.isValid())
      std::cout<<"hey handle hasn't been found"<<std::endl;
       else
      //std::cout<<"valus is       :"<<(*heep70trkIsolMapHandle)[ elPtr ]<<std::endl;
              heep70TrkIso        -> push_back((*heep70trkIsolMapHandle)[ elPtr ]);
      
                     // std::cout<<"valus is       :"<<(*heep70trkIsolMapHandle)[ elPtr ]<<std::endl;

//}//pt cut
}

     rho->push_back(rh);  




    
 

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genPartInputToken_, genParticles);



if(genParticles.isValid()) {//gen level stuff
 //  std::vector<TLorentzVector> hardProcessLeptonLorentzVectors;

     for(reco::GenParticleCollection::const_iterator iPart = genParticles->begin();
         iPart != genParticles->end();
         ++iPart){
   
    if(iPart->isHardProcess()) {
    if( fabs(iPart->pdgId()) >= 11 && fabs(iPart->pdgId()) <= 18){// && iPart->status() == 1 ){

    const Candidate * mom = iPart->mother();
    elstatus->push_back(iPart->status());
    PDGID->push_back( iPart->pdgId() );
    gencharge->push_back( iPart->charge() );
    genPt->push_back(iPart->pt());
    genEta->push_back(iPart->eta());
    genPhi->push_back(iPart->phi());
    genEnergy->push_back(iPart->energy());
    motherPDGID->push_back(mom->pdgId());
    //if(mom->pdgId() == 42)    
    //std::cout<<"mother charge is      "<<mom->threeCharge()<<"                                   "<<mom->pdgId()<<endl;
    
      }
    }
    }

if(genPt->size()==2)

worzsystempt->push_back(genPt->at(0)+genPt->at(1));


}//gen level stuff
    





iEvent.put( worzsystempt , "worzsystempt" );
     // iEvent.put(Eta , "Eta");
      iEvent.put(Et , "Et");
      iEvent.put(ecalDriven , "ecalDriven");
      iEvent.put(DeltaEtain , "DeltaEtain");
      iEvent.put(DeltaPhiin ,"DeltaPhiin");
      iEvent.put(HbE ,"HbE");
     // iEvent.put(SiEtaiEta , "SiEtaiEta");
      iEvent.put(Ecaliso ,"Ecaliso");
      iEvent.put(HD1iso ,"HD1iso");
      //iEvent.put(HD2iso , "HD2iso");
      iEvent.put(trackiso , "trackiso");
      //iEvent.put(e25max , "e25max");
      //iEvent.put(e55 , "e55");
      iEvent.put(e25bye55 , "e25bye55");
      iEvent.put(DeltaEtaSeed , "DeltaEtaSeed");
      iEvent.put(rho , "rho");
      iEvent.put(Charge , "Charge");
     // iEvent.put(ePt , "ePt");
     // iEvent.put(e15 , "e15");
      iEvent.put(ecalEnergy, "ecalEnergy");
      iEvent.put(full55SiEtaiEta , "full55SiEtaiEta");
     // iEvent.put(sce25max , "sce25max");
     // iEvent.put(sce55 , "sce55");
      iEvent.put(sce25bye55 , "sce25bye55");
      iEvent.put(e15bye55 ,"e15bye55");
      iEvent.put(Fullsce25bye55 , "Fullsce25bye55");
      iEvent.put(Fulle15bye55 ,"Fulle15bye55");

      //iEvent.put(DeltaEtaSeedscandTrack , "DeltaEtaSeedscandTrack");
      iEvent.put(Phi , "Phi");
      iEvent.put(eEnergy ,"eEnergy");
      iEvent.put(dxy , "dxy");
      iEvent.put(losthits ,"losthits");
      iEvent.put(ePz, "ePz");
     // iEvent.put(eTheta ,"eTheta");
      iEvent.put(ePx ,"ePx");
      iEvent.put(ePy ,"ePy");
     // iEvent.put(normalizedChi2, "normalizedChi2");
      iEvent.put(PDGID , "PDGID" );
      iEvent.put(gencharge , "gencharge" );
      iEvent.put(genPt, "genPt");
      iEvent.put(genEta, "genEta");
      iEvent.put(genPhi, "genPhi");
      iEvent.put(genEnergy, "genEnergy");
      iEvent.put(motherPDGID, "motherPDGID");
      iEvent.put( elstatus, "elstatus");
      iEvent.put(PtHEEP, "PtHEEP");
     // iEvent.put(scEtaa, "scEtaa");
      iEvent.put(scEta, "scEta");
      iEvent.put(scEnergy, "scEnergy");
      iEvent.put(heep70TrkIso, "heep70TrkIso");
      /*iEvent.put(passShowerShape, "passShowerShape");                   
      iEvent.put(passDeltaEta, "passDeltaEta");      
      iEvent.put(passDeltaPhi, "passDeltaPhi");*/      
      iEvent.put(passEMHD1iso, "passEMHD1iso");      
      /*iEvent.put(passHoverE, "passHoverE");      
      iEvent.put(passDXY, "passDXY");      
      iEvent.put(passMissingHits, "passMissingHits");*/      
      iEvent.put(passEcaldriven, "passEcaldriven");      
      iEvent.put(passN1TrkIso, "passN1TrkIso");      







}

void
HEEPProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------


void
HEEPProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------


void
HEEPProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------


void
HEEPProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------


void
HEEPProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------


void
HEEPProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HEEPProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 //The following says we do not know what parameters are allowed so do no validation
 //  // Please change this to state exactly what you do use, even if it is no parameters


  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in


DEFINE_FWK_MODULE(HEEPProducer);


