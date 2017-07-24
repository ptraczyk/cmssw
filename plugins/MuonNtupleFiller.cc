// -*- C++ -*-
//
// Package:    MuonNtupleFiller
// Class:      MuonNtupleFiller
// 
/**\class MuonNtupleFiller MuonNtupleFiller.cc 

 Description: Fill muon timing information histograms 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: MuonNtupleFiller.cc,v 1.5 2011/04/06 11:44:29 ptraczyk Exp $
//
//

#include "MuonNtupleFiller.h"

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "DataFormats/DTRecHit/interface/DTSLRecSegment2D.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCleaner.h"
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>

//
// constructors and destructor
//
MuonNtupleFiller::MuonNtupleFiller(const edm::ParameterSet& iConfig) 
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  TimeTags_(iConfig.getUntrackedParameter<edm::InputTag>("Timing")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theDebug(iConfig.getParameter<bool>("debug")),
  doSim(iConfig.getParameter<bool>("mctruthMatching")),
  theAngleCut(iConfig.getParameter<double>("angleCut")),
  thePtCut(iConfig.getParameter<double>("PtCut"))
{
  edm::ConsumesCollector collector(consumesCollector());
  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trackToken_ = consumes<reco::TrackCollection>(TKtrackTags_);
  muonToken_ = consumes<reco::MuonCollection>(MuonTags_);
  muons_muonShowerInformation_token_  = consumes<edm::ValueMap<reco::MuonShower>>(edm::InputTag("muons", "muonShowerInformation", "RECO"));
  muCollToken_ = consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"));
  vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  timeMapCmbToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"combined"));
  timeMapDTToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"dt"));
  timeMapCSCToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"csc"));
  genParticleToken_ = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  trackingParticleToken_ = consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
  rpcRecHitToken_ = consumes<RPCRecHitCollection>(edm::InputTag("rpcRecHits")) ;
  l1extraToken_ = consumes<vector<l1extra::L1MuonParticle>>(edm::InputTag("l1extraParticles"));
}


MuonNtupleFiller::~MuonNtupleFiller()
{
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MuonNtupleFiller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;

  bool debug=theDebug;
  bool tpart=false;

  event_run = iEvent.id().run();
  event_lumi = iEvent.id().luminosityBlock();
  event_event = iEvent.id().event();

  if (debug)
    cout << endl << " Event: " << iEvent.id() << "  Orbit: " << iEvent.orbitNumber() << "  BX: " << iEvent.bunchCrossing() << endl;

  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexToken_,recVtxs);
  unsigned int theIndexOfThePrimaryVertex = 999.;
  for (unsigned int ind=0; ind<recVtxs->size(); ++ind) {
    if ( (*recVtxs)[ind].isValid() ) {
      theIndexOfThePrimaryVertex = ind;
      break;
    }
  }
  if (theIndexOfThePrimaryVertex<100) {
    posVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).error();
  }

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
  if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;
    else {
      cout << "No beam spot available from EventSetup." << endl;
      return;
    }
  math::XYZPoint beamspot(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());

  const reco::Vertex pvertex(posVtx,errVtx);
  posVtx = beamSpot.position();
  errVtx(0,0) = beamSpot.BeamWidthX();
  errVtx(1,1) = beamSpot.BeamWidthY();
  errVtx(2,2) = beamSpot.sigmaZ();

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  iEvent.getByToken(trackingParticleToken_, TruthTrackContainer );
  if (!TruthTrackContainer.isValid()) {
    if (debug) cout << "No trackingparticle data in the Event" << endl;
  } else tpart=true;
  
  const TrackingParticleCollection *tPC=0;
  if (tpart) 
    tPC = TruthTrackContainer.product();
  
  edm::Handle<reco::TrackCollection> trackc;
  iEvent.getByToken( trackToken_, trackc);
  const reco::TrackCollection trackC = *(trackc.product());
  if (trackC.size()>2) isCollision=1;
    else isCollision=0;

  // Generated particle collection
  Handle<GenParticleCollection> genParticles;
  if (doSim)   
    iEvent.getByToken(genParticleToken_, genParticles);

  iEvent.getByToken(muonToken_,MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());
  if (debug) cout << " Muon collection size: " << muonC.size() << endl;
  if (!muonC.size()) return;
  MuonCollection::const_iterator imuon,iimuon;

  edm::Handle<edm::ValueMap<reco::MuonShower> > muonShowerInformationValueMapH_;
  iEvent.getByToken(muons_muonShowerInformation_token_, muonShowerInformationValueMapH_);

  // check for back-to-back dimuons
  float angle=0;
  isCosmic=0;
  if (muonC.size()>1) {
    for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon) {
      if ((imuon->isGlobalMuon() || imuon->isTrackerMuon()) && (imuon->track().isNonnull())) {
        for(iimuon = imuon+1; iimuon != muonC.end(); ++iimuon) 
          if ((iimuon->isGlobalMuon() || iimuon->isTrackerMuon()) && (iimuon->track().isNonnull())) {
            double cross = imuon->track()->momentum().Dot(iimuon->track()->momentum());
            angle = acos(-cross/iimuon->track()->p()/imuon->track()->p());
            if (angle<theAngleCut) isCosmic=1;
          }
      }
    }
  }

  // Analyze L1 information
  edm::Handle<l1t::MuonBxCollection> muColl;
  iEvent.getByToken(muCollToken_, muColl);

  // 2015 version
  Handle<View<reco::Candidate> > reco;
  Handle<vector<l1extra::L1MuonParticle> > l1s;

  iEvent.getByToken(l1extraToken_, l1s);

  iEvent.getByToken(timeMapCmbToken_,timeMap1);
//  const reco::MuonTimeExtraMap & timeMapCmb = *timeMap1;
  iEvent.getByToken(timeMapDTToken_,timeMap2);
  const reco::MuonTimeExtraMap & timeMapDT = *timeMap2;
  iEvent.getByToken(timeMapCSCToken_,timeMap3);
  const reco::MuonTimeExtraMap & timeMapCSC = *timeMap3;

  int imucount=0;
  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    reco::MuonRef muonR(MuCollection,imucount);
    imucount++;    

    if (imuon->pt()<thePtCut) continue;
    if (!muon::isLooseMuon(*imuon)) continue;
    if (debug) 
      cout << endl << "   Found muon. Pt: " << imuon->pt() << "   eta: " << imuon->eta() << endl;

    reco::TrackRef glbTrack = imuon->combinedMuon();
    reco::TrackRef trkTrack = imuon->track();
    reco::TrackRef staTrack = imuon->standAloneMuon();
    reco::MuonShower muonShowerInformation = (*muonShowerInformationValueMapH_)[muonR];
    
    // MuonTimeExtra timec = timeMapCmb[muonR];
    MuonTimeExtra timedt = timeMapDT[muonR];
    MuonTimeExtra timecsc = timeMapCSC[muonR];
    reco::MuonTime timerpc = imuon->rpcTime();
    reco::MuonTime timemuon = imuon->time();

    hasSim = 0;
    isSTA = staTrack.isNonnull();
    isGLB = glbTrack.isNonnull();
    isLoose = muon::isLooseMuon(*imuon);
//    isTight = muon::isTightMuon(*imuon, pvertex );
    isTight = muon::isHighPtMuon(*imuon, pvertex );
    isPF = imuon->isPFMuon();

    dxy = imuon->tunePMuonBestTrack()->dxy(pvertex.position());
    dz = imuon->tunePMuonBestTrack()->dz(pvertex.position());

    // fill muon kinematics
    pt = imuon->tunePMuonBestTrack()->pt();
    dPt = imuon->tunePMuonBestTrack()->ptError();
    eta = imuon->tunePMuonBestTrack()->eta();
    phi = imuon->tunePMuonBestTrack()->phi();
    charge = imuon->tunePMuonBestTrack()->charge();

    vector<int> rpchits={0,0,0,0};
    if (isSTA) rpchits=countRPChits(staTrack,iEvent);
    
//    double detaphi=999;
    int l1idx=0;
    for (int i=0;i<10;i++) genPt[i]=0;
    // get L1 information
    for (int ibx=muColl->getFirstBX(); ibx<=muColl->getLastBX(); ibx++)
      for (auto it = muColl->begin(ibx); it != muColl->end(ibx); it++){
        l1t::MuonRef l1muon(muColl, distance(muColl->begin(muColl->getFirstBX()),it) );
        double deta=fabs(l1muon->eta()-eta);
        double dphi=fabs(reco::deltaPhi(l1muon->phi(),phi));
        double dr=sqrt(deta*deta+dphi*dphi);
        
        // insert a protection against cross-matching L1 candidates for close-by muons
        // If the L1 candidate dR w.r.t. a different Loose muon is smaller then don't look at the candidate here
        // (it's blocked by blowing up the dr so that it's not accepted by the dr<0.15 cut)
        for(iimuon = muonC.begin(); iimuon != muonC.end(); ++iimuon){
          if (iimuon!=imuon && iimuon->pt()>thePtCut && muon::isLooseMuon(*iimuon)) {
            double deta2=fabs(l1muon->eta()-iimuon->tunePMuonBestTrack()->eta());
            double dphi2=fabs(reco::deltaPhi(l1muon->phi(),iimuon->tunePMuonBestTrack()->phi()));
            double dr2=sqrt(deta2*deta2+dphi2*dphi2);
            if (dr2<dr) {
              dr=1.;
              break;
            }
          }
        }

        // any L1 match falling inside the 0.2 cone is saved
        if (dr<0.2) {
          hasSim=1;
          genPt[l1idx]=l1muon->pt();
          genEta[l1idx]=l1muon->eta();
          genPhi[l1idx]=l1muon->phi();
          genCharge[l1idx]=l1muon->hwQual();
          genBX[l1idx]=ibx;
          if (l1idx==9) cout << " Too many L1 matches..." << endl;
            else l1idx++;
        }
    }

/*
    int closeidx=0;
    for (int i = 0, n = l1s->size(); i < n; ++i) {
      const l1extra::L1MuonParticle & l1muon = (*l1s)[i];
      double deta=fabs(l1muon.eta()-imuon->bestTrack()->eta());
      double dphi=fabs(reco::deltaPhi(l1muon.phi(),imuon->bestTrack()->phi()));
      double dr=sqrt(deta*deta+dphi*dphi);
      if (dr<detaphi) {
        detaphi=dr;
        closeidx=i;
      }
    }

    if (detaphi<0.15) {
      const l1extra::L1MuonParticle & closestL1muon = (*l1s)[closeidx];
//      cout << " closest L1 pt: " << closestL1muon.pt() << "   dr: " << detaphi << "   BX" << closestL1muon.bx() << endl;
      hasSim=1;
      genPt=closestL1muon.pt();
      genEta=closestL1muon.eta();
      genPhi=closestL1muon.phi();
      genBX=closestL1muon.bx();
      genCharge=closestL1muon.charge();
    }
*/

    muNdof = timemuon.nDof;
    muTime = timemuon.timeAtIpInOut;
    muTimeErr = timemuon.timeAtIpInOutErr;

    rpcNdof = timerpc.nDof;
    rpcTime = timerpc.timeAtIpInOut;
    rpcTimeErr = timerpc.timeAtIpInOutErr;

    dtNdof = timedt.nDof();
    dtTime = timedt.timeAtIpInOut();
    cscNdof = timecsc.nDof();
    cscTime = timecsc.timeAtIpInOut();

    // read muon shower information
    for (int i=0; i<4; i++) { // Loop on stations
      //float sizet    = (muonShowerInformation.stationShowerSizeT).at(i);     // the transverse size of the hit cluster
      nhits[i]  = (muonShowerInformation.nStationHits).at(i);        // number of all the muon RecHits per chamber crossed by a track (1D hits)
      nrpchits[i] = rpchits.at(i);                                   // number of RPC hits in a 0.15 cone around the track
      ssize[i]  = (muonShowerInformation.stationShowerDeltaR).at(i); // the radius of the cone containing all the hits around the track
    }

    bool matched=false;
    if (tpart) {
      for (TrackingParticleCollection::const_iterator iTrack = tPC->begin(); iTrack != tPC->end(); ++iTrack)
        if (fabs(iTrack->pdgId())==13 && iTrack->p4().Pt()>2) {
          if (trkTrack.isNonnull())
            if ((fabs(iTrack->p4().eta()-imuon->track()->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack->p4().phi(),imuon->track()->momentum().phi()))<0.05)) 
              matched=true;
          if (staTrack.isNonnull())
            if ((fabs(iTrack->p4().eta()-staTrack->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack->p4().phi(),staTrack->momentum().phi()))<0.05))
              matched=true;

          if (matched==true) {
            hasSim=1;
            genPt[0]=iTrack->p4().Pt();
            genEta[0]=iTrack->p4().Eta();
            genPhi[0]=iTrack->p4().Phi();
            genBX[0]=iTrack->eventId().bunchCrossing();
            genCharge[0]=iTrack->pdgId()/13;
            break;
          }
        }
    }

    t->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonNtupleFiller::beginJob()
{
   hFile = new TFile( out.c_str(), open.c_str() );
   hFile->cd();

   t = new TTree("MuTree", "MuTree");
   t->Branch("hasSim", &hasSim, "hasSim/O");
   t->Branch("genCharge", &genCharge, "genCharge[10]/I");
   t->Branch("genPt", &genPt, "genPt[10]/F");
   t->Branch("genPhi", &genPhi, "genPhi[10]/F");
   t->Branch("genEta", &genEta, "genEta[10]/F");
   t->Branch("genBX", &genBX, "genBX[10]/I");

   t->Branch("event_run", &event_run, "event_run/i");
   t->Branch("event_lumi", &event_lumi, "event_lumi/i");
   t->Branch("event_event", &event_event, "event_event/i");
   t->Branch("isCosmic", &isCosmic, "isCosmic/O");
   t->Branch("isCollision", &isCollision, "isCollision/O");

   t->Branch("isPF", &isPF, "isPF/O");
   t->Branch("isSTA", &isSTA, "isSTA/O");
   t->Branch("isGLB", &isGLB, "isGLB/O");
   t->Branch("isLoose", &isLoose, "isLoose/O");
   t->Branch("isTight", &isTight, "isTight/O");

   t->Branch("charge", &charge, "charge/I");
   t->Branch("pt", &pt, "pt/F");
   t->Branch("phi", &phi, "phi/F");
   t->Branch("eta", &eta, "eta/F");
   t->Branch("dPt", &dPt, "dPt/F");
   t->Branch("dz", &dz, "dz/F");
   t->Branch("dxy", &dxy, "dxy/F");

   t->Branch("nhits", &nhits, "nhits[4]/I");
   t->Branch("nrpchits", &nrpchits, "nrpchits[4]/I");
   t->Branch("ssize", &ssize, "ssize[4]/F");

//   t->Branch("muNdof", &muNdof, "muNdof/I");
//   t->Branch("muTime", &muTime, "muTime/F");
//   t->Branch("muTimeErr", &muTimeErr, "muTimeErr/F");
   t->Branch("dtNdof", &dtNdof, "dtNdof/I");
   t->Branch("dtTime", &dtTime, "dtTime/F");
//   t->Branch("cscNdof", &cscNdof, "cscNdof/I");
//   t->Branch("cscTime", &cscTime, "cscTime/F");
   t->Branch("rpcNdof", &rpcNdof, "rpcNdof/I");
   t->Branch("rpcTime", &rpcTime, "rpcTime/F");
   t->Branch("rpcTimeErr", &rpcTimeErr, "rpcTimeErr/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonNtupleFiller::endJob() {

  hFile->cd();
  t->Write();
  hFile->Write();
  delete t;  
}

double MuonNtupleFiller::iMass(reco::TrackRef imuon, reco::TrackRef iimuon) {
  double energy1 = sqrt(imuon->p() * imuon->p() + 0.011163691);
  double energy2 = sqrt(iimuon->p() * iimuon->p() + 0.011163691);
  math::XYZTLorentzVector ip4(imuon->px(),imuon->py(),imuon->pz(),energy1);
  math::XYZTLorentzVector iip4(iimuon->px(),iimuon->py(),iimuon->pz(),energy2);
  math::XYZTLorentzVector psum = ip4+iip4;
  double mmumu2 = psum.Dot(psum);
  return sqrt(mmumu2);
}

vector<int> MuonNtupleFiller::countRPChits(reco::TrackRef muon, const edm::Event& iEvent) {
  double RPCCut = 60.;

  int layercount[8]={0,0,0,0,0,0,0,0};
  vector<int> stations;
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitToken_, rpcRecHits);

  for(RPCRecHitCollection::const_iterator hitRPC = rpcRecHits->begin(); hitRPC != rpcRecHits->end(); hitRPC++) {
    if ( !hitRPC->isValid()) continue;
    RPCDetId myChamber((*hitRPC).geographicalId().rawId());
    LocalPoint posLocalRPC = hitRPC->localPosition();

    for(trackingRecHit_iterator hitC = muon->recHitsBegin(); hitC != muon->recHitsEnd(); ++hitC) {
      if (!(*hitC)->isValid()) continue; 
      if ( (*hitC)->geographicalId().det() != DetId::Muon ) continue; 
      if ( (*hitC)->geographicalId().subdetId() != MuonSubdetId::RPC ) continue;

      DetId id = (*hitC)->geographicalId();
      RPCDetId rpcDetIdHit(id.rawId());
      if (rpcDetIdHit!=myChamber) continue;
      LocalPoint posLocalMuon = (*hitC)->localPosition();
      if((fabs(posLocalMuon.x()-posLocalRPC.x())<RPCCut)) 
        layercount[(rpcDetIdHit.station()-1)*2+rpcDetIdHit.layer()-1]++;
    }
  }
  
  stations.push_back(max(layercount[0],layercount[1]));
  stations.push_back(max(layercount[2],layercount[3]));
  stations.push_back(layercount[4]);
  stations.push_back(layercount[6]);
  
//  for (int i=0;i<8;i++) cout << " L:" << layercount[i];
//  for (int i=0;i<4;i++) cout << " S:" << stations[i];
//  cout << endl;
  
  return stations;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonNtupleFiller);
