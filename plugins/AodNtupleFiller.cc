// -*- C++ -*-
//
// Package:    AodNtupleFiller
// Class:      AodNtupleFiller
// 
/**\class AodNtupleFiller AodNtupleFiller.cc 

 Description: Fill muon timing information histograms 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: AodNtupleFiller.cc,v 1.5 2011/04/06 11:44:29 ptraczyk Exp $
//
//

#include "AodNtupleFiller.h"

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
AodNtupleFiller::AodNtupleFiller(const edm::ParameterSet& iConfig) 
  :
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
  muonToken_ = consumes<pat::MuonCollection>(MuonTags_);
  vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  genParticleToken_ = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  trackingParticleToken_ = consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
}


AodNtupleFiller::~AodNtupleFiller()
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
AodNtupleFiller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using pat::MuonCollection;

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
  
  // Generated particle collection
  Handle<GenParticleCollection> genParticles;
  if (doSim)   
    iEvent.getByToken(genParticleToken_, genParticles);

  iEvent.getByToken(muonToken_,MuCollection);
  const pat::MuonCollection muonC = *(MuCollection.product());
  if (debug) cout << " Muon collection size: " << muonC.size() << endl;
  if (!muonC.size()) return;
  MuonCollection::const_iterator imuon,iimuon;

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

  int imucount=0;
  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    bool matched=false;

    reco::TrackRef glbTrack = imuon->combinedMuon();
    reco::TrackRef trkTrack = imuon->track();
    reco::TrackRef staTrack = imuon->standAloneMuon();

    pat::MuonRef muonR(MuCollection,imucount);
    imucount++;    
    
    reco::MuonTime timerpc = imuon->rpcTime();
    reco::MuonTime timemuon = imuon->time();
        
    if (imuon->pt()<thePtCut) continue;
    if (debug) cout << endl << "   Found muon. Pt: " << imuon->pt() << endl;

    hasSim = 0;
    isSTA = staTrack.isNonnull();
    isGLB = glbTrack.isNonnull();
    isLoose = muon::isLooseMuon(*imuon);
    isTight = muon::isTightMuon(*imuon, pvertex );
    isPF = imuon->isPFMuon();

    dxy = imuon->muonBestTrack()->dxy(pvertex.position());
    dz = imuon->muonBestTrack()->dz(pvertex.position());

    pt = imuon->bestTrack()->pt();
    dPt = imuon->bestTrack()->ptError();
    eta = imuon->bestTrack()->eta();
    phi = imuon->bestTrack()->phi();
    charge = imuon->bestTrack()->charge();

    muNdof = timemuon.nDof;
    muTime = timemuon.timeAtIpInOut;
    muTimeErr = timemuon.timeAtIpInOutErr;

    rpcNdof = timerpc.nDof;
    rpcTime = timerpc.timeAtIpInOut;
    rpcTimeErr = timerpc.timeAtIpInOutErr;

    dtNdof = 0;
    dtTime = 0;
    cscNdof = 0;
    cscTime = 0;

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
            genPt=iTrack->p4().Pt();
            genEta=iTrack->p4().Eta();
            genPhi=iTrack->p4().Phi();
            genBX=iTrack->eventId().bunchCrossing();
            genCharge=iTrack->pdgId()/13;
	    break;
          }
        }
    }

    t->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
AodNtupleFiller::beginJob()
{
   hFile = new TFile( out.c_str(), open.c_str() );
   hFile->cd();

   t = new TTree("MuTree", "MuTree");
   t->Branch("hasSim", &hasSim, "hasSim/O");
   t->Branch("genCharge", &genCharge, "genCharge/I");
   t->Branch("genPt", &genPt, "genPt/D");
   t->Branch("genPhi", &genPhi, "genPhi/D");
   t->Branch("genEta", &genEta, "genEta/D");
   t->Branch("genBX", &genBX, "genBX/I");

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
   t->Branch("pt", &pt, "pt/D");
   t->Branch("phi", &phi, "phi/D");
   t->Branch("eta", &eta, "eta/D");
   t->Branch("dPt", &dPt, "dPt/D");
   t->Branch("dz", &dz, "dz/D");
   t->Branch("dxy", &dxy, "dxy/D");

   t->Branch("nhits", &nhits, "nhits[4]/I");
   t->Branch("ssize", &ssize, "ssize[4]/I");

   t->Branch("muNdof", &muNdof, "muNdof/I");
   t->Branch("muTime", &muTime, "muTime/D");
   t->Branch("muTimeErr", &muTimeErr, "muTimeErr/D");
   t->Branch("dtNdof", &dtNdof, "dtNdof/I");
   t->Branch("dtTime", &dtTime, "dtTime/D");
   t->Branch("cscNdof", &cscNdof, "cscNdof/I");
   t->Branch("cscTime", &cscTime, "cscTime/D");
   t->Branch("rpcNdof", &rpcNdof, "rpcNdof/I");
   t->Branch("rpcTime", &rpcTime, "rpcTime/D");
   t->Branch("rpcTimeErr", &rpcTimeErr, "rpcTimeErr/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
AodNtupleFiller::endJob() {

  hFile->cd();
  t->Write();
  hFile->Write();
  delete t;  
}


//define this as a plug-in
DEFINE_FWK_MODULE(AodNtupleFiller);
