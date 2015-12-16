// -*- C++ -*-
//
// Package:    MuonTimingAnalyzer
// Class:      MuonTimingAnalyzer
// 
/**\class MuonTimingAnalyzer MuonTimingAnalyzer.cc 

 Description: Fill muon timing information histograms 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: MuonTimingAnalyzer.cc,v 1.5 2011/04/06 11:44:29 ptraczyk Exp $
//
//

#include "MuonTimingAnalyzer.h"

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
MuonTimingAnalyzer::MuonTimingAnalyzer(const edm::ParameterSet& iConfig) 
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  TimeTags_(iConfig.getUntrackedParameter<edm::InputTag>("Timing")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theDebug(iConfig.getParameter<bool>("debug")),
  doSim(iConfig.getParameter<bool>("mctruthMatching")),
  theIdCut(iConfig.getParameter<string>("requireId")),
  theCollVeto(iConfig.getParameter<bool>("collisionVeto")),
  theVetoCosmics(iConfig.getParameter<bool>("vetoCosmics")),
  theOnlyCosmics(iConfig.getParameter<bool>("onlyCosmics")),
  theAngleCut(iConfig.getParameter<double>("angleCut")),
  theMinEta(iConfig.getParameter<double>("etaMin")),
  theMaxEta(iConfig.getParameter<double>("etaMax")),
  thePtCut(iConfig.getParameter<double>("PtCut")),
  theMinPtres(iConfig.getParameter<double>("PtresMin")),
  theMaxPtres(iConfig.getParameter<double>("PtresMax")),
  theScale(iConfig.getParameter<double>("PlotScale")),
  theDtCut(iConfig.getParameter<int>("DTcut")),
  theCscCut(iConfig.getParameter<int>("CSCcut")),
  theNBins(iConfig.getParameter<int>("nbins"))
{
  edm::ParameterSet matchParameters = iConfig.getParameter<edm::ParameterSet>("MatchParameters");
//  iC=new edm::ConsumesCollector();
  edm::ConsumesCollector collector(consumesCollector());
  theMatcher = new MuonSegmentMatcher(matchParameters, collector);

  cout << " id cut " << theIdCut << endl;

  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trackToken_ = consumes<reco::TrackCollection>(TKtrackTags_);
  muonToken_ = consumes<reco::MuonCollection>(MuonTags_);
  vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  timeMapCmbToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"combined"));
  timeMapDTToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"dt"));
  timeMapCSCToken_ = consumes<reco::MuonTimeExtraMap>(edm::InputTag(TimeTags_.label(),"csc"));
  genParticleToken_ = consumes<GenParticleCollection>(edm::InputTag("genParticles"));

}


MuonTimingAnalyzer::~MuonTimingAnalyzer()
{
  if (theMatcher) delete theMatcher;
//  if (iC) delete iC;
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
MuonTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;

  bool debug=theDebug;
  bool tpart=false;
//  double hiptCut=50.;

  if (debug) {
//    cout << "*** Begin Muon Timing Analyzer " << endl;
    cout << endl << " Event: " << iEvent.id() << "  Orbit: " << iEvent.orbitNumber() << "  BX: " << iEvent.bunchCrossing() << endl;
  }

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);

  if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;
    else {
      cout << "No beam spot available from EventSetup." << endl;
      return;
    }

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  iEvent.getByLabel("mix","MergedTrackTruth",TruthTrackContainer );
  if (!TruthTrackContainer.isValid()) {
    if (debug) cout << "No trackingparticle data in the Event" << endl;
  } else tpart=true;
  
  const TrackingParticleCollection *tPC=0;
  if (tpart) 
    tPC = TruthTrackContainer.product();
  
  edm::Handle<reco::TrackCollection> trackc;
  iEvent.getByToken( trackToken_, trackc);
  const reco::TrackCollection trackC = *(trackc.product());
//  if (debug) cout << " General tracks size: " << trackC.size() << endl;

  // simple "collision event" veto on number of tracker tracks greater than 2
  if (theCollVeto && trackC.size()>2) return;

  iEvent.getByToken(muonToken_,MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());
  if (debug) cout << " Muon collection size: " << muonC.size() << endl;
  if (!muonC.size()) return;

  // Generated particle collection
  Handle<GenParticleCollection> genParticles;
  if (doSim)   
    iEvent.getByToken(genParticleToken_, genParticles);

  MuonCollection::const_iterator imuon,iimuon;

  float angle=0;

  // check for back-to-back dimuons
  if (muonC.size()>1) {
    for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon) {
      if ((imuon->isGlobalMuon() || imuon->isTrackerMuon()) && (imuon->track().isNonnull())) {
        for(iimuon = imuon+1; iimuon != muonC.end(); ++iimuon) 
          if ((iimuon->isGlobalMuon() || iimuon->isTrackerMuon()) && (iimuon->track().isNonnull())) {
            double cross = imuon->track()->momentum().Dot(iimuon->track()->momentum());
//            double cross = imuon->combinedMuon()->momentum().Dot(iimuon->combinedMuon()->momentum());
//            angle = acos(-cross/iimuon->combinedMuon()->p()/imuon->combinedMuon()->p());
            angle = acos(-cross/iimuon->track()->p()/imuon->track()->p());
            // Veto events with a cosmic muon top-bottom pair based on back-to-back angle
            if (theVetoCosmics && angle<theAngleCut) return;
            // Keep only events with a cosmic muon top-bottom pair based on back-to-back angle
            if (theOnlyCosmics && angle>theAngleCut) return;
          }
      }        
    }    
  } 

  // If we want to keep only cosmics - discard events where we couldn't measure the angle
  if (theOnlyCosmics && angle==0) return;

  math::XYZPoint beamspot(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());

  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexToken_,recVtxs);
  unsigned int theIndexOfThePrimaryVertex = 999.;
//  cout << " Vtx size: " << recVtxs->size() << endl;
  for (unsigned int ind=0; ind<recVtxs->size(); ++ind) {
    if ( (*recVtxs)[ind].isValid() ) {
    // && !((*recVtxs)[ind].isFake()) ) {
      theIndexOfThePrimaryVertex = ind;
      break;
    }
  }
  if (theIndexOfThePrimaryVertex<100) {
    posVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).error();
  }

  if (debug) cout << " Pvtx: " << posVtx << " " << theIndexOfThePrimaryVertex << endl;

  const reco::Vertex pvertex(posVtx,errVtx);
  posVtx = beamSpot.position();
  errVtx(0,0) = beamSpot.BeamWidthX();
  errVtx(1,1) = beamSpot.BeamWidthY();
  errVtx(2,2) = beamSpot.sigmaZ();
//  const reco::Vertex pvertex(posVtx,errVtx);

  iEvent.getByToken(timeMapCmbToken_,timeMap1);
  const reco::MuonTimeExtraMap & timeMapCmb = *timeMap1;
  iEvent.getByToken(timeMapDTToken_,timeMap2);
  const reco::MuonTimeExtraMap & timeMapDT = *timeMap2;
  iEvent.getByToken(timeMapCSCToken_,timeMap3);
  const reco::MuonTimeExtraMap & timeMapCSC = *timeMap3;

  double timet=0,timeb=0,phit=0,ptt=0,ptb=0,sptt=0,sptb=0;

  int imucount=0;
  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    double leg=0,outeta=0;
    double genpt=0,stapt=0;
    bool matched=false;

    reco::MuonRef muonR(MuCollection,imucount);
    imucount++;    
    if (imuon->pt()<thePtCut) continue;
    if ((fabs(imuon->eta())<theMinEta) || (fabs(imuon->eta())>theMaxEta)) continue;
    if (theIdCut=="glb"   && !muon::isGoodMuon(*imuon, muon::GlobalMuonPromptTight )) continue;
    if (theIdCut=="loose" && !muon::isLooseMuon(*imuon)) continue;
    if (theIdCut=="tight" && !muon::isTightMuon(*imuon, pvertex )) continue;
//    if (!imuon->track().isNonnull()) continue;

    if (debug) cout << endl << "   Found muon. Pt: " << imuon->pt() << endl;

    if (muon::isGoodMuon(*imuon, muon::GlobalMuonPromptTight )) {
      hi_id_trklay->Fill(imuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      hi_id_trkhit->Fill(imuon->innerTrack()->hitPattern().numberOfValidPixelHits());
      hi_id_statio->Fill(imuon->numberOfMatchedStations());
      hi_id_dxy->Fill(fabs(imuon->muonBestTrack()->dxy(pvertex.position())));
      hi_id_dz->Fill(fabs(imuon->muonBestTrack()->dz(pvertex.position())));
    }

    reco::TrackRef glbTrack = imuon->combinedMuon();
    reco::TrackRef trkTrack = imuon->track();
    reco::TrackRef staTrack = imuon->standAloneMuon();

    double d0=0.;
    if (glbTrack.isNonnull()) d0 = -1.*glbTrack->dxy(beamspot);    
    if (staTrack.isNonnull()) stapt=(*staTrack).pt();

    if (tpart) {
      for (TrackingParticleCollection::const_iterator iTrack = tPC->begin(); iTrack != tPC->end(); ++iTrack)
        if (fabs(iTrack->pdgId())==13 && iTrack->p4().Pt()>1) {
          if (trkTrack.isNonnull())
            if ((fabs(iTrack->p4().eta()-imuon->track()->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack->p4().phi(),imuon->track()->momentum().phi()))<0.05)) 
              matched=true;
          if (staTrack.isNonnull())       
            if ((fabs(iTrack->p4().eta()-staTrack->momentum().eta())<0.05) &&
                (fabs(reco::deltaPhi(iTrack->p4().phi(),staTrack->momentum().phi()))<0.05))
	      matched=true;        
	      
	  if (matched==true) {      
            genpt=iTrack->p4().Pt();
            if (debug) {
              cout << " Matched muon BX: " << iTrack->eventId().bunchCrossing();
              cout << "  hits: " << iTrack->numberOfTrackerLayers();
              cout << "  pT: " << iTrack->p4().Pt() << endl;
	    }
            if (iTrack->eventId().bunchCrossing()==0) matched=false;
	    break;
          }
        }    
    }    

    if (tpart && doSim && !matched) continue;

    if (trkTrack.isNonnull()) { 
      hi_tk_pt->Fill(((*trkTrack).pt()));
      hi_tk_phi->Fill(((*trkTrack).phi()));
      hi_tk_eta->Fill(((*trkTrack).eta()));
      hi_tk_chi2->Fill(((*trkTrack).normalizedChi2()));
//      hi_tk_nhits->Fill((*trkTrack).recHitsSize());
      hi_tk_nvhits->Fill((*trkTrack).found());
      stapt=(*trkTrack).pt();
    }  

    int nrpc=0;
    double trpc=0;
    if (staTrack.isNonnull()) {
      vector<const RPCRecHit*> rpcHits = theMatcher->matchRPC(*staTrack,iEvent);
      for (vector<const RPCRecHit*>::const_iterator hitRPC = rpcHits.begin(); hitRPC != rpcHits.end(); hitRPC++) {
        nrpc++;
        trpc+=(*hitRPC)->BunchX()*25.;
        if (debug) cout << "   RPC hit: " << (*hitRPC)->BunchX()*25. << endl;
      }
      if (nrpc>0) {
        trpc/=(double)nrpc;
        hi_nrpc->Fill(nrpc);
        hi_trpc->Fill(trpc);
        hi_trpc_phi->Fill(trpc,imuon->phi());
        hi_trpc_eta->Fill(trpc,imuon->eta());
      }
      stapt=(*staTrack).pt();
      hi_sta_pt->Fill((*staTrack).pt());
      hi_sta_ptres->Fill((*staTrack).pt()-genpt);
      hi_sta_phi->Fill((*staTrack).phi());
      hi_sta_eta->Fill((*staTrack).eta());
      hi_sta_chi2->Fill((*staTrack).normalizedChi2());
//      hi_sta_nhits->Fill((*staTrack).recHitsSize());
      hi_sta_nvhits->Fill((*staTrack).found());
      math::XYZPoint outerhit = imuon->standAloneMuon()->outerPosition();
      outeta = outerhit.eta();
      leg = outerhit.y();
      if (leg>0) sptt=(*staTrack).pt();
        else sptb=(*staTrack).pt();
    }  
//    cout << " ddd " << endl;

    if (glbTrack.isNonnull()) {
      hi_glb_pt->Fill(imuon->pt());
      hi_glb_ptres->Fill(imuon->pt()-genpt);
      hi_glb_eta->Fill(imuon->eta());
      hi_glb_d0->Fill(d0);
      hi_glb_phi->Fill(imuon->phi());
      hi_glb_chi2->Fill(((*glbTrack).normalizedChi2()));
//      hi_glb_nhits->Fill((*glbTrack).recHitsSize());
      hi_glb_nvhits->Fill((*glbTrack).found());
    }

    if (debug) {
      cout << endl;
      cout << "   Outer point: " << leg << "  eta: " << outeta << endl;
      cout << "|  *Track fit*      |  *pT*  |  *p*  |  *eta*  |  *phi*  |  *chi^2/ndf*  |  *rechits*  |  *valid rechits*  |" << endl;
      if (trkTrack.isNonnull()) {
        cout << "|    Tracker Track ";
        dumpTrack(trkTrack);
      }  
      if (staTrack.isNonnull()) {
        cout << "| StandAlone Track ";
        dumpTrack(staTrack);
      }
      if (glbTrack.isNonnull()) {
        reco::TrackRef fmsTrack = imuon->tpfmsTrack();
        reco::TrackRef pmrTrack = imuon->pickyTrack();
        reco::TrackRef dytTrack = imuon->dytTrack();

        cout << "|     Global Track ";
        dumpTrack(fmsTrack);      
        cout << "|     FMS Track ";
        dumpTrack(pmrTrack);      
        cout << "|     Picky Track ";
        dumpTrack(dytTrack);      
        cout << "|     DYT Track ";
        dumpTrack(glbTrack);      
      }
      cout << endl;
    }
    
    // Analyze the short info stored directly in reco::Muon
    
    reco::MuonTime muonTime;
    if (imuon->isTimeValid()) { 
      muonTime = imuon->time();
      if (debug) cout << "    Time points: " << muonTime.nDof << "  time: " << muonTime.timeAtIpInOut << endl;
      if (muonTime.nDof) {
        hi_mutime_vtx->Fill(muonTime.timeAtIpInOut);
        hi_mutime_vtx_err->Fill(muonTime.timeAtIpInOutErr);
      }
    }
    
    // Analyze the MuonTimeExtra information
    MuonTimeExtra timec = timeMapCmb[muonR];
    MuonTimeExtra timedt = timeMapDT[muonR];
    MuonTimeExtra timecsc = timeMapCSC[muonR];

    if (timec.nDof()) hi_cmbtime_ndof->Fill(timec.nDof());
    if (timedt.nDof()) hi_dttime_ndof->Fill(timedt.nDof());
    if (timecsc.nDof()) hi_csctime_ndof->Fill(timecsc.nDof());
    bool timeok = false;

    if (debug) {
      cout << "          DT nDof: " << timedt.nDof() << endl;
      cout << "         CSC nDof: " << timecsc.nDof() << endl;
      cout << "        Comb nDof: " << timec.nDof() << endl;
    }        
    
    bool idcut = true;
    if (timedt.nDof()>6 && fabs(timedt.timeAtIpInOut())>30) idcut=false;
    if (timecsc.nDof() && fabs(timecsc.timeAtIpInOut())>15) idcut=false;
    if (nrpc && fabs(trpc)>20) idcut=false;     

    if (staTrack.isNonnull()) for (int ii=0;ii<50;ii++) {
      if (nrpc && fabs(trpc)>ii) hi_id_rpccut_sta->Fill(ii);
      if (timecsc.nDof() && fabs(timecsc.timeAtIpInOut())>ii) hi_id_csccut_sta->Fill(ii);
      for (int jj=0;jj<15;jj++)
        if (timedt.nDof()>jj && fabs(timedt.timeAtIpInOut())>ii) hi_id_dtcut_sta->Fill(ii,jj);
    }

    if (glbTrack.isNonnull()) for (int ii=0;ii<50;ii++) {
      if (nrpc && fabs(trpc)>ii) hi_id_rpccut_glb->Fill(ii);
      if (timecsc.nDof() && fabs(timecsc.timeAtIpInOut())>ii) hi_id_csccut_glb->Fill(ii);
      for (int jj=0;jj<15;jj++)
        if (timedt.nDof()>jj && fabs(timedt.timeAtIpInOut())>ii) hi_id_dtcut_glb->Fill(ii,jj);
    }

    if (idcut) {
      if (glbTrack.isNonnull()) hi_glb_pt_cut->Fill(imuon->pt());
      if (staTrack.isNonnull()) hi_sta_pt_cut->Fill((*staTrack).pt());
    }
  
    if (timedt.nDof()>theDtCut) {
      timeok=true;
      if (debug) 
        cout << "          DT Time: " << timedt.timeAtIpInOut() << " +/- " << timedt.inverseBetaErr() << endl;
      hi_dttime_ibt->Fill(timedt.inverseBeta());
      hi_dttime_ibt_pt->Fill(imuon->pt(),timedt.inverseBeta());
      hi_dttime_ibt_err->Fill(timedt.inverseBetaErr());
      hi_dttime_fib->Fill(timedt.freeInverseBeta());
      hi_dttime_fib_err->Fill(timedt.freeInverseBetaErr());
      hi_dttime_vtx->Fill(timedt.timeAtIpInOut());
      hi_dttime_vtxn->Fill(timedt.timeAtIpInOut(),timedt.nDof());
      hi_dttime_vtx_w->Fill(timedt.timeAtIpInOut());
      hi_dttime_vtx_pt->Fill(timedt.timeAtIpInOut(),stapt);
      hi_dttime_vtx_phi->Fill(timedt.timeAtIpInOut(),imuon->phi());
      if (fabs(timedt.timeAtIpInOut())>30.) hi_dttime_etaphi->Fill(imuon->eta(),imuon->phi());
      hi_dttime_vtx_eta->Fill(timedt.timeAtIpInOut(),imuon->eta());
      if (timedt.timeAtIpInOut()<30.) 
        hi_dttime_eeta_lo->Fill(outeta,imuon->eta());
        else
        hi_dttime_eeta_hi->Fill(outeta,imuon->eta());
      hi_dttime_vtx_err->Fill(timedt.timeAtIpInOutErr());
      hi_dttime_vtxr->Fill(timedt.timeAtIpOutIn());
      hi_dttime_vtxr_err->Fill(timedt.timeAtIpOutInErr());
      hi_dttime_errdiff->Fill(timedt.timeAtIpInOutErr()-timedt.timeAtIpOutInErr());

      if (leg>0) {
        timet=timedt.timeAtIpInOut();
        phit=leg;
        ptt=imuon->pt();
        hi_dttime_vtx_t->Fill(timet);
        hi_dttime_vtx_etat->Fill(timedt.timeAtIpInOut(),imuon->eta());
        hi_dttime_vtxp_t->Fill(timet,imuon->phi());
        hi_dttime_fib_t->Fill(timedt.freeInverseBeta());
        hi_dttime_fibp_t->Fill(timedt.freeInverseBeta(),imuon->phi());
        hi_dttime_errdiff_t->Fill(timedt.timeAtIpInOutErr()-timedt.timeAtIpOutInErr());
        if (timecsc.nDof()>theCscCut)
          hi_dtcsc_vtx_t->Fill(timedt.timeAtIpInOut()-timecsc.timeAtIpInOut());
      } else if (leg<0) {
        timeb=timedt.timeAtIpInOut();
        hi_dttime_vtx_b->Fill(timeb);
        hi_dttime_vtx_etab->Fill(timedt.timeAtIpInOut(),imuon->eta());
        hi_dttime_vtxp_b->Fill(timeb,-imuon->phi());
        hi_dttime_fib_b->Fill(timedt.freeInverseBeta());
        hi_dttime_fibp_b->Fill(timedt.freeInverseBeta(),imuon->phi());
        hi_dttime_errdiff_b->Fill(timedt.timeAtIpInOutErr()-timedt.timeAtIpOutInErr());
        if (timecsc.nDof()>theCscCut)
          hi_dtcsc_vtx_b->Fill(timedt.timeAtIpInOut()-timecsc.timeAtIpInOut());
      }

      if (timedt.inverseBetaErr()>0.)
        hi_dttime_ibt_pull->Fill((timedt.inverseBeta()-1.)/timedt.inverseBetaErr());
      if (timedt.freeInverseBetaErr()>0.)    
        hi_dttime_fib_pull->Fill((timedt.freeInverseBeta()-1.)/timedt.freeInverseBetaErr());
      if (timedt.timeAtIpInOutErr()>0.)
        hi_dttime_vtx_pull->Fill(timedt.timeAtIpInOut()/timedt.timeAtIpInOutErr());
      if (timedt.timeAtIpOutInErr()>0.)
        hi_dttime_vtxr_pull->Fill(timedt.timeAtIpOutIn()/timedt.timeAtIpOutInErr());

      if (timecsc.nDof()>theCscCut)
        hi_dtcsc_vtx->Fill(timedt.timeAtIpInOut()-timecsc.timeAtIpInOut());

    }

    if (timecsc.nDof()>theCscCut) {
      timeok=true;
      if (debug) 
        cout << "         CSC Time: " << timecsc.timeAtIpInOut() << " +/- " << timecsc.inverseBetaErr() << endl;
      hi_csctime_ibt->Fill(timecsc.inverseBeta());
      hi_csctime_ibt_pt->Fill(imuon->pt(),timecsc.inverseBeta());
      hi_csctime_ibt_err->Fill(timecsc.inverseBetaErr());
      hi_csctime_fib->Fill(timecsc.freeInverseBeta());
      hi_csctime_fib_err->Fill(timecsc.freeInverseBetaErr());
      hi_csctime_vtx->Fill(timecsc.timeAtIpInOut());
      hi_csctime_vtxn->Fill(timecsc.timeAtIpInOut(),timecsc.nDof());
      hi_csctime_vtx_err->Fill(timecsc.timeAtIpInOutErr());
      hi_csctime_vtx_eta->Fill(timecsc.timeAtIpInOut(),imuon->eta());
      hi_csctime_vtx_phi->Fill(timecsc.timeAtIpInOut(),imuon->phi());
      hi_csctime_vtx_pt->Fill(timecsc.timeAtIpInOut(),stapt);
      if (timecsc.timeAtIpInOut()>-40.) 
        hi_csctime_eeta_lo->Fill(outeta,imuon->eta());
        else
        hi_csctime_eeta_hi->Fill(outeta,imuon->eta());
      hi_csctime_vtxr->Fill(timecsc.timeAtIpOutIn());
      hi_csctime_vtxr_err->Fill(timecsc.timeAtIpOutInErr());

      if (leg>0) {
        hi_csctime_vtx_t->Fill(timecsc.timeAtIpInOut());
        hi_csctime_vtx_etat->Fill(timecsc.timeAtIpInOut(),imuon->eta());
        hi_csctime_fib_t->Fill(timecsc.freeInverseBeta());
      }
      if (leg<0) {
        hi_csctime_vtx_b->Fill(timecsc.timeAtIpInOut());
        hi_csctime_vtx_etab->Fill(timecsc.timeAtIpInOut(),imuon->eta());
        hi_csctime_fib_b->Fill(timecsc.freeInverseBeta());
      }
      
      if (timec.inverseBetaErr()>0.)
        hi_csctime_ibt_pull->Fill((timecsc.inverseBeta()-1.)/timecsc.inverseBetaErr());
      if (timecsc.freeInverseBetaErr()>0.)    
        hi_csctime_fib_pull->Fill((timecsc.freeInverseBeta()-1.)/timecsc.freeInverseBetaErr());
      if (timecsc.timeAtIpInOutErr()>0.)
        hi_csctime_vtx_pull->Fill(timecsc.timeAtIpInOut()/timecsc.timeAtIpInOutErr());
      if (timecsc.timeAtIpOutInErr()>0.)
        hi_csctime_vtxr_pull->Fill(timecsc.timeAtIpOutIn()/timecsc.timeAtIpOutInErr());
    }
    
    if (timec.nDof()>0) {
      if (debug) 
        cout << "        Comb Time: " << timec.timeAtIpInOut() << " +/- " << timec.inverseBetaErr() << endl;
      hi_cmbtime_ibt->Fill(timec.inverseBeta());
      hi_cmbtime_ibt_pt->Fill(imuon->pt(),timec.inverseBeta());
      hi_cmbtime_ibt_err->Fill(timec.inverseBetaErr());
      hi_cmbtime_fib->Fill(timec.freeInverseBeta());
      hi_cmbtime_fib_err->Fill(timec.freeInverseBetaErr());
      hi_cmbtime_vtx->Fill(timec.timeAtIpInOut());
      hi_cmbtime_vtx_err->Fill(timec.timeAtIpInOutErr());
      hi_cmbtime_vtxr->Fill(timec.timeAtIpOutIn());
      hi_cmbtime_vtxr_err->Fill(timec.timeAtIpOutInErr());

      if (timec.inverseBetaErr()>0.)
        hi_cmbtime_ibt_pull->Fill((timec.inverseBeta()-1.)/timec.inverseBetaErr());
      if (timec.freeInverseBetaErr()>0.)    
        hi_cmbtime_fib_pull->Fill((timec.freeInverseBeta()-1.)/timec.freeInverseBetaErr());
      if (timec.timeAtIpInOutErr()>0.)
        hi_cmbtime_vtx_pull->Fill(timec.timeAtIpInOut()/timec.timeAtIpInOutErr());
      if (timec.timeAtIpOutInErr()>0.)
        hi_cmbtime_vtxr_pull->Fill(timec.timeAtIpOutIn()/timec.timeAtIpOutInErr());
    }
    
    if (timeok) {
      if (staTrack.isNonnull()) hi_sta_ptt->Fill((*staTrack).pt());
      if (glbTrack.isNonnull()) hi_glb_ptt->Fill(imuon->pt());
    }
  }  

  double minp = (ptt<ptb)?ptt:ptb;
  
  if (timet!=0 && timeb!=0) {
    hi_dttime_vtx_tb->Fill(timeb-timet);
    hi_dttime_vtxp_tb->Fill(timeb-timet,phit);
    hi_dttime_vtxpt_tb->Fill(timeb-timet,minp);
    hi_dttime_vtx_tb2->Fill(timet,timeb);
  }

  if (ptt!=0 && ptb!=0) 
    hi_glb_ptres_tb->Fill(ptt-ptb);
  if (sptt!=0 && sptb!=0) 
    hi_sta_ptres_tb->Fill(sptt-sptb);
  
  if (timet!=0 && timeb==0) hi_dttime_vtx_to->Fill(timet);
  if (timeb!=0 && timet==0) hi_dttime_vtx_bo->Fill(timeb);
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonTimingAnalyzer::beginJob()
{
   hFile = new TFile( out.c_str(), open.c_str() );
   hFile->cd();

   effStyle = new TStyle("effStyle","Efficiency Study Style");   
   effStyle->SetCanvasBorderMode(0);
   effStyle->SetPadBorderMode(1);
   effStyle->SetOptTitle(0);
   effStyle->SetStatFont(42);
   effStyle->SetTitleFont(22);
   effStyle->SetCanvasColor(10);
   effStyle->SetPadColor(0);
   effStyle->SetLabelFont(42,"x");
   effStyle->SetLabelFont(42,"y");
   effStyle->SetHistFillStyle(1001);
   effStyle->SetHistFillColor(0);
   effStyle->SetOptStat(0);
   effStyle->SetOptFit(0111);
   effStyle->SetStatH(0.05);
   
   hi_id_rpccut_sta = new TH1F("hi_id_rpccut_sta","STA muons rejected by RPC cut",50,0,50);
   hi_id_rpccut_glb = new TH1F("hi_id_rpccut_glb","GLB muons rejected by RPC cut",50,0,50);
   hi_id_csccut_sta = new TH1F("hi_id_csccut_sta","STA muons rejected by CSC cut",50,0,50);
   hi_id_csccut_glb = new TH1F("hi_id_csccut_glb","GLB muons rejected by CSC cut",50,0,50);
   hi_id_dtcut_sta = new TH2F("hi_id_dtcut_sta","STA muons rejected by DT cut",50,0,50,15,0,15);
   hi_id_dtcut_glb = new TH2F("hi_id_dtcut_glb","GLB muons rejected by DT cut",50,0,50,15,0,15);

   hi_id_trklay = new TH1F("hi_id_trklay","Tracker Layers (>5)",18,0.,18);
   hi_id_trkhit = new TH1F("hi_id_trkhit","Pixel hits (>0)",10,0.,10);
   hi_id_statio = new TH1F("hi_id_statio","Matched Stations (>1)",7,0.,7);
   hi_id_dxy = new TH1F("hi_id_dxy","Dxy (<0.2)",theNBins,0.,1);
   hi_id_dz = new TH1F("hi_id_dz","Dz (<0.5)",theNBins,0.,10);

   hi_glb_angle = new TH1F("hi_glb_angle","Dimon global-global opening angle",theNBins,0.,0.1);
   hi_trk_angle = new TH1F("hi_trk_angle","Dimon trk-trk opening angle",theNBins,0.,0.1);
   hi_glb_angle_w = new TH1F("hi_glb_angle_w","Dimon global-global opening angle",theNBins,0.,3.1);
   hi_trk_angle_w = new TH1F("hi_trk_angle_w","Dimon trk-trk opening angle",theNBins,0.,3.1);
   hi_dttime_vtx_tb_angle = new TH2F("hi_dttime_vtx_tb_angle","DT Time at Vertex (BOT-TOP) vs opening angle",60,-100.,80.,theNBins,0.,3.1);

   hi_glb_mass_os = new TH1F("hi_glb_mass_os","Opposite Sign dimuon mass (GLB)",theNBins,50.,130.);
   hi_glb_mass_ss = new TH1F("hi_glb_mass_ss","Same Sign dimuon mass (GLB)",theNBins,0.,200.);
   hi_sta_mass_os = new TH1F("hi_sta_mass_os","Opposite Sign dimuon mass (STA)",theNBins,20.,160.);
   hi_sta_mass_ss = new TH1F("hi_sta_mass_ss","Same Sign dimuon mass (STA)",theNBins,20.,200.);

   hi_sta_pt  = new TH1F("hi_sta_pt","P_{T}^{STA}",theNBins,theMinPtres,theMaxPtres);
   hi_sta_pt_cut  = new TH1F("hi_sta_pt_cut","P_{T}^{STA} after timing cut",theNBins,theMinPtres,theMaxPtres);
   hi_sta_ptres = new TH1F("hi_sta_ptres","P_{T}^{STA} - P_{T}^{gen}",theNBins,-theMaxPtres/10.,theMaxPtres/10.);
   hi_sta_ptg = new TH1F("hi_sta_ptg","P_{T}^{STA} gen matched",theNBins,theMinPtres,theMaxPtres);
   hi_sta_ptt = new TH1F("hi_sta_ptt","P_{T}^{STA} with timing",theNBins,theMinPtres,theMaxPtres);
   hi_sta_ptres_tb= new TH1F("hi_sta_ptres_tb","P_{T}^{TOP} - P_{T}^{BOT}",theNBins+1,-theMaxPtres/10.,theMaxPtres/10.);
   hi_tk_pt   = new TH1F("hi_tk_pt","P_{T}^{TK}",theNBins,theMinPtres,theMaxPtres);
   hi_glb_pt  = new TH1F("hi_glb_pt","Reco muon P_{T}",theNBins,theMinPtres,theMaxPtres);
   hi_glb_pt_cut = new TH1F("hi_glb_pt_cut","Reco muon P_{T}^{GLB} after timing cut",theNBins,theMinPtres,theMaxPtres);
   hi_glb_ptg = new TH1F("hi_glb_ptg","Reco muon P_{T}^{GLB} gen matched",theNBins,theMinPtres,theMaxPtres);
   hi_glb_ptt = new TH1F("hi_glb_ptt","Reco muon P_{T}^{GLB} with timing",theNBins,theMinPtres,theMaxPtres);
   hi_glb_ptres= new TH1F("hi_glb_ptres","P_{T}^{rec} - -P_{T}^{gen}",theNBins+1,-theMaxPtres/10.,theMaxPtres/10.);
   hi_glb_ptresh= new TH1F("hi_glb_ptresh","P_{T}^{rec} - P_{T}^{TK} for P_{T}^{TK}>45 GeV",theNBins+1,-theMaxPtres/40.,theMaxPtres/40.);
   hi_glb_ptres_b= new TH1F("hi_glb_ptres_t","P_{T}^{rec} - P_{T}^{TK} (BOT)",theNBins+1,-theMaxPtres/40.,theMaxPtres/40.);
   hi_glb_ptresh_b= new TH1F("hi_glb_ptresh_t","P_{T}^{rec} - P_{T}^{TK} for P_{T}^{TK}>45 GeV (BOT)",theNBins+1,-theMaxPtres/40.,theMaxPtres/40.);
   hi_glb_ptres_t= new TH1F("hi_glb_ptres_b","P_{T}^{rec} - P_{T}^{TK} (TOP)",theNBins+1,-theMaxPtres/40.,theMaxPtres/40.);
   hi_glb_ptresh_t= new TH1F("hi_glb_ptresh_b","P_{T}^{rec} - P_{T}^{TK} for P_{T}^{TK}>45 GeV (TOP)",theNBins+1,-theMaxPtres/40.,theMaxPtres/40.);
   hi_glb_ptres_tb= new TH1F("hi_glb_ptres_tb","P_{T}^{TOP} - P_{T}^{BOT}",theNBins+1,-theMaxPtres/10.,theMaxPtres/10.);
   hi_glb_d0   = new TH1F("hi_glb_d0","GLB D0",80,-50,50);

   hi_sta_phi = new TH1F("hi_sta_phi","#phi^{STA}",theNBins,-3.0,3.);
   hi_tk_phi  = new TH1F("hi_tk_phi","#phi^{TK}",theNBins,-3.0,3.);
   hi_glb_phi = new TH1F("hi_glb_phi","#phi^{GLB}",theNBins,-3.0,3.);
   hi_sta_eta = new TH1F("hi_sta_eta","#eta^{STA}",theNBins/2,2.5,2.5);
   hi_tk_eta  = new TH1F("hi_tk_eta","#eta^{TK}",theNBins/2,2.5,2.5);
   hi_glb_eta = new TH1F("hi_glb_eta","#eta^{GLB}",theNBins/2,2.5,2.5);
   hi_sta_nhits = new TH1F("hi_sta_nhits","StandAlone number of segments/hits",56,0.,56.0);
   hi_tk_nhits = new TH1F("hi_tk_nhits","Tracker number of hits",30,0.,30.0);
   hi_glb_nhits = new TH1F("hi_glb_nhits","Global number of segments/hits",80,0.,80.0);
   hi_sta_nvhits = new TH1F("hi_sta_nvhits","StandAlone number of valid hits",56,0.,56.0);
   hi_tk_nvhits = new TH1F("hi_tk_nvhits","Tracker number of valid hits",30,0.,30.0);
   hi_glb_nvhits = new TH1F("hi_glb_nvhits","Global number of valid hits",80,0.,80.0);
   hi_sta_chi2 = new TH1F("hi_sta_chi2","StandAlone muon normalized chi2",60,0.,6.0);
   hi_tk_chi2 = new TH1F("hi_tk_chi2","Tracker track normalized chi2",60,0.,6.0);
   hi_glb_chi2 = new TH1F("hi_glb_chi2","Global muon normalized chi2",60,0.,6.0);

   hi_mutime_vtx = new TH1F("hi_mutime_vtx","Time at Vertex (inout)",theNBins,-25.*theScale,25.*theScale);
   hi_mutime_vtx_err = new TH1F("hi_mutime_vtx_err","Time at Vertex Error (inout)",theNBins,0.,25.0);

   hi_dtcsc_vtx = new TH1F("hi_dtcsc_vtx","Time at Vertex (DT-CSC)",theNBins,-25.*theScale,25.*theScale);
   hi_dtcsc_vtx_t = new TH1F("hi_dtcsc_vtx_t","Time at Vertex (TOP DT-CSC)",theNBins,-25.*theScale,25.*theScale);
   hi_dtcsc_vtx_b = new TH1F("hi_dtcsc_vtx_b","Time at Vertex (BOT DT-CSC)",theNBins,-25.*theScale,25.*theScale);

   hi_trpc = new TH1F("hi_trpc","Time at Vertex (RPC)",theNBins/3.,-60.,60.);
   hi_nrpc = new TH1F("hi_nrpc","RPC nHits",8,0,8);
   hi_trpc_phi = new TH2F("hi_trpc_phi","RPC Time vs Phi",theNBins/3.,-60,60,60,-3.14,3.14);
   hi_trpc_eta = new TH2F("hi_trpc_eta","RPC Time vs Eta",theNBins/3.,-60,60,60,-2.1,2.1);

   hi_cmbtime_ibt = new TH1F("hi_cmbtime_ibt","Inverse Beta",theNBins,0.,1.6);
   hi_cmbtime_ibt_pt = new TH2F("hi_cmbtime_ibt_pt","P{T} vs Inverse Beta",theNBins,theMinPtres,theMaxPtres,theNBins,0.7,2.0);
   hi_cmbtime_ibt_err = new TH1F("hi_cmbtime_ibt_err","Inverse Beta Error",theNBins,0.,1.0);
   hi_cmbtime_fib = new TH1F("hi_cmbtime_fib","Free Inverse Beta",theNBins,-5.,5.);
   hi_cmbtime_fib_err = new TH1F("hi_cmbtime_fib_err","Free Inverse Beta Error",theNBins,0,5.);
   hi_cmbtime_vtx = new TH1F("hi_cmbtime_vtx","Time at Vertex (inout)",theNBins,-25.*theScale,25.*theScale);
   hi_cmbtime_vtx_err = new TH1F("hi_cmbtime_vtx_err","Time at Vertex Error (inout)",theNBins,0.,25.0);
   hi_cmbtime_vtxr = new TH1F("hi_cmbtime_vtxR","Time at Vertex (inout)",theNBins,0.,75.*theScale);
   hi_cmbtime_vtxr_err = new TH1F("hi_cmbtime_vtxR_err","Time at Vertex Error (inout)",theNBins,0.,25.0);
   hi_cmbtime_ibt_pull = new TH1F("hi_cmbtime_ibt_pull","Inverse Beta Pull",theNBins,-5.,5.0);
   hi_cmbtime_fib_pull = new TH1F("hi_cmbtime_fib_pull","Free Inverse Beta Pull",theNBins,-5.,5.0);
   hi_cmbtime_vtx_pull = new TH1F("hi_cmbtime_vtx_pull","Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_cmbtime_vtxr_pull = new TH1F("hi_cmbtime_vtxR_pull","Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_cmbtime_ndof = new TH1F("hi_cmbtime_ndof","Number of timing measurements",48,0.,48.0);

   hi_dttime_ibt = new TH1F("hi_dttime_ibt","DT Inverse Beta",theNBins,0.,1.6);
   hi_dttime_ibt_pt = new TH2F("hi_dttime_ibt_pt","P{T} vs DT Inverse Beta",theNBins,theMinPtres,theMaxPtres,theNBins,0.7,2.0);
   hi_dttime_ibt_err = new TH1F("hi_dttime_ibt_err","DT Inverse Beta Error",theNBins,0.,0.3);
   hi_dttime_fib = new TH1F("hi_dttime_fib","DT Free Inverse Beta",theNBins+1,-5.,7.);
   hi_dttime_fib_t = new TH1F("hi_dttime_fib_t","DT Free Inverse Beta (TOP)",theNBins,-5.,5.);
   hi_dttime_fib_b = new TH1F("hi_dttime_fib_b","DT Free Inverse Beta (BOT)",theNBins,-5.,5.);
   hi_dttime_fibp_t = new TH2F("hi_dttime_fibp_t","DT Free Inverse Beta (TOP)",theNBins,-5.,5.,theNBins,0.,3.14);
   hi_dttime_fibp_b = new TH2F("hi_dttime_fibp_b","DT Free Inverse Beta (BOT)",theNBins,-5.,5.,theNBins,0.,3.14);
   hi_dttime_fib_err = new TH1F("hi_dttime_fib_err","DT Free Inverse Beta Error",theNBins,0,5.);
   hi_dttime_vtx = new TH1F("hi_dttime_vtx","DT Time at Vertex",theNBins*2,-60,60);
   hi_dttime_vtxn = new TH2F("hi_dttime_vtxn","DT Time at Vertex",theNBins,-100,100,48,0.,48.0);
   hi_dttime_vtx_w = new TH1F("hi_dttime_vtx_w","DT Time at Vertex (wide)",theNBins*3,-75.*theScale,75.*theScale);
   hi_dttime_vtx_pt = new TH2F("hi_dttime_vtx_pt","Time at Vertex vs STA p_{T}",theNBins,-100,100,theNBins,theMinPtres,theMaxPtres);
   hi_dttime_vtx_phi = new TH2F("hi_dttime_vtx_phi","DT Time at Vertex vs Phi",theNBins,-100,100,60,-3.14,3.14);
   hi_dttime_vtx_eta = new TH2F("hi_dttime_vtx_eta","DT Time at Vertex vs Eta",theNBins,-100,100,60,-2.1,2.1);
   hi_dttime_etaphi = new TH2F("hi_dttime_etaphi","Eta vs Phi of muons with |DT t_{0}|>30ns",60,-2.1,2.1,60,-3.14,3.14);
   hi_dttime_eeta_lo = new TH2F("hi_dttime_eeta_lo","Pt Eta vs Origin Eta for DT in-time",60,-2.1,2.1,60,-2.1,2.1);
   hi_dttime_eeta_hi = new TH2F("hi_dttime_eeta_hi","Pt Eta vs Origin Eta for DT ou-time",60,-2.1,2.1,60,-2.1,2.1);
   hi_dttime_vtx_etat = new TH2F("hi_dttime_vtx_etat","DT Time at Vertex vs Eta (TOP)",theNBins,-100,100,60,-2.1,2.1);
   hi_dttime_vtx_etab = new TH2F("hi_dttime_vtx_etab","DT Time at Vertex vs Eta (BOT)",theNBins,-100,100,60,-2.1,2.1);
   hi_dttime_vtx_t = new TH1F("hi_dttime_vtx_t","DT Time at Vertex (TOP)",theNBins,-100.,140.);
   hi_dttime_vtx_b = new TH1F("hi_dttime_vtx_b","DT Time at Vertex (BOT)",theNBins,-100.,140.);
   hi_dttime_vtx_to = new TH1F("hi_dttime_vtx_to","DT Time at Vertex (TOP only)",theNBins,-100.,140.);
   hi_dttime_vtx_bo = new TH1F("hi_dttime_vtx_bo","DT Time at Vertex (BOT only)",theNBins,-100.,140.);
   hi_dttime_vtx_tb = new TH1F("hi_dttime_vtx_tb","DT Time at Vertex (BOT-TOP)",theNBins,-100.,140.);
   hi_dttime_vtx_tb2 = new TH2F("hi_dttime_vtx_tb2","DT Time at Vertex (BOT-TOP)",theNBins,-100.,140.,60,-100.,140.);
   hi_dttime_vtxp_t = new TH2F("hi_dttime_vtxp_t","DT Time at Vertex (TOP)",theNBins,-100.,140.,theNBins,0.,3.14);
   hi_dttime_vtxp_b = new TH2F("hi_dttime_vtxp_b","DT Time at Vertex (BOT)",theNBins,-100.,140.,theNBins,0.,3.14);
   hi_dttime_vtxp_tb = new TH2F("hi_dttime_vtxp_tb","DT Time at Vertex (BOT-TOP)",theNBins,-100.,140.,theNBins,0.,3.14);
   hi_dttime_vtxpt_tb = new TH2F("hi_dttime_vtxpt_tb","DT Time at Vertex (BOT-TOP)",theNBins,-100.,140.,theNBins,0.,60.);
   hi_dttime_vtx_err = new TH1F("hi_dttime_vtx_err","DT Time at Vertex Error (inout)",theNBins,0.,10.0);
   hi_dttime_vtxr = new TH1F("hi_dttime_vtxR","DT Time at Vertex (inout)",theNBins,0.,75.*theScale);
   hi_dttime_vtxr_err = new TH1F("hi_dttime_vtxR_err","DT Time at Vertex Error (inout)",theNBins,0.,10.0);
   hi_dttime_ibt_pull = new TH1F("hi_dttime_ibt_pull","DT Inverse Beta Pull",theNBins,-5.,5.0);
   hi_dttime_fib_pull = new TH1F("hi_dttime_fib_pull","DT Free Inverse Beta Pull",theNBins,-5.,5.0);
   hi_dttime_vtx_pull = new TH1F("hi_dttime_vtx_pull","DT Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_dttime_vtxr_pull = new TH1F("hi_dttime_vtxR_pull","DT Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_dttime_errdiff = new TH1F("hi_dttime_errdiff","DT Time at Vertex inout-outin error difference",theNBins,-theScale,theScale);
   hi_dttime_errdiff_t = new TH1F("hi_dttime_errdiff_t","DT Time at Vertex inout-outin error difference (Top)",theNBins,-theScale,theScale);
   hi_dttime_errdiff_b = new TH1F("hi_dttime_errdiff_b","DT Time at Vertex inout-outin error difference (Bot)",theNBins,-theScale,theScale);
   hi_dttime_ndof = new TH1F("hi_dttime_ndof","Number of DT timing measurements",48,0.,48.0);

   hi_csctime_ibt = new TH1F("hi_csctime_ibt","CSC Inverse Beta",theNBins,0.,1.6);
   hi_csctime_ibt_pt = new TH2F("hi_csctime_ibt_pt","P{T} vs CSC Inverse Beta",theNBins,theMinPtres,theMaxPtres,theNBins,0.7,2.0);
   hi_csctime_ibt_err = new TH1F("hi_csctime_ibt_err","CSC Inverse Beta Error",theNBins,0.,1.0);
   hi_csctime_fib = new TH1F("hi_csctime_fib","CSC Free Inverse Beta",theNBins,-5.,7.);
   hi_csctime_fib_t = new TH1F("hi_csctime_fib_t","CSC Free Inverse Beta (TOP)",theNBins,-5.,5.);
   hi_csctime_fib_b = new TH1F("hi_csctime_fib_b","CSC Free Inverse Beta (BOT)",theNBins,-5.,5.);
   hi_csctime_fib_err = new TH1F("hi_csctime_fib_err","CSC Free Inverse Beta Error",theNBins,0,5.);
   hi_csctime_vtx = new TH1F("hi_csctime_vtx","CSC Time at Vertex (inout)",theNBins,-100,100);
   hi_csctime_vtxn = new TH2F("hi_csctime_vtxn","CSC Time at Vertex vs nDof",theNBins,-100,100,48,0.,48.0);
   hi_csctime_vtx_pt = new TH2F("hi_csctime_vtx_pt","Time at Vertex vs STA p_{T}",theNBins,-25.*theScale,25.*theScale,theNBins,theMinPtres,theMaxPtres);
   hi_csctime_vtx_t = new TH1F("hi_csctime_vtx_t","CSC Time at Vertex (TOP inout)",theNBins,-25.*theScale,25.*theScale);
   hi_csctime_vtx_b = new TH1F("hi_csctime_vtx_b","CSC Time at Vertex (BOT inout)",theNBins,-25.*theScale,25.*theScale);
   hi_csctime_vtx_phi = new TH2F("hi_csctime_vtx_phi","CSC Time at Vertex vs Phi",theNBins,-100,100,60,-3.14,3.14);
   hi_csctime_vtx_eta = new TH2F("hi_csctime_vtx_eta","CSC Time at Vertex vs Eta",theNBins,-100,100,60,-2.5,2.5);
   hi_csctime_eeta_lo = new TH2F("hi_csctime_eeta_lo","Pt Eta vs Origin Eta for CSC in-time",60,-2.1,2.1,60,-2.1,2.1);
   hi_csctime_eeta_hi = new TH2F("hi_csctime_eeta_hi","Pt Eta vs Origin Eta for CSC ou-time",60,-2.1,2.1,60,-2.1,2.1);
   hi_csctime_vtx_etat = new TH2F("hi_csctime_vtx_etat","CSC Time at Vertex vs Eta (TOP)",theNBins,-100,100,60,-2.1,2.1);
   hi_csctime_vtx_etab = new TH2F("hi_csctime_vtx_etab","CSC Time at Vertex vs Eta (BOT)",theNBins,-100,100,60,-2.1,2.1);
   hi_csctime_vtx_err = new TH1F("hi_csctime_vtx_err","CSC Time at Vertex Error (inout)",theNBins,0.,25.0);
   hi_csctime_vtxr = new TH1F("hi_csctime_vtxR","CSC Time at Vertex (outin)",theNBins,0.,75.*theScale);
   hi_csctime_vtxr_err = new TH1F("hi_csctime_vtxR_err","CSC Time at Vertex Error (outin)",theNBins,0.,25.0);
   hi_csctime_ibt_pull = new TH1F("hi_csctime_ibt_pull","CSC Inverse Beta Pull",theNBins,-5.,5.0);
   hi_csctime_fib_pull = new TH1F("hi_csctime_fib_pull","CSC Free Inverse Beta Pull",theNBins,-5.,5.0);
   hi_csctime_vtx_pull = new TH1F("hi_csctime_vtx_pull","CSC Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_csctime_vtxr_pull = new TH1F("hi_csctime_vtxR_pull","CSC Time at Vertex Pull (inout)",theNBins,-5.,5.0);
   hi_csctime_ndof = new TH1F("hi_csctime_ndof","Number of CSC timing measurements",48,0.,48.0);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonTimingAnalyzer::endJob() {

  hFile->cd();

  gROOT->SetStyle("effStyle");

  hi_id_rpccut_sta->Write();
  hi_id_rpccut_glb->Write();
  hi_id_csccut_sta->Write();
  hi_id_csccut_glb->Write();
  hi_id_dtcut_sta->Write();
  hi_id_dtcut_glb->Write();

  hi_id_trklay->Write();
  hi_id_trkhit->Write();
  hi_id_statio->Write();
  hi_id_dxy->Write();
  hi_id_dz->Write();

  hi_sta_pt->Write();
  hi_sta_pt_cut->Write();
  hi_sta_ptres->Write();
  hi_sta_ptt->Write();
  hi_sta_phi->Write();
  hi_sta_eta->Write();
//  hi_sta_nhits->Write();
  hi_sta_nvhits->Write();
  hi_sta_chi2->Write();

  hi_tk_pt->Write();
  hi_tk_phi->Write();
  hi_tk_eta->Write();
//  hi_tk_nhits->Write();
  hi_tk_nvhits->Write();
  hi_tk_chi2->Write();

  hi_glb_pt->Write();
  hi_glb_pt_cut->Write();
//  hi_glb_ptres->Write();
//  hi_glb_ptt->Write();
  hi_glb_phi->Write();
  hi_glb_eta->Write();
//  hi_glb_nhits->Write();
  hi_glb_nvhits->Write();
//  hi_glb_ptresh->Write();
//  hi_glb_ptres_t->Write();
//  hi_glb_ptresh_t->Write();
//  hi_glb_ptres_b->Write();
//  hi_glb_ptresh_b->Write();
//  hi_glb_ptres_tb->Write();
  hi_glb_d0->Write();
  hi_glb_chi2->Write();

  hi_trpc->Write();
  hi_nrpc->Write();
  hi_trpc_eta->Write();
  hi_trpc_phi->Write();
  hi_mutime_vtx->Write();
  hi_mutime_vtx_err->Write();

  hFile->mkdir("differences");
  hFile->cd("differences");

  hi_dtcsc_vtx->Write();
  hi_dtcsc_vtx_t->Write();
  hi_dtcsc_vtx_b->Write();

  hFile->cd();
  hFile->mkdir("combined");
  hFile->cd("combined");

  hi_cmbtime_ibt->Write();
  hi_cmbtime_ibt_pt->Write();
  hi_cmbtime_ibt_err->Write();
  hi_cmbtime_fib->Write();
  hi_cmbtime_fib_err->Write();
  hi_cmbtime_vtx->Write();
  hi_cmbtime_vtx_err->Write();
  hi_cmbtime_vtxr->Write();
  hi_cmbtime_vtxr_err->Write();
  hi_cmbtime_ibt_pull->Write();
  hi_cmbtime_fib_pull->Write();
  hi_cmbtime_vtx_pull->Write();
  hi_cmbtime_vtxr_pull->Write();
  hi_cmbtime_ndof->Write();

  hFile->cd();
  hFile->mkdir("dt");
  hFile->cd("dt");

  hi_dttime_ibt->Write();
  hi_dttime_ibt_pt->Write();
  hi_dttime_ibt_err->Write();
  hi_dttime_fib->Write();
//  hi_dttime_fib_t->Write();
//  hi_dttime_fib_b->Write();
//  hi_dttime_fibp_t->Write();
//  hi_dttime_fibp_b->Write();
  hi_dttime_fib_err->Write();
  hi_dttime_vtx->Write();
  hi_dttime_vtxn->Write();
  hi_dttime_vtx_w->Write();
  hi_dttime_vtx_pt->Write();
  hi_dttime_vtx_phi->Write();
  hi_dttime_vtx_eta->Write();
  hi_dttime_etaphi->Write();
//  hi_dttime_eeta_lo->Write();
//  hi_dttime_eeta_hi->Write();
//  hi_dttime_vtx_etat->Write();
//  hi_dttime_vtx_etab->Write();
//  hi_dttime_vtx_t->Write();
//  hi_dttime_vtx_b->Write();
//  hi_dttime_vtx_to->Write();
//  hi_dttime_vtx_bo->Write();
  hi_dttime_vtx_tb->Write();
  hi_dttime_vtx_tb2->Write();
//  hi_dttime_vtxp_t->Write();
//  hi_dttime_vtxp_b->Write();
//  hi_dttime_vtxp_tb->Write();
//  hi_dttime_vtxpt_tb->Write();
  hi_dttime_vtx_err->Write();
  hi_dttime_vtxr->Write();
  hi_dttime_vtxr_err->Write();
  hi_dttime_ibt_pull->Write();
  hi_dttime_fib_pull->Write();
  hi_dttime_vtx_pull->Write();
  hi_dttime_vtxr_pull->Write();
  hi_dttime_errdiff->Write();
  hi_dttime_errdiff_t->Write();
  hi_dttime_errdiff_b->Write();
  hi_dttime_ndof->Write();

  hFile->cd();
  hFile->mkdir("csc");
  hFile->cd("csc");

  hi_csctime_ibt->Write();
  hi_csctime_ibt_pt->Write();
  hi_csctime_ibt_err->Write();
  hi_csctime_fib->Write();
//  hi_csctime_fib_t->Write();
//  hi_csctime_fib_b->Write();
  hi_csctime_fib_err->Write();
  hi_csctime_vtx->Write();
  hi_csctime_vtxn->Write();
  hi_csctime_vtx_pt->Write();
//  hi_csctime_vtx_t->Write();
//  hi_csctime_vtx_b->Write();
  hi_csctime_vtx_eta->Write();
  hi_csctime_vtx_phi->Write();
//  hi_csctime_eeta_lo->Write();
//  hi_csctime_eeta_hi->Write();
  hi_csctime_vtx_etat->Write();
  hi_csctime_vtx_etab->Write();
  hi_csctime_vtx_err->Write();
  hi_csctime_vtxr->Write();
  hi_csctime_vtxr_err->Write();
  hi_csctime_ibt_pull->Write();
  hi_csctime_fib_pull->Write();
  hi_csctime_vtx_pull->Write();
  hi_csctime_vtxr_pull->Write();
  hi_csctime_ndof->Write();

  hFile->Write();
}

float 
MuonTimingAnalyzer::calculateDistance(const math::XYZVector& vect1, const math::XYZVector& vect2) {
  float dEta = vect1.eta() - vect2.eta();
  float dPhi = fabs(Geom::Phi<float>(vect1.phi()) - Geom::Phi<float>(vect2.phi()));
  float distance = sqrt(pow(dEta,2) + pow(dPhi,2) );

  return distance;
}

//
// return h1/h2 with recalculated errors
//
TH1F* MuonTimingAnalyzer::divideErr(TH1F* h1, TH1F* h2, TH1F* hout) {

  hout->Reset();
  hout->Divide(h1,h2,1.,1.,"B");

  for (int i = 0; i <= hout->GetNbinsX()+1; i++ ) {
    Float_t tot   = h2->GetBinContent(i) ;
    Float_t tot_e = h2->GetBinError(i);
    Float_t eff = hout->GetBinContent(i) ;
    Float_t Err = 0.;
    if (tot > 0) Err = tot_e / tot * sqrt( eff* (1-eff) );
    if (eff == 1. || isnan(Err) || !isfinite(Err) ) Err=1.e-3;
    hout->SetBinError(i, Err);
  }
  return hout;
}

double MuonTimingAnalyzer::iMass(reco::TrackRef imuon, reco::TrackRef iimuon) {
  double energy1 = sqrt(imuon->p() * imuon->p() + 0.011163691);
  double energy2 = sqrt(iimuon->p() * iimuon->p() + 0.011163691);
  math::XYZTLorentzVector ip4(imuon->px(),imuon->py(),imuon->pz(),energy1);
  math::XYZTLorentzVector iip4(iimuon->px(),iimuon->py(),iimuon->pz(),energy2);
  math::XYZTLorentzVector psum = ip4+iip4;
  double mmumu2 = psum.Dot(psum);
  return sqrt(mmumu2);
}

bool MuonTimingAnalyzer::dumpMuonId(const reco::Muon& muon, const reco::Vertex& vtx, const bool debug){
  if (debug) {
    cout << " isPFmuon " << muon.isPFMuon() << "   isGlobalMuon " << muon.isGlobalMuon() << endl;
    cout << " matched Stations (>1): " << muon.numberOfMatchedStations() << endl;
    cout << " trk Layers (>5): " << muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() << endl;
    cout << " pixel Hits (>0): " << muon.innerTrack()->hitPattern().numberOfValidPixelHits() << endl;
    cout << " dxy (<0.2): " << fabs(muon.muonBestTrack()->dxy(vtx.position())) << "  vtx: " << vtx.position() << endl;
    cout << " dz (<0.5): " << fabs(muon.muonBestTrack()->dz(vtx.position())) << endl;
  }
  if (!muon.isPFMuon() || !muon.isGlobalMuon() || (muon.numberOfMatchedStations()<2)) return(false);
  if ((muon.innerTrack()->hitPattern().trackerLayersWithMeasurement()<6) || 
      (muon.innerTrack()->hitPattern().numberOfValidPixelHits()==0)) return(false);
  return(true);
}

void MuonTimingAnalyzer::dumpTrack(reco::TrackRef track) {
  cout 
       << "|  " << fixed << setprecision(0) << track->pt() << " +/- " << track->ptError()
       << " |  " << track->p() 
       << " |  " << setprecision(2) << track->eta()
       << " |  " << track->phi() 
//       << " |  " << setprecision(3) << track->dxy(beamspot)
       << " |  " << setprecision(2) << track->normalizedChi2() 
//       << " |  " << track->recHitsSize() 
       << " |  " << track->found()
       << " |  " << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTimingAnalyzer);
