// -*- C++ -*-
//
// Package:    AODTimingAnalyzer
// Class:      AODTimingAnalyzer
// 
/**\class AODTimingAnalyzer AODTimingAnalyzer.cc 

 Description: Fill muon timing information histograms 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: AODTimingAnalyzer.cc,v 1.5 2011/04/06 11:44:29 ptraczyk Exp $
//
//

#include "AodTimingAnalyzer.h"

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
#include "DataFormats/PatCandidates/interface/Muon.h"

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
AODTimingAnalyzer::AODTimingAnalyzer(const edm::ParameterSet& iConfig) 
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  VtxTags_(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")),
  TimeTags_(iConfig.getUntrackedParameter<edm::InputTag>("Timing")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theDebug(iConfig.getParameter<bool>("debug")),
  doSim(iConfig.getParameter<bool>("mctruthMatching")),
  theIdCut(iConfig.getParameter<string>("requireId")),
  theCollVeto(iConfig.getParameter<bool>("collisionVeto")),
  theKeepBX(iConfig.getParameter<bool>("keepOnlyBX")),
  theBX(iConfig.getParameter<int>("generatedBX")),
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
  edm::ConsumesCollector collector(consumesCollector());
  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trackToken_ = consumes<reco::TrackCollection>(TKtrackTags_);
  muonToken_ = consumes<pat::MuonCollection>(MuonTags_);
  vertexToken_ = consumes<reco::VertexCollection>(edm::InputTag(VtxTags_));
  genParticleToken_ = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  trackingParticleToken_ = consumes<TrackingParticleCollection>(edm::InputTag("mix","MergedTrackTruth"));
}


AODTimingAnalyzer::~AODTimingAnalyzer()
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
AODTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using pat::MuonCollection;

  bool debug=theDebug;
  bool tpart=false;

  if (debug) {
//    cout << "*** Begin Muon Timing Analyzer " << endl;
    cout << endl << " Event: " << iEvent.id() << "  Orbit: " << iEvent.orbitNumber() << "  BX: " << iEvent.bunchCrossing() << endl;
  }

  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  edm::Handle<reco::VertexCollection> recVtxs;
  cout << vertexToken_.index() << endl;
  iEvent.getByToken(vertexToken_,recVtxs);
  unsigned int theIndexOfThePrimaryVertex = 999.;
  for (unsigned int ind=0; ind<recVtxs->size(); ++ind) {
    if ( (*recVtxs)[ind].isValid() ) {
      theIndexOfThePrimaryVertex = ind;
      break;
    }
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
  iEvent.getByToken(trackingParticleToken_, TruthTrackContainer );
  if (!TruthTrackContainer.isValid()) {
    if (debug) cout << "No trackingparticle data in the Event" << endl;
  } else tpart=true;
  
  const TrackingParticleCollection *tPC=0;
  if (tpart) 
    tPC = TruthTrackContainer.product();
  
  //edm::Handle<reco::TrackCollection> trackc;
  //iEvent.getByToken( trackToken_, trackc);
  //const reco::TrackCollection trackC = *(trackc.product());

  // simple "collision event" veto on number of tracker tracks greater than 2
  //if (theCollVeto && trackC.size()>2) return;

  // Generated particle collection
  Handle<GenParticleCollection> genParticles;
  if (doSim)   
    iEvent.getByToken(genParticleToken_, genParticles);

  // Fill generated muon information
  if (tpart)
    for (TrackingParticleCollection::const_iterator iTrack = tPC->begin(); iTrack != tPC->end(); ++iTrack)
      if (fabs(iTrack->pdgId())==13 && iTrack->p4().Pt()>2 && 
         (fabs(iTrack->p4().eta())<2.5) && ((iTrack->eventId().bunchCrossing()==theBX) || !theKeepBX)) {
        hi_gen_pt->Fill(iTrack->p4().Pt());
        hi_gen_eta->Fill(iTrack->p4().Eta());
        hi_gen_phi->Fill(iTrack->p4().Phi());
      }

  iEvent.getByToken(muonToken_,MuCollection);
  const pat::MuonCollection muonC = *(MuCollection.product());
  if (debug) cout << " Muon collection size: " << muonC.size() << endl;
  if (!muonC.size()) return;
  MuonCollection::const_iterator imuon,iimuon;

  // check for back-to-back dimuons
  float angle=0;
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

  int imucount=0;
  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    double outeta=0;
    double genpt=0;
    bool matched=false;

    reco::TrackRef glbTrack = imuon->combinedMuon();
    reco::TrackRef trkTrack = imuon->track();
    reco::TrackRef staTrack = imuon->standAloneMuon();

    pat::MuonRef muonR(MuCollection,imucount);
    imucount++;    
    
    reco::MuonTime timec = imuon->time();
    reco::MuonTime rpcTime = imuon->rpcTime();
    if (rpcTime.timeAtIpInOut<-60 || rpcTime.timeAtIpInOut>80) rpcTime.nDof=0;
        
    bool idcut = true;
    if (rpcTime.nDof>1 && rpcTime.timeAtIpInOutErr<1) {
      if (fabs(rpcTime.timeAtIpInOut)>20) idcut=false;     
    } else 
      if (timec.nDof>4 && ((timec.timeAtIpInOut>20) || (timec.timeAtIpInOut<-50))) idcut=false;
        
    if (imuon->pt()<thePtCut) continue;
    if ((fabs(imuon->eta())<theMinEta) || (fabs(imuon->eta())>theMaxEta)) continue;
    if (theIdCut=="glb"   && !muon::isGoodMuon(*imuon, muon::GlobalMuonPromptTight )) continue;
    if (theIdCut=="loose" && !muon::isLooseMuon(*imuon)) continue;
    if (theIdCut=="tight" && !muon::isTightMuon(*imuon, pvertex )) continue;
    if (theIdCut=="norpc" && rpcTime.nDof>1) continue;
    if (theIdCut=="norpc3" && rpcTime.nDof>1 && rpcTime.timeAtIpInOutErr==0) continue;
    if (theIdCut=="timeok" && !idcut) continue;

    if (theVetoCosmics && fabs(imuon->muonBestTrack()->dxy(pvertex.position()))>0.1) continue;
    if (theVetoCosmics && fabs(imuon->muonBestTrack()->dz(pvertex.position()))>1.) continue;

    if (debug) cout << endl << "   Found muon. Pt: " << imuon->pt() << endl;

    if (muon::isGoodMuon(*imuon, muon::GlobalMuonPromptTight )) {
      hi_id_trklay->Fill(imuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      hi_id_trkhit->Fill(imuon->innerTrack()->hitPattern().numberOfValidPixelHits());
      hi_id_statio->Fill(imuon->numberOfMatchedStations());
      hi_id_dxy->Fill(fabs(imuon->muonBestTrack()->dxy(pvertex.position())));
      hi_id_dz->Fill(fabs(imuon->muonBestTrack()->dz(pvertex.position())));
    }

    double d0=0.;
    if (glbTrack.isNonnull()) d0 = -1.*glbTrack->dxy(beamspot);    

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
            genpt=iTrack->p4().Pt();
            if (debug) {
              cout << " Matched muon BX: " << iTrack->eventId().bunchCrossing();
              cout << "  hits: " << iTrack->numberOfTrackerLayers();
              cout << "  pT: " << iTrack->p4().Pt() << endl;
	    }
            if ((iTrack->eventId().bunchCrossing()!=theBX) && theKeepBX) matched=false;
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
      hi_tk_nvhits->Fill((*trkTrack).found());
    }  

    if (staTrack.isNonnull()) {
      hi_sta_pt->Fill((*staTrack).pt());
      hi_sta_ptres->Fill((*staTrack).pt()-genpt);
      hi_sta_phi->Fill((*staTrack).phi());
      hi_sta_eta->Fill((*staTrack).eta());
      hi_sta_chi2->Fill((*staTrack).normalizedChi2());
      hi_sta_nvhits->Fill((*staTrack).found());
    }  

    if (glbTrack.isNonnull()) {
      hi_glb_pt->Fill(imuon->pt());
      hi_glb_ptres->Fill(imuon->pt()-genpt);
      hi_glb_eta->Fill(imuon->eta());
      hi_glb_d0->Fill(d0);
      hi_glb_phi->Fill(imuon->phi());
      hi_glb_chi2->Fill(((*glbTrack).normalizedChi2()));
      hi_glb_nvhits->Fill((*glbTrack).found());
    }

    if (debug) {
      cout << endl;
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
/*
        cout << "|     Global Track ";
        dumpTrack(fmsTrack);      
        cout << "|     FMS Track ";
        dumpTrack(pmrTrack);      
        cout << "|     Picky Track ";
        dumpTrack(dytTrack);      
        cout << "|     DYT Track ";
        dumpTrack(glbTrack);      
*/
      }
      cout << endl;
    }
    
    // Analyze the short info stored directly in reco::Muon
    
    reco::MuonTime muonTime;
    if (imuon->isTimeValid()) { 
      muonTime = imuon->time();
      if (debug) cout << "    Time points: " << muonTime.nDof << "  time: " << muonTime.timeAtIpInOut << endl;
      hi_mutime_ndof->Fill(muonTime.nDof);
      if (muonTime.nDof>4) {
        hi_mutime_vtx->Fill(muonTime.timeAtIpInOut);
        hi_mutime_vtx_err->Fill(muonTime.timeAtIpInOutErr);
      }
    }
    
    hi_nrpc->Fill(rpcTime.nDof);
    if (rpcTime.nDof>0) {
      hi_trpc->Fill(rpcTime.timeAtIpInOut);
      hi_trpcerr->Fill(rpcTime.timeAtIpInOutErr);
      if (rpcTime.nDof>1 && rpcTime.timeAtIpInOutErr<1) hi_trpc3->Fill(rpcTime.timeAtIpInOut);
      hi_nrpc_trpc->Fill(rpcTime.nDof,rpcTime.timeAtIpInOut);
      hi_trpc_phi->Fill(rpcTime.timeAtIpInOut,imuon->phi());
      hi_trpc_eta->Fill(rpcTime.timeAtIpInOut,imuon->eta());
      if (debug) cout << "   RPC time: " << rpcTime.timeAtIpInOut << " +/- " << rpcTime.timeAtIpInOutErr << endl;
    }

    bool timeok = false;

    if (timeok) {
      if (staTrack.isNonnull()) hi_sta_ptt->Fill((*staTrack).pt());
      if (glbTrack.isNonnull()) hi_glb_ptt->Fill(imuon->pt());
    }
  }  

}


// ------------ method called once each job just before starting event loop  ------------
void 
AODTimingAnalyzer::beginJob()
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

   hi_gen_pt  = new TH1F("hi_gen_pt","P_{T}^{GEN}",theNBins,theMinPtres,theMaxPtres);
   hi_gen_phi = new TH1F("hi_gen_phi","#phi^{GEN}",theNBins,-3.0,3.);
   hi_gen_eta = new TH1F("hi_gen_eta","#eta^{GEN}",theNBins/2,2.5,2.5);
   
   hi_id_rpccut_sta = new TH1F("hi_id_rpccut_sta","STA muons rejected by RPC cut",50,0,50);
   hi_id_rpccut_glb = new TH1F("hi_id_rpccut_glb","GLB muons rejected by RPC cut",50,0,50);
   hi_id_csccut_sta = new TH1F("hi_id_csccut_sta","STA muons rejected by CSC cut",50,0,50);
   hi_id_csccut_glb = new TH1F("hi_id_csccut_glb","GLB muons rejected by CSC cut",50,0,50);
   hi_id_dtcut_sta = new TH2F("hi_id_dtcut_sta","STA muons rejected by DT cut",50,0,50,15,0,15);
   hi_id_dtcut_glb = new TH2F("hi_id_dtcut_glb","GLB muons rejected by DT cut",50,0,50,15,0,15);
   hi_id_cmbcut_sta = new TH2F("hi_id_cmbcut_sta","STA muons rejected by CMB cut",50,0,50,15,0,15);
   hi_id_cmbcut_glb = new TH2F("hi_id_cmbcut_glb","GLB muons rejected by CMB cut",50,0,50,15,0,15);

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

   hi_mutime_vtx = new TH1F("hi_mutime_vtx","Time at Vertex (inout)",theNBins,-100.,100.);
   hi_mutime_vtx_err = new TH1F("hi_mutime_vtx_err","Time at Vertex Error (inout)",theNBins,0.,25.0);
   hi_mutime_ndof = new TH1F("hi_mutime_ndof","Number of timing measurements",60,0.,60.0);

   hi_trpc = new TH1F("hi_trpc","Time at Vertex (RPC)",theNBins,-100.,100.);
   hi_trpc3= new TH1F("hi_trpc3","Time at Vertex (RPC, RPC nHits>1 RPCerr=0) ",theNBins,-100.,100.);
   hi_nrpc = new TH1F("hi_nrpc","RPC nHits",8,0,8);
   hi_trpcerr = new TH1F("hi_trpcerr","Time at Vertex Error (RPC)",theNBins,0.,25.);
   hi_nrpc_trpc = new TH2F("hi_nrpc_trpc","RPC nHits vs time",8,0,8,theNBins,-100.,100.);
   hi_trpc_phi = new TH2F("hi_trpc_phi","RPC Time vs Phi",theNBins,-100.,100.,60,-3.14,3.14);
   hi_trpc_eta = new TH2F("hi_trpc_eta","RPC Time vs Eta",theNBins,-100.,100.,60,-2.5,2.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AODTimingAnalyzer::endJob() {

  hFile->cd();

  gROOT->SetStyle("effStyle");

  hi_gen_pt->Write();
  hi_gen_phi->Write();
  hi_gen_eta->Write();

  hi_id_rpccut_sta->Write();
  hi_id_rpccut_glb->Write();
  hi_id_csccut_sta->Write();
  hi_id_csccut_glb->Write();
  hi_id_dtcut_sta->Write();
  hi_id_dtcut_glb->Write();
  hi_id_cmbcut_sta->Write();
  hi_id_cmbcut_glb->Write();

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
  hi_trpc3->Write();
  hi_trpcerr->Write();
  hi_nrpc->Write();
  hi_nrpc_trpc->Write();
  hi_trpc_eta->Write();
  hi_trpc_phi->Write();
  hi_mutime_ndof->Write();
  hi_mutime_vtx->Write();
  hi_mutime_vtx_err->Write();

  hFile->Write();
}


double AODTimingAnalyzer::iMass(reco::TrackRef imuon, reco::TrackRef iimuon) {
  double energy1 = sqrt(imuon->p() * imuon->p() + 0.011163691);
  double energy2 = sqrt(iimuon->p() * iimuon->p() + 0.011163691);
  math::XYZTLorentzVector ip4(imuon->px(),imuon->py(),imuon->pz(),energy1);
  math::XYZTLorentzVector iip4(iimuon->px(),iimuon->py(),iimuon->pz(),energy2);
  math::XYZTLorentzVector psum = ip4+iip4;
  double mmumu2 = psum.Dot(psum);
  return sqrt(mmumu2);
}

bool AODTimingAnalyzer::dumpMuonId(const pat::Muon& muon, const reco::Vertex& vtx, const bool debug){
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

void AODTimingAnalyzer::dumpTrack(reco::TrackRef track) {
  if (!track.isNonnull()) return;
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
DEFINE_FWK_MODULE(AODTimingAnalyzer);
