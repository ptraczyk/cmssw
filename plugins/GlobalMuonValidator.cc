// -*- C++ -*-
//
// Package:    GlobalMuonValidator
// Class:      GlobalMuonValidator
// 
/**\class GlobalMuonValidator GlobalMuonValidator.cc AEverett/GlobalMuonValidator/src/GlobalMuonValidator.cc

 Description: Validator tool to study efficiencies and purities.

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: GlobalMuonValidator.cc,v 1.4 2012/08/01 08:13:29 ptraczyk Exp $
//
//

#include "UserCode/HSCPTOF/plugins/GlobalMuonValidator.h"

// system include files
#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

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
#include "RecoLocalMuon/DTSegment/src/DTSegmentCand.h"

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

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
GlobalMuonValidator::GlobalMuonValidator(const edm::ParameterSet& iConfig) 
  :
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theMinEta(iConfig.getParameter<double>("etaMin")),
  theMaxEta(iConfig.getParameter<double>("etaMax")),
  theMinPt(iConfig.getParameter<double>("simPtMin")),
  thePtCut(iConfig.getParameter<double>("PtCut")),
  theMinPtres(iConfig.getParameter<double>("PtresMin")),
  theMaxPtres(iConfig.getParameter<double>("PtresMax")),
  theInvPt(iConfig.getParameter<double>("invPtScale")),
  theNBins(iConfig.getParameter<int>("nbins")),
  theDTRecHitLabel(iConfig.getUntrackedParameter<edm::InputTag>("DTRecHits")),
  theCSCRecHitLabel(iConfig.getUntrackedParameter<edm::InputTag>("CSCRecHits"))
{
  //now do what ever initialization is needed

  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
  // the services
  theService = new MuonServiceProxy(serviceParameters);
  std::ofstream out("Muon_reco.txt",ios::trunc);
  
}


GlobalMuonValidator::~GlobalMuonValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
  if (theService) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GlobalMuonValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::MuonCollection;

  int theHitCut = 2;
  std::ofstream out("Muon_reco.txt",ios::app);

//  cout << "*** Begin Muon Validatior " << endl;

  // Update the services
  theService->update(iSetup);

  iEvent.getByLabel(MuonTags_,MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());

  iEvent.getByLabel(theDTRecHitLabel, theDTRecHits);
  iEvent.getByLabel(theCSCRecHitLabel, theCSCRecHits);

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  MuonCollection::const_iterator imuon;
  if (!muonC.size()) return;

  for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
    
    if (!muon::isGoodMuon(*imuon, muon::GlobalMuonPromptTight )) continue;
    if ((fabs(imuon->eta())<theMinEta) || (fabs(imuon->eta())>theMaxEta)) continue;
    
    hi_glb6_pt->Fill(imuon->pt());
    
    cout << endl << " Found muon with pT: " << imuon->pt();
        
    reco::TrackRef trkTrack = imuon->track();
    if (trkTrack.isNonnull()) { 
      hi_tk_pt->Fill(((*trkTrack).pt()));
      hi_tk_eta->Fill(((*trkTrack).eta()));
    }  

    reco::TrackRef staTrack = imuon->standAloneMuon();
    if (staTrack.isNonnull()) {
      hi_sta_pt->Fill((*staTrack).pt());
      hi_sta_eta->Fill(((*staTrack).eta()));
    }

    reco::TrackRef glbTrack = imuon->combinedMuon();
    if (glbTrack.isNonnull()) {
      hi_glb_pt->Fill((*glbTrack).pt());
      hi_glb_eta->Fill((*glbTrack).eta()); 
    }  
    
    if ((glbTrack.isNonnull()) && (staTrack.isNonnull())) {

      reco::TrackRef fmsTrack = imuon->tpfmsTrack();
      reco::TrackRef pmrTrack = imuon->pickyTrack();
      reco::TrackRef dytTrack = imuon->dytTrack();
   
      if (fmsTrack.isNonnull()) cout << " FMS: " << fmsTrack->pt();
      if (pmrTrack.isNonnull()) cout << " PMR: " << pmrTrack->pt();
      if (dytTrack.isNonnull()) cout << " DYT: " << dytTrack->pt();

      if (fmsTrack.isNonnull()) {
        hi_glb2_pt->Fill(fmsTrack->pt());
        hi_glb2_prob->Fill(trackProbability(*(fmsTrack)));
      }    
      if (pmrTrack.isNonnull()) { 
        hi_glb3_pt->Fill(pmrTrack->pt());
        hi_glb3_prob->Fill(trackProbability(*(pmrTrack)));
      }    
      if (dytTrack.isNonnull()) { 
        hi_glb4_pt->Fill(dytTrack->pt());
        hi_glb4_prob->Fill(trackProbability(*(dytTrack)));
      }    
  
      reco::TrackRef cktTrack = (muon::tevOptimized(*imuon,200.,4.,6.)).first;
      if (cktTrack.isNonnull()) {
        hi_glb5_pt->Fill(cktTrack->pt());
        hi_glb5_prob->Fill(trackProbability(*cktTrack));
      }    
    
      int weird=0;
      if ((glbTrack->pt()<200) || (glbTrack->pt()>1500)) weird=1;

      vector<int> dethits(4,0);
      checkMuonHits(*(imuon->combinedMuon()),dethits);
      int num_sho=0;
      int num_stations=0;
      for (int i=0;i<4;i++) {
        if (dethits[i]>theHitCut) num_sho++;
        if (dethits[i]>0) num_stations++;
        hi_num_hits->Fill(dethits[i]);
      }
      
      hi_num_sho->Fill(num_sho);
      hi_sho_pt->Fill(glbTrack->pt(),num_sho);
      hi_sho_p->Fill(glbTrack->p(),num_sho);
      hi_sho_eta->Fill(glbTrack->eta(),num_sho);
      
      if (num_stations>3) {
      if (!num_sho) {
        hi_glb_pt_0->Fill(glbTrack->pt());
        hi_sta_pt_0->Fill(staTrack->pt());
        if (weird) weird=2;
      } else if (num_sho==1) {
        hi_glb_pt_1->Fill(glbTrack->pt());
        hi_sta_pt_1->Fill(staTrack->pt());
      } else if (num_sho==2) {
        hi_glb_pt_2->Fill(glbTrack->pt());
        hi_sta_pt_2->Fill(staTrack->pt());
      } else if (num_sho==3) {
        hi_glb_pt_3->Fill(glbTrack->pt());
        hi_sta_pt_3->Fill(staTrack->pt());
      } else if (num_sho==4) {
        hi_glb_pt_4->Fill(glbTrack->pt());
        hi_sta_pt_4->Fill(staTrack->pt());
      } 
      }

      if (weird==2) {
        cout << " *** WEIRD EVENT *** " << endl;
        cout << " run: " << iEvent.id().run() << "   event: " << iEvent.id().event() << endl;
        cout << " Showers: " << num_sho << "   Stations: " << num_stations << endl;
        cout << " GLB: " << glbTrack->pt() << endl;
        cout << " PMC: " << cktTrack->pt() << endl;
      }
    }
    
    reco::MuonTime muonTime;
    if (imuon->isTimeValid()) { 
      muonTime = imuon->time();
      if (muonTime.nDof) {
        hi_time_vtx->Fill(muonTime.timeAtIpInOut);
        hi_time_vtx_err->Fill(muonTime.timeAtIpInOutErr);
        hi_time_nstat->Fill(muonTime.nDof);
      }
    }  
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
GlobalMuonValidator::beginJob()
{
   char title[80];

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

   hi_sta_pt   = new TH1F("hi_sta_pt","P_{T}^{STA}",theNBins,0.0,theMaxPtres);
   hi_tk_pt   = new TH1F("hi_tk_pt","P_{T}^{TK}",theNBins,0.0,theMaxPtres);
   hi_glb_pt   = new TH1F("hi_glb_pt","P_{T}^{GLB}",theNBins,0.0,theMaxPtres);
   hi_glb_diff   = new TH1F("hi_glb_diff","P_{T}^{GLB}",theNBins,-1000.0,1000.0);
   hi_glb1_pt   = new TH1F("hi_glb1_pt","P_{T}^{GLB refit}",theNBins,0.0,theMaxPtres);
   hi_glb2_pt   = new TH1F("hi_glb2_pt","P_{T}^{FMS refit}",theNBins,0.0,theMaxPtres);
   hi_glb3_pt   = new TH1F("hi_glb3_pt","P_{T}^{PMR refit}",theNBins,0.0,theMaxPtres);
   hi_glb4_pt   = new TH1F("hi_glb4_pt","P_{T}^{DYT refit}",theNBins,0.0,theMaxPtres);
   hi_glb5_pt   = new TH1F("hi_glb5_pt","P_{T}^{P cocktail}",theNBins,0.0,theMaxPtres);
   hi_glb6_pt   = new TH1F("hi_glb6_pt","P_{T}^{P muon}",theNBins,0.0,theMaxPtres);
   hi_glbu4_pt   = new TH1F("hi_glbu4_pt","P_{T}^{N local cocktail}",theNBins,0.0,theMaxPtres);
   hi_glbu5_pt   = new TH1F("hi_glbu5_pt","P_{T}^{P local cocktail}",theNBins,0.0,theMaxPtres);

   hi_glb1_prob = new TH1F("hi_glb1_prob","#chi^{2} probability {GLB refit}",theNBins,0.0,1.);
   hi_glb2_prob = new TH1F("hi_glb2_prob","#chi^{2} probability {FMS refit}",theNBins,0.0,1.);
   hi_glb3_prob = new TH1F("hi_glb3_prob","#chi^{2} probability {PMR refit}",theNBins,0.0,1.);
   hi_glb4_prob = new TH1F("hi_glb4_prob","#chi^{2} probability {DYT refit}",theNBins,0.0,1.);
   hi_glb5_prob = new TH1F("hi_glb5_prob","#chi^{2} probability {P cocktail}",theNBins,0.0,1.);
   hi_glb6_prob = new TH1F("hi_glb6_prob","#chi^{2} probability {P muon}",theNBins,0.0,1.);

   hi_sta_pt_0 = new TH1F("hi_sta_pt_0","P_{T}^{STA} no showers",theNBins,0.0,theMaxPtres);
   hi_glb_pt_0 = new TH1F("hi_glb_pt_0","P_{T}^{GLB} no showers",theNBins,0.0,theMaxPtres);
   hi_sta_pt_1 = new TH1F("hi_sta_pt_1","P_{T}^{STA} 1 shower",theNBins,0.0,theMaxPtres);
   hi_glb_pt_1 = new TH1F("hi_glb_pt_1","P_{T}^{GLB} 1 shower",theNBins,0.0,theMaxPtres);
   hi_sta_pt_2 = new TH1F("hi_sta_pt_2","P_{T}^{STA} 2 showers",theNBins,0.0,theMaxPtres);
   hi_glb_pt_2 = new TH1F("hi_glb_pt_2","P_{T}^{GLB} 2 showers",theNBins,0.0,theMaxPtres);
   hi_sta_pt_3 = new TH1F("hi_sta_pt_3","P_{T}^{STA} 3 showers",theNBins,0.0,theMaxPtres);
   hi_glb_pt_3 = new TH1F("hi_glb_pt_3","P_{T}^{GLB} 3 showers",theNBins,0.0,theMaxPtres);
   hi_sta_pt_4 = new TH1F("hi_sta_pt_4","P_{T}^{STA} 4 showers",theNBins,0.0,theMaxPtres);
   hi_glb_pt_4 = new TH1F("hi_glb_pt_4","P_{T}^{GLB} 4 showers",theNBins,0.0,theMaxPtres);

   for (int j=0;j<9;j++)
   for (int i=0;i<9;i++) {
     sprintf(title,"hi_tune_%i_%i",i,j);
     hi_tune[i][j] = new TH1F(title,"Cocktail tune",theNBins,0.0,theMaxPtres);
   }

   hi_time_vtx = new TH1F("hi_time_vtx","Time at Vertex (inout)",100,-25.0,25.0);
   hi_time_vtx_err = new TH1F("hi_time_vtx_err","Time at Vertex Error (inout)",100,0.,25.0);
   hi_time_nstat = new TH1F("hi_time_nstat","Number of stations with timing info",8,0.,8.0);

   hi_sta_eta = new TH1F("hi_sta_eta","#eta^{STA}",theNBins/2,theMinEta,theMaxEta);
   hi_tk_eta  = new TH1F("hi_tk_eta","#eta^{TK}",theNBins/2,theMinEta,theMaxEta);
   hi_glb_eta = new TH1F("hi_glb_eta","#eta^{GLB}",theNBins/2,theMinEta,theMaxEta);

  hi_num_sho  = new TH1F("hi_num_sho","Number of showered chambers",5,0,5);
  hi_num_hits = new TH1F("hi_num_hits","Number of hits per chamber (max over layers)",20,0,20);
  hi_sho_pt   = new TH2F("hi_sho_pt","Number of showered chambers vs reco P_{T}",theNBins,0.0,theMaxPtres,5,0,5);
  hi_sho_p    = new TH2F("hi_sho_p","Number of showered chambers vs reco P",theNBins,0.0,theMaxPtres,5,0,5);
  hi_sho_eta  = new TH2F("hi_sho_eta","Number of showered chambers vs reco Eta",theNBins/4,0.,theMaxEta,5,0,5);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlobalMuonValidator::endJob() {

  hFile->cd();

  gROOT->SetStyle("effStyle");

  hi_sta_pt->Write();
  hi_tk_pt->Write();
  hi_glb_pt->Write();
  hi_glb_diff->Write();
  hi_glb1_pt->Write();
  hi_glb2_pt->Write();
  hi_glb3_pt->Write();
  hi_glb4_pt->Write();
  hi_glb5_pt->Write();
  hi_glb6_pt->Write();
  hi_glbu4_pt->Write();
  hi_glbu5_pt->Write();
  hi_glb1_prob->Write();
  hi_glb2_prob->Write();
  hi_glb3_prob->Write();
  hi_glb4_prob->Write();
  hi_glb5_prob->Write();
  hi_glb6_prob->Write();

  hi_sta_pt_0->Write();
  hi_glb_pt_0->Write();
  hi_sta_pt_1->Write();
  hi_glb_pt_1->Write();
  hi_sta_pt_2->Write();
  hi_glb_pt_2->Write();
  hi_sta_pt_3->Write();
  hi_glb_pt_3->Write();
  hi_sta_pt_4->Write();
  hi_glb_pt_4->Write();

//  for (int j=0;j<9;j++)
//    for (int i=0;i<9;i++) 
//      hi_tune[i][j]->Write();

  hi_sta_eta->Write();
  hi_tk_eta->Write();
  hi_glb_eta->Write();

  hi_time_vtx->Write();
  hi_time_vtx_err->Write();
  hi_time_nstat->Write();

  hi_num_sho->Write();
  hi_num_hits->Write();
  hi_sho_pt->Write();
  hi_sho_p->Write();
  hi_sho_eta->Write();

  hFile->Write();
}

float 
GlobalMuonValidator::calculateDistance(const math::XYZVector& vect1, const math::XYZVector& vect2) {
  float dEta = vect1.eta() - vect2.eta();
  float dPhi = fabs(Geom::Phi<float>(vect1.phi()) - Geom::Phi<float>(vect2.phi()));
  float distance = sqrt(pow(dEta,2) + pow(dPhi,2) );

  return distance;
}

//
// return h1/h2 with recalculated errors
//
TH1F* GlobalMuonValidator::divideErr(TH1F* h1, TH1F* h2, TH1F* hout) {

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

//
// choose final trajectory
//
const Track* 
GlobalMuonValidator::chooseTrack(vector<const Track*> t, int muonHitsOption, int p1, int p2) const {

  const Track* result = 0;
  std::ofstream out("Muon_reco.txt",ios::app);
  
  if ( muonHitsOption == 0 ) {
    if (t[0]) result = t[0];
    return result;
  } else if ( muonHitsOption == 1 ) {
    if (t[1]) result = t[1];
    return result;
  } else if ( muonHitsOption == 2 ) {
    if (t[2]) result = t[2];
    return result;
  } else if ( muonHitsOption == 3 ) {
    if (t[3]) result = t[3];
    return result;
  } else if ( muonHitsOption == 4 ) {
    double prob0 = ( t[0] ) ? trackProbability(*t[0]) : 0.0;
    double prob1 = ( t[1] ) ? trackProbability(*t[1]) : 0.0;
    double prob2 = ( t[2] ) ? trackProbability(*t[2]) : 0.0;
    double prob3 = ( t[3] ) ? trackProbability(*t[3]) : 0.0; 
    
    
    if ( t[1] ) result = t[1];
    if ( (t[1] == 0) && t[3] ) result = t[3];
  
    if ( t[1] && t[3] && ( (prob1 - prob3) > 0.05 )  )  result = t[3];

    if ( t[0] && t[2] && fabs(prob2 - prob0) > 30. ) {
      result = t[0];
      return result;
    }

    if ( (t[1] == 0) && (t[3] == 0) && t[2] ) result = t[2];

    const Track* tmin = 0;
    double probmin = 0.0;
    if ( t[1] && t[3] ) {
      probmin = prob3; tmin = t[3];
      if ( prob1 < prob3 ) { probmin = prob1; tmin = t[1]; }
    }
    else if ( (t[3] == 0) && t[1] ) { 
      probmin = prob1; tmin = t[1]; 
    }
    else if ( (t[1] == 0) && t[3] ) {
      probmin = prob3; tmin = t[3]; 
    }

    if ( tmin && t[2] && ( (probmin - prob2) > 3.5 )  ) {
      result = t[2];
    }

  } else if ( muonHitsOption == 5 ) {

    double prob[4];
    double pt;
    int chosen=3;
    for (int i=0;i<4;i++) {
      prob[i] = (t[i]) ? trackProbability(*t[i]) : 0.0; 
      pt=0;
      if (t[i]) {
//        pt = t[i]->lastMeasurement().updatedState().globalMomentum().perp();
        pt = t[i]->pt();
        out << pt << "  " << trackProbability(*t[i]) << "  ";
      } else out << "0.0 0.0 ";
    }    
    out << endl;

    if (!t[3])
      if (t[2]) chosen=2; else
        if (t[1]) chosen=1; else
          if (t[0]) chosen=0;

    if ( t[0] && t[3] && ((prob[3]-prob[0]) > 4) ) chosen=0;
    if ( t[0] && t[1] && ((prob[1]-prob[0]) < p1) ) chosen=1;
    if ( t[2] && ((prob[chosen]-prob[2]) > p2) ) chosen=2;
    
    result=t[chosen];
  }

  return result;

}


//
// calculate the tail probability (-ln(P)) of a fit
//
double 
GlobalMuonValidator::trackProbability(const Track& track) const {

  int nDOF = (int)track.ndof();
  if ( nDOF > 0 && track.chi2()> 0) { 
    return -LnChiSquaredProbability(track.chi2(), nDOF);
  } else { 
    return 0.0;
  }
}


//
//
//
void GlobalMuonValidator::checkMuonHits(const reco::Track& muon, 
				       std::vector<int>& hits) const {

  float coneSize = 20.0;
  int dethits[4];
  for ( int i=0; i<4; i++ ) hits[i]=dethits[i]=0;

  // loop through all muon hits and calculate the maximum # of hits in each chamber
  for (trackingRecHit_iterator imrh = muon.recHitsBegin(); imrh != muon.recHitsEnd(); imrh++ ) {
        
    if (!(*imrh)->isValid()) continue;
  
    int station = 0;
    int detRecHits = 0;
      
    DetId id = (*imrh)->geographicalId();

    // Skip tracker hits
    if (id.det()!=DetId::Muon) continue;
      
    if ( id.subdetId() == MuonSubdetId::DT ) {

      DTChamberId did(id.rawId());
      DTLayerId lid(id.rawId());
      station = did.station();
//      cout << " Station " << station << " DT layer: " << lid.layer() << endl;

      // Get the DT-Segment which relies on this chamber
      DTRecHitCollection::range dRecHits = theDTRecHits->get(lid);
	
      for (DTRecHitCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
	double rhitDistance = fabs(ir->localPosition().x()-(**imrh).localPosition().x());
	if ( rhitDistance < coneSize ) detRecHits++;
//             cout << "       " << (ir)->localPosition() << "  " << (**imrh).localPosition()
//	     << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
      }
    }// end of if DT
    else if ( id.subdetId() == MuonSubdetId::CSC ) {
    
      CSCDetId did(id.rawId());
      // Get the CSC-Segment which relies on this chamber
      CSCRecHit2DCollection::range dRecHits = theCSCRecHits->get(did);      
    
      station = did.station();
//      cout << " Station " << station << " CSC layer: " << did.layer() << endl;

      for (CSCRecHit2DCollection::const_iterator ir = dRecHits.first; ir != dRecHits.second; ir++ ) {
	double rhitDistance = (ir->localPosition()-(**imrh).localPosition()).mag();
	if ( rhitDistance < coneSize ) detRecHits++;
//	cout << ir->localPosition() << "  " << (**imrh).localPosition()
//	     << " Distance: " << rhitDistance << " recHits: " << detRecHits << endl;
      }
    }
    else {
//      cout<<" Wrong Hit Type " << endl;
      continue;      
    }
      
    if ( (station > 0) && (station < 5) ) {
      if ( detRecHits > hits[station-1] ) hits[station-1] = detRecHits;
    }

  } // end of loop over muon rechits

}

//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMuonValidator);
