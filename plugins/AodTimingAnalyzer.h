#ifndef RecoMuon_MuonIdentification_AODTimingAnalyzer_H
#define RecoMuon_MuonIdentification_AODTimingAnalyzer_H

/** \class AODTimingAnalyzer
 *  Analyzer of the timing information in the reco::Muon object
 *
 *  $Date: 2011/04/06 09:56:30 $
 *  $Revision: 1.4 $
 *  \author P. Traczyk    CERN
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "RecoMuon/TrackingTools/interface/MuonSegmentMatcher.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TROOT.h>
#include <TSystem.h>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class InputTag;
}

class TFile;
class TH1F;
class TH2F;

using namespace std;
using namespace edm;
using namespace reco;

class AODTimingAnalyzer : public edm::EDAnalyzer {
//class AODTimingAnalyzer : public edm::stream::EDProducer<> {
public: 

  explicit AODTimingAnalyzer(const edm::ParameterSet&);
  ~AODTimingAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  double iMass(reco::TrackRef imuon, reco::TrackRef iimuon);
  bool dumpMuonId(const pat::Muon& muon, const reco::Vertex& vtx, const bool debug);
  void dumpTrack(reco::TrackRef track);

  // ----------member data ---------------------------

  edm::ConsumesCollector *iC;
  edm::InputTag TKtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag VtxTags_; 
  edm::InputTag TimeTags_; 
  edm::InputTag SIMtrackTags_; 

  string out, open;
  bool theDebug;
  bool doSim;
  string theIdCut;
  bool theCollVeto;
  bool theKeepBX;
  int theBX;
  bool theVetoCosmics;
  bool theOnlyCosmics;
  double theAngleCut;
  double  theMinEta, theMaxEta, thePtCut, theMinPtres, theMaxPtres, theScale;
  int theDtCut, theCscCut;
  int theNBins;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;

  Handle<pat::MuonCollection> MuCollection;
  Handle<pat::MuonCollection> MuCollectionT;
  Handle<reco::TrackCollection> TKTrackCollection;
  Handle<reco::TrackCollection> STATrackCollection;
  Handle<reco::TrackCollection> GLBTrackCollection;
  Handle<reco::TrackCollection> PMRTrackCollection;
  Handle<reco::TrackCollection> GMRTrackCollection;
  Handle<reco::TrackCollection> FMSTrackCollection;
  Handle<reco::TrackCollection> SLOTrackCollection;
  Handle<edm::SimTrackContainer> SIMTrackCollection;

  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;

  TH1F* hi_gen_pt;
  TH1F* hi_gen_eta;
  TH1F* hi_gen_phi;

  TH1F* hi_id_rpccut_sta;
  TH1F* hi_id_rpccut_glb;
  TH1F* hi_id_csccut_sta;
  TH1F* hi_id_csccut_glb;
  TH2F* hi_id_dtcut_sta;
  TH2F* hi_id_dtcut_glb;
  TH2F* hi_id_cmbcut_sta;
  TH2F* hi_id_cmbcut_glb;

  TH1F* hi_id_trklay;
  TH1F* hi_id_trkhit;
  TH1F* hi_id_statio;
  TH1F* hi_id_dxy;
  TH1F* hi_id_dz;

  TH1F* hi_glb_mass_ss  ;
  TH1F* hi_glb_mass_os  ;
  TH1F* hi_sta_mass_ss  ;
  TH1F* hi_sta_mass_os  ;

  TH1F* hi_glb_angle;
  TH1F* hi_trk_angle;
  TH1F* hi_glb_angle_w;
  TH1F* hi_trk_angle_w;
  TH2F* hi_dttime_vtx_tb_angle;

  TH1F* hi_sta_pt  ;
  TH1F* hi_sta_pt_cut;  
  TH1F* hi_sta_ptres  ;
  TH1F* hi_sta_ptg  ;
  TH1F* hi_sta_ptt  ;
  TH1F* hi_sta_ptres_tb  ;
  TH1F* hi_tk_pt   ;
  TH1F* hi_glb_pt  ;
  TH1F* hi_glb_pt_cut  ;
  TH1F* hi_glb_ptg  ;
  TH1F* hi_glb_ptt  ;
  TH1F* hi_glb_ptres  ;
  TH1F* hi_glb_ptresh ;
  TH1F* hi_glb_ptres_t  ;
  TH1F* hi_glb_ptresh_t ;
  TH1F* hi_glb_ptres_b  ;
  TH1F* hi_glb_ptresh_b ;
  TH1F* hi_glb_ptres_tb  ;
  TH1F* hi_glb_d0  ;
  TH1F* hi_sta_phi ;
  TH1F* hi_tk_phi  ;
  TH1F* hi_glb_phi ;
  TH1F* hi_sta_nhits ;
  TH1F* hi_tk_nhits  ;
  TH1F* hi_glb_nhits ;
  TH1F* hi_sta_nvhits ;
  TH1F* hi_tk_nvhits  ;
  TH1F* hi_glb_nvhits ;
  TH1F* hi_sta_chi2 ;
  TH1F* hi_tk_chi2  ;
  TH1F* hi_glb_chi2 ;
  TH1F* hi_tk_eta  ;
  TH1F* hi_sta_eta  ;
  TH1F* hi_glb_eta  ;

  TH1F* hi_mutime_ndof;
  TH1F* hi_mutime_vtx;
  TH1F* hi_mutime_vtx_err;

  TH1F* hi_dtcsc_vtx;
  TH1F* hi_dtcsc_vtx_t;
  TH1F* hi_dtcsc_vtx_b;

  TH2F* hi_dtrpc_vtx;
  TH2F* hi_cscrpc_vtx;
  TH2F* hi_cmbrpc_vtx;
  TH2F* hi_dtrpc3_vtx;
  TH2F* hi_cscrpc3_vtx;
  TH2F* hi_cmbrpc3_vtx;
  TH2F* hi_dtrpc3_vtxw;
  TH2F* hi_cscrpc3_vtxw;
  TH2F* hi_cmbrpc3_vtxw;

  TH1F* hi_trpc;
  TH1F* hi_trpc3;
  TH1F* hi_trpcerr;
  TH1F* hi_nrpc;
  TH2F* hi_nrpc_trpc;
  TH2F* hi_trpc_eta;
  TH2F* hi_trpc_phi;

};
#endif
