#ifndef RecoMuon_MuonIdentification_MuonTimingAnalyzer_H
#define RecoMuon_MuonIdentification_MuonTimingAnalyzer_H

/** \class MuonTimingAnalyzer
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

#include <TROOT.h>
#include <TSystem.h>

namespace edm {
  class ParameterSet;
  //  class Event;
  class EventSetup;
  class InputTag;
}

class TFile;
class TH1F;
class TH2F;
//class TrackRef;
//class SimTrackRef;
//class MuonRef;
class MuonServiceProxy;

using namespace std;
using namespace edm;
using namespace reco;

class MuonTimingAnalyzer : public edm::EDAnalyzer {
//class MuonTimingAnalyzer : public edm::stream::EDProducer<> {
public: 

  explicit MuonTimingAnalyzer(const edm::ParameterSet&);
  ~MuonTimingAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual float calculateDistance(const math::XYZVector&, const math::XYZVector&);
  virtual TH1F* divideErr(TH1F*, TH1F*, TH1F*);
  double iMass(reco::TrackRef imuon, reco::TrackRef iimuon);
  bool dumpMuonId(const reco::Muon& muon, const reco::Vertex& vtx, const bool debug);
  void dumpTrack(reco::TrackRef track);

  // ----------member data ---------------------------

  MuonSegmentMatcher *theMatcher;
  edm::ConsumesCollector *iC;
  edm::InputTag TKtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag TimeTags_; 
  edm::InputTag SIMtrackTags_; 

  string out, open;
  bool theDebug;
  bool doSim;
  bool theOnlyGlb;
  bool theCollVeto;
  bool theVetoCosmics;
  bool theOnlyCosmics;
  double theAngleCut;
  double  theMinEta, theMaxEta, thePtCut, theMinPtres, theMaxPtres, theScale;
  int theDtCut, theCscCut;
  int theNBins;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapCmbToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapDTToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> timeMapCSCToken_;
  edm::EDGetTokenT<GenParticleCollection> genParticleToken_;

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::MuonCollection> MuCollectionT;
  Handle<reco::TrackCollection> TKTrackCollection;
  Handle<reco::TrackCollection> STATrackCollection;
  Handle<reco::TrackCollection> GLBTrackCollection;
  Handle<reco::TrackCollection> PMRTrackCollection;
  Handle<reco::TrackCollection> GMRTrackCollection;
  Handle<reco::TrackCollection> FMSTrackCollection;
  Handle<reco::TrackCollection> SLOTrackCollection;
  Handle<edm::SimTrackContainer> SIMTrackCollection;

  Handle<reco::MuonTimeExtraMap> timeMap1;
  Handle<reco::MuonTimeExtraMap> timeMap2;
  Handle<reco::MuonTimeExtraMap> timeMap3;
  
  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;

  TH1F* hi_id_rpccut_sta;
  TH1F* hi_id_rpccut_glb;
  TH1F* hi_id_csccut_sta;
  TH1F* hi_id_csccut_glb;
  TH2F* hi_id_dtcut_sta;
  TH2F* hi_id_dtcut_glb;

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

  TH1F* hi_mutime_vtx;
  TH1F* hi_mutime_vtx_err;

  TH1F* hi_dtcsc_vtx;
  TH1F* hi_dtcsc_vtx_t;
  TH1F* hi_dtcsc_vtx_b;

  TH1F* hi_trpc;
  TH1F* hi_nrpc;
  TH2F* hi_trpc_eta;
  TH2F* hi_trpc_phi;

  TH1F* hi_cmbtime_ibt;
  TH2F* hi_cmbtime_ibt_pt;
  TH1F* hi_cmbtime_ibt_err;
  TH1F* hi_cmbtime_ibt_pull;
  TH1F* hi_cmbtime_fib;
  TH1F* hi_cmbtime_fib_err;
  TH1F* hi_cmbtime_fib_pull;
  TH1F* hi_cmbtime_vtx;
  TH1F* hi_cmbtime_vtx_err;
  TH1F* hi_cmbtime_vtx_pull;
  TH1F* hi_cmbtime_vtxr;
  TH1F* hi_cmbtime_vtxr_err;
  TH1F* hi_cmbtime_vtxr_pull;
  TH1F* hi_cmbtime_ndof;

  TH1F* hi_dttime_ibt;
  TH2F* hi_dttime_ibt_pt;
  TH1F* hi_dttime_ibt_err;
  TH1F* hi_dttime_ibt_pull;
  TH1F* hi_dttime_fib;
  TH1F* hi_dttime_fib_t;
  TH1F* hi_dttime_fib_b;
  TH2F* hi_dttime_fibp_t;
  TH2F* hi_dttime_fibp_b;
  TH1F* hi_dttime_fib_err;
  TH1F* hi_dttime_fib_pull;
  TH1F* hi_dttime_vtx;
  TH2F* hi_dttime_vtxn;
  TH1F* hi_dttime_vtx_w;
  TH2F* hi_dttime_vtx_pt;
  TH2F* hi_dttime_vtx_phi;
  TH2F* hi_dttime_vtx_eta;
  TH2F* hi_dttime_etaphi;
  TH2F* hi_dttime_eeta_lo;
  TH2F* hi_dttime_eeta_hi;
  TH2F* hi_dttime_vtx_etat;
  TH2F* hi_dttime_vtx_etab;
  TH1F* hi_dttime_vtx_t;
  TH1F* hi_dttime_vtx_b;
  TH1F* hi_dttime_vtx_to;
  TH1F* hi_dttime_vtx_bo;
  TH1F* hi_dttime_vtx_tb;
  TH2F* hi_dttime_vtx_tb2;
  TH2F* hi_dttime_vtxp_t;
  TH2F* hi_dttime_vtxp_b;
  TH2F* hi_dttime_vtxp_tb;
  TH2F* hi_dttime_vtxpt_tb;
  TH1F* hi_dttime_vtx_err;
  TH1F* hi_dttime_vtx_pull;
  TH1F* hi_dttime_vtxr;
  TH1F* hi_dttime_vtxr_err;
  TH1F* hi_dttime_vtxr_pull;
  TH1F* hi_dttime_errdiff;
  TH1F* hi_dttime_errdiff_t;
  TH1F* hi_dttime_errdiff_b;
  TH1F* hi_dttime_ndof;

  TH1F* hi_csctime_ibt;
  TH2F* hi_csctime_ibt_pt;
  TH1F* hi_csctime_ibt_err;
  TH1F* hi_csctime_ibt_pull;
  TH1F* hi_csctime_fib;
  TH1F* hi_csctime_fib_t;
  TH1F* hi_csctime_fib_b;
  TH1F* hi_csctime_fib_err;
  TH1F* hi_csctime_fib_pull;
  TH1F* hi_csctime_vtx;
  TH2F* hi_csctime_vtx_pt;
  TH1F* hi_csctime_vtx_t;
  TH1F* hi_csctime_vtx_b;
  TH2F* hi_csctime_vtx_eta;
  TH2F* hi_csctime_vtx_phi;
  TH2F* hi_csctime_eeta_lo;
  TH2F* hi_csctime_eeta_hi;
  TH2F* hi_csctime_vtx_etat;
  TH2F* hi_csctime_vtx_etab;
  TH1F* hi_csctime_vtx_err;
  TH1F* hi_csctime_vtx_pull;
  TH1F* hi_csctime_vtxr;
  TH1F* hi_csctime_vtxr_err;
  TH1F* hi_csctime_vtxr_pull;
  TH1F* hi_csctime_ndof;


};
#endif
