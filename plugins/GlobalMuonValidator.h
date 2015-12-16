#ifndef RecoMuon_GlobalMuonProducer_GlobalMuonValidator_H
#define RecoMuon_GlobalMuonProducer_GlobalMuonValidator_H

/** \class GlobalMuonValidator
 *  Analyzer of the StandAlone and Global muon tracks
 *
 *  $Date: 2011/09/26 09:46:02 $
 *  $Revision: 1.2 $
 *  \author A. Everett     Purdue University
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
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

class GlobalMuonValidator : public edm::EDAnalyzer {
public:
  explicit GlobalMuonValidator(const edm::ParameterSet&);
  ~GlobalMuonValidator();
  
  typedef std::pair< TrackRef, SimTrackRef> CandToSim;
  typedef std::pair< TrackRef, SimTrackRef> CandStaSim;
  typedef std::pair< TrackRef, SimTrackRef> CandMuonSim;
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual float calculateDistance(const math::XYZVector&, const math::XYZVector&);
  virtual TH1F* divideErr(TH1F*, TH1F*, TH1F*);
  const Track*  chooseTrack(vector<const Track*> t, int muonHitsOption, int p1, int p2) const;
  double trackProbability(const Track& track) const;
  void checkMuonHits(const reco::Track& muon, std::vector<int>& hits) const;


  // ----------member data ---------------------------

  edm::InputTag MuonTags_; 

  string out, open;
  double  theMinEta, theMaxEta, theMinPt, thePtCut, theMinPtres, theMaxPtres, theInvPt;
  int theNBins;

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::TrackCollection> TKTrackCollection;
  
  MuonServiceProxy* theService;

  edm::InputTag theDTRecHitLabel;
  edm::InputTag theCSCRecHitLabel;
    
  // caches that should get filled once per event
  Handle<DTRecHitCollection>    theDTRecHits;
  Handle<CSCRecHit2DCollection> theCSCRecHits;

  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;

  TH1F* hi_sta_pt  ;
  TH1F* hi_tk_pt  ;
  TH1F* hi_glb_pt  ;
  TH1F* hi_glb_diff  ;
  TH1F* hi_glb1_pt  ;
  TH1F* hi_glb2_pt  ;
  TH1F* hi_glb3_pt  ;
  TH1F* hi_glb4_pt  ;
  TH1F* hi_glb5_pt  ;
  TH1F* hi_glb6_pt  ;
  TH1F* hi_glbu4_pt  ;
  TH1F* hi_glbu5_pt  ;
  TH1F* hi_glb1_prob;
  TH1F* hi_glb2_prob;
  TH1F* hi_glb3_prob;
  TH1F* hi_glb4_prob;
  TH1F* hi_glb5_prob;
  TH1F* hi_glb6_prob;

  TH1F* hi_sta_pt_0;
  TH1F* hi_glb_pt_0;
  TH1F* hi_sta_pt_1;
  TH1F* hi_glb_pt_1;
  TH1F* hi_sta_pt_2;
  TH1F* hi_glb_pt_2;
  TH1F* hi_sta_pt_3;
  TH1F* hi_glb_pt_3;
  TH1F* hi_sta_pt_4;
  TH1F* hi_glb_pt_4;
  
  TH1F* hi_tune[10][10];

  TH1F* hi_time_vtx;
  TH1F* hi_time_nstat;
  TH1F* hi_time_vtx_err;

  TH1F* hi_tk_eta  ;
  TH1F* hi_sta_eta  ;
  TH1F* hi_glb_eta  ;

  TH1F* hi_num_sho;
  TH1F* hi_num_hits;
  TH2F* hi_sho_pt;
  TH2F* hi_sho_p;
  TH2F* hi_sho_eta;

};
#endif
