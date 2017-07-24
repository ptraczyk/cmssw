#ifndef RecoMuon_MuonIdentification_AodNtupleFiller_H
#define RecoMuon_MuonIdentification_AodNtupleFiller_H

/** \class AodNtupleFiller
 *  Dump muon information into a simple ntuple
 *
 *  $Date: 2011/04/06 09:56:30 $
 *  $Revision: 1.4 $
 *  \author P. Traczyk   
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "RecoMuon/TrackingTools/interface/MuonSegmentMatcher.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include <TROOT.h>
#include <TSystem.h>
#include "Riostream.h"    
#include "TObject.h"
#include <TTree.h>

namespace edm {
  class ParameterSet;
  //  class Event;
  class EventSetup;
  class InputTag;
}

class TFile;
class TTree;

using namespace std;
using namespace edm;
using namespace reco;

class AodNtupleFiller : public edm::EDAnalyzer {
//class AodNtupleFiller : public edm::stream::EDProducer<> {
public: 

  explicit AodNtupleFiller(const edm::ParameterSet&);
  ~AodNtupleFiller();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  edm::ConsumesCollector *iC;
  edm::InputTag TKtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag TimeTags_; 
  edm::InputTag SIMtrackTags_; 

  string out, open;
  bool theDebug;
  bool doSim;
  double theAngleCut;
  double thePtCut;

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

  Handle<reco::MuonTimeExtraMap> timeMap1;
  Handle<reco::MuonTimeExtraMap> timeMap2;
  Handle<reco::MuonTimeExtraMap> timeMap3;
  
  // ROOT Pointers
  TFile* hFile;
  TTree* t;

  // ***** Tree structure *******
// generator info (if available)
  bool hasSim;
  int genCharge;
  double genPt, genPhi, genEta;
  int genBX;

// event info
  unsigned int event_run;
  unsigned int event_lumi;
  unsigned int event_event;
  bool isCosmic;
  bool isCollision;

// muon ID
  bool isPF;
  bool isSTA;
  bool isGLB;
  bool isLoose;
  bool isTight;

// muon kinematics and track fit parameters
  int charge;
  double pt, phi, eta;
  double dPt;
  double dz;
  double dxy;
  
  int nhits[4];
  int ssize[4];
  
// muon timing
  int muNdof;
  double muTime;
  double muTimeErr;
  int dtNdof;
  double dtTime;
  int cscNdof;
  double cscTime;
  int rpcNdof;
  double rpcTime;
  double rpcTimeErr;

};
#endif
