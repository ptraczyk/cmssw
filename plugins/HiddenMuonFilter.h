#ifndef RecoMuon_MuonIdentification_HiddenMuonFilter_H
#define RecoMuon_MuonIdentification_HiddenMuonFilter_H

/** \class HiddenMuonFilter
 *  Filter out Hidden muon events
 *
 *  $Date: 2010/12/15 11:01:59 $
 *  $Revision: 1.1 $
 *  \author P. Traczyk    CERN
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


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

class HiddenMuonFilter : public edm::EDFilter {
public:
  explicit HiddenMuonFilter(const edm::ParameterSet&);
  ~HiddenMuonFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::MuonCollection> MuCollectionC;
  Handle<reco::VertexCollection> vtxCollectionC;
};
#endif
