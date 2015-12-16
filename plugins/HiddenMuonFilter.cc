// -*- C++ -*-
//
// Package:    HiddenMuonFilter
// Class:      HiddenMuonFilter
// 
/**\class HiddenMuonFilter HiddenMuonFilter.cc 

 Description: An EDFilter selecting events with "Hidden" muons

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Piotr Traczyk
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: HiddenMuonFilter.cc,v 1.1 2010/12/15 11:01:59 ptraczyk Exp $
//
//

#include "HiddenMuonFilter.h"

// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDFilter.h"

//#include "FWCore/Framework/interface/Event.h"
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

#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"


//
// constructors and destructor
//
HiddenMuonFilter::HiddenMuonFilter(const edm::ParameterSet& iConfig) {
}

HiddenMuonFilter::~HiddenMuonFilter() {
}

//
// member functions
//

// ------------ method called to for each event  ------------
bool
HiddenMuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::MuonCollection;

  iEvent.getByLabel("muons",MuCollection);
  const reco::MuonCollection muon = *(MuCollection.product());

  iEvent.getByLabel("muonsFromCosmics",MuCollectionC);
  const reco::MuonCollection muonC = *(MuCollectionC.product());

  iEvent.getByLabel("offlinePrimaryVertices",vtxCollectionC);
  const reco::VertexCollection vtxC = *(vtxCollectionC.product());
  
  cout << " First vertex has ndof = " << vtxC[0].ndof() << endl;
  
  if (vtxC[0].ndof()<4) return true;

  if (muonC.size()) {
//    cout << " Cosmic muons: " <<;
//    if (!muon.size()) return(true);
  }

  return(false);
}


void 
HiddenMuonFilter::beginJob() {
}

void 
HiddenMuonFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiddenMuonFilter);
