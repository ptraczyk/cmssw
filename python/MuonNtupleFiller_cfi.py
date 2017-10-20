import FWCore.ParameterSet.Config as cms

muonNtupleFiller = cms.EDAnalyzer("MuonNtupleFiller",

    mctruthMatching = cms.bool(True),

    Muons = cms.untracked.InputTag("muons"),
    TKtracks = cms.untracked.InputTag("generalTracks"),
    Timing = cms.untracked.InputTag("muons"),

    # cosmic ID back-to-back angle cut 
    angleCut = cms.double(0.02),
    PtCut = cms.double(200.0),

    open = cms.string('recreate'),
    out = cms.string('muonNtuple.root'),
    debug= cms.bool(False)
)
