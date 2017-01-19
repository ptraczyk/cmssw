import FWCore.ParameterSet.Config as cms

aodTimingAnalyzer = cms.EDAnalyzer("AODTimingAnalyzer",

    mctruthMatching = cms.bool(True),

# Event input tags
    Muons = cms.untracked.InputTag("muons"),
    TKtracks = cms.untracked.InputTag("generalTracks"),
    Timing = cms.untracked.InputTag("muons"),
#    PrimaryVertex = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),

# Event-level cuts
    collisionVeto = cms.bool(False),
    keepOnlyBX = cms.bool(False),
    generatedBX = cms.int32(0),
    vetoCosmics = cms.bool(False),
    onlyCosmics = cms.bool(False),
    # cosmic ID back-to-back angle cut 
    angleCut = cms.double(0.02),

# Muon-level cuts
    requireId = cms.string(""),
    PtCut = cms.double(5.0),
    etaMin = cms.double(0.0),
    etaMax = cms.double(2.5),
    DTcut  = cms.int32(6),
    CSCcut = cms.int32(4),

# Output plot parameters
    PtresMax = cms.double(400.0),
    PtresMin = cms.double(0.0),
    PlotScale = cms.double(1.0),
    nbins = cms.int32(100),
    open = cms.string('recreate'),
    out = cms.string('aodTimingAnalyzer.root'),
    debug= cms.bool(False)
)
