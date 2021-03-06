import FWCore.ParameterSet.Config as cms

hltResults = cms.EDAnalyzer("FourVectorHLT",
    plotAll = cms.untracked.bool(True),
    ptMax = cms.untracked.double(100.0),
    ptMin = cms.untracked.double(0.0),
    filters = cms.VPSet(),
    # data best guess
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD","","HLT")
)


