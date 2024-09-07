import FWCore.ParameterSet.Config as cms

scoutingPFMET = cms.EDProducer("Run3ScoutingPFMETProducer",
  pt = cms.InputTag("hltScoutingPFPacker", "pfMetPt", "HLT"),
  phi = cms.InputTag("hltScoutingPFPacker", "pfMetPhi", "HLT")
)
