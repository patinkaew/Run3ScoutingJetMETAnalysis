import FWCore.ParameterSet.Config as cms

offlinePuppiJetTightLeptonVetoIdCriteria = cms.VPSet(
  cms.PSet(
    min_abs_eta = cms.double(0.0),
    max_abs_eta = cms.double(2.6),
    max_neutral_hadron_energy_fraction = cms.double(0.99),
    max_neutral_em_energy_fraction = cms.double(0.90),
    min_num_constituents = cms.int32(1),
    max_muon_energy_fraction = cms.double(0.80), # for LepVeto
    min_charged_hadron_energy_fraction = cms.double(0.01),
    min_charged_multiplicity = cms.int32(0),
    max_charged_em_energy_fraction = cms.double(0.80), # for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
  cms.PSet( 
    min_abs_eta = cms.double(2.6),
    max_abs_eta = cms.double(2.7),
    max_neutral_hadron_energy_fraction = cms.double(0.90),
    max_neutral_em_energy_fraction = cms.double(0.99),
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(0.80), # for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not used
    min_charged_multiplicity = cms.int32(-1), # not used
    max_charged_em_energy_fraction = cms.double(0.80), # for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
  cms.PSet( 
    min_abs_eta = cms.double(2.7),
    max_abs_eta = cms.double(3.0),
    max_neutral_hadron_energy_fraction = cms.double(0.99),
    max_neutral_em_energy_fraction = cms.double(99), # not used
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(99), # not used, for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not used
    min_charged_multiplicity = cms.int32(-1), # not used
    max_charged_em_energy_fraction = cms.double(99), # not used, for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
  cms.PSet( 
    min_abs_eta = cms.double(3.0),
    max_abs_eta = cms.double(5.0),
    max_neutral_hadron_energy_fraction = cms.double(99), # not used
    max_neutral_em_energy_fraction = cms.double(0.4),
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(99), # not used, for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not used
    min_charged_multiplicity = cms.int32(-1), # not used
    max_charged_em_energy_fraction = cms.double(99), # not used, for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
)

offlineCHSJetTightLeptonVetoIdCriteria = cms.VPSet(
  cms.PSet(
    min_abs_eta = cms.double(0.0),
    max_abs_eta = cms.double(2.6),
    max_neutral_hadron_energy_fraction = cms.double(0.99),
    max_neutral_em_energy_fraction = cms.double(0.90),
    min_num_constituents = cms.int32(1),
    max_muon_energy_fraction = cms.double(0.80), # for LepVeto
    min_charged_hadron_energy_fraction = cms.double(0.01),
    min_charged_multiplicity = cms.int32(0),
    max_charged_em_energy_fraction = cms.double(0.80), # for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
  cms.PSet( 
    min_abs_eta = cms.double(2.6),
    max_abs_eta = cms.double(2.7),
    max_neutral_hadron_energy_fraction = cms.double(0.90),
    max_neutral_em_energy_fraction = cms.double(0.99),
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(0.80), # for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not use
    min_charged_multiplicity = cms.int32(0),
    max_charged_em_energy_fraction = cms.double(0.80), # for LepVeto
    min_neutral_multiplicity = cms.int32(-1), # not used
  ),
  cms.PSet( 
    min_abs_eta = cms.double(2.7),
    max_abs_eta = cms.double(3.0),
    max_neutral_hadron_energy_fraction = cms.double(0.99),
    max_neutral_em_energy_fraction = cms.double(0.99),
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(99), # not used, for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not used
    min_charged_multiplicity = cms.int32(-1), # not used
    max_charged_em_energy_fraction = cms.double(99), # not used, for LepVeto
    min_neutral_multiplicity = cms.int32(2), # >1
  ),
  cms.PSet( 
    min_abs_eta = cms.double(3.0),
    max_abs_eta = cms.double(5.0),
    max_neutral_hadron_energy_fraction = cms.double(99), # not used
    max_neutral_em_energy_fraction = cms.double(0.4),
    min_num_constituents = cms.int32(-1), # not used
    max_muon_energy_fraction = cms.double(99), # not used, for LepVeto
    min_charged_hadron_energy_fraction = cms.double(-1), # not used
    min_charged_multiplicity = cms.int32(-1), # not used
    max_charged_em_energy_fraction = cms.double(99), # not used, for LepVeto
    min_neutral_multiplicity = cms.int32(11), # >10
  ),
)

offlinePuppiJetTightLeptonVetoId = cms.EDProducer("PATJetIdProducer",
  src = cms.InputTag("slimmedJetsPuppi", "", "RECO"),
  criteria = offlinePuppiJetTightLeptonVetoIdCriteria,
)

offlineCHSJetTightLeptonVetoId = cms.EDProducer("PATJetIdProducer",
  src = cms.InputTag("slimmedJets", "", "RECO"),
  criteria = offlineCHSJetTightLeptonVetoIdCriteria,
)

offlinePFJetTightLeptonVetoId = cms.EDProducer("RecoPFJetIdProducer",
  src = cms.InputTag("offlinePFJet"),
  criteria = offlineCHSJetTightLeptonVetoIdCriteria,
)

scoutingPFJetTightLeptonVetoId = cms.EDProducer("ScoutingPFJetIdProducer",
  src = cms.InputTag("hltScoutingPFPacker", "", "HLT"),
  criteria = offlineCHSJetTightLeptonVetoIdCriteria,
)
