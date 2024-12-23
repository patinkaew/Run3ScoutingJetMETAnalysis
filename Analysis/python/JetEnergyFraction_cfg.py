import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:ac2fde62-8866-4e50-8f36-77d604bfd4e1.root", # ScoutingPFRun3
        "file:e4029e0e-090a-4096-8295-1e6e4faee6f9.root", # ScoutingPFMonitor 2024F
        #"file:bfafe96c-4a6c-43a2-a711-762d7fec8530.root", # ScoutingPFMonitor 2023D
    ),
    #skipEvents = cms.untracked.uint32(9),
)

#import FWCore.PythonUtilities.LumiList as LumiList
#json_filepath = "/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_384380_Golden.json"
#process.source.lumisToProcess = LumiList.LumiList(filename = json_filepath).getVLuminosityBlockRange()

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("jet_energy_fraction.root")
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "140X_dataRun3_Prompt_v4", "")

# L1 triggger
process.load("L1Trigger.Configuration.L1TRawToDigi_cff")

from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
#process.gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
process.gtStage2DigisScouting = gtStage2Digis.clone()

from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits
#process.l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting")
process.l1bitsScouting = l1bits.clone()

process.l1bit_sequence = cms.Sequence(process.L1TRawToDigi + cms.Sequence(cms.Task(process.gtStage2DigisScouting, process.l1bitsScouting)))

process.load("Run3ScoutingJetMETAnalysis.Analysis.JetIdProducer_cff")
#process.load("Run3ScoutingJetMETAnalysis.Analysis.Run3ScoutingPFMETProducer_cff")

trigger_task_matrix = {
  "2023_early" : { # until 2023C-v2, Run367620
    "triggers" : cms.VPSet(
      cms.PSet(
          name=cms.string("DST_PFScouting_JetHT"), # this expr is taken from confdb, l1 seeds for path: DST_Run3_PFScoutingPixelTracking_v* (only hadronic, drop muon)
          expr=cms.vstring("L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_HTT400er", "L1_HTT450er", "L1_ETT2000", 
            "L1_SingleJet180", "L1_SingleJet200", "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5")
      ),
      #cms.PSet(expr=cms.vstring("L1_HTT200er")),
      #cms.PSet(expr=cms.vstring("L1_HTT255er")),
      #cms.PSet(expr=cms.vstring("L1_HTT280er")),
      #cms.PSet(expr=cms.vstring("L1_HTT320er")),
      #cms.PSet(expr=cms.vstring("L1_HTT360er")),
      #cms.PSet(expr=cms.vstring("L1_HTT400er")),
      #cms.PSet(expr=cms.vstring("L1_HTT450er")),
      #cms.PSet(expr=cms.vstring("L1_ETT2000")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet180")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet200")),
    ),
  },
  "2023_late" : { # from 2023C-v3, Run367621
    "triggers" : cms.VPSet(
      cms.PSet(
        name=cms.string(""), # the actual path is DST_Run3_JetHT_PFScoutingPixelTracking. One l1 seed: L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 is added
        expr=cms.vstring("DST_Run3_JetHT_PFScoutingPixelTracking")
      ),
      #cms.PSet(expr=cms.vstring("L1_HTT200er")),
      #cms.PSet(expr=cms.vstring("L1_HTT255er")),
      #cms.PSet(expr=cms.vstring("L1_HTT280er")),
      #cms.PSet(expr=cms.vstring("L1_HTT320er")),
      #cms.PSet(expr=cms.vstring("L1_HTT360er")),
      #cms.PSet(expr=cms.vstring("L1_HTT400er")),
      #cms.PSet(expr=cms.vstring("L1_HTT450er")),
      #cms.PSet(expr=cms.vstring("L1_ETT2000")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet180")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet200")),
    ),
  },
  "2024" : {
    "triggers" : cms.VPSet(
      cms.PSet(expr=cms.vstring("DST_PFScouting_JetHT")),
      #cms.PSet(expr=cms.vstring("L1_HTT200er")),
      #cms.PSet(expr=cms.vstring("L1_HTT255er")),
      #cms.PSet(expr=cms.vstring("L1_HTT280er")),
      #cms.PSet(expr=cms.vstring("L1_HTT320er")),
      #cms.PSet(expr=cms.vstring("L1_HTT360er")),
      #cms.PSet(expr=cms.vstring("L1_HTT400er")),
      #cms.PSet(expr=cms.vstring("L1_HTT450er")),
      #cms.PSet(expr=cms.vstring("L1_ETT2000")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet180")),
      #cms.PSet(expr=cms.vstring("L1_SingleJet200")),
      cms.PSet(expr=cms.vstring("L1_ZeroBias")),
      cms.PSet(expr=cms.vstring("DST_PFScouting_ZeroBias")),
    ),
  }
}

# select task to run
trigger_task = trigger_task_matrix["2024"]

process.ScoutingPFJet = cms.EDAnalyzer("ScoutingPFJetEnergyFractionAnalyzer",
  src = cms.InputTag("scoutingPFJetTightLeptonVetoId"),
  L1TriggerResults = cms.InputTag("l1bitsScouting"),
  HLTTriggerResults = cms.InputTag("TriggerResults", "",  "HLT"),
  num_jets_to_fill = cms.int32(-1),
  energy_fractions = cms.VPSet( #neutral hadron contains HFHadron and neutral EM contains HF EM
    cms.PSet(name=cms.string("Charged_Hadron"), label=cms.string("Charged Hadron Energy Fraction"), func=cms.string("chargedHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Neutral_Hadron"), label=cms.string("Neutral Hadron Energy Fraction"), func=cms.string("neutralHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Neutral_Hadron_noHF"), label=cms.string("Neutral Hadron Energy Fraction (no HF)"), func=cms.string("(neutralHadronEnergy()-HFHadronEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Charged_EM"), label=cms.string("Charged EM Energy Fraction"), func=cms.string("electronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Neutral_EM"), label=cms.string("Neutral EM Energy Fraction"), func=cms.string("(photonEnergy()+HFEMEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Electron"), label=cms.string("Electron Energy Fraction"), func=cms.string("electronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Photon"), label=cms.string("Photon Energy Fraction"), func=cms.string("photonEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Muon"), label=cms.string("Muon Energy Fraction"), func=cms.string("muonEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("HF_Hadron"), label=cms.string("HF Hadron Energy Fraction"), func=cms.string("HFHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("HF_EM"), label=cms.string("HF EM Energy Fraction"), func=cms.string("HFEMEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Hadron"), label=cms.string("Hadron Energy Fraction"), func=cms.string("(chargedHadronEnergy()+neutralHadronEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Hadron_noHF"), label=cms.string("Hadron Energy Fraction (noHF)"), func=cms.string("(chargedHadronEnergy()+neutralHadronEnergy()-HFHadronEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
    cms.PSet(name=cms.string("Lepton"), label=cms.string("Lepton Energy Fraction"), func=cms.string("(electronEnergy()+muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy()+HFEMEnergy())")),
  ),
  **trigger_task
)
process.path0 = cms.Path(process.l1bit_sequence + process.scoutingPFJetTightLeptonVetoId + process.ScoutingPFJet)

process.OfflinePuppiJet = cms.EDAnalyzer("PATJetEnergyFractionAnalyzer",
  src = cms.InputTag("offlinePuppiJetTightLeptonVetoId"),
  L1TriggerResults = cms.InputTag("l1bitsScouting"),
  HLTTriggerResults = cms.InputTag("TriggerResults", "",  "HLT"),
  cut = cms.string("pt()>20"),
  jet_pt_func = cms.string("correctedP4('Uncorrected').Pt()"),
  jet_eta_func = cms.string("correctedP4('Uncorrected').Eta()"),
  jet_phi_func = cms.string("correctedP4('Uncorrected').Phi()"),
  num_jets_to_fill = cms.int32(-1),
  energy_fractions = cms.VPSet( # undo corrections
    cms.PSet(name=cms.string("Charged_Hadron"), label=cms.string("Charged Hadron Energy Fraction"), func=cms.string("correctedJet('Uncorrected').chargedHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_Hadron"), label=cms.string("Neutral Hadron Energy Fraction"), func=cms.string("correctedJet('Uncorrected').neutralHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_Hadron_noHF"), label=cms.string("Neutral Hadron Energy Fraction (no HF)"), func=cms.string("correctedJet('Uncorrected').neutralHadronEnergyFraction()-correctedJet('Uncorrected').HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Charged_EM"), label=cms.string("Charged EM Energy Fraction"), func=cms.string("correctedJet('Uncorrected').chargedEmEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_EM"), label=cms.string("Neutral EM Energy Fraction"), func=cms.string("correctedJet('Uncorrected').neutralEmEnergyFraction()")),
    cms.PSet(name=cms.string("Electron"), label=cms.string("Electron Energy Fraction"), func=cms.string("correctedJet('Uncorrected').electronEnergyFraction()")),
    cms.PSet(name=cms.string("Photon"), label=cms.string("Photon Energy Fraction"), func=cms.string("correctedJet('Uncorrected').photonEnergyFraction()")),
    cms.PSet(name=cms.string("Muon"), label=cms.string("Muon Energy Fraction"), func=cms.string("correctedJet('Uncorrected').muonEnergyFraction()")),
    cms.PSet(name=cms.string("HF_Hadron"), label=cms.string("HF Hadron Energy Fraction"), func=cms.string("correctedJet('Uncorrected').HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("HF_EM"), label=cms.string("HF EM Energy Fraction"), func=cms.string("correctedJet('Uncorrected').HFEMEnergyFraction()")),
    cms.PSet(name=cms.string("Hadron"), label=cms.string("Hadron Energy Fraction"), func=cms.string("correctedJet('Uncorrected').chargedHadronEnergyFraction()+correctedJet('Uncorrected').neutralHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Hadron_noHF"), label=cms.string("Hadron Energy Fraction (noHF)"), func=cms.string("correctedJet('Uncorrected').chargedHadronEnergyFraction()+correctedJet('Uncorrected').neutralHadronEnergyFraction()-correctedJet('Uncorrected').HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Lepton"), label=cms.string("Lepton Energy Fraction"), func=cms.string("correctedJet('Uncorrected').photonEnergyFraction()+correctedJet('Uncorrected').muonEnergyFraction()")),
  ),
  **trigger_task
)
process.path1 = cms.Path(process.l1bit_sequence + process.offlinePuppiJetTightLeptonVetoId + process.OfflinePuppiJet)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.offlinePFJet = ak4PFJets.clone(
  src = cms.InputTag("packedPFCandidates", "", "RECO"),
  jetPtMin = 20,
)

process.OfflinePFJet = cms.EDAnalyzer("RecoPFJetEnergyFractionAnalyzer",
  src = cms.InputTag("offlinePFJetTightLeptonVetoId"),
  L1TriggerResults = cms.InputTag("l1bitsScouting"),
  HLTTriggerResults = cms.InputTag("TriggerResults", "",  "HLT"),
  cut = cms.string("pt()>20"),
  num_jets_to_fill = cms.int32(-1),
  energy_fractions = cms.VPSet( # no corrections
    cms.PSet(name=cms.string("Charged_Hadron"), label=cms.string("Charged Hadron Energy Fraction"), func=cms.string("chargedHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_Hadron"), label=cms.string("Neutral Hadron Energy Fraction"), func=cms.string("neutralHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_Hadron_noHF"), label=cms.string("Neutral Hadron Energy Fraction (no HF)"), func=cms.string("neutralHadronEnergyFraction()-HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Charged_EM"), label=cms.string("Charged EM Energy Fraction"), func=cms.string("chargedEmEnergyFraction()")),
    cms.PSet(name=cms.string("Neutral_EM"), label=cms.string("Neutral EM Energy Fraction"), func=cms.string("neutralEmEnergyFraction()")),
    cms.PSet(name=cms.string("Electron"), label=cms.string("Electron Energy Fraction"), func=cms.string("electronEnergyFraction()")),
    cms.PSet(name=cms.string("Photon"), label=cms.string("Photon Energy Fraction"), func=cms.string("photonEnergyFraction()")),
    cms.PSet(name=cms.string("Muon"), label=cms.string("Muon Energy Fraction"), func=cms.string("muonEnergyFraction()")),
    cms.PSet(name=cms.string("HF_Hadron"), label=cms.string("HF Hadron Energy Fraction"), func=cms.string("HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("HF_EM"), label=cms.string("HF EM Energy Fraction"), func=cms.string("HFEMEnergyFraction()")),
    cms.PSet(name=cms.string("Hadron"), label=cms.string("Hadron Energy Fraction"), func=cms.string("chargedHadronEnergyFraction()+neutralHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Hadron_noHF"), label=cms.string("Hadron Energy Fraction (noHF)"), func=cms.string("chargedHadronEnergyFraction()+neutralHadronEnergyFraction()-HFHadronEnergyFraction()")),
    cms.PSet(name=cms.string("Lepton"), label=cms.string("Lepton Energy Fraction"), func=cms.string("photonEnergyFraction()+muonEnergyFraction()")),
  ),
  **trigger_task
)
process.path2 = cms.Path(process.l1bit_sequence + process.offlinePFJet + process.offlinePFJetTightLeptonVetoId + process.OfflinePFJet)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletio
