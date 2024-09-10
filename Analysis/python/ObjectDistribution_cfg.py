import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
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
    fileName = cms.string("object_distribution.root")
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

#process.load("Run3ScoutingJetMETAnalysis.Analysis.JetIdProducer_cff")
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
    ),
  }
}

# select task to run
trigger_task = trigger_task_matrix["2024"]

process.ScoutingPFCandidate = cms.EDAnalyzer("ScoutingPFCandidateDistributionAnalyzer",
  src = cms.InputTag("hltScoutingPFPacker", "", "HLT"),
  L1TriggerResults = cms.InputTag("l1bitsScouting"),
  HLTTriggerResults = cms.InputTag("TriggerResults", "",  "HLT"),
  categories = cms.VPSet(
    cms.PSet(name=cms.string("Electron"), cut=cms.string("pdgId()==11 || pdgId()==-11")),
    cms.PSet(name=cms.string("Muon"), cut=cms.string("pdgId()==13 || pdgId()==-13")),
    cms.PSet(name=cms.string("Photon"), cut=cms.string("pdgId()==22 || pdgId()==-22")),
    cms.PSet(name=cms.string("Charged Hadron"), cut=cms.string("pdgId()==211 || pdgId()==-211")),
    cms.PSet(name=cms.string("Neutral Hadron"), cut=cms.string("pdgId()==130 || pdgId()==-130")),
  ),
  num_objects_to_fill = cms.int32(-1),
  **trigger_task
)
process.path0 = cms.Path(process.l1bit_sequence + process.ScoutingPFCandidate)

process.OfflinePFCandidate = cms.EDAnalyzer("PatPFCandidateDistributionAnalyzer",
  src = cms.InputTag("packedPFCandidates", "", "RECO"),
  L1TriggerResults = cms.InputTag("l1bitsScouting"),
  HLTTriggerResults = cms.InputTag("TriggerResults", "",  "HLT"),
  cut = cms.string("pt()>0.6 && eta()>-3.0 && eta()<3.0"),
  categories = cms.VPSet(
    cms.PSet(name=cms.string("Electron"), cut=cms.string("pdgId()==11 || pdgId()==-11")),
    cms.PSet(name=cms.string("Muon"), cut=cms.string("pdgId()==13 || pdgId()==-13")),
    cms.PSet(name=cms.string("Photon"), cut=cms.string("pdgId()==22 || pdgId()==-22")),
    cms.PSet(name=cms.string("Charged Hadron"), cut=cms.string("pdgId()==211 || pdgId()==-211")),
    cms.PSet(name=cms.string("Neutral Hadron"), cut=cms.string("pdgId()==130 || pdgId()==-130")),
  ),
  num_objects_to_fill = cms.int32(-1),
  **trigger_task
)
process.path1 = cms.Path(process.l1bit_sequence + process.OfflinePFCandidate)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletio
