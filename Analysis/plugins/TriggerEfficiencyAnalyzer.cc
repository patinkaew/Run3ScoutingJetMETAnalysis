// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// File IO include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

// Trigger bits
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

// Math
#include "DataFormats/Math/interface/deltaR.h"

// cut and function
//#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

// ROOT include files
#include "TH1.h"

// C++ include files
#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

template <typename JetType, typename METType, typename MuonType>
class TriggerEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using JetCollection = std::vector<JetType>;
  using METCollection = std::vector<METType>;
  using MuonCollection = std::vector<MuonType>;

  TriggerEfficiencyAnalyzer(const edm::ParameterSet&);
  ~TriggerEfficiencyAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  bool isGoodMuon(MuonType const &);
  bool isGoodJet(JetType const &, MuonCollection const &);
  
  // input tokens 
  edm::EDGetTokenT<JetCollection> jet_collection_token_;
  edm::EDGetTokenT<METCollection> met_collection_token_;
  edm::EDGetTokenT<MuonCollection> muon_collection_token_;

  // trigger token
  edm::EDGetTokenT<edm::TriggerResults> l1TriggerResults_token_;
  edm::EDGetTokenT<edm::TriggerResults> hltTriggerResults_token_;

  // functor for pt
  StringObjectFunction<JetType> jet_pt_func_;
  StringObjectFunction<METType> met_pt_func_;

  // file IO
  edm::Service<TFileService> fs_;
  
  // signal triggers
  std::vector<std::vector<std::string>> signal_triggers_;

  // reference trigger
  std::vector<std::string> reference_trigger_;

  // histograms
  std::vector<TH1D*> jet_pt_histograms_;
  std::vector<TH1D*> event_HT_histograms_;
  std::vector<TH1D*> met_pt_histograms_;
};

//
// constructors and destructor
//
template <typename JetType, typename METType, typename MuonType>
TriggerEfficiencyAnalyzer<JetType, METType, MuonType>::TriggerEfficiencyAnalyzer(const edm::ParameterSet& iConfig)
    : jet_collection_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("jet"))),
      met_collection_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("met"))),
      muon_collection_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("muon"))),
      l1TriggerResults_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("L1TriggerResults"))),
      hltTriggerResults_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("HLTTriggerResults"))),
      jet_pt_func_(iConfig.getUntrackedParameter<std::string>("jet_pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
      met_pt_func_(iConfig.getUntrackedParameter<std::string>("met_pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")) {

  auto reference_trigger_pset = iConfig.getParameter<edm::ParameterSet>("reference_trigger");
  reference_trigger_ = reference_trigger_pset.getParameter<std::vector<std::string>>("expr");
  std::string reference_trigger_name = reference_trigger_pset.getParameter<std::string>("name");
  if (reference_trigger_name.empty()) reference_trigger_name = boost::algorithm::join(reference_trigger_, "-or-");

  auto signal_trigger_vpset = iConfig.getParameter<std::vector<edm::ParameterSet>>("signal_triggers");
  std::vector<std::string> signal_trigger_names;
  for (auto const &signal_trigger_pset : signal_trigger_vpset) {
    std::vector<std::string> signal_trigger = signal_trigger_pset.getParameter<std::vector<std::string>>("expr");
    std::string signal_trigger_name = signal_trigger_pset.getParameter<std::string>("name");
    if (signal_trigger_name.empty()) signal_trigger_name = boost::algorithm::join(signal_trigger, "-or-");
    signal_triggers_.push_back(signal_trigger);
    signal_trigger_names.push_back(signal_trigger_name);
  }
        
  double pt_bins[] = {0., 20., 40., 50, 60, 70, 80., 90, 100., 110, 120.,  130, 140., 150, 160.,
         170, 180.,  190, 200.,  210, 220., 230, 240., 250, 260., 280., 300.,  320., 340.,
         360., 380.,  400., 450., 500., 600., 800, 1000, 1200, 1400, 1800};
  double ht_bins[] = {0., 20.,   40.,   60.,   80.,  100.,  120.,  140.,  160., 
         180.,  200.,  220.,  240.,  260.,  280.,  300.,  320.,  340.,
         360.,  380.,  400., 450., 500., 600., 700, 800, 900, 1000, 1200, 1400, 1800, 2000, 2500, 3000};
  double met_bins[] = {0., 20.,   40.,   60.,   80.,  100.,  120.,  140.,  160., 
         180.,  200.,  220.,  240.,  260.,  280.,  300.,  320.,  340.,
         360.,  380.,  400., 450., 500., 600., 700, 800, 900, 1000, 1200, 1400, 1800, 2000, 
         2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 4000, 4500, 5000, 6000};

  // initialise histograms
  // directories
  TFileDirectory jet_pt_dir = fs_->mkdir("jet_pt");
  TFileDirectory event_HT_dir = fs_->mkdir("event_HT");
  TFileDirectory met_pt_dir = fs_->mkdir("met_pt");

  // no trigger at index 0
  std::string no_trigger_name = "No_Trigger";
  jet_pt_histograms_.push_back(jet_pt_dir.make<TH1D>(no_trigger_name.c_str(), ";pT (GeV); Events", sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins));
  event_HT_histograms_.push_back(event_HT_dir.make<TH1D>(no_trigger_name.c_str(), ";HT (GeV); Events", sizeof(ht_bins)/sizeof(ht_bins[0])-1, ht_bins));
  met_pt_histograms_.push_back(met_pt_dir.make<TH1D>(no_trigger_name.c_str(), "; pT (GeV); Events", sizeof(met_bins)/sizeof(met_bins[0])-1, met_bins));

  // reference trigger at index 1
  jet_pt_histograms_.push_back(jet_pt_dir.make<TH1D>(reference_trigger_name.c_str(), ";pT (GeV); Events", sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins));
  event_HT_histograms_.push_back(event_HT_dir.make<TH1D>(reference_trigger_name.c_str(), ";HT (GeV); Events", sizeof(ht_bins)/sizeof(ht_bins[0])-1, ht_bins));
  met_pt_histograms_.push_back(met_pt_dir.make<TH1D>(reference_trigger_name.c_str(), "; pT (GeV); Events", sizeof(met_bins)/sizeof(met_bins[0])-1, met_bins));

  // signal triggers from index 2
  for (const auto &signal_trigger_name : signal_trigger_names) {
    // signal trigger
    jet_pt_histograms_.push_back(jet_pt_dir.make<TH1D>(signal_trigger_name.c_str(), ";pT (GeV); Events", sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins));
    event_HT_histograms_.push_back(event_HT_dir.make<TH1D>(signal_trigger_name.c_str(), ";HT (GeV); Events", sizeof(ht_bins)/sizeof(ht_bins[0])-1, ht_bins));
    met_pt_histograms_.push_back(met_pt_dir.make<TH1D>(signal_trigger_name.c_str(), "; pT (GeV); Events", sizeof(met_bins)/sizeof(met_bins[0])-1, met_bins));

    // signal AND reference trigger
    std::string intersection_trigger_name = "[" + signal_trigger_name + "]-and-[" + reference_trigger_name + "]";
    jet_pt_histograms_.push_back(jet_pt_dir.make<TH1D>(intersection_trigger_name.c_str(), ";pT (GeV); Events", sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins));
    event_HT_histograms_.push_back(event_HT_dir.make<TH1D>(intersection_trigger_name.c_str(), ";HT (GeV); Events", sizeof(ht_bins)/sizeof(ht_bins[0])-1, ht_bins));
    met_pt_histograms_.push_back(met_pt_dir.make<TH1D>(intersection_trigger_name.c_str(), "; pT (GeV); Events", sizeof(met_bins)/sizeof(met_bins[0])-1, met_bins));
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename JetType, typename METType, typename MuonType>
void TriggerEfficiencyAnalyzer<JetType, METType, MuonType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("jet");
  desc.addUntracked<edm::InputTag>("met");
  desc.addUntracked<edm::InputTag>("muon");
  desc.addUntracked<edm::InputTag>("L1TriggerResults", edm::InputTag("l1bits"));
  desc.addUntracked<edm::InputTag>("HLTTriggerResults", edm::InputTag("TriggerResults", "", "HLT"));

  desc.addUntracked<std::string>("jet_pt_func", "pt()");
  desc.addUntracked<std::string>("met_pt_func", "pt()");
  desc.addUntracked<bool>("lazy_eval", false);
    
  edm::ParameterSetDescription trigger_desc;
  trigger_desc.add<std::string>("name", "");
  trigger_desc.add<std::vector<std::string>>("expr");
  desc.add<edm::ParameterSetDescription>("reference_trigger", trigger_desc);
  
  edm::ParameterSet trigger_pset;
  trigger_pset.addParameter<std::string>("name", "");
  trigger_pset.addParameter<std::vector<std::string>>("expr", {"DST_PFScouting_JetHT"});
  std::vector<edm::ParameterSet> trigger_vpset;
  trigger_vpset.push_back(trigger_pset);
  desc.addVPSet("signal_triggers", trigger_desc, trigger_vpset);
  descriptions.addWithDefaultLabel(desc);
}

template <typename JetType, typename METType, typename MuonType>
bool TriggerEfficiencyAnalyzer<JetType, METType, MuonType>::isGoodMuon(MuonType const &muon) {
  return true;
}

template <typename JetType, typename METType, typename MuonType>
bool TriggerEfficiencyAnalyzer<JetType, METType, MuonType>::isGoodJet(JetType const &jet, MuonCollection const &muon_collection) {
  return true;
}

// ------------ method called for each event  ------------
template <typename JetType, typename METType, typename MuonType>
void TriggerEfficiencyAnalyzer<JetType, METType, MuonType>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // retrieve trigger decisions
  auto l1TriggerResults_handle = iEvent.getHandle(l1TriggerResults_token_);
  if (!l1TriggerResults_handle.isValid()) {
    edm::LogWarning ("Handle") << "L1 TriggerResults is invalid";
    return;
  }

  auto hltTriggerResults_handle = iEvent.getHandle(hltTriggerResults_token_);
  if (!hltTriggerResults_handle.isValid()) {
    edm::LogWarning ("Handle") << "HLT TriggerResults is invalid";
    return;
  }

  auto hltTriggerResultsByName = iEvent.triggerResultsByName(*hltTriggerResults_handle);

  // retrieve objects
  auto muon_collection_handle = iEvent.getHandle(muon_collection_token_);
  if (!muon_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Muon is invalid";
    return;
  }

  auto jet_collection_handle = iEvent.getHandle(jet_collection_token_);
  if (!jet_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Jet is invalid";
    return;
  }

  auto met_collection_handle = iEvent.getHandle(met_collection_token_);
  if (!met_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "MET is invalid";
  }
  if (met_collection_handle->size() != 1) {
    edm::LogError ("MET") << "MET collection must contain exactly one object, but get " << met_collection_handle->size();
    return; 
  }

  // select muons
  auto muon_collection_good_ptr = std::make_unique<MuonCollection>();
  for (const auto &muon : *muon_collection_handle) {
    if (isGoodMuon(muon)) {
      muon_collection_good_ptr->push_back(muon);
    }
  }

  //std::cout << "num muon: " << muon_collection_handle->size() << " good: " << muon_collection_good_ptr->size() << std::endl;
  if (muon_collection_good_ptr->empty()) {
    return;
  }

  // select jets
  auto jet_collection_good_ptr = std::make_unique<JetCollection>();
  for (const auto &jet : *jet_collection_handle) {
    if (isGoodJet(jet, *muon_collection_good_ptr)) {
      jet_collection_good_ptr->push_back(jet);
    }
  }
  if (jet_collection_good_ptr->empty()) {
    return;
  }

  //std::cout << "finish selection" << std::endl;

  // compute HT
  double HT = 0.0;
  for (const auto &jet : *jet_collection_good_ptr) HT += jet.pt();

  // get leading jet pt
  double leading_jet_pt = jet_pt_func_(jet_collection_good_ptr->at(0));

  // get met pt
  double met_pt = met_pt_func_(met_collection_handle->at(0));

  // fill histograms
  int i_histogram = 0;

  // no trigger
  jet_pt_histograms_[i_histogram]->Fill(leading_jet_pt);
  event_HT_histograms_[i_histogram]->Fill(HT);
  met_pt_histograms_[i_histogram]->Fill(met_pt);
  i_histogram++;

  // reference trigger
  bool reference_trigger_accept = util::isAnyTriggerAccept(reference_trigger_, *l1TriggerResults_handle, hltTriggerResultsByName);
  //std::cout << "reference trigger accept: " << reference_trigger_accept << std::endl;
  if (reference_trigger_accept) {
    jet_pt_histograms_[i_histogram]->Fill(leading_jet_pt);
    event_HT_histograms_[i_histogram]->Fill(HT);
    met_pt_histograms_[i_histogram]->Fill(met_pt);
  }
  i_histogram++;

  // signal triggers
  for (const auto &signal_trigger : signal_triggers_) {
    bool signal_trigger_accept = util::isAnyTriggerAccept(signal_trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
    if (signal_trigger_accept) {
      jet_pt_histograms_[i_histogram]->Fill(leading_jet_pt);
      event_HT_histograms_[i_histogram]->Fill(HT);
      met_pt_histograms_[i_histogram]->Fill(met_pt);
    }
    i_histogram++;

    if (signal_trigger_accept && reference_trigger_accept) {
      jet_pt_histograms_[i_histogram]->Fill(leading_jet_pt);
      event_HT_histograms_[i_histogram]->Fill(HT);
      met_pt_histograms_[i_histogram]->Fill(met_pt);
    }
    i_histogram++;
  }

 /*
  auto hltTriggerResults_handle = iEvent.getHandle(hltTriggerResults_token_);
  if (hltTriggerResults_handle.isValid()) {
    //auto trigger_names = iEvent.triggerNames(*triggerResults_handle);
    auto trigger_results_by_name = iEvent.triggerResultsByName(*hltTriggerResults_handle);
    for (size_t itrigger = 0; itrigger < trigger_results_by_name.size(); itrigger++) {
      std::string trigger_name = trigger_results_by_name.triggerName(itrigger);
      if (trigger_name.compare(0, 3, "DST") == 0 && trigger_name.find("Scouting") != std::string::npos) {
        //std::cout << "path idx: " << itrigger << " name: " << trigger_name << " accept: " << trigger_results_by_name.accept(itrigger) << std::endl;
      }
    }
  }

  auto l1TriggerResults_handle = iEvent.getHandle(l1TriggerResults_token_);
  if (l1TriggerResults_handle.isValid()) {
    auto trigger_names = l1TriggerResults_handle->getTriggerNames();
    for (size_t itrigger = 0; itrigger < trigger_names.size(); itrigger++) {
      std::string trigger_name = trigger_names[itrigger];
      bool trigger_result = l1TriggerResults_handle->accept(itrigger);
      if (trigger_result){
        //std::cout << "path idx: " << itrigger << " name: " << trigger_name << " accept: " << l1TriggerResults_handle->accept(itrigger) << std::endl;
      }
    }
  }
  */
}

/*
float dr2(float eta1, float eta2, float phi1, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  return deta * deta + dphi * dphi
}

double dr2(double eta1, double eta2, double phi1, double phi2) {
  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;
  return deta * deta + dphi * dphi
}
*/

// Run3ScoutingPFJet / Run3ScoutingMuon TriggerEfficiency
// Run3Scouting dataformat include files
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "AnalysisDataFormats/Scouting/interface/Run3ScoutingPFMET.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

using ScoutingPFJetTriggerEfficiencyAnalyzer = TriggerEfficiencyAnalyzer<Run3ScoutingPFJet, Run3ScoutingPFMET, Run3ScoutingMuon>;
DEFINE_FWK_MODULE(ScoutingPFJetTriggerEfficiencyAnalyzer);

bool isGoodScoutingMuon(Run3ScoutingMuon const &scoutingMuon) {
  if (scoutingMuon.pt() > 10
      && abs(scoutingMuon.eta()) < 2.4
      && abs(scoutingMuon.trk_dxy()) < 0.2
      && abs(scoutingMuon.trackIso()) < 0.15
      && abs(scoutingMuon.trk_dz()) < 0.5
      && scoutingMuon.normalizedChi2() < 10
      && scoutingMuon.nValidRecoMuonHits() > 0
      && scoutingMuon.nRecoMuonMatchedStations() > 1
      && scoutingMuon.nValidPixelHits() > 0
      && scoutingMuon.nTrackerLayersWithMeasurement() > 5
      ) {
    return true;
  } else {
    return false;
  }
}

template <>
bool ScoutingPFJetTriggerEfficiencyAnalyzer::isGoodMuon(Run3ScoutingMuon const &scoutingMuon) {
  return isGoodScoutingMuon(scoutingMuon);
}

using Run3ScoutingMuonCollection = std::vector<Run3ScoutingMuon>;
template <>
bool ScoutingPFJetTriggerEfficiencyAnalyzer::isGoodJet(Run3ScoutingPFJet const &scoutingPFJet, Run3ScoutingMuonCollection const &scoutingMuon_collection) {
  if (scoutingPFJet.pt() > 30 && abs(scoutingPFJet.eta()) < 2.5) {
    for (const auto &scoutingMuon : scoutingMuon_collection) { // dr(jet, muon) >= 0.4 for all muons
      if (deltaR2<Run3ScoutingPFJet, Run3ScoutingMuon>(scoutingPFJet, scoutingMuon) <= 0.4 * 0.4) {
        return false;
      }
    }
    return true;
  }
  return false;
}

// PAT jets / PAT Muon
// pat dataformat include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

using PATJetTriggerEfficiencyAnalyzer = TriggerEfficiencyAnalyzer<pat::Jet, pat::MET, pat::Muon>;
DEFINE_FWK_MODULE(PATJetTriggerEfficiencyAnalyzer);

bool isGoodPATMuon(pat::Muon const &patMuon) {
  if (patMuon.pt() > 10 && abs(patMuon.eta()) < 2.4) {
    double pfRelIso04_all = (patMuon.pfIsolationR04().sumChargedHadronPt + std::max(patMuon.pfIsolationR04().sumNeutralHadronEt + patMuon.pfIsolationR04().sumPhotonEt - patMuon.pfIsolationR04().sumPUPt/2, 0.0f))/patMuon.pt();
    bool looseId = patMuon.passed(reco::Muon::CutBasedIdLoose);
    if (pfRelIso04_all < 0.25 && looseId) {
      return true;
    }
  }
  return false;
}

template <>
bool PATJetTriggerEfficiencyAnalyzer::isGoodMuon(pat::Muon const &patMuon) {
  return isGoodPATMuon(patMuon);
}

using PATMuonCollection = std::vector<pat::Muon>;
template <>
bool PATJetTriggerEfficiencyAnalyzer::isGoodJet(pat::Jet const &patJet, PATMuonCollection const &patMuon_collection) {
  if (patJet.pt() > 30 && abs(patJet.eta()) < 2.5) {
      for (const auto &patMuon : patMuon_collection) { // dr(jet, muon) >= 0.4 for all muons
        if (deltaR2<pat::Jet, pat::Muon>(patJet, patMuon) <= 0.4 * 0.4) {
          return false;
        }
      }
    return true;
  }
  return false;
}

// reco::PFJet / PAT Muon
// reco dataformat include files
#include "DataFormats/JetReco/interface/PFJet.h"

using RecoPFJetTriggerEfficiencyAnalyzer = TriggerEfficiencyAnalyzer<reco::PFJet, pat::MET, pat::Muon>;
DEFINE_FWK_MODULE(RecoPFJetTriggerEfficiencyAnalyzer);

template <>
bool RecoPFJetTriggerEfficiencyAnalyzer::isGoodMuon(pat::Muon const &patMuon) {
  return isGoodPATMuon(patMuon);
}

template <>
bool RecoPFJetTriggerEfficiencyAnalyzer::isGoodJet(reco::PFJet const& recoPFJet, PATMuonCollection const &patMuon_collection) {
  if (recoPFJet.pt() > 30 && abs(recoPFJet.eta()) < 2.5) {
      for (const auto &patMuon : patMuon_collection) { // dr(jet, muon) >= 0.4 for all muons
        if (deltaR2<reco::PFJet, pat::Muon>(recoPFJet, patMuon) <= 0.4 * 0.4) {
          return false;
        }
      }
    return true;
  }
  return false;
}
