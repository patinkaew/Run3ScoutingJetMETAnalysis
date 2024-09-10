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
#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// cut and function
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

// ROOT include files
#include "TH1.h"

// C++ include files
#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

template <typename JetType>
class DijetTriggerEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using JetCollection = std::vector<JetType>; 

  DijetTriggerEfficiencyAnalyzer(const edm::ParameterSet&);
  ~DijetTriggerEfficiencyAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  bool isGoodTagAndProbeDijetPair(const JetType &, const JetType &);

  // input tokens
  edm::EDGetTokenT<JetCollection> jet_collection_token_;

  // trigger tokens
  edm::EDGetTokenT<edm::TriggerResults> l1TriggerResults_token_;
  edm::EDGetTokenT<edm::TriggerResults> hltTriggerResults_token_;

  // functor
  const StringCutObjectSelector<JetType> tag_jet_cut_; // cut applied to tag jet
  const StringCutObjectSelector<JetType> probe_jet_cut_; // cut applied to probe jet
  const StringObjectFunction<JetType> jet_pt_func_; // function to get jet pt
  const StringObjectFunction<JetType> jet_eta_func_; // function to get jet eta
  const StringObjectFunction<JetType> jet_phi_func_; // function to get jet phi
  const StringObjectFunction<JetType> jet_mass_func_; // function to get jet mass

  // file IO
  edm::Service<TFileService> fs_;

  // triggers
  std::vector<std::vector<std::string>> triggers_; // trigger expressions

  // histograms
  std::vector<TH1D*> dijet_mass_histograms_;
};

template <typename JetType>
DijetTriggerEfficiencyAnalyzer<JetType>::DijetTriggerEfficiencyAnalyzer(const edm::ParameterSet& iConfig)
  : jet_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("jet"))),
    l1TriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("L1TriggerResults"))),
    hltTriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("HLTTriggerResults"))),
    tag_jet_cut_(iConfig.getParameter<std::string>("tag_jet_cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    probe_jet_cut_(iConfig.getParameter<std::string>("probe_jet_cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_pt_func_(iConfig.getParameter<std::string>("jet_pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")), 
    jet_eta_func_(iConfig.getParameter<std::string>("jet_eta_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")), 
    jet_phi_func_(iConfig.getParameter<std::string>("jet_phi_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")), 
    jet_mass_func_(iConfig.getParameter<std::string>("jet_mass_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")) {

  std::vector<std::string> trigger_names;

  std::string no_trigger_name = "No_Trigger";
  trigger_names.push_back(no_trigger_name);

  auto trigger_vpset = iConfig.getParameter<std::vector<edm::ParameterSet>>("signal_triggers");
  for (auto const &trigger_pset : trigger_vpset) {
    std::vector<std::string> trigger = trigger_pset.getParameter<std::vector<std::string>>("expr");
    std::string trigger_name = trigger_pset.getParameter<std::string>("name");
    if (trigger_name.empty()) trigger_name = boost::algorithm::join(trigger, "-or-");
    triggers_.push_back(trigger);
    trigger_names.push_back(trigger_name);
  }
  
  // initialize histograms
  TFileDirectory dijet_mass_dir = fs_->mkdir("dijet_mass");
  for (const auto &trigger_name : trigger_names) {
    dijet_mass_histograms_.push_back(dijet_mass_dir.make<TH1D>(trigger_name.c_str(), ";m_{jj} (GeV); Events", 100, 0., 1000.));
  }
}

template <typename JetType>
void DijetTriggerEfficiencyAnalyzer<JetType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jet");
  desc.add<edm::InputTag>("L1TriggerResults", edm::InputTag("l1bits"));
  desc.add<edm::InputTag>("HLTTriggerResults", edm::InputTag("TriggerResults", "", "HLT"));

  desc.addUntracked<bool>("lazy_eval", false);
  desc.add<std::string>("tag_jet_cut", "pt()>=30 && eta()>=-2.5 && eta()<=2.5");
  desc.add<std::string>("probe_jet_cut", "pt()>=30 && eta()>=-2.5 && eta()<=2.5");
  desc.add<std::string>("jet_pt_func", "pt()");
  desc.add<std::string>("jet_eta_func", "eta()");
  desc.add<std::string>("jet_phi_func", "phi()");
  desc.add<std::string>("jet_mass_func", "mass()");

  edm::ParameterSetDescription trigger_desc;
  trigger_desc.add<std::string>("name", "");
  trigger_desc.add<std::vector<std::string>>("expr");
  edm::ParameterSet trigger_pset;
  trigger_pset.addParameter<std::string>("name", "");
  trigger_pset.addParameter<std::vector<std::string>>("expr", {"L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"});
  std::vector<edm::ParameterSet> trigger_vpset;
  trigger_vpset.push_back(trigger_pset);
  desc.addVPSet("triggers", trigger_desc, trigger_vpset);

  descriptions.addWithDefaultLabel(desc);
}

template <typename JetType>
void DijetTriggerEfficiencyAnalyzer<JetType>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  // retrieve jet collection
  auto jet_collection_handle = iEvent.getHandle(jet_collection_token_);
  if (!jet_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Jet is invalid";
    return;
  }
  
  size_t num_jets = jet_collection_handle->size();
  // at least two jets for dijet event
  if (num_jets < 2) {
    return;
  }
  
  // loop jets to find max dijet mass
  double max_dijet_mass = -1.;
  for (unsigned int itag = 0; itag < num_jets; itag++) {
    auto tag_jet = jet_collection_handle->at(itag);
    if (!tag_jet_cut_(tag_jet)) continue;
    for (unsigned int iprobe = itag+1; iprobe < num_jets; iprobe++) {
      auto probe_jet = jet_collection_handle->at(iprobe);
      if (!probe_jet_cut_(probe_jet)) continue;
      if (isGoodTagAndProbeDijetPair(tag_jet, probe_jet)) {
        // compute dijet mass
        math::PtEtaPhiMLorentzVector tag_jet_p4(jet_pt_func_(tag_jet), jet_eta_func_(tag_jet), jet_phi_func_(tag_jet), jet_mass_func_(tag_jet));
        math::PtEtaPhiMLorentzVector probe_jet_p4(jet_pt_func_(probe_jet), jet_eta_func_(probe_jet), jet_phi_func_(probe_jet), jet_mass_func_(probe_jet));
        double dijet_mass = (tag_jet_p4 + probe_jet_p4).M();
        if (dijet_mass > max_dijet_mass) {
          max_dijet_mass = dijet_mass;
        }
      }
    }
  }

  // no valid dijet pair
  if (max_dijet_mass < 0) return;
  
  // fill histograms
  //unsigned int ihistogram = 0;
  
  // No trigger first
  dijet_mass_histograms_[0]->Fill(max_dijet_mass);

  // with trigger
  for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
    auto trigger = triggers_[itrigger];
    bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
    if (trigger_accept) {
      dijet_mass_histograms_[itrigger+1]->Fill(max_dijet_mass);
    }
  }
}

template <typename JetType>
bool DijetTriggerEfficiencyAnalyzer<JetType>::isGoodTagAndProbeDijetPair(JetType const &tag_jet, JetType const &probe_jet) {
  double deta = abs(jet_eta_func_(tag_jet) - jet_eta_func_(probe_jet));
  double dphi = deltaPhi(jet_phi_func_(tag_jet), jet_phi_func_(probe_jet));
  return (deta < 1.5) && (dphi > 2.7); // deta from L1 and dphi from usual dijet selection (should drop dphi > 2.7 for trigger efficiency measurement?)
}

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
using ScoutingPFJetDijetTriggerEfficiencyAnalyzer = DijetTriggerEfficiencyAnalyzer<Run3ScoutingPFJet>;
DEFINE_FWK_MODULE(ScoutingPFJetDijetTriggerEfficiencyAnalyzer);

#include "DataFormats/PatCandidates/interface/Jet.h"
using PATJetDijetTriggerEfficiencyAnalyzer = DijetTriggerEfficiencyAnalyzer<pat::Jet>;
DEFINE_FWK_MODULE(PATJetDijetTriggerEfficiencyAnalyzer);

#include "DataFormats/JetReco/interface/PFJet.h"
using RecoPFJetDijetTriggerEfficiencyAnalyzer = DijetTriggerEfficiencyAnalyzer<reco::PFJet>;
DEFINE_FWK_MODULE(RecoPFJetDijetTriggerEfficiencyAnalyzer);
