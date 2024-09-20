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
//#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"

// cut and function
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

// ROOT include files
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

// C++ include files
#include <iostream>
#include <cmath>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

template <typename JetType>
class JetEnergyFractionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using JetCollection = std::vector<JetType>; 

  JetEnergyFractionAnalyzer(const edm::ParameterSet&);
  ~JetEnergyFractionAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  // input tokens
  edm::EDGetTokenT<JetCollection> jet_collection_token_;

  // trigger tokens
  edm::EDGetTokenT<edm::TriggerResults> l1TriggerResults_token_;
  edm::EDGetTokenT<edm::TriggerResults> hltTriggerResults_token_;

  // file IO
  edm::Service<TFileService> fs_;

  // functor
  const StringCutObjectSelector<JetType> cut_; // cut applied to select jets
  const StringObjectFunction<JetType> jet_pt_func_;
  const StringObjectFunction<JetType> jet_eta_func_;
  const StringObjectFunction<JetType> jet_phi_func_;
  /*const StringObjectFunction<JetType> jet_chf_func_; // charged hadron energy fraction
  const StringObjectFunction<JetType> jet_nhf_func_; // neutral hadron energy fraction
  const StringObjectFunction<JetType> jet_cef_func_; // charged em energy fraction 
  const StringObjectFunction<JetType> jet_nef_func_; // neutral em energy fraction 
  const StringObjectFunction<JetType> jet_muf_func_; // muon energy fraction 
  const StringObjectFunction<JetType> jet_fhf_func_; // hf hadron energy fraction 
  const StringObjectFunction<JetType> jet_fef_func_; // hf em energy fraction
  */
  std::vector<StringObjectFunction<JetType>> jet_fraction_funcs_; // all energy fractions in a vector

  // triggers
  std::vector<std::vector<std::string>> triggers_; // trigger expressions

  // other parameters
  int num_jets_to_fill_;
  std::vector<double> eta_ranges_;
  size_t num_fractions_;

  // histograms
  std::vector<TH1D*> ef_histograms_; // for rough eta range, inclusive phi, inclusive pt
  std::vector<TH2D*> pt_ef_histograms_; // for rough eta range, inclusive phi
  std::vector<TProfile*> pt_ef_profiles_; // for rough eta range, inclusive phi , profile on ef
  /*std::vector<TH3D*> eta_pt_ef_histograms_; // inclusive phi
  std::vector<TProfile2D*> eta_pt_ef_profiles_; // inclusive phi, profile on ef
  std::vector<THnD*> eta_phi_pt_ef_histograms_; // 4D
  std::vector<TProfile3D*> eta_phi_pt_ef_profiles_; // profile on ef */
};

template <typename JetType>
JetEnergyFractionAnalyzer<JetType>::JetEnergyFractionAnalyzer(const edm::ParameterSet& iConfig)
  : jet_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    l1TriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("L1TriggerResults"))),
    hltTriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("HLTTriggerResults"))),
    cut_(iConfig.getParameter<std::string>("cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_pt_func_(iConfig.getParameter<std::string>("pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_eta_func_(iConfig.getParameter<std::string>("eta_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_phi_func_(iConfig.getParameter<std::string>("phi_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    /*jet_chf_func_(iConfig.getParameter<std::string>("chf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_nhf_func_(iConfig.getParameter<std::string>("nhf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_cef_func_(iConfig.getParameter<std::string>("cef_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_nef_func_(iConfig.getParameter<std::string>("nef_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_elf_func_(iConfig.getParameter<std::string>("elf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_phf_func_(iConfig.getParameter<std::string>("phf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_muf_func_(iConfig.getParameter<std::string>("muf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_fhf_func_(iConfig.getParameter<std::string>("fhf_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_fef_func_(iConfig.getParameter<std::string>("fef_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),*/
    num_jets_to_fill_(iConfig.getParameter<int>("num_jets_to_fill")) {

  // parse energy fraction functions
  bool lazy_eval = iConfig.getUntrackedParameter<bool>("lazy_eval");
  std::vector<std::string> jet_fraction_names = {}; //{"Charged-Hadron", "Neutral-Hadron", "Charged-EM", "Neutral-EM", "Electron", "Photon", "Muon", "HF-Hadron", "HF-EM"};
  std::vector<std::string> jet_fraction_labels = {};
  jet_fraction_funcs_ = {};
  auto fraction_vpset = iConfig.getParameter<std::vector<edm::ParameterSet>>("energy_fractions");
  for (auto const &fraction_pset : fraction_vpset) {
    std::string fraction_func = fraction_pset.getParameter<std::string>("func");
    jet_fraction_funcs_.push_back(StringObjectFunction<JetType>(fraction_func, lazy_eval));
    std::string fraction_name = fraction_pset.getParameter<std::string>("name"); // name in root file
    if (fraction_name.empty()) fraction_name = fraction_func;
    std::string fraction_label = fraction_pset.getParameter<std::string>("label"); // label on histogram
    if (fraction_label.empty()) fraction_label = fraction_name;
    jet_fraction_names.push_back(fraction_name);
    jet_fraction_labels.push_back(fraction_label);
  }

  // parse triggers
  std::vector<std::string> trigger_names;

  std::string no_trigger_name = "No_Trigger";
  trigger_names.push_back(no_trigger_name);

  auto trigger_vpset = iConfig.getParameter<std::vector<edm::ParameterSet>>("triggers");
  for (auto const &trigger_pset : trigger_vpset) {
    std::vector<std::string> trigger = trigger_pset.getParameter<std::vector<std::string>>("expr");
    std::string trigger_name = trigger_pset.getParameter<std::string>("name");
    if (trigger_name.empty()) trigger_name = boost::algorithm::join(trigger, "-or-");
    triggers_.push_back(trigger);
    trigger_names.push_back(trigger_name);
  }

  eta_ranges_ = {0.0, 1.3, 2.5, 3.0, 5.0};
  double pt_bins[] = {15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 272, 330, 395, 468, 548, 686, 846, 1032, 1248, 1588, 2000, 2500, 3103};
  // initialize directories and histograms
  for (auto const &trigger_name : trigger_names) {
    TFileDirectory trigger_dir = fs_->mkdir(trigger_name.c_str());
    for (unsigned int ieta = 0; ieta < eta_ranges_.size()-1; ieta++) {
      std::string eta_range = fmt::format("{:.2f}",eta_ranges_[ieta])+"<|eta|<"+fmt::format("{:.2f}",eta_ranges_[ieta+1]);
      TFileDirectory eta_range_dir = trigger_dir.mkdir(eta_range);
      TFileDirectory ef_histogram_dir = eta_range_dir.mkdir("Energy_Fraction");
      TFileDirectory pt_ef_histogram_dir = eta_range_dir.mkdir("Energy_Fraction_vs_pT");
      TFileDirectory pt_ef_profile_dir = eta_range_dir.mkdir("mean_Energy_Fraction_vs_pT");
      for (unsigned int ifraction = 0; ifraction < jet_fraction_funcs_.size(); ifraction++) {
      //for (auto const &jet_fraction_name : jet_fraction_names) {
        std::string jet_fraction_name = jet_fraction_names[ifraction];
        std::string jet_fraction_label = jet_fraction_labels[ifraction];
        ef_histograms_.push_back(ef_histogram_dir.make<TH1D>(jet_fraction_name.c_str(), fmt::format(";{};Events", jet_fraction_label).c_str(), 40, 0., 1.));
        pt_ef_histograms_.push_back(pt_ef_histogram_dir.make<TH2D>(jet_fraction_name.c_str(), (";p_{T}"+fmt::format(" (GeV);{};Events", jet_fraction_label)).c_str(), sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins, 40, 0., 1.));
        pt_ef_profiles_.push_back(pt_ef_profile_dir.make<TProfile>(jet_fraction_name.c_str(), (";p_{T}"+fmt::format(" (GeV);{}",jet_fraction_label)).c_str(), sizeof(pt_bins)/sizeof(pt_bins[0])-1, pt_bins));
      }
    }
  }
}

template <typename JetType>
void JetEnergyFractionAnalyzer<JetType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<edm::InputTag>("L1TriggerResults", edm::InputTag("l1bits"));
  desc.add<edm::InputTag>("HLTTriggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<int>("num_jets_to_fill", -1); // -1 for all objects
  desc.addUntracked<bool>("lazy_eval", false);
  desc.add<std::string>("cut", "");

  desc.add<std::string>("pt_func", "pt()");
  desc.add<std::string>("eta_func", "eta()");
  desc.add<std::string>("phi_func", "phi()");
  /*desc.add<std::string>("chf_func", "chargedHadronEnergyFraction()");
  desc.add<std::string>("nhf_func", "neutralHadronEnergyFraction()");
  desc.add<std::string>("cef_func", "chargedEmEnergyFraction()");
  desc.add<std::string>("nef_func", "neutralEmEnergyFraction()");
  desc.add<std::string>("elf_func", "electronEnergyFraction()");
  desc.add<std::string>("phf_func", "photonEnergyFraction()");
  desc.add<std::string>("muf_func", "muonEnergyFraction()");
  desc.add<std::string>("fhf_func", "HFHadronEnergyFraction()");
  desc.add<std::string>("fef_func", "HFEMEnergyFraction()");*/
  
  edm::ParameterSetDescription fraction_desc;
  fraction_desc.add<std::string>("name", "");
  fraction_desc.add<std::string>("label", "");
  fraction_desc.add<std::string>("func");

  std::vector<edm::ParameterSet> fraction_vpset;
  std::vector<std::string> fraction_default_names = { 
                                                      "Charged_Hadron", "Neutral_Hadron",
                                                      "Charged_EM", "Neutral_EM",
                                                      "Electron", "Photon", "Muon",
                                                      "HF_Hadron", "HF_EM"
                                                    };
  std::vector<std::string> fraction_default_labels = { 
                                                       "Charged Hadron Energy Fraction", "Neutral Hadron Energy Fraction",
                                                       "Charged EM Energy Fraction", "Neutral EM Energy Fraction",
                                                       "Electron Energy Fraction", "Photon Energy Fraction", "Muon Energy Fraction",
                                                       "HF Hadron Energy Fraction", "HF EM Energy Fraction"
                                                    };
  std::vector<std::string> fraction_default_funcs = {
                                                      "chargedHadronEnergyFraction()", "neutralHadronEnergyFraction()", "neutralHadronEnergyFraction()-HFHadronEnergyFraction()",
                                                      "chargedEmEnergyFraction()", "neutralEmEnergyFraction()",
                                                      "electronEnergyFraction()", "photonEnergyFraction()", "muonEnergyFraction()",
                                                      "HFHadronEnergyFraction()", "HFEMEnergyFraction()"
                                                    };
  for (unsigned int ifraction = 0; ifraction < fraction_default_names.size(); ifraction++) {
    edm::ParameterSet fraction_pset;
    fraction_pset.addParameter<std::string>("name", fraction_default_names[ifraction]);
    fraction_pset.addParameter<std::string>("label", fraction_default_labels[ifraction]);
    fraction_pset.addParameter<std::string>("func", fraction_default_funcs[ifraction]);
    fraction_vpset.push_back(fraction_pset);
  }
  desc.addVPSet("energy_fractions", fraction_desc, fraction_vpset);

  edm::ParameterSetDescription trigger_desc;
  trigger_desc.add<std::string>("name", "");
  trigger_desc.add<std::vector<std::string>>("expr");
  edm::ParameterSet trigger_pset;
  trigger_pset.addParameter<std::string>("name", "");
  trigger_pset.addParameter<std::vector<std::string>>("expr", {"DST_PFScouting_JetHT"});
  std::vector<edm::ParameterSet> trigger_vpset;
  trigger_vpset.push_back(trigger_pset);
  desc.addVPSet("triggers", trigger_desc, trigger_vpset);

  descriptions.addWithDefaultLabel(desc);
}

template <typename JetType>
void JetEnergyFractionAnalyzer<JetType>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  // retrieve object
  auto jet_collection_handle = iEvent.getHandle(jet_collection_token_);
  if (!jet_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Jet is invalid";
    return;
  }
  
  int num_jets_filled = 0;
  size_t num_jet_fractions = jet_fraction_funcs_.size();
  size_t num_eta_ranges = eta_ranges_.size()-1;
  for (unsigned int ijet = 0; ijet < jet_collection_handle->size(); ijet++) {
    auto jet = jet_collection_handle->at(ijet);
    if (cut_(jet)) {
      // search for ieta
      double abs_eta = abs(jet_eta_func_(jet));
      int jet_ieta = -1;
      for (unsigned int ieta = 0; ieta < num_eta_ranges; ieta++) { 
        if (abs_eta > eta_ranges_[ieta] && abs_eta <= eta_ranges_[ieta+1]) {
          jet_ieta = ieta;
          break;
        }
      }
      if (jet_ieta <= -1) continue; // out of eta range

      // no trigger
      for (unsigned int ifraction = 0; ifraction <  num_jet_fractions; ifraction++) {
        auto jet_fraction_func = jet_fraction_funcs_[ifraction];
        ef_histograms_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_fraction_func(jet));
        pt_ef_histograms_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt_func_(jet), jet_fraction_func(jet));
        pt_ef_profiles_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt_func_(jet), jet_fraction_func(jet));
      }
      // with trigger
      for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
        auto trigger = triggers_[itrigger];
        bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
        if (trigger_accept) {
          unsigned int offset = (1+itrigger)*num_eta_ranges*num_jet_fractions;
          for (unsigned int ifraction = 0; ifraction < num_jet_fractions; ifraction++) {
            auto jet_fraction_func = jet_fraction_funcs_[ifraction];
            ef_histograms_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_fraction_func(jet));
            pt_ef_histograms_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt_func_(jet), jet_fraction_func(jet));
            pt_ef_profiles_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt_func_(jet), jet_fraction_func(jet));
          }
        }
      }
      num_jets_filled++; // filled one more object
    }
    if ((num_jets_to_fill_ != -1) && (num_jets_filled >= num_jets_to_fill_)) {
      break;
    }
  }
}

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
using ScoutingPFJetEnergyFractionAnalyzer = JetEnergyFractionAnalyzer<Run3ScoutingPFJet>;
DEFINE_FWK_MODULE(ScoutingPFJetEnergyFractionAnalyzer);

#include "DataFormats/PatCandidates/interface/Jet.h"
using PATJetEnergyFractionAnalyzer = JetEnergyFractionAnalyzer<pat::Jet>;
DEFINE_FWK_MODULE(PATJetEnergyFractionAnalyzer);

#include "DataFormats/JetReco/interface/PFJet.h"
using RecoPFJetEnergyFractionAnalyzer = JetEnergyFractionAnalyzer<reco::PFJet>;
DEFINE_FWK_MODULE(RecoPFJetEnergyFractionAnalyzer);
