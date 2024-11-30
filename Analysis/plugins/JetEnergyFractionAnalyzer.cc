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
//#include "THn.h"
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
  std::vector<StringObjectFunction<JetType>> jet_fraction_funcs_; // all energy fractions in a vector

  // triggers
  std::vector<std::vector<std::string>> triggers_; // trigger expressions

  // other parameters
  int num_jets_to_fill_;
  std::vector<double> eta_ranges_;
  std::vector<double> pt_ranges_;
  //const int num_phi_bins_ = 50;
  //const int num_ef_bins_ = 40;
  size_t num_fractions_;

  // histograms
  // rough eta range
  std::vector<TH1D*> ef_histograms_; // for rough eta range, inclusive phi, inclusive pt
  std::vector<TH2D*> pt_ef_histograms_; // for rough eta range, inclusive phi
  std::vector<TProfile*> pt_ef_profiles_; // for rough eta range, inclusive phi , profile on ef
  // rough pt range
  std::vector<TH2D*> eta_ef_histograms_; // inclusive phi
  std::vector<TProfile*> eta_ef_profiles_; // inclusive phi, profile on ef
  std::vector<TH3D*> eta_phi_ef_histograms_;
  std::vector<TProfile2D*> eta_phi_ef_profiles_; // profile on ef
  // all
  std::vector<TH3D*> eta_pt_ef_histograms_; // inclusive phi
  std::vector<TProfile2D*> eta_pt_ef_profiles_; // inclusive phi, profile on ef
  //std::vector<THnD*> eta_phi_pt_ef_histograms_; // 4D
  std::vector<TProfile3D*> eta_phi_pt_ef_profiles_; // profile on ef
};

template <typename JetType>
JetEnergyFractionAnalyzer<JetType>::JetEnergyFractionAnalyzer(const edm::ParameterSet& iConfig)
  : jet_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    l1TriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("L1TriggerResults"))),
    hltTriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("HLTTriggerResults"))),
    cut_(iConfig.getParameter<std::string>("cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_pt_func_(iConfig.getParameter<std::string>("jet_pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_eta_func_(iConfig.getParameter<std::string>("jet_eta_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    jet_phi_func_(iConfig.getParameter<std::string>("jet_phi_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
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
    if (fraction_name.empty())
      fraction_name = fraction_func;
    std::string fraction_label = fraction_pset.getParameter<std::string>("label"); // label on histogram
    if (fraction_label.empty())
      fraction_label = fraction_name;
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
  
  // eta binning: BB, EC1, EC2, HF 
  eta_ranges_ = {0.0, 1.3, 2.5, 3.0, 5.0};
  
  // eta binning for L2Relative
  double eta_bins[] = {
    -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139,
    -2.964, -2.853, -2.65, -2.5,
    -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305,
    -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
    0,
    0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218,
    1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322,
    2.5, 2.65, 2.853, 2.964,
    3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

  // eta binning for L2Residual
  //double eta_bins[] = {0., 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.839, 4.013, 4.583, 5.191};

  int num_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0])-1;

  // pt binning
  pt_ranges_ = {15, 50, 110, 500, 1000, 4500};

  // pt binning for energy fraction shape studies from Mikael
  double pt_bins[] = {15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 272, 330, 395, 468, 548, 686, 846, 1032, 1248, 1588, 2000, 2500, 3103};
  int num_pt_bins = sizeof(pt_bins)/sizeof(pt_bins[0])-1;

  // phi binning
  constexpr int num_phi_bins = 50;
  double phi_bins[num_phi_bins+1];
  for (int i=0; i<=num_phi_bins; i++)
    phi_bins[i] = -M_PI + 2*M_PI/num_phi_bins*i;

  // energy fraction binning
  constexpr int num_ef_bins = 40;
  double ef_bins[num_ef_bins+1];
  for (int i=0; i<=num_ef_bins; i++)
    ef_bins[i] = 0. + 1./num_ef_bins*i;
  
  // initialize directories and histograms
  for (auto const &trigger_name : trigger_names) {
    std::cout << trigger_name << std::endl;
    TFileDirectory trigger_dir = fs_->mkdir(trigger_name.c_str());
    // rough eta ranges
    for (unsigned int ieta = 0; ieta < eta_ranges_.size()-1; ieta++) {
      std::string eta_range = fmt::format("{:.2f}",eta_ranges_[ieta])+"<|eta|<"+fmt::format("{:.2f}",eta_ranges_[ieta+1]);
      TFileDirectory eta_range_dir = trigger_dir.mkdir(eta_range);
      TFileDirectory ef_histogram_dir = eta_range_dir.mkdir("Energy_Fraction");
      TFileDirectory pt_ef_histogram_dir = eta_range_dir.mkdir("Energy_Fraction_vs_pT");
      TFileDirectory pt_ef_profile_dir = eta_range_dir.mkdir("mean_Energy_Fraction_vs_pT");
      for (unsigned int ifraction = 0; ifraction < jet_fraction_funcs_.size(); ifraction++) {
        std::string jet_fraction_name = jet_fraction_names[ifraction];
        std::string jet_fraction_label = jet_fraction_labels[ifraction];
        ef_histograms_.push_back(ef_histogram_dir.make<TH1D>(jet_fraction_name.c_str(), fmt::format(";{};Events", jet_fraction_label).c_str(), 40, 0., 1.));
        pt_ef_histograms_.push_back(pt_ef_histogram_dir.make<TH2D>(jet_fraction_name.c_str(), (";p_{T}"+fmt::format(" (GeV);{};Events", jet_fraction_label)).c_str(), num_pt_bins, pt_bins, 40, 0., 1.));
        pt_ef_profiles_.push_back(pt_ef_profile_dir.make<TProfile>(jet_fraction_name.c_str(), (";p_{T}"+fmt::format(" (GeV);{}", jet_fraction_label)).c_str(), num_pt_bins, pt_bins));
      }
    }

    // rough pt ranges
    for (unsigned int ipt = 0; ipt < pt_ranges_.size()-1; ipt++){
      std::string pt_range = fmt::format("{:.2f}",pt_ranges_[ipt])+"<pt<"+fmt::format("{:.2f}",pt_ranges_[ipt+1]);
      TFileDirectory pt_range_dir = trigger_dir.mkdir(pt_range);
      TFileDirectory eta_ef_histogram_dir = pt_range_dir.mkdir("Energy_Fraction_vs_eta");
      TFileDirectory eta_ef_profile_dir = pt_range_dir.mkdir("mean_Energy_Fraction_vs_eta");
      TFileDirectory eta_phi_ef_histogram_dir = pt_range_dir.mkdir("Energy_Fraction_vs_eta_phi");
      TFileDirectory eta_phi_ef_profile_dir = pt_range_dir.mkdir("mean_Energy_Fraction_vs_eta_phi");
      for (unsigned int ifraction = 0; ifraction < jet_fraction_funcs_.size(); ifraction++) {
        std::string jet_fraction_name = jet_fraction_names[ifraction];
        std::string jet_fraction_label = jet_fraction_labels[ifraction];
        eta_ef_histograms_.push_back(eta_ef_histogram_dir.make<TH2D>(jet_fraction_name.c_str(), (";#eta"+fmt::format(";{};Events", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, 40, 0., 1.));
        eta_ef_profiles_.push_back(eta_ef_profile_dir.make<TProfile>(jet_fraction_name.c_str(), (";#eta"+fmt::format(";{}", jet_fraction_label)).c_str(), num_eta_bins, eta_bins));
        eta_phi_ef_histograms_.push_back(eta_phi_ef_histogram_dir.make<TH3D>(jet_fraction_name.c_str(), (";#eta;#phi"+fmt::format(";{};Events", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, num_phi_bins, phi_bins, num_ef_bins, ef_bins));
        eta_phi_ef_profiles_.push_back(eta_phi_ef_profile_dir.make<TProfile2D>(jet_fraction_name.c_str(), (";#eta;#phi"+fmt::format(";{}", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, num_phi_bins, phi_bins));
      }
    }

    // all bins
    TFileDirectory eta_pt_ef_histogram_dir = trigger_dir.mkdir("Energy_Fraction_vs_eta_pT");
    TFileDirectory eta_pt_ef_profile_dir = trigger_dir.mkdir("mean_Energy_Fraction_vs_eta_pT");
    //TFileDirectory eta_phi_pt_ef_histogram_dir = trigger_dir.mkdir("Energy_Fraction_vs_eta_phi_pT");
    TFileDirectory eta_phi_pt_ef_profile_dir = trigger_dir.mkdir("mean_Energy_Fraction_vs_eta_phi_pT");
    for (unsigned int ifraction = 0; ifraction < jet_fraction_funcs_.size(); ifraction++) {
      std::string jet_fraction_name = jet_fraction_names[ifraction];
      std::string jet_fraction_label = jet_fraction_labels[ifraction];
      eta_pt_ef_histograms_.push_back(eta_pt_ef_histogram_dir.make<TH3D>(jet_fraction_name.c_str(), (";#eta;p_{T}"+fmt::format(" (GeV);{};Events", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, num_pt_bins, pt_bins, num_ef_bins, ef_bins));
      eta_pt_ef_profiles_.push_back(eta_pt_ef_profile_dir.make<TProfile2D>(jet_fraction_name.c_str(), (";#eta;p_{T}"+fmt::format(" (GeV);{}", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, num_pt_bins, pt_bins));
      //eta_phi_pt_ef_histograms_.push_back(eta_phi_pt_ef_histogram_dir.make<THnD>(jet_fraction_name.c_str(), (";#eta;p_{T}"+fmt::format(" (GeV);{};Events")).c_str(), 4, num_eta_bins, eta_bins, num_pt_bins, pt_bins, 40, 0., 1.));
      eta_phi_pt_ef_profiles_.push_back(eta_phi_pt_ef_profile_dir.make<TProfile3D>(jet_fraction_name.c_str(), (";#eta;#phi;p_{T}"+fmt::format(" (GeV);{}", jet_fraction_label)).c_str(), num_eta_bins, eta_bins, num_phi_bins, phi_bins, num_pt_bins, pt_bins));
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

  desc.add<std::string>("jet_pt_func", "pt()");
  desc.add<std::string>("jet_eta_func", "eta()");
  desc.add<std::string>("jet_phi_func", "phi()");
  
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
  
  //int num_jets_filled = 0;
  unsigned int num_jets = std::max((int)jet_collection_handle->size(), num_jets_to_fill_);
  size_t num_jet_fractions = jet_fraction_funcs_.size();
  size_t num_eta_ranges = eta_ranges_.size()-1;
  size_t num_pt_ranges = pt_ranges_.size()-1;
  for (unsigned int ijet = 0; ijet < num_jets; ijet++) {
    auto jet = jet_collection_handle->at(ijet);
    if (cut_(jet)) {
      // retrieve variables
      double jet_eta = jet_eta_func_(jet);
      double jet_pt = jet_pt_func_(jet);
      double jet_phi = jet_phi_func_(jet);

      // search for ieta
      double abs_eta = abs(jet_eta);
      int jet_ieta = -1;
      for (unsigned int ieta = 0; ieta < num_eta_ranges; ieta++) { 
        if (abs_eta > eta_ranges_[ieta] && abs_eta <= eta_ranges_[ieta+1]) {
          jet_ieta = ieta;
          break;
        }
      }
      //if (jet_ieta <= -1) continue; // out of eta range
      
      // seach for ipt 
      int jet_ipt = -1;
      for (unsigned int ipt = 0; ipt < num_pt_ranges; ipt++) { 
        if (jet_pt > pt_ranges_[ipt] && jet_pt <= pt_ranges_[ipt+1]) {
          jet_ipt = ipt;
          break;
        }
      }

      // no trigger
      for (unsigned int ifraction = 0; ifraction <  num_jet_fractions; ifraction++) {
        //auto jet_fraction_func = jet_fraction_funcs_[ifraction];
        double jet_ef = jet_fraction_funcs_[ifraction](jet);
        if (jet_ieta > -1) {
          ef_histograms_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_ef);
          pt_ef_histograms_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt, jet_ef);
          pt_ef_profiles_[jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt, jet_ef);
        }

        if (jet_ipt > -1) {
          eta_ef_histograms_[jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_ef);
          eta_ef_profiles_[jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_ef);
          eta_phi_ef_histograms_[jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_phi, jet_ef);
          eta_phi_ef_profiles_[jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_phi, jet_ef);
        }
        
        eta_pt_ef_histograms_[ifraction]->Fill(jet_eta, jet_pt, jet_ef);
        eta_pt_ef_profiles_[ifraction]->Fill(jet_eta, jet_pt, jet_ef);
        eta_phi_pt_ef_profiles_[ifraction]->Fill(jet_eta, jet_phi, jet_pt, jet_ef);
      }
      // with trigger
      for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
        auto trigger = triggers_[itrigger];
        bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
        if (trigger_accept) {
          for (unsigned int ifraction = 0; ifraction < num_jet_fractions; ifraction++) {
            //auto jet_fraction_func = jet_fraction_funcs_[ifraction];
            double jet_ef = jet_fraction_funcs_[ifraction](jet);
            if (jet_ieta > -1) {
              unsigned int offset = (1+itrigger)*num_eta_ranges*num_jet_fractions;
              ef_histograms_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_ef);
              pt_ef_histograms_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt, jet_ef);
              pt_ef_profiles_[offset + jet_ieta*num_jet_fractions + ifraction]->Fill(jet_pt, jet_ef);
            }

            if (jet_ipt > -1) {
              unsigned int offset = (1+itrigger)*num_pt_ranges*num_jet_fractions;
              eta_ef_histograms_[offset + jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_ef);
              eta_ef_profiles_[offset + jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_ef);
              eta_phi_ef_histograms_[offset + jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_phi, jet_ef);
              eta_phi_ef_profiles_[offset + jet_ipt*num_jet_fractions + ifraction]->Fill(jet_eta, jet_phi, jet_ef);
            }

            unsigned int offset = (1+itrigger)*num_jet_fractions;
            eta_pt_ef_histograms_[offset + ifraction]->Fill(jet_eta, jet_pt, jet_ef);
            eta_pt_ef_profiles_[offset + ifraction]->Fill(jet_eta, jet_pt, jet_ef);
            eta_phi_pt_ef_profiles_[offset + ifraction]->Fill(jet_eta, jet_phi, jet_pt, jet_ef);
          
          } // for jet_fractions
        } // if trigger_accept
      } // for triggers 
      //num_jets_filled++; // filled one more object
    }
    /*
    if ((num_jets_to_fill_ != -1) && (num_jets_filled >= num_jets_to_fill_)) {
      break;
    }
    */
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
