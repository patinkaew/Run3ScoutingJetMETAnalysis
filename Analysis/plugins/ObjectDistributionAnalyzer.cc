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

// C++ include files
#include <iostream>
#include <cmath>
#include <boost/algorithm/string/join.hpp>

#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

template <typename ObjectType>
class ObjectDistributionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using ObjectCollection = std::vector<ObjectType>; 

  ObjectDistributionAnalyzer(const edm::ParameterSet&);
  ~ObjectDistributionAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};

  // input tokens
  edm::EDGetTokenT<ObjectCollection> object_collection_token_;

  // trigger tokens
  edm::EDGetTokenT<edm::TriggerResults> l1TriggerResults_token_;
  edm::EDGetTokenT<edm::TriggerResults> hltTriggerResults_token_;

  // file IO
  edm::Service<TFileService> fs_;

  // functor
  StringCutObjectSelector<ObjectType> cut_; // general cut applied to all object
  std::vector<StringCutObjectSelector<ObjectType>> category_cuts_; // cut to each category
  const StringObjectFunction<ObjectType> object_pt_func_;
  const StringObjectFunction<ObjectType> object_eta_func_;
  const StringObjectFunction<ObjectType> object_phi_func_;

  // triggers
  std::vector<std::vector<std::string>> triggers_; // trigger expressions

  // other parameters
  int num_objects_to_fill_;

  // histograms
  std::vector<TH2D*> eta_phi_histograms_;
  std::vector<TH2D*> eta_pt_histograms_;
};

template <typename ObjectType>
ObjectDistributionAnalyzer<ObjectType>::ObjectDistributionAnalyzer(const edm::ParameterSet& iConfig)
  : object_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))),
    l1TriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("L1TriggerResults"))),
    hltTriggerResults_token_(consumes(iConfig.getParameter<edm::InputTag>("HLTTriggerResults"))),
    cut_(iConfig.getParameter<std::string>("cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    object_pt_func_(iConfig.getParameter<std::string>("pt_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    object_eta_func_(iConfig.getParameter<std::string>("eta_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    object_phi_func_(iConfig.getParameter<std::string>("phi_func"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    num_objects_to_fill_(iConfig.getParameter<int>("num_objects_to_fill")) {
  
  // parse categories
  bool lazy_eval = iConfig.getUntrackedParameter<bool>("lazy_eval");
  std::vector<std::string> category_names = {};
  category_cuts_ = {};
  if (iConfig.existsAs<std::vector<edm::ParameterSet>>("categories")) {
    auto category_vpset = iConfig.getParameter<std::vector<edm::ParameterSet>>("categories");
    for (auto const &category_pset : category_vpset) {
      std::string category_cut = category_pset.getParameter<std::string>("cut");
      category_cuts_.push_back(StringCutObjectSelector<ObjectType>(category_cut, lazy_eval));
      std::string category_name = category_pset.getParameter<std::string>("name");
      if (!category_name.empty()) {
        category_names.push_back(category_name);
      } else {
        category_names.push_back(category_cut);
      }
    }
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
  
  // initialize directories and histograms
  if (category_cuts_.size() == 0) { // no categorization
    for (auto const &trigger_name : trigger_names) {
      TFileDirectory trigger_dir = fs_->mkdir(trigger_name.c_str());

      eta_phi_histograms_.push_back(trigger_dir.make<TH2D>("eta-phi", ";#eta; #phi; Events", 50, -3., 3., 50, -M_PI, M_PI));
      eta_pt_histograms_.push_back(trigger_dir.make<TH2D>("eta-pt", ";#eta; pT (GeV); Events", 50, -3., 3., 50, 0., 500.));
    }
  } else {
    for (auto const &category_name : category_names) {
      TFileDirectory category_dir = fs_->mkdir(category_name.c_str());
      for (auto const &trigger_name : trigger_names) {
        TFileDirectory trigger_dir = category_dir.mkdir(trigger_name.c_str());
        eta_phi_histograms_.push_back(trigger_dir.make<TH2D>("eta-phi", ";#eta; #phi; Events", 50, -3., 3., 50, -M_PI, M_PI));
        eta_pt_histograms_.push_back(trigger_dir.make<TH2D>("eta-pt", ";#eta; pT (GeV); Events", 50, -3., 3., 50, 0., 500.));
      }
    }
  } 
}

template <typename ObjectType>
void ObjectDistributionAnalyzer<ObjectType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<edm::InputTag>("L1TriggerResults", edm::InputTag("l1bits"));
  desc.add<edm::InputTag>("HLTTriggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<int>("num_objects_to_fill", -1); // -1 for all objects
  desc.addUntracked<bool>("lazy_eval", false);
  desc.add<std::string>("cut", "");

  desc.add<std::string>("pt_func", "pt()");
  desc.add<std::string>("eta_func", "eta()");
  desc.add<std::string>("phi_func", "phi()");

  edm::ParameterSetDescription category_desc;
  category_desc.add<std::string>("name", "");
  category_desc.add<std::string>("cut");
  edm::ParameterSet category_pset;
  category_pset.addParameter<std::string>("name", "");
  category_pset.addParameter<std::string>("cut", "");
  std::vector<edm::ParameterSet> category_vpset;
  category_vpset.push_back(category_pset);
  desc.addVPSetOptional("categories", category_desc, category_vpset);

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

template <typename ObjectType>
void ObjectDistributionAnalyzer<ObjectType>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  auto object_collection_handle = iEvent.getHandle(object_collection_token_);
  if (!object_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Object is invalid";
    return;
  }
  
  int num_objects_filled = 0;
  for (unsigned int iobj = 0; iobj < object_collection_handle->size(); iobj++) {
    auto object = object_collection_handle->at(iobj);
    if (cut_(object)) {
      if (category_cuts_.size() == 0) { // no categorization
        // no trigger
        eta_phi_histograms_[0]->Fill(object_eta_func_(object), object_phi_func_(object));
        eta_pt_histograms_[0]->Fill(object_eta_func_(object), object_pt_func_(object));
        // with trigger
        for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
          auto trigger = triggers_[itrigger];
          bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
          if (trigger_accept) {
            eta_phi_histograms_[itrigger+1]->Fill(object_eta_func_(object), object_phi_func_(object));
            eta_pt_histograms_[itrigger+1]->Fill(object_eta_func_(object), object_pt_func_(object));
          }
        }
      } else { // categorization
        for (unsigned int icategory = 0; icategory < category_cuts_.size(); icategory++) {
          unsigned int offset = icategory * (triggers_.size() + 1); // each category has (num_trigger+1) subdirectories
          if (category_cuts_[icategory](object)) {
            // no trigger
            eta_phi_histograms_[offset]->Fill(object_eta_func_(object), object_phi_func_(object));
            eta_pt_histograms_[offset]->Fill(object_eta_func_(object), object_pt_func_(object));
            // with trigger
            for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
              auto trigger = triggers_[itrigger];
              bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
              if (trigger_accept) {
                eta_phi_histograms_[offset+itrigger+1]->Fill(object_eta_func_(object), object_phi_func_(object));
                eta_pt_histograms_[offset+itrigger+1]->Fill(object_eta_func_(object), object_pt_func_(object));
              }
            }
          }
        }
      }
      num_objects_filled++; // filled one more object
    }
    if ((num_objects_to_fill_ != -1) && (num_objects_filled >= num_objects_to_fill_)) {
      break;
    }
  }
}

#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
using ScoutingPFCandidateDistributionAnalyzer = ObjectDistributionAnalyzer<Run3ScoutingParticle>;
DEFINE_FWK_MODULE(ScoutingPFCandidateDistributionAnalyzer);

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
using PatPFCandidateDistributionAnalyzer = ObjectDistributionAnalyzer<pat::PackedCandidate>;
DEFINE_FWK_MODULE(PatPFCandidateDistributionAnalyzer);
