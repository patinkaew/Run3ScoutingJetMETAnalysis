// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysis/Run3ScoutingJetMETAnalysis
// Class:      MatchingEfficiencyAnalyzer
//
/**\class MatchingEfficiencyAnalyzer MatchingEfficiencyAnalyzer.cc Run3ScoutingAnalysis/Run3ScoutingJetMETAnalysis/plugins/MatchingEfficiencyAnalyzer.cc

 Description: Fill histograms necessary for matching efficiency calculation

 Implementation:
   This analyzer fill four histograms necessary for calculating matching efficiency
     1) object1
     2) object2
     3) object1-to-object2
     4) object2-to-object1
   Then, for example, the matching efficiency of object1-to-object2 is calculated as object1-to-object2 / object1

   The matching efficiency can be used both for efficiency and purity
   For example:
     If object1 is offline and object2 is scouting,
     Efficiency of Scouting is matching efficiency of object1-to-object2 (fraction of offline objects with scouting counterpart)
     Purity of Scouting is matching efficiency of object2-to-object1 (fraction of scouting objects with offline counterpart)
*/
//
// Original Author:  Patin Inkaew

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

// object selection include files
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

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

// ROOT include files
#include "TH1.h"
#include "TH2.h"

// C++ include files
#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <fmt/core.h>

#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

template <typename ObjType1, typename ObjType2>
class MatchingEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  using Obj1Collection = std::vector<ObjType1>;
  using Obj2Collection = std::vector<ObjType2>;

  MatchingEfficiencyAnalyzer(const edm::ParameterSet&);
  ~MatchingEfficiencyAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {};
  
  // define criteria whether object1 is matched with object2
  // currently use symmetric criteria
  // this is the only function ones need to define
  bool isMatched(ObjType1 const &, ObjType2 const &);
  
  // for given object1 (object2), does it have a matched in collection of object2 (object1)
  bool isObj1Matched(ObjType1 const &object1, Obj2Collection const &object2_collection, std::set<unsigned int> &already_matched_indices);
  bool isObj2Matched(ObjType2 const &object2, Obj1Collection const &object1_collection, std::set<unsigned int> &already_matched_indices);

  // input tokens
  edm::EDGetTokenT<Obj1Collection> obj1_collection_token_;
  edm::EDGetTokenT<Obj2Collection> obj2_collection_token_;

  // trigger token
  edm::EDGetTokenT<edm::TriggerResults> l1TriggerResults_token_;
  edm::EDGetTokenT<edm::TriggerResults> hltTriggerResults_token_;

  // file IO
  edm::Service<TFileService> fs_;
  
  // triggers
  std::vector<std::vector<std::string>> triggers_; // trigger expressions
  
  // number of objects to match, i.e. if this is 1, study matching efficiency of leading, if 2, study matching efficiency of leading and subleading
  unsigned int num_objects_to_match_;

  const StringCutObjectSelector<ObjType1> object1_cut_;
  const StringCutObjectSelector<ObjType2> object2_cut_;
  
  std::vector<double> eta_ranges_;

  // histograms
  std::vector<TH1D*> pt_histograms;
  //std::vector<TH1D*> eta_histograms;
  //std::vector<TH2D*> pt_eta_histograms;
};

template <typename ObjType1, typename ObjType2>
MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::MatchingEfficiencyAnalyzer(const edm::ParameterSet& iConfig)
  : obj1_collection_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("object1"))),
    obj2_collection_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("object2"))),
    l1TriggerResults_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("L1TriggerResults"))),
    hltTriggerResults_token_(consumes(iConfig.getUntrackedParameter<edm::InputTag>("HLTTriggerResults"))),
    object1_cut_(iConfig.getUntrackedParameter<std::string>("object1_cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")),
    object2_cut_(iConfig.getUntrackedParameter<std::string>("object2_cut"), iConfig.getUntrackedParameter<bool>("lazy_eval")) {
  
  std::string object1_name = iConfig.getUntrackedParameter<std::string>("object1_name");
  if (object1_name.empty()) object1_name = "Object1";
  std::string object2_name = iConfig.getUntrackedParameter<std::string>("object2_name");
  if (object2_name.empty()) object2_name = "Object2";
  
  num_objects_to_match_ = iConfig.getUntrackedParameter<unsigned int>("num_objects_to_match");
  
  /*
  bool lazy_eval = iConfig.getUntrackedParameter<bool>("lazy_eval");
  object1_cut_ = StringCutObjectSelector<ObjType1>(iConfig.getUntrackedParameter<std::string>("object1_cut"), lazy_eval);
  object2_cut_ = StringCutObjectSelector<ObjType2>(iConfig.getUntrackedParameter<std::string>("object2_cut"), lazy_eval);
  */
  
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
  
  // prepare directories in files and histograms
  eta_ranges_ = {0.0, 0.5}; //{0.0, 0.5, 1.0, 1.5, 2.0, 2.6, 2.7, 3.0, 5.0};

  for (auto const &trigger_name : trigger_names) {
    TFileDirectory trigger_dir = fs_->mkdir(trigger_name.c_str());
    for (size_t iobject_to_match = 0; iobject_to_match < num_objects_to_match_; iobject_to_match++) {
      TFileDirectory object_to_match_dir = trigger_dir.mkdir((std::string("Match")+std::to_string(iobject_to_match)).c_str());
      for (size_t ieta = 0; ieta < eta_ranges_.size()-1; ieta++) {
        TFileDirectory eta_range_dir = object_to_match_dir.mkdir((fmt::format("{:.2f}",eta_ranges_[ieta])+"<|eta|<"+fmt::format("{:.2f}",eta_ranges_[ieta+1])).c_str());
        pt_histograms.push_back(eta_range_dir.make<TH1D>(object1_name.c_str(), ";pT (GeV); Events", 30, 0., 3000.));
        pt_histograms.push_back(eta_range_dir.make<TH1D>((object1_name+"-to-"+object2_name).c_str(), ";pT (GeV); Events", 30, 0., 3000.));
        pt_histograms.push_back(eta_range_dir.make<TH1D>(object2_name.c_str(), ";pT (GeV); Events", 30, 0., 3000.));
        pt_histograms.push_back(eta_range_dir.make<TH1D>((object2_name+"-to-"+object1_name).c_str(), ";pT (GeV); Events", 30, 0., 3000.));
      }
    }
  }
}

template <typename ObjType1, typename ObjType2>
void MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("object1");
  desc.addUntracked<edm::InputTag>("object2");
  desc.addUntracked<edm::InputTag>("L1TriggerResults", edm::InputTag("l1bits"));
  desc.addUntracked<edm::InputTag>("HLTTriggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.addUntracked<std::string>("object1_name", "");
  desc.addUntracked<std::string>("object2_name", "");
  desc.addUntracked<unsigned int>("num_objects_to_match", 1);
  desc.addUntracked<bool>("lazy_eval", false);
  desc.addUntracked<std::string>("object1_cut", "");
  desc.addUntracked<std::string>("object2_cut", "");

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

template <typename ObjType1, typename ObjType2>
void MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  auto obj1_collection_handle = iEvent.getHandle(obj1_collection_token_);
  if (!obj1_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Object1 is invalid";
    return;
  }

  auto obj2_collection_handle = iEvent.getHandle(obj2_collection_token_);
  if (!obj2_collection_handle.isValid()) {
    edm::LogWarning ("Handle") << "Object2 is invalid";
    return;
  }

  // select object1
  auto obj1_collection_good_ptr = std::make_unique<Obj1Collection>();
  for (const auto &object1 : *obj1_collection_handle) {
    if (object1_cut_(object1)) {
      obj1_collection_good_ptr->push_back(object1);
    }
  }

  // select object2
  auto obj2_collection_good_ptr = std::make_unique<Obj2Collection>();
  for (const auto &object2 : *obj2_collection_handle) {
    if (object2_cut_(object2)) {
      obj2_collection_good_ptr->push_back(object2);
    }
  }

  // build matching flags for object1
  std::vector<bool> obj1_match_flags(num_objects_to_match_, false);
  std::set<unsigned int> obj2_already_matched_indices;
  for (unsigned int iobj1 = 0; iobj1 < num_objects_to_match_; iobj1++) {
    if (iobj1 < obj1_collection_good_ptr->size()) {
      auto object1 = obj1_collection_good_ptr->at(iobj1);
      obj1_match_flags[iobj1] = isObj1Matched(object1, *obj2_collection_good_ptr, obj2_already_matched_indices);
    } else {
      // if number of objects are less than requested number to match, set flag to false
      // the vector is already initialized with false filled, so we don't need to do anything
      // we do it anyway for clarity
      obj1_match_flags[iobj1] = false;
    }
  }

  // build matching flags for object2
  std::vector<bool> obj2_match_flags(num_objects_to_match_, false);
  std::set<unsigned int> obj1_already_matched_indices;
  for (unsigned int iobj2 = 0; iobj2 < num_objects_to_match_; iobj2++) {
    if (iobj2 < obj2_collection_good_ptr->size()) {
      auto object2 = obj2_collection_good_ptr->at(iobj2);
      obj2_match_flags[iobj2] = isObj2Matched(object2, *obj1_collection_good_ptr, obj1_already_matched_indices);
    } else {
      // if number of objects are less than requested number to match, set flag to false
      // the vector is already initialized with false filled, so we don't need to do anything
      // we do it anyway for clarity
      obj2_match_flags[iobj2] = false;
    }
  }

  // fill histograms  
  size_t num_eta_ranges = eta_ranges_.size() - 1;
  for (unsigned int iobj = 0; iobj < num_objects_to_match_; iobj++) {
    bool has_object1 = iobj < obj1_collection_good_ptr->size();
    bool has_object2 = iobj < obj2_collection_good_ptr->size();

    double object1_pt = -1;
    double object1_abs_eta = -1;
    double object2_pt = -1;
    double object2_abs_eta = -1;

    // eta bin index
    int object1_ieta = -1;
    int object2_ieta = -1;

    if (has_object1) {
      auto object1 = obj1_collection_good_ptr->at(iobj);
      object1_pt = object1.pt();
      object1_abs_eta = abs(object1.eta());
      for (unsigned int ieta = 0; ieta < num_eta_ranges; ieta++) {
        if (object1_abs_eta > eta_ranges_[ieta] && object1_abs_eta <= eta_ranges_[ieta+1]) {
          object1_ieta = ieta;
          break;
        }
      }
      if (object1_ieta == -1) has_object1 = false; // out of eta range
    }

    if (has_object2) {
      auto object2 = obj2_collection_good_ptr->at(iobj);
      object2_pt = object2.pt();
      object2_abs_eta = abs(object2.eta());
      for (unsigned int ieta = 0; ieta < num_eta_ranges; ieta++) {
        if (object2_abs_eta > eta_ranges_[ieta] && object2_abs_eta <= eta_ranges_[ieta+1]) {
          object2_ieta = ieta;
          break;
        }
      }
      if (object2_ieta == -1) has_object2 = false; // out of eta range
    }

    // no trigger
    unsigned int offset = iobj*num_eta_ranges*4;
    if (has_object1) {
      pt_histograms[offset + object1_ieta*4]->Fill(object1_pt);
      if (obj1_match_flags[iobj]) {
        pt_histograms[offset + object1_ieta*4 + 1]->Fill(object1_pt);
      }
    }
    if (has_object2) {
      pt_histograms[offset + object2_ieta*4 + 2]->Fill(object2_pt);
      if (obj2_match_flags[iobj]) {
        pt_histograms[offset + object2_ieta*4 + 3]->Fill(object2_pt);
      }
    }
    // with trigger
    for (unsigned int itrigger = 0; itrigger < triggers_.size(); itrigger++) {
      auto trigger = triggers_[itrigger];
      bool trigger_accept = util::isAnyTriggerAccept(trigger, *l1TriggerResults_handle, hltTriggerResultsByName);
      offset = (itrigger+1)*num_objects_to_match_*num_eta_ranges*4 + iobj*num_eta_ranges*4;
      if (trigger_accept) {
        if (has_object1) {
          pt_histograms[offset + object1_ieta*4]->Fill(object1_pt);
          if (obj1_match_flags[iobj]) {
            pt_histograms[offset + object1_ieta*4 + 1]->Fill(object1_pt);
          }
        }
        if (has_object2) {
          pt_histograms[offset + object2_ieta*4 + 2]->Fill(object2_pt);
          if (obj2_match_flags[iobj]) {
            pt_histograms[offset + object2_ieta*4 + 3]->Fill(object2_pt);
          }
        }
      } // if (trigger_accept)
    } // for loop over triggers
  } // object to match
}

template <typename ObjType1, typename ObjType2>
bool MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::isMatched(ObjType1 const &object1, ObjType2 const& object2) {
  return true;
}

template <typename ObjType1, typename ObjType2>
bool MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::isObj1Matched(ObjType1 const &object1, Obj2Collection const &object2_collection, std::set<unsigned int> &already_matched_indices) {
  for (size_t iobj2 = 0; iobj2 < object2_collection.size(); iobj2++) {
    if (!already_matched_indices.empty() && already_matched_indices.find(iobj2) != already_matched_indices.end()) continue;
    ObjType2 object2 = object2_collection[iobj2];
    if (isMatched(object1, object2)) {
      already_matched_indices.insert(iobj2);
      return true;
    }
  }
  return false;
}

template <typename ObjType1, typename ObjType2>
bool MatchingEfficiencyAnalyzer<ObjType1, ObjType2>::isObj2Matched(ObjType2 const &object2, Obj1Collection const &object1_collection, std::set<unsigned int> &already_matched_indices) {
  for (size_t iobj1 = 0; iobj1 < object1_collection.size(); iobj1++) {
    if (!already_matched_indices.empty() && already_matched_indices.find(iobj1) != already_matched_indices.end()) continue;
    ObjType1 object1 = object1_collection[iobj1];
    if (isMatched(object1, object2)) {
      already_matched_indices.insert(iobj1);
      return true;
    }
  }
  return false;
}

// generic match functions
template <typename ObjType1, typename ObjType2>
bool isDeltaRMatched(ObjType1 object1, ObjType2 object2, double dr_max = 0.2) {
  return deltaR2<ObjType1, ObjType2>(object1, object2) <= dr_max;
}

// specialization
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

using ScoutingPFJetToPATJetMatchingEfficiencyAnalyzer = MatchingEfficiencyAnalyzer<Run3ScoutingPFJet, pat::Jet>;
DEFINE_FWK_MODULE(ScoutingPFJetToPATJetMatchingEfficiencyAnalyzer);

template <>
bool ScoutingPFJetToPATJetMatchingEfficiencyAnalyzer::isMatched(Run3ScoutingPFJet const &scoutingPFJet, pat::Jet const &patJet) {
  return isDeltaRMatched<Run3ScoutingPFJet, pat::Jet>(scoutingPFJet, patJet, 0.2);
}

using ScoutingPFJetToRecoPFJetMatchingEfficiencyAnalyzer = MatchingEfficiencyAnalyzer<Run3ScoutingPFJet, reco::PFJet>;
DEFINE_FWK_MODULE(ScoutingPFJetToRecoPFJetMatchingEfficiencyAnalyzer);

template <>
bool ScoutingPFJetToRecoPFJetMatchingEfficiencyAnalyzer::isMatched(Run3ScoutingPFJet const &scoutingPFJet, reco::PFJet const &recoPFJet) {
  return isDeltaRMatched<Run3ScoutingPFJet, reco::PFJet>(scoutingPFJet, recoPFJet, 0.2);
}
