// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"

struct JetIdCriterion {
  double min_abs_eta;
  double max_abs_eta;
  double max_neutral_hadron_energy_fraction;
  double max_neutral_em_energy_fraction;
  double max_muon_energy_fraction;
  double min_charged_hadron_energy_fraction;
  double max_charged_em_energy_fraction;
  int min_num_constituents;
  int min_charged_multiplicity;
  int min_neutral_multiplicity;
};

template <typename JetType>
class JetIdProducer : public edm::stream::EDProducer<> {
public:
  using JetCollection = std::vector<JetType>; 

  JetIdProducer(const edm::ParameterSet&);
  ~JetIdProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  
  bool passJetId(JetType const &);

  edm::EDGetTokenT<JetCollection> jet_collection_token_;

  std::vector<JetIdCriterion> jet_id_criteria_;
};

template <typename JetType>
JetIdProducer<JetType>::JetIdProducer(const edm::ParameterSet& iConfig)
  : jet_collection_token_(consumes(iConfig.getParameter<edm::InputTag>("src"))) {
  std::vector<edm::ParameterSet> criteria = iConfig.getParameter<std::vector<edm::ParameterSet>>("criteria");
  for (auto const &criterion : criteria) {
    JetIdCriterion jet_id_criterion;
    jet_id_criterion.min_abs_eta = criterion.getParameter<double>("min_abs_eta");
    jet_id_criterion.max_abs_eta = criterion.getParameter<double>("max_abs_eta");
    jet_id_criterion.max_neutral_hadron_energy_fraction = criterion.getParameter<double>("max_neutral_hadron_energy_fraction");
    jet_id_criterion.max_neutral_em_energy_fraction = criterion.getParameter<double>("max_neutral_em_energy_fraction");
    jet_id_criterion.max_muon_energy_fraction = criterion.getParameter<double>("max_muon_energy_fraction");
    jet_id_criterion.min_charged_hadron_energy_fraction = criterion.getParameter<double>("min_charged_hadron_energy_fraction");
    jet_id_criterion.max_charged_em_energy_fraction = criterion.getParameter<double>("max_charged_em_energy_fraction");
    jet_id_criterion.min_num_constituents = criterion.getParameter<int>("min_num_constituents");
    jet_id_criterion.min_charged_multiplicity = criterion.getParameter<int>("min_charged_multiplicity");
    jet_id_criterion.min_neutral_multiplicity = criterion.getParameter<int>("min_neutral_multiplicity");
    jet_id_criteria_.push_back(jet_id_criterion);
  }

  produces<JetCollection>();
}

template <typename JetType>
void JetIdProducer<JetType>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto jet_collection_handle = iEvent.getHandle(jet_collection_token_);
  if (!jet_collection_handle.isValid()) {
    return;
  }
  auto jet_collection_good = std::make_unique<JetCollection>();
  for (auto const &jet : *jet_collection_handle) {
    if (passJetId(jet)) {
      jet_collection_good->push_back(jet);
    }
  }
  iEvent.put(std::move(jet_collection_good));
}

template <typename JetType>
bool passJetId(JetType const &jet) {
  return true;
}

template <typename JetType>
void JetIdProducer<JetType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  
  edm::ParameterSetDescription criterion_desc;
  criterion_desc.add<double>("min_abs_eta");
  criterion_desc.add<double>("max_abs_eta");
  criterion_desc.add<double>("max_neutral_hadron_energy_fraction", 99);
  criterion_desc.add<double>("max_neutral_em_energy_fraction", 99);
  criterion_desc.add<double>("max_muon_energy_fraction", 99);
  criterion_desc.add<double>("min_charged_hadron_energy_fraction", -1);
  criterion_desc.add<double>("max_charged_em_energy_fraction", -1);
  criterion_desc.add<int>("min_num_constituents", -1);
  criterion_desc.add<int>("min_charged_multiplicity", -1);
  criterion_desc.add<int>("min_neutral_multiplicity", -1);
  std::vector<edm::ParameterSet> criteria;
  edm::ParameterSet criterion;
  criterion.addParameter<double>("min_abs_eta", 0);
  criterion.addParameter<double>("max_abs_eta", 5);
  criterion.addParameter<double>("max_neutral_hadron_energy_fraction", 99);
  criterion.addParameter<double>("max_neutral_em_energy_fraction", 99);
  criterion.addParameter<double>("max_muon_energy_fraction", 99);
  criterion.addParameter<double>("min_charged_hadron_energy_fraction", -1);
  criterion.addParameter<double>("max_charged_em_energy_fraction", -1);
  criterion.addParameter<int>("min_num_constituents", -1);
  criterion.addParameter<int>("min_charged_multiplicity", -1);
  criterion.addParameter<int>("min_neutral_multiplicity", -1);
  criteria.push_back(criterion);
  
  desc.addVPSet("criteria", criterion_desc, criteria);
  descriptions.addWithDefaultLabel(desc);
}

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
using ScoutingPFJetIdProducer = JetIdProducer<Run3ScoutingPFJet>;
DEFINE_FWK_MODULE(ScoutingPFJetIdProducer);

template <>
bool ScoutingPFJetIdProducer::passJetId(Run3ScoutingPFJet const & scoutingPFJet) {
  float abs_eta = abs(scoutingPFJet.eta());

  int num_charged_hadrons = scoutingPFJet.chargedHadronMultiplicity();
  int num_neutral_hadrons = scoutingPFJet.neutralHadronMultiplicity();
  int num_photons = scoutingPFJet.photonMultiplicity();
  int num_electrons = scoutingPFJet.electronMultiplicity();
  int num_muons = scoutingPFJet.muonMultiplicity();
  int num_hf_hadrons = scoutingPFJet.HFHadronMultiplicity();
  int num_hf_ems = scoutingPFJet.HFEMMultiplicity();

  int num_charged = num_charged_hadrons + num_electrons + num_muons;
  int num_neutral = num_neutral_hadrons + num_photons + num_hf_hadrons + num_hf_ems;
  int num_constituents = num_charged + num_neutral;

  float charged_hadron_energy = scoutingPFJet.chargedHadronEnergy();
  float neutral_hadron_energy = scoutingPFJet.neutralHadronEnergy() + scoutingPFJet.HFHadronEnergy();
  float neutral_em_energy = scoutingPFJet.photonEnergy() + scoutingPFJet.HFEMEnergy();
  float charged_em_energy = scoutingPFJet.electronEnergy() + scoutingPFJet.muonEnergy();
  float energy = charged_hadron_energy + neutral_hadron_energy + charged_em_energy + neutral_em_energy;
  double charged_hadron_energy_fraction = charged_hadron_energy / energy;
  double neutral_hadron_energy_fraction = neutral_hadron_energy / energy;
  double charged_em_energy_fraction = charged_em_energy / energy;
  double neutral_em_energy_fraction = neutral_em_energy / energy;
  double muon_energy_fraction = scoutingPFJet.muonEnergy() / energy;

  // this is from offline CHS
  
  for (auto const &jet_id_criterion : jet_id_criteria_) {
    if (abs_eta > jet_id_criterion.min_abs_eta
        && abs_eta <= jet_id_criterion.max_abs_eta
        && neutral_hadron_energy_fraction < jet_id_criterion.max_neutral_hadron_energy_fraction
        && neutral_em_energy_fraction < jet_id_criterion.max_neutral_hadron_energy_fraction
        && num_constituents > jet_id_criterion.min_num_constituents
        && muon_energy_fraction < jet_id_criterion.max_muon_energy_fraction
        && charged_hadron_energy_fraction > jet_id_criterion.min_charged_hadron_energy_fraction
        && num_charged > jet_id_criterion.min_charged_multiplicity
        && charged_em_energy_fraction < jet_id_criterion.max_charged_em_energy_fraction
        && num_neutral >= jet_id_criterion.min_neutral_multiplicity
       ) {
      return true;
    }
  }
  return false;
}

#include "DataFormats/PatCandidates/interface/Jet.h"
using PATJetIdProducer = JetIdProducer<pat::Jet>;
DEFINE_FWK_MODULE(PATJetIdProducer);

template<>
bool PATJetIdProducer::passJetId(pat::Jet const &patJet) {
  auto abs_eta = abs(patJet.eta());
  for (auto const &jet_id_criterion : jet_id_criteria_) {
    if (abs_eta > jet_id_criterion.min_abs_eta
        && abs_eta <= jet_id_criterion.max_abs_eta
        && patJet.neutralHadronEnergyFraction() < jet_id_criterion.max_neutral_hadron_energy_fraction
        && patJet.neutralEmEnergyFraction() < jet_id_criterion.max_neutral_em_energy_fraction
        && patJet.chargedMultiplicity()+patJet.neutralMultiplicity() > jet_id_criterion.min_num_constituents
        && patJet.muonEnergyFraction() < jet_id_criterion.max_muon_energy_fraction
        && patJet.chargedHadronEnergyFraction() > jet_id_criterion.min_charged_hadron_energy_fraction
        && patJet.chargedMultiplicity() > jet_id_criterion.min_charged_multiplicity
        && patJet.chargedEmEnergyFraction() < jet_id_criterion.max_charged_em_energy_fraction
        && patJet.neutralMultiplicity() >= jet_id_criterion.min_neutral_multiplicity
       ) {
      return true;
    }
  }
  return false;
}

#include "DataFormats/JetReco/interface/PFJet.h"
using RecoPFJetIdProducer = JetIdProducer<reco::PFJet>;
DEFINE_FWK_MODULE(RecoPFJetIdProducer);

template<>
bool RecoPFJetIdProducer::passJetId(reco::PFJet const &recoPFJet) {
  auto abs_eta = abs(recoPFJet.eta());
  for (auto const &jet_id_criterion : jet_id_criteria_) {
    if (abs_eta > jet_id_criterion.min_abs_eta
        && abs_eta <= jet_id_criterion.max_abs_eta
        && recoPFJet.neutralHadronEnergyFraction() < jet_id_criterion.max_neutral_hadron_energy_fraction
        && recoPFJet.neutralEmEnergyFraction() < jet_id_criterion.max_neutral_em_energy_fraction
        && recoPFJet.chargedMultiplicity()+recoPFJet.neutralMultiplicity() > jet_id_criterion.min_num_constituents
        && recoPFJet.muonEnergyFraction() < jet_id_criterion.max_muon_energy_fraction
        && recoPFJet.chargedHadronEnergyFraction() > jet_id_criterion.min_charged_hadron_energy_fraction
        && recoPFJet.chargedMultiplicity() > jet_id_criterion.min_charged_multiplicity
        && recoPFJet.chargedEmEnergyFraction() < jet_id_criterion.max_charged_em_energy_fraction
        && recoPFJet.neutralMultiplicity() >= jet_id_criterion.min_neutral_multiplicity
       ) {
      return true;
    }
  }
  return false;
}
