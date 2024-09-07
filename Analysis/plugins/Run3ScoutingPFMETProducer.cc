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

#include "AnalysisDataFormats/Scouting/interface/Run3ScoutingPFMET.h"

class Run3ScoutingPFMETProducer : public edm::stream::EDProducer<> {
public:
  Run3ScoutingPFMETProducer(const edm::ParameterSet&);
  ~Run3ScoutingPFMETProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<double> pt_token_;
  edm::EDGetTokenT<double> phi_token_;
};

Run3ScoutingPFMETProducer::Run3ScoutingPFMETProducer(const edm::ParameterSet& iConfig) :
  pt_token_(consumes(iConfig.getParameter<edm::InputTag>("pt"))),
  phi_token_(consumes(iConfig.getParameter<edm::InputTag>("phi"))) {

  produces<std::vector<Run3ScoutingPFMET>>(); // need to make this a vector similar to pat::MET in MINIAOD
}

void Run3ScoutingPFMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pt");
  desc.add<edm::InputTag>("phi");

  descriptions.addWithDefaultLabel(desc);
}

void Run3ScoutingPFMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto pt_handle = iEvent.getHandle(pt_token_);
  auto phi_handle = iEvent.getHandle(phi_token_);
  auto pfmet_collection = std::make_unique<std::vector<Run3ScoutingPFMET>>();
  if (pt_handle.isValid() && phi_handle.isValid()) {
    pfmet_collection->push_back(Run3ScoutingPFMET(*pt_handle, *phi_handle));
  }
  iEvent.put(std::move(pfmet_collection));
}

DEFINE_FWK_MODULE(Run3ScoutingPFMETProducer);
