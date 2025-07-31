#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

class Run3ScoutingPFJetToRecoPFJetProducer : public edm::stream::EDProducer<> {
  public:
    explicit Run3ScoutingPFJetToRecoPFJetProducer(edm::ParameterSet const& params);
    ~Run3ScoutingPFJetToRecoPFJetProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

    reco::PFJet createPFJet(Run3ScoutingPFJet const& scoutingPFJet);

    edm::EDGetTokenT<Run3ScoutingPFJetCollection> scoutingPFJet_collection_token_;
};

Run3ScoutingPFJetToRecoPFJetProducer::Run3ScoutingPFJetToRecoPFJetProducer(edm::ParameterSet const& params)
  : scoutingPFJet_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingPFJet"))) {
  produces<reco::PFJetCollection>();
}

void Run3ScoutingPFJetToRecoPFJetProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
    // produce reco::PFJet
    edm::Handle<Run3ScoutingPFJetCollection> scoutingPFJet_collection_handle = iEvent.getHandle(scoutingPFJet_collection_token_);
    auto recoPFJet_collection_ptr = std::make_unique<reco::PFJetCollection>();
    //auto scoutingPFJetRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingPFJet>>();
    if (scoutingPFJet_collection_handle.isValid()) {
        for (size_t scoutingPFJet_index = 0; scoutingPFJet_index < scoutingPFJet_collection_handle->size(); scoutingPFJet_index++) {
            auto &scoutingPFJet = scoutingPFJet_collection_handle->at(scoutingPFJet_index);
            recoPFJet_collection_ptr->push_back(createPFJet(scoutingPFJet));
            //scoutingPFJetRef_collection_ptr->push_back(edm::Ref<Run3ScoutingPFJetCollection>(scoutingPFJet_collection_handle, scoutingPFJet_index));
        }
    }

    iEvent.put(std::move(recoPFJet_collection_ptr));
}

reco::PFJet Run3ScoutingPFJetToRecoPFJetProducer::createPFJet(Run3ScoutingPFJet const& scoutingPFJet) {
    // fill LorentzVector P4
    float px = scoutingPFJet.pt() * cos(scoutingPFJet.phi());
    float py = scoutingPFJet.pt() * sin(scoutingPFJet.phi());
    float pz = scoutingPFJet.pt() * sinh(scoutingPFJet.eta());
    float p = scoutingPFJet.pt() * cosh(scoutingPFJet.eta());
    float energy = std::sqrt(p * p + scoutingPFJet.m() * scoutingPFJet.m());
    reco::Particle::LorentzVector p4(px, py, pz, energy);
    
    // fill vertex with default (0, 0, 0)
    reco::Particle::Point vertex(0, 0, 0);
    
    // fill specific
    reco::PFJet::Specific specific;
    specific.mChargedHadronEnergy = scoutingPFJet.chargedHadronEnergy();
    specific.mNeutralHadronEnergy = scoutingPFJet.neutralHadronEnergy();
    specific.mPhotonEnergy = scoutingPFJet.photonEnergy();
    specific.mElectronEnergy = scoutingPFJet.electronEnergy();
    specific.mMuonEnergy = scoutingPFJet.muonEnergy();
    specific.mHFHadronEnergy = scoutingPFJet.HFHadronEnergy();
    specific.mHFEMEnergy = scoutingPFJet.HFEMEnergy();
    
    specific.mChargedHadronMultiplicity = scoutingPFJet.chargedHadronMultiplicity();
    specific.mNeutralHadronMultiplicity = scoutingPFJet.neutralHadronMultiplicity();
    specific.mPhotonMultiplicity = scoutingPFJet.photonMultiplicity();
    specific.mElectronMultiplicity = scoutingPFJet.electronMultiplicity();
    specific.mMuonMultiplicity = scoutingPFJet.muonMultiplicity();
    specific.mHFHadronMultiplicity = scoutingPFJet.HFHadronMultiplicity();
    specific.mHFEMMultiplicity = scoutingPFJet.HFEMMultiplicity();

    specific.mChargedEmEnergy = scoutingPFJet.electronEnergy();
    specific.mChargedMuEnergy = scoutingPFJet.muonEnergy();
    specific.mNeutralEmEnergy = scoutingPFJet.photonEnergy() + + scoutingPFJet.HFEMEnergy();
    specific.mChargedMultiplicity = scoutingPFJet.chargedHadronMultiplicity() + scoutingPFJet.electronMultiplicity() + scoutingPFJet.muonMultiplicity();
    specific.mNeutralMultiplicity = scoutingPFJet.neutralHadronMultiplicity() + scoutingPFJet.photonMultiplicity() + scoutingPFJet.HFHadronMultiplicity() + scoutingPFJet.HFEMMultiplicity();

    specific.mHOEnergy = scoutingPFJet.HOEnergy();
    
    // create reco::PFJet output
    reco::PFJet recoPFJet(p4, vertex, specific);

    // set jetArea
    recoPFJet.setJetArea(scoutingPFJet.jetArea());
    
    return recoPFJet;
}

void Run3ScoutingPFJetToRecoPFJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("scoutingPFJet", edm::InputTag("hltScoutingPFPacker"));

    descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Run3ScoutingPFJetToRecoPFJetProducer);
