#include "Run3ScoutingJetMETAnalysis/Utils/interface/util.h"

bool util::isAnyTriggerAccept(const std::vector<std::string> &triggers, const edm::TriggerResults &l1TriggerResults, const edm::TriggerResultsByName &hltTriggerResultsByName) {
  for (auto const &trigger: triggers) {
    if (trigger.compare(0, 2, "L1") == 0) {
      auto trigger_names = l1TriggerResults.getTriggerNames();
      for (size_t itrigger = 0; itrigger < trigger_names.size(); itrigger++) {
        std::string trigger_name = trigger_names[itrigger];
        if (trigger_name.compare(0, trigger.length(), trigger) == 0) {
          if (l1TriggerResults.accept(itrigger)) return true;
        }
      }
    } else {
      for (size_t itrigger = 0; itrigger < hltTriggerResultsByName.size(); itrigger++) {
        std::string trigger_name = hltTriggerResultsByName.triggerName(itrigger);
        if (trigger_name.compare(0, trigger.length(), trigger) == 0) { // ignore _v*
          if (hltTriggerResultsByName.accept(itrigger)) return true;
        }
      }
    }
  }
  return false;
}
