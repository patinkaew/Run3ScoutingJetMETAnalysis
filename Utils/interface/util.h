#ifndef Run3ScoutingJetMETAnalysis_Utils_util_h
#define Run3ScoutingJetMETAnalysis_Utils_util_h

// Trigger bits
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

// C++ include files
#include <string>
#include <vector>

namespace util {
  bool isAnyTriggerAccept(const std::vector<std::string> &triggers, const edm::TriggerResults &l1TriggerResults, const edm::TriggerResultsByName &hltTriggerResultsByName);

} // namespace util

#endif // Run3ScoutingJetMETAnalysis_Utils_util_h
