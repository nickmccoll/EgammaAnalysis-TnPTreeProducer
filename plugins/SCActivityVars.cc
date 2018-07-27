#include "SCActivityVars.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"


typedef SCActivityVars<reco::RecoEcalCandidate> SCActivityVarHelper;
DEFINE_FWK_MODULE(SCActivityVarHelper);
