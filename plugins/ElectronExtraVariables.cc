#include "ElectronExtraVariables.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/Framework/interface/MakerMacros.h"

typedef ElectronExtraVariables<pat::Electron> PatElectronExtraVariables;
DEFINE_FWK_MODULE(PatElectronExtraVariables);
