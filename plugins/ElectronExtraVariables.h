#ifndef _ELECTRONVARIABLEHELPER_H
#define _ELECTRONVARIABLEHELPER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaR.h"

template <class T>
class ElectronExtraVariables : public edm::EDProducer {
 public:
  explicit ElectronExtraVariables(const edm::ParameterSet & iConfig);
  virtual ~ElectronExtraVariables() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  
private:
  edm::EDGetTokenT<std::vector<T> > probesToken_;
//  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection>     token_pfCand     ;
  const edm::EDGetTokenT<double>                             token_rho      ;
};

namespace Isolations{

double electronEA(double eta);
double muonEA(double eta);

double getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
              const reco::Candidate* ptcl,
              double r_iso_min, double r_iso_max, double kt_scale,
              bool charged_only, bool use_EA_corr, double EA, double rho);
}


template<class T>
ElectronExtraVariables<T>::ElectronExtraVariables(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
//  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  token_pfCand      (consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfParticles"))),
  token_rho       (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))){

  produces<edm::ValueMap<float> >("sip3D");
  produces<edm::ValueMap<float> >("miniiso");
}

template<class T>
ElectronExtraVariables<T>::~ElectronExtraVariables()
{}

template<class T>
void ElectronExtraVariables<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  edm::Handle<std::vector<T> > probes;
  edm::Handle<reco::VertexCollection> vtxH;
  edm::Handle<double>                             han_rho      ;
  edm::Handle<pat::PackedCandidateCollection>     han_pfCand     ;
  
  iEvent.getByToken(probesToken_, probes);
//  iEvent.getByToken(vtxToken_, vtxH);
  iEvent.getByToken(token_rho    ,han_rho    );
  iEvent.getByToken(token_pfCand   ,han_pfCand   );

  // prepare vector for output
  std::vector<float> sip3DVals;sip3DVals.reserve(probes->size());
  std::vector<float> miniIsoVals;miniIsoVals.reserve(probes->size());

  for(const auto& probe: *probes){
      const float sip3d=std::fabs(probe.dB(pat::Electron::PV3D) / probe.edB(pat::Electron::PV3D));
      const float eA = Isolations::electronEA(probe.superCluster()->eta());
      const float eISO = Isolations::getPFMiniIsolation(han_pfCand, dynamic_cast<const reco::Candidate *>
      (&probe), 0.05, 0.2, 10., false, true, eA, *han_rho);
      sip3DVals.push_back(sip3d);
      miniIsoVals.push_back(eISO);

  }

  auto fill = [&](const std::string& name, const std::vector<float>& vals) {
      // convert into ValueMap and store
      std::auto_ptr<edm::ValueMap<float> > chi2ValMap(new edm::ValueMap<float>());
      edm::ValueMap<float>::Filler chi2Filler(*chi2ValMap);
      chi2Filler.insert(probes, vals.begin(), vals.end());
      chi2Filler.fill();
      iEvent.put(chi2ValMap, name);
  };

  fill("sip3D",sip3DVals);
  fill("miniiso",miniIsoVals);

  
}


#include "DataFormats/Math/interface/deltaR.h"

// Version 2.0

// Source:
// V1.0
// https://hypernews.cern.ch/HyperNews/CMS/get/susy/1991.html
// https://github.com/manuelfs/CfANtupler/blob/master/minicfa/interface/miniAdHocNTupler.h#L54
// V2.0:
// Added EA correction option for pile-up
// https://hypernews.cern.ch/HyperNews/CMS/get/b2g-selections/259.html

double Isolations::getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
              const reco::Candidate* ptcl,
              double r_iso_min, double r_iso_max, double kt_scale,
              bool charged_only, bool use_EA_corr=false, double EA_03=0, double rho=0) {

  if (ptcl->pt()<5.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if(ptcl->isElectron()) {
    if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  } else if(ptcl->isMuon()) {
    deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
  } else {
    //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }

  double iso_nh(0.); double iso_ch(0.);
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh(0.5);
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;

    double dr = reco::deltaR(pfc, *ptcl);
    if (dr > r_iso) continue;

    //////////////////  NEUTRALS  /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {
    /////////// PHOTONS ////////////
    if (abs(pfc.pdgId())==22) {
      if(dr < deadcone_ph) continue;
      iso_ph += pfc.pt();
      /////////// NEUTRAL HADRONS ////////////
        } else if (abs(pfc.pdgId())==130) {
      if(dr < deadcone_nh) continue;
      iso_nh += pfc.pt();
    }
      }
      //////////////////  CHARGED from PV  /////////////////////////
    } else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211) {
    if(dr < deadcone_ch) continue;
    iso_ch += pfc.pt();
      }
      //////////////////  CHARGED from PU  /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
    if(dr < deadcone_pu) continue;
    iso_pu += pfc.pt();
      }
    }
  }
  double iso(0.);
  if (charged_only){
    iso = iso_ch;
  } else {
    iso = iso_ph + iso_nh;
    if (use_EA_corr) {
      double EA_miniIso = EA_03 * (r_iso/0.3)*(r_iso/0.3);
      iso -= rho * EA_miniIso;
    } else iso -= 0.5*iso_pu;
    if (iso>0) iso += iso_ch;
    else iso = iso_ch;
  }
  iso = iso/ptcl->pt();

  return iso;
}

double Isolations::electronEA(double eta) {
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
//https://github.com/cmsb2g/B2GAnaFW/blob/v8.0.x_v3.2/src/ElectronUserData.cc
      //These are Effective areas suitable for 80X samples post ICHEP
      float effArea = 0.;
      if(std::fabs(eta)>0.0 &&   std::fabs(eta)<=1.0)   effArea = 0.1703;
      if(std::fabs(eta)>1.0 &&   std::fabs(eta)<=1.479) effArea = 0.1715;
      if(std::fabs(eta)>1.479 && std::fabs(eta)<=2.0) effArea = 0.1213;
      if(std::fabs(eta)>2.0 &&   std::fabs(eta)<=2.2)   effArea = 0.1230;
      if(std::fabs(eta)>2.2 &&   std::fabs(eta)<=2.3)   effArea = 0.1635;
      if(std::fabs(eta)>2.3 &&   std::fabs(eta)<=2.4)   effArea = 0.1937;
      if(std::fabs(eta)>2.4 &&   std::fabs(eta)<=5.0)   effArea = 0.2393;
      return effArea;
}

double Isolations::muonEA(double eta){
//https://github.com/cmsb2g/B2GAnaFW/blob/v8.0.x_v3.2/src/MuonUserData.cc
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=0.8) effArea = 0.0735;
  if(abs(eta)>0.8 && abs(eta)<=1.3) effArea = 0.0619;
  if(abs(eta)>1.3 && abs(eta)<=2.0) effArea = 0.0465;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.0433;
  if(abs(eta)>2.2 && abs(eta)<=2.5) effArea = 0.0577;
  return effArea;
}






#endif
