#ifndef _SCACTIVITYVARS_H
#define _SCACTIVITYVARS_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTTrackIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include <TVector2.h>


template <class T>
class SCActivityVars : public edm::EDProducer {
public:
    explicit SCActivityVars(const edm::ParameterSet & iConfig);
    virtual ~SCActivityVars() ;

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >  CylLorentzVectorF;
    float electronEA(float eta);
private:
    const edm::EDGetTokenT<std::vector<T> > probesToken_;
    const edm::EDGetTokenT<pat::PackedCandidateCollection > token_pfCand;
    const edm::EDGetTokenT<reco::VertexCollection>  token_vtx;
    const edm::EDGetTokenT<double>                             token_rho      ;

    const bool applyVertex_ ;
    const double trkIsoConeSize_;
    const double trkIsoDeltaEtaVeto_;
    const double trkIsoDeltaPhiVeto_;
    const double trkIsoConeVeto_;
};

template<class T>
SCActivityVars<T>::SCActivityVars(const edm::ParameterSet & iConfig) :
probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
token_pfCand      (consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfParticles"))),
token_vtx         (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
token_rho       (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
applyVertex_       (iConfig.getParameter<bool>("applyVertex")),
trkIsoConeSize_     (iConfig.getParameter<double>("trkIsoConeSize")),
trkIsoDeltaEtaVeto_    (iConfig.getParameter<double>("trkIsoDeltaEtaVeto")),
trkIsoDeltaPhiVeto_    (iConfig.getParameter<double>("trkIsoDeltaPhiVeto")),
trkIsoConeVeto_    (iConfig.getParameter<double>("trkIsoConeVeto"))
{

    produces<edm::ValueMap<float> >("scIsoTrkPT");
    produces<edm::ValueMap<float> >("scIsoTrkDR");
    produces<edm::ValueMap<float> >("scIsoPT");
    produces<edm::ValueMap<float> >("scIsoDR");
}

template<class T>
SCActivityVars<T>::~SCActivityVars()
{}

template<class T>
void SCActivityVars<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

    // read input
    edm::Handle<std::vector<T> > probes;
    iEvent.getByToken(probesToken_,  probes);

    edm::Handle<reco::VertexCollection>             han_vtx        ;
    edm::Handle<pat::PackedCandidateCollection>     han_pfCand     ;
    edm::Handle<double>                             han_rho      ;

    iEvent.getByToken(token_vtx      ,han_vtx      );
    iEvent.getByToken(token_pfCand   ,han_pfCand   );
    iEvent.getByToken(token_rho    ,han_rho    );

    // prepare vector for output
    std::vector<float> scIsoPT;scIsoPT.reserve(probes->size());
    std::vector<float> scIsoDR;scIsoDR.reserve(probes->size());
    std::vector<float> scIsoTrkPT;scIsoTrkPT.reserve(probes->size());
    std::vector<float> scIsoTrkDR;scIsoTrkDR.reserve(probes->size());

    auto vtx_pt =  (applyVertex_ && han_vtx->size() > 0) ?  (*han_vtx)[0].position() : reco::Vertex::Point();
    const float drc2 = trkIsoConeSize_*trkIsoConeSize_;
    const float drv2 = trkIsoConeVeto_*trkIsoConeVeto_;
    const float eaCorrRho = (*han_rho)*drc2/(0.3*0.3);

    for(const auto& probe : *probes){
        const auto& pos = ((const reco::RecoCandidate *)&probe)->superCluster()->position();
        GlobalVector mom(pos.x()-vtx_pt.x(),pos.y()-vtx_pt.y(),pos.z()-vtx_pt.z());
        CylLorentzVectorF pIso;
        CylLorentzVectorF pTkIso;
        float iso_charged=0;
        float iso_neutral=0;

        for(const auto& can : *han_pfCand){
            float deta = can.eta()-mom.eta();
            float dphi = TVector2::Phi_mpi_pi(can.phi()-mom.phi()) ;
            float dr2 = (deta*deta + dphi*dphi);
            if(dr2 > drc2) continue;
            if (can.charge()!=0 && can.fromPV()>1){
                iso_charged += can.pt();
                pIso += can.p4();
            } else if(can.charge() == 0){
                iso_neutral += can.pt();
                pIso += can.p4();
            }
            if(can.pdgId() != 211) continue;
            if(std::fabs(dphi) < trkIsoDeltaPhiVeto_) continue;
            if(std::fabs(deta) < trkIsoDeltaEtaVeto_) continue;
            if(trkIsoConeVeto_ > 0 &&  dr2  < drv2) continue;
            pTkIso += can.p4();
        }
        auto getDR = [&](const CylLorentzVectorF& isoMom ) -> float{
            float deta = isoMom.eta()-mom    .eta();
            float dphi = TVector2::Phi_mpi_pi(isoMom.phi()-mom.phi()) ;
            return std::sqrt((deta*deta + dphi*dphi));
        };
        float eACorrectedIso = iso_charged + std::max(0.0f, iso_neutral - electronEA(pos.eta())*eaCorrRho);
        scIsoPT.push_back(eACorrectedIso);
        scIsoDR.push_back(getDR(pIso));
        scIsoTrkPT.push_back(pTkIso.pt());
        scIsoTrkDR.push_back(getDR(pTkIso));
    }

    auto fill = [&](const std::string& name, const std::vector<float>& vs) {
        // convert into ValueMap and store
        std::auto_ptr<edm::ValueMap<float> > scIsoValMap(new edm::ValueMap<float>());
        edm::ValueMap<float>::Filler scIsoFiller(*scIsoValMap);
        scIsoFiller.insert(probes, vs.begin(), vs.end());
        scIsoFiller.fill();
        iEvent.put(scIsoValMap, name);
    };

    fill("scIsoPT",scIsoPT);
    fill("scIsoDR",scIsoDR);
    fill("scIsoTrkPT",scIsoTrkPT);
    fill("scIsoTrkDR",scIsoTrkDR);



}

template<class T>
float SCActivityVars<T>::electronEA(float eta) {
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

#endif
