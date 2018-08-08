#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrack.h"

namespace TMTT {

// Digitize track and degrade helix parameter resolution according to effect of digitisation.

void L1fittedTrack::digitizeTrack(const string& fitterName){
  if (settings_->enableDigitize()) {
    if (! digitizedTrack_ ) {
      digitizedTrack_ = true;
     
      bool consistent = this->consistentHTcell();
      int  mbinhelix  = int(this->getCellLocationFit().first) - floor(settings_->houghNbinsPt()/2);
      int  cbinhelix  = int(this->getCellLocationFit().second) - floor(settings_->houghNbinsPhi()/2);
      int  mBinHT     = int(this->getCellLocationHT().first) - floor(settings_->houghNbinsPt()/2);
      int  cBinHT     = int(this->getCellLocationHT().second) - floor(settings_->houghNbinsPhi()/2);

      if(matchedTP_ != nullptr){
        digitalTrack_.init(fitterName, nHelixParam_,
         iPhiSec_, iEtaReg_, mBinHT, cBinHT, mbinhelix, cbinhelix, 
	 qOverPt_, d0_, phi0_,tanLambda_, z0_, chi2_,
	 qOverPt_bcon_, phi0_bcon_, chi2_bcon_,
	 nLayers_, consistent, this->accepted(),
	 matchedTP_->qOverPt(), matchedTP_->d0(), matchedTP_->phi0(), matchedTP_->tanLambda(), matchedTP_->z0(), matchedTP_->eta(), 
	 matchedTP_->index(), matchedTP_->useForAlgEff(), matchedTP_->useForEff(), matchedTP_->pdgId());
      } else {
        digitalTrack_.init(fitterName, nHelixParam_,
         iPhiSec_, iEtaReg_, mBinHT, cBinHT, mbinhelix, cbinhelix,       
	 qOverPt_, d0_, phi0_, tanLambda_, z0_, chi2_, 
	 qOverPt_bcon_, phi0_bcon_, chi2_bcon_,
	 nLayers_, consistent, this->accepted(), 
	 0, 0, 0, 0, 0, 0, 
	 -1, 0, 0, 0);
      }

      // Digitize track
      digitalTrack_.makeDigitalTrack();

      // Convert digitized track params back to floating point with degraded resolution.
      qOverPt_   = digitalTrack_.qOverPt();
      d0_        = digitalTrack_.d0();
      phi0_      = digitalTrack_.phi0();
      z0_        = digitalTrack_.z0();
      tanLambda_ = digitalTrack_.tanLambda();
      chi2_      = digitalTrack_.chisquared();

      // Ditto for beam-spot constrained values.
      qOverPt_bcon_   = digitalTrack_.qOverPt_bcon();
      phi0_bcon_      = digitalTrack_.phi0_bcon();
      chi2_bcon_      = digitalTrack_.chisquared_bcon();
    }
  }
}

}
