#include "L1Trigger/TrackFindingTMTT/interface/L1track2D.h"

namespace TMTT {

// Digitize track and degrade helix parameter resolution according to effect of digitisation.

void L1track2D::digitizeTrack(){
  if (settings_->enableDigitize()) {
    if (! digitizedTrack_ ) {
      digitizedTrack_ = true;
     
      int  mBinHT     = int(this->getCellLocationHT().first) - floor(settings_->houghNbinsPt()/2);
      int  cBinHT     = int(this->getCellLocationHT().second) - floor(settings_->houghNbinsPhi()/2);

  // void init(const string& fitterName, unsigned int nHelixParams,
  //     unsigned int iPhiSec, unsigned int iEtaReg, int mbin, int cbin, int mBinhelix, int cBinhelix, 
  //     unsigned int hitPattern,
  //     float qOverPt_orig, float d0_orig, float phi0_orig, float tanLambda_orig, float z0_orig, float chisquared_orig, 
  //     float qOverPt_bcon_orig, float phi0_bcon_orig, float chisquared_bcon_orig, // beam-spot constrained values. 
  //     unsigned int nLayers, bool consistent, bool accepted, 
  //     float tp_qOverPt, float tp_d0, float tp_phi0, float tp_tanLambda, float tp_z0, float tp_eta, 
  //     int tp_index, bool tp_useForAlgEff, bool tp_useForEff, int tp_pdgId);


      if(matchedTP_ != nullptr){
        digitalTrack_.init("MiniHT", 0,
	 iPhiSec_, iEtaReg_, mBinHT, cBinHT, 0, 0, 0,
	 this->qOverPt(), 0, this->phi0(),estTanLambda_, estZ0_, 0,
	 0, 0, 0,
	 nLayers_, true, true,
	 matchedTP_->qOverPt(), matchedTP_->d0(), matchedTP_->phi0(), matchedTP_->tanLambda(), matchedTP_->z0(), matchedTP_->eta(), 
	 matchedTP_->index(), matchedTP_->useForAlgEff(), matchedTP_->useForEff(), matchedTP_->pdgId());
      } else {
        digitalTrack_.init("MiniHT", 0,
	 iPhiSec_, iEtaReg_, mBinHT, cBinHT, -1, -1, 0,     
	 this->qOverPt(), 0, this->phi0(), estTanLambda_, estZ0_, 0, 
	 0, 0, 0,
	 nLayers_, true, true, 
	 0, 0, 0, 0, 0, 0, 
	 -1, 0, 0, 0);
      }

      // Digitize track
      digitalTrack_.makeDigitalTrack();

      // // Convert digitized track params back to floating point with degraded resolution.
      // qOverPt_   = digitalTrack_.qOverPt();
      // d0_        = digitalTrack_.d0();
      // phi0_      = digitalTrack_.phi0();
      // z0_        = digitalTrack_.z0();
      // tanLambda_ = digitalTrack_.tanLambda();
      // chi2_      = digitalTrack_.chisquared();

      // // Ditto for beam-spot constrained values.
      // qOverPt_bcon_   = digitalTrack_.qOverPt_bcon();
      // phi0_bcon_      = digitalTrack_.phi0_bcon();
      // chi2_bcon_      = digitalTrack_.chisquared_bcon();
    }
  }
}

}
