#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/DigitalStub.h"

#include "DataFormats/Math/interface/deltaPhi.h"

namespace TMTT {

//=== Note configuration parameters.

DigitalStub::DigitalStub(const Settings* settings) :

  // To checkn that DigitalStub is correctly initialized.
  ranInit_(false),
  ranMakeGPinput_(false),
  ranMakeHTinput_(false),
  ranMakeSForTFinput_(""),
  ranMakeHybridInput_(false),

  // Digitization configuration parameters
  phiSectorBits_ (settings->phiSectorBits()), // No. of bits to store phi sector number
  //--- Parameters available in HT board.
  phiSBits_      (settings->phiSBits()),      // No. of bits to store phiS coord.
  phiSRange_     (settings->phiSRange()),     // Range of phiS coord. in radians.
  rtBits_        (settings->rtBits()),        // No. of bits to store rT coord.
  rtRange_       (settings->rtRange()),       // Range of rT coord. in cm.
  zBits_         (settings->zBits()),         // No. of bits to store z coord.
  zRange_        (settings->zRange()),        // Range of z coord in cm.
  //--- Parameters available in GP board (excluding any in common with HT specified above).
  phiOBits_      (settings->phiOBits()),      // No. of bits to store phiO parameter.
  phiORange_     (settings->phiORange()),     // Range of phiO parameter
  bendBits_      (settings->bendBits()),      // No. of bits to store stub bend.

  // Note if using reduced layer ID, so tracker layer can be encoded in 3 bits.
  reduceLayerID_ (settings->reduceLayerID()),

  // Number of phi sectors and phi nonants.
  numPhiSectors_ (settings->numPhiSectors()),
  numPhiNonants_ (settings->numPhiNonants()),
  // Phi sector and phi nonant width (radians)
  phiSectorWidth_(2.*M_PI / float(numPhiSectors_)), 
  phiNonantWidth_(2.*M_PI / float(numPhiNonants_)), 
  // Radius from beamline with respect to which stub r coord. is measured.
  chosenRofPhi_  (settings->chosenRofPhi()),

  // Number of q/Pt bins in Hough  transform array.
  nbinsPt_       ((int) settings->houghNbinsPt())
{
  // Calculate multipliers to digitize the floating point numbers.
  phiSMult_ = pow(2, phiSBits_)/phiSRange_;
  rtMult_   = pow(2, rtBits_  )/rtRange_;
  zMult_    = pow(2, zBits_   )/zRange_;
  phiOMult_ = pow(2, phiOBits_)/phiORange_;

  bendMult_ = 4.; // No precision lost by digitization, since original bend (after encoding) has steps of 0.25 (in units of pitch).
  bendRange_ = round(pow(2, bendBits_)/bendMult_); // discrete values, so digitisation different 
}

//=== Initialize stub with original, floating point stub coords,
//=== range of m bin (= q/Pt bin) values allowed by bend filter, 
//=== normal & "reduced" tracker layer of stub, stub bend, and pitch & seperation of module,
//=== and half-length of strip or pixel in r and in z, and if it's in barrel, tilted barrel and/or PS module.

void DigitalStub::init(float phi_orig, float r_orig, float z_orig,
           float alpha,
           unsigned int min_qOverPt_bin_orig, unsigned int max_qOverPt_bin_orig, 
           unsigned int layerID, unsigned int layerIDreduced, float bend_orig,
	   float pitch, float sep, float rErr, float zErr, bool barrel, bool tiltedBarrel, bool psModule) {

  ranInit_ = true; // Note we ran init().
  // Variables in HT.
  phi_orig_             = phi_orig; 
  r_orig_               = r_orig;
  z_orig_               = z_orig;
  min_qOverPt_bin_orig_ = min_qOverPt_bin_orig;
  max_qOverPt_bin_orig_ = max_qOverPt_bin_orig;
  layerID_              = layerID;
  layerIDreduced_       = layerIDreduced;
  // Variables exclusively in GP.
  bend_orig_            = bend_orig;
  rErr_orig_            = rErr;
  zErr_orig_            = zErr;
  alpha_ = alpha;

  isPSModule_ = psModule;

  // Calculate unique module type ID, allowing pitch/sep of module to be determined.
  // (N.B. Module types 0 & 1 have identical pitch & sep, but one is tilted & one is flat barrel module).

  moduleType_ = 999;
  // EJC new (additional) module type for tilted geometry
  const vector<float> pitchVsType  = {0.0099, 0.0099, 0.0099, 0.0099, 0.0089, 0.0099, 0.0089, 0.0089};
  const vector<float> sepVsType    = {0.26  , 0.26  , 0.16  , 0.4   , 0.18  , 0.4   , 0.18  , 0.4   };
  const vector<bool>  barrelVsType = {true  , true  , true  , true  , true  , false , false , false };
  const vector<bool>  psVsType     = {true  , true  , true  , true  , false , true  , false , false };
  const vector<bool>  tiltedVsType = {false , true  , false , true  , false,  false , false , false };
  if (pitchVsType.size() != sepVsType.size()) throw cms::Exception("DigitalStub: module type array size wrong");
  const float tol = 0.001; // Tolerance
  for (unsigned int i = 0; i < pitchVsType.size(); i++) {
    if (fabs(pitch - pitchVsType[i]) < tol && fabs(sep - sepVsType[i]) < tol && barrel == barrelVsType[i] && tiltedBarrel == tiltedVsType[i] && psModule == psVsType[i]) {
      moduleType_ = i;
    }
  }
  if (moduleType_ == 999) throw cms::Exception("DigitalStub: unknown module type")<<"pitch="<<pitch<<" separation="<<sep<<" barrel="<<barrel<<" tilted="<<tiltedBarrel<<" PS="<<psModule<<endl;
}

//=== Digitize stub for input to Geographic Processor, with stub phi coord. measured relative to phi nonant that contains specified phi sector.

void DigitalStub::makeGPinput(unsigned int iPhiSec) {

  if (! ranInit_) throw cms::Exception("DigitalStub: You forgot to call init() before makeGPinput()!");

  unsigned int iPhiNon = floor(iPhiSec*numPhiNonants_/numPhiSectors_); // Find nonant corresponding to this sector.

  // If this stub was already digitized, we don't have to redo all the work again. Save CPU.
  if (ranMakeGPinput_) {
    if (iPhiNon == iDigi_Nonant_) {
      return; // Work already done.
    } else {
      this->quickMakeGPinput(iPhiSec);
      return;
    }
  }

  ranMakeGPinput_ = true; // Note we ran make().

  //--- Shift axes of coords. if required.

  // r coordinate relative to specified point.
  rt_orig_ = r_orig_ - chosenRofPhi_;

  // Phi coord. of stub relative to centre of nonant.
  double phiNonantCentre = phiNonantWidth_ * double(iPhiNon);

  phiO_orig_ = reco::deltaPhi(phi_orig_, phiNonantCentre);

  //--- Digitize variables used exclusively in GP.
  iDigi_Nonant_ = iPhiNon;
  iDigi_PhiO_   = floor(phiO_orig_*phiOMult_);
  iDigi_Bend_   = round(bend_orig_*bendMult_);   // discrete values, so digitisation different 
  //--- Digitize variables used in both GP & HT.
  iDigi_Rt_     = floor(rt_orig_*rtMult_);
  iDigi_Z_      = floor(z_orig_*zMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used exclusively in GP.
  phiO_      = (iDigi_PhiO_ + 0.5)/phiOMult_;
  bend_      = iDigi_Bend_/bendMult_;  // discrete values, so digitisation different
  phi_       = reco::deltaPhi(phiO_, -phiNonantCentre);
  //--- Then for variables used in both GP & HT.
  rt_        = (iDigi_Rt_ + 0.5)/rtMult_;  
  r_         = rt_ + chosenRofPhi_;
  z_         = (iDigi_Z_ + 0.5)/zMult_; 
}

//=== Digitize stub for input to Hough transform, with stub phi coord. measured relative to specified phi sector.

void DigitalStub::makeHTinput(unsigned int iPhiSec) {

  if (! ranInit_) throw cms::Exception("DigitalStub: You forgot to call init() before makeHTinput()!");

  // Digitize for GP input if not already done, since some variables are shared by GP & HT.
  this->makeGPinput(iPhiSec);

  // If this stub was already digitized, we don't have to redo all the work again. Save CPU.
  if (ranMakeHTinput_) {
    if (iPhiSec == iDigi_PhiSec_) {
      return; // Work already done.
    } else {
      this->quickMakeHTinput(iPhiSec);
      return;
    }
  }

  ranMakeHTinput_ = true; // Note we ran makeHTinput().

  //--- Shift axes of coords. if required.

  // Centre of this sector in phi
  double phiCentreSec0 = -M_PI/float(numPhiNonants_) + M_PI/float(numPhiSectors_);
  double phiSectorCentre = phiSectorWidth_ * float(iPhiSec) + phiCentreSec0; 

  // Point in sector from which stub phiS should be measured.
  double phiSectorRef = phiSectorCentre;

  // Phi coord. of stub relative to centre of sector.
  phiS_orig_ = reco::deltaPhi(phi_orig_, phiSectorRef); 

  //--- Digitize variables used exclusively by HT.
  iDigi_PhiSec_ = iPhiSec;
  iDigi_PhiS_   = floor(phiS_orig_*phiSMult_);
  // Don't bother digitising here variables used by both GP & HT, as makeGPinput() will already have digitized them.

  // N.B. If using daisy-chain firmware, then should logically recalculate m bin range here, since it depends on the now
  // digitized r and z coordinates. But too lazy to move code here from Stub::digitizeForHTinput(), where calculation
  // actually done.
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used exclusively by HT
  phiS_      = (iDigi_PhiS_ + 0.5)/phiSMult_;
  phi_       = reco::deltaPhi(phiS_, -phiSectorRef); // N.B. phi_ measured w.r.t sector here, but w.r.t. nonant in makeGPinput()
  // Don't bother with  variables used by both GP & HT, as makeGPinput() will already have already calculated them.

  //--- Do next two checks here rather than in makeGPinput(), as the latter may be called for stubs in the wrong nonant, so out of range.

  // Check that stub coords. are within assumed digitization range.
  this->checkInRange();

  // Check that digitization followed by undigitization doesn't change results too much.
  this->checkAccuracy();

  // Adjust m-bin range, as hardware counts q/Pt bins in HT array using a signed integer in a symmetric range about zero.
  const int min_array_bin = (nbinsPt_%2 == 0)  ?  -(nbinsPt_/2)      :  -(nbinsPt_ - 1)/2;
  m_min_ = min_qOverPt_bin_orig_ + min_array_bin;
  m_max_ = max_qOverPt_bin_orig_ + min_array_bin;

  //--- Produce tracker layer identifier, encoded as it is sent along the optical link.
  
  if (reduceLayerID_) {
    // Firmware is using "reduced" layer ID, which can be packed into 3 bits in range 1-7.
    iDigi_LayerID_ = layerIDreduced_;

  } else {
    // Firmware is using normal layer ID, which needs more than 3 bits to store it.
    // Encode barrel layers as 0 to 5.
    iDigi_LayerID_ = layerID_ - 1;
    // Endcode endcap layers as 6 to 10, not bothering to distinguish the two endcaps.
    if        (iDigi_LayerID_ == 10 || iDigi_LayerID_ == 20) {
      iDigi_LayerID_ = 6;
    } else if (iDigi_LayerID_ == 11 || iDigi_LayerID_ == 21) {
      iDigi_LayerID_ = 7;
    } else if (iDigi_LayerID_ == 12 || iDigi_LayerID_ == 22) {
      iDigi_LayerID_ = 8;
    } else if (iDigi_LayerID_ == 13 || iDigi_LayerID_ == 23) {
      iDigi_LayerID_ = 9;
    } else if (iDigi_LayerID_ == 14 || iDigi_LayerID_ == 24) {
      iDigi_LayerID_ = 10;
    }
  }
}

//=== Digitize stub for input to r-z Seed Filter or Track Fitter.
//=== Argument is "SeedFilter" or name of Track Fitter.
//=== N.B. This digitisation is done internally within the FPGA not for transmission along opto-links.

void DigitalStub::makeSForTFinput(string SForTF) {
  if (! ranInit_) throw cms::Exception("DigitalStub: You forgot to call init() before makeSForTFinput()!");

  // Save CPU by not digitizing stub again if already done.
  if (ranMakeSForTFinput_ != SForTF) {
    ranMakeSForTFinput_ = SForTF; // Note we ran makeSForTFinput().

    // The stub r coordinate is calculated inside the seed filter from the value of rT.
    // Best to choose a ref. radius such that digitising then undigitising it leaves it unchanged.
    iDigi_R_      = iDigi_Rt_ + std::round(chosenRofPhi_*rtMult_); 

    //--- Determine floating point variables from digtized numbers.
    r_            = (iDigi_R_ + 0.5)/rtMult_;  
    rt_           = r_ - chosenRofPhi_;

    if (SForTF.find("KF") != string::npos) { 
      // Digitize variables that are exclusive to Kalman filter.
      // Data format uses z directly from HT, as its multiplier is exactly a factor 2
      // smaller than for r, which is easy to fix in VHDL-HLS interface.
      iDigi_Z_KF_ = iDigi_Z_;
      // Determine floating point variable from digitized one, using z multiplier.
      z_          = (iDigi_Z_KF_ + 0.5)/zMult_;
    } else {
      // If not using KF, then restore z value from HT.
      iDigi_Z_KF_ = 0;  // Shoudln't be used in this case.
      z_          = (iDigi_Z_ + 0.5)/zMult_; 
    }

    if (SForTF == "SeedFilter") {
      //--- Digitize variables used exclusively in seed filter.
      // The stub (r,z) uncertainties are actually calculated inside the seed filter and not passed to it along optical links.
      iDigi_rErr_   = ceil(rErr_orig_*rtMult_); // Round up to avoid zero uncertainty ...
      iDigi_zErr_   = ceil(zErr_orig_*zMult_);
      //--- Determine floating point variables from digtized numbers
      rErr_         = (iDigi_rErr_ - 0.5)/rtMult_;
      zErr_         = (iDigi_zErr_ - 0.5)/zMult_;
    } else {
      // Restore original value. Although perhaps better to leave?
      rErr_         = rErr_orig_;
      zErr_         = zErr_orig_;
    }
  }
}

void DigitalStub::makeDRinput(unsigned int stubId){
  if (! ranInit_) throw cms::Exception("DigitalStub: You forgot to call init() before makeDRinput()!");

  ranMakeDRinput_ = true; // Note we ran makeDRinput().
  stubId_ = stubId;  
}

void DigitalStub::makeHybridInput(unsigned int iPhiSec, const Settings* settings) {

  this->makeHTinput( iPhiSec );

  unsigned int iPhiNon = floor(iPhiSec*numPhiNonants_/numPhiSectors_); // Find nonant corresponding to this sector.

  ranMakeHybridInput_ = true; // Note we ran makeHybridInput().

  iDigi_Hybrid_Layer_ = layerID_ - 1;

  // Hybrid endcap layer ids are lower for +ve z, larger for -ve z
  if ( iDigi_Hybrid_Layer_ >= 20 ) iDigi_Hybrid_Layer_ -= 14;
  else if ( iDigi_Hybrid_Layer_ >= 10 ) iDigi_Hybrid_Layer_ += 1;

  //--- Shift axes of coords. if required.

  // Get reference value of z/r to use
  unsigned int layerIndex = iDigi_Hybrid_Layer_;
  if ( layerIndex > 10 ) layerIndex -= 5;
  double r_ref = settings->getLayerRReference()[layerIndex];
  double z_ref = settings->getLayerZReference()[layerIndex];

  float rMultiplier = 0;
  float zMultiplier = 0;
  float phiMultiplier = 0;
  double zRange = 0;
  if ( iDigi_Hybrid_Layer_ <= 5 ) {
    if ( isPSModule() ) {
      rMultiplier = pow(2, settings->hybrid_rBits_BarrelPS()  ) / settings->hybrid_rRange_BarrelPS();
      zMultiplier = pow(2, settings->hybrid_zBits_BarrelPS()  ) / settings->hybrid_zRange_BarrelPS();
      phiMultiplier = pow(2, settings->hybrid_phiBits_BarrelPS()  ) / settings->hybrid_phiRange_BarrelPS();
    }
    else {
      rMultiplier = pow(2, settings->hybrid_rBits_Barrel2S()  ) / settings->hybrid_rRange_Barrel2S();
      zMultiplier = pow(2, settings->hybrid_zBits_Barrel2S()  ) / settings->hybrid_zRange_Barrel2S();
      phiMultiplier = pow(2, settings->hybrid_phiBits_Barrel2S()  ) / settings->hybrid_phiRange_Barrel2S();
    }
  }
  else {
    if ( isPSModule() ) {
      rMultiplier = pow(2, settings->hybrid_rBits_EndcapPS()  ) / settings->hybrid_rRange_EndcapPS();
      zMultiplier = pow(2, settings->hybrid_zBits_EndcapPS()  ) / settings->hybrid_zRange_EndcapPS();
      phiMultiplier = pow(2, settings->hybrid_phiBits_EndcapPS()  ) / settings->hybrid_phiRange_EndcapPS();
    }
    else {
      zMultiplier = pow(2, settings->hybrid_zBits_Endcap2S()  ) / settings->hybrid_zRange_Endcap2S();
      phiMultiplier = pow(2, settings->hybrid_phiBits_Endcap2S()  ) / settings->hybrid_phiRange_Endcap2S();
    }
  }

  int sign = (z_orig_ >= 0) ? 1 : -1;
  iDigi_Hybrid_Z_ = floor( ( z_orig_ - sign*z_ref ) * zMultiplier );
  iDigi_Hybrid_Phi_ = floor( phiO_orig_ * phiMultiplier );

  if ( iDigi_Hybrid_Layer_ <= 5 || this->isPSModule() ) {
    // Barrel and PS endcap modules
    iDigi_Hybrid_R_ = floor( ( r_orig_ - r_ref ) * rMultiplier );
    iDigi_Hybrid_Alpha_ = 0;
  }
  else {
    // 2S endcap modules
    vector<double> moduleR_2S = settings->getDisk1To2_2S_R();

    if ( ( iDigi_Hybrid_Layer_ >= 8 && iDigi_Hybrid_Layer_ <= 10 ) || iDigi_Hybrid_Layer_ >= 13 ) {
      moduleR_2S = settings->getDisk3To5_2S_R();
    }

    int r_index = -1;
    for ( unsigned int i = 0; i < moduleR_2S.size(); ++i ) {
      if ( fabs( r_orig_ - moduleR_2S[i] ) < 0.2 ) {
        r_index = i;
        break;
      }
    }

    if ( r_index < 0 ) {
       throw cms::Exception("Couldn't find radius of 2S endcap module corresponding to this stub")<<" Stub r = " <<r_orig_<<endl;  
    }
    iDigi_Hybrid_R_ = r_index;


    // alpha
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HybridDataFormat
    // "alpha is simply defined as the strip number in the module (with coarse granularity)."
    float alphaMultiplier = pow(2, settings->hybrid_alphaBits_Endcap2S()  ) / settings->hybrid_alphaRange_Endcap2S();
    // iDigi_Hybrid_Alpha_ = floor( iPhi_ * alphaMultiplier );
    // But it's actually scaled between -1 and 1 : https://gitlab.cern.ch/cms-tracker-phase2-backend-development/BE_software/L1Tracking/blob/master/TrackFindingTracklet/interface/L1TStub.h#L223-232
    // and then digitised.
    iDigi_Hybrid_Alpha_ = floor( alpha_ * alphaMultiplier );

  }


  // Bend
  int hybridBend=2.0*bend_orig_;
  if ( this->isPSModule() ) {
    if ( hybridBend >= 0 ) {
        iDigi_Hybrid_Bend_ = std::min( floor( hybridBend / 2 ), 3. );
    }
    else {
        iDigi_Hybrid_Bend_ = std::min( fabs( floor( ( hybridBend - 1 ) / 2 ) ) + 3, 7. );
    }
    // Above is more compact
    // Below is implementation from hybrid code
    // if (hybridBend==0||hybridBend==1) iDigi_Hybrid_Bend_ = 0;
    // else if (hybridBend==2||hybridBend==3) iDigi_Hybrid_Bend_ = 1;
    // else if (hybridBend==4||hybridBend==5) iDigi_Hybrid_Bend_ = 2;
    // else if (hybridBend>=6) iDigi_Hybrid_Bend_ = 3;
    // else if (hybridBend==-1||hybridBend==-2) iDigi_Hybrid_Bend_ = 4;
    // else if (hybridBend==-3||hybridBend==-4) iDigi_Hybrid_Bend_ = 5;
    // else if (hybridBend==-5||hybridBend==-6) iDigi_Hybrid_Bend_ = 6;
    // else if (hybridBend<=-7) iDigi_Hybrid_Bend_ = 7;
  }
  else {

    if ( hybridBend >= 0 ) {
        iDigi_Hybrid_Bend_ = std::min( floor( hybridBend / 2 ), 7. );
    }
    else {
        iDigi_Hybrid_Bend_ = std::min( fabs( floor( ( hybridBend - 1 ) / 2 ) ) + 7, 15. );
    }
    // Above is more compact
    // Below is implementation from hybrid code
    // if (hybridBend==0||hybridBend==1) iDigi_Hybrid_Bend_ = 0;
    // else if (hybridBend==2||hybridBend==3) iDigi_Hybrid_Bend_ = 1;
    // else if (hybridBend==4||hybridBend==5) iDigi_Hybrid_Bend_ = 2;
    // else if (hybridBend==6||hybridBend==7) iDigi_Hybrid_Bend_ = 3;
    // else if (hybridBend==8||hybridBend==9) iDigi_Hybrid_Bend_ = 4;
    // else if (hybridBend==10||hybridBend==11) iDigi_Hybrid_Bend_ = 5;
    // else if (hybridBend==12||hybridBend==13) iDigi_Hybrid_Bend_ = 6;
    // else if (hybridBend>=14) hybridBend = 7;
    // else if (hybridBend==-1||hybridBend==-2) iDigi_Hybrid_Bend_ = 8;
    // else if (hybridBend==-3||hybridBend==-4) iDigi_Hybrid_Bend_ = 9;
    // else if (hybridBend==-5||hybridBend==-6) iDigi_Hybrid_Bend_ = 10;
    // else if (hybridBend==-7||hybridBend==-8) iDigi_Hybrid_Bend_ = 11;
    // else if (hybridBend==-9||hybridBend==-10) iDigi_Hybrid_Bend_ = 12;
    // else if (hybridBend==-10||hybridBend==-11) iDigi_Hybrid_Bend_ = 13;
    // else if (hybridBend==-12||hybridBend==-13) iDigi_Hybrid_Bend_ = 14;
    // else if (hybridBend<=-14) hybridBend = 15;
  }
}


//=== Redigitize stub for input to Geographic Processor, if it was previously digitized for a different phi sector.

void DigitalStub::quickMakeGPinput(int iPhiSec) {

  //--- Shift axes of coords. if required.

  // Phi coord. of stub relative to centre of nonant.
  unsigned int iPhiNon = floor(iPhiSec*numPhiNonants_/numPhiSectors_);
  double phiNonantCentre = phiNonantWidth_ * double(iPhiNon);
  phiO_orig_ = reco::deltaPhi(phi_orig_, phiNonantCentre);

  //--- Digitize variables used exclusively in GP.
  iDigi_Nonant_ = iPhiNon;
  iDigi_PhiO_   = floor(phiO_orig_*phiOMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- Variables used exclusively in GP.
  phiO_      = (iDigi_PhiO_)/phiOMult_;
  phi_       = reco::deltaPhi(phiO_, -phiNonantCentre);
}

//=== Redigitize stub for input to Hough transform, if it was previously digitized for a different phi sector.

void DigitalStub::quickMakeHTinput(int iPhiSec) {

  //--- Shift axes of coords. if required.

  // Centre of this sector in phi
  double phiCentreSec0 = -M_PI/float(numPhiNonants_) + M_PI/float(numPhiSectors_);
  double phiSectorCentre = phiSectorWidth_ * float(iPhiSec) + phiCentreSec0; 

  // Point in sector from which stub phiS should be measured.
  double phiSectorRef = phiSectorCentre;

  // Phi coord. of stub relative to centre of sector.
  phiS_orig_ = reco::deltaPhi(phi_orig_, phiSectorRef); 

  // Check that stub coords. are within assumed digitization range.
  this->checkInRange();

  //--- Digitize variables used in HT.
  iDigi_PhiSec_ = iPhiSec;
  iDigi_PhiS_   = floor(phiS_orig_*phiSMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used in HT.
  phiS_      = (iDigi_PhiS_)/phiSMult_;
  phi_       = reco::deltaPhi(phiS_, -phiSectorRef);
}

//=== Check that stub coords. are within assumed digitization range.

void DigitalStub::checkInRange() const {
  // All ranges are centred at zero, except for rho, which is +ve-definate.
  if (fabs(phiS_orig_) >= 0.5*phiSRange_)   throw cms::Exception("DigitalStub: Stub phiS is out of assumed digitization range.")<<" |phiS| = " <<fabs(phiS_orig_) <<" > "<<0.5*phiSRange_<<endl;  
  if (fabs(rt_orig_)   >= 0.5*rtRange_)     throw cms::Exception("DigitalStub: Stub rT is out of assumed digitization range.")  <<" |rt| = "   <<fabs(rt_orig_)   <<" > "<<0.5*rtRange_  <<endl; 
  if (fabs(z_orig_)    >= 0.5*zRange_)      throw cms::Exception("DigitalStub: Stub z is out of assumed digitization range.")   <<" |z| = "    <<fabs(z_orig_)   <<" > "<<0.5*zRange_  <<endl;  
  if (fabs(phiO_orig_)   >= 0.5*phiORange_) throw cms::Exception("DigitalStub: Stub phiO is out of assumed digitization range.")<<" |phiO| = " <<fabs(phiO_orig_)<<" > "<<0.5*phiORange_  <<endl;  
  if (fabs(bend_orig_)   >= 0.5*bendRange_) throw cms::Exception("DigitalStub: Stub bend is out of assumed digitization range.")<<" |bend| = " <<fabs(bend_orig_)<<" > "<<0.5*bendRange_  <<endl;  
}

//=== Check that digitisation followed by undigitisation doesn't change significantly the stub coordinates.

void DigitalStub::checkAccuracy() const {
  float TA = reco::deltaPhi(phi_, phi_orig_);
  float TB = r_    - r_orig_;
  float TC = z_    - z_orig_;
  float TD = phiO_ - phiO_orig_;
  float TE = bend_ - bend_orig_;

  static unsigned int nErr = 0;
  const  unsigned int maxErr = 20; // Print error message only this number of times.
  if (nErr < maxErr) {
    if (fabs(TA) > 0.001 || fabs(TB) > 0.3 || fabs(TC) > 0.25 || fabs(TD) > 0.005 || fabs(TE) > 0.01) {
      nErr++;
      cout<<"WARNING: DigitalStub lost precision: "<<TA<<" "<<TB<<" "<<TC<<" "<<TD<<" "<<TE<<endl;
    }
  }
}

}
