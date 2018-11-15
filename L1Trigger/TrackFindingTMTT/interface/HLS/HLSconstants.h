#ifndef __HLSconstants__
#define __HLSconstants__

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#endif

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

//-- Copied from Maxeller code Constants.maxj

// Digitisation multipliers (from data format doc).
// KF uses same multiplier for r as for stubs in DTC, but one extra bit to accomodate larger range,
// since KF measures r w.r.t. beamline. And it uses r multiplier for z too.
static const float rMult = pow(2.,BSR-1)/103.1103;
static const float phiMult = pow(2.,BSP)/0.698131700;
static const float rphiMult = rMult*phiMult;
static const float inv2R_Mult = (phiMult/rMult);

// Beam spot length & reference radii w.r.t. beamline.
static const float beamSpotLength= 15.0;
static const float chosenRofPhi = 61.273;

static const StubHLS::TR chosenRofPhi_digi = chosenRofPhi*rMult;
static const float chosenRofZ = 50.0;
static const StubHLS::TR chosenRofZ_digi = chosenRofZ*rMult;

static const float bField = 3.8112;
static const float invPtToInvR = bField*(2.9979e8/1.0e11);
static const float minPt_HT = 3.; // Range of Hough transform
static const float invRmin_HT = invPtToInvR*(1./minPt_HT);
static const float ptCut_4lay = 2.95; // Slightly smaller cut to allow for resolution during KF fit.
static const float ptCut_2lay = 2.90; // Slightly smaller cut to allow for resolution during KF fit.
static const float inv2Rcut_4lay = 0.5*invPtToInvR*(1./ptCut_4lay);
static const float inv2Rcut_2lay = 0.5*invPtToInvR*(1./ptCut_2lay);

static const float kalmanMultScatTerm = 0.00075; // Same as cfg param of same name in CMSSW TMTT code.

// Phi sectors
static const float TWO_PI = 2*3.14159265;
static const int numPhiSectors = 18;
static const float phiSectorWidth = TWO_PI / numPhiSectors;

// Bit shift *_bitShift to calculate HT cell from digitized (phi, invR) of helix params.
// Chosen such that pow(2,+-shift) = (dcBin_digi, dmBin_digi) calculated below.
// (where sign diff is because in KalmanUpdate.cc, one is used to shift right & one to shift left).
// Adjust if you change bits assigned to stubs.
enum {phiToCbin_bitShift = 7, invRToMbin_bitShift = 4};
enum {BCH=BH1-phiToCbin_bitShift, BMH=BH0+invRToMbin_bitShift};

// Size of HT cells
static const int numPhiBins = 64; // Bins in HT
static const int numPtBins = 36;  // Nonants
//static const int numPtBins = 32;  // Octants
static const AP_INT(BCH) minPhiBin = -numPhiBins/2; // BCH & BMH should be larger than BHT_C & BHT_M to monitor overflow.
static const AP_INT(BCH) maxPhiBin =  numPhiBins/2 - 1;
static const AP_INT(BMH) minPtBin  = -numPtBins/2;
static const AP_INT(BMH) maxPtBin  =  numPtBins/2 - 1;

/*
static const float dcBin = numPhiBins / phiSectorWidth; 
static const float dmBin = numPtBins / (invRmin_HT); 
static const float dcBin_digi = dcBin/phiMult; // = 1 / pow(2,7)
static const float dmBin_digi = dmBin/inv2R_Mult; // = pow(2,2+BEX)
*/

// Eta sector boundaries in z at reference radius (assumed symmetric).
// (As this is complex, ROM initialization fails unless stored in a class ...)

class EtaBoundaries {
public:
  enum {nSec=9};

  EtaBoundaries() {
    static const float eta[nSec+1] = {0.0, 0.31, 0.61, 0.89, 1.16, 1.43, 1.7, 1.95, 2.16, 2.4};
    for (unsigned int i = 0; i <= nSec; i++) {
      float zAtRefR = chosenRofZ/tan(2 * atan(exp(-eta[i])));
      z_[i] = rMult*zAtRefR;
    }
  }

public:
  StubHLS::TZ  z_[nSec+1];
};

// Also failed in VHDL
//static const EtaBoundaries etaBoundaries;

// Cuts to select acceptable track states.
static const KFstateHLS::TZ z0_digi_cut = beamSpotLength*rMult;  // r multiplier used for z in KF.
static const KFstateHLS::TR inv2R_digi_cut_4lay = inv2Rcut_4lay*inv2R_Mult;
static const KFstateHLS::TR inv2R_digi_cut_2lay = inv2Rcut_2lay*inv2R_Mult;
// 1/2R cut for different #stubs on track. Element 0 is never used.
static const KFstateHLS::TR inv2R_digi_cut[] = {0, 0, inv2R_digi_cut_2lay, inv2R_digi_cut_2lay, inv2R_digi_cut_4lay, inv2R_digi_cut_4lay, inv2R_digi_cut_4lay, inv2R_digi_cut_4lay};

// Chi2 cut for different #stubs on track. Taken from KF*ParamsComb::isGoodState(). Element 0 is never used.
//static const KFstateHLS::TCHI chi2_digi_cut[] = {1023, 1023, 15, 100, 320, 1023, 1023, 1023}; 
static const KFstateHLS::TCHI chi2_digi_cut[] = {1023, 1023, 10, 30, 80, 120, 160, 1023}; 

#ifdef CMSSW_GIT_HASH
}
}
#endif

#endif
