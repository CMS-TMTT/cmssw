#ifndef __KF4interfaceHLS__
#define __KF4interfaceHLS__

// If defined, then add 2 extra bits to stub r & z coords relative to Dec. 2016 format.
#define EXTRA_BITS

// If defined, use 25 bits for the off-diagonal elements of the helix param covariance matrix.
// Ian thinks this is needed to avoid bit overflows messing up some tracks (although the effect on tracking
// performance is small). And it would require changing the KF interface in the firmware ...
// It is unnecessary if the hit errors are inflated to allow for scattering.
//#define COV_EXTRA_BITS

// If defined, set optimum numbers of bits for Ultrascale instead of Virtex7 FPGAs.
//#define ULTRASCALE

// Must use AP_UINT(1) instead of bool, due to bug in HLS IP export.

// Copied from /opt/ppd/tools/xilinx/Vivado_HLS/2016.4/include/
#include "ap_int.h"
#include "ap_fixed.h"

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#else
#include "HLSutilities.h"
#endif

///=== Hard-wired configuration parameters.
///=== WARNING: Since this code must be used in Vivado HLS, it can't use the Settings class.
///=== Therefore all constants are hard-wired here, so may break if you change the configuration parameters.

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

///=== Data format of Stubs & KF helix state passed to KF updator as assumed within Maxeller implementation.
///=== Numbers of bits hard-wired, since same code also used within Vivado HLS, where access to Settings class not possible.

// Format of Stub taken from existing KF Maxeller firmware, keeping only useful elements.
//https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/stubs.vhd

// Virtex7 DSP = (18 bits * 25 bits = 48 bits); Ultrascale DSP = (18 bits * 27 bits = 48 bits).
// Though if multiplying unsigned variables, must use 1 bit less than this?

enum B_DSP {
  // Number of bits used by DSP for multiplication of signed numbers in FPGA.
#ifdef ULTRASCALE
  B18=18, B25=25+2, B35=2*B18-1, B48=48,
#else
  B18=18, B25=25  , B35=2*B18-1, B48=48,
#endif
  // Number of bits used by DSP for multiplication of unsigned numbers in FPGA.
  B17=B18-1, B24=B25-1, B34=B35-1
};

// Extra bits for stub (r,z) coords. -- Also influences helix digitisation.
#ifdef EXTRA_BITS
enum B_EXTRA {BEX = 2};
#else
enum B_EXTRA {BEX = 0};
#endif

// Total number of bits taken from https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/stubs.vhd

// Number of bits needed for integer part of stub info.
// (For r & z, this is 1 larger than HT output format, since KF VHDL internally redigitizes them, measuring
//  r w.r.t. beamline & uses r digi multipler for z. Except that in nonant data format, KF VHDL actually
// uses z multiplier that is factor 2 smaller than r multiplier. By using BSZ1 = BSZ - 1 below,
// the additional factor 2 is applied at the VHDL-HLS interface. ).
// For octant format
//enum B_STUB {BSR = 10+1+BEX, BSZ = 12+1+BEX, BSP=14, BSZ1 = BSZ};
// For nonant format
enum B_STUB {BSR = 10+1+BEX, BSZ = 12+1+BEX, BSP=14, BSZ1 = BSZ - 1};

// Total number of bits from  https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/KFstates.vhd .
// Fractional number of bits from dfeFixMax() or dfeFix() in Maxeller code https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/cgn/src/formats/State.maxj
// Can change if modify seed covariance in VHDL accordingly.

// Number of bits needed for integer part of helix parameters & chi2 & their sign
//enum B_HELIX {BH0 = 7-BEX, BH1 = 14, BH2 = 5, BH3 = 9+BEX, BCHI = 10};
// IRT make compatible with data format doc.
enum B_HELIX {BH0 = 5-BEX, BH1 = 15, BH2 = 5, BH3 = 9+BEX, BCHI = 10};
// Number of bits needed for integer part of helix covariance matrix & their sign.
//enum B_HCOV {BHC00 = -1-2*BEX, BHC11 = 17, BHC22 = -2, BHC33=13+2*BEX, BHC01=8-BEX, BHC23=6+BEX};
// IRT: scaled up tanL uncertainty from seed by factor 4, requiring this change.
enum B_HCOV {BHC00 = -1-2*BEX, BHC11 = 17, BHC22 = 0, BHC33=13+2*BEX, BHC01=8-BEX, BHC23=6+BEX};
// Number of bits needed for integer part of track fit chi3.

// Bits used for Hough (m,c) bins.
enum B_HT {BHT_C = 6, BHT_M = 6};   // Nonants
//enum B_HT {BHT_C = 6, BHT_M = 5};  // Octants

struct StubHLS {
  typedef AP_UFIXED(BSR,BSR) TR;
  //  typedef AP_UFIXED(BSR,BSR) TR;
  typedef AP_FIXED(BSZ1,BSZ)  TZ;
  typedef AP_FIXED(BSP,BSP)  TP;
  // The digitized stub parameters here differ slightly from those used in the HT, to simplify the KF maths.
  TR r;      // Note this is r, not rT, so is unsigned. (Misnamed in Maxeller) -- IRT: I suspect only 10 bits are needed.
  TZ z;      // This is (rMult/zMult) larger than digitised z used by HT, so needs one more bit.
  TP phiS;
  // IRT: Maxeller transmits stub layerId (not reducedLayerId). I assume this can't be used for anything, so don't do it here.
  AP_UINT(1)       valid; // Used by external code to indicate if input data is valid.
};

// Format of KF helix state existing Maxeller firmware, keeping only useful elements.

struct KFstateHLS {
  typedef AP_FIXED(B18,BH0) TR;   
  typedef AP_FIXED(B18,BH1) TP;   
  typedef AP_FIXED(B18,BH2) TT;   
  typedef AP_FIXED(B18,BH3) TZ;   

  typedef AP_FIXED(B25,BHC00) TC00; // Seems silly that this is signed with 25 bits, rather than unsigned with 24.
  typedef AP_FIXED(B25,BHC11) TC11;
  typedef AP_FIXED(B25,BHC22) TC22;
  typedef AP_FIXED(B25,BHC33) TC33;
#ifdef COV_EXTRA_BITS 
  typedef AP_FIXED(B25,BHC01) TC01;
  typedef AP_FIXED(B25,BHC23) TC23;
#else
  typedef AP_FIXED(B18,BHC01) TC01;
  typedef AP_FIXED(B18,BHC23) TC23;
#endif

  typedef AP_UFIXED(B17,BCHI) TCHI;

  // The digitized helix & covariance parameters specified here are scaled relative to the floating 
  // point ones by factors appearing in KF4ParamsCombHLS::getDigiState().

  AP_INT(BHT_C) cBin;     // The HT cell (cbin, mbin) are centred on zero here.
  AP_INT(BHT_M) mBin;     
  TR  inv2R; // This is misnamed as rInv in Maxeller. Integer bits = 1+ceil(log2(51));
  TP  phi0;  // Measured with respect to centre of phi sector. Integer bits = 1+ceil(log2(8191));
  TT  tanL;  // This is misnamed as tanTheta in Maxeller. Integer bits = 1+ceil(log2(12));
  TZ  z0;  // Integer bits = 1+ceil(log2(150));
  TC00 cov_00; 
  TC11 cov_11;
  TC22 cov_22;
  TC33 cov_33;
  TC01 cov_01; // (inv2R, phi0) -- other off-diagonal elements assumed negligible.
  TC23 cov_23; // (tanL,  z0)   -- other off-diagonal elements assumed negligible.
  TCHI chiSquared;    // No idea why Maxeller doesn't use 18 bits for this.
  // This is the KF layer that the KF updator is currently looking for stubs in, encoded by L1KalmanComb::doKF(), which in any eta region increases from 0-7 as a particle goes through each layer in turn. The KF updator in HLS/Maxeller does not incremement it.
  AP_UINT(3)  layerID;  
  // This is the number of skipped layers assuming we find a stub in the layer the KF updator is currently searched. The KF updator in HLS/Maxeller does not incremement it.
  AP_UINT(2)  nSkippedLayers;
  AP_UINT(6)  candidateID;    // Not used by KF updator. Just helps VHDL keep track of which state this is. 
  AP_UINT(4)  eventID;        // Not used by KF updator. Just helps VHDL keep track of which event this is. - 2 bits more than Andy's VHDL for no good reason?
  AP_UINT(4)  etaSectorID; // Eta sector ID, but counting away from 0 near theta=PI/2 & increasing to 8 near endcap. (Named SectorID in Maxeller).
  AP_UINT(1)  etaSectorZsign;  // True if eta sector is in +ve z side of tracker; False otherwise. (Named zSign in Maxeller).
  AP_UINT(1)  valid; // Used by external code when calculation finished on valid input state & stub.
};

// Additional output parameters returned by KF updated.
//https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/KalmanFilter/KalmanWorker.vhd?peg=4914

struct ExtraOutHLS {
  AP_UINT(1)    z0Cut; // Did updated state pass cut on z0 etc.
  AP_UINT(1)    ptCut;
  AP_UINT(1)    chiSquaredCut;
  AP_UINT(1)    sufficientPScut; // Enough PS layers
  AP_UINT(1)    mBinInRange;  // Not needed. Must be true if track passes Pt cut.
  AP_INT(BHT_C) mBinHelix;    // HT bin that fitted helix params lie within.
  AP_INT(BHT_M) cBinHelix;
  AP_UINT(1)    sectorCut;   // Helix parameters lie within Sector.
  AP_UINT(1)    consistent;  // Duplicate removal -- helix parameters lie within original HT cell.
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
