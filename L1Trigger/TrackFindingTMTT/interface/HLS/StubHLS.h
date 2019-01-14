#ifndef __StubHLS__
#define __StubHLS__

/**
 * This defines Stubs for the Kalman Filter HLS code.
 * N.B. It therefore can't use the Settings class or any external libraries! Nor can it be a C++ class.
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */

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

// Total number of bits taken from https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/stubs.vhd

// But updated to include 2 extra bits for stub (r,z) & factor 4 increase in seed tanL uncertainty.

// Number of bits needed for integer part of stub info.
// (For r & z, this is 1 larger than HT output format, since KF VHDL internally redigitizes them, measuring
//  r w.r.t. beamline & uses r digi multipler for z. Except that in nonant data format, KF VHDL actually
// uses z multiplier that is factor 2 smaller than r multiplier. By using BSZ1 = BSZ - 1 below,
// the additional factor 2 is applied at the VHDL-HLS interface. ).

enum B_STUB {BSR = 12+1, BSZ = 14+1, BSP=14, BSZ1 = BSZ - 1}; // nonants

class StubHLS {
public:
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

#ifdef PRINT_HLSARGS
public:
  void print(const char* text) const {
    std::cout<<text<<" r="<<r<<" phiS="<<phiS<<" z="<<z/2<<std::endl;
  }
#endif
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
