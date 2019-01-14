#ifndef __KFstateHLS__
#define __KFstateHLS__

/**
 * This defines Helix States for the Kalman Filter HLS code.
 * N.B. It therefore can't use the Settings class or any external libraries! Nor can it be a C++ class.
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */

// Copied from /opt/ppd/tools/xilinx/Vivado_HLS/2016.4/include/
#include "ap_int.h"
#include "ap_fixed.h"

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSpragmaOpts.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#else
#include "HLSpragmaOpts.h"
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
// Though if multiplying unsigned variables, must use 1 bit less than this.

enum B_DSP {
  // Number of bits used by DSP for multiplication of signed numbers in FPGA.
#ifdef ULTRASCALE
  B18=18, B27=27, B35=2*B18-1, B48=48,
#else
  B18=18, B27=27-2, B35=2*B18-1, B48=48,
#endif
  // Number of bits used by DSP for multiplication of unsigned numbers in FPGA.
  B17=B18-1, B26=B27-1, B34=B35-1,
  // Number of bits used for interface to VHDL (historic, but increasing may increase VHDL BRAM use).
  B25=25, B24=B25-1
};

// Total number of bits from  https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/KFstates.vhd .
// Fractional number of bits from dfeFixMax() or dfeFix() in Maxeller code https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/cgn/src/formats/State.maxj
// Can change if modify seed covariance in VHDL accordingly.

// But updated to include 2 extra bits for stub (r,z) & factor 4 increase in seed tanL uncertainty.

// Number of bits for helix params & chi2.
enum B_HELIX {BH0 = 3, BH1 = 15, BH2 = 5, BH3 = 11, BH4=25, BCHI = 10};
// Number of bits needed for integer part of helix covariance matrix & their sign.
enum B_HCOV {BHC00 = -5, BHC11 = 17, BHC22 = 0, BHC33=17, BHC44=42+2, BHC01=6, BHC23=8, BHC04=18, BHC14=20+8+2};
// Total number of bits needed for off-diagonal elements of helix covariance matrix.
#ifdef COV_EXTRA_BITS 
enum B_HCOV_OPT {B18or25 = B25};
#else
enum B_HCOV_OPT {B18or25 = B18};
#endif

// Bits used for Hough (m,c) bins.
enum B_HT {BHT_C = 6, BHT_M = 6};   // Nonants

// Format of KF helix state to match VHDL, for both 4 & 5 param helix states.

template <unsigned int NPAR> class KFstateHLS;

template <> class KFstateHLS<4> {
public:
  typedef AP_FIXED(B18,BH0) TR;   
  typedef AP_FIXED(B18,BH1) TP;   
  typedef AP_FIXED(B18,BH2) TT;   
  typedef AP_FIXED(B18,BH3) TZ;   

  typedef AP_FIXED(B25,BHC00) TC00; // Seems silly that this is signed with 25 bits, rather than unsigned with 24.
  typedef AP_FIXED(B25,BHC11) TC11;
  typedef AP_FIXED(B25,BHC22) TC22;
  typedef AP_FIXED(B25,BHC33) TC33;
  typedef AP_FIXED(B18or25,BHC01) TC01;
  typedef AP_FIXED(B18or25,BHC23) TC23;

  typedef AP_UFIXED(B17,BCHI) TCHI;

  // The digitized helix & covariance parameters specified here are scaled relative to the floating 
  // point ones by factors appearing in KF4ParamsCombHLS::getDigiState().

  AP_INT(BHT_C) cBin;     // The HT cell (cbin, mbin) are centred on zero here.
  AP_INT(BHT_M) mBin;     
  TR  inv2R; // This is misnamed as rInv in Maxeller. Integer bits = 1+ceil(log2(51));
  TP  phi0;  // Measured with respect to centre of phi sector. Integer bits = 1+ceil(log2(8191));
  TT  tanL;  // This is misnamed as tanTheta in Maxeller. Integer bits = 1+ceil(log2(12));
  TZ  z0;    // Integer bits = 1+ceil(log2(150));

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
  AP_UINT(3)  eventID;        // Not used by KF updator. Just helps VHDL keep track of which event this is.
  AP_UINT(4)  etaSectorID; // Eta sector ID, but counting away from 0 near theta=PI/2 & increasing to 8 near endcap. (Named SectorID in Maxeller).
  AP_UINT(1)  etaSectorZsign;  // True if eta sector is in +ve z side of tracker; False otherwise. (Named zSign in Maxeller).
  AP_UINT(1)  valid; // Used by external code when calculation finished on valid input state & stub.

#ifdef PRINT_HLSARGS
public:
  void print(const char* text) const {
    std::cout<<text
	   <<" HT (m,c)=("<<mBin<<","<<cBin<<")"
           <<" layers (ID, skip)=("<<layerID<<","<<nSkippedLayers<<")"
           <<" 1/2R="<<ap_fixed<B18,B18>(inv2R.range( B18 - 1, 0))
	   <<" phi0="<<ap_fixed<B18,B18>(phi0.range( B18 - 1, 0))
	   <<" tanL="<<ap_fixed<B18,B18>(tanL.range( B18 - 1, 0))
	   <<" z0="  <<ap_fixed<B18,B18>(z0.range( B18 - 1, 0))
	   <<" chi2="<<ap_ufixed<B17,B17>(chiSquared.range( B17 - 1, 0))
	   <<std::endl;
    std::cout<<text
           <<" cov00="<<ap_fixed<B25,B25>(cov_00.range( B25 - 1, 0))
           <<" cov11="<<ap_fixed<B25,B25>(cov_11.range( B25 - 1, 0))
           <<" cov22="<<ap_fixed<B25,B25>(cov_22.range( B25 - 1, 0))
           <<" cov33="<<ap_fixed<B25,B25>(cov_33.range( B25 - 1, 0))
           <<" cov01="<<ap_fixed<B18or25,B18or25>(cov_01.range( B18or25 - 1, 0))
           <<" cov23="<<ap_fixed<B18or25,B18or25>(cov_23.range( B18or25 - 1, 0))
	   <<std::endl;
  }
#endif
};


template <> class KFstateHLS<5> : public KFstateHLS<4> {
public:
  typedef AP_FIXED(B18,BH4) TD;   
  typedef AP_FIXED(B25,BHC44) TC44;
  typedef AP_FIXED(B18or25,BHC04) TC04;
  typedef AP_FIXED(B18or25,BHC14) TC14;

public:
  TD  d0;

  TC44 cov_44; // (d0,    d0)   
  TC04 cov_04; // (inv2R, d0)   -- other off-diagonal elements assumed negligible.
  TC14 cov_14; // (phi0,  d0)   -- other off-diagonal elements assumed negligible.

#ifdef PRINT_HLSARGS
public:
  void print(const char* text) const {
    this->KFstateHLS<4>::print(text);
    std::cout<<text
             <<" d0="<<ap_fixed<B18,B18>(d0.range( B18 - 1, 0))
             <<" cov44="<<ap_fixed<B25,B25>(cov_44.range( B25 - 1, 0))
             <<" cov04="<<ap_fixed<B18or25,B18or25>(cov_04.range( B18or25 - 1, 0))
             <<" cov14="<<ap_fixed<B18or25,B18or25>(cov_14.range( B18or25 - 1, 0))
             <<std::endl;
  }
#endif
};

// Additional output parameters returned by KF updated, for both 4 & 5 param helix fits.
//https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/KalmanFilter/KalmanWorker.vhd?peg=4914

template <unsigned int NPAR> class ExtraOutHLS;

template <> class ExtraOutHLS<4> {
public:
  // Must use AP_UINT(1) instead of bool, due to bug in HLS IP export.
  AP_UINT(1)    z0Cut; // Did updated state pass cut on z0 etc.
  AP_UINT(1)    ptCut;
  AP_UINT(1)    chiSquaredCut;
  AP_UINT(1)    sufficientPScut; // Enough PS layers
  AP_UINT(1)    htBinWithin1Cut;  
  AP_INT(BHT_C) mBinHelix;    // HT bin that fitted helix params lie within.
  AP_INT(BHT_M) cBinHelix;
  AP_UINT(1)    sectorCut;   // Helix parameters lie within Sector.
  AP_UINT(1)    consistent;  // Duplicate removal -- helix parameters lie within original HT cell.

#ifdef PRINT_HLSARGS
public:
  void print(const char* text) const {
    std::cout<<"HLS OUTPUT EXTRA:"
             <<" Helix (m,c)=("<<mBinHelix<<","<<cBinHelix<<")"
    	     <<std::endl;
  }
#endif
};

template <> class ExtraOutHLS<5> : public ExtraOutHLS<4> {
public:
  AP_UINT(1)    d0Cut;
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
