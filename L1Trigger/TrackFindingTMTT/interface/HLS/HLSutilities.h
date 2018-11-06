/**
 * General HLS utilities, not specific to KF.
 *
 * Author: Ian Tomalin
 */

#ifndef __HLSutilities__
#define __HLSutilities__

// Define macros to easily change between fixed bit & floating representations.
#define USE_FIXED // Comment out to use floating point representation.

#ifndef __SYNTHESIS__
// Switch on short debug summary printout at end of job.
#define PRINT_SUMMARY
// Switch on verbose debug printout (but never during Vivado synthesis).
//#define PRINT
#endif

#ifdef PRINT
#define PRINT_SUMMARY
#endif

#ifdef PRINT_SUMMARY
#include <iostream>
#include <string>
#include <typeinfo>
#include <map>
#endif

// Copied from /opt/ppd/tools/xilinx/Vivado_HLS/2016.4/include/
#include "ap_int.h"
#include "ap_fixed.h"
#include "hls_half.h"
// Gives access to fp_struct used to cast half to ap_fixed.
//#include "hls/utils/x_hls_utils.h"

///=== Hard-wired configuration parameters.
///=== WARNING: Since this code must be used in Vivado HLS, it can't use the Settings class.
///=== Therefore all constants are hard-wired here, so may break if you change the configuration parameters.

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Specify which type to use in floating point representation (float, double ...)
#define SW_FLOAT       float

// Define macros to easily change default options for fixed bit floats.
#define AP_OPTS        SC_TRN,SC_WRAP
//#define AP_OPTS        SC_TRN,SC_SAT_SYM

// Used where one definately wants to used fixed bit types.
#define AP_FIXED(N,I)  ap_fixed<N,I,AP_OPTS>
#define AP_UFIXED(N,I) ap_ufixed<N,I,AP_OPTS>
#define AP_INT(N)      ap_int<N>   // Not needed, but pretty to have same style as for AP_FIXED.
#define AP_UINT(N)     ap_uint<N>

// Used where one is unsure whether to use fixed bit or float types.
#ifdef USE_FIXED
#define SW_FIXED(N,I)  ap_fixed<N,I,AP_OPTS>
#define SW_UFIXED(N,I) ap_ufixed<N,I,AP_OPTS>
#define SW_INT(N)      ap_int<N>
#define SW_UINT(N)     ap_uint<N>
#else
#define SW_FIXED(N,I)  SW_FLOAT
#define SW_UFIXED(N,I) SW_FLOAT
#define SW_INT(N)      int
#define SW_UINT(N)     unsigned int
#endif

// Casting between half & ap_fixed is slow. Use ap_fixed::to_half() to convert to half.
// Use the following functions to convert from half. Only works in HLS 2017.X.
/*
template<int W, int I>
ap_ufixed<W,I,AP_OPTS> to_ap_ufixed(half x_in) {
  // Split input number into sign, mantissa & power of 2 exponent.
  fp_struct<half> x_str(x_in);
  int mantX   = x_str.mantissa();
  int iExpX   = x_str.expv();
  ap_ufixed<W,I,AP_OPTS> x_out_abs = ap_ufixed<W+I,I,AP_OPTS>(mantX) << iExpX;
  return x_out_abs;
}

template<int W, int I>
ap_fixed<W,I,AP_OPTS> to_ap_fixed(half x_in) {
  // Split input number into sign, mantissa & power of 2 exponent.
  fp_struct<half> x_str(x_in);
  bool iSignX = x_str.__signbit();
  int mantX   = x_str.mantissa();
  int iExpX   = x_str.expv();
  ap_ufixed<W-1,I-1,AP_OPTS> x_out_abs = ap_ufixed<W+I-1,I-1,AP_OPTS>(mantX) << iExpX;
  ap_fixed<W,I,AP_OPTS> x_out;
  if (iSignX) {
	x_out = - x_out_abs;
  } else {
    x_out =   x_out_abs;
  }
  return x_out;
}
*/

// Used to calculate number of required bits, (which are constants downloaded into FPGA).
#define MAX2(a,b)     ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c)   (MAX2(MAX2(a,b),c))
#define MAX4(a,b,c,d) (MAX2(MAX2(a,b),MAX2(c,d)))

#ifdef PRINT_SUMMARY 

// Detect poor precision in fixed bit calculations.

namespace CHECK_AP {

class INFO {
public:
  INFO() {}
  INFO(std::string className, int intBitsCfg, int intBitsSeenHigh, int intBitsSeenLow) : className_(className), intBitsCfg_(intBitsCfg), intBitsSeenHigh_(intBitsSeenHigh), intBitsSeenLow_(intBitsSeenLow)  {}
public:
  std::string className_;
  int         intBitsCfg_;
  int         intBitsSeenHigh_;
  int         intBitsSeenLow_;
};

class INFO_int {
public:
  INFO_int() {}
  INFO_int(int intCfgHigh, int intCfgLow, int intSeenHigh, int intSeenLow) : intCfgHigh_(intCfgHigh), intCfgLow_(intCfgLow), intSeenHigh_(intSeenHigh), intSeenLow_(intSeenLow) {}
public:
  int         intCfgHigh_;
  int         intCfgLow_;
  int         intSeenHigh_;
  int         intSeenLow_;
};

// Map containing info about variable ranges in terms of numbers of bits required to represent them.
extern std::map<std::string, INFO> apCheckMap_;

// Map containing info about integer ranges expressed as numbers,
extern std::map<std::string, INFO_int> intRangeMap_;

// Print contents of map about variable ranges at end of job.
void printCheckMap();

// Fill info for summary table & check if fixed bit calculation suffered precision loss.
// (Two versions of this function, one for calculations in float which merely notes max & min
// values of each variable for summary printout, & one for finite bit which also checks precision).

bool checkCalc(std::string varName, SW_FLOAT res_fix, double res_float, double reltol = 0.1, double tol = 0.0);

template <class C>
bool checkCalc(std::string varName, C res_fix, double res_float, double reltol = 0.1, double tol = 0.0) {
  bool OK = true;

  double res_float_abs = fabs(res_float);
  bool   res_float_sign = (res_float >= 0);
  int  intBitsSeen = (res_float_abs > 0)  ?  std::ceil(log(res_float_abs)/log(2.))  :  -99;

  std::string cNameUgly = typeid(C).name();
  std::string cName = "unknown";
  int intBitsCfg = 0;
  if (cNameUgly.find("ap_int") != std::string::npos) {
    cName = "ap_int   ";
    intBitsCfg = res_fix.width - 1; // One bit less for sign
  } else if (cNameUgly.find("ap_uint") != std::string::npos) {
    cName = "ap_uint  ";
    intBitsCfg = res_fix.width;
  } else if (cNameUgly.find("ap_fixed") != std::string::npos) {
    cName = "ap_fixed ";
    intBitsCfg = res_fix.iwidth - 1;
  } else if (cNameUgly.find("ap_ufixed") != std::string::npos) {
    cName = "ap_ufixed";
    intBitsCfg = res_fix.iwidth;
  }

  if (apCheckMap_.find(varName) == apCheckMap_.end()) {
    apCheckMap_[varName] = INFO(cName, intBitsCfg, intBitsSeen, intBitsSeen); 
  } else {
    int intHighOld = apCheckMap_[varName].intBitsSeenHigh_;
    int intLowOld  = apCheckMap_[varName].intBitsSeenLow_;
    if (intHighOld < intBitsSeen) {
	apCheckMap_[varName] = INFO(cName, intBitsCfg, intBitsSeen, intLowOld); 
    } else if (intLowOld > intBitsSeen) {
      apCheckMap_[varName] = INFO(cName, intBitsCfg, intHighOld, intBitsSeen); 
    }
  }

#ifdef PRINT_SUMMARY

#ifdef PRINT
#define NPRINTMAX 99999
#else
#define NPRINTMAX 100
#endif

  static unsigned int nErrors = 0;

  // Check -ve numbers aren't stored in unsigned variables.
  if (res_float < 0) {
    if (cName == "ap_uint" || cName == "ap_ufixed") {
      OK = false;
      if (nErrors < NPRINTMAX) std::cout<<"checkCalc SIGN ERROR ("<<cName<<"): "<<varName<<" "<<res_float<<std::endl;
    }
  }

  if (varName == "RES0" && intBitsSeen == 13) {
    if (nErrors < NPRINTMAX) std::cout<<"checkCalc WIERD RES0 ERROR "<<res_float<<" "<<res_fix<<std::endl;  
  }

  // Check precision is OK.
  if (res_float != 0.) {
    double err = fabs(double(res_fix) - res_float);
    double relerr = err/fabs(res_float);
    if (relerr > reltol && err > tol) {
      OK = false;
      if (nErrors < NPRINTMAX) std::cout<<"checkCalc PRECISION LOSS ("<<cName<<"): "<<varName<<" relerr="<<relerr<<" err="<<err<<" fix="<<res_fix<<" flt="<<res_float<<std::endl;
      if (intBitsSeen > intBitsCfg) if (nErrors < NPRINTMAX) std::cout<<"checkCalc TOO FEW INT BITS ("<<cName<<"): "<<varName<<" "<<intBitsSeen<<" > "<<intBitsCfg<<std::endl;
    }
  }
  if (not OK) nErrors++;
#endif
  return OK;
} 

// Fill info for integer range summary table.
bool checkIntRange(std::string varName, int intCfgHigh, int intCfgLow, int intValue);

// Check covariance matrix determinants are positive.
bool checkDet(std::string matrixName, double m00, double m11, double m01);

}
#endif

#ifdef CMSSW_GIT_HASH
}
}
#endif

#endif
