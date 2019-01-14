/**
 * This defines the KF matrices and the operations performance on them.
 *
 *  All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 *
 * Author: Ian Tomalin
 */
 
#ifndef __KalmanMatricesHLS__
#define __KalmanMatricesHLS__

// Defines StateHLS & KFstateHLS. Also defines finite bit integers & floats.
#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSconstants.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#include "HLSconstants.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Forward declarations
class InvR;
class VectorM;
class MatrixV;
template<unsigned int NPAR> class BODGE;
template<unsigned int NPAR> class MatrixInverseR;
template<unsigned int NPAR> class MatrixH;
template<unsigned int NPAR> class MatrixS;
template<unsigned int NPAR> class MatrixS_transpose;
template<unsigned int NPAR> class MatrixR;
template<unsigned int NPAR> class MatrixK;
template<unsigned int NPAR> class VectorRes;
template<unsigned int NPAR> class VectorX;
template<unsigned int NPAR> class MatrixC;

// These hard-wired constants reduce the number of bits relative to that predicted by bit counting
// based on variable ranges observed in CMSSW. So they should be checked if the code is changed.
// (The dependence of one bodge on another is chosen to ensure that changing each affects the declaration
// of only one variable).

enum BIT_ADJUST_PLAY {BODGE_V=5, BODGE_CHI2=12};

template<>
class BODGE<4> {
public:
  enum BIT_ADJUST {V=BODGE_V, S=2+1, R=5-S, IR=10+R+S, DET=15-2*R-2*S, K=23-IR+R, RES=2, CHI2=BODGE_CHI2};
};

template<>
class BODGE<5> {
public:
  enum BIT_ADJUST {V=BODGE_V, S=2+2, R=5-S+2, IR=10+R+S-3, DET=15-2*R-2*S+2, K=23-IR+R-2, RES=2, CHI2=BODGE_CHI2};
};

// Allow some bits for correlation in helix params between r-phi and r-z planes.
enum RZ_PHI_CORR {BCORR=1}; // Set this to 1, so this correlation is almost neglected.

// Vector of stub coords.

class VectorM {
public:
  typedef AP_FIXED(BSP+1,BSP) TMP;
  typedef AP_FIXED(BSZ1+1,BSZ) TMZ;

  VectorM(const StubHLS::TP& phiS, const StubHLS::TZ& z) : _0(phiS), _1(z) {
    // Compensate for stubs being rounded down when digitized, by adding extra LSB to coords set to '1'.
    _0[0] = 1;
    _1[0] = 1;
  } 
public:
  //  StubHLS::TP _0;
  //  StubHLS::TZ _1;  
  TMP _0;
  TMZ _1;  
};

// Covariance matrix of stub coords.

class MatrixV {
public:
  // Use same granularity for resolution as for residuals.

  // But with integer range reduced by BODGE_V, as hit resolution much smaller than max. stub coordinate.
  enum {BVP=BSP-BODGE_V, BVZ=BSZ-BODGE_V, BVPP=2*BVP, BVZZ=2*BVZ};
  typedef AP_UFIXED(B34,BVPP) TVPP;
  typedef AP_UFIXED(B34,BVZZ) TVZZ;
  typedef AP_UFIXED(1,1)       TV0;

  enum {BM=12}; // Used for pitch. May need increasing if larger r or phi multipliers used.
  typedef AP_UFIXED(B17,BM) TM;

public:
  MatrixV(const StubHLS::TR& r, const StubHLS::TZ& z, const KFstateHLS<5>::TR& inv2R, const KFstateHLS<5>::TT& tanL, const AP_INT(BHT_M)& mBin);

public:
  TVPP _00;
  TVZZ _11;
  const TV0  _01;
  const TV0& _10; // Matrix symmetric so use reference to top half.

  // Record if stub is in 2S module or not. (Not pretty to include it in this class, but convenient ...)
  bool _2Smodule;
};

// Utility for calculating pow(pitch/r, 2) from ROM. 

class PitchOverR_2 {
public:
  enum {BRED = 4}; // To save ROM resources, reduce granularity in r by this number of bits. 
  enum {MAXN = 1 << (BSR - BRED)}; // pow(2,BSR) // Max. value of [r / pow(2,BRED)].
  // Number of bits chosen based on CalcCheck job summary.
  typedef AP_UFIXED(12,5)   TPOR;
public:

  PitchOverR_2(const MatrixV::TM& pitch) {
    for (unsigned int n = 2; n < MAXN; n++) { // Don't bother initializing first two, as would overflow bit range.
      float pitchOverR = float(pitch)/(float((n << BRED)) + 0.5*(1 << BRED));
      get[n] = pitchOverR * pitchOverR; // Round to nearest half integer
    }
  }
public:
  TPOR get[MAXN];
};

// Utility for estimating pow(const/Pt, 2) from ROM, where 1/Pt is taken from HT m-bin.

class InvPt2 {
public:
  InvPt2(const AP_UINT(1)& iOpt = 0, float scaleFactor = 1.) {
    float theConst;
    if (iOpt == 0) {
      // Used to estimate increase in phi stub error due to scattering.
      theConst = kalmanMultScatTerm*phiMult;
    } else {
      // Used to estimate increase in phi stub error due to conversion from (r,phi) to (z,phi) in endcap.
      theConst = 0.5*invPtToInvR*scaleFactor*inv2R_Mult;
    }
    for (int m = minPtBin; m <= maxPtBin; m++) { 
      // Estimate Pt from Hough m-bin cell.
      float constOverPt = theConst * (1./minPt_HT)*(float(m) + 0.5)/float(numPtBins/2);
      get[m - minPtBin] = constOverPt*constOverPt;
    }
  }

  const MatrixV::TVPP& getIt(const AP_INT(BHT_M)& m) const {return this->get[m - minPtBin];} 

public:
  MatrixV::TVPP get[numPtBins];
};

// Utility for calculating 1/(unsigned int), where unsigned int has MaxInverseR::BDET bits and the
// most significant of these bits has value 1. (Loaded into ROM).

class OneOverInt {
public:
  enum {BDET=9}; // Number of significant bits used to calculate 1/determinant. Keep small to save resources.
  enum {MINN = (1 << (BDET - 1)), MAXN = (1 << BDET)}; // pow(2,BSR) // Min & max. value of r
  enum {BOI=2-BDET}; // Large enough to contain reciprocal.
  typedef SW_UFIXED(BDET,BOI)   TOI;
public:
  OneOverInt() {
    for (unsigned int n = MINN; n < MAXN; n++) { // Don't bother initializing first two, as would overflow bit range.
      get[n - MINN] = 1./float(n); // Round to nearest half integer
    }
  }

  const TOI& getIt(const AP_UINT(BDET)& i) const {return get[i - MINN];}

private:
  TOI get[MAXN - MINN + 1];
};


// Inverse of matrix R. 

template <unsigned int NPAR>
class MatrixInverseR {
public:
  // Calculate number of integer bits required for elements of R (assuming matrix approximately diagonal).
  enum {BIR00=B34 - MatrixR<NPAR>::BR00 - BODGE<NPAR>::IR,
        BIR11=B34 - MatrixR<NPAR>::BR11 - BODGE<NPAR>::IR,
        BIR01=0};   // correlation between r-phi & r-z planes not used.
  typedef SW_UFIXED(B34,   BIR00)  TRI00;
  typedef SW_UFIXED(B34,   BIR11)  TRI11;
  typedef SW_UFIXED(BCORR, BIR01)  TRI01;

  // Additional types used to cast this matrix to a lower precision one for use in chi2 calculation.
  typedef SW_UFIXED(B26,   BIR00) TRI00_short;
  typedef SW_UFIXED(B26,   BIR11) TRI11_short;
  typedef SW_UFIXED(BCORR, BIR01) TRI01_short;

  enum {BDET=OneOverInt::BDET}; // Number of significant bits used to calculate 1/determinant. Keep small to save resources.

public:
  MatrixInverseR(const MatrixR<NPAR>& R);

public:
  // R(inverse)
  TRI00  _00;
  TRI11  _11;
  TRI01  _01;
  TRI01& _10;  // Matrix symmetric, so can use reference.
};

// Since chi2 can be large, use more bits for internal calculation than for external number.
typedef AP_UFIXED(B17+BODGE_CHI2,BCHI+BODGE_CHI2) TCHI_INT;

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif




