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
class VectorX;
class MatrixC;
class MatrixS; 
class MatrixR;
class MatrixInverseR;
class MatrixK;
class VectorRes;

// These hard-wired constants reduce the number of bits relative to that predicted by bit counting
// based on variable ranges observed in CMSSW. So they should be checked if the code is changed.
// (The dependence of one bodge on another is chosen to ensure that changing each affects the declaration
// of only one variable).
//enum BIT_ADJUST {BODGE_V=5, BODGE_S=2+1, BODGE_R=5-BODGE_S, BODGE_IR=10+BODGE_R+BODGE_S, BODGE_DET=14-2*BODGE_R-2*BODGE_S, BODGE_K=24-BODGE_IR+BODGE_R, BODGE_RES=2, BODGE_CHI2=12};
// Change to allow increased seed uncertainty in tanL.
//enum BIT_ADJUST {BODGE_V=5, BODGE_S=2+1, BODGE_R=5-BODGE_S, BODGE_IR=10+BODGE_R+BODGE_S, BODGE_DET=14-2*BODGE_R-2*BODGE_S, BODGE_K=23-BODGE_IR+BODGE_R, BODGE_RES=2, BODGE_CHI2=12};
// Change to avoid MSB=0 error during matrix inversion.
enum BIT_ADJUST {BODGE_V=5, BODGE_S=2+1, BODGE_R=5-BODGE_S, BODGE_IR=10+BODGE_R+BODGE_S, BODGE_DET=15-2*BODGE_R-2*BODGE_S, BODGE_K=23-BODGE_IR+BODGE_R, BODGE_RES=2, BODGE_CHI2=12};

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
  MatrixV(const StubHLS::TR& r, const StubHLS::TZ& z, const KFstateHLS::TR& inv2R, const KFstateHLS::TT& tanL, const AP_INT(BHT_M)& mBin);

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

// Calculate matrix of derivatives of predicted stub coords w.r.t. helix params.

class MatrixH {
public:
  enum {BH=BSR+1};
  typedef AP_FIXED(BH,BH)  TH;  // One extra bit, since "-r" can be -ve.
  typedef AP_UFIXED(1,1)   T1;
  MatrixH(const StubHLS::TR& r) : _00(-r), _12(r),
                                           _01(1), _02(0), _03(0),
                                  _10(0),  _11(0),         _13(1) {}
public:
  TH       _00, _12;
  const T1      _01, _02, _03, 
           _10, _11,      _13;
};

// S = H * C

class MatrixS {
public:
  enum {BH=MatrixH::BH,
  // Calculate number of integer bits required for all non-zero elements of S.
  // (Assumes that some elements of C & H are zero and that all non-trivial elements of H have BH bits).
        BS00=MAX2(BH+BHC00, BHC01) - BODGE_S,  // H00*C00 + H01*C10 + (H02*C20 + H03*C30 = zero).
        BS01=MAX2(BH+BHC01, BHC11) - BODGE_S,  // H00*C01 + H01*C11 + (H02*C21 + H03*C31 = zero).
        BS12=MAX2(BH+BHC22, BHC23) - BODGE_S,  // (H00*C02 + H01*C12 = zero) + H02*C22 + H03*C32.
	BS13=MAX2(BH+BHC23, BHC33) - BODGE_S,  // (H00*C03 + H01*C13 = zero) + H02*C23 + H03*C33.
        BS=0}; // Neglect correlation between r-phi & r-z planes for now.
  typedef AP_FIXED(B25,BS00)  TS00;
  typedef AP_FIXED(B25,BS01)  TS01;
  typedef AP_FIXED(B25,BS12)  TS12;
  typedef AP_FIXED(B25,BS13)  TS13;
  typedef AP_FIXED(BCORR,BS)  TS; 

public:
  MatrixS(const MatrixH& H, const MatrixC& C);

public:
  
  TS00 _00;
  TS01 _01;
  TS12 _12;
  TS13 _13;
  TS            _02, _03,
      _10, _11          ;
};

// S(transpose) = C*H(transpose)

class MatrixS_transpose {
public:
  typedef MatrixS::TS00  TS00;
  typedef MatrixS::TS01  TS01;
  typedef MatrixS::TS12  TS12;
  typedef MatrixS::TS13  TS13;
  typedef MatrixS::TS    TS;
  MatrixS_transpose(const MatrixS& S) : _00(S._00), _10(S._01), _21(S._12), _31(S._13),
					_01(S._10), _11(S._11), _20(S._02), _30(S._03) {}
public:
  const TS00&  _00;
  const TS01&  _10;
  const TS12&  _21;
  const TS13&  _31;
  const TS&       _01,
                  _11,
             _20,
             _30     ;
};

// Covariance matrix of predicted residuals R = V + H*C*Ht = V + H*St.

class MatrixR {
public:
  enum {BH=MatrixH::BH,
	BS00=MatrixS::BS00,
	BS01=MatrixS::BS01,
	BS12=MatrixS::BS12,
	BS13=MatrixS::BS13,
	BS  =MatrixS::BS,
        // Calculate number of integer bits required for elements of R.
        BR00 = MAX2(MatrixV::BVPP, MAX2(BH+BS00, BS01)) - BODGE_R, // H00*St00 + H01*St10 + (H02*St20 + H03*St30 = zero)
	BR11 = MAX2(MatrixV::BVZZ, MAX2(BH+BS12, BS13)) - BODGE_R, // (H10*St01 + H11*St11 = zero) + H12*St21 + H13*St31
	BR01 = MAX2(BH+BS, BS)                                     // H00*St01 + H01*St11 + (H02*St21 + H03*St31 = zero)
       };  
  typedef SW_UFIXED(B34,BR00) TR00;
  typedef SW_UFIXED(B34,BR11) TR11;
  typedef SW_UFIXED(B34,BR01) TR01;

public:
  MatrixR(const MatrixV& V, const MatrixH& H, const MatrixS_transpose& St);

public:
  TR00  _00;
  TR11  _11;
  TR01  _01;
  TR01& _10; // Matrix symmetric, so can use reference.
};

// Inverse of matrix R. 

class MatrixInverseR {
public:
  // Calculate number of integer bits required for elements of R (assuming matrix approximately diagonal).
  
  enum {BIR00=B34 - MatrixR::BR00 - BODGE_IR,
	BIR11=B34 - MatrixR::BR11 - BODGE_IR,
        BIR01=0};   // correlation between r-phi & r-z planes not used.
  typedef SW_UFIXED(B34,BIR00)  TRI00;
  typedef SW_UFIXED(B34,BIR11)  TRI11;
  typedef SW_UFIXED(B34,BIR01)  TRI01;

  // Additional types used to cast this matrix to a lower precision one for use in chi2 calculation.
  typedef SW_UFIXED(B25, BIR00) TRI00_short;
  typedef SW_UFIXED(B25, BIR11) TRI11_short;
  typedef SW_UFIXED(B25, BIR01) TRI01_short;

  // IRT - increase to cope with larger tanL seed uncertainty.
  //enum {BDET=8}; // Number of significant bits used to calculate 1/determinant. Keep small to save resources.
  enum {BDET=9}; // Number of significant bits used to calculate 1/determinant. Keep small to save resources.

public:
  MatrixInverseR(const MatrixR& R);

public:
  // R(inverse)
  TRI00  _00;
  TRI11  _11;
  TRI01  _01;
  TRI01& _10;  // Matrix symmetric, so can use reference.
};

// Utility for calculating 1/(unsigned int), where unsigned int has MaxInverseR::BDET bits and the
// most significant of these bits has value 1. (Loaded into ROM).

class OneOverInt {
public:
  enum {BDET=MatrixInverseR::BDET};
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

// Kalman gain matrix K = S(transpose)*R(inverse).

class MatrixK {
public:
  enum {BS00=MatrixS::BS00,
	BS01=MatrixS::BS01,
	BS12=MatrixS::BS12,
	BS13=MatrixS::BS13,
        BIR00=MatrixInverseR::BIR00,
        BIR11=MatrixInverseR::BIR11,
        BIR01=MatrixInverseR::BIR01,
        BK00=(BS00+BIR00) - BODGE_K,  // St00*Rinv00 (+ St01*Rinv10 = zero)
        BK10=(BS01+BIR00) - BODGE_K,  // St10*Rinv00 (+ St11*Rinv10 = zero)
        BK21=(BS12+BIR11) - BODGE_K,  // (St20*Rinv01 = zero) + St21*Rinv11
        BK31=(BS13+BIR11) - BODGE_K}; // (St30*Rinv01 = zero) + St31*Rinv11
  typedef SW_FIXED(B35,BK00)  TK00;
  typedef SW_FIXED(B35,BK10)  TK10;
  typedef SW_FIXED(B35,BK21)  TK21;
  typedef SW_FIXED(B35,BK31)  TK31;
  typedef SW_FIXED(BCORR,0)   TK; // Neglect correlation between r-phi & r-z
  MatrixK(const MatrixS_transpose& St, const MatrixInverseR& RmatInv);
public:
  // Additional types used to cast this matrix to a lower precision one for updated helix param calculation.
  typedef SW_FIXED(B25,BK00)  TK00_short;
  typedef SW_FIXED(B25,BK10)  TK10_short;
  typedef SW_FIXED(B25,BK21)  TK21_short;
  typedef SW_FIXED(B25,BK31)  TK31_short;
public:
  TK00  _00;
  TK10  _10;
  TK21      _21;
  TK31      _31;
  TK        _01,
            _11,
        _20,
        _30    ;
};

// Hit residuals: res = m - H*x. 

class VectorRes {
public:
  // Use higher granularity for residuals than for stubs.
  // BODGE_RES should be slightly larger than BODGE_V as hits can be several sigma wrong.
  // Add one extra fractional bit relative to stub, to avoid additional precision loss.
  typedef AP_FIXED(B18-BODGE_RES+1,BSP-BODGE_RES) TRP;
  typedef AP_FIXED(B18-BODGE_RES+1,BSZ-BODGE_RES) TRZ;

public:
  VectorRes(const VectorM& m, const MatrixH& H, const VectorX& x);

public:
  TRP _0;
  TRZ _1;
};


// Vector of helix params.

class VectorX {
public:
  // Determine input helix params.
  VectorX(const KFstateHLS::TR& inv2R, const KFstateHLS::TP& phi0, const KFstateHLS::TT& tanL, const KFstateHLS::TZ& z0) : _0(inv2R), _1(phi0), _2(tanL), _3(z0) {} 

  // Calculate output helix params: x' = x + K*res
  VectorX(const VectorX& x, const MatrixK& K, const VectorRes& res);

public:
  KFstateHLS::TR _0;
  KFstateHLS::TP _1;  
  KFstateHLS::TT _2;  
  KFstateHLS::TZ _3;  
};

// Covariance matrix of helix params.

class MatrixC {
public:
  typedef AP_UFIXED(1,1)      T0; // HLS doesn't like zero bit variables.

  // Determine input helix coviaraiance matrix.
  MatrixC(const KFstateHLS::TC00& c00, const KFstateHLS::TC11& c11, const KFstateHLS::TC22& c22, 
          const KFstateHLS::TC33& c33, const KFstateHLS::TC01& c01, const KFstateHLS::TC23& c23) :
             _00(c00), _11(c11), _22(c22), _33(c33), _01(c01), _23(c23), 
             _02(0), _03(0), _12(0), _13(0),
             _10(_01), _32(_23), _20(_02), _30(_03), _21(_12), _31(_13) {}

  // Calculate output helix covariance matrix: C' = C - K*H*C = C - K*S.
  MatrixC(const MatrixC& C, const MatrixK& K, const MatrixS& S);

public:
  // Elements that are finite
// Maxeller wierdly uses signed 25 bits for these, so use unsigned 24 instead to match the DSP abilities.
//  KFstateHLS::TC00 _00;
//  KFstateHLS::TC11 _11;
//  KFstateHLS::TC22 _22;
//  KFstateHLS::TC33 _33;
  AP_UFIXED(B24,BHC00-1) _00; // One less integer bit as no sign required.
  AP_UFIXED(B24,BHC11-1) _11;
  AP_UFIXED(B24,BHC22-1) _22;
  AP_UFIXED(B24,BHC33-1) _33;
  KFstateHLS::TC01 _01; // (inv2R, phi0) -- other off-diagonal elements assumed negligeable.
  KFstateHLS::TC23 _23; // (tanL,  z0)   -- other off-diagonal elements assumed negligeable.
  // Elements that are zero.
  const T0 _02, _03, _12, _13;
  // Elements below the diagonal of this symmetric matrix.
  const KFstateHLS::TC01 &_10;
  const KFstateHLS::TC23 &_32;
  const T0 &_20, &_30, &_21, &_31;
};

// Since chi2 can be large, use more bits for internal calculation than for external number.
typedef AP_UFIXED(B17+BODGE_CHI2,BCHI+BODGE_CHI2) TCHI_INT;

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif




