/**
 * This defines the KF matrices and the operations performance on them.
 *
 *  All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 *
 * Author: Ian Tomalin
 */
 
#ifndef __KalmanMatricesHLS5__
#define __KalmanMatricesHLS5__

// Defines StateHLS & KFstateHLS. Also defines finite bit integers & floats.
#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSconstants.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS4.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#include "HLSconstants.h"
#include "KalmanMatricesHLS.h"
#include "KalmanMatricesHLS5.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Calculate matrix of derivatives of predicted stub coords w.r.t. helix params.

template <>
class MatrixH<5> {
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

template <>
class MatrixS<5> {
public:
  enum {BH=MatrixH<5>::BH,
  // Calculate number of integer bits required for all non-zero elements of S.
  // (Assumes that some elements of C & H are zero and that all non-trivial elements of H have BH bits).
        BS00=MAX2(BH+BHC00, BHC01) - BODGE<5>::S,  // H00*C00 + H01*C10 + (H02*C20 + H03*C30 = zero).
        BS01=MAX2(BH+BHC01, BHC11) - BODGE<5>::S,  // H00*C01 + H01*C11 + (H02*C21 + H03*C31 = zero).
        BS12=MAX2(BH+BHC22, BHC23) - BODGE<5>::S,  // (H00*C02 + H01*C12 = zero) + H02*C22 + H03*C32.
	BS13=MAX2(BH+BHC23, BHC33) - BODGE<5>::S,  // (H00*C03 + H01*C13 = zero) + H02*C23 + H03*C33.
        BS=0}; // Neglect correlation between r-phi & r-z planes for now.
  typedef AP_FIXED(B25,BS00)  TS00;
  typedef AP_FIXED(B25,BS01)  TS01;
  typedef AP_FIXED(B25,BS12)  TS12;
  typedef AP_FIXED(B25,BS13)  TS13;
  typedef AP_FIXED(BCORR,BS)  TS; 

public:

  MatrixS(const MatrixH<5>& H, const MatrixC<5>& C);

public:
  
  TS00 _00;
  TS01 _01;
  TS12 _12;
  TS13 _13;
  TS            _02, _03,
      _10, _11          ;
};

// Covariance matrix of helix params.

template <>
class MatrixC<5> {
public:
  typedef AP_UFIXED(1,1)      T0; // HLS doesn't like zero bit variables.

  // Determine input helix coviaraiance matrix.
  MatrixC(const KFstateHLS<5>& stateIn) :
             _00(stateIn.cov_00), _11(stateIn.cov_11), _22(stateIn.cov_22), _33(stateIn.cov_33), _01(stateIn.cov_01), _23(stateIn.cov_23), 
             _02(0), _03(0), _12(0), _13(0),
             _10(_01), _32(_23), _20(_02), _30(_03), _21(_12), _31(_13) {}

  // Calculate output helix covariance matrix: C' = C - K*H*C = C - K*S.
  MatrixC(const MatrixC<5>& C, const MatrixK<5>& K, const MatrixS<5>& S);

public:
  // Elements that are finite
// Maxeller wierdly uses signed 25 bits for these, so use unsigned 24 instead to match the DSP abilities.
//  KFstateHLS<5>::TC00 _00;
//  KFstateHLS<5>::TC11 _11;
//  KFstateHLS<5>::TC22 _22;
//  KFstateHLS<5>::TC33 _33;
  AP_UFIXED(B24,BHC00-1) _00; // One less integer bit as no sign required.
  AP_UFIXED(B24,BHC11-1) _11;
  AP_UFIXED(B24,BHC22-1) _22;
  AP_UFIXED(B24,BHC33-1) _33;
  KFstateHLS<5>::TC01 _01; // (inv2R, phi0) -- other off-diagonal elements assumed negligeable.
  KFstateHLS<5>::TC23 _23; // (tanL,  z0)   -- other off-diagonal elements assumed negligeable.
  // Elements that are zero.
  const T0 _02, _03, _12, _13;
  // Elements below the diagonal of this symmetric matrix.
  const KFstateHLS<5>::TC01 &_10;
  const KFstateHLS<5>::TC23 &_32;
  const T0 &_20, &_30, &_21, &_31;
};

// S(transpose) = C*H(transpose)

template <>
class MatrixS_transpose<5> {
public:
  typedef MatrixS<5>::TS00  TS00;
  typedef MatrixS<5>::TS01  TS01;
  typedef MatrixS<5>::TS12  TS12;
  typedef MatrixS<5>::TS13  TS13;
  typedef MatrixS<5>::TS    TS;
  MatrixS_transpose(const MatrixS<5>& S) : _00(S._00), _10(S._01), _21(S._12), _31(S._13),
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

template <>
class MatrixR<5> {
public:
  enum {BH=MatrixH<5>::BH,
	BS00=MatrixS<5>::BS00,
	BS01=MatrixS<5>::BS01,
	BS12=MatrixS<5>::BS12,
	BS13=MatrixS<5>::BS13,
	BS  =MatrixS<5>::BS,
        // Calculate number of integer bits required for elements of R.
        BR00 = MAX2(MatrixV::BVPP, MAX2(BH+BS00, BS01)) - BODGE<5>::R, // H00*St00 + H01*St10 + (H02*St20 + H03*St30 = zero)
	BR11 = MAX2(MatrixV::BVZZ, MAX2(BH+BS12, BS13)) - BODGE<5>::R, // (H10*St01 + H11*St11 = zero) + H12*St21 + H13*St31
	BR01 = MAX2(BH+BS, BS)                                     // H00*St01 + H01*St11 + (H02*St21 + H03*St31 = zero)
       };  
  typedef SW_UFIXED(B34,BR00) TR00;
  typedef SW_UFIXED(B34,BR11) TR11;
  typedef SW_UFIXED(B34,BR01) TR01;

public:
  MatrixR(const MatrixV& V, const MatrixH<5>& H, const MatrixS_transpose<5>& St);

public:
  TR00  _00;
  TR11  _11;
  TR01  _01;
  TR01& _10; // Matrix symmetric, so can use reference.
};

// Kalman gain matrix K = S(transpose)*R(inverse).

template <>
class MatrixK<5> {
public:
  enum {BS00=MatrixS<5>::BS00,
	BS01=MatrixS<5>::BS01,
	BS12=MatrixS<5>::BS12,
	BS13=MatrixS<5>::BS13,
        BIR00=MatrixInverseR<5>::BIR00,
        BIR11=MatrixInverseR<5>::BIR11,
        BIR01=MatrixInverseR<5>::BIR01,
        BK00=(BS00+BIR00) - BODGE<5>::K,  // St00*Rinv00 (+ St01*Rinv10 = zero)
        BK10=(BS01+BIR00) - BODGE<5>::K,  // St10*Rinv00 (+ St11*Rinv10 = zero)
        BK21=(BS12+BIR11) - BODGE<5>::K,  // (St20*Rinv01 = zero) + St21*Rinv11
        BK31=(BS13+BIR11) - BODGE<5>::K}; // (St30*Rinv01 = zero) + St31*Rinv11
  typedef SW_FIXED(B35,BK00)  TK00;
  typedef SW_FIXED(B35,BK10)  TK10;
  typedef SW_FIXED(B35,BK21)  TK21;
  typedef SW_FIXED(B35,BK31)  TK31;
  typedef SW_FIXED(BCORR,0)   TK; // Neglect correlation between r-phi & r-z
  MatrixK(const MatrixS_transpose<5>& St, const MatrixInverseR<5>& RmatInv);
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

template <>
class VectorRes<5> {
public:
  // Use higher granularity for residuals than for stubs.
  // BODGE<5>::RES should be slightly larger than BODGE_V as hits can be several sigma wrong.
  // Add one extra fractional bit relative to stub, to avoid additional precision loss.
  typedef AP_FIXED(B18-BODGE<5>::RES+1,BSP-BODGE<5>::RES) TRP;
  typedef AP_FIXED(B18-BODGE<5>::RES+1,BSZ-BODGE<5>::RES) TRZ;

public:
  VectorRes(const VectorM& m, const MatrixH<5>& H, const VectorX<5>& x);

public:
  TRP _0;
  TRZ _1;
};

// Vector of helix params.

template <>
class VectorX<5> {
public:
  // Determine input helix params.
  VectorX(const KFstateHLS<5>::TR& inv2R, const KFstateHLS<5>::TP& phi0, const KFstateHLS<5>::TT& tanL, const KFstateHLS<5>::TZ& z0) : _0(inv2R), _1(phi0), _2(tanL), _3(z0) {} 

  // Calculate output helix params: x' = x + K*res
  VectorX(const VectorX<5>& x, const MatrixK<5>& K, const VectorRes<5>& res);

public:
  KFstateHLS<5>::TR _0;
  KFstateHLS<5>::TP _1;  
  KFstateHLS<5>::TT _2;  
  KFstateHLS<5>::TZ _3;  
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
