/**
 * This defines the KF matrices and the operations performance on them.
 *
 *  All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 *
 * Author: Ian Tomalin
 */

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#else
#include "KalmanMatricesHLS.h"
#endif

#ifdef PRINT_SUMMARY
#include <iostream>
#endif

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Covariance matrix of stub coords.

MatrixV::MatrixV(const StubHLS::TR& r, const StubHLS::TZ& z, const KFstateHLS::TR& inv2R, const KFstateHLS::TT& tanL, const AP_INT(BHT_M)& mBin) : _01(0), _10(_01) {

  // Numbers from http://ghugo.web.cern.ch/ghugo/layouts/July/OT613_200_IT4025/layout.html
  // Module pitch for tilted tracker geometry divided by sqrt(12) to get resolution.
  //static const float invRoot12  = 1./sqrt(12.);
  static const float invRoot12  = 0.288675; // 1/sqrt(12)
  // Declaring these static causes cosimulation to set them to zero. Why?
  // But it is OK if invRoot12 is defined as 0.288 instead of 1/sqrt(12.)
  static const TM pitchPS    = rphiMult*invRoot12*0.0099;  // < 1024
  static const TM pitch2S    = rphiMult*invRoot12*0.0089;
  // Factor 1.41 allows for digitisation granularity (Sioni).
  //static const TM lengthPS   = rMult*invRoot12*0.15;
  static const TM lengthPS   = 1.41*rMult*invRoot12*0.15; 
  static const TM length2S   = rMult*invRoot12*5.02;       // < 16

  // These boundaries are for tilted tracker geometry.
  static const StubHLS::TZ zBarrel     = rMult*125.; // Largest z in barrel.
  static const StubHLS::TZ zWheel12    = rMult*170.; // Largest z of endcap wheels 1 or 2.
  static const StubHLS::TR rPSbarrel   = rMult*60.0;  // r of PS-2S transition in barrel.
  static const StubHLS::TR rPSwheel12  = rMult*66.4; // r below which stub certain to be PS if in endcap wheels 1 or 2.
  static const StubHLS::TR rPSwheel345 = rMult*64.6; // r below which stub certain to be PS if in endcap wheels 3, 4 or 5.

  // Initialise pitch/r ROMs.
  static const PitchOverR_2 calcPitchPSoverR_2(pitchPS);
  static const PitchOverR_2 calcPitch2SoverR_2(pitch2S);

  PitchOverR_2::TPOR pitchPSoverR_2 = calcPitchPSoverR_2.get[r.to_uint() >> PitchOverR_2::BRED];
  PitchOverR_2::TPOR pitch2SoverR_2 = calcPitch2SoverR_2.get[r.to_uint() >> PitchOverR_2::BRED];

#ifdef PRINT_SUMMARY
  CHECK_AP::checkCalc("p*p/r*r", pitch2SoverR_2, double(pitch2S*pitch2S)/double(r*r));
#endif
#ifdef PRINT
  std::cout<<"p/r check "<<pitchPSoverR_2<<" vs "<<double(pitchPS*pitchPS)/double(r*r)<<std::endl;
#endif

  StubHLS::TZ absZ = z;   // HLS has no abs() function.
  if (z < 0) {absZ = -z;};

  // Use same granularity for resolution as for residuals.
  // (Define as signed, so dont have to worry if tanL or inv2R are -ve).
  MatrixV::TVPP sigmaPhi2;  // Uncertainty squared in phi.
  AP_FIXED(B18,BVZ) sigmaZ; // Uncertainty in z. 
  MatrixV::TVPP sigmaPhiExtra2;

  // Initialize ROM used to calculate contribution to phi uncertainty from (r,phi) to (z,phi) conversion in endcap.
  static const InvPt2 calcPhiExtra2_PS(1, lengthPS);
  static const InvPt2 calcPhiExtra2_2S(1, length2S);

  if (absZ < zBarrel) {
    // Barrel
#ifdef PRINT
    std::cout<<"BARREL "<<absZ<<" "<<zBarrel<<std::endl;
#endif
    if (r < rPSbarrel) {
      _2Smodule = false;
      sigmaPhi2 = pitchPSoverR_2;	
      sigmaZ    = lengthPS;	
    } else {
      _2Smodule = true;
      sigmaPhi2 = pitch2SoverR_2;	
      sigmaZ    = length2S;	
    }
    sigmaPhiExtra2 = 0.;

  } else {

    // Endcap
    bool psEndcap = (absZ < zWheel12)  ?  (r < rPSwheel12)  :  (r < rPSwheel345);
    if (psEndcap) {
#ifdef PRINT
      std::cout<<"ENDCAP1 "<<tanL<<std::endl;
#endif
      _2Smodule = false;
      sigmaPhi2 = pitchPSoverR_2;	
      sigmaZ    = lengthPS*tanL;
      sigmaPhiExtra2 = calcPhiExtra2_PS.getIt(mBin);
    } else { 
#ifdef PRINT
      std::cout<<"ENDCAP2 "<<tanL<<std::endl;
#endif
      _2Smodule = true;
      sigmaPhi2 = pitch2SoverR_2;	
      sigmaZ    = length2S*tanL;
      sigmaPhiExtra2 = calcPhiExtra2_2S.getIt(mBin);
    }
  }

  // Allow for scattering in r-phi plane (since hit resolution is best there)
  static const InvPt2 calcScatTerm2;
  TVPP sigmaPhiScat2 = calcScatTerm2.getIt(mBin);

  // IRT - check if using DSPs gives better accuracy that LUT
  //static const float junk = 2*rMult*kalmanMultScatTerm/invPtToInvR;
  //static const float junk2 = junk*junk;
  //TVPP sigmaPhiScat2 = junk2*float(inv2R*inv2R); 

  _00 = sigmaPhi2 + sigmaPhiExtra2 + sigmaPhiScat2;
  _11 = sigmaZ*sigmaZ;

#ifdef PRINT_SUMMARY
  CHECK_AP::checkCalc("sigmaPhiScat2", sigmaPhiScat2,
    pow(kalmanMultScatTerm*2.*double(rMult)*double(inv2R)/invPtToInvR, 2), 0.2, pow(0.0002*phiMult,2));
#endif
#ifdef PRINT
  std::cout<<"2S="<<_2Smodule<<" ENDCAP="<<(absZ > zBarrel)<<" (r,z)=("<<r<<", "<<z<<")"<<std::endl;
  std::cout<<"SIGMA RPHI="<<sqrt(double(_00))/double(phiMult)<<" SIGMA_RZ="<<sqrt(double(_11))/double(rMult)<<" EXTRA="<<double(sigmaPhiExtra)/double(phiMult)<<" SCAT="<<sqrt(double(sigmaPhiScat2))/double(phiMult)<<std::endl;
  std::cout<<"SCAT CHECK: "<<mBin<<" "<<double(inv2R)/double(inv2R_digi_cut)<<" RESULT: DIGI="<<double(sigmaPhiScat2)<<" FLOAT="<<pow(kalmanMultScatTerm*2.*double(rMult)*double(inv2R)/invPtToInvR, 2)<<std::endl;
  std::cout<<"  V00="<<_00<<"   V11="<<_11<<std::endl;
#endif

  //static const float rats = sqrt(2.9713); // FAIL
  //static const AP_UFIXED(40,20) length2ST   = rats;
  //_11 = length2ST;

  //static const float rats = 3./2.; // GOOD
  //static const AP_UFIXED(40,20) length2ST   = rats;
  //_11 = length2ST;

  //static const float rats = sqrt(2.9713); // GOOD
  //const AP_UFIXED(40,20) length2ST   = rats;
  //_11 = length2ST;

  //static const float rats = sqrt(2.9713); // GOOD
  //_11 = rats;
}

// Calculate S = H * C

MatrixS::MatrixS(const MatrixH& H, const MatrixC& C) {
  _00 = H._00 * C._00 + H._01 * C._10 + H._02 * C._20 + H._03 * C._30;  
  _01 = H._00 * C._01 + H._01 * C._11 + H._02 * C._21 + H._03 * C._31;  
  _02 = H._00 * C._02 + H._01 * C._12 + H._02 * C._22 + H._03 * C._32;  
  _03 = H._00 * C._03 + H._01 * C._13 + H._02 * C._23 + H._03 * C._33;  
  _10 = H._10 * C._00 + H._11 * C._10 + H._12 * C._20 + H._13 * C._30;  
  _11 = H._10 * C._01 + H._11 * C._11 + H._12 * C._21 + H._13 * C._31;  
  _12 = H._10 * C._02 + H._11 * C._12 + H._12 * C._22 + H._13 * C._32;  
  _13 = H._10 * C._03 + H._11 * C._13 + H._12 * C._23 + H._13 * C._33;  

#ifdef PRINT_SUMMARY
  double s00 = H._00 * C._00 + H._01 * C._10 + H._02 * C._20 + H._03 * C._30;
  double s01 = H._00 * C._01 + H._01 * C._11 + H._02 * C._21 + H._03 * C._31;
  double s02 = H._00 * C._02 + H._01 * C._12 + H._02 * C._22 + H._03 * C._32;
  double s03 = H._00 * C._03 + H._01 * C._13 + H._02 * C._23 + H._03 * C._33;
  double s10 = H._10 * C._00 + H._11 * C._10 + H._12 * C._20 + H._13 * C._30;
  double s11 = H._10 * C._01 + H._11 * C._11 + H._12 * C._21 + H._13 * C._31;
  double s12 = H._10 * C._02 + H._11 * C._12 + H._12 * C._22 + H._13 * C._32;
  double s13 = H._10 * C._03 + H._11 * C._13 + H._12 * C._23 + H._13 * C._33;
  CHECK_AP::checkCalc("S00", _00, s00, 0.03);
  CHECK_AP::checkCalc("S01", _01, s01, 0.03);
  CHECK_AP::checkCalc("S02", _02, s02, 0.03);
  CHECK_AP::checkCalc("S03", _03, s03, 0.03);
  CHECK_AP::checkCalc("S10", _10, s10, 0.03);
  CHECK_AP::checkCalc("S11", _11, s11, 0.03);
  CHECK_AP::checkCalc("S12", _12, s12, 0.03);
  CHECK_AP::checkCalc("S13", _13, s13, 0.03);
#endif
}

// Calculate covariance matrix of predicted residuals R = V + H*C*Ht = V + H*St.

MatrixR::MatrixR(const MatrixV& V, const MatrixH& H, const MatrixS_transpose& St) : 
  _10(_01) 
{
  _00 = V._00 + (H._00*St._00 + H._01*St._10 + H._02*St._20 + H._03*St._30);
  _01 = V._01 + (H._00*St._01 + H._01*St._11 + H._02*St._21 + H._03*St._31);
  // R._10 // Matrix symmetric so don't need to calculate this element.
  _11 = V._11 + (H._10*St._01 + H._11*St._11 + H._12*St._21 + H._13*St._31);

#ifdef PRINT_SUMMARY
  double r00 = V._00 + (H._00*St._00 + H._01*St._10 + H._02*St._20 + H._03*St._30);
  double r01 = V._01 + (H._00*St._01 + H._01*St._11 + H._02*St._21 + H._03*St._31);
  double r11 = V._11 + (H._10*St._01 + H._11*St._11 + H._12*St._21 + H._13*St._31);
  CHECK_AP::checkCalc("R00", _00, r00);
  CHECK_AP::checkCalc("R01", _01, r01);
  CHECK_AP::checkCalc("R11", _11, r11);
#endif
}

// Inverse of matrix R. 

MatrixInverseR::MatrixInverseR(const MatrixR& R) : _10(_01) 
{
  // The determinant is positive, as the covariance matrix is almost diagonal.
  enum {BDW = B48, BDI = (MatrixR::BR00 + MatrixR::BR11) - BODGE_DET};
  enum {B7 = 7}; // Number of bits needed to describe a number in range 0 to B48 safely.
  const AP_UFIXED(BDW,BDI) det = (R._00 * R._11 - R._01 * R._10);

  //--- Asking HLS to calculate 1/det is slow & expensive. So calculate it with home-made algorithm,
  //--- involving finding the leading non-zero bit, and then using a look up table to take the reciprocal.

  // Find leading non-zero bit.
  enum {iMIN = 2*BDET}; // Increasing this reduces FPGA resources. But don't make so big that most significant bit of det is below this value. If reduced below 2*BDET, code below must be modified to allow for det_short being less than 2*BDET bits.

  AP_UINT(B7) msb = iMIN; // most-significant bit counting from zero at the right-most bit.

  // This takes 5 clock cycles & uses 4 BRAM.
  for (AP_UINT(B7) i = iMIN+1; i < BDW; i++) {
    if (det[i]) msb = i;
  }

  // // This uses built-in C function to save 1 clock cycle, but at cost of 2 extra BRAM.
  // AP_UINT(B7) lzero = __builtin_clz((unsigned int)(det.range(BDW-1,iMIN+1))); // Finds location of leading non-zero bit.
  // AP_UINT(B7) msb = (32+iMIN)-lzero;

  AP_UINT(B7) lsb = msb - BDET + 1;
  const AP_UINT(BDET) det_veryshort = det.range(msb, lsb);
  SW_UFIXED(2*BDET,BDET) det_short;
  det_short.range(2*BDET-1, 0) = det.range(msb, lsb - BDET);

  // // This saves 2 clock cycles, at cost of 1 extra DSP. But it fails timing when compiled in Vivado.
  // AP_FIXED(BDW,BDI) det_noLeadZeros = det;
  // for (AP_UINT(B7) i = BDW - 1; i > iMIN; i--) {
  //   if (det_noLeadZeros[BDW-1]) {
  //     msb = i;
  //     break;
  //   } else {
  //     det_noLeadZeros = det_noLeadZeros << 1;
  //   }
  // }
  // AP_UINT(B7) lsb = msb - BDET + 1;
  // const AP_UINT(BDET) det_veryshort = det_noLeadZeros.range(BDW-1, BDW-BDET);
  // SW_UFIXED(2*BDET,BDET) det_short;
  // det_short.range(2*BDET-1, 0) = det_noLeadZeros.range(BDW-1, BDW-2*BDET);

  // Take reciprocal with lookup table.
  static const OneOverInt calcOneOverInt;
  SW_UFIXED(BDET,OneOverInt::BOI) invDet_veryshort = calcOneOverInt.getIt(det_veryshort);

  // Determine higher order corrections to reciprocal.
  // (If det = a + b, where b small, this solves for x, where x is small, (a + b) * (A + x) = 1,
  // where A is a LUT approximation to 1/a. So x = A*(1 - A*a - A*b), and (A+x) = A*(2 - A*a - A*b). 

  // This is inverse determinant, aside from shift factor SHIFT.
  // First line uses one less clock cycle than second line, but one more DSP ...
  SW_UFIXED(B18, OneOverInt::BOI) invDet_short = invDet_veryshort * (SW_UFIXED(1,2)(2.) - invDet_veryshort * det_short); 
  //SW_UFIXED(B18, OneOverInt::BOI) invDet_short = 2 * invDet_veryshort - (invDet_veryshort * invDet_veryshort) * det_short;
 
  AP_INT(B7) SHIFT = lsb - (BDW - BDI); // Bit shift to be applied to inverse determinant.

  // Calculate max & min values that SHIFT can take over all events.
  enum {MAX_LSB = (BDW - 1) - BDET + 1, MAX_SHIFT = MAX_LSB - (BDW - BDI),
	MIN_LSB = iMIN      - BDET + 1, MIN_SHIFT = MIN_LSB - (BDW - BDI)};

  // Invert matrix.
  _00 =  SW_UFIXED(B34-MIN_SHIFT+MAX_SHIFT, BIR11+MAX_SHIFT) (invDet_short*R._11) >> SHIFT;
  _11 =  SW_UFIXED(B34-MIN_SHIFT+MAX_SHIFT, BIR00+MAX_SHIFT) (invDet_short*R._00) >> SHIFT;
  _01 =  SW_FIXED(BCORR-MIN_SHIFT+MAX_SHIFT, BIR01-B34+BCORR+MAX_SHIFT) (-(invDet_short*R._10)) >> SHIFT;

#ifdef PRINT
  std::cout<<"MatrixInverseR: Det="<<det<<" det_veryshort="<<det_veryshort<<" invDet_veryshort="<<invDet_veryshort<<" det_range2="<<det_range2<<" invDet_short="<<invDet_short<<" det*invDet_short="<<double(det)*double(invDet_short)/double(1 << SHIFT)<<std::endl;
#endif

#ifdef PRINT_SUMMARY
  // Check assumed bit ranges are OK.
  CHECK_AP::checkIntRange("MSB", BDW-1, iMIN, msb);
  CHECK_AP::checkIntRange("SHIFT", MAX_SHIFT, MIN_SHIFT, SHIFT);
  CHECK_AP::checkIntRange("Det[MSB]", 1, 1, det_veryshort[AP_UINT(B7)(BDET-1)]);
  double trueDet = double(R._00)*double(R._11)-double(R._01)*double(R._10);
  double trueInvDet = 1./trueDet;
  double true_ri00 =  double(R._11)*trueInvDet;
  double true_ri11 =  double(R._00)*trueInvDet;
  double true_ri01 = -double(R._10)*trueInvDet;
  double invDet = double(invDet_short)*double(AP_UFIXED(1-MIN_SHIFT+MAX_SHIFT, 1-MIN_SHIFT)(1) >>  SHIFT);
  CHECK_AP::checkCalc("DET", det, trueDet, 0.00001);
  // Precision of this (controlled by BDET) is critical.
  CHECK_AP::checkCalc("INVDET", invDet_short,
                      trueInvDet/double(AP_UFIXED(1-MIN_SHIFT+MAX_SHIFT, 1-MIN_SHIFT)(1) >>  SHIFT), 0.0001);
  CHECK_AP::checkCalc("INVR00", _00, true_ri00, 0.001);
  CHECK_AP::checkCalc("INVR01", _01, true_ri01, 0.001);
  CHECK_AP::checkCalc("INVR11", _11, true_ri11, 0.001);
#endif
}

// Kalman gain matrix K = S*R(inverse). (Actually this *Det(R) for precision reasons).

MatrixK::MatrixK(const MatrixS_transpose& St, const MatrixInverseR& RmatInv) {
#ifdef USE_FIXED
  _00 =  St._00 * RmatInv._00 + St._01 * RmatInv._10;
  _10 =  St._10 * RmatInv._00 + St._11 * RmatInv._10;
  _20 =  St._20 * RmatInv._00 + St._21 * RmatInv._10;
  _30 =  St._30 * RmatInv._00 + St._31 * RmatInv._10;
  _01 =  St._00 * RmatInv._01 + St._01 * RmatInv._11;
  _11 =  St._10 * RmatInv._01 + St._11 * RmatInv._11;
  _21 =  St._20 * RmatInv._01 + St._21 * RmatInv._11;
  _31 =  St._30 * RmatInv._01 + St._31 * RmatInv._11;
#else
  _00 =  SW_FLOAT(St._00) * RmatInv._00 + SW_FLOAT(St._01) * RmatInv._10;
  _10 =  SW_FLOAT(St._10) * RmatInv._00 + SW_FLOAT(St._11) * RmatInv._10;
  _20 =  SW_FLOAT(St._20) * RmatInv._00 + SW_FLOAT(St._21) * RmatInv._10;
  _30 =  SW_FLOAT(St._30) * RmatInv._00 + SW_FLOAT(St._31) * RmatInv._10;
  _01 =  SW_FLOAT(St._00) * RmatInv._01 + SW_FLOAT(St._01) * RmatInv._11;
  _11 =  SW_FLOAT(St._10) * RmatInv._01 + SW_FLOAT(St._11) * RmatInv._11;
  _21 =  SW_FLOAT(St._20) * RmatInv._01 + SW_FLOAT(St._21) * RmatInv._11;
  _31 =  SW_FLOAT(St._30) * RmatInv._01 + SW_FLOAT(St._31) * RmatInv._11;
#endif

#ifdef PRINT_SUMMARY
  double k00 =  double(St._00) * double(RmatInv._00) + double(St._01) * double(RmatInv._10);
  double k10 =  double(St._10) * double(RmatInv._00) + double(St._11) * double(RmatInv._10);
  double k20 =  double(St._20) * double(RmatInv._00) + double(St._21) * double(RmatInv._10);
  double k30 =  double(St._30) * double(RmatInv._00) + double(St._31) * double(RmatInv._10);
  double k01 =  double(St._00) * double(RmatInv._01) + double(St._01) * double(RmatInv._11);
  double k11 =  double(St._10) * double(RmatInv._01) + double(St._11) * double(RmatInv._11);
  double k21 =  double(St._20) * double(RmatInv._01) + double(St._21) * double(RmatInv._11);
  double k31 =  double(St._30) * double(RmatInv._01) + double(St._31) * double(RmatInv._11);
  CHECK_AP::checkCalc("K00", _00, k00, 0.001);
  CHECK_AP::checkCalc("K10", _10, k10, 0.001);
  CHECK_AP::checkCalc("K20", _20, k20, 0.001);
  CHECK_AP::checkCalc("K30", _30, k30, 0.001);
  CHECK_AP::checkCalc("K01", _01, k01, 0.001);
  CHECK_AP::checkCalc("K11", _11, k11, 0.001);
  CHECK_AP::checkCalc("K21", _21, k21, 0.001);
  CHECK_AP::checkCalc("K31", _31, k31, 0.001);
#endif
}

// Hit residuals: res = m - H*x. 

VectorRes::VectorRes(const VectorM& m, const MatrixH& H, const VectorX& x) {
  _0 = m._0 - (H._00 * x._0 + H._01 * x._1 + H._02 * x._2 + H._03 * x._3);  
  _1 = m._1 - (H._10 * x._1 + H._11 * x._1 + H._12 * x._2 + H._13 * x._3);  
#ifdef PRINT_SUMMARY
  double r0 =  double(m._0) - (double(H._00) * double(x._0) + double(H._01) * double(x._1) + 
                               double(H._02) * double(x._2) + double(H._03) * double(x._3));
  double r1 =  double(m._1) - (double(H._10) * double(x._0) + double(H._11) * double(x._1) + 
                               double(H._12) * double(x._2) + double(H._13) * double(x._3));
  CHECK_AP::checkCalc("RES0", _0, r0, 0.1, 0.1);
  CHECK_AP::checkCalc("RES1", _1, r1, 0.1, 0.1);
#endif
}

// Calculate output helix params: x' = x + K*res

VectorX::VectorX(const VectorX& x, const MatrixK& K, const VectorRes& res) {
#ifdef USE_FIXED
  typedef MatrixK::TK00_short TK00_short;
  typedef MatrixK::TK10_short TK10_short;
  typedef MatrixK::TK21_short TK21_short;
  typedef MatrixK::TK31_short TK31_short;
  typedef MatrixK::TK         TK;
  _0 = x._0 + KFstateHLS::TR(TK00_short(K._00) * res._0 + TK        (K._01) * res._1);
  _1 = x._1 + KFstateHLS::TP(TK10_short(K._10) * res._0 + TK        (K._11) * res._1);
  _2 = x._2 + KFstateHLS::TT(TK        (K._20) * res._0 + TK21_short(K._21) * res._1);
  _3 = x._3 + KFstateHLS::TZ(TK        (K._30) * res._0 + TK31_short(K._31) * res._1);
#else
  _0 = x._0 + KFstateHLS::TR(K._00 * SW_FLOAT(res._0) + K._01 * SW_FLOAT(res._1)); 
  _1 = x._1 + KFstateHLS::TP(K._10 * SW_FLOAT(res._0) + K._11 * SW_FLOAT(res._1)); 
  _2 = x._2 + KFstateHLS::TT(K._20 * SW_FLOAT(res._0) + K._21 * SW_FLOAT(res._1)); 
  _3 = x._3 + KFstateHLS::TZ(K._30 * SW_FLOAT(res._0) + K._31 * SW_FLOAT(res._1)); 
#endif
}

// Calculate output helix covariance matrix: C' = C - K*H*C = C - K*S.

MatrixC::MatrixC(const MatrixC& C, const MatrixK& K, const MatrixS& S) :
  _02(0), _03(0), _12(0), _13(0),
  _10(_01), _32(_23), _20(_02), _30(_03), _21(_12), _31(_13)
{
  // Covariance matrix is symmetric & some elements can be neglected.
#ifdef USE_FIXED
  _00 =  C._00 - KFstateHLS::TC00(K._00 * S._00 + K._01 * S._10);
  _11 =  C._11 - KFstateHLS::TC11(K._10 * S._01 + K._11 * S._11);
  _22 =  C._22 - KFstateHLS::TC22(K._20 * S._02 + K._21 * S._12);
  _33 =  C._33 - KFstateHLS::TC33(K._30 * S._03 + K._31 * S._13);
  _01 =  C._01 - KFstateHLS::TC01(K._00 * S._01 + K._01 * S._11);
  _23 =  C._23 - KFstateHLS::TC23(K._20 * S._03 + K._21 * S._13);
#else
  _00 =  C._00 - KFstateHLS::TC00(K._00 * SW_FLOAT(S._00) + K._01 * SW_FLOAT(S._10));
  _11 =  C._11 - KFstateHLS::TC11(K._10 * SW_FLOAT(S._01) + K._11 * SW_FLOAT(S._11));
  _22 =  C._22 - KFstateHLS::TC22(K._20 * SW_FLOAT(S._02) + K._21 * SW_FLOAT(S._12));
  _33 =  C._33 - KFstateHLS::TC33(K._30 * SW_FLOAT(S._03) + K._31 * SW_FLOAT(S._13));
  _01 =  C._01 - KFstateHLS::TC01(K._00 * SW_FLOAT(S._01) + K._01 * SW_FLOAT(S._11));
  _23 =  C._23 - KFstateHLS::TC23(K._20 * SW_FLOAT(S._03) + K._21 * SW_FLOAT(S._13));
#endif

#ifdef PRINT_SUMMARY
  double c00new = double(C._00) - (double(K._00) * double(S._00) + double(K._01) * double(S._10));
  double c11new = double(C._11) - (double(K._10) * double(S._01) + double(K._11) * double(S._11));
  double c22new = double(C._22) - (double(K._20) * double(S._02) + double(K._21) * double(S._12));
  double c33new = double(C._33) - (double(K._30) * double(S._03) + double(K._31) * double(S._13));
  double c01new = double(C._01) - (double(K._00) * double(S._01) + double(K._01) * double(S._11));
  double c23new = double(C._23) - (double(K._20) * double(S._03) + double(K._21) * double(S._13));
  CHECK_AP::checkCalc("C00_new", _00, c00new, 0.01);
  CHECK_AP::checkCalc("C11_new", _11, c11new, 0.01);
  CHECK_AP::checkCalc("C22_new", _22, c22new, 0.01);
  CHECK_AP::checkCalc("C33_new", _33, c33new, 0.01);
  CHECK_AP::checkCalc("C01_new", _01, c01new, 0.01);
  CHECK_AP::checkCalc("C23_new", _23, c23new, 0.01);
  CHECK_AP::checkDet("C_new(rphi)",_00,_11,_01);
  CHECK_AP::checkDet("C_new(rz)"  ,_22,_33,_23);
#endif
}

#ifdef CMSSW_GIT_HASH
}

}
#endif

