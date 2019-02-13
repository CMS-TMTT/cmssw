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
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS4.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS5.h"
#else
#include "KalmanMatricesHLS.h"
#include "KalmanMatricesHLS4.h"
#include "KalmanMatricesHLS5.h"
#endif

#ifdef PRINT_SUMMARY
#include <iostream>
#endif

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

//--- Explicit instantiation required for all non-specialized templates, to allow them to be implemented 
//--- in .cc files.

template MatrixInverseR<4>::MatrixInverseR(const MatrixR<4>& R);

template MatrixInverseR<5>::MatrixInverseR(const MatrixR<5>& R);

// Covariance matrix of stub coords.

MatrixV::MatrixV(const StubHLS::TR& r, const StubHLS::TZ& z, const KFstateHLS<4>::TR& inv2R, const KFstateHLS<4>::TT& tanL, const AP_INT(BHT_M)& mBin) : _01(0), _10(_01) {

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

#ifdef CMSSW_GIT_HASH
  StubHLS::TZ absZ = fabs(float(z));  
#else
  StubHLS::TZ absZ = hls::abs(z);  
#endif

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
  CHECK_AP::checkCalc("V00", _00,
		      double(sigmaPhi2) + double(sigmaPhiExtra2) + pow(kalmanMultScatTerm*2.*double(rMult)*double(inv2R)/invPtToInvR, 2), 0.4, 30.); // Very inaccurate as mBin not identical to fitted inv2R?
  CHECK_AP::checkCalc("V11", _11, pow(double(sigmaZ), 2));
#endif
#ifdef PRINT
  std::cout<<"2S="<<_2Smodule<<" ENDCAP="<<(absZ > zBarrel)<<" (r,z)=("<<r<<", "<<z<<")"<<std::endl;
  std::cout<<"SIGMA RPHI="<<sqrt(double(_00))/double(phiMult)<<" SIGMA_RZ="<<sqrt(double(_11))/double(rMult)<<" EXTRA="<<sqrt(double(sigmaPhiExtra2))/double(phiMult)<<" SCAT="<<sqrt(double(sigmaPhiScat2))/double(phiMult)<<std::endl;
  std::cout<<"SCAT CHECK: "<<mBin<<" "<<double(inv2R)/double(inv2Rcut)<<" RESULT: DIGI="<<double(sigmaPhiScat2)<<" FLOAT="<<pow(kalmanMultScatTerm*2.*double(rMult)*double(inv2R)/invPtToInvR, 2)<<std::endl;
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

// Inverse of matrix R. 

template <unsigned int NPAR>
MatrixInverseR<NPAR>::MatrixInverseR(const MatrixR<NPAR>& R) : _10(_01) 
{
  // The determinant is positive, as the covariance matrix is almost diagonal.
  enum {BDW = B48, BDI = (MatrixR<NPAR>::BR00 + MatrixR<NPAR>::BR11) - BODGE<NPAR>::DET};
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
  // If ap_fixed<N,I> is shifted right by SHIFT, which can be signed, no truncation occurs if it is 
  // first cast to AP<N+|S|, I+Max(0,-S)>
  // If SHIFT can take any value between MIN_SHIFT & MAX_SHIFT, it must instead be cast to <N+nr,I+ir>
  // where ir = Max(0, -MIN_SHIFT), so nr = Max(MAX_SHIFT + ir, -MIN_SHIFT ) = MAX_SHIFT + Max(0, -MIN_SHIFT).
  // If ap_fixed<N,I> is shifted left by SHIFT, one replaces SHIFT --> -SHIFT and 
  // MAX_SHIFT <--> -MIN_SHIFT in previous formula, so no 
  // truncation occurs if it is first cast to <N+nl,I+il>
  // where il = Max(0, MAX_SHIFT) and nr = -MIN_SHIFT + Max(0, MAX_SHIFT).
  enum {NR_EXTRA =  MAX_SHIFT + MAX2(0, -MIN_SHIFT), IR_EXTRA = MAX2(0, -MIN_SHIFT),
        NL_EXTRA = -MIN_SHIFT + MAX2(0,  MAX_SHIFT), IL_EXTRA = MAX2(0,  MAX_SHIFT)};

  // Invert matrix.
  //  _00 =  SW_UFIXED(B34-MIN_SHIFT+MAX_SHIFT, BIR11+MAX_SHIFT) (invDet_short*R._11) >> SHIFT;
  //  _11 =  SW_UFIXED(B34-MIN_SHIFT+MAX_SHIFT, BIR00+MAX_SHIFT) (invDet_short*R._00) >> SHIFT;
  //  _01 =  SW_FIXED(BCORR-MIN_SHIFT+MAX_SHIFT, BIR01-B34+BCORR+MAX_SHIFT) (-(invDet_short*R._10)) >> SHIFT;
  _00 =  SW_UFIXED(B34+NL_EXTRA, BIR11+IL_EXTRA) (invDet_short*R._11) >> SHIFT;
  _11 =  SW_UFIXED(B34+NL_EXTRA, BIR00+IL_EXTRA) (invDet_short*R._00) >> SHIFT;
  _01 =  SW_FIXED(BCORR+NL_EXTRA, BIR01+IL_EXTRA) (-(invDet_short*R._10)) >> SHIFT;

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
  //  double invDet = double(invDet_short)*double(AP_UFIXED(1-MIN_SHIFT+MAX_SHIFT, 1-MIN_SHIFT)(1) >>  SHIFT);
  double invDet = double(invDet_short)*double(AP_UFIXED(1+NR_EXTRA, 1+IR_EXTRA)(1) >>  SHIFT);
  CHECK_AP::checkCalc("DET", det, trueDet, 0.00001);
  // Precision of this (controlled by BDET) is critical.
  //  CHECK_AP::checkCalc("INVDET", invDet_short,
  //                  trueInvDet/double(AP_UFIXED(1-MIN_SHIFT+MAX_SHIFT, 1-MIN_SHIFT)(1) >>  SHIFT), 0.0001);
  CHECK_AP::checkCalc("INVDET", invDet_short,
                      trueInvDet/double(AP_UFIXED(1+NR_EXTRA, 1+IR_EXTRA)(1) >>  SHIFT), 0.0001);
  CHECK_AP::checkCalc("INVR00", _00, true_ri00, 0.001);
  CHECK_AP::checkCalc("INVR01", _01, true_ri01, 0.001);
  CHECK_AP::checkCalc("INVR11", _11, true_ri11, 0.001);
#endif
}

#ifdef CMSSW_GIT_HASH
}

}
#endif

