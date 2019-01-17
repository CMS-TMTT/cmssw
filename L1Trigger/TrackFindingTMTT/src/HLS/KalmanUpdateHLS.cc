/**
 * This is the top-level HLS function, which updates a helix state by adding a stub to it.
 * N.B. It therefore can't use the Settings class or any external libraries! Nor can it be a C++ class.
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanUpdateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSconstants.h"
#else
#include "KalmanUpdateHLS.h"
#include "KalmanMatricesHLS.h"
#include "HLSutilities.h"
#include "HLSconstants.h"
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

template void kalmanUpdateHLS(const StubHLS& stub, const KFstateHLS<4>& stateIn, KFstateHLS<4>& stateOut, ExtraOutHLS<4>& extraOut);

template void kalmanUpdateHLS(const StubHLS& stub, const KFstateHLS<5>& stateIn, KFstateHLS<5>& stateOut, ExtraOutHLS<5>& extraOut);

template TCHI_INT calcDeltaChi2(const VectorRes<4>& res, const MatrixInverseR<4>& Rinv);

template TCHI_INT calcDeltaChi2(const VectorRes<5>& res, const MatrixInverseR<5>& Rinv);

//=== Add stub to old KF helix state to get new KF helix state.

template <unsigned int NPAR>
void kalmanUpdateHLS(const StubHLS& stub, const KFstateHLS<NPAR>& stateIn, KFstateHLS<NPAR>& stateOut, ExtraOutHLS<NPAR>& extraOut) {
  stateOut.cBin = stateIn.cBin;
  stateOut.mBin = stateIn.mBin;
  stateOut.layerID = stateIn.layerID;
  stateOut.nSkippedLayers = stateIn.nSkippedLayers;
  stateOut.candidateID = stateIn.candidateID;
  stateOut.eventID = stateIn.eventID;
  stateOut.etaSectorID = stateIn.etaSectorID;
  stateOut.etaSectorZsign = stateIn.etaSectorZsign;
  stateOut.valid = (stateIn.valid && stateIn.valid);

#ifdef PRINT_SUMMARY
  static bool first = true;
  if (first) {
    first = false;
    std::cout<<std::endl<<"KF HLS bodge bits: V="<<BODGE<NPAR>::V<<" S="<<BODGE<NPAR>::S<<" R="<<BODGE<NPAR>::R<<" IR="<<BODGE<NPAR>::IR<<" DET="<<BODGE<NPAR>::DET<<" K="<<BODGE<NPAR>::K<<" RES="<<BODGE<NPAR>::RES<<" CHI2="<<BODGE<NPAR>::CHI2<<std::endl<<std::endl;
  }
#endif

#ifdef PRINT
  std::cout<<"KalmanUpdate call: layerID="<<stateIn.layerID<<" nSkipped="<<stateIn.nSkippedLayers<<std::endl;
#endif

  // Store vector of stub coords.
  VectorM m(stub.phiS, stub.z);

  // Store covariance matrix of stub coords.
  MatrixV V(stub.r, stub.z, stateIn.inv2R, stateIn.tanL, stateIn.mBin);

  // Store vector of input helix params.
  VectorX<NPAR> x(stateIn);

  // Store covariance matrix of input helix params.
  MatrixC<NPAR> C(stateIn);

  // Calculate matrix of derivatives of predicted stub coords w.r.t. helix params.
  MatrixH<NPAR> H(stub.r);

  // Calculate S = H*C, and its transpose St, which is equal to C*H(transpose).
  MatrixS<NPAR>           S(H, C);
  MatrixS_transpose<NPAR> St(S);

  // Calculate covariance matrix of predicted residuals R = V + H*C*Ht = V + H*St, and its inverse.
  // (Call this Rmat instead of R to avoid confusion with track radius).
  MatrixR<NPAR>  Rmat(V, H, St);
  MatrixInverseR<NPAR> RmatInv(Rmat);

  // Calculate Kalman gain matrix * determinant(R): K = S*R(inverse)
  MatrixK<NPAR> K(St, RmatInv);

  // Calculate hit residuals.
  VectorRes<NPAR> res(m, H, x); 

  // Calculate output helix params & their covariance matrix.
  VectorX<NPAR> x_new(x, K, res);
  MatrixC<NPAR> C_new(C, K, S);
 
  /*
  // Useful to debug C matrices with negative determinants, by fully recalculating them in double precision.
  double s00 = double(H._00) * double(C._00) + double(H._01) * double(C._10) + double(H._02) * double(C._20) + double(H._03) * double(C._30);
  double s01 = double(H._00) * double(C._01) + double(H._01) * double(C._11) + double(H._02) * double(C._21) + double(H._03) * double(C._31);
  double s02 = double(H._00) * double(C._02) + double(H._01) * double(C._12) + double(H._02) * double(C._22) + double(H._03) * double(C._32);
  double s03 = double(H._00) * double(C._03) + double(H._01) * double(C._13) + double(H._02) * double(C._23) + double(H._03) * double(C._33);
  double s10 = double(H._10) * double(C._00) + double(H._11) * double(C._10) + double(H._12) * double(C._20) + double(H._13) * double(C._30);
  double s11 = double(H._10) * double(C._01) + double(H._11) * double(C._11) + double(H._12) * double(C._21) + double(H._13) * double(C._31);
  double s12 = double(H._10) * double(C._02) + double(H._11) * double(C._12) + double(H._12) * double(C._22) + double(H._13) * double(C._32);
  double s13 = double(H._10) * double(C._03) + double(H._11) * double(C._13) + double(H._12) * double(C._23) + double(H._13) * double(C._33);
  double st00 = s00;
  double st10 = s01;
  double st20 = s02;
  double st30 = s03;
  double st01 = s10;
  double st11 = s11;
  double st21 = s12;
  double st31 = s13;  
  double r00 = double(V._00) + double(H._00)*(st00) + double(H._01)*(st10) + double(H._02)*(st20) + double(H._03)*(st30);
  double r11 = double(V._11) + double(H._10)*(st01) + double(H._11)*(st11) + double(H._12)*(st21) + double(H._13)*(st31);
  double rinv00 = 1./r00;
  double rinv11 = 1./r11;
  double k00 =  (st00)*rinv00; 
  double k10 =  (st10)*rinv00;
  double k20 =  (st20)*rinv00;
  double k30 =  (st30)*rinv00;
  double k01 =  (st01)*rinv11;
  double k11 =  (st11)*rinv11;
  double k21 =  (st21)*rinv11;
  double k31 =  (st31)*rinv11;
  double c22 =  double(C._22) - (k20 * (s02) + k21 * (s12));
  double c33 =  double(C._33) - (k30 * (s03) + k31 * (s13));
  double c23 =  double(C._23) - (k20 * (s03) + k21 * (s13));
  std::cout<<"recalc C new: TT="<<c22<<" ZZ="<<c33<<" TZ="<<c23<<std::endl;
  CHECK_AP::checkDet("recalc_rz",c22,c33,c23); 
 */

  // Calculate increase in chi2 from adding new stub: delta(chi2) = res(transpose) * R(inverse) * res
  TCHI_INT deltaChi2 = calcDeltaChi2(res, RmatInv);
  TCHI_INT chi2 = stateIn.chiSquared + deltaChi2;
  // Truncate chi2 to avoid overflow.
  static const TCHI_INT MAX_CHI2 = (1 << BCHI) - 1;
  if (chi2 > MAX_CHI2) chi2 = MAX_CHI2;
  stateOut.chiSquared = chi2;
  
  stateOut.inv2R = x_new._0;
  stateOut.phi0  = x_new._1;
  stateOut.tanL  = x_new._2;
  stateOut.z0    = x_new._3;
  stateOut.cov_00 = C_new._00;
  stateOut.cov_11 = C_new._11;
  stateOut.cov_22 = C_new._22;
  stateOut.cov_33 = C_new._33;
  stateOut.cov_01 = C_new._01;
  stateOut.cov_23 = C_new._23;

  // Check if output helix passes cuts.
  // (Copied from Maxeller code KFWorker.maxj)
  AP_UINT(3) nStubs = stateIn.layerID - stateIn.nSkippedLayers; // Number of stubs on state including current one.

  // IRT - feed in test helix params to debug cases seen in QuestaSim. (1/2r, phi, tanl, z0)
  //x_new._0 = float(-8163)/float(1 << (B18 - BH0));
  //x_new._1 = float(-57543)/float(1 << (B18 - BH1));
  //x_new._2 = float(4285)/float(1 << (B18 - BH2));
  //x_new._3 = float(-7652)/float(1 << (B18 - BH3));

  KFstateHLS<5>::TZ   cut_z0          = z0Cut[nStubs];
  KFstateHLS<5>::TZ   cut_z0_minus    = z0CutMinus[nStubs];
  KFstateHLS<5>::TR   cut_inv2R       = inv2Rcut[nStubs];
  KFstateHLS<5>::TR   cut_inv2R_minus = inv2RcutMinus[nStubs];
  KFstateHLS<5>::TCHI cut_chi2        = chi2Cut[nStubs];
  // Don't do "hls::abs(x_new._3) <= cut_z0)" as this wastes 2 clk cycles.
  // Also, don't do "cut_z0_minus = - cut_z0" or this fails Vivado implementation with timing errors.
  extraOut.z0Cut = ((x_new._3 >= cut_z0_minus    && x_new._3 <= cut_z0)    || (cut_z0 == 0));  // cut = 0 means cut not applied.
  extraOut.ptCut = ((x_new._0 >= cut_inv2R_minus && x_new._0 <= cut_inv2R) || (cut_inv2R == 0)); 
  extraOut.chiSquaredCut = ((chi2 <= cut_chi2)        || (cut_chi2 == 0));
  extraOut.sufficientPScut = not (nStubs <= 2 && V._2Smodule);
  // IRT -- very useful whilst optimising variable bit ranges, to skip all but first iteration.
  //extraOut.ptCut = false;

  //=== Set output helix params & associated cov matrix related to d0, & check if d0 passes cut.
  //=== (Relevant only to 5-param helix fit) 
  setOutputsD0(x_new, C_new, nStubs, stateOut, extraOut);

  typename KFstateHLS<NPAR>::TP phiAtRefR = x_new._1 - chosenRofPhi * x_new._0;
  StubHLS::TZ                   zAtRefR = x_new._3 + chosenRofZ * x_new._2; // Intentional use of StubHLS::TZ type
  // Constants BMH & BCH below set in HLSconstants.h
  // Casting from ap_fixed to ap_int rounds to zero, not -ve infinity, so cast to ap_fixed with no fractional part first.
  AP_INT(BMH) mBinHelix_tmp = AP_FIXED(BMH,BMH)( 
				AP_FIXED(B18+invRToMbin_bitShift,BH0+invRToMbin_bitShift)(x_new._0) << invRToMbin_bitShift
					       );
  AP_INT(BCH) cBinHelix_tmp = AP_FIXED(BCH,BCH)(
						AP_FIXED(B18+phiToCbin_bitShift,BH1)(phiAtRefR) >> phiToCbin_bitShift
					       );
  bool cBinInRange = (cBinHelix_tmp >= minPhiBin && cBinHelix_tmp <= maxPhiBin);

  // Duplicate removal works best in mBinHelix is forced back into HT array if it lies just outside.
  AP_INT(BHT_M) mBinHelix_tmp_trunc;
  if (mBinHelix_tmp < minPtBin) {
    mBinHelix_tmp_trunc = minPtBin;
  } else if (mBinHelix_tmp > maxPtBin) {
    mBinHelix_tmp_trunc = maxPtBin;
  } else {
    mBinHelix_tmp_trunc = mBinHelix_tmp;
  }
  AP_INT(BHT_C) cBinHelix_tmp_trunc = cBinHelix_tmp;
  extraOut.mBinHelix = mBinHelix_tmp_trunc;
  extraOut.cBinHelix = cBinHelix_tmp_trunc;
  //std::cout<<"MBIN helix "<<extraOut.mBinHelix<<" tmp "<<mBinHelix_tmp<<" ht "<<stateIn.mBin<<std::endl;
  //std::cout<<"CBIN helix "<<extraOut.cBinHelix<<" tmp "<<cBinHelix_tmp<<" ht "<<stateIn.cBin<<std::endl;

  static const EtaBoundaries etaBounds;

  // IRT -- feed in test helix params to debug cases seen in QuestaSim.
  //bool TMPSIGN = false;
  //if (TMPSIGN) zAtRefR = -zAtRefR;
  //unsigned int TMPS = 1;
  //bool inEtaSector = (zAtRefR > etaBounds.z_[TMPS] && zAtRefR < etaBounds.z_[TMPS+1]);

  if (stateIn.etaSectorZsign == 1) zAtRefR = -zAtRefR;
  bool inEtaSector = (zAtRefR > etaBounds.z_[stateIn.etaSectorID] && zAtRefR < etaBounds.z_[stateIn.etaSectorID+1]);

  extraOut.sectorCut = (cBinInRange && inEtaSector);
  extraOut.consistent = (mBinHelix_tmp_trunc == stateIn.mBin && cBinHelix_tmp_trunc == stateIn.cBin);
  // The following long-winded calc. saves a clock cycle.
  AP_INT(BHT_M) mPlus1  = stateIn.mBin + AP_INT(BHT_M)(1);
  AP_INT(BHT_M) mMinus1 = stateIn.mBin - AP_INT(BHT_M)(1);
  AP_INT(BHT_C) cPlus1  = stateIn.cBin + AP_INT(BHT_C)(1);
  AP_INT(BHT_C) cMinus1 = stateIn.cBin - AP_INT(BHT_C)(1);
  // This cut is not needed within KF HLS. Will be done in KF or DR VDHL. 
  //extraOut.htBinWithin1Cut = ((mBinHelix_tmp_trunc == stateIn.mBin || mBinHelix_tmp_trunc == mPlus1 || mBinHelix_tmp_trunc == mMinus1) && (cBinHelix_tmp_trunc == stateIn.cBin || cBinHelix_tmp_trunc == cPlus1 || cBinHelix_tmp_trunc == cMinus1));
  extraOut.htBinWithin1Cut = true;

  //std::cout<<"ZCALC "<<x_new._3<<" "<<chosenRofZ<<" "<<x_new._2<<std::endl;

  // IRT -- feed in test helix params to debug cases seen in QuestaSim.
  //std::cout<<"ZZZ RANGE TMP "<<etaBounds.z_[TMPS]<<" < "<<zAtRefR<<" < "<<etaBounds.z_[TMPS+1]<<" sec="<<TMPS<<" zsign="<<TMPSIGN<<std::endl;

  //std::cout<<"ZZZ RANGE "<<etaBounds.z_[stateIn.etaSectorID]<<" < "<<zAtRefR<<" < "<<etaBounds.z_[stateIn.etaSectorID+1]<<" sec="<<stateIn.etaSectorID<<" zsign="<<stateIn.etaSectorZsign<<std::endl;

  //std::cout<<"CHECK HT WITHIN 1 BIN: "<<extraOut.htBinWithin1Cut<<std::endl;
  //std::cout<<"CHECK IN RANGE: c"<<cBinInRange<<" sec "<<inEtaSector<<std::endl;
  
  //std::cout<<"EXTRA: z0Cut="<<extraOut.z0Cut<<" ptCut="<<extraOut.ptCut<<" chi2Cut="<<extraOut.chiSquaredCut<<" PScut="<<extraOut.sufficientPScut<<std::endl;
  //std::cout<<"EXTRA: mBin="<<int(stateIn.mBin)<<" "<<int(mBinHelix_tmp)<<" cBin="<<int(stateIn.cBin)<<" "<<int(cBinHelix_tmp)<<" consistent="<<extraOut.consistent<<std::endl;
  //std::cout<<"EXTRA: in sector="<<extraOut.sectorCut<<" in eta="<<inEtaSector<<" phiAtR="<<phiAtRefR<<" zAtR="<<zAtRefR<<std::endl;
  
#ifdef PRINT_HLSARGS
  stub.print("HLS INPUT stub:");
  stateIn.print("HLS INPUT state:");
  stateOut.print("HLS OUTPUT state:");
  extraOut.print("HLS OUTPUT extra:");
#endif
}

//=== Calculate increase in chi2 from adding new stub: delta(chi2) = res(transpose) * R(inverse) * res

template <unsigned int NPAR>
TCHI_INT calcDeltaChi2(const VectorRes<NPAR>& res, const MatrixInverseR<NPAR>& Rinv) {
  // Simplify calculation by noting that Rinv is symmetric.
#ifdef USE_FIXED
  typedef typename MatrixInverseR<NPAR>::TRI00_short TRI00_short;
  typedef typename MatrixInverseR<NPAR>::TRI11_short TRI11_short;
  typedef typename MatrixInverseR<NPAR>::TRI01_short TRI01_short;
  TCHI_INT dChi2 = (res._0 * res._0) * TRI00_short(Rinv._00) +
                   (res._1 * res._1) * TRI11_short(Rinv._11) +
               2 * (res._0 * res._1) * TRI01_short(Rinv._01);
#else
  TCHI_INT dChi2 = SW_FLOAT(res._0) * Rinv._00 * SW_FLOAT(res._0) + 
                   SW_FLOAT(res._1) * Rinv._11 * SW_FLOAT(res._1) +
	      2 * (SW_FLOAT(res._0) * Rinv._01 * SW_FLOAT(res._1));
#endif
#ifdef PRINT_SUMMARY
  double chi2_phi = double(res._0) * double(res._0) * double(Rinv._00);
  double chi2_z   = double(res._1) * double(res._1) * double(Rinv._11);
  double chi2_c   = double(res._0) * double(res._1) * double(Rinv._01);
  CHECK_AP::checkCalc("dchi2", dChi2, chi2_phi + chi2_z + 2*chi2_c, 0.1, 0.1);
#ifdef PRINT
  std::cout<<"Delta chi2 = "<<dChi2<<" res (phi,z) = "<<res._0<<" "<<res._1<<" chi2 (phi,z) = "<<chi2_phi<<" "<<chi2_z<<std::endl;
#endif
#endif
  return dChi2;
}

//=== Set output helix params & associated cov matrix related to d0, & check if d0 passes cut.
//=== (Relevant only to 5-param helix fit)

void setOutputsD0(const VectorX<4>& x_new, const MatrixC<4>& C_new, const AP_UINT(3)& nStubs, KFstateHLS<4>& stateOut, ExtraOutHLS<4>& extraOut) {}

void setOutputsD0(const VectorX<5>& x_new, const MatrixC<5>& C_new, const AP_UINT(3)& nStubs, KFstateHLS<5>& stateOut, ExtraOutHLS<5>& extraOut) {
  stateOut.d0 = x_new._4;
  stateOut.cov_44 = C_new._44;
  stateOut.cov_04 = C_new._04;
  stateOut.cov_14 = C_new._14;
  KFstateHLS<5>::TD cut_d0        = d0Cut[nStubs];
  KFstateHLS<5>::TD cut_d0_minus  = d0CutMinus[nStubs];
  extraOut.d0Cut = ((x_new._4 >= cut_d0_minus && x_new._4 <= cut_d0) || (cut_d0 == 0));
}

#ifdef CMSSW_GIT_HASH
}

}
#endif

