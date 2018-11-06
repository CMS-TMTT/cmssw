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

#ifdef PRINT_SUMMARY
//#define PRINT_HLSARGS
#endif

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

//=== Add stub to old KF helix state to get new KF helix state.

void KalmanUpdateHLS(const StubHLS& stub, const KFstateHLS& stateIn, KFstateHLS& stateOut, ExtraOutHLS& extraOut) {
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
    std::cout<<std::endl<<"KF HLS bodge bits: V="<<BODGE_V<<" S="<<BODGE_S<<" R="<<BODGE_R<<" IR="<<BODGE_IR<<" DET="<<BODGE_DET<<" K="<<BODGE_K<<" RES="<<BODGE_RES<<" CHI2="<<BODGE_CHI2<<std::endl<<std::endl;
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
  VectorX x(stateIn.inv2R, stateIn.phi0, stateIn.tanL, stateIn.z0);

  // Store covariance matrix of input helix params.
  MatrixC C(stateIn.cov_00, stateIn.cov_11, stateIn.cov_22, stateIn.cov_33, stateIn.cov_01, stateIn.cov_23);

  // Calculate matrix of derivatives of predicted stub coords w.r.t. helix params.
  MatrixH H(stub.r);

  // Calculate S = H*C, and its transpose St, which is equal to C*H(transpose).
  MatrixS           S(H, C);
  MatrixS_transpose St(S);

  // Calculate covariance matrix of predicted residuals R = V + H*C*Ht = V + H*St, and its inverse.
  // (Call this Rmat instead of R to avoid confusion with track radius).
  MatrixR        Rmat(V, H, St);
  MatrixInverseR RmatInv(Rmat);

  // Calculate Kalman gain matrix * determinant(R): K = S*R(inverse)
  MatrixK K(St, RmatInv);

  // Calculate hit residuals.
  VectorRes res(m, H, x); 

  // Calculate output helix params & their covariance matrix.
  VectorX x_new(x, K, res);
  MatrixC C_new(C, K, S);
 
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
  bool innerLayer = (nStubs < 2);

  // IRT - feed in test helix params to debug cases seen in QuestaSim. (1/2r, phi, tanl, z0)
  //x_new._0 = float(-8163)/float(1 << (B18 - BH0));
  //x_new._1 = float(-57543)/float(1 << (B18 - BH1));
  //x_new._2 = float(4285)/float(1 << (B18 - BH2));
  //x_new._3 = float(-7652)/float(1 << (B18 - BH3));

  extraOut.z0Cut = (x_new._3 >= (- z0_digi_cut)   && x_new._3 <= z0_digi_cut)    || innerLayer;
  extraOut.ptCut = (x_new._0 >= (-inv2R_digi_cut[nStubs]) && x_new._0 <= inv2R_digi_cut[nStubs]) || innerLayer;
  extraOut.chiSquaredCut = (chi2 <= chi2_digi_cut[nStubs]) || innerLayer;
  extraOut.sufficientPScut = not (nStubs <= 2 && V._2Smodule); 

  KFstateHLS::TP phiAtRefR = x_new._1 - chosenRofPhi_digi * x_new._0;
  StubHLS::TZ      zAtRefR = x_new._3 + chosenRofZ_digi * x_new._2; // Intentional use of StubHLS::TZ type
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
  extraOut.htBinWithin1Cut = ((mBinHelix_tmp_trunc == stateIn.mBin || mBinHelix_tmp_trunc == mPlus1 || mBinHelix_tmp_trunc == mMinus1) && (cBinHelix_tmp_trunc == stateIn.cBin || cBinHelix_tmp_trunc == cPlus1 || cBinHelix_tmp_trunc == cMinus1));

  //std::cout<<"ZCALC "<<x_new._3<<" "<<chosenRofZ_digi<<" "<<x_new._2<<std::endl;

  // IRT -- feed in test helix params to debug cases seen in QuestaSim.
  //std::cout<<"ZZZ RANGE TMP "<<etaBounds.z_[TMPS]<<" < "<<zAtRefR<<" < "<<etaBounds.z_[TMPS+1]<<" sec="<<TMPS<<" zsign="<<TMPSIGN<<std::endl;

  //std::cout<<"ZZZ RANGE "<<etaBounds.z_[stateIn.etaSectorID]<<" < "<<zAtRefR<<" < "<<etaBounds.z_[stateIn.etaSectorID+1]<<" sec="<<stateIn.etaSectorID<<" zsign="<<stateIn.etaSectorZsign<<std::endl;

  //std::cout<<"CHECK HT WITHIN 1 BIN: "<<extraOut.htBinWithin1Cut<<std::endl;
  //std::cout<<"CHECK IN RANGE: c"<<cBinInRange<<" sec "<<inEtaSector<<std::endl;
  
  //std::cout<<"EXTRA: z0Cut="<<extraOut.z0Cut<<" ptCut="<<extraOut.ptCut<<" chi2Cut="<<extraOut.chiSquaredCut<<" PScut="<<extraOut.sufficientPScut<<std::endl;
  //std::cout<<"EXTRA: mBin="<<int(stateIn.mBin)<<" "<<int(mBinHelix_tmp)<<" cBin="<<int(stateIn.cBin)<<" "<<int(cBinHelix_tmp)<<" consistent="<<extraOut.consistent<<std::endl;
  //std::cout<<"EXTRA: in sector="<<extraOut.sectorCut<<" in eta="<<inEtaSector<<" phiAtR="<<phiAtRefR<<" zAtR="<<zAtRefR<<std::endl;
  
#ifdef PRINT_HLSARGS
  std::cout<<"HLS INPUT stub: r="<<stub.r<<" phiS="<<stub.phiS<<" z="<<stub.z/2<<std::endl;
  std::cout<<"HLS INPUT: HT (m,c)=("<<stateIn.mBin<<","<<stateIn.cBin<<")"
           <<" layers (ID, skip)=("<<stateIn.layerID<<","<<stateIn.nSkippedLayers<<")"
           <<" 1/2R="<<ap_fixed<B18,B18>(stateIn.inv2R.range( B18 - 1, 0))
	   <<" phi0="<<ap_fixed<B18,B18>(stateIn.phi0.range( B18 - 1, 0))
	   <<" tanL="<<ap_fixed<B18,B18>(stateIn.tanL.range( B18 - 1, 0))
	   <<" z0="  <<ap_fixed<B18,B18>(stateIn.z0.range( B18 - 1, 0))
	   <<" chi2="<<ap_ufixed<B17,B17>(stateIn.chiSquared.range( B17 - 1, 0))
	   <<std::endl;
  std::cout<<"HLS INPUT cov:"
           <<" cov00="<<ap_fixed<B25,B25>(stateIn.cov_00.range( B25 - 1, 0))
           <<" cov11="<<ap_fixed<B25,B25>(stateIn.cov_11.range( B25 - 1, 0))
           <<" cov22="<<ap_fixed<B25,B25>(stateIn.cov_22.range( B25 - 1, 0))
           <<" cov33="<<ap_fixed<B25,B25>(stateIn.cov_33.range( B25 - 1, 0))
           <<" cov01="<<ap_fixed<B18,B18>(stateIn.cov_01.range( B18 - 1, 0))
           <<" cov23="<<ap_fixed<B18,B18>(stateIn.cov_23.range( B18 - 1, 0))
	   <<std::endl;
  std::cout<<"HLS OUTPUT: HT (m,c)=("<<stateOut.mBin<<","<<stateOut.cBin<<")"
           <<" layers (ID, skip)=("<<stateOut.layerID<<","<<stateOut.nSkippedLayers<<")"
           <<" 1/2R="<<ap_fixed<B18,B18>(stateOut.inv2R.range( B18 - 1, 0))
	   <<" phi0="<<ap_fixed<B18,B18>(stateOut.phi0.range( B18 - 1, 0))
	   <<" tanL="<<ap_fixed<B18,B18>(stateOut.tanL.range( B18 - 1, 0))
	   <<" z0="  <<ap_fixed<B18,B18>(stateOut.z0.range( B18 - 1, 0))
	   <<" chi2="<<ap_ufixed<B17,B17>(stateOut.chiSquared.range( B17 - 1, 0))
	   <<std::endl;
  std::cout<<"HLS OUTPUT cov:"
           <<" cov00="<<ap_fixed<B25,B25>(stateOut.cov_00.range( B25 - 1, 0))
           <<" cov11="<<ap_fixed<B25,B25>(stateOut.cov_11.range( B25 - 1, 0))
           <<" cov22="<<ap_fixed<B25,B25>(stateOut.cov_22.range( B25 - 1, 0))
           <<" cov33="<<ap_fixed<B25,B25>(stateOut.cov_33.range( B25 - 1, 0))
           <<" cov01="<<ap_fixed<B18,B18>(stateOut.cov_01.range( B18 - 1, 0))
           <<" cov23="<<ap_fixed<B18,B18>(stateOut.cov_23.range( B18 - 1, 0))
	   <<std::endl;
  std::cout<<"HLS OUTPUT EXTRA:"
           <<" Helix (m,c)=("<<extraOut.mBinHelix<<","<<extraOut.cBinHelix<<")"
	   <<std::endl;
#endif

}

// Calculate increase in chi2 from adding new stub: delta(chi2) = res(transpose) * R(inverse) * res
TCHI_INT calcDeltaChi2(const VectorRes& res, const MatrixInverseR& Rinv) {
  // Simplify calculation by noting that Rinv is symmetric.
#ifdef USE_FIXED
  typedef MatrixInverseR::TRI00_short TRI00_short;
  typedef MatrixInverseR::TRI11_short TRI11_short;
  typedef MatrixInverseR::TRI01_short TRI01_short;
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

#ifdef CMSSW_GIT_HASH
}

}
#endif

