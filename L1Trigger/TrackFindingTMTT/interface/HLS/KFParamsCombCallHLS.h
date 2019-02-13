/**
 * This is the interface between the C++ KF framework in CMSSW and the HLS code.
 * It is identical to KF4ParamsComb, except that the state updator is modified 
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */
 
#ifndef __KFPARAMSCOMBCALLHLS__
#define __KFPARAMSCOMBCALLHLS__
 
#include "L1Trigger/TrackFindingTMTT/interface/KFParamsComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSpragmaOpts.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSconstants.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"

using namespace std;

namespace TMTT {

class TP; 
class kalmanState;
class StubCluster;

class KFParamsCombCallHLS : public KFParamsComb {

public:

  KFParamsCombCallHLS(const Settings* settings, const uint nPar, const string &fitterName ) : KFParamsComb(settings, nPar, fitterName) {
#ifdef PT_2GEV
    if (settings->houghMinPt() > 2.5) throw cms::Exception("KFParamsConmbCallHLS: Edit HLSpragmaOpts.h to undefine PT_2GEV");
#else
    if (settings->houghMinPt() < 2.5) throw cms::Exception("KFParamsConmbCallHLS: Edit HLSpragmaOpts.h to define PT_2GEV");
#endif
  }

  virtual ~KFParamsCombCallHLS(){}

  // Update KF helix params with this stub.
  // (Override KF state updator in L1KalmanComb with version suitable for HLS).
  const kalmanState *kalmanUpdate( unsigned skipped, unsigned layer, const StubCluster* stubCluster, const kalmanState &stateIn, const TP *);

  // Print summary info to help tune bit ranges of HLS calculations.
  virtual void endJob() {KalmanHLS::CHECK_AP::printCheckMap();}

protected:

  // Is this HLS code?
  virtual bool isHLS() {return true;};

private:

  // Get digital stub info that the KF VHDL injects into the KF state updater (Maxeller/HLS)
  KalmanHLS::StubHLS    getDigiStub(const StubCluster* stubCluster, const kalmanState* state);

  // Get digitised KF state info that the KF VHDL injects into the KF state updater (Maxeller/HLS),
  // both for NPAR = 4 & 5 param helix states.
  template <unsigned int NPAR> 
  KalmanHLS::KFstateHLS<NPAR> getDigiStateIn(unsigned int skipped, unsigned int layer, const kalmanState* state) const;

  // Implement NPAR-specific code called by getDigiStateIn(...).
  template <unsigned int NPAR>
  void getDigiStateInUtil(const vector<double>& helixParams, const TMatrixD& cov, KalmanHLS::KFstateHLS<NPAR>& stateDigi) const;

  // Convert digitized ourput KF state to floating point,
  // both for NPAR = 4 & 5 param helix states.
  template <unsigned int NPAR>
  const kalmanState* getStateOut(const kalmanState* stateIn, const StubCluster* stubCluster, const KalmanHLS::KFstateHLS<NPAR>& stateOutDigi, const KalmanHLS::ExtraOutHLS<NPAR>& extraOut);

  // Implement NPAR-specific code call by getStateOut(...).
  template <unsigned int NPAR>
  void getStateOutUtil(const KalmanHLS::KFstateHLS<NPAR>& stateOutDigi, const KalmanHLS::ExtraOutHLS<NPAR>& extraOut, vector<double>& x, TMatrixD& pxx);

  // This is identical to version in KFParamsComb, deciding if a state passes cuts,
  // except that it also checks the cut decisions produced by the HLS KalmanUpdate.
  bool isGoodState( const kalmanState &state ) const;

private:
  // Digitisation multipliers
  double rMult_;
  double zMult_;
  double phiMult_;
  double inv2R_Mult_;
  double d0_Mult_;
  // Reference radius in r-phi plane.
  double chosenRofPhi_;
  // Number of eta sectors.
  unsigned int numEtaRegions_;

  // Store the extra info provided by the HLS updator about whether the state passes cuts.
  KalmanHLS::ExtraOutHLS<4> extraOut4_;
  KalmanHLS::ExtraOutHLS<5> extraOut5_;
};

// Fully specialized templates must be declared outside class, but inside .h file, to ensure they are found.

template <>
void KFParamsCombCallHLS::getDigiStateInUtil<4>(const vector<double>& helixParams, const TMatrixD& cov, KalmanHLS::KFstateHLS<4>& stateDigi) const;

template <>
void KFParamsCombCallHLS::getDigiStateInUtil<5>(const vector<double>& helixParams, const TMatrixD& cov, KalmanHLS::KFstateHLS<5>& stateDigi) const;

template <>
void KFParamsCombCallHLS::getStateOutUtil<4>(const KalmanHLS::KFstateHLS<4>& stateOutDigi, const KalmanHLS::ExtraOutHLS<4>& extraOut, vector<double>& x, TMatrixD& pxx);

template <>
void KFParamsCombCallHLS::getStateOutUtil<5>(const KalmanHLS::KFstateHLS<5>& stateOutDigi, const KalmanHLS::ExtraOutHLS<5>& extraOut, vector<double>& x, TMatrixD& pxx);
}

#endif
