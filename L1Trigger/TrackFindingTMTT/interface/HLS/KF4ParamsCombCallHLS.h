/**
 * This is the interface between the C++ KF framework in CMSSW and the HLS code.
 * It is identical to KF4ParamsComb, except that the state updator is modified 
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */
 
#ifndef __KF4PARAMSCOMBCALLHLS__
#define __KF4PARAMSCOMBCALLHLS__
 
#include "L1Trigger/TrackFindingTMTT/interface/KF4ParamsComb.h"
// Defines StubHLS & KFstateHLS & ExtraOutHLS.
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"

using namespace std;

namespace TMTT {

class TP; 
class kalmanState;
class StubCluster;

class KF4ParamsCombCallHLS : public KF4ParamsComb {

public:

  KF4ParamsCombCallHLS(const Settings* settings, const uint nPar, const string &fitterName ) : KF4ParamsComb(settings, nPar, fitterName) {}

  virtual ~KF4ParamsCombCallHLS(){}

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

  // Get digitised KF state info that the KF VHDL injects into the KF state updater (Maxeller/HLS)
  KalmanHLS::KFstateHLS getDigiStateIn(unsigned int skipped, unsigned int layer, const kalmanState* state);

  // Convert digitized ourput KF state to floating point.
  const kalmanState* getStateOut(const kalmanState* stateIn, const StubCluster* stubCluster, const KalmanHLS::KFstateHLS& stateOutDigi, const KalmanHLS::ExtraOutHLS& extraOut);

  // This is identical to version in KF4ParamsComb, deciding if a state passes cuts,
  // except that it also checks the cut decisions produced by the HLS KalmanUpdate.
  bool isGoodState( const kalmanState &state ) const;

private:
  // Digitisation multipliers
  double rMult_;
  double zMult_;
  double phiMult_;
  double inv2R_Mult_;
  // Reference radius in r-phi plane.
  double chosenRofPhi_;
  // Number of eta sectors.
  unsigned int numEtaRegions_;

  // Store the extra info provided by the HLS updator about whether the state passes cuts.
  KalmanHLS::ExtraOutHLS extraOut_;
};

}

#endif
