/**
 * This is the top-level HLS function, which updates a helix state by adding a stub to it.
 * N.B. It therefore can't use the Settings class or any external libraries! Nor can it be a C++ class.
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */

 
#ifndef __KalmanUpdateHLS__
#define __KalmanUpdateHLS__

// Defines StateHLS & KFstateHLS. Also defines finite bit integers & floats.
#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#include "KalmanMatricesHLS.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Add stub to old KF helix state to get new KF helix state.
void KalmanUpdateHLS(const StubHLS& stub, const KFstateHLS& stateIn, KFstateHLS& stateOut, ExtraOutHLS& extraOut);

// Calculate increase in chi2 from adding new stub: delta(chi2) = res(transpose) * R(inverse) * res
TCHI_INT calcDeltaChi2(const VectorRes& res, const MatrixInverseR& Rinv);

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif




