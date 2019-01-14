/**
 * This is the top-level HLS function within CMSSW, which updates a helix state by adding a stub to it.
 * N.B. It therefore can't use the Settings class or any external libraries! Nor can it be a C++ class.
 *
 * All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 * 
 * Author: Ian Tomalin
 */

 
#ifndef __KalmanUpdateHLS__
#define __KalmanUpdateHLS__

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS4.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS5.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#include "KalmanMatricesHLS.h"
#include "KalmanMatricesHLS4.h"
#include "KalmanMatricesHLS5.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Internal interface.
// Add stub to old KF helix state to get new KF helix state for NPAR = 4 or 5 param helix fits.
template <unsigned int NPAR>
void kalmanUpdateHLS(const StubHLS& stub, const KFstateHLS<NPAR>& stateIn, KFstateHLS<NPAR>& stateOut, ExtraOutHLS<NPAR>& extraOut);

// Calculate increase in chi2 from adding new stub: delta(chi2) = res(transpose) * R(inverse) * res
template <unsigned int NPAR>
TCHI_INT calcDeltaChi2(const VectorRes<NPAR>& res, const MatrixInverseR<NPAR>& Rinv);

// Set output helix params & associated cov matrix related to d0, & check if d0 passes cut.
// (Relevant only to 5-param helix fit)
template <unsigned int NPAR>
void setOutputsD0(const VectorX<NPAR>& x_new, const MatrixC<NPAR>& C_new, const AP_UINT(3)& nStubs, KFstateHLS<NPAR>& stateOut, ExtraOutHLS<NPAR>& extraOut);

// Fully specialized function templates must also be declared to ensure they are found.

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif




