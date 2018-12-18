/**
 * This is the top-level function for Vivado HLS compilation. 
 * It is not used by CMSSW.
 * 
 * It is required because HLS does not allow the top-level function to be templated.
 * 
 * Author: Ian Tomalin
 */

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanUpdateHLS_top.h"
#else
#include "KalmanUpdateHLS_top.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

void kalmanUpdateHLS_top(const StubHLS& stub, const KFstateHLS<N_HELIX_PAR>& stateIn, KFstateHLS<N_HELIX_PAR>& stateOut, ExtraOutHLS<N_HELIX_PAR>& extraOut) {
  kalmanUpdateHLS<N_HELIX_PAR>(stub, stateIn, stateOut, extraOut);
}

#ifdef CMSSW_GIT_HASH
}

}
#endif
