#ifndef __MINIHTSTAGE_H__
#define __MINIHTSTAGE_H__

#include "L1Trigger/TrackFindingTMTT/interface/HTrphi.h"

using namespace std;
using boost::numeric::ublas::matrix;

namespace TMTT {

class Settings;

class MiniHTstage {

public:

  MiniHTstage( const Settings* settings );

  ~MiniHTstage(){};

  void exec( matrix< HTrphi >& mHtRphis ) const;

private:

  const Settings* settings_; // Configuration parameters
  bool miniHTstage_;
  unsigned int muxOutputsHT_;
  unsigned int miniHoughNbinsPt_;
  unsigned int miniHoughNbinsPhi_;
  float miniHoughMinPt_;
  bool miniHoughDontKill_;
  float miniHoughDontKillMinPt_;
  unsigned int numSubSecsEta_;
  unsigned int numPhiOctants_;
  unsigned int numPhiSecPerOct_;
  unsigned int numEtaRegions_;
  bool busySectorKill_;
  unsigned int busySectorNumStubs_;
  vector< unsigned int > busySectorMbinRanges_;
  float chosenRofPhi_;
  float binSizeQoverPtAxis_;
  float binSizePhiTrkAxis_;
  float invPtToDphi_;
};

}

#endif

