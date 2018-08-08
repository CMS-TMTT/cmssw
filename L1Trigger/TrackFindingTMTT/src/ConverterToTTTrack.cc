#include "L1Trigger/TrackFindingTMTT/interface/ConverterToTTTrack.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrk4and5.h"
#include "FWCore/Utilities/interface/Exception.h"

//=== Convert our non-persistent L1track3D object (track candidate found by Hough transform prior to fit)
//=== to the official persistent CMSSW EDM TTrack format.

namespace TMTT {


TTTrack< Ref_Phase2TrackerDigi_ > ConverterToTTTrack::makeTTTrack(const L1track3D& trk, unsigned int iPhiSec, unsigned int iEtaReg) const {

  // Get references to stubs on this track.
  std::vector<TTStubRef> ttstubrefs = this->getStubRefs(trk); 

  // Create TTTrack object using these stubs.
  TTTrack< Ref_Phase2TrackerDigi_ > track(ttstubrefs);

  // Note which (eta,phi) sector this track was reconstructed in by HT.
  track.setSector(iPhiSec);
  track.setWedge(iEtaReg);

  // Set helix parameters.
  const unsigned int nPar4 = 4; // Number of helix parameters determined by HT.
  GlobalPoint bsPosition(0.0, 0.0, trk.z0()); // Point of closest approach of track to beam-line.
  track.setPOCA(bsPosition, nPar4);
  float pt = 1./trk.invPt(); // pt
  track.setMomentum(
		    GlobalVector(
				 GlobalVector::Cylindrical(
							   pt, 
							   trk.phi0(), // phi
							   pt*trk.tanLambda()  // pz
							   )
				 ),
		    nPar4
		    );
  track.setRInv(invPtToInvR_ * trk.qOverPt(), nPar4);
  track.setChi2(-1, nPar4); // meaningless since no track fit was done.
  track.setStubPtConsistency(-1, nPar4); // not filled.

  return track;
}

//=== Convert our non-persistent L1fittedTrack object (fitted track candidate)
//=== to the official persistent CMSSW EDM TTrack format.

TTTrack< Ref_Phase2TrackerDigi_ > ConverterToTTTrack::makeTTTrack(const L1fittedTrack& trk, unsigned int iPhiSec, unsigned int iEtaReg) const{

  // Check that this track is valid.
  if (! trk.accepted()) throw cms::Exception("ConverterToTTTrack ERROR: requested to convert invalid L1fittedTrack.");

  // Get references to stubs on this track.
  std::vector<TTStubRef> ttstubrefs = this->getStubRefs(trk);

  // Create TTTrack object using these stubs. 
  TTTrack< Ref_Phase2TrackerDigi_ > track(ttstubrefs);

  // Note which (eta,phi) sector this track was reconstructed in by HT.
  track.setSector(iPhiSec);
  track.setWedge(iEtaReg);

  // Set helix parameters.
  const unsigned int nPar = trk.nHelixParam(); // Number of helix parameters determined by HT.
  // Get and store fitted track parameters 
  GlobalPoint bsPosition(-trk.d0()*sin(trk.phi0()),
			  trk.d0()*cos(trk.phi0()),
			  trk.z0()); // Point of closest approach of track to beam-line.
  track.setPOCA(bsPosition, nPar);
  float pt = trk.pt(); // pt
  track.setMomentum(
		    GlobalVector(
				 GlobalVector::Cylindrical(
							   pt,
							   trk.phi0(), // phi
							   pt*trk.tanLambda()  // pz
							   )
				 ),
		    nPar
		    );
  track.setRInv(invPtToInvR_ * trk.qOverPt(), nPar);
  track.setChi2(trk.chi2(), nPar);
  track.setStubPtConsistency(-1, nPar);

  return track;
}

//=== Convert our non-persistent L1fittedTrk4and5 object (track candidate with both 4 and 5 parameter helix fits)
//=== to the official persistent CMSSW EDM TTrack format.

TTTrack< Ref_Phase2TrackerDigi_ > ConverterToTTTrack::makeTTTrack(const L1fittedTrk4and5& trk4and5, unsigned int iPhiSec, unsigned int iEtaReg) const {

  // Check if 4 and 5 parameter track fits are valid
  bool valid4par = trk4and5.validL1fittedTrack(4);
  bool valid5par = trk4and5.validL1fittedTrack(5);

  // Check that at least one is valid.
  if (! (valid4par || valid5par) ) throw cms::Exception("ConverterToTTTrack ERROR: requested to convert invalid L1fittedTrk4and5.");

  // Get references to stubs on fitted track.
  // The official EDM TTrack object doesn't allow the 4 & 5 parameter fits to use different stubs, so use either of them.
  std::vector<TTStubRef> ttstubrefs;
  if (valid4par) {
    ttstubrefs = this->getStubRefs(trk4and5.getL1fittedTrack(4)); 
  } else {
    ttstubrefs = this->getStubRefs(trk4and5.getL1fittedTrack(5)); 
  }

  // Create TTTrack object using these stubs.
  TTTrack< Ref_Phase2TrackerDigi_ > track(ttstubrefs);

  // Note which (eta,phi) sector this track was reconstructed in by HT.
  track.setSector(iPhiSec);
  track.setWedge(iEtaReg);

  //--- Copy helix parameters for 4 & 5 parameter fits.

  for (unsigned int iPar = 4; iPar <= 5; iPar++) {
    // Check fit with this number of parameters gave valid track.
    if (trk4and5.validL1fittedTrack(iPar)) {
      // Get track obtained with this number of helix parameters in fit.
      const L1fittedTrack& trk = trk4and5.getL1fittedTrack(iPar); 
      // Get and store track parameters.
      GlobalPoint bsPosition(-trk.d0()*sin(trk.phi0()),
			      trk.d0()*cos(trk.phi0()),
			      trk.z0()); // Point of closest approach of track to beam-line.
      track.setPOCA(bsPosition, iPar);
      float pt = 1./trk.invPt(); // pt
      track.setMomentum(
			GlobalVector(
				     GlobalVector::Cylindrical(
							       pt, 
							       trk.phi0(), // phi
							       pt*trk.tanLambda()  // pz
							       )
				     ),
			iPar
			);
      track.setRInv(invPtToInvR_ * trk.qOverPt(), iPar);
      track.setChi2(trk.chi2(), iPar);
      track.setStubPtConsistency(-1, iPar);
    }
  }

  return track;
}

}