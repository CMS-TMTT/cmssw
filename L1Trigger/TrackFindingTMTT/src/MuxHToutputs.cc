
//--- Note that the word "link" appearing in the C++ or comments in this class actually corresponds 
//--- to a pair of links in the hardware.

#include "L1Trigger/TrackFindingTMTT/interface/MuxHToutputs.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"

#include "FWCore/Utilities/interface/Exception.h"

namespace TMTT {

//=== Initialize constants from configuration parameters.

MuxHToutputs::MuxHToutputs(const Settings* settings) : 
  settings_(settings),
  muxOutputsHT_( settings_->muxOutputsHT() ),
  numPhiOctants_( settings_->numPhiOctants() ),
  numPhiSectors_( settings_->numPhiSectors() ),
  numPhiSecPerOct_( numPhiSectors_ / numPhiOctants_ ),
  numEtaRegions_( settings_->numEtaRegions() ),
  busySectorKill_( settings_->busySectorKill() ),              // Kill excess tracks flowing out of HT?
  busySectorNumStubs_( settings_->busySectorNumStubs()),       // Max. num. of stubs that can be sent within TM period
  busySectorMbinRanges_( settings_->busySectorMbinRanges() ),  // Individual m bin (=q/Pt) ranges to be output to opto-links. 
  busySectorUseMbinRanges_( busySectorMbinRanges_.size() > 0)  // m bin ranges option disabled if vector empty. 
{
  // Implemented MUX algorithm relies on same number of sectors per octant.
  if (numPhiSectors_%numPhiOctants_ != 0) throw cms::Exception("MuxHToutputs: Number of phi sectors is not a multiple of number of octants!");

  if ( ! busySectorUseMbinRanges_) throw cms::Exception("MuxHToutputs: The implemented MUX algorithm requires you to be using the busySectorMbinRanges cfg option!");

  // Check that the MUX algorithm implemented in linkID() is not obviously wrong.
  this->sanityCheck();

  bool static first = true;
  if (first) {
    first = false;
    cout<<"=== The R-PHI HT output is multiplexed onto "<<this->numLinksPerOctant()<<" pairs of opto-links per octant."<<endl;
  }
}

//=== Determine which tracks are transmitted on each HT output optical link, taking into account the multiplexing
//=== of multiple (eta,phi) sectors onto single links and the truncation of the tracks caused by the requirement
//=== to output all the tracks within the time-multiplexed period.
//=== This function replaces the 2D track collection in the r-phi HT with the subset surviving the TM cut.

void MuxHToutputs::exec(matrix<HTrphi>& mHtRphis) const {

  // As this loops over sectors in order of increasing sector number, this MUX algorithm always transmits tracks
  // from the lowest sector numbers on each link first. So the highest sector numbers are more likely to be
  // truncated by the TM period. The algorithm assumes that two or more m-bin ranges from the same sector will never
  // be transmitted down the same link, as if this happens, it does not predict the order in which they will be 
  // transmitted.

  for (unsigned int iPhiOct = 0; iPhiOct < numPhiOctants_; iPhiOct++) {
    vector<unsigned int> numStubsPerLink( this->numLinksPerOctant(), 0 );

    for (unsigned int iSecInOct = 0; iSecInOct < numPhiSecPerOct_; iSecInOct++) {
      unsigned int iPhiSec = iPhiOct * numPhiSecPerOct_ + iSecInOct;

      for (unsigned int iEtaReg = 0; iEtaReg < numEtaRegions_; iEtaReg++) {

  	HTrphi& htRphi = mHtRphis(iPhiSec, iEtaReg); // Get a mutable version of the r-phi HT.

	vector<const L1track2D*> keptTracks;
	vector<L1track2D> tracks = htRphi.trackCands2D();

	for (L1track2D& trk : tracks) {
	  unsigned int nStubs = trk.getNumStubs(); // #stubs on this track.
	  unsigned int mBinRange = htRphi.getMbinRange(trk); // Which m bin range is this track in?
	  // Get the output optical link corresponding to this sector & m-bin range.
	  unsigned int link = this->linkID(iSecInOct, iEtaReg, mBinRange);
	  // Make a note of opto-link number inside track object.
	  trk.setOptoLinkID(link);

	  numStubsPerLink[ link ] += nStubs;
	  // Check if this track can be output within the time-multiplexed period.
          bool keep = ( (not busySectorKill_) || (numStubsPerLink[ link ] <= busySectorNumStubs_));
	  if (keep) keptTracks.push_back(&trk);
	}

	// Replace the collection of 2D tracks in the r-phi HT with the subset of them surviving the TM cut.
	htRphi.replaceTrackCands2D(keptTracks);
      }
    } 
  }
}

//=== Define the number of (eta,phi) sectors that each output opto-link takes tracks from. (Depends on MUX scheme).

unsigned int MuxHToutputs::muxFactor() const {
  if (muxOutputsHT_ == 1) {
    return 6;
  } else {
    return numEtaRegions_;
  }
}

//=== Define the MUX algorithm by which tracks from the specified m-bin range in the HT for a given (phi,eta)
//=== sector within a phi octant are multiplexed onto a single output optical link.

unsigned int MuxHToutputs::linkID(unsigned int iSecInOct, unsigned int iEtaReg, unsigned int mBinRange) const {
  unsigned int link;

  // This algorithm multiplexes tracks from different eta sectors onto the a single optical link.

  if (muxOutputsHT_ == 1) {

    //--- This is Mux used for the Dec. 2016 demonstrator.

    // Link 0 contains eta sectors 0, 3, 6, 9, 12 & 15, whilst Link 1 contains eta sectors 1, 4, 7, 10, 13 & 16 etc.
    // The multiplexing is independent of the phi sector or m-bin range, except in that two tracks that
    // differ in either of these quantities will always be sent to different Links.

    if (numEtaRegions_ == 18) {
      link = iEtaReg%3;          // In range 0 to 2
      link += 3*iSecInOct;       // In range 0 to (3*numPhiSecPerOct - 1)
      link += 3*numPhiSecPerOct_ * mBinRange; // In range 0 to (3*numPhiSecsPerOct*numMbinRanges - 1)
    } else {
      throw cms::Exception("MuxHToutputs: MUX algorithm only implemented for 18 eta sectors!");
    }

  } else if (muxOutputsHT_ == 2) {

      //--- This is the Mux for the transverse HT readout organised by m-bin. (Each phi sector & m bin range go to a different link).

      link = 0;          
      link += iSecInOct;     
      link += numPhiSecPerOct_ * mBinRange; 

    } else { 

    throw cms::Exception("MuxHToutputs: Unknown MuxOutputsHT configuration option!");

  }

  if (link >= this->numLinksPerOctant() ) throw cms::Exception("MuxHToutputs: Calculated link ID exceeded expected number of links! ")<<link<<" "<<this->numLinksPerOctant()<<endl;
  return link;
}

//=== Do sanity check of the MUX algorithm implemented in linkID().

void MuxHToutputs::sanityCheck() {
  if ( numPhiSecPerOct_ * numEtaRegions_ % this->muxFactor() != 0) throw cms::Exception("MuxHToutputs: Number of sectors per phi octant is not a multiple of muxFactor().");

  vector<unsigned int> nObsElementsPerLink ( this->numLinksPerOctant(), 0 );
  for (unsigned int iSecInOct = 0; iSecInOct < numPhiSecPerOct_; iSecInOct++) {
    for (unsigned int iEtaReg = 0; iEtaReg < numEtaRegions_; iEtaReg++) {
      for (unsigned int mBinRange = 0; mBinRange < busySectorMbinRanges_.size(); mBinRange++) {
	unsigned int link = this->linkID(iSecInOct, iEtaReg, mBinRange);
	nObsElementsPerLink[ link ] += 1;
      }
    }
  }
  for (const unsigned int& n : nObsElementsPerLink) {
    // Assume good algorithms will distribute sectors & m-bin ranges equally across links.
    if (n != this->muxFactor()) throw cms::Exception("MuxHToutputs: MUX algorithm is not assigning equal numbers of elements per link! ")<<n<<" "<<this->muxFactor()<<endl;
  }
}

}
