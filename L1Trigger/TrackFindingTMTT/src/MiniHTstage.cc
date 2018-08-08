#include "L1Trigger/TrackFindingTMTT/interface/MiniHTstage.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"

#include <L1Trigger/TrackFindingTMTT/interface/Sector.h>

namespace TMTT {

MiniHTstage::MiniHTstage( const Settings* settings ) :
        settings_( settings ),
        miniHTstage_( settings_->miniHTstage() ),
        muxOutputsHT_( settings_->muxOutputsHT() ),
        miniHoughNbinsPt_( settings_->miniHoughNbinsPt() ),
        miniHoughNbinsPhi_( settings_->miniHoughNbinsPhi() ),
        miniHoughMinPt_( settings_->miniHoughMinPt() ),
        miniHoughDontKill_( settings_->miniHoughDontKill() ),
        miniHoughDontKillMinPt_( settings_->miniHoughDontKillMinPt() ),
        numSubSecsEta_( settings_->numSubSecsEta() ),
        numPhiOctants_( settings_->numPhiOctants() ),
        numPhiSecPerOct_( settings_->numPhiSectors() / numPhiOctants_ ),
        numEtaRegions_( settings_->numEtaRegions() ),
        busySectorKill_( settings_->busySectorKill() ),
        busySectorNumStubs_( settings_->busySectorNumStubs()),
        busySectorMbinRanges_( settings_->busySectorMbinRanges() ),
        chosenRofPhi_( settings_->chosenRofPhi() ),
	// Get size of 1st stage HT cells.
        binSizeQoverPtAxis_( miniHoughNbinsPt_ * 2. / (float)settings->houghMinPt() / (float)settings->houghNbinsPt() ),
        binSizePhiTrkAxis_( miniHoughNbinsPhi_ * 2. * M_PI / (float)settings->numPhiSectors() / (float)settings->houghNbinsPhi() ),
        invPtToDphi_( settings_->invPtToDphi() ) {}

void MiniHTstage::exec( matrix< HTrphi >& mHtRphis ) const {

  for ( unsigned int iPhiOct = 0; iPhiOct < numPhiOctants_; iPhiOct++ ) {
    map< unsigned int, unsigned int> numStubsPerLink; // Indices are (link ID, #stubs).
    for ( unsigned int iSecInOct = 0; iSecInOct < numPhiSecPerOct_; iSecInOct++ ) {
      unsigned int iPhiSec = iPhiOct * numPhiSecPerOct_ + iSecInOct;
      for ( unsigned int iEtaReg = 0; iEtaReg < numEtaRegions_; iEtaReg++ ) {

        Sector sector;
        sector.init( settings_, iPhiSec, iEtaReg );
        const float& phiCentre = sector.phiCentre();
        HTrphi& htRphi = mHtRphis( iPhiSec, iEtaReg );
        const vector< L1track2D >& roughTracks = htRphi.trackCands2D();
        vector< L1track2D > fineTracks;

        for (const L1track2D& roughTrk : roughTracks) {

          float roughTrkPhi = reco::deltaPhi( roughTrk.phi0() - chosenRofPhi_ * invPtToDphi_ * roughTrk.qOverPt() - phiCentre, 0. );
          const pair< unsigned int, unsigned int >& cell = roughTrk.getCellLocationHT();
          const vector< const Stub* >& stubs = roughTrk.getStubs();
          bool fineTrksFound = false;
	  bool storeCoarseTrack = false;
          const unsigned int& link = roughTrk.optoLinkID();

	  if (fabs(roughTrk.qOverPt()) < 1./miniHoughMinPt_) { // Not worth using mini-HT at low Pt due to scattering.

	    for ( unsigned int mBin = 0; mBin < miniHoughNbinsPt_; mBin++ ) {
	      float qOverPtBin = roughTrk.qOverPt() - binSizeQoverPtAxis_ / 2. + ( mBin + .5 ) * binSizeQoverPtAxis_ / settings_->miniHoughNbinsPt();
	      for ( unsigned int cBin = 0; cBin < miniHoughNbinsPhi_; cBin++ ){
		float phiBin = reco::deltaPhi( roughTrkPhi - binSizePhiTrkAxis_ / 2. + ( cBin + .5 ) * binSizePhiTrkAxis_ / settings_->miniHoughNbinsPhi(), 0. );
		HTcell htCell;
		const bool mergedCell = false; // This represents mini cell.
		htCell.init( settings_, iPhiSec, iEtaReg, sector.etaMin(), sector.etaMax(), qOverPtBin, cell.first + mBin, mergedCell );

		for ( auto& stub: stubs ) {
		  // Ensure stubs are digitized with respect to the current phi sector.
		  if ( settings_->enableDigitize() )
		    ( const_cast< Stub* >( stub ) )->digitizeForHTinput( iPhiSec );
		  float phiStub = reco::deltaPhi( stub->phi() + invPtToDphi_ * qOverPtBin * ( stub->r() - chosenRofPhi_ ) - phiCentre, 0. );
		  float dPhi = reco::deltaPhi( phiBin - phiStub, 0. );
		  float dPhiMax = binSizePhiTrkAxis_ / miniHoughNbinsPhi_ / 2. + invPtToDphi_ * binSizeQoverPtAxis_ / (float)miniHoughNbinsPt_ * fabs( stub->r() - chosenRofPhi_ ) / 2.;
		  if ( fabs( dPhi ) <= fabs( reco::deltaPhi( dPhiMax, 0. ) ) ) htCell.store( stub, sector.insideEtaSubSecs( stub ) );
		}
		htCell.end();
		if ( htCell.trackCandFound() ) {
		  pair< unsigned int, unsigned int > cellLocation( cell.first + mBin, cell.second + cBin );
		  pair< float, float > helix2D( qOverPtBin, reco::deltaPhi( phiBin + chosenRofPhi_ * invPtToDphi_ * qOverPtBin + phiCentre, 0. ) );
		  L1track2D fineTrk( settings_, htCell.stubs(), cellLocation, helix2D, iPhiSec, iEtaReg, roughTrk.optoLinkID(), mergedCell );
		  // Truncation due to output opto-link bandwidth.
		  bool keep( true );
		  if ( settings_->muxOutputsHT() ) {
		    numStubsPerLink[ link ] += htCell.numStubs();
		    if ( busySectorKill_ && numStubsPerLink[ link ] > busySectorNumStubs_ ) keep = false;
		  }
		  if ( keep ) {
		    fineTracks.push_back( fineTrk );
		    fineTrksFound = true;
		  }
		}
	      }
	    }

	  } else { 

	    // Keep rough track if below Pt threshold where mini-HT in use.
	    storeCoarseTrack = true;
	    
	  }

          if ( storeCoarseTrack || ( (not fineTrksFound) && miniHoughDontKill_ && fabs(roughTrk.qOverPt()) < 1./miniHoughDontKillMinPt_ ) ) {
            bool keep( true );
            if ( settings_->muxOutputsHT() ) {
              numStubsPerLink[ link ] += roughTrk.getNumStubs();
              if ( busySectorKill_ && numStubsPerLink[ link ] > busySectorNumStubs_ )
                keep = false;
            }
            if ( keep ) {
              fineTracks.push_back( roughTrk );
	    }
          }
        }
        // Replace all existing tracks inside HT array with new ones.
        htRphi.replaceTrackCands2D( fineTracks );
      }
    }
  }

}

}
