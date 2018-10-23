#ifndef __KALMAN_STATE__
#define __KALMAN_STATE__
 
#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1KalmanComb.h"
#include <map>

namespace TMTT {

class L1KalmanComb;
class kalmanState;
class StubCluster;

typedef std::map<std::string, double> (*GET_TRACK_PARAMS)( const L1KalmanComb *p, const kalmanState *state );
 
class kalmanState{
    public:
	kalmanState();
	kalmanState( const L1track3D& candidate, unsigned n_skipped, unsigned kLayer_next, unsigned layerId, const kalmanState *last_state, 
		const std::vector<double> &x, const TMatrixD &pxx, const TMatrixD &K, const TMatrixD &dcov, 
		const StubCluster* stubcl, double chi2, 
		L1KalmanComb *fitter, GET_TRACK_PARAMS f );
	kalmanState(const kalmanState &p);
	~kalmanState(){}

	kalmanState & operator=( const kalmanState &other );

	unsigned             nextLayer()const{ return      kLayerNext_; }
	unsigned               layerId()const{ return         layerId_; }
	unsigned            endcapRing()const{ return      endcapRing_; }
	bool                    barrel()const{ return          barrel_; }
	unsigned        nSkippedLayers()const{ return       n_skipped_; }
        // Hit coordinates.
	double                       r()const{ return               r_; }
	double                       z()const{ return               z_; }
	const kalmanState  *last_state()const{ return      last_state_; }
        // Helix parameters (1/2R, phi relative to sector, z0, tanLambda) 
	std::vector<double>         xa()const{ return              xa_; }
        // Covariance matrix on helix params.
	TMatrixD                  pxxa()const{ return            pxxa_; }
        // Kalman Gain matrix 
	TMatrixD                     K()const{ return               K_; }
        // Hit position covariance matrix.
	TMatrixD                  dcov()const{ return            dcov_; }
        // Hit
	const StubCluster* stubCluster()const{ return     stubCluster_; }
	double                    chi2()const{ return            chi2_; }
	unsigned           nStubLayers()const{ return         n_stubs_; }
        L1track3D            candidate()const{ return       l1track3D_; }

	bool                            good( const TP *tp )const;
	double                   reducedChi2()const;
	const kalmanState *last_update_state()const;
	std::vector<const Stub *>      stubs()const;
	L1KalmanComb                 *fitter()const{ return fitter_; }
	GET_TRACK_PARAMS     fXtoTrackParams()const{ return fXtoTrackParams_; };


	static bool orderReducedChi2(const kalmanState *left, const kalmanState *right);
	static bool orderChi2(const kalmanState *left, const kalmanState *right);
	static bool orderMinSkipChi2(const kalmanState *left, const kalmanState *right);

	static bool order(const kalmanState *left, const kalmanState *right);
	void dump( ostream &os, const TP *tp=0, bool all=0 )const;
	void setChi2( double p ){ chi2_ = p; }

        // If using HLS, note/get additional output produced by HLS core.
        void setHLSextra(unsigned int mBinHelix, unsigned int cBinHelix, bool consistent) { mBinHelixHLS_ = mBinHelix; cBinHelixHLS_ = cBinHelix; consistentHLS_ = consistent;}
        void getHLSextra(unsigned int& mBinHelix, unsigned int& cBinHelix, bool& consistent) const { mBinHelix = mBinHelixHLS_; cBinHelix = cBinHelixHLS_; consistent = consistentHLS_;}

    private:
	unsigned              kLayerNext_;
	unsigned                 layerId_;
	unsigned              endcapRing_;
	double                         r_;
	const kalmanState    *last_state_;
	std::vector<double>           xa_;
	TMatrixD                    pxxa_;
	TMatrixD                       K_;
	TMatrixD                    dcov_;
	const StubCluster   *stubCluster_;
	double                      chi2_;
	unsigned                 n_stubs_;
	L1KalmanComb             *fitter_;
	GET_TRACK_PARAMS fXtoTrackParams_;
	bool                      barrel_;
	unsigned               n_skipped_;
	double                         z_;
        L1track3D              l1track3D_;

       // Additional output from HLS if using it.
       unsigned int mBinHelixHLS_; 
       unsigned int cBinHelixHLS_; 
       bool consistentHLS_; 
};

}

#endif


