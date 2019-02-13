#ifndef __HLSpragmaOpts__
#define __HLSpragmaOpts__

/**
 * Define options specified by pragma statements.
 *
 * Author: Ian Tomalin
 */

// OPTION 1:
// If defined, use 25 bits for the off-diagonal elements of the helix param covariance matrix.
// Ian thinks this is needed to avoid bit overflows messing up some tracks (although the effect on tracking
// performance is small). And it would require changing the KF interface in the firmware ...
// It is unnecessary if the hit errors are inflated to allow for scattering.

//#define COV_EXTRA_BITS

// OPTION 2:
// If defined, set optimum numbers of bits for Ultrascale instead of Virtex7 FPGAs. 

#define ULTRASCALE

// OPTION 3:
// If defined, HLS KF will cope with tracking down to 2 GeV Pt instead of 3 GeV.
//#define PT_2GEV

#endif
