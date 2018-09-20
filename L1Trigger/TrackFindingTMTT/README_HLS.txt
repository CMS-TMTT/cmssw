To run the TMTT tracking chain in CMSSW with the HLS version of the KF, you must have 
Vivado on your computer, and have set it up by doing something like:

setenv VIVADO_DIR /opt/ppd/tools/xilinx/Vivado/2018.1
source $VIVADO_DIR/settings.csh

Then checkout the TMTT CMSSW software following instructions in https://github.com/CMS-TMTT/cmssw/blob/TMTT_1025/L1Trigger/TrackFindingTMTT/README.md .

You should then:

1) cd L1Trigger/TrackFindingTMTT/
1) mv  BuildFile.xml      BuildFile_original.xml
2) mv  BuildFile_HLS.xml  BuildFile.xml
3) Edit hls.xml , changing variable HLS_BASE to point to your local Vivado directory.
4) scram setup hls.xml
5) Set the TMTT cfg parameter:
    TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring("KF4ParamsCombHLS")


Notes:

a) The HLS code is in src/HLS/ & interface/HLS/

b) BuildFile_HLS.xml defines a pragma variable USE_HLS , which is used to switch on a few lines of 
C++ in TrackFitGeneric.cc that call the HLS code.

c) In addition to the usual track performance summary, the HLS code prints a summary recording 
if any finite bit HLS variables overflowed. This would print lots of warnings because we use slightly 
too few bits for the off-diagonal elements of the helix covariance matrix (forced on us by Andy's KF
VHDL). These warnings do not results in any significant loss of tracking performance, so are unimportant.
So to silence them, it is convenient to add extra bits to these off-diagonal elements, which 
BuildFile_HLS.xml does by defining pragma COV_EXTRA_BITS. You can disable this line to check that they
aren't needed to improve the tracking performance.
