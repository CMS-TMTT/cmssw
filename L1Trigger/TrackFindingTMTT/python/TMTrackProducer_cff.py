import FWCore.ParameterSet.Config as cms

#---------------------------------------------------------------------------------------------------------
# This describes the full TMTT track reconstruction chain with 3 GeV threshold, where:
# the GP divides the tracker into 18 eta sectors (each sub-divided into 2 virtual eta subsectors);  
# the HT uses a  32x18 array followed by 2x2 Mini-HT array, with transverese HT readout & multiplexing, 
# followed by the KF (or optionally SF+SLR) track fit; duplicate track removal (Algo50) is run.
#
# This usually corresponds to the current firmware.
#---------------------------------------------------------------------------------------------------------

#=== Import default values for all parameters & define EDProducer.

from L1Trigger.TrackFindingTMTT.TMTrackProducer_Defaults_cfi import TMTrackProducer_params

TMTrackProducer = cms.EDProducer('TMTrackProducer',
  # Load cfg parameters from TMTrackProducer_Defaults_cfi.py
  TMTrackProducer_params
)

#===================================================================================================
#=== Parameters changed from their default values.
#===================================================================================================

#--- Disable internal digitisation of SimpleLR fitter, as it was never retuned for nonants.
TMTrackProducer.TrackFitSettings.DigitizeSLR = cms.bool(False)

#===================================================================================================
#=== All the following parameters already have identical values in TMTrackProducer_Defaults_cfi .
#=== They are listed here just to remind you of the most interesting parameters to play with.
#===================================================================================================

#--- Configure track fitting

# Use only 4 or 5 parameter helix fit Kalman Filter (which automatically runs on tracks produced with no r-z track filter)
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring("KF5ParamsComb", "KF4ParamsComb")
# Use only Linear Regression Fitter (which automatically runs on tracks produced by r-z track filter).
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring("SimpleLR")

# Allow KF to assign stubs in up to this many layers to fitted tracks.
#TMTrackProducer.TrackFitSettings.KalmanMaxNumStubs  = cms.uint32(6)
# Enable more sophisticated fit mathematics in KF.
#TMTrackProducer.TrackFitSettings.KalmanHOtilted     = cms.bool(True)
#TMTrackProducer.TrackFitSettings.KalmanHOhelixExp   = cms.bool(True)
#TMTrackProducer.TrackFitSettings.KalmanHOalpha      = cms.uint32(2)
#TMTrackProducer.TrackFitSettings.KalmanHOprojZcorr  = cms.uint32(2)
#TMTrackProducer.TrackFitSettings.KalmanHOdodgy      = cms.bool(False)

#--- Switch off parts of the track reconstruction chain.

#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi   = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlg3D     = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit    = cms.uint32(0)
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()

#--- Keep Pt threshold at 3 GeV, with coarse HT, but switch off Mini-HT.

#TMTrackProducer.HTArraySpecRphi.MiniHTstage         = cms.bool(False)  
#TMTrackProducer.HTFillingRphi.MaxStubsInCell        = cms.uint32(16) 
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt        = cms.uint32(18)  # Mini cells in whole HT array.
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi       = cms.uint32(32) 
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges  = cms.vuint32(2,2,2,2,2,2,2,2,2)   
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder   = cms.vuint32(0,9, 1,10, 2,11, 3,12, 4,13, 5,14, 6,15, 7,16, 8,17)

#--- Reduce Pt threshold to 2 GeV, with coarse HT, and switch off Mini-HT.

#TMTrackProducer.HTArraySpecRphi.MiniHTstage        = cms.bool(False)  
#TMTrackProducer.HTFillingRphi.MaxStubsInCell       = cms.uint32(16) 
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt       = cms.uint32(27)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(32) 
#TMTrackProducer.GenCuts.GenMinPt                   = cms.double(2.0)
#TMTrackProducer.HTArraySpecRphi.HoughMinPt         = cms.double(2.0)
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(2,2,2,2,2,2,2,2,2,2,2,2,2,1)   
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder  = cms.vuint32(0,14, 1,15, 2,16, 3,17, 4,18, 5,19, 6,20, 7,21, 8,22, 9,23, 10,24, 11,25, 12,26, 13)

#--- Reduce Pt threshold to 2 GeV, with coarse HT, followed  by Mini-HT.

#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt        = cms.uint32(54)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi       = cms.uint32(64) 
#TMTrackProducer.GenCuts.GenMinPt                    = cms.double(2.0)
#TMTrackProducer.HTArraySpecRphi.HoughMinPt          = cms.double(2.0)
#TMTrackProducer.HTArraySpecRphi.MiniHoughMinPt      = cms.double(3.0) # Mini-HT not used below this Pt, to reduce sensitivity to scattering.
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges  = cms.vuint32(2,2,2,2,2,2,2,2,2,2,2,2,2,1, 27)   
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder   = cms.vuint32(0,28, 2,30, 4,32, 6,34, 8,36, 10,38, 12,40, 14,42, 16,44, 18,46, 20,48, 22,50, 24,52, 26, 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53)

#--- Additional Mini-HT options to improve electron/displaced tracking.

# Next 2 lines cause tracks found by 1st stage HT to be output if above specified Pt threshold & mini-HT found no tracks.
# Improves electron tracking. Setting Pt threshold to 0 improves displaced tracking.
#TMTrackProducer.HTArraySpecRphi.MiniHoughDontKill   = cms.bool(True)  
#TMTrackProducer.HTArraySpecRphi.MiniHoughDontKillMinPt   = cms.double(8.)  
# Extreme displaced tracking also benefits from following.
#TMTrackProducer.L1TrackDef.MinStubLayers            = cms.uint32(4)  # HT accepts tracks with >= 4 layers
#TMTrackProducer.TrackFitSettings.KalmanRemove2PScut = cms.bool(True)
#To study displaced tracking, include non-prompt particles in efficiency definition.
#TMTrackProducer.GenCuts.GenMaxVertR                 = cms.uint32(30.) 

#--- Unusual HT cell shapes

# Simplify HT MUX to allow easy playing with the number of m bins.
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder  = cms.vuint32() 

# Diamond shaped cells: (64,62), (34,32) or (46,44) sized array interesting.
#TMTrackProducer.HTArraySpecRphi.Shape         = cms.uint32(1)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt  = cms.uint32(38)  
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(32)  

# Hexagonal shaped cells: (64,42), (50,32) or (56,36) sized array interesting.
#TMTrackProducer.HTArraySpecRphi.Shape         = cms.uint32(2)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt  = cms.uint32(56)  
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(32)  

# Brick-wall arranged cells: (64,30) or (66,32) sized array interesting.
#TMTrackProducer.HTArraySpecRphi.Shape         = cms.uint32(3)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt  = cms.uint32(64)  
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(27)  

#--- Older cfg giving similar tracking performance with slightly larger resource use.

#TMTrackProducer.PhiSectors.NumPhiSectors           = cms.uint32(36)
#TMTrackProducer.EtaSectors.EtaRegions = cms.vdouble(-2.4, -2.0, -1.53, -0.98, -0.37, 0.37, 0.98, 1.53, 2.0, 2.4)
#TMTrackProducer.EtaSectors.ChosenRofZ              = cms.double(45.)     
#TMTrackProducer.EtaSectors.AllowOver2EtaSecs       = cms.bool(False)
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi      = cms.uint32(32)
#TMTrackProducer.HTArraySpecRphi.NumSubSecsEta      = cms.uint32(1)

#--- Stub digitization (switch on/off and/or change defaults).

#TMTrackProducer.StubDigitize.EnableDigitize  = cms.bool(True)

#--- Reduce requirement on number of layers a track must have stubs in, either globally or in specific eta regions.

#TMTrackProducer.L1TrackDef.MinStubLayers       = cms.uint32(4)  # Reduce it globally
#TMTrackProducer.L1TrackDef.EtaSecsReduceLayers = cms.vuint32(5,12) # barrel-endcap transition region

#--- If globally reducing number of layers cut, best to also use just one HT output opto-link per m-bin.
# For 3 GeV threshold with no mini-HT.
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)   
# For 2 GeV threshold with mini-HT.
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 27)   

#--- Change TP to track matching criteria.

#TMTrackProducer.GenCuts.GenMinStubLayers               = cms.uint32(4)
#TMTrackProducer.TrackMatchDef.MinNumMatchLayers        = cms.uint32(4)

#--- Switch off data truncation due to finite band-width.

#TMTrackProducer.HTFillingRphi.BusySectorKill       = cms.bool(False)
#TMTrackProducer.HTFillingRphi.BusyInputSectorKill  = cms.bool(False)

# Don't order stubs by bend in DTC, such that highest Pt stubs are transmitted first.
#TMTrackProducer.StubCuts.OrderStubsByBend = cms.bool(False)

#--- Switch on FPGA-friendly approximation to B parameter in GP - will be used in future GP firmware.
#--- (used to relate track angle dphi to stub bend) 
#TMTrackProducer.GeometricProc.UseApproxB           = cms.bool(True)

#--- Use octants instead of nonants. 

#TMTrackProducer.PhiSectors.NumPhiOctants      = cms.uint32(8)   
#TMTrackProducer.PhiSectors.NumPhiSectors      = cms.uint32(16)   
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt  = cms.uint32(32)   
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(16)  
#TMTrackProducer.HTFillingRphi.BusySectorMbinRanges = cms.vuint32(2,2,2,2,2,2,2,2)
#TMTrackProducer.HTFillingRphi.BusySectorMbinOrder  = cms.vuint32(0,8, 1,9, 2,10, 3,11, 4,12, 5,13, 6,14, 7,15)
#
#TMTrackProducer.StubDigitize = cms.PSet(
#     EnableDigitize  = cms.bool(True),  # Digitize stub coords? If not, use floating point coords.
#     FirmwareType    = cms.uint32(1),    # 0 = Old Thomas 2-cbin data format, 1 = new Thomas data format used for daisy chain, 2-98 = reserved for demonstrator daisy chain use, 99 = Systolic array data format.
#     #
#     #--- Parameters available in MP board.
#     #
#     PhiSectorBits   = cms.uint32(6),    # Bits used to store phi sector number
#     PhiSBits        = cms.uint32(14),   # Bits used to store phiS coord. (13 enough?)
#     PhiSRange       = cms.double(0.78539816340),  # Range phiS coord. covers in radians.
#     RtBits          = cms.uint32(12),   # Bits used to store Rt coord.
#     RtRange         = cms.double(103.0382), # Range Rt coord. covers in units of cm.
#     ZBits           = cms.uint32(14),   # Bits used to store z coord.
#     ZRange          = cms.double(640.), # Range z coord. covers in units of cm.
#     # The following four parameters do not need to be specified if FirmwareType = 1 (i.e., daisy-chain firmware) 
#     DPhiBits        = cms.untracked.uint32(8),    # Bits used to store Delta(phi) track angle.
#     DPhiRange       = cms.untracked.double(1.),   # Range Delta(phi) covers in radians.
#     RhoBits         = cms.untracked.uint32(6),    # Bits used to store rho parameter.
#     RhoRange        = cms.untracked.double(0.25), # Range rho parameter covers.
#     #
#     #--- Parameters available in GP board (excluding any in common with MP specified above).
#     #
#     PhiOBits        = cms.uint32(15),      # Bits used to store PhiO parameter.
#     PhiORange       = cms.double(1.5707963268), # Range PhiO parameter covers.
#     BendBits        = cms.uint32(6)        # Bits used to store stub bend.
#)

#TMTrackProducer.TrackDigi=cms.PSet(
#    # For firmware reasons, can't use common digitisation cfg for all fitters.
#
#    #======= SimpleLR digi parameters ========
#    SLR_skipTrackDigi = cms.bool( False ), # Optionally skip track digitisation if done internally inside fitting code.
#    SLR_oneOver2rBits = cms.uint32(13),
#    SLR_oneOver2rRange = cms.double(0.0076223979397932),
#    SLR_d0Bits = cms.uint32(12), # Made up by Ian as never yet discussed.
#    SLR_d0Range  = cms.double(10.),
#    SLR_phi0Bits = cms.uint32(18),
#    SLR_phi0Range = cms.double(0.78539816340), # phi0 is actually only digitised relative to centre of sector.
#    SLR_z0Bits = cms.uint32(12),
#    SLR_z0Range  = cms.double(160),
#    SLR_tanlambdaBits = cms.uint32(15),
#    SLR_tanlambdaRange = cms.double(49.6903090310196),
#    SLR_chisquaredBits = cms.uint32(8),
#    SLR_chisquaredRange = cms.double(128.),
#    
#    #====== Kalman Filter Digi parameters ========
#    KF_skipTrackDigi = cms.bool( False ), # Optionally skip track digitisation if done internally inside fitting code.
#    KF_oneOver2rBits = cms.uint32(18),
#    KF_oneOver2rRange = cms.double(0.06097882386778998),
#    KF_d0Bits = cms.uint32(12), # Made up by Ian as never yet discussed.
#    KF_d0Range  = cms.double(10.),
#    KF_phi0Bits = cms.uint32(18),
#    KF_phi0Range = cms.double(0.7855158),  # phi0 is actually only digitised relative to centre of sector.
#    KF_z0Bits = cms.uint32(18),
#    KF_z0Range  = cms.double(51.5194204),
#    KF_tanlambdaBits = cms.uint32(18),
#    KF_tanlambdaRange = cms.double(32.),
#    KF_chisquaredBits = cms.uint32(17),
#    KF_chisquaredRange = cms.double(1024.),
#
#    #====== Other track fitter Digi params.
#    # Currently equal to those for KF, although you can skip track digitisation for them with following.
#    Other_skipTrackDigi = cms.bool( True ) 
#)

