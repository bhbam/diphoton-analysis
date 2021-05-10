import FWCore.ParameterSet.Config as cms

# from os.path import basename
# import os
# import sys
# import importlib
# submit_utils = importlib.import_module("diphoton-analysis.CommonClasses.submit_utils")

# JEC = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
# process.load("FWCore.MessageService.MessageLogger_cfi")

process = cms.Process("exoefficiency")
# process.load("FWCore.MessageService.MessageLogger_cfi")
# process.MessageLogger.cerr.FwkReport.reportEvery = 100
# process.MessageLogger.suppressWarning.append('ExoEfficiencyAnalyzer')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#
# process.load("Configuration.StandardSequences.GeometryDB_cff")

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/uscms/home/bbbam/mc_RunIIAutumn18MiniAOD_ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8_MINIAODSIM_102X_upgrade2018_realistic_v15-v1_70000_5C1E393E-C20D-E84F-9CD0-6399FA05DE2D.root'

    )
)


process.TFileService = cms.Service("TFileService",
                fileName = cms.string("ExoEfficiency.root")
                            )

process.exoeff = cms.EDAnalyzer('ExoEfficiencyAnalyzer',
    genparticles = cms.InputTag("prunedGenParticles"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    slimmedAddPileupInfo = cms.InputTag("slimmedAddPileupInfo", "", "PAT"),
    BeamHaloSummary = cms.InputTag("BeamHaloSummary", "", "RECO"),
    photonsMiniAOD = cms.InputTag("slimmedPhotons"),
    rho = cms.InputTag("fixedGridRhoAll"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO")
    #minPhotonPt = cms.double(125.),
)


process.p = cms.Path(process.exoeff)
