import FWCore.ParameterSet.Config as cms

process = cms.Process("exoefficiency")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscms/home/bbbam/mc_RunIIAutumn18MiniAOD_ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8_MINIAODSIM_102X_upgrade2018_realistic_v15-v1_70000_5C1E393E-C20D-E84F-9CD0-6399FA05DE2D.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.exoeff = cms.EDAnalyzer("ExoEfficiencyAnalyzer",
    BeamHaloSummary = cms.InputTag("BeamHaloSummary","","RECO"),
    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
    genInfo = cms.InputTag("generator","","SIM"),
    genparticles = cms.InputTag("prunedGenParticles"),
    photonsMiniAOD = cms.InputTag("slimmedPhotons"),
    rho = cms.InputTag("fixedGridRhoAll"),
    slimmedAddPileupInfo = cms.InputTag("slimmedAddPileupInfo","","PAT"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ExoEfficiency.root')
)


process.p = cms.Path(process.exoeff)


