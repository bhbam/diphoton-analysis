import FWCore.ParameterSet.Config as cms



process = cms.Process("exoefficiency")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    ##for Photons----------------------------------------------------------------------------------------
    # fileNames = cms.untracked.vstring(
    #     'file:/uscms/home/bbbam/mc_RunIIAutumn18MiniAOD_ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8_MINIAODSIM_102X_upgrade2018_realistic_v15-v1_70000_5C1E393E-C20D-E84F-9CD0-6399FA05DE2D.root'
    ##for Electrons------------------------------------------------------------------------------------
    fileNames = cms.untracked.vstring(
        'file:/uscms/home/bbbam/signal_EB_pass.root'
    ##---------------------------------------------------------------------------------------------------
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
    # photonsMiniAOD = cms.InputTag("slimmedElectrons"),
    rho = cms.InputTag("fixedGridRhoAll"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    electronSource = cms.InputTag("slimmedElectrons", "", "PAT"),
    conversionSource = cms.InputTag( "reducedEgamma","reducedConversions","PAT" )


)


process.p = cms.Path(process.exoeff)
