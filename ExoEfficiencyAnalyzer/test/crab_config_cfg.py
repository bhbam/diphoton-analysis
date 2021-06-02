#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

# config.General.requestName = 'TESTFirst-'
# config.General.requestName = 'TESTFirstANa'
# config.General.requestName = 'LambdaT_13000'
# config.General.requestName = 'LambdaT_13000_1'
# config.General.requestName = 'LambdaT_13000_2'
# config.General.requestName = 'LambdaT_8000'
# config.General.requestName = 'LambdaT_8000_1'
# config.General.requestName = 'LambdaT_8000_2'
# config.General.requestName = 'LambdaT_6000'
# config.General.requestName = 'LambdaT_6000_1'
# config.General.requestName = 'LambdaT_6000_2'


# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000'
# config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000'
config.General.requestName = '2017_ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000'


config.General.workArea = 'outcrab'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/uscms/home/bbbam/nobackup/CMSSW_10_2_26/src/diphoton-analysis/ExoEfficiencyAnalyzer/python/confFile_eff_cfg.py'
# config.JobType.maxMemoryMB = 2800
config.JobType.allowUndistributedCMSSW = True
config.section_("Data")
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
# config.Data.splitting = 'Automatic'
# config.Data.unitsPerJob = 200

## 2018 dataset

# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'


## 2017 data sets

# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
# config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDataset = '/ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'




# config.Data.totalUnits = 500
# config.Data.publication = False
# config.Data.outputDatasetTag = '%s_%s_%s'

config.Data.outLFNDirBase = '/store/user/bhbam/ExoEfficiency017'  # put cern user name instead of fnal in order work proprely
config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
# config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T1_RU_JINR', 'T2_US_Vanderbilt']

#crab submit -c crabConfig_MC.py