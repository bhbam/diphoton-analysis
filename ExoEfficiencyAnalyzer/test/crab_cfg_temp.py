#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = '%s_%s_%s'%(CFG,PU,EVTCONT)

dataset = 'DATASETNAME'
taskname = dataset[1:].replace('/','__')
if(len(taskname)>100): taskname = taskname[0:99]
config.General.requestName = taskname



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

config.Data.inputDataset = 'DATASETNAME'



# config.Data.totalUnits = 500
# config.Data.publication = False
# config.Data.outputDatasetTag = '%s_%s_%s'
# config.Data.outputPrimaryDataset = 'ExoEfficiency'
config.Data.outLFNDirBase = '/store/user/bhbam/ExoEfficiencyBackground018'  # put cern user name instead of fnal in order work proprely
# config.Data.outLFNDirBase = '/store/user/bhbam/16_Exoeff'  # put cern user name instead of fnal in order work proprely
config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
# config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T1_RU_JINR', 'T2_US_Vanderbilt']

#crab submit -c crabConfig_MC.py
