import os
DATASETS = [[]]

## 2018 MC sets

# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])


## 2017 MC sets

# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'])

## 2016 MC sets

# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-500To1000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-6000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-500To1000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-8000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-500To1000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-10000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-500To1000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-1000To2000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])
# DATASETS.append(['/ADDGravToGG_NegInt-0_LambdaT-13000_M-2000To4000_TuneCP2_13TeV-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'])

## 2018 Background MC DATASETS
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_50_120/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_120_200/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_200_400/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_400_800/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_800_1400/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])
DATASETS.append(['/ZToEE_NNPDF30_13TeV-powheg_M_3500_4500/RunIIAutumn18MiniAOD-NZSPU0to70_102X_upgrade2018_realistic_v15-v1/MINIAODSIM'])


for ilist in DATASETS:

  for ids in ilist:
    cmssw_base = os.getenv("CMSSW_BASE")
    datasetID = ids.replace('/', '', 1).replace('/', '_', 1)
    datasetID = datasetID[0:datasetID.find('/')]
    inputfile = cmssw_base + "/src/diphoton-analysis/ExoEfficiencyAnalyzer/test/crab_cfg_temp.py"
    outputfile = "crab_cfg_" + datasetID + ".py"

    s = open(inputfile).read()
    s = s.replace('DATASETNAME', ids)
    # s = s.replace('NEVENTS', str(nevents))
    f = open(outputfile, 'w')
    f.write(s)
    f.close()
    print "Wrote crab configuration file " + outputfile

    cmd = "crab submit -c " + outputfile
    os.system(cmd)
    print "Submitted ", ids
