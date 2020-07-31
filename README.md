# EXO DiPhoton Analysis Code

## To build. Note that the recipe is different for 2016/2017 and 2018 data/MC.

```
# select cmssw environment
## 2015 (bash syntax)
export SCRAM_ARCH=slc6_amd64_gcc493  
export CMSSW_VERSION=CMSSW_7_6_4 
## 2016 + 2017 (tcsh syntax)
setenv SCRAM_ARCH slc6_amd64_gcc630
setenv CMSSW_VERSION CMSSW_9_4_13
## 2018 (bash syntax)
export SCRAM_ARCH=slc6_amd64_gcc700
export CMSSW_VERSION=CMSSW_10_2_16

# setup cmssw release
cmsrel $CMSSW_VERSION  
cd $CMSSW_BASE/src  
cmsenv  

# EGM recipe from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2
# only for 2018 data/MC
# reverts the changes to EgammaPostRecoTools needed to run in 10_6_X
git cms-merge-topic -u christopheralanwest:revert-41-EgammaPostRecoTools_dev
#git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029 #optional but speeds up the photon ID value module so things fun faster
#now to add the scale and smearing for 2018 (eventually this will not be necessary in later releases but is harmless to do regardless)
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data

# only for 2016/2017 data/MC
# reverts the changes to EgammaPostRecoTools needed to run in 10_6_X
git cms-merge-topic -u christopheralanwest:revert-41-EgammaPostRecoTools_dev
#git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier

# clone repository
git clone git@github.com:cms-exotica-diphotons/diphoton-analysis  

# compile code
cd diphoton-analysis  
scram b -j 8
```
