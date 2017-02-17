#!/bin/tcsh  

# delete all pervious ntuples, so only one CRAB3 date directory (e.g. 161025_234618) exists
# change input_base (and possibly ntuple_path) and out_path accordingly

if ( "$1" != "QCD" && "$1" != "GJets") then
    echo " "
    echo "Set argument 1 as: QCD or GJets"
    echo " "
    exit 1
endif

set sample = $1

echo " "
echo "Merging samples:" $sample
echo " "

set input_base = "/store/user/abuccill/DiPhotonAnalysis/FakeRateClosureTest"

foreach dir (`eosls $input_base | grep $sample`)
    echo $dir
    set ntuple_path = "$input_base/$dir"
    set ntuple_path = "$ntuple_path/`eosls $ntuple_path | grep $sample`"
    set ntuple_path = "$ntuple_path/`eosls $ntuple_path | grep 1`"
    set ntuple_path = "$ntuple_path/`eosls $ntuple_path | grep 00`"
    set out_path = "/store/user/abuccill/diphoton-analysis/fake_rate_closure_test/${sample}_merged_ntuples/diphoton_fake_rate_closure_test_${dir}_76X_MiniAOD_merged.root"
    hadd -f root://cmseos.fnal.gov/$out_path `xrdfsls -u $ntuple_path/ | grep \.root`
    echo " "
end
