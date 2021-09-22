FilePath="/Users/adamjackson/Documents/summer_project/general_model/"
ExecPath="/Users/adamjackson/Documents/summer_project/general_model/QOMiC_poisson_through_substates_multiv_adaptive/QOMIC/"

simul_seed=1
simulName="2"

inputDates="Original_HMM2014/dates_simul_${simulName}_${simul_seed}.txt"
inputExp="Original_HMM2014/mat_sim_exposure_std.txt"
inputSpecs="Original_HMM2014/details1_simul_${simulName}_${simul_seed}.txt"
inputMortFemale="Original_HMM2014/VBT_mortality_female.txt"
inputMortMale="Original_HMM2014/VBT_mortality_male.txt"

OutFilePath="Results_simulation_${simulName}_${simul_seed}"
mkdir $OutFilePath
seed=1
Name=${simulName}"_"${simul_seed}
NameRun=$Name"_"$seed

output1=$(printf "%s%s" $OutFilePath/  $NameRun)
output2=$(printf "%s%s%s" $OutFilePath/ "log_" $NameRun ".txt")

mydate=$(date)
echo Start time: $mydate

niter=100000
burn_in=10000

$ExecPath"QOMiC" -dates $FilePath$inputDates -exp $FilePath$inputExp -specs $FilePath$inputSpecs \
-mort $FilePath/$inputMortFemale $FilePath/$inputMortMale \
-seed $seed -k 2 -iter $niter -burn_in $burn_in -adaptive 1 -iter_adaptive 100 -update_frequency 100 \
-mu -2 0.01 -lambda1 0.1 0.0004 -lambda2 0.0 0.0 -lambda3 0.0 0.0 -lambda4 0.0 0.0 -gamma 0.5 0.04 \
-out $output1 > $output2

mydate=$(date)
echo End time: $mydate

Rscript $ExecPath"R_monitor.R" -execpath $ExecPath -namefolder $OutFilePath -name $Name -iter $niter -burn_in $burn_in -seed $seed

