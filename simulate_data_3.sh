FilePath="/Users/adamjackson/Documents/summer_project/general_model/"
ExecPath="/Users/adamjackson/Documents/summer_project/general_model/QOMiC_poisson_through_substates_multiv_adaptive/QOMiC_simul/"

simul_seed=1
simulName="2"

inputDates="Original_HMM2014/dates_simul_${simulName}_${simul_seed}.txt"
inputExp="Original_HMM2014/mat_sim_exposure_std.txt"
inputSpecs="Original_HMM2014/details1_simul_${simulName}_${simul_seed}.txt"
inputMortFemale="Original_HMM2014/VBT_mortality_female.txt"
inputMortMale="Original_HMM2014/VBT_mortality_male.txt"
inputMCMC="Results_simulation_2_1/2_1_1_History_iter_100000_k_2.txt"

OutFilePath="./Results"
mkdir $OutFilePath
seed=1
Name=${simulName}"_"${simul_seed}
NameRun=$Name"_"$seed

scenario="simulation_"$simulName"_"$simul_seed

mydate=$(date)
echo Start time: $mydate

niter=10
burn_in=60000
sweeps=100000

$ExecPath"QOMiC_simul" -dates $FilePath$inputDates -exp $FilePath$inputExp -specs $FilePath$inputSpecs \
-mort $FilePath/$inputMortFemale $FilePath/$inputMortMale \
-seed $seed -k 2 -iter $niter -burn_in $burn_in \
-MCMC $inputMCMC -scenario $scenario -sweeps $sweeps \
> ./Results/log.txt

mydate=$(date)
echo End time: $mydate
