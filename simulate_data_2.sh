FilePath="/Users/adamjackson/Documents/summer_project/general_model/"
ExecPath="/Users/adamjackson/Documents/summer_project/general_model/QOMiC_poisson_through_substates_multiv_adaptive/Data_simulation/"

inputDates="Original_HMM2014/dates_pooled.txt"
inputExp="Original_HMM2014/mat_sim_exposure_std.txt"
inputSpecs="Original_HMM2014/misc_details_pooled.txt"
inputMortFemale="Original_HMM2014/VBT_mortality_female.txt"
inputMortMale="Original_HMM2014/VBT_mortality_male.txt"

seed=1
myname="2"
output1="simulation_"$myname"_"$seed

mkdir Simulations

mydate=$(date)
echo Start time: $mydate

$ExecPath"Data_simul" -dates $FilePath$inputDates -exp $FilePath$inputExp -specs $FilePath$inputSpecs \
-mort $FilePath/$inputMortFemale $FilePath/$inputMortMale \
-k 2 -seed $seed -scenario $output1 \
-mu -3 -lambda1 0.05 -lambda2 0.0 -lambda3 0.0 -lambda4 0.0 -gamma 1 > ./Simulations/log_$output1.txt

Rscript $ExecPath"create_input_files.R" -dates $FilePath$inputDates -specs $FilePath$inputSpecs \
-execpath $ExecPath \
-trajectories "./Simulations/Trajectories_"$output1".txt" -output $myname \
-filepath $FilePath"Original_HMM2014/dates_pooled.txt" -logpath "./Simulations" -seed $seed

mydate=$(date)
echo End time: $mydate

