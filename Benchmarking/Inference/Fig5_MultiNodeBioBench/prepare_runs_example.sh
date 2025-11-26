

run_dir="MediumModel_lmax-3"
#/ must be escaped:
model_path="\/path\/to\/Model\/"
#. must be escaped:
declare -a model_array=("MediumModel_lmax-3_kernel\.nequip\.pt2" "MediumModel_lmax-3\.nequip\.pt2" "MediumModel_lmax-3\.nequip\.pth")
declare -a seed_array=("401")
#Compute information
declare -a numnodes_array=("1" "2" "4" "8" "16" "32" "64" "128" "256" "512" "1024")
declare -a numtaskspernode_array=("8")
partition="batch"
gpu_name="MI250X" #A100 or H100

#Optional, additional inputs to the lammps scripts, Could have them either indexed together or in a nested loop, depending on later on treatment in the script.

#declare -a l_array=("1" "2" "3" "4" "5" "6" "7" "8")


mkdir $run_dir
cd $run_dir

for model in "${model_array[@]}" ; do
for numnodes in  "${numnodes_array[@]}" ; do
for numtaskspernode in "${numtaskspernode_array[@]}"; do
for seed in "${seed_array[@]}"; do


target=model-${model//\\/}_N-${numnodes}_n-${numtaskspernode}_seed-${seed}


cp -r ../Template_DHFR ./${target//\'/}
#cp -r ../Template_Cellulose ./${target//\'/}

#submit script variables
sed -i -e "s/NUMNODESBLANK/${numnodes}/g" ${target}/submit_md.job
sed -i -e "s/NUMTASKSPERNODEBLANK/${numtaskspernode}/g" ${target}/submit_md.job
#sed -i -e "s/GPUBLANK/${gpu}/g" ${target}/submit_md.job
sed -i -e "s/GPUNAMEBLANK/${gpu_name}/g" ${target}/submit_md.job
sed -i -e "s/PARTITIONBLANK/${partition}/g" ${target}/submit_md.job

#lammps updates
sed -i -e "s/MODELBLANK/${model_path}\/${model}/g" ${target}/lammps.in
sed -i -e "s/SEEDBLANK/${seed}/g" ${target}/lammps.in


done; done; done; done;
