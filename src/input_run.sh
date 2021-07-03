# the name of the input file that we are going to run with it
input_file_name=$1
# where the postprocessing data will be exported to
save_location=$2

# start by making the required dirs
mkdir "$save_location"

# clear out the old data from the local output directory
rm output/*
rm output/velocity/*
rm output/energy/*
rm output/velocity_field/*

# add gitignore files so they will not be overwritten
touch output/.gitignore
touch output/velocity/.gitignore
touch output/velocity_field/.gitignore
touch output/energy/.gitignore

# run the solver
nproc=16
mpirun -np "$nproc" ./hit3d.x "$input_file_name" "nosplit" &&\

# copy the input file to the new location
mv "$input_file_name.in" "$save_location/$input_file_name.in" &&\

# run postprocessing stuff
cd ~/github/hit3d-utils/ &&\
sh run.sh "/home/brooks/github/hit3d/src/output/" "$save_location"

cd /home/brooks/github/hit3d/src/
mv output/energy.csv "$save_location"
