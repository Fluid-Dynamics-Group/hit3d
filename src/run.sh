# clear out the old data
rm output/*
rm output/velocity/*
rm output/energy/*
rm output/velocity_field/*

# add gitignore files so they will not be overwritten
touch output/.gitignore
touch output/velocity/.gitignore
touch output/velocity_field/.gitignore
touch output/energy/.gitignore

make
nproc=16
mpirun -np "$nproc" ./hit3d.x sample_inp "nosplit"

# run postprocessing stuff
cd ~/github/hit3d-utils/
sh run.sh
