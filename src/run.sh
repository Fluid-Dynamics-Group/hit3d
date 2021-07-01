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

base_output="/home/brooks/data/hit3d/"
base_save="/home/brooks/sync/hit3d"

# clear the written data from the last runs
# sudo rm -rf "$base_output"
# sudo rm -rf "$base_save"
mkdir "$base_output"
mkdir "$base_save"

#parameters that are constant
order="first"
dt="const"

# parameters that change every line
diffusion="no"
size="64"

# make sure we are working with a fresh executable
make

sh input_run.sh "64_no_diff" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="yes"

sh input_run.sh "64_ye_diff" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

size="128"
diffusion="no"

sh input_run.sh "128no_diff" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="yes"

sh input_run.sh "128ye_diff" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="no"
dt="const-0.00025"

sh input_run.sh "128shortdt" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

dt="const-0.0005"

sh input_run.sh "128long_dt" "$base_output/$size-$dt-dt-$diffusion-diffusion-$order-order" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
