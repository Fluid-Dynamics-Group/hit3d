base_save="/home/brooks/sync/hit3d"

# clear the written data from the last runs
mkdir "$base_save"

#parameters that are constant
order="first"
dt="const"

# parameters that change every line
diffusion="no"
size="64"

# make sure we are working with a fresh executable
make

# generate the config files that we are going to run
hit3d-config --n 64 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 1 64_no_diff.in &&\
hit3d-config --n 64 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 64_ye_diff.in &&\

hit3d-config --n 128 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 1 128no_diff.in &&\
hit3d-config --n 128 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 128ye_diff.in &&\

# todo - this is basically a repitition of the previous one
hit3d-config --n 128 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 --dt -0.0005 128shortdt.in &&\
hit3d-config --n 128 --steps 80000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 --dt -0.00025 128long_dt.in

sh input_run.sh "64_no_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="yes"

sh input_run.sh "64_ye_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

size="128"
diffusion="no"

sh input_run.sh "128no_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="yes"

sh input_run.sh "128ye_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

diffusion="no"
dt="const-0.00025"

sh input_run.sh "128shortdt" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

dt="const-0.0005"

sh input_run.sh "128long_dt" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
