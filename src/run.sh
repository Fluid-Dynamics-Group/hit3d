base_save="/home/brooks/sync/hit3d"

# clear the written data from the last runs
mkdir "$base_save"

#parameters that are constant
order="first"
dt="const-0.001"

# parameters that change every line
diffusion="no"
size="64"

# make sure we are working with a fresh executable
make

# generate the config files that we are going to run
hit3d-config --n 64 --steps 20000 --steps-between-io 100 --flow-type 0 --skip-diffusion 1 --dt -0.001 64_no_diff.in --restarts 0&&\
hit3d-config --n 64 --steps 20000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 --dt -0.001 64_ye_diff.in --restarts 0&&\

#hit3d-config --n 128 --steps 20000 --steps-between-io 100 --flow-type 0 --skip-diffusion 1 --dt -0.001 128no_diff.in &&\
#hit3d-config --n 128 --steps 20000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 --dt -0.001 128ye_diff.in &&\
#
## todo - this is basically a repitition of the previous one
#hit3d-config --n 128 --steps 40000 --steps-between-io 100 --flow-type 0 --skip-diffusion 0 --dt -0.0005 128long_dt.in

sh input_run.sh "64_no_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
echo "finished 1/n"

diffusion="yes"

sh input_run.sh "64_ye_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"

#sh input_run.sh "64_ye_diff" "$base_save/mgm_compare"
#echo "finished 2/n"
#
#size="128"
#diffusion="no"
#
#sh input_run.sh "128no_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
#echo "finished 3/n"
#
#diffusion="yes"
#
#sh input_run.sh "128ye_diff" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
#echo "finished 4/n"
#
#dt="const-0.0005"
#
#sh input_run.sh "128long_dt" "$base_save/$size-$dt-dt-$diffusion-diffusion-$order-order"
#echo "finished 5/n"
