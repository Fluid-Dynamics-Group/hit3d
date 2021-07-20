import subprocess
import shutil
import os
import traceback

BASE_SAVE = "/home/brooks/sync/hit3d"
MPI_PROC = 1

class RunCase():
    def __init__(self,skip_diffusion, size, dt, steps, restarts, reynolds_number,path):
        # if there are restarts, find the number of steps spent in that those restarts
        # and add them to the current number of steps
        simulation_time_restart = restarts * 1.0
        steps_restarts = int(simulation_time_restart / dt)
        steps += steps_restarts

        self.skip_diffusion =  skip_diffusion
        self.size =            size
        self.dt =              dt
        self.steps =           steps
        self.restarts =        restarts
        self.reynolds_number = reynolds_number
        self.path =            path

    def run(self, iteration):
        if self.dt == 0.001:
            io_steps = 100
        elif self.dt == 0.0005:
            io_steps = 200
        elif self.dt == 0.00025:
            io_steps = 400
        else:
            io_steps = 100

        run_case(
            self.skip_diffusion, 
            self.size, 
            self.steps, 
            self.dt, 
            self.restarts, 
            self.reynolds_number, 
            self.path,
            iteration,
            io_steps
        )

    def __repr__(self):
        return f"RunCase N={self.size} | {skip_diffusion_to_str(self.skip_diffusion)} | dt={self.dt} | restarts = {self.restarts} | steps = {self.steps} | Re = {self.reynolds_number}"

def main():
    if os.path.exists(BASE_SAVE):
        shutil.rmtree(BASE_SAVE)
    os.mkdir(BASE_SAVE)
    create_file(BASE_SAVE + "/errors.txt")

    run_shell_command("make")
    i = 0

    cases = [
        # reynolds
        RunCase(skip_diffusion=0,size=64, dt=0.0005, steps=20000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/compare_reynolds/re40'),
        RunCase(skip_diffusion=0,size=64, dt=0.0005, steps=20000, restarts=3, reynolds_number=80, path= BASE_SAVE + '/compare_reynolds/re80'),
        # restarts
        RunCase(skip_diffusion=0,size=128, dt=0.001, steps=10000, restarts=0, reynolds_number=40, path= BASE_SAVE + '/restarts/0restarts'),
        RunCase(skip_diffusion=0,size=128, dt=0.001, steps=10000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/restarts/3restarts'),
        RunCase(skip_diffusion=0,size=128, dt=0.001, steps=10000, restarts=6, reynolds_number=40, path= BASE_SAVE + '/restarts/6restarts'),
        # size
        RunCase(skip_diffusion=0,size=64, dt=0.001, steps=10000, restarts=3, reynolds_number=40,  path= BASE_SAVE + '/size/N64'),
        RunCase(skip_diffusion=0,size=128, dt=0.001, steps=10000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/size/N128'),
        # exploding reynolds
        RunCase(skip_diffusion=0,size=128, dt=0.0005, steps=20000, restarts=3, reynolds_number=40,  path= BASE_SAVE + '/explode_reynolds/re40'),
        RunCase(skip_diffusion=0,size=128, dt=0.0005, steps=20000, restarts=3, reynolds_number=120, path= BASE_SAVE + '/explode_reynolds/re120'),
        RunCase(skip_diffusion=0,size=128, dt=0.0005, steps=20000, restarts=3, reynolds_number=200, path= BASE_SAVE + '/explode_reynolds/re200'),
        RunCase(skip_diffusion=0,size=128, dt=0.0005, steps=20000, restarts=3, reynolds_number=280, path= BASE_SAVE + '/explode_reynolds/re280'),
        RunCase(skip_diffusion=0,size=128, dt=0.0005, steps=20000, restarts=3, reynolds_number=320, path= BASE_SAVE + '/explode_reynolds/re360'),
        # different dt
        RunCase(skip_diffusion=0,size=128, dt=0.001  , steps=10000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/dt/001'),
        RunCase(skip_diffusion=0,size=128, dt=0.0005 , steps=20000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/dt/0005'),
        RunCase(skip_diffusion=0,size=128, dt=0.00025, steps=40000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/dt/00025'),
    ]
    
    for case in cases:
        try:
            case.run(i)
        except Exception as e:

            traceback_str = ''.join(traceback.format_tb(e.__traceback__))
            with open(BASE_SAVE + "/errors.txt", 'a') as f:
                print(f"failed for case {case}")
                f.write(f"failed for case {case} {case.path}\ntraceback:\n{traceback_str}")

        i += 1
        print(f"{i/len(cases)*100}% done with cases")

# skip diffusion - 0 / 1 - whether or not to skip the diffusion calculations in rhs_velocity
# size param - the number of divisions in each dimension (size_param=64 means 64x64x64 slices)
# steps - the number of steps the solver will take
# dt - the constant time step used in the solver
# restarts - the number of restarts the solver has, each restart lasts 1 second
# reynolds_number - the reynolds number of the simulation (float)
# save_folder - a hard coded save folder. if none is provided then one will be generated based on the parameters
def run_case(skip_diffusion_param, size_param, steps,dt, restarts, reynolds_number, save_folder=None, iteration=1, steps_between_io=None):
    if save_folder is None:
        save_folder = BASE_SAVE + f"/{size_param}N-dt{dt}-{skip_diffusion_to_str(skip_diffusion_param)}-{restarts}-restarts-re{reynolds_number}-steps{steps}"

    if steps_between_io is None:
        steps_between_io = 100

    # delete the entire folder and remake it
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    os.makedirs(save_folder)

    print(f"creating a config for N={size_param} | {skip_diffusion_to_str(skip_diffusion_param)} | dt={dt} | restarts = {restarts} | steps = {steps}")
    run_shell_command(f"hit3d-config --n {size_param} --steps {steps} --steps-between-io {steps_between_io} --flow-type 0 --skip-diffusion {skip_diffusion_param} --dt -{dt} --restarts {restarts} --reynolds {reynolds_number} input_file.in ")

    restart_time_slice = restarts * 1.

    if restart_time_slice > steps*dt:
        print(steps*dt, restart_time_slice)
        raise ValueError("restarts last longer than the simulation - this will cause a crash down the line")

    run_hit3d()

    # save simple run info to a json file

    # start perforing postprocessessing stuff on another thread
    postprocessing("output/", save_folder, restart_time_slice, steps, dt, iteration)
    #p = Process(target=postprocessing, args= (solver_folder, save_folder, restart_time_slice, steps, dt, iteration,))
    #p.start()

def run_hit3d():
    clean_output_dir()

    run_shell_command(f'mpirun -np {MPI_PROC} --mca opal_warn_on_missing_libcuda 0 ./hit3d.x "input_file" "nosplit"')

def postprocessing(solver_folder, output_folder, restart_time_slice, steps, dt, iteration):
    shutil.move("input_file.in", output_folder + "/input_file.in")

    os.mkdir(output_folder + '/flowfield')

    flowfield_files = [i for i in os.listdir(f"{solver_folder}/velocity_field") if i !=".gitignore"]
    print(flowfield_files)

    # group each of the flowfield files by the timestep that they belong to
    # so that we can combine all of the files that are from the same t
    groupings = {}
    for filename in flowfield_files:
        _mpi_id, timestep = parse_filename(filename)

        if groupings.get(timestep) is None:
            groupings[timestep] = []

        groupings[timestep].append(filename)

    # move all of the files in each group to a tmp directory and process all of the
    # files in that directory with hit3d-utils
    for filegroup in groupings.values():

        #
        # Combine all partial flowfield csv files written by each mpi process into a
        # singular csv file - add q criterion information to the csv file and then
        # re-export that csv to a vtk file for viewing in paraview
        #

        # clear the tmp file
        clean_and_create_folder(f"{solver_folder}/tmp_velo")

        # move all of the current files into the tmp folder to compress them
        for current_file in filegroup:
            shutil.move(f"{solver_folder}/velocity_field/" + current_file, f"{solver_folder}/tmp_velo/" + current_file)

        _, timestep = parse_filename(filegroup[0])

        concat_csv = f"{solver_folder}/velocity_field/{timestep}.csv"
        combined_csv = f"{solver_folder}/combined_csv.csv"
        # file is in a directory accessable by the paraview headless renderer
        vtk_save = output_folder + f"/flowfield/timestep_{timestep}.vtk" 

        # concat all of the data together
        run_shell_command(f'hit3d-utils concat {solver_folder}/tmp_velo {solver_folder}/size.csv "{concat_csv}" {combined_csv}')

        # add qcriterion data to the csv
        run_shell_command(f'python3 /home/brooks/github/hit3d-utils/src/q_criterion.py {concat_csv} {combined_csv}')

        # re-export the csv file to a vtk for viewing
        run_shell_command(f'hit3d-utils vtk {concat_csv} {combined_csv} {vtk_save}')

    #
    # run the paraview headless renderer with our new vtk files 
    # and render write the videos to the results folder
    #

    # # TODO: better dynamic framerate here
    # framerate = 20
    # if steps == 20000:
    #     framerate = 20
    # else:
    #     framerate = 10

    # # run paraview over the files
    # run_shell_command("sh /home/brooks/github/hit3d-utils/paraview/run.sh")
    # # animate the frames into a video
    # run_shell_command(f"sh /home/brooks/github/hit3d-utils/paraview/video.sh {render_folder} {vtk_render_results_folder} {framerate}")

    # # cleanup files used by paraview rendering
    # shutil.rmtree(vtk_data_folder)
    # shutil.rmtree(vtk_render_results_folder)

    #
    # Handle time step energy files
    #

    run_shell_command(f"hit3d-utils add {solver_folder}/energy/ {solver_folder}/energy.csv")

    #
    # parse and re-export spectral information
    #

    # the number of points on the spectra that we want
    NUM_DATAPOINTS = 5
    steps_for_restarts = int(restart_time_slice / dt) + 1
    effective_steps = steps - steps_for_restarts
    # how many actual values from es.gp are available after restarts are complete
    datapoints = int(effective_steps / 10)

    stepby = max(int(datapoints / NUM_DATAPOINTS) - 5,1)

    run_shell_command(f"hit3d-utils spectral {solver_folder}/es.gp {solver_folder}/spectra.json --step-by {stepby} --skip-before-time {restart_time_slice}")

    #
    # run plotting for energy / helicity information from energy.csv
    # and for the spectra
    #

    run_shell_command(f'python3 /home/brooks/github/hit3d-utils/plots/energy_helicity.py {solver_folder}/energy.csv {output_folder} "{restart_time_slice}"')
    run_shell_command(f'python3 /home/brooks/github/hit3d-utils/plots/spectra.py {solver_folder}/spectra.json {output_folder}')

    # move some of the important files to the save folder so they do not get purged
    shutil.move(f"{solver_folder}/energy.csv", output_folder + '/energy.csv')
    shutil.move(f"{solver_folder}/es.gp", output_folder + '/es.gp')

# parse csv files for flowfield output by fortran
def parse_filename(filename):
    mpi_id = filename[0:2]
    timestep = filename[3:8]
    return mpi_id,timestep

def clean_and_create_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)

# clean out the output directory and  place gitignore files in it
def clean_output_dir():
    shutil.rmtree("output")
    os.mkdir("output")
    os.mkdir("output/velocity")
    os.mkdir("output/energy")
    os.mkdir("output/velocity_field")

    create_file("output/.gitignore")
    create_file("output/velocity/.gitignore")
    create_file("output/velocity_field/.gitignore")
    create_file("output/energy/.gitignore")

def create_file(path):
    with open(path,'w') as _:
        pass

def run_shell_command(command):
    print(f"running {command}")
    output = subprocess.run(command,shell=True, check=True)
    if not output.stdout is None:
        print(output.stdout)

def skip_diffusion_to_str(skip_diffusion):
    if skip_diffusion == 1:
        return "no-diffusion"
    elif skip_diffusion == 0:
        return "yes-diffusion"
    else:
        raise ValueError("skip diffusion must be either 0 or 1")

# helpful function for runnning one-off cases
def one_case():
    run_shell_command("make")
    case =  RunCase(skip_diffusion=0,size=128, dt=0.001, steps=100, restarts=0, reynolds_number=40, path= BASE_SAVE + '/testcase_1proc')
    #case =  RunCase(skip_diffusion=0,size=64, dt=0.001, steps=10000, restarts=3, reynolds_number=40, path= BASE_SAVE + '/12proc_10000')
    case.run(0)

if __name__ == "__main__":
    #main()
    one_case()
