import subprocess
import shutil
import os
import traceback
import csv
import json

UNR = True 
IS_DISTRIBUTED = True

if UNR:
    BASE_SAVE = "/home/brooks/sync/hit3d"
else:
    BASE_SAVE = "/home/brooks/lab/hit3d-cases"

if IS_DISTRIBUTED:
    HIT3D_UTILS_BASE = "../../hit3d-utils"
    BASE_SAVE = "/home/brooks/distribute-hit3d"
else:
    if UNR:
        HIT3D_UTILS_BASE = "/home/brooks/github/hit3d-utils"
    else:
        HIT3D_UTILS_BASE = "/home/brooks/github/fluids/hit3d-utils"

IC_SPEC_NAME = "initial_condition_espec.pkg"
IC_WRK_NAME = "initial_condition_wrk.pkg"
IC_JSON_NAME = "initial_condition_vars.json"

class RunCase():
    def __init__(self,skip_diffusion, size, dt, steps, restarts, reynolds_number,path, load_initial_data=0, nprocs=16, export_vtk=False, epsilon1=0.0, epsilon2=0.0, restart_time=1.0, skip_steps=0 ):
        # if there are restarts, find the number of steps spent in that those restarts
        # and add them to the current number of steps
        simulation_time_restart = restarts * restart_time
        steps_restarts = int(simulation_time_restart / dt)
        steps += steps_restarts + skip_steps

        self.skip_diffusion =  skip_diffusion
        self.size =            size
        self.dt =              dt
        self.steps =           steps
        self.restarts =        restarts
        self.reynolds_number = reynolds_number
        self.path =            path
        self.load_initial_data = load_initial_data
        self.nprocs            = nprocs
        self.export_vtk        = export_vtk
        self.epsilon1          = epsilon1
        self.epsilon2          = epsilon2 
        self.restart_time      = restart_time
        self.skip_steps        = skip_steps

    def run(self, iteration):
        # automatically calculate a reasonable number of steps between io 
        # operations (~ 100 every 10_000 steps)

        if self.steps > 1000:
            #io_steps = 100
            io_steps = int(self.steps * 100 / 10_000)
        else:
            io_steps = int(self.steps * 100 / 10_000)

        run_case(
            self.skip_diffusion, 
            self.size, 
            self.steps, 
            self.dt, 
            self.restarts, 
            self.reynolds_number, 
            self.load_initial_data,
            self.nprocs,
            self.path,
            iteration,
            io_steps,
            self.export_vtk,
            self.epsilon1,
            self.epsilon2,
            self.restart_time
        )

    def __repr__(self):
        return f"RunCase N={self.size} | {skip_diffusion_to_str(self.skip_diffusion)} | dt={self.dt} | restarts = {self.restarts} | steps = {self.steps} | Re = {self.reynolds_number} | load-initial-data = {self.load_initial_data}"

    def effective_restart_time(self):
        self.restart_time + (self.skip_steps / self.dt)

    def write_to_json(self, job_name, folder_path):
        file_name = f"{folder_path}/{job_name}.json"

        print(f"writing json file for job `{file_name}`")
        json_data = {
            "skip_diffusion":    self.skip_diffusion,
            "size":              self.size,
            "dt":                self.dt,
            "steps":             self.steps,
            "restarts":          self.restarts,
            "reynolds_number":   self.reynolds_number,
            "path":              self.path,
            "load_initial_data": self.load_initial_data,
            "nprocs":            self.nprocs,
            "export_vtk":        self.export_vtk,
            "epsilon1":          self.epsilon1,
            "epsilon2":          self.epsilon2,
            "restart_time":      self.restart_time,
            "skip_steps":        self.skip_steps,
            "job_name": job_name
        }

        with open(file_name, "w", encoding="utf-8") as file:
            json.dump(json_data, file)

# skip diffusion - 0 / 1 - whether or not to skip the diffusion calculations in rhs_velocity
# size param - the number of divisions in each dimension (size_param=64 means 64x64x64 slices)
# steps - the number of steps the solver will take
# dt - the constant time step used in the solver
# restarts - the number of restarts the solver has, each restart lasts 1 second
# reynolds_number - the reynolds number of the simulation (float)
# save_folder - a hard coded save folder. if none is provided then one will be generated based on the parameters
def run_case(
    skip_diffusion_param, size_param, steps,dt, 
    restarts, reynolds_number, load_initial_data, 
    nprocs, save_folder, iteration, 
    steps_between_io, export_vtk, epsilon1, epsilon2, 
    restart_time):
    if save_folder is None:
        save_folder = BASE_SAVE + f"/{size_param}N-dt{dt}-{skip_diffusion_to_str(skip_diffusion_param)}-{restarts}-restarts-re{reynolds_number}-steps{steps}"

    if steps_between_io is None:
        steps_between_io = 100

    # delete the entire folder and remake it
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)

    os.makedirs(save_folder)

    print(f"creating a config for N={size_param} | {skip_diffusion_to_str(skip_diffusion_param)} | dt={dt} | restarts = {restarts} | steps = {steps}")
    run_shell_command(f"hit3d-utils config-generator --n {size_param} --steps {steps} --steps-between-io {steps_between_io} --flow-type 0 --skip-diffusion {skip_diffusion_param} --dt -{dt} --restarts {restarts} --reynolds {reynolds_number} --initial-condition {load_initial_data} --epsilon1 {epsilon1} --epsilon2 {epsilon2} --restart-time {restart_time} input_file.in ")

    restart_time_slice = restarts * 1.

    if restart_time_slice > steps*dt:
        print(steps*dt, restart_time_slice)
        raise ValueError("restarts last longer than the simulation - this will cause a crash down the line")

    run_hit3d(nprocs)

    postprocessing("output/", save_folder, restart_time_slice, steps, dt, export_vtk)

    if load_initial_data == 1:
        organize_initial_condition(save_folder)

def run_hit3d(nprocs):
    clean_output_dir()

    run_shell_command(f'mpirun -np {nprocs} --mca opal_warn_on_missing_libcuda 0 ./hit3d.x "input_file" "nosplit"')

    print("finished hit3d, back in python")

# we have just run a process to generate an initial condition for the other datasets
# this function organizes the outputs so that they may be used by other processes
def organize_initial_condition(save_folder):
    #run_shell_command(f"hit3d-utils add {solver_folder}/energy/ {solver_folder}/energy.csv")

    with open(save_folder + "/energy.csv", "r") as file:
        reader = csv.reader(file)

        first = True

        for row in reader:
            if first:
                first=False
                continue

            energy = float(row[1])
            helicity = float(row[3])
            fdot_u = float(row[5])
            fdot_h = float(row[6])

    if first==True:
        raise ValueError("there were no rows in the energy.csv outputted by the solver")

    EpsilonControl.to_json(energy, helicity, fdot_u, fdot_h)

    if IS_DISTRIBUTED:
        file = IC_JSON_NAME
        shutil.copy(file, save_folder + "/" + file)
        file = IC_SPEC_NAME
        shutil.copy(file, save_folder + "/" + file)
        file = IC_WRK_NAME
        shutil.copy(file, save_folder + "/" + file)


class EpsilonControl():
    def __init__(self):
        with open(IC_JSON_NAME, "r") as file:
            data = json.load(file)

        self.energy = data["energy"]
        self.helicity = data["helicity"]
        self.fdot_u = data["fdot_u"]
        self.fdot_h = data["fdot_h"]

    @staticmethod
    def to_json(energy, helicity, fdot_u, fdot_h):
        data = {
            "energy": energy,
            "helicity": helicity,
            "fdot_u": fdot_u,
            "fdot_h": fdot_h,
        }

        print("IC json data is")
        print(data)

        with open(IC_JSON_NAME, "w") as file:
            json.dump(data, file)

    @staticmethod
    def load_json():
        return EpsilonControl()

    def epsilon_1(self,delta):
        return delta * self.energy / self.fdot_u

    def epsilon_2(self,delta):
        return delta * self.helicity/ self.fdot_h

def postprocessing(solver_folder, output_folder, restart_time_slice, steps, dt, save_vtk):
    shutil.move("input_file.in", output_folder + "/input_file.in")
    shutil.move(f"{solver_folder}/slice", output_folder + "/fortran_slice_data")

    os.mkdir(output_folder + '/flowfield')

    flowfield_files = [i for i in os.listdir(f"{solver_folder}/velocity_field") if i !=".gitignore"]
    print(flowfield_files)

    if save_vtk:
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
            run_shell_command(f'python3 {HIT3D_UTILS_BASE}/src/q_criterion.py {concat_csv} {combined_csv}')

            # re-export the csv file to a vtk for viewing
            run_shell_command(f'hit3d-utils vtk {concat_csv} {combined_csv} {vtk_save}')

            #run_shell_command(f'python3 {HIT3D_UTILS_BASE}/src/slice_vtk.py vorticity {output_folder}')
            #run_shell_command(f'python3 {HIT3D_UTILS_BASE}/src/slice_vtk.py forcing {output_folder}')

    #
    # Handle time step energy files
    #

    #run_shell_command(f"hit3d-utils add {solver_folder}/energy/ {solver_folder}/energy.csv")

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

    run_shell_command(f'python3 {HIT3D_UTILS_BASE}/plots/energy_helicity.py {solver_folder}/energy.csv {output_folder} "{restart_time_slice}"')

    # move some of the important files to the save folder so they do not get purged
    shutil.move(f"{solver_folder}/energy.csv", output_folder + '/energy.csv')
    shutil.move(f"{solver_folder}/es.gp", output_folder + '/es.gp')
    shutil.move(f"{solver_folder}/spectra.json", output_folder + '/spectra.json')

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
    if os.path.exists("output"):
        shutil.rmtree("output")
    os.mkdir("output")
    os.mkdir("output/velocity")
    os.mkdir("output/energy")
    os.mkdir("output/velocity_field")
    os.mkdir("output/slice")

    create_file("output/.gitignore")
    create_file("output/velocity/.gitignore")
    create_file("output/velocity_field/.gitignore")
    create_file("output/energy/.gitignore")
    create_file("output/slice/.gitignore")

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

def mpi_routine_study():
    run_shell_command("make")

    skip_diffusion = 1

    case =  RunCase(skip_diffusion=skip_diffusion,size=64, dt=0.001, steps=1, restarts=3, reynolds_number=40, path= BASE_SAVE + '/vary_proc/initial_field', load_initial_data=1, nprocs=1)
    case.run(0)

    for i in range(1,17):
        if 64 %i ==0:
            print(i)
    
            case =  RunCase(skip_diffusion=skip_diffusion,size=64, dt=0.001, steps=10000, restarts=0, reynolds_number=40, path= BASE_SAVE + f'/vary_proc/{i}proc_10000', nprocs=i)
            case.run(0)

def load_spectra_initial_condition():
    run_shell_command("make")

    skip_diffusion = 0
    case =  RunCase(skip_diffusion=skip_diffusion, size=64, dt=0.001, steps=1, restarts=3, reynolds_number=40, path= BASE_SAVE + '/spectra_ic/initial_field', load_initial_data=1)
    case.run(0)

    case =  RunCase(skip_diffusion=skip_diffusion, size=64, dt=0.001, steps=10_000, restarts=0, reynolds_number=40, path= BASE_SAVE + '/spectra_ic/derived_field', load_initial_data=0)
    case.run(1)

    case =  RunCase(skip_diffusion=skip_diffusion, size=64, dt=0.001, steps=10_000, restarts=0, reynolds_number=40, path= BASE_SAVE + '/spectra_ic/no_restarts_no_ic', load_initial_data=2)
    case.run(1)

def wrap_error_case(case, filepath):
        try:
            case.run(0)
        except Exception as e:
            traceback_str = ''.join(traceback.format_tb(e.__traceback__))
            with open(filepath, 'a') as f:
                print(f"failed for case {case}")
                f.write(f"failed for case {case} {case.path}\ntraceback:\n{traceback_str}")

def initial_condition():
    dt = 0.0005
    size = 128
    re = 40
    forcing_folder = "forcing_0005_dt_longer_steps_40_000"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"
    output_folder = f"../../distribute_save"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    case =  RunCase(skip_diffusion=0, size=size, dt=dt, steps=15_000*2, restarts=0, reynolds_number=re, path=f'{output_folder}/initial_field' ,load_initial_data=1, restart_time=1.)
    case.write_to_json("initial_condition", save_json_folder)

    if UNR:
        build_location= "/home/brooks/github/hit3d-utils/build.py"
        nodes_location = "/home/brooks/distribute/distribute-nodes.yaml"
        run_py = "/home/brooks/github/hit3d/src/run.py"
    else:
        build_location = "/home/brooks/github/fluids/hit3d-utils/build.py"
        nodes_location = "/home/brooks/github/fluids/distribute/run/server/distribute-nodes.yaml"
        run_py = "/home/brooks/github/fluids/hit3d/src/run.py"

    run_shell_command(f"hit3d-utils distribute-gen --output-folder {save_json_folder} --library {run_py} --library-save-name hit3d_helpers.py --batch-name hit3d_initial_condition {save_json_folder}/*.json")

    shutil.copy(build_location, f"{save_json_folder}/build.py")
    shutil.copy(nodes_location, f"{save_json_folder}/distribute-nodes.yaml")

# in order to calculate forcing cases we need to have an initial condition file
def forcing_cases():
    run_shell_command("make")
    forcing_folder = "forcing_cases"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    dt = 0.0005
    size = 128
    re = 40
    steps = 20_000 * 4
    save_vtk = True
    batch_name = "hit3d_forcing"

    epsilon_generator = EpsilonControl.load_json()

    delta_1 = 0.01
    delta_2 = 0.02

    cases = [
        [0., 0., "baseline"],

        #epsilon 1 cases
        [delta_1, 0., "ep1-pos"],
        [-1*delta_1, 0., "ep1-neg"],

        # epsilon 2  cases
        [ 0., delta_2, "ep2-pos"],
        [ 0., -1*delta_2, "ep2-neg"],

        # both ep1 and ep2 cases
        [ delta_1, delta_2, "epboth-pos"],
        [ -1*delta_1, -1 * delta_1, "epboth-neg"],
    ]

    for skip_diffusion in [0,1]:
        for delta_1, delta_2, folder in cases:
            diffusion_str = skip_diffusion_to_str(skip_diffusion)
            epsilon1 = epsilon_generator.epsilon_1(delta_1)
            epsilon2 = epsilon_generator.epsilon_2(delta_2)

            if epsilon1 == -0.0:
                epsilon1 = 0.0
            if epsilon2 == -0.0:
                epsilon2 = 0.0

            output_folder = f"../../distribute_save/{diffusion_str}/{folder}"

            case =  RunCase(skip_diffusion=skip_diffusion, 
                size=size,
                dt=dt,
                steps=steps,
                restarts=0,
                restart_time=1.,
                reynolds_number=re,
                path= output_folder,
                load_initial_data=0,
                epsilon1=epsilon1,
                epsilon2=epsilon2,
                export_vtk=save_vtk
            )

            case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    if UNR:
        build_location= "/home/brooks/github/hit3d-utils/build.py"
        nodes_location = "/home/brooks/distribute/distribute-nodes.yaml"
        run_py = "/home/brooks/github/hit3d/src/run.py"
    else:
        build_location = "/home/brooks/github/fluids/hit3d-utils/build.py"
        nodes_location = "/home/brooks/github/fluids/distribute/run/server/distribute-nodes.yaml"
        run_py = "/home/brooks/github/fluids/hit3d/src/run.py"

    run_shell_command(f"hit3d-utils distribute-gen --output-folder {save_json_folder} --library {run_py} --library-save-name hit3d_helpers.py --batch-name {batch_name} --required-files {IC_SPEC_NAME} --required-files {IC_WRK_NAME} {save_json_folder}/*.json")

    shutil.copy(build_location, f"{save_json_folder}/build.py")
    shutil.copy(nodes_location, f"{save_json_folder}/distribute-nodes.yaml")

    shutil.copy(IC_SPEC_NAME, f"{save_json_folder}/{IC_SPEC_NAME}")
    shutil.copy(IC_WRK_NAME, f"{save_json_folder}/{IC_WRK_NAME}")

def resolution_study():
    run_shell_command("make")
    N = 256

    num_stencil = [ int(i*N) for i in [
        1/2,
        1,
        2,
    ]]

    dt = 0.0005
    re = 80
    steps = 10_000 
    initial_steps = 25_000
    resolution_folder = "resolution_study_3"

    for n in num_stencil:
        for skip_diffusion in [0,1]:
            diff_str = skip_diffusion_to_str(skip_diffusion)
            folder_name = f"{n}"

            case = RunCase(
                skip_diffusion=skip_diffusion,
                size = n,
                dt = dt,
                steps = steps,
                restarts = 0,
                reynolds_number=re,
                path = f"{BASE_SAVE}/{resolution_folder}/{diff_str}/{folder_name}",
                load_initial_data=2,
                export_vtk=True,
                skip_steps=initial_steps
            )

            wrap_error_case(case, f"{BASE_SAVE}/{resolution_folder}/errors.txt")

def temporal_study():
    run_shell_command("make")

    DT = 0.0005
    STEPS = 10_000 
    INITIAL_STEPS = 25_000

    timesteps = [ [i*DT, int(STEPS / i), int(INITIAL_STEPS / i)] for i in [
        1/2,
        1,
        2,
    ]]

    N = 256
    re = 80
    temporal_folder = "temporal_study_3"

    print("timesteps are", timesteps)

    for dt, steps, initial_steps in timesteps:
        for skip_diffusion in [0,1]:
            diff_str = skip_diffusion_to_str(skip_diffusion)
            folder_name = f"{dt}"

            case = RunCase(
                skip_diffusion=skip_diffusion,
                size = N,
                dt = dt,
                steps = steps,
                restarts = 0,
                reynolds_number=re,
                path = f"{BASE_SAVE}/{temporal_folder}/{diff_str}/{folder_name}",
                load_initial_data=2,
                export_vtk=True,
                skip_steps=initial_steps
            )

            wrap_error_case(case, f"{BASE_SAVE}/{temporal_folder}/errors.txt")

# helpful function for runnning one-off cases
def one_case():

    # 25_000 for 64
    #
    run_shell_command("make")
    case =  RunCase(
        skip_diffusion=1,
        size=128,
        dt=0.001,
        steps=35000,
        restarts=0,
        reynolds_number=80,
        path= BASE_SAVE + f'/spectra_study_128',
        load_initial_data=2,
        epsilon1=0.0000,
        epsilon2=0.0000,
    )

    case.run(0)

def remove_restart_files():
    for i in ["initial_condition_espec.pkg", "initial_condition_wrk.pkg", "initial_condition_vars.json"]:
        if os.path.exists(i):
            os.remove(i) 

if __name__ == "__main__":
    initial_condition()
    forcing_cases()
