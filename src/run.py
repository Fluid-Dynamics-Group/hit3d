import subprocess
import shutil
import os
import traceback
import csv
import json
from pprint import pprint
from glob import glob

UNR = True
#IS_DISTRIBUTED = True
IS_DISTRIBUTED = False
IS_SINGULARITY = False

if UNR:
    BASE_SAVE = "/home/brooks/sync/hit3d"
    INITIAL_COND_FOLDER = "/home/brooks/distribute-cases/initial_consitions_1/initial_consitions_1"
else:
    BASE_SAVE = "/home/brooks/lab/hit3d-cases"
    # TODO: pull these cases locally
    #INITIAL_COND_FOLDER = "/home/brooks/lab/sync/initial_conditions"

if IS_DISTRIBUTED:
    HIT3D_UTILS_BASE = "../../hit3d-utils"
    if UNR:
        BASE_SAVE = "/home/brooks/distribute-cases"
    else:
        BASE_SAVE = "/home/brooks/lab/distribute-cases"
else:
    if UNR:
        HIT3D_UTILS_BASE = "/home/brooks/github/hit3d-utils"
    else:
        HIT3D_UTILS_BASE = "/home/brooks/github/fluids/hit3d-utils"

IC_SPEC_NAME = "initial_condition_espec.pkg"
IC_WRK_NAME = "initial_condition_wrk.pkg"
IC_JSON_NAME = "initial_condition_vars.json"


class RunCase:
    def __init__(
        self,
        skip_diffusion,
        size,
        dt,
        steps,
        restarts,
        reynolds_number,
        path,
        load_initial_data=0,
        nprocs=16,
        export_vtk=False,
        epsilon1=0.0,
        epsilon2=0.0,
        restart_time=1.0,
        skip_steps=0,
        scalar_type=0,
        validate_viscous_compensation=0,
        viscous_compensation=0,
        require_forcing=0,
        io_steps=None,
        export_divergence=0,
        nscalar=0,
    ):

        if steps >= 100_000:
            raise ValueError(
                "steps cannot currently be a number with more than 5 digits due to how data is outputted"
            )

        # if there are restarts, find the number of steps spent in that those restarts
        # and add them to the current number of steps
        simulation_time_restart = restarts * restart_time
        steps_restarts = int(simulation_time_restart / dt)
        steps += steps_restarts + skip_steps

        self.skip_diffusion = skip_diffusion
        self.size = size
        self.dt = dt
        self.steps = steps
        self.restarts = restarts
        self.reynolds_number = reynolds_number
        self.path = path
        self.load_initial_data = load_initial_data
        self.nprocs = nprocs
        self.export_vtk = export_vtk
        self.epsilon1 = epsilon1
        self.epsilon2 = epsilon2
        self.restart_time = restart_time
        self.skip_steps = skip_steps
        self.scalar_type = scalar_type
        self.export_divergence = export_divergence

        self.validate_viscous_compensation = validate_viscous_compensation
        self.viscous_compensation = viscous_compensation

        self.require_forcing = require_forcing
        self.io_steps = io_steps
        self.nscalar = nscalar

    def validate_params(self):
        if (
            self.viscous_compensation == 1 or self.validate_viscous_compensation == 1
        ) and self.skip_diffusion == 1:
            raise ValueError("cant have viscous compensation for the inviscid flow")

    def run(self, iteration):
        # automatically calculate a reasonable number of steps between io
        # operations (~ 100 every 10_000 steps)

        # keep the amount of data produced constant
        # higher number = more steps saved
        if self.io_steps is None:
            io_steps = int(self.steps * 600 / 80_000)
        else:
            io_steps = self.io_steps

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
            self.restart_time,
            self.scalar_type,
            self.validate_viscous_compensation,
            self.viscous_compensation,
            self.require_forcing,
            self.export_divergence,
            self.nscalar,
        )

    def __repr__(self):
        return f"RunCase N={self.size} | {skip_diffusion_to_str(self.skip_diffusion)} | dt={self.dt} | restarts = {self.restarts} | steps = {self.steps} | Re = {self.reynolds_number} | load-initial-data = {self.load_initial_data}"

    def effective_restart_time(self):
        self.restart_time + (self.skip_steps / self.dt)

    def write_to_json(self, job_name, folder_path):
        file_name = f"{folder_path}/{job_name}.json"

        print(f"writing json file for job `{file_name}`")
        json_data = {
            "skip_diffusion": self.skip_diffusion,
            "size": self.size,
            "dt": self.dt,
            "steps": self.steps,
            "restarts": self.restarts,
            "reynolds_number": self.reynolds_number,
            "path": self.path,
            "load_initial_data": self.load_initial_data,
            "nprocs": self.nprocs,
            "export_vtk": self.export_vtk,
            "epsilon1": self.epsilon1,
            "epsilon2": self.epsilon2,
            "restart_time": self.restart_time,
            "skip_steps": self.skip_steps,
            "scalar_type": self.scalar_type,
            "job_name": job_name,
            "validate_viscous_compensation": self.validate_viscous_compensation,
            "viscous_compensation": self.viscous_compensation,
            "require_forcing": self.require_forcing,
            "io_steps": self.io_steps,
            "export_divergence": self.export_divergence,
            "nscalar": self.nscalar,
        }

        with open(file_name, "w", encoding="utf-8") as file:
            json.dump(json_data, file)


# define some schema so that we can more dynamically build
# the correct branches of code
class Build:
    def __init__(self, hit3d_branch, hit3d_utils_branch):

        self.hit3d_branch = hit3d_branch
        self.hit3d_utils_branch = hit3d_utils_branch

    def to_json(self, output_folder):
        json_data = {
            "hit3d_branch": self.hit3d_branch,
            "hit3d_utils_branch": self.hit3d_utils_branch,
        }

        with open(output_folder + "/build.json", "w", encoding="utf-8") as file:
            json.dump(json_data, file)

    @staticmethod
    def load_json(from_folder):
        json_path = from_folder + "/build.json"
        print("loading json from ", json_path)

        with open(json_path, "r", encoding="utf-8") as file:
            json_data = json.load(file)

            hit3d_branch = json_data["hit3d_branch"]
            hit3d_utils_branch = json_data["hit3d_utils_branch"]

            return Build(hit3d_branch, hit3d_utils_branch)


# skip diffusion - 0 / 1 - whether or not to skip the diffusion calculations in rhs_velocity
# size param - the number of divisions in each dimension (size_param=64 means 64x64x64 slices)
# steps - the number of steps the solver will take
# dt - the constant time step used in the solver
# restarts - the number of restarts the solver has, each restart lasts 1 second
# reynolds_number - the reynolds number of the simulation (float)
# save_folder - a hard coded save folder. if none is provided then one will be generated based on the parameters
def run_case(
    skip_diffusion_param,
    size_param,
    steps,
    dt,
    restarts,
    reynolds_number,
    load_initial_data,
    nprocs,
    save_folder,
    iteration,
    steps_between_io,
    export_vtk,
    epsilon1,
    epsilon2,
    restart_time,
    scalar_type,
    viscous_compensation_validation,
    viscous_compensation,
    require_forcing,
    export_divergence,
    nscalar,
):

    if IS_DISTRIBUTED:
        print("running in distributed mode - singularity: ", IS_SINGULARITY)

    if save_folder is None:
        save_folder = (
            BASE_SAVE
            + f"/{size_param}N-dt{dt}-{skip_diffusion_to_str(skip_diffusion_param)}-{restarts}-restarts-re{reynolds_number}-steps{steps}"
        )

    if steps_between_io is None:
        steps_between_io = 100

    # we can only remove the directories if we are not
    # using singularity to run the jobs
    if not IS_SINGULARITY:
        # delete the entire folder and remake it
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)

        os.makedirs(save_folder)

    print(
        f"creating a config for N={size_param} | {skip_diffusion_to_str(skip_diffusion_param)} | dt={dt} | restarts = {restarts} | steps = {steps}"
    )
    run_shell_command(
        f"hit3d-utils config-generator \
            --n {size_param} \
            --steps {steps} \
            --steps-between-io {steps_between_io} \
            --flow-type 0 \
            --skip-diffusion {skip_diffusion_param} \
            --dt -{dt} --restarts {restarts} \
            --reynolds {reynolds_number} \
            --initial-condition {load_initial_data} \
            --epsilon1 {epsilon1} \
            --epsilon2 {epsilon2} \
            --restart-time {restart_time} \
            --tscalar -0.1 \
            --nscalar {nscalar} \
            --scalar-type {scalar_type} \
            --validate-viscous-compensation {viscous_compensation_validation} \
            --viscous-compensation {viscous_compensation} \
            --require-forcing {require_forcing} \
            --export-divergence {export_divergence} \
            input_file.in "
    )

    restart_time_slice = restarts * 1.0

    if restart_time_slice > steps * dt:
        print(steps * dt, restart_time_slice)
        raise ValueError(
            "restarts last longer than the simulation - this will cause a crash down the line"
        )

    try:
        run_hit3d(nprocs)
    except Exception as e:
        print(
            "exception running hit3d: {}\n ----> calling postprocessessing anyway".format(
                e
            )
        )

    postprocessing(
        "output/",
        save_folder,
        restart_time_slice,
        steps,
        dt,
        export_vtk,
        size_param,
        epsilon1,
        epsilon2,
        export_divergence,
        steps_between_io,
    )

    if load_initial_data == 1:
        organize_initial_condition(save_folder, epsilon1, epsilon2)


def run_hit3d(nprocs):
    # if we are not using singularity then we are responsible for cleaning
    # up our own mess
    if not IS_SINGULARITY:
        clean_output_dir()

    run_shell_command(
        f'mpirun -np {nprocs} --mca opal_warn_on_missing_libcuda 0 ./hit3d.x "input_file" "nosplit"'
    )

    print("finished hit3d, back in python")


# we have just run a process to generate an initial condition for the other datasets
# this function organizes the outputs so that they may be used by other processes
def organize_initial_condition(save_folder, epsilon_1, epsilon_2):
    # run_shell_command(f"hit3d-utils add {solver_folder}/energy/ {solver_folder}/energy.csv")

    with open(save_folder + "/energy.csv", "r") as file:
        reader = csv.reader(file)

        first = True

        for row in reader:
            if first:
                first = False
                continue

            energy = float(row[1])
            helicity = float(row[3])

            # these need to change if the CSV every changes
            fdot_u_l = float(row[5])
            fdot_u_r = float(row[6])

            fdot_h_l = float(row[7])
            fdot_h_r = float(row[8])

    if first == True:
        raise ValueError("there were no rows in the energy.csv outputted by the solver")

    fdot_u = fdot_u_l + fdot_u_r
    fdot_h = fdot_h_l + fdot_h_r

    EpsilonControl.to_json(energy, helicity, fdot_u, fdot_h)

    if IS_DISTRIBUTED:
        file = IC_JSON_NAME
        shutil.copy(file, save_folder + "/" + file)
        file = IC_SPEC_NAME
        shutil.copy(file, save_folder + "/" + file)
        file = IC_WRK_NAME
        shutil.copy(file, save_folder + "/" + file)


class EpsilonControl:
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

    def epsilon_1(self, delta):
        return adjust_negative_epsilon(delta * self.energy / self.fdot_u)

    def epsilon_2(self, delta):
        return adjust_negative_epsilon(delta * self.helicity / self.fdot_h)


def postprocessing(
    solver_folder,
    output_folder,
    restart_time_slice,
    steps,
    dt,
    save_vtk,
    size,
    epsilon_1,
    epsilon_2,
    export_divergence,
    io_steps,
):

    # move viscous compensation files if they are present
    if os.path.exists(f"{solver_folder}/dE_dt_history.binary"):
        dest = f"{output_folder}/dE_dt_history.binary"
        print("renaming dE_dt binary to", dest)
        os.rename(f"{solver_folder}/dE_dt_history.binary", dest)
    if os.path.exists(f"{solver_folder}/dh_dt_history.binary"):
        dest = f"{output_folder}/dh_dt_history.binary"
        print("renaming dh_dt binary")
        os.rename(f"{solver_folder}/dh_dt_history.binary", dest)

    shutil.move("input_file.in", output_folder + "/input_file.in")
    shutil.copy(f"{solver_folder}/energy.csv", output_folder + "/energy.csv")

    os.mkdir(output_folder + "/flowfield")
    clean_and_create_folder(f"{output_folder}/scalars")
    clean_and_create_folder(f"{output_folder}/fortran_slice_data")

    if export_divergence == 1:
        clean_and_create_folder(f"{output_folder}/divergences")
        div_files = [
            i for i in os.listdir(f"{solver_folder}/divergences") if i != ".gitignore"
        ]

        for filename in div_files:
            timestep = parse_filename(filename)
            run_shell_command(
                f"hit3d-utils vtk-div {solver_folder}/divergences/{filename} {size} {output_folder}/divergences/div_{timestep:05}.vtr"
            )

    flowfield_files = [
        i for i in os.listdir(f"{solver_folder}/velocity_field") if i != ".gitignore"
    ]
    scalar_files = [
        i for i in os.listdir(f"{solver_folder}/scalars") if i != ".gitignore"
    ]

    # copy the fortran logging files
    logs_dir = f"{output_folder}/logs/"
    clean_and_create_folder(logs_dir)
    for file in glob(f"{solver_folder}/d*.txt"):
        shutil.move(file, logs_dir)

    # plots from energy.csv
    run_shell_command(
        f'python3 {HIT3D_UTILS_BASE}/plots/energy_helicity.py {solver_folder}/energy.csv {output_folder} "{restart_time_slice}" {epsilon_1} {epsilon_2}'
    )

    print(flowfield_files)

    # if we are saving the vtk data then rip it from the csv files, and delete the csv file once donw
    if save_vtk:
        for filename in flowfield_files:
            timestep = parse_filename(filename)

            csv_file = f"{solver_folder}/velocity_field/{filename}"
            vtk_file = f"{output_folder}/flowfield/flow_{timestep:05}.vtk"
            print_csv_file_header(csv_file, 5)

            run_shell_command(f"hit3d-utils vtk  {csv_file} {size} {vtk_file}")

            # os.remove(csv_file)

    # othewise just remove the csv files anyway
    else:
        for filename in flowfield_files:
            os.remove(f"{solver_folder}/velocity_field/{filename}")

    print("now printing the VTK files that have been generated")
    pprint(list(os.listdir(f"{output_folder}/flowfield/")))

    for filename in scalar_files:
        timestep = parse_scalar_name(filename)
        file = f"{solver_folder}/scalars/{filename}"
        run_shell_command(
            f"hit3d-utils scalars {file} {size} {output_folder}/scalars/sc_{timestep:06}.vtk"
        )
        os.remove(file)

    #
    # parse and re-export spectral information
    #

    # the number of points on the spectra that we want
    NUM_DATAPOINTS = 5
    steps_for_restarts = int(restart_time_slice / dt) + 1
    effective_steps = steps - steps_for_restarts
    # how many actual values from es.gp are available after restarts are complete
    # TODO: this value can change depending on how often write_data is called in main.f90
    datapoints = int(effective_steps / (io_steps * 4))

    stepby = max(int(datapoints / NUM_DATAPOINTS) - 5, 1)

    run_shell_command(
        f"hit3d-utils spectral {solver_folder}/es.gp {solver_folder}/spectra.json --step-by {stepby} --skip-before-time {restart_time_slice}"
    )

    #
    # run plotting for energy / helicity information from energy.csv
    # and for the spectra
    #

    run_shell_command(
        f"python3 {HIT3D_UTILS_BASE}/plots/spectra.py {solver_folder}/spectra.json {output_folder}"
    )

    # move some of the important files to the save folder so they do not get purged
    shutil.move(f"{solver_folder}/es.gp", output_folder + "/es.gp")
    shutil.move(f"{solver_folder}/spectra.json", output_folder + "/spectra.json")

    # now animate the VTK files together
    if save_vtk:
        animation_dir = f"{solver_folder}/animation/"
        rendering_dir = f"{output_folder}/renders/"
        os.mkdir(animation_dir)
        os.mkdir(rendering_dir)

        q_criterion_colored_velocity = (
            f"{rendering_dir}/q_criterion/colored_velocity_mag"
        )
        q_criterion_colored_forcing = f"{rendering_dir}/q_criterion/colored_forcing_mag"
        q_criterion_colored_strain = f"{rendering_dir}/q_criterion/colored_strain_mag"
        q_criterion_colored_omega = f"{rendering_dir}/q_criterion/colored_omega_mag"

        os.makedirs(q_criterion_colored_velocity)
        os.makedirs(q_criterion_colored_forcing)
        os.makedirs(q_criterion_colored_strain)
        os.makedirs(q_criterion_colored_omega)

        run_shell_command(f"ip addr | grep '\\.'")
        run_shell_command(f"which pvpython")

        run_shell_command(
            f"pvpython {HIT3D_UTILS_BASE}/paraview/main.py velocity {output_folder}/flowfield {animation_dir} {q_criterion_colored_velocity}"
        )
        run_shell_command(
            f"pvpython {HIT3D_UTILS_BASE}/paraview/main.py forcing {output_folder}/flowfield {animation_dir} {q_criterion_colored_forcing}"
        )
        run_shell_command(
            f"pvpython {HIT3D_UTILS_BASE}/paraview/main.py strain {output_folder}/flowfield {animation_dir} {q_criterion_colored_strain}"
        )
        run_shell_command(
            f"pvpython {HIT3D_UTILS_BASE}/paraview/main.py omega {output_folder}/flowfield {animation_dir} {q_criterion_colored_omega}"
        )


# parse csv files for flowfield output by fortran
# should be the same for both flowfield files and divergence files
def parse_filename(filename):
    timestep = filename[0:5]
    return int(timestep)


def parse_scalar_name(filename):
    timestep = filename[5 : 5 + 6]
    return int(timestep)


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
    os.mkdir("output/divergences")
    os.mkdir("output/slice")
    os.mkdir("output/scalars")

    create_file("output/.gitignore")
    create_file("output/velocity/.gitignore")
    create_file("output/velocity_field/.gitignore")
    create_file("output/energy/.gitignore")
    create_file("output/slice/.gitignore")
    create_file("output/divergences/.gitignore")


def create_file(path):
    with open(path, "w") as _:
        pass

def run_shell_command(command):
    print(f"running {command}")
    output = subprocess.run(
        command,
        shell=True,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    if (not output.stdout is None) and len(output.stdout) != 0:
        print(str(output.stdout.decode()))

    if (not output.stderr is None) and len(output.stderr) != 0:
        print(len(output.stderr.decode()))


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

    case = RunCase(
        skip_diffusion=skip_diffusion,
        size=64,
        dt=0.001,
        steps=1,
        restarts=3,
        reynolds_number=40,
        path=BASE_SAVE + "/vary_proc/initial_field",
        load_initial_data=1,
        nprocs=1,
    )
    case.run(0)

    for i in range(1, 17):
        if 64 % i == 0:
            print(i)

            case = RunCase(
                skip_diffusion=skip_diffusion,
                size=64,
                dt=0.001,
                steps=10000,
                restarts=0,
                reynolds_number=40,
                path=BASE_SAVE + f"/vary_proc/{i}proc_10000",
                nprocs=i,
            )
            case.run(0)


def load_spectra_initial_condition():
    run_shell_command("make")

    skip_diffusion = 0
    case = RunCase(
        skip_diffusion=skip_diffusion,
        size=64,
        dt=0.001,
        steps=1,
        restarts=3,
        reynolds_number=40,
        path=BASE_SAVE + "/spectra_ic/initial_field",
        load_initial_data=1,
    )
    case.run(0)

    case = RunCase(
        skip_diffusion=skip_diffusion,
        size=64,
        dt=0.001,
        steps=10_000,
        restarts=0,
        reynolds_number=40,
        path=BASE_SAVE + "/spectra_ic/derived_field",
        load_initial_data=0,
    )
    case.run(1)

    case = RunCase(
        skip_diffusion=skip_diffusion,
        size=64,
        dt=0.001,
        steps=10_000,
        restarts=0,
        reynolds_number=40,
        path=BASE_SAVE + "/spectra_ic/no_restarts_no_ic",
        load_initial_data=2,
    )
    case.run(1)


def wrap_error_case(case, filepath):
    try:
        case.run(0)
    except Exception as e:
        traceback_str = "".join(traceback.format_tb(e.__traceback__))
        with open(filepath, "a") as f:
            print(f"failed for case {case}")
            f.write(f"failed for case {case} {case.path}\ntraceback:\n{traceback_str}")


def copy_init_files(size):
    if size == 128:
        prefix = (
            INITIAL_COND_FOLDER
            + "/ic_128/"
        )
    elif size == 256:
        prefix = (
            INITIAL_COND_FOLDER
            + "/ic_256/"
        )
    elif size == 64:
        prefix = (
            INITIAL_COND_FOLDER
            + "/ic_64/"
        )
    else:
        raise ValueError("unknown initial condition file formats")

    shutil.copy(prefix + "/" + IC_JSON_NAME, IC_JSON_NAME)
    shutil.copy(prefix + "/" + IC_SPEC_NAME, IC_SPEC_NAME)
    shutil.copy(prefix + "/" + IC_WRK_NAME, IC_WRK_NAME)


def define_output_folder(error_on_local=False):
    if IS_DISTRIBUTED:
        if IS_SINGULARITY:
            return f"/distribute_save/"
        else:
            return "../../distribute_save/"
    else:

        if error_on_local:
            raise ValueError("Cant run this set of cases locally")
        else:
            folder = "./current_results"

            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.mkdir(folder)

            return folder


def ep1_temporal_cases():
    delta_1 = 0.01
    delta_2 = 0.02

    run_shell_command("make")
    forcing_folder = f"temporal_check_{delta_1},{delta_2}"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    dt = 0.0001
    size = 128
    re = 40
    steps = 10_000 * 5
    save_vtk = True

    cases = [
        [dt, steps, "dt_1E-4"],
        [dt / 10, steps * 10, "dt_1E-5"],
        [dt / 100, steps * 100, "dt_1E-6"],
    ]

    epsilon_generator = EpsilonControl.load_json()
    epsilon1 = epsilon_generator.epsilon_1(delta_1)
    epsilon2 = epsilon_generator.epsilon_1(delta_2)

    for dt, steps, folder in cases:
        case = RunCase(
            skip_diffusion=1,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            restart_time=1.0,
            reynolds_number=re,
            path=folder,
            load_initial_data=0,
            epsilon1=epsilon1,
            epsilon2=epsilon2,
            export_vtk=save_vtk,
            scalar_type=14,
        )
        case.write_to_json(folder, save_json_folder)

    copy_distribute_files(save_json_folder, forcing_folder)


def copy_distribute_files(target_folder, batch_name, extra_caps):
    if UNR:
        build_location = "/home/brooks/github/hit3d-utils/build.py"
        run_py = "/home/brooks/github/hit3d/src/run.py"
    else:
        build_location = "/home/brooks/github/fluids/hit3d-utils/build.py"
        run_py = "/home/brooks/github/fluids/hit3d/src/run.py"

    # read in all the json config files (skip build.json)
    files = [i for i in glob(f"{target_folder}/*.json") if not "build.json" in i]
    files = " ".join(files)

    # format the capabilities to a string
    caps_str = ""
    first_cap = True
    for cap in extra_caps:
        if first_cap:
            first_cap = False
        else:
            caps_str += ","

        caps_str += cap

    if len(extra_caps) == 0:
        caps_str = ","

    print("files are \n", files)

    run_shell_command(
        f"hit3d-utils distribute-gen \
            --output-folder {target_folder} \
            --library {run_py} \
            --library-save-name hit3d_helpers.py \
            --batch-name {batch_name} \
            --extra-caps {caps_str} \
            --required-files {IC_SPEC_NAME} \
            --required-files {IC_WRK_NAME} \
            {files}"
    )

    shutil.copy(build_location, f"{target_folder}/build.py")

    shutil.copy(IC_SPEC_NAME, f"{target_folder}/{IC_SPEC_NAME}")
    shutil.copy(IC_WRK_NAME, f"{target_folder}/{IC_WRK_NAME}")

    shutil.copy(f"{HIT3D_UTILS_BASE}/generic_run.py", target_folder)
    shutil.copy(f"{HIT3D_UTILS_BASE}/build.py", target_folder)


def proposal_figures():
    delta_1 = 0.0001
    delta_2 = 0.002

    run_shell_command("make")
    forcing_folder = f"proposal_figures"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"
    extra_caps = ["lab1-2"]

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    dt = 0.0001
    size = 256
    re = 40
    steps = 99_000
    save_vtk = True
    batch_name = forcing_folder

    copy_init_files(size)

    epsilon_generator = EpsilonControl.load_json()

    cases = [
        [-1 * delta_1, 0.0, "ep1-neg"],
        [0.0, -1 * delta_2, "ep2-neg"],
    ]

    output_folder = define_output_folder()

    for skip_diffusion in [0, 1]:

        # TODO: remove this after generating
        if skip_diffusion == 0:
            continue

        for delta_1, delta_2, folder in cases:
            diffusion_str = skip_diffusion_to_str(skip_diffusion)
            epsilon1 = epsilon_generator.epsilon_1(delta_1)
            epsilon2 = epsilon_generator.epsilon_2(delta_2)

            # normalize for -0.0
            epsilon1 = adjust_negative_epsilon(epsilon1)
            epsilon2 = adjust_negative_epsilon(epsilon2)

            case = RunCase(
                skip_diffusion=skip_diffusion,
                size=size,
                dt=dt,
                steps=steps,
                restarts=0,
                restart_time=1.0,
                reynolds_number=re,
                path=output_folder,
                load_initial_data=0,
                epsilon1=epsilon1,
                epsilon2=epsilon2,
                export_vtk=save_vtk,
                scalar_type=14,
            )

            case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)


# remaps a -0.0 value for epsilon to 0.0; otherwise just return the original epsilon value
def adjust_negative_epsilon(epsilon_value: float) -> float:
    if epsilon_value == -0.0:
        return 0.0
    else:
        return epsilon_value


def validate_viscous_compensation():
    batch_name = "viscous_compensation_short_128_final_validation_8"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    size = 128
    steps = 100
    extra_caps = ["lab3"]
    io_steps = int(steps / 10)

    output_folder = define_output_folder(error_on_local=True)

    # if the directory exists remove any older files from the dir
    if os.path.exists(save_json_folder):
        for f in os.listdir(save_json_folder):
            os.remove(os.path.join(save_json_folder, f))

    os.makedirs(save_json_folder, exist_ok=True)

    run_shell_command("make")
    dt = 0.0005
    restarts = 0
    re = 40

    visc_params = [
        # no viscous compensation, use MGM array
        [0, "mgm-forcing"],
        # use viscous compensation
        [1, "visc-compensation"],
    ]

    for skip_diffusion in [0, 1]:
        if skip_diffusion == 0:
            diffusion = "yes-diffusion"
        else:
            diffusion = "no-diffusion"

        for viscous_compensation, case_name in visc_params:
            case_name = f"{case_name}_{diffusion}"

            case = RunCase(
                skip_diffusion=skip_diffusion,
                size=size,
                dt=dt,
                steps=steps,
                restarts=restarts,
                reynolds_number=re,
                path=output_folder,
                load_initial_data=2,
                epsilon1=0.0,
                epsilon2=0.0,
                export_vtk=False,
                scalar_type=14,
                validate_viscous_compensation=viscous_compensation,
                viscous_compensation=viscous_compensation,
                require_forcing=1,
                io_steps=io_steps,
            )

            case.write_to_json(case_name, save_json_folder)

    base_name = "unforced"

    for skip_diffusion in [0, 1]:
        if skip_diffusion == 0:
            diffusion = "yes-diffusion"
        else:
            diffusion = "no-diffusion"

        case_name = f"{base_name}_{diffusion}"

        case = RunCase(
            skip_diffusion=skip_diffusion,
            size=size,
            dt=dt,
            steps=steps,
            restarts=restarts,
            reynolds_number=re,
            path=output_folder,
            load_initial_data=2,
            epsilon1=0.0,
            epsilon2=0.0,
            export_vtk=False,
            scalar_type=14,
            require_forcing=0,
            io_steps=io_steps,
        )

        case.write_to_json(case_name, save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    # write a build.json file so that our code pulls the correct versions that we want

    build = Build("master", "master")
    build.to_json(save_json_folder)


def remove_restart_files():
    for i in [
        "initial_condition_espec.pkg",
        "initial_condition_wrk.pkg",
        "initial_condition_vars.json",
    ]:
        if os.path.exists(i):
            os.remove(i)


# print the first `n` number of lines from a csv file at `path`
# n must be > 0
def print_csv_file_header(path: str, n: int) -> None:
    with open(path, "r") as f:
        lines = f.readlines()

        if len(lines) > n:
            print("csv header:\n")
            for i in range(0, n):
                print(lines[i])


if __name__ == "__main__":
    from cases import forcing_sweep
    from cases import forcing_cases
    from cases import full_system_test
    from cases import figure2
    from cases import one_case
    from cases import one_case
    from cases import track_inviscid_compensation_local
    from cases import generate_initial_conditions

    # forcing_cases()
    # one_case()
    # proposal_figures()
    # forcing_sweep()
    # forcing_cases()
    # full_system_test()
    # figure2()
    track_inviscid_compensation_local()
    # generate_initial_conditions()
    pass
