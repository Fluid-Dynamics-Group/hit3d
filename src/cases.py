import os
import shutil
from pprint import pprint
from run import run_shell_command
from run import define_output_folder
from run import BASE_SAVE
from run import IS_DISTRIBUTED
from run import IS_SINGULARITY
from run import UNR
from run import skip_diffusion_to_str
from run import copy_init_files
from run import copy_distribute_files
from run import EpsilonControl
from run import RunCase
from run import Build


def forcing_sweep():
    # Conclusion:
    # for viscous cases, the max I values are
    # ep1 : 1 -> 10^(i-2) -> 0.1
    # ep2 : 3 -> 10^(i-2) -> 10

    run_shell_command("make")
    n = 5
    forcing_folder = f"epsilon_parameter_sweep_{n}"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for _ in os.listdir(save_json_folder):
        raise ValueError("directory is not empty")
        # os.remove(os.path.join(save_json_folder, f))

    END_TIME = 10
    dt = 0.0005
    size = 128
    re = 40
    steps = int(END_TIME / dt)
    save_vtk = False
    batch_name = forcing_folder
    extra_caps = []

    copy_init_files(size)

    epsilon_generator = EpsilonControl.load_json()

    cases = []

    NUM_CASES = 10

    delta_1_max = 0.5
    delta_1_min = 0.01
    delta_1_inc = (delta_1_max - delta_1_min) / (NUM_CASES - 1)

    delta_2_max = 100
    delta_2_min = 0.1
    delta_2_inc = (delta_2_max - delta_2_min) / (NUM_CASES - 1)

    for epsilon_value in [1, 2]:
        for i in range(0, NUM_CASES):

            if epsilon_value == 1:
                delta = delta_1_min + (delta_1_inc * i)
                newcase = [delta, 0, f"ep1-{i}"]
            else:
                delta = delta_2_min + (delta_2_inc * i)
                newcase = [0, delta, f"ep2-{i}"]

            cases.append(newcase)

    delta_1_cases = [i for i, _, _ in cases[0:NUM_CASES]]
    epsilon_1_cases = [epsilon_generator.epsilon_1(i) for i in delta_1_cases]
    z1 = list(zip(delta_1_cases, epsilon_1_cases))

    delta_2_cases = [i for _, i, _ in cases[NUM_CASES:]]
    epsilon_2_cases = [epsilon_generator.epsilon_2(i) for i in delta_2_cases]
    z2 = list(zip(delta_2_cases, epsilon_2_cases))

    print("asdf")

    print(f"forced energy cases:")
    pprint(z1)
    print(f"forced helicity cases:")
    pprint(z2)

    if IS_SINGULARITY and IS_DISTRIBUTED:
        output_folder = f"/distribute_save/"
    else:
        output_folder = f"../../distribute_save/"

    for skip_diffusion in [0, 1]:

        # TODO: remove this after generating
        if skip_diffusion == 1:
            continue

        for delta_1, delta_2, folder in cases:
            diffusion_str = skip_diffusion_to_str(skip_diffusion)
            epsilon1 = epsilon_generator.epsilon_1(delta_1)
            epsilon2 = epsilon_generator.epsilon_2(delta_2)

            if epsilon1 == -0.0:
                epsilon1 = 0.0
            if epsilon2 == -0.0:
                epsilon2 = 0.0

            print("ep1: {:.5E} ep2: {:.5E}".format(epsilon1, epsilon2))

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


def generate_initial_conditions():
    re = 40

    n = 3
    batch_name = f"initial_conditions_{n}"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    output_folder = define_output_folder()
    extra_caps = []

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)
    else: 
        raise ValueError(f"json output folder {save_json_folder} already exists. Stopping here")

    input_cases = [
        [0.0005, 64, 5000],
        [0.0005, 128, 5000],
        [0.0001, 256, 25000],
    ]

    for dt, size, steps in input_cases:
        case = RunCase(
            skip_diffusion=0,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            reynolds_number=re,
            path=output_folder,
            load_initial_data=1,
        )

        case.write_to_json(f"ic_{case.size}", save_json_folder)

    build = Build("master", "master")
    build.to_json(save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps, False)


# in order to calculate forcing cases we need to have an initial condition file
def forcing_cases():
    delta_1 = 0.1
    delta_2 = 0.1

    run_shell_command("make")
    batch_name = f"forcing_viscous_7"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    END_TIME = 5
    dt = 0.0005
    size = 128
    re = 40
    steps = int(END_TIME / dt)
    save_vtk = True
    extra_caps = []

    io_steps = 38

    copy_init_files(size)

    epsilon_generator = EpsilonControl.load_json()

    cases = [
        [0.0, 0.0, "baseline"],
        [-1 * delta_1, 0.0, "ep1-neg"],
        [0.0, -1 * delta_2, "ep2-neg"],
    ]

    output_folder = define_output_folder()

    for skip_diffusion in [0, 1]:

        # TODO: remove this after generating
        if skip_diffusion == 1:
            continue

        for delta_1, delta_2, folder in cases:
            diffusion_str = skip_diffusion_to_str(skip_diffusion)
            epsilon1 = epsilon_generator.epsilon_1(delta_1)
            epsilon2 = epsilon_generator.epsilon_2(delta_2)

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
                # scalar_type=14,
                scalar_type=0,
                io_steps=io_steps,
            )

            case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)


# designed to easily tell if the full codebase is functioning as it should be
def full_system_test():
    delta_1 = 0.1
    delta_2 = 0.1

    n = 20
    run_shell_command("make")
    batch_name = f"system_test_{n}"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for _ in os.listdir(save_json_folder):
        raise ValueError("folder should be empty")

    dt = 0.0005
    size = 128
    re = 40
    steps = 500
    save_vtk = True
    extra_caps = []

    io_steps = 10

    copy_init_files(size)

    epsilon_generator = EpsilonControl.load_json()

    cases = [
        # [0., 0., "baseline"],
        [-1 * delta_1, 0.0, "ep1-neg"],
        [0.0, -1 * delta_2, "ep2-neg"],
    ]

    output_folder = define_output_folder()

    for delta_1, delta_2, folder in cases:
        diffusion_str = skip_diffusion_to_str(0)

        epsilon1 = epsilon_generator.epsilon_1(delta_1)
        epsilon2 = epsilon_generator.epsilon_2(delta_2)

        case = RunCase(
            skip_diffusion=0,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            reynolds_number=re,
            path=output_folder,
            load_initial_data=0,
            epsilon1=epsilon1,
            epsilon2=epsilon2,
            export_vtk=save_vtk,
            scalar_type=14,
            # scalar_type=0,
            io_steps=io_steps,
        )

        case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)


# plots and data specifically for generating figure 2 of aditya's paper
def figure2():
    # n=1..7 had delta_2 as .1     : delta_1 = 0.1
    # n= 8..9 had delta_2 as 10.0  : delta_1 = 0.1
    # n= 10 had delta_2=1          : delta_1 = 0.1
    # n= 11..17 had delta_2=0.7          : delta_1 = 0.05
    delta_1 = -0.05
    delta_2 = -0.7

    include_viscous_compensation = True

    run_shell_command("make")
    n = 24

    if include_viscous_compensation:
        batch_name = f"figure2_viscous_comp_{n}"
    else:
        batch_name = f"figure2_{n}"

    save_json_folder = f"{BASE_SAVE}/{batch_name}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    END_TIME = 5.0
    dt = 0.0005
    size = 128
    re = 40
    steps = int(END_TIME / dt)
    save_vtk = True
    extra_caps = []

    io_steps = 38

    copy_init_files(size)

    epsilon_generator = EpsilonControl.load_json()

    if include_viscous_compensation:
        # load the data
        print("generating figure 2 data with viscous compensation being read")
        visc_comp_param = 1
    else:
        # dont do any viscous compensation 
        print("no viscous_compensation_in_effect")
        visc_comp_param = 0

    cases = [
        # baseline
        # n=7 has baseline
        # always generate a viscous compensation case when we are creating baseline information
        #[0.0, 0.0, "baseline_viscous", 0, 2],
        # viscous stuff
        [delta_1, 0.0, "energy_modification_viscous", 0, visc_comp_param],
        [0.0, 1 * delta_2, "helicity_modification_viscous", 0, visc_comp_param],
        [delta_1, 1 * delta_2, "both_modification_viscous", 0, visc_comp_param],
        # inviscid stuff
        # [delta_1, 0., "energy_modification_inviscid", 1, visc_comp_param],
        # [ 0., delta_2, "helicity_modification_inviscid", 1, visc_comp_param],
        # [delta_1, delta_2, "both_modification_inviscid", 1, visc_comp_param],
    ]

    output_folder = define_output_folder()

    for delta_1, delta_2, folder, skip_diffusion, viscous_compensation in cases:
        epsilon1 = epsilon_generator.epsilon_1(delta_1)
        epsilon2 = epsilon_generator.epsilon_2(delta_2)

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
            # scalar_type=14,
            scalar_type=0,
            io_steps=io_steps,
            viscous_compensation=viscous_compensation
        )

        case.write_to_json(folder, save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps, include_viscous_compensation, size)

    build = Build("visc-compensation-tracking", "master")
    build.to_json(save_json_folder)

# helpful function for runnning one-off cases
def one_case():
    TIME_END = 0.1
    batch_name = "viscous_baseline_2"
    job_name = "single-case"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    size = 128
    dt = 0.0005
    steps = int(TIME_END / dt)
    nprocs = 16
    extra_caps = []
    io_steps = int(steps / 10)
    load_initial_data = 0
    export_divergence = 0

    ep1 = 0.001
    ep2 = 0.001

    skip_diffusion = 0

    output_folder = define_output_folder()

    if not (load_initial_data == 2):
        copy_init_files(size)

    # if the directory exists remove any older files from the dir
    if os.path.exists(save_json_folder):
        for f in os.listdir(save_json_folder):
            os.remove(os.path.join(save_json_folder, f))

    os.makedirs(save_json_folder, exist_ok=True)

    run_shell_command("make")

    case = RunCase(
        skip_diffusion=skip_diffusion,
        size=size,
        dt=dt,
        steps=steps,
        restarts=0,
        reynolds_number=40,
        path=output_folder,
        load_initial_data=load_initial_data,
        epsilon1=ep1,
        # delta2 is negative, this will decrease helicity
        epsilon2=ep2,
        export_vtk=True,
        scalar_type=14,
        require_forcing=1,
        viscous_compensation=0,
        validate_viscous_compensation=0,
        io_steps=io_steps,
        nprocs=nprocs,
        export_divergence=export_divergence,
    )

    if IS_DISTRIBUTED:
        print("creating files to run on distributed compute")
        case.write_to_json(job_name, save_json_folder)

        copy_distribute_files(save_json_folder, batch_name, extra_caps, False, size)

        build = Build("master", "master")
        build.to_json(save_json_folder)

    else:
        print("running the case locally")
        # init files have already been copied above
        case.run(1)

# helpful function for runnning one-off cases
def track_inviscid_compensation_local():
    TIME_END = 5
    size = 128
    dt = 0.0005
    steps = int(TIME_END / dt)
    nprocs = 16
    load_initial_data = 0
    export_divergence = 0

    copy_init_files(size)
    epsilon_generator = EpsilonControl.load_json()
    delta_1 = -0.05
    delta_2 = -0.7
    ep1 = epsilon_generator.epsilon_1(delta_1)
    ep2 = epsilon_generator.epsilon_2(delta_2)

    output_folder = define_output_folder()

    run_shell_command("make")

    def make_case(skip_diffusion: int, viscous_compensation: int, ep1: float, ep2: float) -> RunCase:
        return RunCase(
            skip_diffusion=skip_diffusion,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            reynolds_number=40,
            path=output_folder,
            load_initial_data=load_initial_data,
            epsilon1=ep1,
            # delta2 is negative, this will decrease helicity
            epsilon2=ep2,
            export_vtk=False,
            scalar_type=14,
            require_forcing=1,
            viscous_compensation=viscous_compensation,
            validate_viscous_compensation=0,
            io_steps=10,
            nprocs=nprocs,
            export_divergence=export_divergence,
        )

    if IS_DISTRIBUTED == False:
        # first, an inviscid case

        print("running baseline case")
        baseline_case = make_case(0, 2, 0.0, 0.0)
        baseline_case.run(1)

        # move the binary files so we can execute them here
        os.rename("current_results/dE_dt_history.binary","./dE_dt_history.binary")
        os.rename("current_results/dh_dt_history.binary","./dh_dt_history.binary")

        shutil.move("current_results", "visc_baseline_results");
        
        print("running viscous case to track")
        print("WARNING: make sure the .binary files are in the current directory")
        viscous_case = make_case(0, 1, ep1, ep2)
        viscous_case.run(2)
    else:
        raise ValueError("local viscous compensation validation must be run locally")

def track_viscous_compensation_remote(is_inviscid: bool):
    # THESE SHOULD BE THE SAME FROM figure_2()
    delta_1 = -0.05
    delta_2 = -0.7


    END_TIME = 5.0
    dt = 0.0005
    size = 128
    re = 40
    steps = int(END_TIME / dt)
    extra_caps = []

    copy_init_files(size)
    epsilon_generator = EpsilonControl.load_json()
    ep1 = epsilon_generator.epsilon_1(delta_1)
    ep2 = epsilon_generator.epsilon_1(delta_2)

    n = 3

    if is_inviscid:
        batch_name = f"viscous_compensation_tracking_inviscid_{n}"
    else:
        batch_name = f"viscous_compensation_tracking_viscous_{n}"

    save_json_folder = f"{BASE_SAVE}/{batch_name}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)
    else: 
        raise ValueError(f"json output folder {save_json_folder} already exists. Stopping here")

    def make_case(skip_diffusion: int, viscous_compensation: int, ep1: float, ep2: float) -> RunCase:
        output_folder = define_output_folder()

        return RunCase(
            skip_diffusion=skip_diffusion,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            reynolds_number=re,
            path=output_folder,
            epsilon1=ep1,
            # delta2 is negative, this will decrease helicity
            epsilon2=ep2,
            require_forcing=1,
            viscous_compensation=viscous_compensation,
        )

    build = Build("visc-compensation-tracking", "master")
    build.to_json(save_json_folder)

    if IS_DISTRIBUTED:
        if is_inviscid:
            case = make_case(0, 2, ep1, ep2)
            case.write_to_json(batch_name, save_json_folder)
            copy_distribute_files(save_json_folder, batch_name, extra_caps, False)
        else:
            case = make_case(0, 1, ep1, ep2)
            case.write_to_json(batch_name, save_json_folder)
            copy_distribute_files(save_json_folder, batch_name, extra_caps, True)

    else:
        raise ValueError("remote run must be run with IS_DISTRIBUTED = True")
