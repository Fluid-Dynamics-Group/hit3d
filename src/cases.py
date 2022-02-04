import os
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


def initial_condition():
    re = 40

    batch_name = f"initial_conditions"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    output_folder = f"../../distribute_save"
    extra_caps = []

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
            restart_time=1.0,
        )

        case.write_to_json(batch_name, save_json_folder)

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    copy_distribute_files(save_json_folder, batch_name, extra_caps)


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
    delta_1 = 0.05
    delta_2 = 0.7

    run_shell_command("make")
    n = 20
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

    cases = [
        # baseline
        # n=7 has baseline
        [0.0, 0.0, "baseline_viscous", 0],
        # viscous stuff
        #[-1 * delta_1, 0.0, "energy_modification_viscous", 0],
        #[0.0, -1 * delta_2, "helicity_modification_viscous", 0],
        #[-1 * delta_1, -1 * delta_2, "both_modification_viscous", 0],
        # inviscid stuff
        # [-1*delta_1, 0., "energy_modification_inviscid", 1],
        # [ 0., -1*delta_2, "helicity_modification_inviscid", 1],
        # [-1*delta_1, -1*delta_2, "both_modification_inviscid", 1],
    ]

    output_folder = define_output_folder()

    for delta_1, delta_2, folder, skip_diffusion in cases:
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

        case.write_to_json(folder, save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)


# helpful function for runnning one-off cases
def generate_initial_conditions():
    TIME_END = 2.5
    batch_name = "initial_consitions_1"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    nprocs = 16
    extra_caps = []

    output_folder = define_output_folder()

    # if the directory exists remove any older files from the dir
    if os.path.exists(save_json_folder):
        for f in os.listdir(save_json_folder):
            os.remove(os.path.join(save_json_folder, f))

    os.makedirs(save_json_folder, exist_ok=True)

    run_shell_command("make")

    def make_case(dt: float, size: int):
        steps = int(TIME_END / dt)

        return RunCase(
            # initial conditions are calculated from inviscid data
            skip_diffusion=1,
            size=size,
            dt=dt,
            steps=steps,
            restarts=0,
            reynolds_number=40,
            path=output_folder,
            # write an initial condition from the default condition
            load_initial_data=1,
            epsilon1=0.0,
            # delta2 is negative, this will decrease helicity
            epsilon2=0.0,
            export_vtk=False,
            require_forcing=1,
            viscous_compensation=0,
            validate_viscous_compensation=0,
            nprocs=nprocs,
            export_divergence=0,
        )

    cases = [
        [64, 0.0005],
        [128, 0.0005],
        [256, 0.0001],
    ]

    cases = [make_case(dt, size) for size, dt in cases]

    if IS_DISTRIBUTED:
        print("creating files to run on distributed compute")
        for case in cases:
            job_name = f"ic_{case.size}"

            case.write_to_json(job_name, save_json_folder)

        copy_distribute_files(save_json_folder, batch_name, extra_caps)

        build = Build("master", "master")
        build.to_json(save_json_folder)

    else:
        raise ValueError("cannot generate an initial condition locally")

# helpful function for runnning one-off cases
def one_case():
    TIME_END = 0.5
    batch_name = "viscous_baseline_2"
    job_name = "single-case"
    save_json_folder = f"{BASE_SAVE}/{batch_name}"
    size = 128
    dt = 0.0005
    steps = int(TIME_END / dt)
    nprocs = 16
    extra_caps = []
    io_steps = None
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

        copy_distribute_files(save_json_folder, batch_name, extra_caps)

        build = Build("master", "master")
        build.to_json(save_json_folder)

    else:
        print("running the case locally")
        # init files have already been copied above
        case.run(1)

# helpful function for runnning one-off cases
def track_inviscid_compensation_local():
    TIME_END = 0.1
    size = 64
    dt = 0.0005
    steps = int(TIME_END / dt)
    nprocs = 16
    load_initial_data = 0
    export_divergence = 0

    epsilon_generator = EpsilonControl.load_json()
    ep1 = epsilon_generator.epsilon_1(0.01)
    ep2 = epsilon_generator.epsilon_1(0.01)

    output_folder = define_output_folder()

    copy_init_files(size)

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

        # print("running inviscid case")
        # inviscid_case = make_case(1, 2, ep1, ep2)
        # inviscid_case.run(1)
        
        print("running viscous case to track")
        print("WARNING: make sure the .binary files are in the current directory")
        viscous_case = make_case(0, 1, 0.0, 0.0)
        viscous_case.run(2)
    else:
        raise ValueError("local viscous compensation validation must be run locally")
