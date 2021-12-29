import os
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
    run_shell_command("make")
    forcing_folder = f"epsilon_parameter_sweep_4"
    save_json_folder = f"{BASE_SAVE}/{forcing_folder}"

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for _ in os.listdir(save_json_folder):
        raise ValueError("directory is not empty") 
        #os.remove(os.path.join(save_json_folder, f))

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

    for epsilon_value in [1,2]:
        for i in range(0,7):

            if epsilon_value == 1:
                epsilon = pow(10, i-2)
                newcase = [epsilon, 0, f"ep1-{i}"]
            else:
                epsilon = pow(10, i-2)
                newcase = [0, epsilon, f"ep2-{i}"]

            cases.append(newcase)

    if IS_SINGULARITY and IS_DISTRIBUTED:
        output_folder = f"/distribute_save/"
    else:
        output_folder = f"../../distribute_save/"

    for skip_diffusion in [0,1]:

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
        case =  RunCase(
            skip_diffusion=0, 
            size=size, 
            dt=dt, 
            steps=steps, 
            restarts=0, 
            reynolds_number=re, 
            path=output_folder,
            load_initial_data=1, 
            restart_time=1.
        )

        case.write_to_json(batch_name, save_json_folder)

    if not os.path.exists(save_json_folder):
        os.mkdir(save_json_folder)

    for f in os.listdir(save_json_folder):
        os.remove(os.path.join(save_json_folder, f))

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

# in order to calculate forcing cases we need to have an initial condition file
def forcing_cases():
    delta_1 = .1
    delta_2 = .1

    run_shell_command("make")
    batch_name = f"forcing_viscous_5"
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
        #[0., 0., "baseline"],
        [-1*delta_1, 0., "ep1-neg"],
        [ 0., -1*delta_2, "ep2-neg"],
    ]

    output_folder = define_output_folder()

    for skip_diffusion in [0,1]:

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
                export_vtk=save_vtk,
                #scalar_type=14,
                scalar_type=0,
                io_steps=io_steps
            )

            case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)

# designed to easily tell if the full codebase is functioning as it should be
def full_system_test():
    delta_1 = .1
    delta_2 = .1

    run_shell_command("make")
    batch_name = f"system_test_14"
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
        #[0., 0., "baseline"],
        [-1*delta_1, 0., "ep1-neg"],
        [ 0., -1*delta_2, "ep2-neg"],
    ]

    output_folder = define_output_folder()

    for delta_1, delta_2, folder in cases:
        diffusion_str = skip_diffusion_to_str(0)

        epsilon1 = epsilon_generator.epsilon_1(delta_1)
        epsilon2 = epsilon_generator.epsilon_2(delta_2)

        if epsilon1 == -0.0:
            epsilon1 = 0.0
        if epsilon2 == -0.0:
            epsilon2 = 0.0

        case =  RunCase(skip_diffusion=0, 
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
            #scalar_type=0,
            io_steps=io_steps
        )

        case.write_to_json(f"{folder}_{diffusion_str}", save_json_folder)

    copy_distribute_files(save_json_folder, batch_name, extra_caps)

    build = Build("master", "master")
    build.to_json(save_json_folder)
