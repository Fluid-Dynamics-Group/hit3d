from run import run_shell_command
from run import BASE_SAVE
from run import RunCase
from run import skip_diffusion_to_str
from run import wrap_error_case

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
