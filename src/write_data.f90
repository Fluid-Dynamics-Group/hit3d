module m_data
    real*8, dimension(:, :, :), allocatable :: dUdX, dUdY, dUdZ, dVdX, dVdY, dVdZ, dWdX, dWdY, dWdZ
    real*8, dimension(:, :, :), allocatable :: dPxdX, dPxdY, dPxdZ, dPydX, dPydY, dPydZ, dPzdX, dPzdY, dPzdZ
    real*8, dimension(:, :, :), allocatable :: OmgX, OmgY, OmgZ
    ! second order derivatives for third term in dE/dt
    real*8, dimension(:, :, :), allocatable :: dUdX2, dUdY2, dUdZ2, dVdX2, dVdY2, dVdZ2, dWdX2, dWdY2, dWdZ2
    ! curl ( RHS) terms
    real*8, dimension(:, :, :), allocatable :: dRHS1dX, dRHS1dY, dRHS1dZ, dRHS2dX, dRHS2dY, dRHS2dZ, dRHS3dX, dRHS3dY, dRHS3dZ
    real*8, dimension(:,:,:), allocatable :: scalars_global 
    real*8, dimension(:,:,:,:), allocatable :: wrk_global

    ! io variables
    real*8 :: energy, helicity, solver_energy, solver_helicity, udotw, umag, wmag
    real*8 :: fcomp_u_left, fcomp_u_right
    real*8 :: fcomp_omega_left, fcomp_omega_right
end module m_data

subroutine write_scalars(current_timestep)
    use m_work
    use m_parameters
    use x_fftw
    implicit none

    integer :: current_timestep, n

    if (int_scalars) then
        do n = 1, n_scalars
            wrk(1:nx, 1:ny, 1:nz, 3+n) = fields(1:nx, 1:ny, 1:nz, 3+n)
            call xFFT3d(-1, 3 + n)
            
            ! write all the scalar data to an output file
            call send_scalars(n, current_timestep)
        end do
    end if

end subroutine write_scalars

! collect all of the scalar information for a single scalar at the 
! current timestep and write it to a csv file
subroutine send_scalars(wrk_idx, current_timestep)
    use m_work
    use m_parameters
    use m_data

    implicit none

    integer :: wrk_idx, current_timestep
    integer :: i, j, k
    character(len=80) :: filename
    real*8 :: scalar_value

    count = nx * ny * nz

    if (myid .ne. master) then
        id_to = master
        tag = myid
        call MPI_SEND(wrk(1:nx, 1:ny, 1:nz, wrk_idx), count, MPI_REAL8, master, tag, MPI_COMM_TASK, mpi_err)
    else
        ! allocate a global array - since we cant run any of the other modes at the same time as a write-mode
        ! we can be sure that this array has not already been allocated

        ! copy the master proc's wrk array in
        scalars_global(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, wrk_idx)

        do id_from = 1, numprocs - 1
            tag = id_from
            call MPI_RECV(wrk(1:nx, 1:ny, 1:nz, wrk_idx), count, MPI_REAL8, id_from, tag, MPI_COMM_TASK, mpi_status, mpi_err)

            scalars_global(1:nx, 1:ny, id_from*nz + 1:(id_from + 1)*nz) = wrk(1:nx, 1:ny, 1:nz, wrk_idx)
        end do

        ! we now have a giant array of all the data points from each mpi process, we write them to a file now
        write (filename, "('output/scalars/sc',I2.2,'.',I6.6,'.csv')") wrk_idx, current_timestep

        open (1015, file=filename, status="new")
        write(1015, "('scalar')")

        do k= 1,nz_all
            do j = 1,ny
                do i = 1,ny
                    scalar_value = scalars_global(i,j,k)
                    write(1015, "(E16.10)") scalar_value
                    
                end do
            end do
        end do


        call flush (1015)
        close(1015)

        ! in a previous
    end if
end subroutine send_scalars

subroutine write_velocity_field(current_timestep)
    use m_work ! wrk()
    use m_data
    use m_parameters ! nx, ny, nz, pertamp1, pertamp2

    implicit none

    integer :: filenumber, meta_file, i, j, k, current_timestep
    !real*8, allocatable :: fields(:,:,:,:)
    real*8 :: u, v, w, omgx_, omgy_, omgz_, omg_mag
    real*8 :: fcomp_left, fcomp_right, fcomp_total, epsilon_1, epsilon_2

    character(len=60) :: filename, sizefile

    if (load_initial_condition == 1) then
        epsilon_1 = 1.0
        epsilon_2 = 1.0
    else
        epsilon_1 = pertamp1
        epsilon_2 = pertamp2
    end if

    filenumber = 619

    call send_wrk_global()

    if (myid == master) then 
        write (filename, "('output/velocity_field/', i5.5, '.csv')") current_timestep
        open (filenumber, file=filename, status="new")
        write (filenumber, "('u,v,w,omega_mag,forcing')") 

        do k = 1, nz_all
            do j = 1, ny
                do i = 1, nx
                    u = wrk_global(i, j, k, 1)
                    v = wrk_global(i, j, k, 2)
                    w = wrk_global(i, j, k, 3)

                    omgx_ = wrk_global(i, j, k, 4)
                    omgy_ = wrk_global(i, j, k, 5)
                    omgz_ = wrk_global(i, j, k, 6)

                    omg_mag = omgx_**2 + omgy_**2 + omgz_**2

                    udotw = (u*omgx_) + (v*omgy_) + (w*omgz_)
                    umag = u**2 + v**2 + w**2
                    wmag = omgx_**2 + omgy_**2 + omgz_**2

                    fcomp_left = epsilon_1*(udotw*omgx_ - wmag*u) + &
                              epsilon_1*(udotw*omgy_ - wmag*v) + &
                              epsilon_1*(udotw*omgz_ - wmag*w)

                    fcomp_right = epsilon_2*(udotw*u - umag*omgx_) + &
                              epsilon_2*(udotw*v - umag*omgy_) + &
                              epsilon_2*(udotw*w - umag*omgz_)

                    fcomp_total = fcomp_left + fcomp_right

                    write (filenumber, "(E16.10, ',', E16.10, ',', E16.10 ',', E16.10, ',' E16.10)") &
                        u, v, w, omg_mag, fcomp_total
                end do
            end do
        end do

        flush (filenumber)

        meta_file = 620

        !write(sizefile, "('output/size_', i2.2, '.csv')") myid_world

        open (meta_file, file="output/size.csv")

        write (meta_file, "(A8)") "nx,ny,nz"
        write (meta_file, "(I5.5, A1, I5.5, A1, I5.5)") nx, ",", ny, ",", nz_all
    end if

end subroutine write_velocity_field

subroutine send_wrk_global
    use m_parameters
    use m_data
    use m_work

    implicit none

    count = nx * ny * nz * 6

    if (myid .ne. master) then
        id_to = master
        tag = myid
        call MPI_SEND(wrk(1:nx, 1:ny, 1:nz, 1:6), count, MPI_REAL8, master, tag, MPI_COMM_TASK, mpi_err)
    else
        ! allocate a global array - since we cant run any of the other modes at the same time as a write-mode
        ! we can be sure that this array has not already been allocated

        ! copy the master proc's wrk array in
        wrk_global(1:nx, 1:ny, 1:nz, 1:6) = wrk(1:nx, 1:ny, 1:nz, 1:6)

        do id_from = 1, numprocs - 1
            tag = id_from
            call MPI_RECV(wrk(1:nx, 1:ny, 1:nz, 1:6), count, MPI_REAL8, id_from, tag, MPI_COMM_TASK, mpi_status, mpi_err)

            wrk_global(1:nx, 1:ny, id_from*nz + 1:(id_from + 1)*nz, 1:6) = wrk(1:nx, 1:ny, 1:nz, 1:6)
        end do

    end if

end subroutine send_wrk_global

! create the file name and open the CSV files that we are using to write E(t) and h(t)
subroutine init_write_energy
    use m_parameters
    use m_data
    implicit none

    integer:: filenumber
    character(len=40) filename

    filenumber = 621
    call create_energy_filename(filename)

    if (master == myid) then
        open (filenumber, file="output/energy.csv")
        write (filenumber, "('current_time,', 'energy,', 'solver_energy,', 'helicity,', 'solver_helicity,', &
              'fdot_u_1,', 'fdot_u_2,', 'fdot_omega_1,', 'fdot_omega_2,', 'f_rate_e,', 'f_rate_h,', &
              're_lambda' &
          )")
    end if

    filenumber = 622
    open (filenumber, file="output/velocity.csv")
    write (filenumber, "('u,v,w')")

    allocate (dUdX(nx, ny, nz)); allocate (dUdY(nx, ny, nz)); allocate (dUdZ(nx, ny, nz));
    allocate (dVdX(nx, ny, nz)); allocate (dVdY(nx, ny, nz)); allocate (dVdZ(nx, ny, nz));
    allocate (dWdX(nx, ny, nz)); allocate (dWdY(nx, ny, nz)); allocate (dWdZ(nx, ny, nz));
    allocate (OmgX(nx, ny, nz)); allocate (OmgY(nx, ny, nz)); allocate (OmgZ(nx, ny, nz));
    ! second order vars for Del^2(u)
    allocate (dUdX2(nx, ny, nz)); allocate (dUdY2(nx, ny, nz)); allocate (dUdZ2(nx, ny, nz));
    allocate (dVdX2(nx, ny, nz)); allocate (dVdY2(nx, ny, nz)); allocate (dVdZ2(nx, ny, nz));
    allocate (dWdX2(nx, ny, nz)); allocate (dWdY2(nx, ny, nz)); allocate (dWdZ2(nx, ny, nz));
    ! pressure terms
    allocate (dPxdX(nx, ny, nz)); allocate (dPxdY(nx, ny, nz)); allocate (dPxdZ(nx, ny, nz)); 
    allocate (dPydX(nx, ny, nz)); allocate (dPydY(nx, ny, nz)); allocate (dPydZ(nx, ny, nz)); 
    allocate (dPzdX(nx, ny, nz)); allocate (dPzdY(nx, ny, nz)); allocate (dPzdZ(nx, ny, nz)); 
    ! RHS curl terms
    allocate (dRHS1dX(nx, ny, nz)); allocate (dRHS1dY(nx, ny, nz)); allocate (dRHS1dZ(nx, ny, nz)); 
    allocate (dRHS2dX(nx, ny, nz)); allocate (dRHS2dY(nx, ny, nz)); allocate (dRHS2dZ(nx, ny, nz)); 
    allocate (dRHS3dX(nx, ny, nz)); allocate (dRHS3dY(nx, ny, nz)); allocate (dRHS3dZ(nx, ny, nz)); 

    ! for writing scalar stuff
    allocate (scalars_global(nx, ny, nz_all))
    allocate (wrk_global(nx, ny, nz_all, 6))
end subroutine init_write_energy

! initalize the name of the file that we are writing to with E(t) and h(t) values
subroutine create_energy_filename(filename)
    use m_parameters
    character(len=40) :: filename
    write (filename, "('output/energy/energy_', i2.2, '_.csv')") myid_world
end subroutine create_energy_filename

! calculate and write the energy and helicity values for the current time step
! you must ensure that you call init_write_energy to open the correct files
! before calling this subroutine
subroutine write_energy(current_time)
    use m_work ! wrk
    use m_parameters
    use x_fftw
    use m_data
    implicit none

    integer :: filenumber
    real*8 :: current_time
    ! tmp variables for calculations
    real*8 :: u, v, w, u_rhs, v_rhs, w_rhs, omg_x, omg_y, omg_z
    real*8 :: f_rate, f_rate_e, f_rate_h
    real*8 :: energy_derivative, helicity_derivative
    ! loop variables
    integer :: i, j, k

    real*8 :: epsilon_1, epsilon_2
    real*8 :: tmp_val, frac
    real*8 :: re_lambda


    ! ifft the RHS variables so we are in x-space
    call ifft_rhs

    ! load variables from fields into wrk
    ! wrk(1:nx,1:ny,1:nz,1) = fields(1:nx,1:ny,1:nz,1)
    ! wrk(1:nx,1:ny,1:nz,2) = fields(1:nx,1:ny,1:nz,2)
    ! wrk(1:nx,1:ny,1:nz,3) = fields(1:nx,1:ny,1:nz,3)

    wrk(:, :, :, 1) = fields(:, :, :, 1)
    wrk(:, :, :, 2) = fields(:, :, :, 2)
    wrk(:, :, :, 3) = fields(:, :, :, 3)

    ! copy the velocity data into wrk (1,2,3)
    tmp_wrk(:, :, :, 1:3) = wrk(:, :, :, 1:3)

    !
    ! calculate derivatives of wrk varaibles to get omegas
    ! into (4,5,6)
    !

    call x_derivative(3, 'y', 6) ! dw/dy
    call x_derivative(3, 'x', 5) ! dw/dx

    call x_derivative(2, 'z', 4) ! dv/dz
    call x_derivative(2, 'x', 3) ! dv/dx

    call x_derivative(1, 'z', 2) ! du/dz
    call x_derivative(1, 'y', 1) ! du/dy

    wrk(:, :, :, 4) = wrk(:, :, :, 6) - wrk(:, :, :, 4)  ! omega_1 = w_y - v_z
    wrk(:, :, :, 5) = wrk(:, :, :, 2) - wrk(:, :, :, 5)  ! omega_2 = u_z - w_x
    wrk(:, :, :, 6) = wrk(:, :, :, 3) - wrk(:, :, :, 1)  ! omega_3 = v_x - u_y

    ! copy velocities back into wrk (1,2,3)
    wrk(:, :, :, 1:3) = tmp_wrk(:, :, :, 1:3)

    ! truncate and convert data to x-space
    do i = 1, 6
        !call truncate_and_inverse_wrk_idx(i)
        call xFFT3d(-1,i)
    end do

    ! wrk( 1-3) contains velocity information
    ! wrk( 4-6) contains omega information

    !
    ! Main calculation loop - perform required integrals
    !

    energy = 0.
    helicity = 0.
    solver_energy = 0.
    solver_helicity = 0.
    fcomp_omega_left = 0.
    fcomp_omega_right = 0.
    fcomp_u_left = 0.
    fcomp_u_right = 0.

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u = wrk(i, j, k, 1)
                v = wrk(i, j, k, 2)
                w = wrk(i, j, k, 3)

                omg_x = wrk(i, j, k, 4)
                omg_y = wrk(i, j, k, 5)
                omg_z = wrk(i, j, k, 6)

                u_rhs = rhs_saved(i, j, k, 1)
                v_rhs = rhs_saved(i, j, k, 2)
                w_rhs = rhs_saved(i, j, k, 3)

                ! define some common variables through the calculations
                udotw = (u*omg_x) + (v*omg_y) + (w*omg_z)
                umag = u**2 + v**2 + w**2
                wmag = omg_x**2 + omg_y**2 + omg_z**2

                ! start calculating variables that are ouput to IO

                helicity = helicity + udotw
                solver_helicity = solver_helicity + ((u_rhs*omg_x) + (v_rhs*omg_y) + (w_rhs*omg_z))

                energy = energy + umag
                solver_energy = solver_energy + (u_rhs*u + v_rhs*v + w_rhs*w)

                ! calculate forcing stuff
                fcomp_u_left = fcomp_u_left + &
                               u*(udotw*omg_x - wmag*u) + &
                               v*(udotw*omg_y - wmag*v) + &
                               w*(udotw*omg_z - wmag*w)

                fcomp_u_right = fcomp_u_right + &
                                u*(udotw*u - umag*omg_x) + &
                                v*(udotw*v - umag*omg_y) + &
                                w*(udotw*w - umag*omg_z)

                fcomp_omega_left = fcomp_omega_left + &
                                   omg_x*(udotw*omg_x - wmag*u) + &
                                   omg_y*(udotw*omg_y - wmag*v) + &
                                   omg_z*(udotw*omg_z - wmag*w)

                fcomp_omega_right = fcomp_omega_right + &
                                    omg_x*(udotw*u - umag*omg_x) + &
                                    omg_y*(udotw*v - umag*omg_y) + &
                                    omg_z*(udotw*w - umag*omg_z)

                f_rate = f_rate + &
                    (v *omg_z - w *omg_y)**2 + &
                    (u *omg_z - w *omg_x)**2 + &
                    (u *omg_y - v *omg_x)**2
            end do
        end do
        !write(*,*) "energy", energy
    end do

    frac = (2*3.1415)**3/(nx*ny*nz_all)

    helicity = helicity*frac/2
    solver_helicity = solver_helicity*frac
    energy = energy*frac/2.
    solver_energy = solver_energy*frac

    fcomp_u_left = fcomp_u_left * frac
    fcomp_u_right = fcomp_u_right * frac
    fcomp_omega_left = fcomp_omega_left * frac
    fcomp_omega_right = fcomp_omega_right * frac

    ! we initialize these values for the rate specifically
    epsilon_1 = PERTamp1
    epsilon_2 = PERTamp2

    f_rate_e = (-1* f_rate * epsilon_1 * frac)
    f_rate_h = (-1* f_rate * epsilon_2 * frac)

    !
    ! check that all of the variables are not NAN
    !

    call error_on_nan(helicity, "helicity")
    call error_on_nan(solver_helicity, "solver helicity")
    call error_on_nan(energy, "energy")
    call error_on_nan(solver_energy, "solver energy")

    ! check the nan on the forcing terms
    call error_on_nan(fcomp_u_left, "f_u 1")
    call error_on_nan(fcomp_u_right, "f_u 2")
    call error_on_nan(fcomp_omega_left, "f_omega 1")
    call error_on_nan(fcomp_omega_right, "f_omega 2")

    call error_on_nan(f_rate_e, "f_rate_e")
    call error_on_nan(f_rate_h, "f_rate_h")

    !
    ! sum the values through mpi
    !

    call add_through_mpi(helicity)
    call add_through_mpi(solver_helicity)
    call add_through_mpi(energy)
    call add_through_mpi(solver_energy)

    ! add the forcing terms
    call add_through_mpi(fcomp_u_left)
    call add_through_mpi(fcomp_u_right)
    call add_through_mpi(fcomp_omega_left)
    call add_through_mpi(fcomp_omega_right)

    call add_through_mpi(f_rate_e)
    call add_through_mpi(f_rate_h)

    ! calculate Re_lambda (taylor reynolds number) according to MGM's 
    ! formulation in m_stats.f90 subroutine stat_velocity
    call calculate_re_lambda(re_lambda)

    !
    ! write the summary data to file if we are the master process
    !

    if (myid == master) then
        filenumber = 621
        ! initialize the name of the csv that this mpi process will write to
        open (filenumber, file="output/energy.csv")
        write (filenumber, "(E16.10, ',', E16.10, ',', E16.10, ',', &
  &            E16.10, ',', E16.10, ',', E16.10, ',', E16.10, ',' E16.10, ',', E16.10, ',' &
               E16.10, ',', E16.10, ',', E16.10 &
              )") &
            current_time, energy, solver_energy, helicity, solver_helicity, fcomp_u_left, &
            fcomp_u_right, fcomp_omega_left, fcomp_omega_right, f_rate_e, f_rate_h, re_lambda

        flush (filenumber)
    end if

end subroutine

subroutine calculate_re_lambda(re_lambda_us)
    use m_stats

    implicit none

    real * 8 :: re_lambda_us

    call stat_velocity

    re_lambda_us = re_lambda
end 

subroutine ifft_rhs
    use m_work !tmp_wrk + wrk + rhs_saved
    use x_fftw ! fft stuff
    use m_parameters ! kmax

    ! save the current wrk array
    tmp_wrk(:, :, :, 1:3) = wrk(:, :, :, 1:3)

    ! copy the RHS variables into wrk
    wrk(:, :, :, 1:3) = rhs_saved(:, :, :, 1:3)

    ! truncate all the variables + perform ifft
    do i = 1, 3
        !call truncate_and_inverse_wrk_idx(i)
        call xFFT3d(-1, i)
    end do

    ! copy the data out of wrk and into the rhs
    rhs_saved(:, :, :, 1:3) = wrk(:, :, :, 1:3)

    ! copy the original data back into wrk
    wrk(:, :, :, 1:3) = tmp_wrk(:, :, :, 1:3)
end subroutine ifft_rhs

! runs through an index of the wrk variable and truncates the values
! before inverting the array to x-space
subroutine truncate_and_inverse_wrk_idx(idx)
    use m_work
    use m_parameters
    use x_fftw
    implicit none

    integer :: idx, i, j, k
    real*8  :: wmag2, rkmax2

    rkmax2 = real(kmax, 8)**2
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx + 2
                ! magnitude of the wave number
                wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2

                if (wmag2 .gt. rkmax2) then
                    wrk(i, j, k, idx) = zip
                else
                    ! do nothing, we already have wrk loaded in
                end if

            end do
        end do
    end do

    call xFFT3d(-1, idx)
end subroutine truncate_and_inverse_wrk_idx

! sum `variable_to_add` across every mpi process
! and return the result
subroutine add_through_mpi(variable_to_add)
    use m_parameters
    implicit none
    real*8 :: variable_to_add, tmp_val

    tmp_val = variable_to_add
    count = 1
    call MPI_REDUCE(tmp_val, variable_to_add, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
end subroutine add_through_mpi

subroutine error_on_nan(variable_to_check, variable_name)
    use m_openmpi
    use m_io
    implicit none
    real*8 :: variable_to_check
    character(len=*) :: variable_name

    if (variable_to_check /= variable_to_check) then
        write(out, *) variable_name, " was NAN - killing the simulation"
        write(*, *) variable_name, " was NAN - killing the simulation"
        call my_exit(1)
        call m_openmpi_exit
    end if
end subroutine error_on_nan

subroutine write_slice(current_timestep)
    use m_work
    use m_parameters
    implicit none

    integer :: current_timestep, i,j, filenumber,k, desired_proc
    character(len=200) :: filename
    real*8 :: u, v, w, omgx, omgy, omgz, omg_mag
    real*8 :: fcomp_left, fcomp_right, fcomp_total, epsilon_1, epsilon_2
    real*8 :: umag, udotw, wmag


    filenumber = 995
    k = nz

    ! the mpi process number that will be responsible
    ! for output in this subroutine
    desired_proc = MAX(INT(CEILING(REAL(numprocs)/2.0)) - 1, 0)

    if (myid == desired_proc) then
        if (load_initial_condition == 1) then
            epsilon_1 = 1.0
            epsilon_2 = 1.0
        else
            epsilon_1 = pertamp1
            epsilon_2 = pertamp2
        end if

        write(filename, "('output/slice/', i5.5, '.csv')") current_timestep
        open(filenumber, file=filename, status="new")
        write(filenumber, "('u,v,w,omega_mag,forcing')")

        do j = 1,ny
            do i = 1,ny
                u = wrk(i, j, k, 1)
                v = wrk(i, j, k, 2)
                w = wrk(i, j, k, 3)

                omgx = wrk(i, j, k, 4)
                omgy = wrk(i, j, k, 5)
                omgz = wrk(i, j, k, 6)

                omg_mag = omgx**2 + omgy**2 + omgz**2

                udotw = (u*omgx) + (v*omgy) + (w*omgz)
                umag = u**2 + v**2 + w**2
                wmag = omgx**2 + omgy**2 + omgz**2

                fcomp_left = epsilon_1*(udotw*omgx - wmag*u) + &
                          epsilon_1*(udotw*omgy - wmag*v) + &
                          epsilon_1*(udotw*omgz - wmag*w)

                fcomp_right = epsilon_2*(udotw*u - umag*omgx) + &
                          epsilon_2*(udotw*v - umag*omgy) + &
                          epsilon_2*(udotw*w - umag*omgz)

                fcomp_total = fcomp_left + fcomp_right
                write (filenumber, "(E16.10, ',', E16.10, ',', E16.10 ',', E16.10, ',' E16.10)") &
                    u, v, w, omg_mag, fcomp_total

                
            end do
        end do

    end if
end subroutine write_slice
