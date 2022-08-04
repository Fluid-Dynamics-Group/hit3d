module write_data_m
    real*8 :: energy
    integer:: all_energy_handle
    integer:: after_turbulence_energy_handle
    integer:: u_handle, v_handle, w_handle
    logical :: wrote_velocity_field
    integer :: force_turbulence_step
end module

subroutine calculate_energy()
    use m_work
    !use m_fields
    use m_parameters
    use x_fftw
    use write_data_m
    implicit none

    real*8 :: u, v, w
    integer :: i, j, k
    real*8 :: frac

    !
    ! invert fields back to x-space
    !
    wrk(1:nx, 1:ny, 1:nz, 1:3) = fields(1:nx, 1:ny, 1:nz, 1:3)

    do i = 1, 3
        call xFFT3d(-1, i)
    end do

    ! calculate the total energy in the field for each mpi node
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u = wrk(i, j, k, 1)
                v = wrk(i, j, k, 2)
                w = wrk(i, j, k, 3)

                energy = energy + (u**2 + v**2 + w**2)
            enddo
        enddo
    enddo

    ! normalize the summation to the dx, dy, dz
    frac = (2*3.1415)**3/(nx*ny*nz_all)

    energy = (energy / 2.) * frac

    ! master node will have the the energy from all procs
    call add_through_mpi(energy)

    !if (iammaster) write(*,*) "energy is ", energy

end subroutine calculate_energy

subroutine add_through_mpi(variable_to_add)
    use m_parameters
    implicit none
    real*8 :: variable_to_add, tmp_val

    tmp_val = variable_to_add
    count = 1
    call MPI_REDUCE(tmp_val, variable_to_add, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
end subroutine add_through_mpi

subroutine open_files()
    use write_data_m
    use m_parameters
    implicit none

    all_energy_handle = 99995
    after_turbulence_energy_handle = 99996
    u_handle = 99997
    v_handle = 99998
    w_handle = 99999

    force_turbulence_step = 20000

    wrote_velocity_field = .false.

    if (iammaster) then
        open(all_energy_handle, file="output/all_energy.binary", status="new", access="stream")
        open(after_turbulence_energy_handle, file="output/after_turbulence_energy.binary", status="new", access="stream")

        open(u_handle, file="output/u_ic.binary", status="new", access="stream")
        open(v_handle, file="output/v_ic.binary", status="new", access="stream")
        open(w_handle, file="output/w_ic.binary", status="new", access="stream")
    endif

end subroutine open_files

subroutine close_files()
    use write_data_m
    implicit none

end subroutine close_files

subroutine write_energy(t, step)
    use write_data_m
    use m_parameters
    implicit none
    real*8 :: t;
    integer :: step

    if (iammaster) then
        write(all_energy_handle) energy
        write(all_energy_handle) t

        if (step > force_turbulence_step) then
            write (after_turbulence_energy_handle) energy
            write (after_turbulence_energy_handle) time 
        endif
    endif

end subroutine write_energy

subroutine write_velocity_field(t, step)
    use write_data_m
    use m_parameters
    implicit none

    ! the current time in the solver
    real*8 :: t
    integer :: step

    ! if the > condition applies, then we already started writing the energy trajectories out to the file,
    ! and so we will be writing an initial condition at a data point that we should not be!
    if ((step > force_turbulence_step) .and. (wrote_velocity_field .eqv. .false.)) then
        write(*,*) "failed to write the velocity field at the correct time, this should probably not happen!"
        call my_exit(-1)
    endif

    ! do the writing
    if (step == force_turbulence_step) then
        call write_field_to_file(1, u_handle)
        call write_field_to_file(2, v_handle)
        call write_field_to_file(3, w_handle)

        wrote_velocity_field = .true.
    endif

end subroutine 

subroutine write_field_to_file(velocity_comp, filehandle)
    use x_fftw
    use m_parameters
    implicit none
    integer :: velocity_comp, filehandle

    ! invert the required velocity field to the wrk index
    wrk(:,:,:,velocity_comp) = fields(:,:,:, velocity_comp)
    call xFFT3d(-1, velocity_comp)

    ! wrk(:,:,:, velocity_comp) now contains x-space velocity information that we need to write
    
    count = nx * ny * nz


    if (iammaster) then
        ! we are the master node, we are in charge of RX all arrays and writing them to a file

        ! first write our own buffer before receiving anything that would overwrite the data
        call write_wrk_buffer_to_file(velocity_comp, filehandle)

        ! read in all the buffers sent from the master proc
        do id_from = 1, numprocs - 1
            tag = id_from
            call MPI_RECV(wrk(1:nx, 1:ny, 1:nz, velocity_comp), count, MPI_REAL8, id_from, tag, MPI_COMM_TASK, mpi_status, mpi_err)

            call write_wrk_buffer_to_file(velocity_comp, filehandle)
        end do

    else
        ! we are a child node, we are in charge of  TX our velocity arrays
        id_to = master
        tag = myid
        call MPI_SEND(wrk(1:nx, 1:ny, 1:nz, velocity_comp), count, MPI_REAL8, master, tag, MPI_COMM_TASK, mpi_err)

    endif

end subroutine write_field_to_file

subroutine write_wrk_buffer_to_file(velocity_comp, filehandle)
    use m_parameters
    use m_work
    implicit none
    integer :: velocity_comp, filehandle
    integer :: i, j, k

    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                write(filehandle) wrk(i,j,k, velocity_comp)
            enddo
        enddo
    enddo

endsubroutine write_wrk_buffer_to_file
