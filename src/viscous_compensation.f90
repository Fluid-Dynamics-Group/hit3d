module m_viscous_compensation 
    integer :: energy_filehandle
    integer :: helicity_filehandle
    real*8, allocatable :: dE_dt_tracking(:)
    real*8, allocatable :: dh_dt_tracking(:)
end module

! open the required files for viscous compensation data
subroutine init_viscous_compensation_files()
    use m_viscous_compensation
    use m_parameters

    implicit none

    integer :: i
    logical :: is_writing_files, is_reading_files

    ! set the integers for the files we will read and write to
    energy_filehandle = 990
    helicity_filehandle = 991

    call check_read_viscous_compensation(is_reading_files, .true.)
    call check_write_viscous_compensation(is_writing_files, .true.)


    if (is_writing_files) then
        if (iammaster) then 
            allocate(dE_dt_tracking(ITMAX))
            allocate(dh_dt_tracking(ITMAX))
        endif
    elseif (is_reading_files) then
        ! allocate the arrays that we will store the target values in
        allocate(dE_dt_tracking(ITMAX))
        allocate(dh_dt_tracking(ITMAX))

        open(energy_filehandle, file="dE_dt_history.binary", status="old", access="STREAM")
        open(helicity_filehandle, file="dh_dt_history.binary", status="old", access="STREAM")

        read(energy_filehandle) (dE_dt_tracking(i), i=1, ITMAX)
        read(helicity_filehandle) (dh_dt_tracking(i), i=1, ITMAX)
    endif

end subroutine

! finish up all viscous compensation related activities
subroutine finish_viscous_compensation_files()
    use m_viscous_compensation
    use m_parameters

    implicit none

    integer :: i
    logical :: is_writing_files, is_reading_files

    call check_read_viscous_compensation(is_reading_files, .false.)
    call check_write_viscous_compensation(is_writing_files, .false.)

    if (is_writing_files) then
        ! only the master node will have access to the arrays
        ! required for writing the binary files
        if (iammaster) then
            write(out, *) "writing binary files"
            write(*, *) "writing binary files"

            open(energy_filehandle, file="output/dE_dt_history.binary", status="new", access="STREAM")
            open(helicity_filehandle, file="output/dh_dt_history.binary", status="new", access="STREAM")

            ! now we write the contents of all the arrays to a file
            write(energy_filehandle) (dE_dt_tracking(i), i=1,ITMAX)
            write(helicity_filehandle) (dh_dt_tracking(i), i=1,ITMAX)

            write(out, *) "printing dh_dt data points"
            write(out, *) (dh_dt_tracking(i), i=1,max(ITMAX,10))

            flush(energy_filehandle)
            flush(helicity_filehandle)

            close(energy_filehandle)
            close(helicity_filehandle)
        endif

    elseif (is_reading_files) then
        ! nothing to do here
    endif

end subroutine

subroutine check_write_viscous_compensation(check, is_verbose)
    use m_parameters

    implicit none

    logical :: check
    logical :: is_verbose

    if (viscous_compensation == 2 .and. viscous_compensation_validation == 0) then
        if (is_verbose .and. iammaster) write(out, *) "we are writing viscous compensation files"
        if (is_verbose .and. iammaster) write(*, *) "we are writing viscous compensation files"
        check = .true.
    else
        if (is_verbose) write(out, *) "we are NOT writing viscous compensation files"
        check = .false.
    endif

end subroutine

! see if we are meant to read in viscous compensation files
subroutine check_read_viscous_compensation(check, is_verbose)
    use m_parameters

    implicit none

    logical :: check
    logical :: is_verbose

    if (viscous_compensation == 1 .and. viscous_compensation_validation == 0) then
        if (is_verbose .and. iammaster) write(out, *)"we are reading viscous compensation files"
        if (is_verbose .and. iammaster) write(*, *) "we are reading viscous compensation files"
        check = .true.
    else 
        if (is_verbose) write(out, *) "we are NOT reading viscous compensation files"
        check = .false.
    endif

end subroutine

! write to the arrays so they can be finalized later
! subroutine does nothing if we are not in write-mode
subroutine write_visc_comp_rates()
    use m_viscous_compensation
    use m_parameters
    use m_work
    use m_fields
    use x_fftw

    implicit none

    logical :: is_writing_files, is_reading_files

    real*8 :: dE_dt
    real*8 :: dh_dt
    real*8 :: frac ! scalaing factor for the energy at a given point in time

    ! inner loop variables for calculation
    real*8:: u_rhs, v_rhs, w_rhs
    real*8:: omg_x, omg_y, omg_z
    real*8:: u, v, w
    integer :: i, j, k

    call check_read_viscous_compensation(is_reading_files, .false.)
    call check_write_viscous_compensation(is_writing_files, .false.)

    if (is_writing_files) then 
        ! after this call
        ! wrk( 1-3) contains velocity information (x-space)
        ! wrk( 4-6) contains omega information (x-space)
        call calculate_vorticity_velocity_x_space()

        ! make sure RHS is in X-space
        call ifft_rhs

        ! initialize the WRK array to the velocity values from fields
        wrk(:,:,:,1:3) = fields(:,:,:,1:3)

        ! invert the wrk values to X-SPACE
        do i = 1, 3
            call xFFT3d(-1,i)
        end do

        !
        ! loop through the coordinates and calculate dE / dt at every time step
        !
        do i = 1, nx
            do j = 1, nx
                do k = 1, nz
                    
                    u_rhs = rhs_saved(i, j, k, 1)
                    v_rhs = rhs_saved(i, j, k, 2)
                    w_rhs = rhs_saved(i, j, k, 3)

                    u = wrk(i, j, k, 1)
                    v = wrk(i, j, k, 2)
                    w = wrk(i, j, k, 3)

                    omg_x = wrk(i, j, k, 4)
                    omg_y = wrk(i, j, k, 5)
                    omg_z = wrk(i, j, k, 6)

                    dE_dt = dE_dt+ (u_rhs*u + v_rhs*v + w_rhs*w)
                    dh_dt = dh_dt + ((u_rhs*omg_x) + (v_rhs*omg_y) + (w_rhs*omg_z))
                end do
            end do
        end do

        frac = (2*3.1415)**3/(nx*ny*nz_all)

        ! TODO: these might need to be / 2
        dE_dt = dE_dt * frac
        dh_dt = dh_dt * frac

        call error_on_nan(dE_dt, "dE_dt")
        call error_on_nan(dh_dt, "dh_dt")

        call add_through_mpi(dE_dt)
        call add_through_mpi(dh_dt)

        ! store the values in an array for later
        ! only the master node has these arrays allocated
        if (iammaster) then
            dE_dt_tracking(ITIME) = dE_dt
            dh_dt_tracking(ITIME) = dh_dt
        endif

        ! invert the RHS variables back to where they should be for future analysis
        call fft_rhs

        ! we dont really care what exists in the WRK array here
        ! since it will be overwritten later
    end if
end subroutine

! fetch the current derivatives that we neet to track
! this subroutine must only be called if you have validated that 
! we have vectors to read from
subroutine read_visc_dQ_dt_values(dq1, dq2)
    use m_parameters
    use m_viscous_compensation
    
    implicit none

    real*8 dq1, dq2

    dq1 = dE_dt_tracking(ITIME)
    dq2 = dh_dt_tracking(ITIME)

    !if (iammaster) then
    !    write(*, *) "dq1, dq2 are", dq1, dq2
    !endif

end subroutine
