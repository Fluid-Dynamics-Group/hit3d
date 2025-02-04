program x_code

    use m_openmpi
    use m_io
    use m_parameters
    use m_fields
    use m_work
    use x_fftw
    use m_stats
    use m_timing
    use m_force
    use m_rand_knuth
    use m_particles
    use m_filter_xfftw
    use m_les

    implicit none

    integer :: n
    integer :: itemp, jtemp, ktemp
    integer :: filenumber_temp
    real* 8 :: expectedtemp
    character :: sym
    logical :: finished_restarts
    finished_restarts = .false.

    call m_timing_init   ! Setting the time zero
    call m_openmpi_init
    call m_io_init
    ! read the input file (subroutine read_input_file)
    call m_parameters_init
    call m_les_init
    call m_fields_init
    call m_work_init
    call init_write_energy

    write(*,*) "finished calling first init functions"

    ! allocating and initializing FFTW arrays
    call x_fftw_allocate(1)
    call x_fftw_init

    call m_stats_init
    call m_force_init

    ! allocating and initializing particles
    if (task .eq. 'parts') then
        call particles_init
    end if

    if (NRESCALE == 0) finished_restarts = .true.

    write (out, *) "IN THE PROGRAM."
    call flush (out)

    ! initializing the random number generator
    ! call rand_knuth_init

    ! getting the wallclock runlimit for the job
    call get_job_runlimit

!-----------------------------------------------------------------------
!     Starting from the beginning or from the saved flowfield
!-----------------------------------------------------------------------
    if (load_initial_condition == 2 .or. load_initial_condition == 1) then
        call begin_new
    elseif (load_initial_condition == 0) then
        call begin_new
        call load_initial_spectral_data()
    else
        write (*, *) "invalid value for load_initial_condition"
        stop 1
    end if

    ! Initializing the LES stuff
    if (les) call m_les_begin

    ! checking divergence
    if (task .eq. 'hydro') call divergence

    ! indicators whether to use first-order in time
    ! for velocities and scalars
    fov = .true.
    fos = .true.

    ! need to dealias the fields at the beginning
    if (task .eq. 'hydro') call dealias_all


!********************************************************************************
!  call benchmark
!********************************************************************************

!================================================================================
!  MAIN CYCLE
!================================================================================

    do 100 ITIME = ITMIN + 1, ITMAX

        ! output some runtime statistics on how far into the simulation we are
        if (mod(ITIME, 1000) == 0) then
            if (myid == master) then
                write(*,*) 100* ITIME / ITMAX, " percent done"
            endif
        endif

        ! getting the file extension for current iteration
        call get_file_ext
!--------------------------------------------------------------------------------
!  now performing the core of the cycle.
!  This is done with "if" rather than "select case" because if we're not
!  splitting tasks then we want everything to be done consequently by the
!  same set of processors.
!
!  All the syncronization calls (fields_to_parts, fields_to_stats) will be
!  called only if (task_split).
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!                             HYDRO PART
!   note that even in the case where there is no task splitting,
!   'hydro' part is still there.  all processors will have task = 'hydro'
!--------------------------------------------------------------------------------
        hydro: if (task .eq. 'hydro') then

            ! ------------------------------------------------------------
            ! taking care of rescaling when running decaying turbulence
            ! if the time just was divisible by TRESCALE
            ! ------------------------------------------------------------
            if (flow_type .eq. 0 .and. floor((time - dt)/TRESCALE) .lt. floor(time/TRESCALE)) then
                ! ...and if we haven't rescaled NRESCALE times
                if (floor(time/TRESCALE) .le. NRESCALE .and. itime .ne. 1) then
                    write (out, *) "MAIN: Rescaling velocities"
                    call flush (out)
                    call velocity_rescale
                    ! after rescaling, the time-stepping needs to be first order
                    fov = .true.; fos = .true.
                    if (.not. task_split .and. mod(itime, iprint1) .eq. 0) call stat_main

                    if (floor(time/TRESCALE) == NRESCALE) then
                        finished_restarts = .true.
                        write (*, *) "finishing restarts at time step = ", itime
                    end if
                end if
            end if

            ! RHS for scalars
            call rhs_scalars

            ! now the velocities in x-space are contained in wrk1...3
            ! if we are moving particles, then we want to send the velocity field
            ! to the "parts" part of the code
            if (task_split) call fields_to_parts

            ! advance scalars - either Euler or Adams-Bashforth
            if (int_scalars .or. n_les > 0) then
                call flush (out)
                n = 3 + n_scalars + n_les

                if (fos) then
                    rhs_old(:, :, :, 4:n) = wrk(:, :, :, 4:n)
                    fields(:, :, :, 4:n) = fields(:, :, :, 4:n) + dt*rhs_old(:, :, :, 4:n)
                    fos = .false.
                else
                    fields(:, :, :, 4:n) = fields(:, :, :, 4:n) + &
                                           dt*(1.5d0*wrk(:, :, :, 4:n) - 0.5d0*rhs_old(:, :, :, 4:n))
                    rhs_old(:, :, :, 4:n) = wrk(:, :, :, 4:n)
                end if

            end if

            ! RHS for velocities
            call rhs_velocity

            ! adding forcing, if computing forced flow
            if (flow_type .eq. 1) call force_velocity

            ! advance velocity - either Euler or Adams-Bashforth
            if (fov) then
                ! copy to this new variable to work with later
                rhs_saved(:, :, :, :) = wrk(:, :, :, 1:3)

                rhs_old(:, :, :, 1:3) = wrk(:, :, :, 1:3)
                fields(:, :, :, 1:3) = fields(:, :, :, 1:3) + dt*rhs_old(:, :, :, 1:3)
                fov = .false.
            else
                ! dot this varaible with u and omega to get helicity
                rhs_saved(:, :, :, :) = 1.5d0*wrk(:, :, :, 1:3) - 0.5d0*rhs_old(:, :, :, 1:3)

                fields(:, :, :, 1:3) = fields(:, :, :, 1:3) + &
                                       dt*(rhs_saved)
                rhs_old(:, :, :, 1:3) = wrk(:, :, :, 1:3)
            end if

            ! solve for pressure and update velocities so they are incompressible
            ! we solve for the pressure after finding the pressure
            ! fractionary step method
            call pressure

            ! advance the time
            TIME = TIME + DT

            ! write the restart file if it's the time
            ! Brooks - dont use his restart files
            !if (mod(itime, IPRINT2) .eq. 0) call restart_write_parallel

            ! change the timestep in case we're running with variable timestep
            if (variable_dt) call my_dt

            ! CPU usage statistics
            if (mod(itime, iprint1) .eq. 0) then
                call m_timing_check
                if (mod(itime, iwrite4) .eq. 0) then
                    sym = "*"
                else
                    sym = " "
                end if
                write (out, 9000) itime, time, dt, courant, cpu_hrs, cpu_min, cpu_sec, &
                    sym, les_model_name
                call flush (out)
            end if

            if (mod(itime, iprint1) .eq. 0 .or. mod(itime, iwrite4) .eq. 0) then

                ! send the velocities to the "stats" part of the code for statistics
                if (task_split) call fields_to_stats
                ! checking if we need to stop the calculations due to simulation time
                ! Brooks - removed a stop to simulation time
                ! if (TIME .gt. TMAX) call my_exit(1)

                ! checking if we need to start advancing scalars
                if (n_scalars .gt. 0 .and. .not. int_scalars .and. time .gt. TSCALAR) then
                    int_scalars = .true.
                    call init_scalars
                    write (out, *) "Starting to move the scalars."
                    call flush (out)
                end if

            end if
        end if hydro

        ! if we are at the specified timestep and our job is to write a restart file ...
        if (ITIME >= ITMAX .and. load_initial_condition == 1) then
            ! write a full flowfield file
            call write_velocity_field(int(itime))
            ! write the energy info at the final timestep
            call write_energy(time)
            ! create a file to restart from
            call write_initial_data()
            call my_exit(0)
            call m_openmpi_exit
            stop
        ! write data at the first iteration of any IWRITE4 (config file) point after
        elseif (mod(itime, iwrite4) == 0 .or. itime == ITMIN+1) then
            ! write energy and helicity to a csv
            call write_energy(time)
            call write_slice(int(itime))
            call write_scalars(int(itime))


            ! we only write vtk files once every *8 time steps because the post processessing
            ! is very slow
            ! also only write them after the restarts are done
            if ( &
                ( &
                    finished_restarts .AND. mod(itime, iwrite4*4) .eq. 0) &
                    .or. &
                    itime == ITMIN + 1 &
                ) then
                write(out, *) "writing flowfield at itime ", itime
                call write_velocity_field(int(itime))
            end if

            if (export_divergence == 1) then
                call write_derivatives(int(itime));
            end if
        end if

!--------------------------------------------------------------------------------
!                             STATISTICS PART
!--------------------------------------------------------------------------------
        stats: if (task .eq. 'stats' .or. .not. task_split) then

            if (mod(itime, iprint1) .eq. 0 .or. mod(itime, iwrite4) .eq. 0) then
                ! if this is a separate set of processors, then...
                stats_task_split: if (task_split) then
                    ! checking if we need to stop the calculations due to simulation time
                    ! Brooks - remove any stops due to simulation time
                    ! if (TIME .gt. TMAX) call my_exit(1)
                end if stats_task_split

                ! these are executed regardless of the processor configuration
                if (task_split) call fields_to_stats
                if (mod(itime, iprint1) .eq. 0) call stat_main
                ! we are at a write timestep and we are also not writing initial condition data to a file
                if (mod(itime, iwrite4) .eq. 0 .and. load_initial_condition /= 1) then
                    !call io_write_4
                end if

            end if
        end if stats

!--------------------------------------------------------------------------------
!                             PARTICLE PARTS
!  NOTE: This is not enabled to work when not task_split.
!  Need to return to it later.
!  Currently the particles can be calculated only if we split the tasks due to
!  requirements on the wrk array sizes in the particle interpolation routines.
!--------------------------------------------------------------------------------
        particles: if (task .eq. 'parts') then

            call fields_to_parts

            if (int_particles) then
                call particles_move
                if (mod(itime, iwrite4) .eq. 0) call particles_restart_write_binary
            end if

            ! remove any exit due to time overrun
            !if (mod(itime, iprint1) .eq. 0 .or. mod(itime, iwrite4) .eq. 0) then
            !    if (TIME .gt. TMAX) call my_exit(1)
            !end if
        end if particles
!!$!--------------------------------------------------------------------------------
!!$!                             OTHER PARTS
!!$!--------------------------------------------------------------------------------
!!$
!!$
!!$     write(out,*) "skipping the time step",ITIME
!!$     call flush(out)
!!$
!!$
!--------------------------------------------------------------------------------
!                             COMMON PARTS
!--------------------------------------------------------------------------------

        ! every 10 iterations checking
        ! 1) for the run time: are we getting close to the job_runlimit?
        ! 2) for the user termination: is there a file "stop" in directory?
        if (mod(ITIME, 10) .eq. 0) then

            ! synchronize all processors, hard
!!$     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)

            if (myid_world .eq. 0) call m_timing_check
            count = 1
            call MPI_BCAST(cpu_min_total, count, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpi_err)

            ! allowing 5 extra minutes for writing restart file
            ! note that for large-scale calculations (e.g. 1024^3)
            ! the restart writing time can be long (up to 20 minutes or so).
            ! this should be taken care of in the job submission script
            ! via file job_parameters.txt
            
            ! Brooks - no early exit to the simulation
            ! if (cpu_min_total + 5 .gt. job_runlimit) call my_exit(2)

            ! user termination.  If the file "stop" is in the directory, stop
            inquire (file='stop', exist=there)
            if (there) call my_exit(3)
        end if

100     continue
!================================================================================

!--------------------------------------------------------------------------------
!  In a case when we've gone to ITMAX, write the restart file
!--------------------------------------------------------------------------------

        ITIME = ITIME - 1
        if (task .eq. 'hydro') call restart_write_parallel
        call my_exit(0)
        call m_openmpi_exit

        stop
9000    format('ITIME=', i6, 3x, 'TIME=', f8.4, 4x, 'DT=', f8.5, 3x, 'Courn= ', f6.4, &
               2x, 'CPU:(', i4.4, ':', i2.2, ':', i2.2, ')', x, a1, x, a3)

        end program x_code

!=============================================================================
        subroutine benchmark

            use m_openmpi
            use m_io
            use m_parameters
            use m_fields
            use m_work
            use x_fftw
            use m_stats
            use m_timing
            use m_force
            use m_rand_knuth
            use m_particles
            use m_filter_xfftw
            use m_les

            implicit none

            integer :: n, nmax

            call m_timing_init
            benchmarking = .true.
            bm = 0

            wrk(:, :, :, 1) = fields(:, :, :, 1)

            do n = 1, nmax
                call xfft3d(1, 1)
                call xfft3d(-1, 1)
            end do

            if (myid .eq. 0) then
                write (out, *) "BENF: statistics on forward transform"
                write (out, *) "BENF: R2C: ", bm(1)/nmax
                write (out, *) "BENF: T13: ", bm(2)/nmax
                write (out, *) "BENF: C2C: ", bm(3)/nmax
                write (out, *) "BENF: ====="
                write (out, *) "BENF: TOT: ", bm(11)/nmax
                write (out, *) "BENB: statistics on backward transform"
                write (out, *) "BENB: C2C: ", bm(4)/nmax
                write (out, *) "BENB: T13: ", bm(5)/nmax
                write (out, *) "BENB: C2R: ", bm(6)/nmax
                write (out, *) "BENB: ====="
                write (out, *) "BENB: TOT: ", bm(12)/nmax
            end if

            close (out)

            stop
        end subroutine benchmark
