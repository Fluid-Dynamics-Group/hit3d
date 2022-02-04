module m_ic
    real*8, allocatable, dimension(:, :, :, :) :: wrk_global
end module m_ic

subroutine write_initial_data
    use m_work
    use m_ic
    use m_stats
    use m_fields
    use m_parameters
    use x_fftw
    implicit none

    integer :: i, j, k, v
    real*8  :: wmag2, rkmax2

    count = nx*ny*nz*3

    if (myid .ne. master) then
        id_to = master
        tag = myid
        call MPI_SEND(fields(1:nx, 1:ny, 1:nz, 1:3), count, MPI_REAL8, master, tag, MPI_COMM_TASK, mpi_err)
    else
        write (*, *) "writing an initial data file"

        ! allocate a global array - since we cant run any of the other modes at the same time as a write-mode
        ! we can be sure that this array has not already been allocated
        allocate (wrk_global(nx, ny, nz_all, 3))

        ! copy the master proc's wrk array in
        wrk_global(1:nx, 1:ny, 1:nz, 1:3) = fields(1:nx, 1:ny, 1:nz, 1:3)

        do id_from = 1, numprocs - 1
            tag = id_from
            call MPI_RECV(fields(1:nx, 1:ny, 1:nz, 1:3), count, MPI_REAL8, id_from, tag, MPI_COMM_TASK, mpi_status, mpi_err)
            wrk_global(1:nx, 1:ny, id_from*nz + 1:(id_from + 1)*nz, 1:3) = fields(1:nx, 1:ny, 1:nz, 1:3)
        end do

        ! we now have a giant array of all the data points from each mpi process, we write them to a file now
        open (994, file="initial_condition_wrk.pkg", form='unformatted', access="stream")
        write (994) ((((wrk_global(i, j, k, v), i=1, nx), j=1, ny), k=1, nz_all), v=1, 3)
        call flush (994)
        close(994)

        ! in a previous

    end if

    ! writing espec to file
    ! we dont use any mpi here because the number of processes can change between restarts 
    ! and it is not clear how the processes should be adjusted to only write 1:kmax number 
    ! of points

    if (myid == master) then 
        open(994, file="initial_condition_espec.pkg", form="unformatted", access="stream")
        write(994) (e_spec(i), i=1,kmax)
        flush(994)
        close(994)
    end if

end subroutine write_initial_data

subroutine load_initial_velocity_data
    use m_work
    use m_ic
    use m_parameters
    use x_fftw
    implicit none

    integer :: i, j, k, v

    if (myid == master) then 
        write(*,*) "reading velocity initial condition"
    end if

    ! read in the data from the file
    ! there will be an EOF error here probably if the input field was generated for
    ! a different N sized field
    allocate (wrk_global(nx, ny, nz_all, 3))
    open (994, file="initial_condition_wrk.pkg", form='unformatted', access="stream")
    read (994) ((((wrk_global(i, j, k, v), i=1, nx), j=1, ny), k=1, nz_all), v=1, 3)
    close (994)

    ! fetch the wrk data that is relevant to us
    fields(1:nx, 1:ny, 1:nz, 1:3) = wrk_global(1:nx, 1:ny, myid*nz + 1:(myid + 1)*nz, 1:3)

    deallocate (wrk_global)

end subroutine load_initial_velocity_data

subroutine load_initial_spectral_data
    use m_ic
    use m_parameters
    use m_stats
    implicit none

    integer :: i

    if (myid == master) then 
        write(*,*) "reading spectral initial condition"
    end if

    open (994, file="initial_condition_espec.pkg", form='unformatted', access="stream")
    read (994) (e_spec(i), i=1,kmax)
    close (994)

end subroutine load_initial_spectral_data

