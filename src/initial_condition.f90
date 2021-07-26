module m_ic
    integer :: nztotal
    real*8, allocatable, dimension(:,:,:,:) :: wrk_global
end module m_ic


subroutine write_initial_data
    use m_work
    use m_ic
    use m_fields
    use m_parameters
    use x_fftw
    implicit none

    integer :: i,j,k,v
    real*8  :: wmag2, rkmax2

    nztotal = nz * numprocs
    count = nx * ny * nz * 3

    ! ! write the summary statistics for the first time in this run
    ! ! this will include information on E(t_0), h(t_0), fe(t_0), fh(t_0)
    ! ! that we later use to calculate epsilons
    ! call write_energy(time)

    if (myid.ne.master) then
        id_to = master
        tag = myid
        call MPI_SEND( wrk(1:nx, 1:ny, 1:nz, 1:3), count, MPI_REAL8,master,tag, MPI_COMM_TASK,mpi_err )
    else
        write(*,*) "writing an initial data file"

        ! allocate a global array - since we cant run any of the other modes at the same time as a write-mode
        ! we can be sure that this array has not already been allocated
        allocate(wrk_global(nx,ny,nztotal,3))

        ! copy the master proc's wrk array in
        wrk_global(1:nx, 1:ny, 1:nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)

        do id_from=1,numprocs-1
            tag = id_from
            call MPI_RECV(wrk(1:nx, 1:ny, 1:nz, 1:3),count,MPI_REAL8,id_from,tag,MPI_COMM_TASK,mpi_status,mpi_err)
            wrk_global(1:nx, 1:ny, id_from*nz+1:(id_from+1)*nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)
        end do

        ! we now have a giant array of all the data points from each mpi process, we write them to a file now
        open(994, file="initial_condition_wrk.pkg", form='unformatted', access="stream")
        write(994) ((((wrk_global(i,j,k,v),i=1,nx),j=1,ny),k=1,nztotal),v=1,3)
        call flush(994)

    end if

end subroutine write_initial_data

subroutine load_initial_data
    use m_work
    use m_ic
    use m_parameters
    use x_fftw
    implicit none

    integer :: i,j,k,v

    nztotal = nz * numprocs
    write(*,*) " reading initial data file ", myid

    ! read in the data from the file
    ! there will be an EOF error here probably if the input field was generated for 
    ! a different N sized field
    allocate(wrk_global(nx,ny,nztotal,3))
    open(994, file="initial_condition_wrk.pkg", form='unformatted', access="stream")
    read(994) ((((wrk_global(i,j,k,v),i=1,nx),j=1,ny),k=1,nz),v=1,3)

    ! fetch the wrk data that is relevant to us
    wrk(1:nx, 1:ny, 1:nz, 1:3) = wrk_global(1:nx, 1:ny, myid*nz+1:(myid+1)*nz, 1:3)

    deallocate(wrk_global)

end subroutine load_initial_data
