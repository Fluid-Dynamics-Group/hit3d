!================================================================================
! M_IO - module that contains the I/O related subroutines.
!
! Time-stamp: <2008-11-04 15:31:31 (chumakov)>
!================================================================================

module m_io

    use m_openmpi
    implicit none

    character*6  :: file_ext
    character*80 :: fname, fname_u, fname_v, fname_w
    ! output file handle
    integer :: in = 10, out = 11

!================================================================================
contains
!================================================================================
!================================================================================
    subroutine m_io_init

        implicit none

        write (fname, "('output/d',i4.4,'.txt')") myid_world
        open (out, file=fname, position="append")
        write (out, "('-------------------------------------')")
        write (out, "('Process ',i4,' of ',i4,'(',i4.4,') is alive.')") &
            myid_world, numprocs_world, numprocs_world - 1
        write (out, "('My task is ""',a5,'"", my id is',i4)") task, myid
        call flush (out)

        return
    end subroutine m_io_init

!================================================================================
!================================================================================
    subroutine m_io_exit
        implicit none
        write (out, "('Done.')")
        close (out)
        return
    end subroutine m_io_exit

!================================================================================
!================================================================================

!================================================================================
!================================================================================
end module m_io
