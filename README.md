# Possible UB in HIT3d

`akx` `aky` and `akz` are defined in `x_fftw.f90`:

```fortran
    do ix = 1, nx + 1, 2
       akx(ix) = real((ix - 1) / 2, 8)		! <--- here
       akx(ix + 1) = akx(ix)
       coskx2(ix) = dcos(half * akx(ix))
       sinkx2(ix) = dsin(half * akx(ix))
       coskx2(ix + 1) = coskx2(ix)
       sinkx2(ix + 1) = sinkx2(ix)
       rezkax(ix) = 0
       if (dabs(akx(ix)) > (real(nz_all, 8)) / 3.0D0) rezkax(ix) = 1
    end do
    ! in Fourier space ky-axis is distributed among the processors
    do iy = 1, nz
       aky(iy) = real(myid * nz + iy - 1, 8)	! <--- here
       if (aky(iy) > (0.5D0 * real(ny, 8))) aky(iy) = aky(iy) - real(ny, 8)
       cosky2(iy) = dcos(half * aky(iy))
       sinky2(iy) = dsin(half * aky(iy))
       rezkay(iy) = 0
       if (dabs(aky(iy)) > (real(ny, 8)) / 3.0D0) rezkay(iy) = 1
    end do
    ! in Fourier space the z wavenumbers are aligned along the second index
    do iz = 1, ny
       akz(iz) = real(iz - 1, 8)		! <--- here
       if (akz(iz) > (0.5D0 * real(nz_all, 8))) akz(iz) = akz(iz) - real(nz_all, 8)
       coskz2(iz) = dcos(half * akz(iz))
       sinkz2(iz) = dsin(half * akz(iz))
       rezkaz(iz) = 0
       if (dabs(akz(iz)) > (real(nz_all, 8)) / 3.0D0) rezkaz(iz) = 1
    end do
```

If you notice, if the mpi rank `myid` is zero (for instance, on the master process)
all of `akx(1)`, `aky(1)` and `akz(1)` will be zero. This routine is called in `main.f90`:


```
  implicit none

  integer :: n
  character :: sym

  call m_timing_init   ! Setting the time zero
  call m_openmpi_init
  call m_io_init
  call m_parameters_init
  call m_les_init
  call m_fields_init
  call m_work_init


  ! allocating and initializing FFTW arrays
  call x_fftw_allocate(1)
  call x_fftw_init				! < --------  called here

  call m_stats_init
  call m_force_init

  ! allocating and initializing particles
  if (task.eq.'parts') then
     call particles_init
  end if

  write(out,*) "IN THE PROGRAM."
  call flush(out)

  ! initializing the random number generator
  ! call rand_knuth_init

  ! getting the wallclock runlimit for the job
  call get_job_runlimit


!-----------------------------------------------------------------------
!     Starting from the beginning or from the saved flowfield
!-----------------------------------------------------------------------
  if(ITMIN.eq.0) then
     call begin_new				! <---- where we are going next
  else
     call begin_restart
  endif
```

Then, in `begin_new`:

```
subroutine begin_new

  use m_openmpi
  use m_parameters
  use m_fields
  use m_stats
  implicit none


  ! defining time
  TIME = zip

  ! deciding if we advance scalars or not
  if (TSCALAR.le.zip .and. n_scalars.gt.0) int_scalars = .true.

  ! defining the iteration number
  ITIME = 0
  file_ext = '000000'



  if (task.eq.'hydro') then
    call init_velocity						! <-----  we use "hydro", so this routine gets called
    if (n_scalars.gt.0) call init_scalars
    call io_write_4
  end if

  if (task_split) then
     call fields_to_stats
     if (task.eq.'stats') call stat_main
  end if

  return
end subroutine begin_new
```

The problem with `ak_(1)` arises in `init_velocity.f90` deep inside the `init_velocity` subroutine.
At line 222:

```
    write(*,*) "akx aky akz", akx(1), aky(1), akz(1)
    write(*,*) "nshell before indexing is", nint(sqrt(real(akx(1)**2 + aky(1)**2 + akz(1)**2, 4)))
    write(*,*) "espec at nshell", e_spec(nint(sqrt(real(akx(1)**2 + aky(1)**2 + akz(1)**2, 4))))
  ! now go over all Fourier shells and multiply the velocities in a shell by
  ! the sqrt of ratio of the resired to the current spectrum
  do k = 1,nz
     do j = 1,ny
        do i = 1,nx+2

           n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))
           if (n_shell .gt. 0 .and. n_shell .le. kmax .and. e_spec(n_shell) .gt. zip) then
              fields(i,j,k,1:3) = fields(i,j,k,1:3) * sqrt(e_spec1(n_shell)/e_spec(n_shell))
           else
              fields(i,j,k,1:3) = zip
           end if

        end do
     end do
  end do
```

You will notice that `n_shell` is defined based on the magnitude of `ak_( )`. At `i,j,k = 1`, so since
`akx(1) = 0.`, `aky(1) = 0.`, `akz(1) = 0.`, the magnitude on the master node is `n_shell = 0`. 
On the next line, we see that `n_shell` is used to index `e_spec`, and `e_spec` is allocated
for indicied `e_spec(1:kmax)`.

Running the code with compiler optimizations removes this issue (read: undefined behavior):

```
FCFLAGS = -O4 -c -fallow-argument-mismatch
LDFLAGS = -O4 -L/usr/local/lib -lfftw3
```

Running: 

```
# from src/
make clean && make && mpirun -np 1 ./hit3d.x "input_file.in" "nosplit"
```

The output from the `write` statements is:

```
 akx aky akz   0.0000000000000000        0.0000000000000000        0.0000000000000000
 nshell before indexing is           0
 espec at nshell   8.6955553668059392E-322
```

Then, if you run the code without compiler optimizations:

```
FCFLAGS = -O -c -fallow-argument-mismatch -ffree-line-length-none -ffixed-line-length-none -Wall -fcheck=all -g -fbacktrace 
LDFLAGS = -O -L/usr/local/lib -lfftw3
```

```
# from src/
make clean && make && mpirun -np 1 ./hit3d.x "input_file.in" "nosplit"
```

```
akx aky akz   0.0000000000000000        0.0000000000000000        0.0000000000000000
 nshell before indexing is           0
At line 224 of file init_velocity.f90
Fortran runtime error: Index '0' of dimension 1 of array 'e_spec' below lower bound of 1

Error termination. Backtrace:
#0  0x55d131750cbc in init_velocity_
        at /home/brooks/github/hit3d/src/init_velocity.f90:224
#1  0x55d13174c027 in begin_new_
        at /home/brooks/github/hit3d/src/begin_new.f90:23
#2  0x55d1317486a3 in x_code
        at /home/brooks/github/hit3d/src/main.f90:57
#3  0x55d13174b6a3 in main
        at /home/brooks/github/hit3d/src/main.f90:3
--------------------------------------------------------------------------
Primary job  terminated normally, but 1 process returned
a non-zero exit code. Per user-direction, the job has been aborted.
--------------------------------------------------------------------------
--------------------------------------------------------------------------
mpirun detected that one or more processes exited with non-zero status, thus causing
the job to be terminated. The first process to do so was:

  Process name: [[42808,1],0]
  Exit code:    2
--------------------------------------------------------------------------
```
