module m_stats

    use m_parameters

    real*8, allocatable :: e_spec(:), e_spec1(:), moments(:, :), sc_diss(:), sc_min(:), sc_max(:)
    integer*8, allocatable :: hits(:), hits1(:)

    real*8 :: energy, eps_v, eta, etakmax, enstrophy, helicity, re_lambda, uvar, x_length
    real*8 :: lambda, re_lambda1, tau_e
    real*8 :: sctmp
    real*8 :: uf1, uf2, omgf1, omgf2

contains

!================================================================================
!  This is a small module so the arrays get allocated in all parts of the code
!================================================================================

    subroutine m_stats_init
        use m_parameters

        implicit none

        allocate (e_spec(kmax), e_spec1(kmax), moments(3 + n_scalars, 4), hits(kmax), hits1(kmax), &
                  sc_diss(n_scalars), sc_min(n_scalars), sc_max(n_scalars), stat=ierr)

        write (out, *) "Stat allocated.", ierr

        e_spec = zip
        e_spec1 = zip
        hits = 0
        hits1 = 0
        sc_diss = zip
        sc_min = zip
        sc_max = zip

        ! load e_spec and e_spec1 arrays from storage if we are loading initial conditions
        if (load_initial_condition == 0) then
            call load_initial_spectral_data()
        end if


        return
    end subroutine m_stats_init

!================================================================================
!================================================================================

    subroutine stat_main

        use m_parameters
        use m_fields
        use m_work
        use x_fftw

        implicit none

        integer :: n

        logical :: there2

        real*8 :: fac, fac2

        call stat_velocity

        if (int_scalars) then
            call stat_scalars

            ! now outputting the scalar statistics
            do n = 1, n_scalars
                if (myid .eq. 0) then

                    write (fname, "('sc',i2.2,'.gp')") n
                    inquire (file=fname, exist=there, opened=there2)
                    if (.not. there) then
                        open (105 + n, file=fname, form='formatted')
                       write (105 + n, '(A)') '# 1.itime 2.time          3.sc.diss       4. mean      5.variance    6.min     7.max'
                    end if
                    if (there .and. .not. there2) then
                        open (105 + n, file=fname, position='append')
                    end if
                    write (105 + n, "(i7,10e15.6)") itime, time, sc_diss(n), moments(3 + n, 1:2), sc_min(n), sc_max(n)
                    call flush (105 + n)
                end if
            end do

        end if

        if (task_split) then
            write (out, 9000) ITIME, TIME
            call flush (out)
        end if

        return
9000    format('ITIME=', i7, ' TIME=', e15.7, ' Stat files are written.')
    end subroutine stat_main

!================================================================================

    subroutine stat_velocity

        use m_parameters

        implicit none

        logical :: there2

        integer :: k, n

        ! getting the enstrophy
        call get_gradient_statistics

        ! getting the helicity
        ! ==========================MGM=================================
        call get_helicity

        ! getting the forcing components
        !--------------------------------------------------------------------------------
        ! ===========================MGM-Forcing=====================
        !--------------------------------------------------------------------------------
        if (PERT .eq. 1) then
            call get_forcing
        end if

        ! getting the energy spectrum e_spec to the main process
        call get_e_spec

        ! outputting the statistics into files
        if (myid .eq. 0) then

            ! getting the total energy
            energy = sum(e_spec(1:kmax))

            ! finding dissipation spectrum and total dissipation
            do k = 1, kmax
                e_spec1(k) = e_spec(k)*real(k**2, 8)*two*nu
            end do
            eps_v = sum(e_spec1(1:kmax))

            ! finding Kolmogorov scale
            eta = (nu**3/eps_v)**0.25
            etakmax = eta*real(kmax, 8)

            ! variance
            uvar = two/three*energy
            ! integral length scale
            sctmp = zip
            do k = 1, kmax
                sctmp = sctmp + e_spec(k)/real(k, 8)
            end do
            x_length = PI/two*sctmp/uvar

            ! Taylor microscale
            lambda = sqrt(15.d0*uvar*nu/eps_v)

            ! Taylor-Reynolds number
            re_lambda = uvar*sqrt(15.d0/eps_v*RE)
            re_lambda1 = sqrt(uvar)*lambda/nu

            ! Eddy turnover time
            tau_e = x_length/sqrt(uvar)

            ! outputting all necessary in the stat0 file

            inquire (file='output/stat0.gp', exist=there, opened=there2)
            if (.not. there) then
                open (69, file='output/stat0.gp', form='formatted')
                write (69, *) '# 1.itime 2.time         3.energy       4.enstrophy    5.helicity'
            end if
            if (there .and. .not. there2) then
                open (69, file='output/stat0.gp', position='append')
            end if
            write (69, "(i8,20e15.6)") itime, time, energy, enstrophy, helicity
            call flush (69)

            ! outputting all necessary in the statF file
            !--------------------------------------------------------------------------------
            ! ===========================MGM-Forcing=====================
            !--------------------------------------------------------------------------------
            if (PERT .eq. 1) then
                inquire (file='output/statF.gp', exist=there, opened=there2)
                if (.not. there) then
                    open (69, file='output/statF.gp', form='formatted')
                    write (69, *) '# 1.itime 2.time         3.uF1         4.uF2         5.omgF1         6.omgF2'
                end if
                if (there .and. .not. there2) then
                    open (69, file='output/statF.gp', position='append')
                end if
                write (69, "(i8,20e15.6)") itime, time, uf1, uf2, omgf1, omgf2
                call flush (69)
            end if

            ! outputting all this in the stat1 file

            inquire (file='output/stat1.gp', exist=there, opened=there2)
            if (.not. there) then
                open (69, file='output/stat1.gp', form='formatted')
                write (69, *) '# 1.itime 2.time         3.energy       4.diss         5.eta          6.enstrophy    7.R_lambda'
            end if
            if (there .and. .not. there2) then
                open (69, file='output/stat1.gp', position='append')
            end if
            write (69, "(i8,20e15.6)") itime, time, energy, eps_v, eta, enstrophy, re_lambda
            call flush (69)

            ! outputting all this in the stat2 file

            inquire (file='output/stat2.gp', exist=there, opened=there2)
            if (.not. there) then
                open (70, file='output/stat2.gp', form='formatted')
                write (70, '(A)') '# 1.itime  2.time         3.int LS       4. lambda      5.R_lambda1    6.tau_e        7.etakmax'
            end if
            if (there .and. .not. there2) then
                open (70, file='output/stat2.gp', position='append')
            end if
            write (70, "(i8,20e15.6)") itime, time, x_length, lambda, re_lambda1, tau_e, etakmax
            call flush (70)

            ! outputting the energy spectrum
            open (900, file='output/es.gp', position='append')
            write (900, "()")
            write (900, "()")
            write (900, "('# ITIME=',i7,' TIME=',e17.8)") ITIME, TIME
            do k = 1, kmax !min(kmax,nx/3)
                write (900, "(i4,2e15.6)") k, e_spec(k), e_spec1(k)
            end do
            close (900)

        end if
        return
    end subroutine stat_velocity

!================================================================================
!================================================================================
!================================================================================

    subroutine get_e_spec

        use m_io
        use m_fields
        use x_fftw
        implicit none

        real*8    :: sc_rad1, sc_rad2, fac, fac2
        integer :: i, j, k, n_shell

        real*8 :: energy2

        ! need this normalization factor because the FFT is unnormalized
        fac = one/real(nx*ny*nz_all)**2

        e_spec1 = zip
        e_spec = zip

        ! assembling the total energy in each shell and number of hits in each shell
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))

                    if (n_shell .gt. 0 .and. n_shell .le. kmax) then
                        fac2 = fac*(fields(i, j, k, 1)**2 + fields(i, j, k, 2)**2 + fields(i, j, k, 3)**2)
                        if (akx(i) .eq. 0.d0) fac2 = 0.5d0*fac2
                        e_spec1(n_shell) = e_spec1(n_shell) + fac2
                    end if

                end do
            end do
        end do

        ! reducing the number of hits and energy to two arrays on master node
        count = kmax
        call MPI_REDUCE(e_spec1, e_spec, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)

        return
    end subroutine get_e_spec

!================================================================================-
!================================================================================-

    subroutine get_gradient_statistics

        use m_fields
        use m_work
        use x_fftw

        implicit none
        integer :: i
        real*8  :: fac

        ! normalization  factor
        fac = one/real(nx*ny*nz_all, 8)

        do i = 1, 3
            wrk(:, :, :, i) = fields(:, :, :, i)
        end do

        ! Taking derivatives

        call x_derivative(3, 'y', 6)
        call x_derivative(3, 'x', 5)

        call x_derivative(2, 'z', 4)
        call x_derivative(2, 'x', 3)

        call x_derivative(1, 'z', 2)
        call x_derivative(1, 'y', 1)

!------------------------------------------------------------
!   getting vorticity and enstrophy
!------------------------------------------------------------
        wrk(:, :, :, 3) = wrk(:, :, :, 3) - wrk(:, :, :, 1)  ! omega_3 = v_x - u_y
        wrk(:, :, :, 2) = wrk(:, :, :, 2) - wrk(:, :, :, 5)  ! omega_2 = u_z - w_x
        wrk(:, :, :, 1) = wrk(:, :, :, 6) - wrk(:, :, :, 4)  ! omega_1 = w_y - v_z

        call xFFT3d(-1, 1)
        call xFFT3d(-1, 2)
        call xFFT3d(-1, 3)

        ! getting mean enstrophy
        wrk(:, :, :, 0) = wrk(:, :, :, 1)**2 + wrk(:, :, :, 2)**2 + wrk(:, :, :, 3)**2

        sctmp = sum(wrk(1:nx, :, :, 0))*fac
        count = 1
        call MPI_REDUCE(sctmp, enstrophy, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)

!------------------------------------------------------------
        return
    end subroutine get_gradient_statistics

    !================================================================================-
    !================================================================================-

    subroutine get_helicity
        !------------------------------------------------------------
        !   ===========================MGM===========================
        !   getting helicity
        !------------------------------------------------------------

        use m_fields
        use m_work
        use x_fftw

        implicit none
        integer :: i
        real*8  :: fac

        ! normalization  factor
        fac = one/real(nx*ny*nz_all, 8)

        do i = 1, 3
            wrk(:, :, :, i) = fields(:, :, :, i)
        end do

        ! Taking derivatives

        call x_derivative(3, 'y', 6)
        call x_derivative(3, 'x', 5)

        call x_derivative(2, 'z', 4)
        call x_derivative(2, 'x', 3)

        call x_derivative(1, 'z', 2)
        call x_derivative(1, 'y', 1)

        !------------------------------------------------------------
        !   getting vorticity and helicity
        !------------------------------------------------------------
        wrk(:, :, :, 3) = wrk(:, :, :, 3) - wrk(:, :, :, 1)  ! omega_3 = v_x - u_y
        wrk(:, :, :, 2) = wrk(:, :, :, 2) - wrk(:, :, :, 5)  ! omega_2 = u_z - w_x
        wrk(:, :, :, 1) = wrk(:, :, :, 6) - wrk(:, :, :, 4)  ! omega_1 = w_y - v_z

        call xFFT3d(-1, 1)
        call xFFT3d(-1, 2)
        call xFFT3d(-1, 3)

        do i = 4, 6
            wrk(:, :, :, i) = fields(:, :, :, i - 3)
        end do

        call xFFT3d(-1, 4)
        call xFFT3d(-1, 5)
        call xFFT3d(-1, 6)

        ! getting mean helicity
        wrk(:, :, :, 0) = wrk(:, :, :, 1)*wrk(:, :, :, 4) + wrk(:, :, :, 2)*wrk(:, :, :, 5) + wrk(:, :, :, 3)*wrk(:, :, :, 6)

        sctmp = sum(wrk(1:nx, :, :, 0))*fac
        count = 1
        call MPI_REDUCE(sctmp, helicity, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)

        !------------------------------------------------------------
        return
    end subroutine get_helicity

    !================================================================================-
    !================================================================================-

    subroutine get_forcing
        !------------------------------------------------------------
        !   ===========================MGM===========================
        !   getting forcing stats (u . f) and (omega . f)
        !------------------------------------------------------------

        use m_fields
        use m_work
        use x_fftw
        use m_parameters

        implicit none
        integer :: i
        real*8  :: fac, epsilon1, epsilon2

        ! Forcing amplitude
        epsilon1 = PERTamp1
        epsilon2 = PERTamp2

        ! normalization  factor
        fac = one/real(nx*ny*nz_all, 8)

        !do i = 1, 3
        !    wrk(:, :, :, i) = fields(:, :, :, i)
        !end do

        !! Taking derivatives

        !call x_derivative(3, 'y', 6)
        !call x_derivative(3, 'x', 5)

        !call x_derivative(2, 'z', 4)
        !call x_derivative(2, 'x', 3)

        !call x_derivative(1, 'z', 2)
        !call x_derivative(1, 'y', 1)

        !!------------------------------------------------------------
        !!   getting vorticity and helicity
        !!------------------------------------------------------------
        !wrk(:, :, :, 3) = wrk(:, :, :, 3) - wrk(:, :, :, 1)  ! omega_3 = v_x - u_y
        !wrk(:, :, :, 2) = wrk(:, :, :, 2) - wrk(:, :, :, 5)  ! omega_2 = u_z - w_x
        !wrk(:, :, :, 1) = wrk(:, :, :, 6) - wrk(:, :, :, 4)  ! omega_1 = w_y - v_z
        !! velocities
        !wrk(:, :, :, 4:6) = fields(:, :, :, 1:3)
        !! convert to X-space
        !do i = 1, 6
        !    call xFFT3d(-1, i)
        !end do

        !! ******2. Compute forcing terms & take FFTs
        !! u \cdot omg
        !fcomp(:, :, :, 0) = wrk(:, :, :, 1)*wrk(:, :, :, 4) + wrk(:, :, :, 2)*wrk(:, :, :, 5) + wrk(:, :, :, 3)*wrk(:, :, :, 6)
        !! ||omg||^2
        !fcomp(:, :, :, 1) = wrk(:, :, :, 1)**2 + wrk(:, :, :, 2)**2 + wrk(:, :, :, 3)**2
        !! ||u||^2
        !fcomp(:, :, :, 2) = wrk(:, :, :, 4)**2 + wrk(:, :, :, 5)**2 + wrk(:, :, :, 6)**2
        !! components of (u . f)
        !fcomp(:, :, :, 3) = 0.0
        !fcomp(:, :, :, 4) = 0.0
        !do i = 1, 3
        !    fcomp(:, :, :, 3) = fcomp(:, :, :, 3) + &
        !                  (wrk(:, :, :, 3 + i)*epsilon1*(fcomp(:, :, :, 0)*wrk(:, :, :, i) - fcomp(:, :, :, 1)*wrk(:, :, :, 3 + i)))
        !    !                 velocity_i * epsilon _1 * ( fcomp_0 * omega_i - fcomp_1 * velocity_i)
        !    !                 velocity_i * epsilon _1 * ( [u \cdot omg] * omega_i - [||omg||^2] * velocity_i )
        !    fcomp(:, :, :, 4) = fcomp(:, :, :, 4) + &
        !                  (wrk(:, :, :, 3 + i)*epsilon2*(fcomp(:, :, :, 0)*wrk(:, :, :, 3 + i) - fcomp(:, :, :, 2)*wrk(:, :, :, i)))
        !    !                 velocity_i * epsilon_2 * ( fcomp_0 * velocity_i - fcomp_2 * omega_i)
        !    !                 velocity_i * epsilon_2 * ( [u \cdot omg] * velocity_i - [||u||^2] * omega_i)
        !end do
        !! components of (omega . f)
        !fcomp(:, :, :, 5) = 0.0
        !fcomp(:, :, :, 6) = 0.0
        !do i = 1, 3
        !    fcomp(:, :, :, 5) = fcomp(:, :, :, 5) + &
        !                      (wrk(:, :, :, i)*epsilon1*(fcomp(:, :, :, 0)*wrk(:, :, :, i) - fcomp(:, :, :, 1)*wrk(:, :, :, 3 + i)))
        !    fcomp(:, :, :, 6) = fcomp(:, :, :, 6) + &
        !                      (wrk(:, :, :, i)*epsilon2*(fcomp(:, :, :, 0)*wrk(:, :, :, 3 + i) - fcomp(:, :, :, 2)*wrk(:, :, :, i)))
        !end do

        ! getting sum
        sctmp = sum(fcomp(1:nx, :, :, 3))*fac
        count = 1
        call MPI_REDUCE(sctmp, uf1, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
        sctmp = sum(fcomp(1:nx, :, :, 4))*fac
        count = 1
        call MPI_REDUCE(sctmp, uf2, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
        sctmp = sum(fcomp(1:nx, :, :, 5))*fac
        count = 1
        call MPI_REDUCE(sctmp, omgf1, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
        sctmp = sum(fcomp(1:nx, :, :, 6))*fac
        count = 1
        call MPI_REDUCE(sctmp, omgf2, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)

        !------------------------------------------------------------
        return
    end subroutine get_forcing

!================================================================================
!================================================================================
!================================================================================

    subroutine stat_scalars

        use m_openmpi
        use m_fields
        use m_work
        use x_fftw
        implicit none

        integer :: i, j, k, n
        real*8  :: q1, q2, fac

        if (.not. int_scalars) return

        ! getting the spectra of the scalar variances
        call get_scalar_spectra

        ! scaling factor
        fac = one/real(nx*ny*nz_all)

        ! --- Calculating moments of scalars
        do n = 1, n_scalars

            ! putting the scalar in wrk0
            wrk(:, :, :, 0) = fields(:, :, :, 3 + n)

            ! taking derivatives
            call x_derivative(0, 'x', 1)
            call x_derivative(0, 'y', 2)
            call x_derivative(0, 'z', 3)

            ! converting the derivatives to X-space
            call xFFT3d(-1, 1)
            call xFFT3d(-1, 2)
            call xFFT3d(-1, 3)

            ! getting the dissipation rate of the variance
            wrk(:, :, :, 4) = wrk(:, :, :, 1)**2 + wrk(:, :, :, 2)**2 + wrk(:, :, :, 3)**2
            q1 = two*pe(n)*sum(wrk(1:nx, :, :, 4))*fac
            count = 1
            call MPI_REDUCE(q1, q2, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
            if (myid .eq. 0) sc_diss(n) = q2

            ! converting the scalar itself to X-space
            call xFFT3d(-1, 0)

            ! First moment - mean
            q1 = sum(wrk(1:nx, :, :, 0))*fac
            count = 1
            call MPI_REDUCE(q1, q2, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
            if (myid .eq. 0) moments(3 + n, 1) = q2

            ! Second moment - variance
            q1 = sum(wrk(1:nx, :, :, 0)**2)*fac
            count = 1
            call MPI_REDUCE(q1, q2, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)
            if (myid .eq. 0) moments(3 + n, 2) = q2 - moments(3 + n, 1)**2

            ! Min and max of the scalar
            q1 = minval(wrk(1:nx, :, :, 0))
            q2 = maxval(wrk(1:nx, :, :, 0))
            count = 1
            call MPI_REDUCE(q1, sc_min(n), count, MPI_REAL8, MPI_MIN, 0, MPI_COMM_TASK, mpi_err)
            call MPI_REDUCE(q2, sc_max(n), count, MPI_REAL8, MPI_MAX, 0, MPI_COMM_TASK, mpi_err)

        end do

        return
    end subroutine stat_scalars

!================================================================================
!================================================================================

    subroutine get_scalar_spectra

        use m_io
        use m_fields
        use x_fftw
        implicit none

        real*8    :: sc_rad1, sc_rad2, fac, fac2
        integer :: i, j, k, n, n_shell

        real*8 :: energy2

        ! cycle over the scalars
        do n = 1, n_scalars

            ! need this normalization factor because the FFT is unnormalized
            fac = one/real(nx*ny*nz_all)**2

            ! using the neergy spectra arrays to keep the scalar spectra
            e_spec1 = zip
            e_spec = zip
            hits = 0
            hits1 = 0

            ! assembling the total scalar energy in each shell and number of hits in each shell
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx

                        n_shell = nint(sqrt(real(akx(i)**2 + aky(k)**2 + akz(j)**2, 4)))

                        if (n_shell .gt. 0 .and. n_shell .le. kmax) then
                            fac2 = fac*fields(i, j, k, 3 + n)**2
                            if (akx(i) .eq. 0.d0) fac2 = 0.5d0*fac2
                            e_spec1(n_shell) = e_spec1(n_shell) + fac2
                        end if

                    end do
                end do
            end do

            ! reducing the energy to two arrays on master node
            count = kmax
            call MPI_REDUCE(e_spec1, e_spec, count, MPI_REAL8, MPI_SUM, 0, MPI_COMM_TASK, mpi_err)

            ! now the master node counts the energy density in each shell
            master_node: if (myid .eq. 0) then

                ! multiplying by two because we summed only half of the scalar energy
                e_spec = 2.d0*e_spec

                ! now the master node puts the scalar energy in the file es_sc##.gp
                write (fname, "('es_',i2.2,'.gp')") n
                open (900, file=fname, position='append')
                write (900, "()")
                write (900, "()")
                write (900, "('# ITIME=',i7,' TIME=',e17.8)") ITIME, TIME
                do k = 1, kmax
                    write (900, "(i4,4e15.6)") k, e_spec(k)
                end do
                close (900)

            end if master_node

        end do

        return
    end subroutine get_scalar_spectra

!================================================================================
!================================================================================
end module m_stats
