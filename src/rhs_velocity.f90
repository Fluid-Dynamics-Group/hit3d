module forcing_vaues
    real*8 F_1, F_2
    real*8 D_1, D_2
end module

subroutine rhs_velocity

    use m_openmpi
    use m_io
    use m_parameters
    use m_fields
    use m_work
    use x_fftw
    use m_les

    implicit none

    integer :: i, j, k, n, nx3
    real*8  :: t1(0:6), rtmp, wnum2, rnx3, epsilon1, epsilon2

    ! The IFFT of velocities has been done earlier in rhs_scalars
    ! the velocities were kept in wrk1...wrk3, intact.

!!$  ! putting the velocity field in the wrk array
!!$  wrk(:,:,:,1:3) = fields(:,:,:,1:3)
!!$  ! performing IFFT to convert them to the X-space
!!$  call xFFT3d(-1,1)
!!$  call xFFT3d(-1,2)
!!$  call xFFT3d(-1,3)

    ! right now the velocities should be in x space

    ! getting the Courant number (on the master process only)
    wrk(:, :, :, 4) = abs(wrk(:, :, :, 1)) + abs(wrk(:, :, :, 2)) + abs(wrk(:, :, :, 3))
    rtmp = maxval(wrk(1:nx, :, :, 4))
    call MPI_REDUCE(rtmp, courant, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_TASK, mpi_err)
    if (variable_dt) then
        count = 1
        call MPI_BCAST(courant, count, MPI_REAL8, 0, MPI_COMM_TASK, mpi_err)
    end if
    courant = courant*dt/dx

!-------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!  Calculating the right-hand side for the velocities
!
!  There are two options available: the standard 2/3 rule (dealias=0) and
!  combination of phase shift and truncation (dealias=1).  The latter retains
!  more modes but requires more calculations thus slowing down the simulation.
!  These are treated separately in two different "if" blocks.  This is done in
!  order not to complicate the logic.  Also this way both blocks can be
!  optimized separately.
!--------------------------------------------------------------------------------

    two_thirds_rule: if (dealias .eq. 0) then

        ! getting all 6 products of velocities
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    t1(1) = wrk(i, j, k, 1)*wrk(i, j, k, 1)
                    t1(2) = wrk(i, j, k, 1)*wrk(i, j, k, 2)
                    t1(3) = wrk(i, j, k, 1)*wrk(i, j, k, 3)
                    t1(4) = wrk(i, j, k, 2)*wrk(i, j, k, 2)
                    t1(5) = wrk(i, j, k, 2)*wrk(i, j, k, 3)
                    t1(6) = wrk(i, j, k, 3)*wrk(i, j, k, 3)

                    do n = 1, 6
                        wrk(i, j, k, n) = t1(n)
                    end do

                end do
            end do
        end do

        ! converting the products to the Fourier space
        do n = 1, 6
            call xFFT3d(1, n)
        end do

        ! Building the RHS.
        ! First, put into wrk arrays the convective terms (that will be multiplyed by "i"
        ! later) and the factor that corresponds to the diffusion

        ! Do not forget that in Fourier space the indicies are (ix, iz, iy)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 2
                    t1(1) = -(akx(i)*wrk(i, j, k, 1) + aky(k)*wrk(i, j, k, 2) + akz(j)*wrk(i, j, k, 3))
                    t1(2) = -(akx(i)*wrk(i, j, k, 2) + aky(k)*wrk(i, j, k, 4) + akz(j)*wrk(i, j, k, 5))
                    t1(3) = -(akx(i)*wrk(i, j, k, 3) + aky(k)*wrk(i, j, k, 5) + akz(j)*wrk(i, j, k, 6))

                    ! this looks like a diffusion term
                    t1(4) = -nu*(akx(i)**2 + aky(k)**2 + akz(j)**2)

                    do n = 1, 4
                        wrk(i, j, k, n) = t1(n)
                    end do

                end do
            end do
        end do

        ! call a subroutine to handle both viscous compensation and regular compensation
        ! if the viscous compensation parameter was not used
        if (PERT == 1) then
            call update_forcing_viscous_compensation(PERTamp1, PERTamp2)
        end if

        ! now take the actual fields from fields(:,:,:,:) and calculate the RHSs

        ! at this moment the contains of wrk(:,:,:,1:3) are the convective terms in the RHS
        ! which are not yet multiplied by "i"
        ! wrk(:,:,:,4) contains the Laplace operator in Fourier space.  To get the diffusion term
        ! we need to take wrk(:,:,:,4) and multiply it by the velocity

        t1(6) = real(kmax, 8)

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    ! If the dealiasing option is 2/3-rule (dealias=0) then we retain the modes
                    ! inside the cube described by $| k_i | \leq  k_{max}$, $i=1,2,3$.
                    ! The rest of the modes is purged

                    if (ialias(i, j, k) .gt. 0) then ! run for dealias=0
                        ! setting the Fourier components to zero
                        wrk(i, j, k, 1:3) = zip
                        wrk(i + 1, j, k, 1:3) = zip

                    else ! run for dealias=1

                        ! RHS for u, v and w
                        do n = 1, 3

                            ! this is what the code was like before MGM
                            ! https://github.com/Fluid-Dynamics-Group/hit3d/blob/92db6d9e864bdbf2248b025a71f315cba21d72fd/rhs_velocity.f90#L129-L135
                            if (skip_diffusion == 1) then ! skip diffusion calculation
                                ! ===========================MGM-Forcing=====================
                                ! Validation by inviscid flow
                                rtmp = -wrk(i + 1, j, k, n) + fcomp(i, j, k, n)
                                wrk(i + 1, j, k, n) = wrk(i, j, k, n) + fcomp(i + 1, j, k, n)
                                wrk(i, j, k, n) = rtmp
                            else ! dont skip the diffusion calculation
                                ! ==========================================================
                                ! BROOKS: original hit3d code with MGM forcing is here (with diffusion, i think)
                                ! ==========================================================

                                ! taking the convective term, multiply it by "i"
                                ! (see how it's done in x_fftw.f90)
                                ! and adding the diffusion term
                                rtmp = -wrk(i + 1, j, k, n) + wrk(i, j, k, 4)*fields(i, j, k, n) + fcomp(i, j, k, n)
                                wrk(i + 1, j, k, n) = wrk(i, j, k, n) + wrk(i + 1, j, k, 4)*fields(i + 1, j, k, n) + fcomp(i + 1, j, k, n)
                                wrk(i, j, k, n) = rtmp

                            end if

                        end do

                    end if

                end do
            end do
        end do

    end if two_thirds_rule

!--------------------------------------------------------------------------------
!  The second option (dealias=1).  All pairwise products of velocities are
!  dealiased using one phase shift of (dx/2,dy/2,dz/2).
!--------------------------------------------------------------------------------
    phase_shifting: if (dealias .eq. 1) then

        ! work parameters
        wrk(:, :, :, 0) = zip

        ! getting all 6 products of velocities
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    t1(1) = wrk(i, j, k, 1)*wrk(i, j, k, 1)
                    t1(2) = wrk(i, j, k, 1)*wrk(i, j, k, 2)
                    t1(3) = wrk(i, j, k, 1)*wrk(i, j, k, 3)
                    t1(4) = wrk(i, j, k, 2)*wrk(i, j, k, 2)
                    t1(5) = wrk(i, j, k, 2)*wrk(i, j, k, 3)
                    t1(6) = wrk(i, j, k, 3)*wrk(i, j, k, 3)
                    do n = 1, 6
                        wrk(i, j, k, n) = t1(n)
                    end do
                end do
            end do
        end do

        ! converting the products to the Fourier space
        do n = 1, 6
            call xFFT3d(1, n)
        end do

        ! Building the RHS.
        ! First, put into wrk arrays the convectove terms (that will be multiplyed by "i"
        ! later) and the factor that corresponds to the diffusion

        ! Do not forget that in Fourier space the indicies are (ix, iz, iy)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 2

                    t1(1) = -(akx(i)*wrk(i, j, k, 1) + aky(k)*wrk(i, j, k, 2) + akz(j)*wrk(i, j, k, 3))
                    t1(2) = -(akx(i)*wrk(i, j, k, 2) + aky(k)*wrk(i, j, k, 4) + akz(j)*wrk(i, j, k, 5))
                    t1(3) = -(akx(i)*wrk(i, j, k, 3) + aky(k)*wrk(i, j, k, 5) + akz(j)*wrk(i, j, k, 6))

                    t1(4) = -nu*(akx(i)**2 + aky(k)**2 + akz(j)**2)

                    do n = 1, 4
                        wrk(i, j, k, n) = t1(n)
                    end do
                end do
            end do
        end do

        ! now use the actual fields from fields(:,:,:,:) to calculate the RHSs

        ! at this moment the contains of wrk(:,:,:,1:3) are the convective terms in the RHS
        ! which are not yet multiplied by "i"

        ! wrk(:,:,:,4) contains the Laplace operator in Fourier space.  To get the diffusion term
        ! we need to take wrk(:,:,:,4) and multiply it by the velocity

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2

                    ! If the dealiasing option is (dealias=1) then we retain the modes
                    ! for which no more than one component of the k-vector is larger than nx/3.
                    ! The rest of the modes is purged.

                    if (ialias(i, j, k) .gt. 1) then ! run for dealias = 1
                        ! setting the Fourier components to zero
                        wrk(i, j, k, 1:3) = zip
                        wrk(i + 1, j, k, 1:3) = zip
                    else ! run for dealias = 0
                        if (skip_diffusion == 1) then ! we _ARE NOT_ doing diffusion calculations
                            ! RHS for u, v and w
                            do n = 1, 3
                                rtmp = -0.5d0*wrk(i + 1, j, k, n) !+ wrk(i  ,j,k,4) * fields(i  ,j,k,n)
                                wrk(i + 1, j, k, n) = 0.5d0*wrk(i, j, k, n) !+ wrk(i+1,j,k,4) * fields(i+1,j,k,n)
                                wrk(i, j, k, n) = rtmp
                            end do
                        else ! we _ARE_ doing diffusion calculations
                            ! RHS for u, v and w
                            do n = 1, 3
                                ! taking the HALF of the convective term, multiply it by "i"
                                ! and adding the diffusion term
                                rtmp = -0.5d0*wrk(i + 1, j, k, n) + wrk(i, j, k, 4)*fields(i, j, k, n)
                                wrk(i + 1, j, k, n) = 0.5d0*wrk(i, j, k, n) + wrk(i + 1, j, k, 4)*fields(i + 1, j, k, n)
                                wrk(i, j, k, n) = rtmp
                            end do

                        end if
                    end if

                end do
            end do
        end do

!--------------------------------------------------------------------------------
!  Second part of the phase shifting technique
!--------------------------------------------------------------------------------

        ! since wrk1...3 are taken by parts of RHS constructed earlier, we can use
        ! only wrk0 and wrk4...6.

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2

                    ! computing sines and cosines for the phase shift of dx/2,dy/2,dz/2
                    ! and putting them into wrk0
                    wrk(i, j, k, 0) = cos(half*(akx(i) + aky(k) + akz(j))*dx)
                    wrk(i + 1, j, k, 0) = sin(half*(akx(i + 1) + aky(k) + akz(j))*dx)

                    ! wrk4 will have phase-shifted u
                    wrk(i, j, k, 4) = fields(i, j, k, 1)*wrk(i, j, k, 0) - fields(i + 1, j, k, 1)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 4) = fields(i + 1, j, k, 1)*wrk(i, j, k, 0) + fields(i, j, k, 1)*wrk(i + 1, j, k, 0)

                    ! wrk5 will have phase-shifted v
                    wrk(i, j, k, 5) = fields(i, j, k, 2)*wrk(i, j, k, 0) - fields(i + 1, j, k, 2)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 5) = fields(i + 1, j, k, 2)*wrk(i, j, k, 0) + fields(i, j, k, 2)*wrk(i + 1, j, k, 0)

                end do
            end do
        end do

        ! transforming u+ and v+ into X-space
        call xFFT3d(-1, 4)
        call xFFT3d(-1, 5)

        ! now wrk4 and wrk5 contain u+ and v+

        ! getting (u+)*(u+) in real space, converting it to Fourier space,
        ! phase shifting back and adding -0.5*(the results)  to the RHS for u
        wrk(:, :, :, 6) = wrk(:, :, :, 4)**2
        call xFFT3d(1, 6)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    rtmp = wrk(i, j, k, 6)*wrk(i, j, k, 0) + wrk(i + 1, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 6) = wrk(i + 1, j, k, 6)*wrk(i, j, k, 0) - wrk(i, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i, j, k, 6) = rtmp
                end do
            end do
        end do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    wrk(i, j, k, 1) = wrk(i, j, k, 1) + 0.5d0*akx(i + 1)*wrk(i + 1, j, k, 6)
                    wrk(i + 1, j, k, 1) = wrk(i + 1, j, k, 1) - 0.5d0*akx(i)*wrk(i, j, k, 6)
                end do
            end do
        end do

        ! getting (u+)*(v+) in real space, converting it to Fourier space,
        ! phase shifting back and adding -0.5*(the results)  to the RHSs for u and v
        wrk(:, :, :, 6) = wrk(:, :, :, 4)*wrk(:, :, :, 5)
        call xFFT3d(1, 6)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    rtmp = wrk(i, j, k, 6)*wrk(i, j, k, 0) + wrk(i + 1, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 6) = wrk(i + 1, j, k, 6)*wrk(i, j, k, 0) - wrk(i, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i, j, k, 6) = rtmp
                end do
            end do
        end do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    wrk(i, j, k, 1) = wrk(i, j, k, 1) + 0.5d0*aky(k)*wrk(i + 1, j, k, 6)
                    wrk(i + 1, j, k, 1) = wrk(i + 1, j, k, 1) - 0.5d0*aky(k)*wrk(i, j, k, 6)

                    wrk(i, j, k, 2) = wrk(i, j, k, 2) + 0.5d0*akx(i + 1)*wrk(i + 1, j, k, 6)
                    wrk(i + 1, j, k, 2) = wrk(i + 1, j, k, 2) - 0.5d0*akx(i)*wrk(i, j, k, 6)
                end do
            end do
        end do

        ! getting (v+)*(v+) in real space, converting it to Fourier space,
        ! phase shifting back and adding -0.5*(the results)  to the RHS for v
        wrk(:, :, :, 6) = wrk(:, :, :, 5)**2
        call xFFT3d(1, 6)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    rtmp = wrk(i, j, k, 6)*wrk(i, j, k, 0) + wrk(i + 1, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 6) = wrk(i + 1, j, k, 6)*wrk(i, j, k, 0) - wrk(i, j, k, 6)*wrk(i + 1, j, k, 0)
                    wrk(i, j, k, 6) = rtmp
                end do
            end do
        end do
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    wrk(i, j, k, 2) = wrk(i, j, k, 2) + 0.5d0*aky(k)*wrk(i + 1, j, k, 6)
                    wrk(i + 1, j, k, 2) = wrk(i + 1, j, k, 2) - 0.5d0*aky(k)*wrk(i, j, k, 6)
                end do
            end do
        end do

        ! now get the (w+) in wrk6
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2
                    ! wrk6 will have phase-shifted w
                    wrk(i, j, k, 6) = fields(i, j, k, 3)*wrk(i, j, k, 0) - fields(i + 1, j, k, 3)*wrk(i + 1, j, k, 0)
                    wrk(i + 1, j, k, 6) = fields(i + 1, j, k, 3)*wrk(i, j, k, 0) + fields(i, j, k, 3)*wrk(i + 1, j, k, 0)
                end do
            end do
        end do
        ! transforming w+ into X-space
        call xFFT3d(-1, 6)

        ! at this point wrk4..6 contain (u+), (v+) and (w+) in real space.
        ! the combinations that we have not dealt with are: uw, vw and ww.
        ! we'll deal with all three of them at once.

        ! first get all three of these in wrk4...6 and
        wrk(:, :, :, 4) = wrk(:, :, :, 4)*wrk(:, :, :, 6)
        wrk(:, :, :, 5) = wrk(:, :, :, 5)*wrk(:, :, :, 6)
        wrk(:, :, :, 6) = wrk(:, :, :, 6)**2

        ! transform them into Fourier space
        call xFFT3d(1, 4)
        call xFFT3d(1, 5)
        call xFFT3d(1, 6)

        ! phase shift back to origianl grid and add to corresponding RHSs
        do n = 4, 6
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx + 1, 2
                        rtmp = wrk(i, j, k, n)*wrk(i, j, k, 0) + wrk(i + 1, j, k, n)*wrk(i + 1, j, k, 0)
                        wrk(i + 1, j, k, n) = wrk(i + 1, j, k, n)*wrk(i, j, k, 0) - wrk(i, j, k, n)*wrk(i + 1, j, k, 0)
                        wrk(i, j, k, n) = rtmp
                    end do
                end do
            end do
        end do

        ! adding to corresponding RHSs
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx + 1, 2

                    ! If the dealiasing option is (dealias=1) then we retain the modes
                    ! for which no more than one component of the k-vector is larger than nx/3.
                    ! The rest of the modes is purged.

                    if (ialias(i, j, k) .lt. 2) then

                        wrk(i, j, k, 1) = wrk(i, j, k, 1) + 0.5d0*akz(j)*wrk(i + 1, j, k, 4)
                        wrk(i + 1, j, k, 1) = wrk(i + 1, j, k, 1) - 0.5d0*akz(j)*wrk(i, j, k, 4)

                        wrk(i, j, k, 2) = wrk(i, j, k, 2) + 0.5d0*akz(j)*wrk(i + 1, j, k, 5)
                        wrk(i + 1, j, k, 2) = wrk(i + 1, j, k, 2) - 0.5d0*akz(j)*wrk(i, j, k, 5)

                        wrk(i, j, k, 3) = wrk(i, j, k, 3) + 0.5d0* &
                                          (akx(i + 1)*wrk(i + 1, j, k, 4) + aky(k)*wrk(i + 1, j, k, 5) + akz(j)*wrk(i + 1, j, k, 6))
                        wrk(i + 1, j, k, 3) = wrk(i + 1, j, k, 3) - 0.5d0* &
                                              (akx(i)*wrk(i, j, k, 4) + aky(k)*wrk(i, j, k, 5) + akz(j)*wrk(i, j, k, 6))

                    else
                        wrk(i:i + 1, j, k, 1) = zip
                        wrk(i:i + 1, j, k, 2) = zip
                        wrk(i:i + 1, j, k, 3) = zip
                    end if

                end do
            end do
        end do

    end if phase_shifting

    ! if performing large eddy simulations, call LES subroutine to augment
    ! the right hand side for velocioties
    les_active: if (les) then
        call les_rhs_velocity
    end if les_active

    return
end subroutine rhs_velocity

!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================

subroutine test_rhs_velocity

    use m_openmpi
    use m_io
    use m_parameters
    use m_fields
    use m_work
    use x_fftw

    implicit none

    integer :: i, j, k, n
    real*8 :: a, b, c, x, y, z

    ! defining very particular velocities so the RHS can be computed analytically

    if (task .eq. 'hydro') then

        write (out, *) 'inside.'
        call flush (out)

        a = 1.d0
        b = 1.d0
        c = 1.d0

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    x = dx*real(i - 1)
                    y = dx*real(j - 1)
                    z = dx*real(myid*nz + k - 1)

                    wrk(i, j, k, 1) = sin(a*x)
                    wrk(i, j, k, 2) = sin(b*y)
                    wrk(i, j, k, 3) = sin(c*z)
                end do
            end do
        end do

        write (out, *) 'did work'
        call flush (out)

        do n = 1, 3
            call xFFT3d(1, n)
            fields(:, :, :, n) = wrk(:, :, :, n)
        end do

        write (out, *) 'did FFTs'
        call flush (out)

        nu = 0.d0

        call rhs_velocity

        write (out, *) 'got rhs'
        call flush (out)

        do n = 1, 3
            call xFFT3d(-1, n)
        end do

        write (out, *) 'did FFTs'
        call flush (out)

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx

                    x = dx*real(i - 1)
                    y = dx*real(j - 1)
                    z = dx*real(myid*nz + k - 1)

                    ! checking u
                    wrk(i, j, k, 4) = -sin(a*x)*(two*a*cos(a*x) + b*cos(b*y) + c*cos(c*z) + nu*a**2)

                    ! checking v
                    wrk(i, j, k, 5) = -sin(b*y)*(two*b*cos(b*y) + a*cos(a*x) + c*cos(c*z) + nu*b**2)

                    ! checking w
                    wrk(i, j, k, 6) = -sin(c*z)*(two*c*cos(c*z) + b*cos(b*y) + a*cos(a*x) + nu*c**2)

                end do
            end do
        end do

!!$    do k = 1,nz
!!$      write(out,"(3e15.6)") wrk(1,1,k,3),wrk(1,1,k,5),wrk(1,1,k,4)
!!$    end do

        wrk(:, :, :, 0) = &
            abs(wrk(:, :, :, 1) - wrk(:, :, :, 4)) + &
            abs(wrk(:, :, :, 2) - wrk(:, :, :, 5)) + &
            abs(wrk(:, :, :, 3) - wrk(:, :, :, 6))

        print *, 'Maximum error is ', maxval(wrk(1:nx, :, :, 0))

        tmp4(:, :, :) = wrk(1:nx, :, :, 1) - wrk(1:nx, :, :, 4)
        fname = 'e1.arr'
        call write_tmp4

        tmp4(:, :, :) = wrk(1:nx, :, :, 2) - wrk(1:nx, :, :, 5)
        fname = 'e2.arr'
        call write_tmp4

        tmp4(:, :, :) = wrk(1:nx, :, :, 3) - wrk(1:nx, :, :, 6)
        fname = 'e3.arr'
        call write_tmp4

    end if
    return
end subroutine test_rhs_velocity

! calculates vorticity using fields (in fourier space) and sets wrk(:,:,:1-3) = omg x omgy omgz
! WARNING: This function will overwrite any values currently in wrk
subroutine calculate_vorticity()
    use m_fields
    use m_work
    use x_fftw
    use m_parameters ! dx

    implicit none

    integer :: i, j, k

    ! velocities in Fourier space
    wrk(:, :, :, 1:3) = fields(:, :, :, 1:3)
    ! Taking derivatives
    ! derivative of (1) WRT (2) -> store it in (3)
    call x_derivative(3, 'y', 6) ! dw/dy -> 6
    call x_derivative(3, 'x', 5) ! dw/dx -> 5
    call x_derivative(2, 'z', 4) ! dv/dz -> 4
    call x_derivative(2, 'x', 3) ! dv/dx -> 3
    call x_derivative(1, 'z', 2) ! du/dz -> 2
    call x_derivative(1, 'y', 1) ! du/dy -> 1

    ! getting vorticity
    wrk(:, :, :, 3) = wrk(:, :, :, 3) - wrk(:, :, :, 1)  ! omega_3 = v_x - u_y
    wrk(:, :, :, 2) = wrk(:, :, :, 2) - wrk(:, :, :, 5)  ! omega_2 = u_z - w_x
    wrk(:, :, :, 1) = wrk(:, :, :, 6) - wrk(:, :, :, 4)  ! omega_1 = w_y - v_z
end

SUBROUTINE gradient3D(m, n, o, f, hx, hy, hz, dfdx, dfdy, dfdz)
    INTEGER, INTENT(IN) :: m, n, o
    REAL*8, DIMENSION(m, n, o), INTENT(IN) :: f
    REAL*8, INTENT(IN) :: hx, hy, hz
    REAL*8, DIMENSION(m, n, o), INTENT(OUT) :: dfdx, dfdy, dfdz
    INTEGER :: i, j, k

    dfdx = 0.0; dfdy = 0.0; dfdz = 0.0; 
    !dfdx
    do k = 1, o
        do j = 1, n
            !        forward difference at start point
            dfdx(1, j, k) = (f(2, j, k) - f(1, j, k))/hx; 
            !        central difference at middle region
            dfdx(2:m - 1, j, k) = (f(3:m, j, k) - f(1:m - 2, j, k))/(2.*hx); 
            !        backward difference at end point
            dfdx(m, j, k) = (f(m, j, k) - f(m - 1, j, k))/hx; 
        end do
    end do

    !dfdy
    do k = 1, o
        do i = 1, m
            !        forward difference at start point
            dfdy(i, 1, k) = (f(i, 2, k) - f(i, 1, k))/hy; 
            !        central difference at middle region
            dfdy(i, 2:n - 1, k) = (f(i, 3:n, k) - f(i, 1:n - 2, k))/(2.*hy); 
            !        backward difference at end point
            dfdy(i, n, k) = (f(i, n, k) - f(i, n - 1, k))/hy; 
        end do
    end do

    !dfdz
    do j = 1, n
        do i = 1, m
!        forward difference at start point
            dfdz(i, j, 1) = (f(i, j, 2) - f(i, j, 1))/hz; 
!        central difference at middle region
            dfdz(i, j, 2:o - 1) = (f(i, j, 3:o) - f(i, j, 1:o - 2))/(2.*hz); 
!        backward difference at end point
            dfdz(i, j, o) = (f(i, j, o) - f(i, j, o - 1))/hz; 
        end do
    end do

END SUBROUTINE gradient3D

! the input to this subroutine has the wrk array in fourier space
! the output arrays should also be in fourier space
! wrk(:,:,:,1-3) contains the convective terms
! wrk(:,:,:,4) contains the diffusion term
subroutine update_forcing_viscous_compensation(epsilon_1, epsilon_2)
    use m_work
    use m_fields
    use x_fftw
    use forcing_vaues
    use m_parameters

    implicit none

    ! F_1, D_i and d/dt (Q_1)
    real*8 :: dQ_1, dQ_2
    real*8 :: diffusion_x, diffusion_y, diffusion_z
    real*8 :: epsilon_1, epsilon_2
    real*8 :: new_epsilon_1, new_epsilon_2
    real*8 :: f_left, f_right, f_total, diffusion
    integer :: i, j, k, n
    real*8 :: frac

    logical :: is_reading_file

    call check_read_viscous_compensation(is_reading_file, .false.)

    if (ITIME == 1) then
        write (out, *) "calculating forcing components with visc compensation routine"
    end if

    F_1 = 0.
    D_1 = 0.
    dQ_1 = 0.

    F_2 = 0.
    D_2 = 0.
    dQ_2 = 0.

    ! copy the first 3 terms (used for derivatives or something) to the temp array
    ! so that we can use the first three indicies to store the results of the diffusion term
    tmp_wrk(:, :, :, 1:3) = wrk(:, :, :, 1:3)

    ! copy the diffusion term to the 7th slot
    tmp_wrk(:, :, :, 7) = wrk(:, :, :, 4)

    ! FIRST:
    ! we evaluate the diffusion term from the RHS in fourier space
    !
    ! this code is replicated (read: Copied) from the for-loop that is directly AFTER this subroutine is called
    ! copy all the old varaibles to tmp_work
    if (viscous_compensation == 1 .or. viscous_compensation == 2) then
        call calculate_diffusion_term()
    end if

    ! after calling calculate_diffusion_term the diffusion term WRT u,v,w is now populated in wrk(:,:,:,1-3)
    ! so we can now use it to calculate whatever derivatives we may need. Swap the diffusion terms to X-space
    ! so we can use them when calculating the forcing

    do n = 1, 3
        call xFFT3d(-1, n)
    end do

    ! Lets copy these diffusion terms
    ! to tmp_wrk so that we can use all 6 values of wrk() when calculating vorticity

    tmp_wrk(:, :, :, 4:6) = wrk(:, :, :, 1:3)

    ! so now:
    ! tmp_wrk(:,:,:,1) -> some fourier derivative term we need to return to wrk() in this position
    ! tmp_wrk(:,:,:,2) -> some fourier derivative term we need to return to wrk() in this position
    ! tmp_wrk(:,:,:,3) -> some fourier derivative term we need to return to wrk() in this position
    ! tmp_wrk(:,:,:,4) -> diffusion WRT x (X-space)
    ! tmp_wrk(:,:,:,5) -> diffusion WRT y (X-space)
    ! tmp_wrk(:,:,:,6) -> diffusion WRT z (X-space)

    ! now - we dont have anything useful in the wrk array so we can use wrk to calculate the vorticity
    ! which requires 6 slots in wrk

    call calculate_vorticity()

    ! now:
    ! wrk(:,:,:,1) -> omg x (fourier-space)
    ! wrk(:,:,:,2) -> omg y (fourier-space)
    ! wrk(:,:,:,3) -> omg z (fourier-space)

    ! copy the velocities into wrk(:,:,:,4-6) so that we can shamelessly copy-paste more code from murali

    wrk(:, :, :, 4:6) = fields(:, :, :, 1:3)

    ! convert all wrk variables to X-space
    do n = 1, 6
        call xFFT3d(-1, n)
    end do

    ! now:
    ! wrk(:,:,:,1) -> omg x (x-space)
    ! wrk(:,:,:,2) -> omg y (x-space)
    ! wrk(:,:,:,3) -> omg z (x-space)
    ! wrk(:,:,:,4) -> u (x-space)
    ! wrk(:,:,:,5) -> v (x-space)
    ! wrk(:,:,:,6) -> w (x-space)

    ! ******2. Compute forcing terms & take FFTs
    ! u \cdot omg
    fcomp(:, :, :, 0) = wrk(:, :, :, 1)*wrk(:, :, :, 4) + wrk(:, :, :, 2)*wrk(:, :, :, 5) + wrk(:, :, :, 3)*wrk(:, :, :, 6)
    ! ||omg||^2
    fcomp(:, :, :, 1) = wrk(:, :, :, 1)**2 + wrk(:, :, :, 2)**2 + wrk(:, :, :, 3)**2
    ! ||u||^2
    fcomp(:, :, :, 2) = wrk(:, :, :, 4)**2 + wrk(:, :, :, 5)**2 + wrk(:, :, :, 6)**2

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                ! for each n direction
                do n = 1, 3
                    !
                    ! Calculate the forcing components for 1 and 2
                    !

                    ! u \cdot omg * omg_i - |omg|^2 * u_i
                    f_left = fcomp(i, j, k, 0)*wrk(i, j, k, n) - fcomp(i, j, k, 1)*wrk(i, j, k, 3 + n)
                    ! u \cdot omg * u_i - |u|^2 * omg_i
                    f_right = fcomp(i, j, k, 0)*wrk(i, j, k, 3 + n) - fcomp(i, j, k, 2)*wrk(i, j, k, n)

                    ! store fcomp left and right so they can be exported to a VTK file later
                    ! in a 3D flowfield

                    fcomp_individual(i, j, k, n) = f_left
                    fcomp_individual(i, j, k, n + 3) = f_right

                    ! the total forcing
                    f_total = (f_left*epsilon_1) + (f_right*epsilon_2)

                    ! integrate all the forcing components without the corresponding
                    ! epsilons into each F term
                    !
                    ! these values are outside the if statement because we want to know F_1 and F_2
                    ! for all the time steps so that we can (sometimes) output them into the energy.csv file
                    F_1 = F_1 + (f_left*wrk(i, j, k, 3 + n))
                    F_2 = F_2 + (f_right*wrk(i, j, k, n))

                    diffusion = tmp_wrk(i, j, k, n + 3)

                    ! D_1 = u \cdot d_u
                    D_1 = D_1 + wrk(i, j, k, n + 3)*diffusion

                    ! D_2 = \omega \cdot d_u
                    D_2 = D_2 + wrk(i, j, k, n)*diffusion

                    if (viscous_compensation == 1) then

                        ! d/dt Q_1 = u \cdot (d_u + f_u)
                        dQ_1 = dQ_1 + &
                               wrk(i, j, k, n + 3)*(diffusion + f_total)

                        ! d/dt Q_2 = \omega \cdot (d_u + f_u)
                        dQ_2 = dQ_2 + &
                               wrk(i, j, k, n)*(diffusion + f_total)
                    else
                        ! forcing results if no viscous compensation go in 3:5
                        fcomp(i, j, k, n + 2) = f_total
                    end if
                end do

                ! end
            end do
        end do
    end do

    ! now we have evaluated the integral so we can set the forcing components to their true value
    if (viscous_compensation == 1) then
        frac = (2*3.1415)**3/(nx*ny*nz_all)

        D_1 = D_1 * frac
        D_2 = D_2 * frac

        F_1 = F_1 * frac
        F_2 = F_2 * frac

        call add_broadcast_mpi(D_1)
        call add_broadcast_mpi(D_2)

        call add_broadcast_mpi(F_1)
        call add_broadcast_mpi(F_2)

        F_1 = -1*F_1
        F_2 = -1*F_2

        if (viscous_compensation_validation == 1) then
            dQ_1 = 0.0
            dQ_2 = 0.0
        elseif (is_reading_file) then
            ! load in the values from a file
            call read_visc_dQ_dt_values(dQ_1, dQ_2)
        end if

        if (skip_diffusion == 1) then
            D_1 = 0.0
            D_2 = 0.0
        end if

        new_epsilon_1 = epsilon_1 + ((D_1 - dQ_1)/F_1)
        new_epsilon_2 = epsilon_2 + ((D_2 - dQ_2)/F_2)

        write(*, *) "D1, DQ_1, F_1", D_1, dQ_1, F_1
        write(*, *) "ep1 inc, ep2 inc", ((D_1 - dQ_1)/F_1), ((D_2 - dQ_2)/F_2)

        ! recalculate the forcing components with the new values
        do n = 1, 3
            fcomp(:, :, :, 2 + n) = new_epsilon_1*(fcomp(:, :, :, 0)*wrk(:, :, :, n) - fcomp(:, :, :, 1)*wrk(:, :, :, 3 + n)) + &
                                    new_epsilon_2*(fcomp(:, :, :, 0)*wrk(:, :, :, 3 + n) - fcomp(:, :, :, 2)*wrk(:, :, :, n))
        end do
    end if

    ! copy either forcing result to the wrk array (we dont care about any of the data in it anymore)
    ! so that we can go back to fourier space

    wrk(:, :, :, 1:3) = fcomp(:, :, :, 3:5)

    ! Do the FFT over the 3 forcing components
    do n = 1, 3
        call xFFT3d(1, n)
    end do

    fcomp(:, :, :, 1:3) = wrk(:, :, :, 1:3)

    ! copy back whatever fourier derivatives used to be at this position
    wrk(:, :, :, 1:3) = tmp_wrk(:, :, :, 1:3)

    ! copy back the original diffusion term (fourier space) and just recalculate diffusion components as normal
    ! This is a missed optimization in the parent routine

    wrk(:, :, :, 4) = tmp_wrk(:, :, :, 7)
end

! calculate and store the result of the diffusion term
! in wrk(:,:,:,4)
! in order to call this subroutine you must have ensured that the laplace
! operation has been evaluated in fourier space in the main rhs_velocity
! initialization and that it is stored in wrk(:,:,:)
!
! This subroutine does _NOT_ take care of any forcing
subroutine calculate_diffusion_term()
    use m_work      !wrk
    use m_fields  ! fields
    use x_fftw ! ialias

    implicit none

    integer :: k, j, i, n

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx + 1, 2

                if (ialias(i, j, k) .gt. 0) then ! run for dealias=0
                    ! setting the Fourier components to zero
                    wrk(i, j, k, 1:3) = zip
                    wrk(i + 1, j, k, 1:3) = zip

                else
                    ! RHS for u, v and w
                    do n = 1, 3
                        ! for the inviscid case
                        if (skip_diffusion == 1) then
                            wrk(i + 1, j, k, n) = 0.
                            wrk(i, j, k, n) = 0.
                            ! for the viscouse case - we evaluate the diffusion term
                        else
                            !write(*,*) "calculating diffusion terms"
                            ! the reason these indicies look like this is the convective term has
                            ! a negative in front of it in the impt in the main rhs_velocity subroutine
                            ! that we never fixed.
                            !wrk(i+1,j,k,n)  =   wrk(i  ,j,k,n) + wrk(i+1,j,k,4) * fields(i+1,j,k,n)
                            !wrk(i, j, k, n) =- wrk(i+1,j,k,n)  + wrk(i  ,j,k,4) * fields(i  ,j,k,n)
                            wrk(i + 1, j, k, n) = wrk(i + 1, j, k, 4)*fields(i + 1, j, k, n)
                            wrk(i, j, k, n) = wrk(i, j, k, 4)*fields(i, j, k, n)
                        end if
                    end do

                end if

            end do
        end do
    end do

end
