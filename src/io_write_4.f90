subroutine io_write_4

    ! Writing out the velocities and scalars in X-space
    ! to the real*4 file

    use m_parameters
    use m_fields
    use m_work
    use x_fftw
    use m_les

    implicit none

    integer :: n_out, n, i, j, k
    real*8  :: wmag2, rkmax2

    ! every variable will undergo a mode truncation for all modes
    ! that are higher than kmax.  This will ensure that the written
    ! variables are isotropic

    rkmax2 = real(kmax, 8)**2

    ! number of variables to write out
    ! =============================================================================
    ! brooks: used to be =3 but make =6 so that pressure could be modified as well
    ! =============================================================================
    n_out = 3
    if (int_scalars) n_out = n_out + n_scalars
    if (les .and. n_les > 0) n_out = n_out + n_les

    ! putting all variables in wrk array
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx + 2
                ! magnitude of the wave number
                wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2

                if (wmag2 .gt. rkmax2) then
                    wrk(i, j, k, 1:n_out) = zip
                else
                    wrk(i, j, k, 1:n_out) = fields(i, j, k, 1:n_out)
                end if

            end do
        end do
    end do

    ! ! velocities
    call xFFT3d(-1, 1)
    fname = 'output/velocity/u.'//file_ext
    tmp4(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, 1)
    call write_tmp4

    !write(*,*) "u velocity"
    !do i=1,nx
    !  write(*,*) tmp4(i,1,1)
    !end do

    call xFFT3d(-1, 2)
    fname = 'output/velocity/v.'//file_ext
    tmp4(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, 2)
    call write_tmp4

    !write(*,*) "v velocity"
    !do i=1,nx
    !  write(*,*) tmp4(i,1,1)
    !end do

    call xFFT3d(-1, 3)
    fname = 'output/velocity/w.'//file_ext
    tmp4(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, 3)
    call write_tmp4

    !write(*,*) "w velocity"
    !do i=1,nx
    !  write(*,*) tmp4(i,1,1)
    !end do

    ! copy 4:6 from wrk so we can save the state for later
    tmp_wrk(:, :, :, 4:6) = wrk(:, :, :, 4:6)

    ! putting all variables in wrk array
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx + 2
                wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2

                if (wmag2 .gt. rkmax2) then
                    wrk(i, j, k, 4:6) = zip
                else
                    wrk(i, j, k, 4:6) = pressure_field(i, j, k, 1:3)
                end if

            end do
        end do
    end do

    ! invert all the pressures
    call xFFT3d(-1, 4)
    call xFFT3d(-1, 5)
    call xFFT3d(-1, 6)

    ! copy wrk over to pressure_field
    pressure_field(:, :, :, 1:3) = wrk(:, :, :, 4:6)

    ! set the values back to what they were before the copy
    wrk(:, :, :, 4:6) = tmp_wrk(:, :, :, 4:6)

    ! scalars
    if (int_scalars) then
        do n = 1, n_scalars
            call xFFT3d(-1, 3 + n)
            write (fname, "('output/sc',i2.2,'.',a6)") n, file_ext
            tmp4(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, 3 + n)
            call write_tmp4

        end do
    end if

    ! LES quantities
    if (les) then
        ! turbulent viscosity
        if (allocated(turb_visc)) then
            write (fname, "('output/nu_t.',a6)") file_ext
            tmp4 = turb_visc
            call write_tmp4
        end if

        if (n_les > 0) then
            do n = 1, n_les
                call xFFT3d(-1, 3 + n_scalars + n)
                write (fname, "('output/les',i1,'.',a6)") n, file_ext
                tmp4(1:nx, 1:ny, 1:nz) = wrk(1:nx, 1:ny, 1:nz, 3 + n_scalars + n)
                call write_tmp4
            end do
        end if

    end if

  !!------------------------------------------------------------
!!   getting vorticity (MGM; 07/19/2018)
!!------------------------------------------------------------
!  ! putting all variables in wrk array
!  do k = 1,nz
!     do j = 1,ny
!        do i = 1,nx+2
!           wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2
!
!           if (wmag2 .gt. rkmax2) then
!              wrk(i,j,k,1:n_out) = zip
!           else
!              wrk(i,j,k,1:n_out) = fields(i,j,k,1:n_out)
!           end if
!
!        end do
!     end do
!  end do
!
! ! Taking derivatives
!
!    call x_derivative(3,'y',6)
!    call x_derivative(3,'x',5)
!
!    call x_derivative(2,'z',4)
!    call x_derivative(2,'x',3)
!
!    call x_derivative(1,'z',2)
!    call x_derivative(1,'y',1)
!
!    wrk(:,:,:,3) = wrk(:,:,:,3) - wrk(:,:,:,1)  ! omega_3 = v_x - u_y
!    wrk(:,:,:,2) = wrk(:,:,:,2) - wrk(:,:,:,5)  ! omega_2 = u_z - w_x
!    wrk(:,:,:,1) = wrk(:,:,:,6) - wrk(:,:,:,4)  ! omega_1 = w_y - v_z
!
!    call xFFT3d(-1,1)
!    call xFFT3d(-1,2)
!    call xFFT3d(-1,3)
!
!  fname = 'output/vorticity/omgx.'//file_ext
!  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,1)
!  call write_tmp4
!
!  fname = 'output/vorticity/omgy.'//file_ext
!  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,2)
!  call write_tmp4
!
!  fname = 'output/vorticity/omgz.'//file_ext
!  tmp4(1:nx,1:ny,1:nz) = wrk(1:nx,1:ny,1:nz,3)
!  call write_tmp4

    return
end subroutine io_write_4
