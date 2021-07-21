module m_data
    integer :: nz_global
    real*8, dimension(:,:,:), allocatable :: dUdX, dUdY, dUdZ, dVdX, dVdY, dVdZ, dWdX, dWdY, dWdZ
    real*8, dimension(:,:,:), allocatable :: dPxdX, dPxdY, dPxdZ, dPydX, dPydY, dPydZ, dPzdX, dPzdY, dPzdZ
    real*8, dimension(:,:,:), allocatable :: OmgX, OmgY, OmgZ
    ! second order derivatives for third term in dE/dt
    real*8, dimension(:,:,:), allocatable :: dUdX2, dUdY2, dUdZ2, dVdX2, dVdY2, dVdZ2, dWdX2, dWdY2, dWdZ2
    ! curl ( RHS) terms 
    real*8, dimension(:,:,:), allocatable :: dRHS1dX, dRHS1dY, dRHS1dZ, dRHS2dX, dRHS2dY, dRHS2dZ, dRHS3dX, dRHS3dY, dRHS3dZ

    ! an array of ALL the data globally [i,j,k, (u/v/w)]
    real*8, dimension(:,:,:,:), allocatable :: velocity_all, pressure_all

end module m_data

subroutine write_velocity_field(current_timestep)
    use m_work ! wrk()
    use m_parameters ! nx, ny, nz

    implicit none

    integer :: filenumber, meta_file, i, j, k, current_timestep
    !real*8, allocatable :: fields(:,:,:,:)
    real*8 :: u, v,w

    character(len=60) :: filename, sizefile

    filenumber = 619

    write(filename, "('output/velocity_field/', i2.2, '_', i5.5, '.csv')") myid_world, current_timestep
    open(filenumber, file=filename, status="new")

    write(filenumber, "(A5)") "u,v,w"

    do k = 1,nz
      do j = 1,ny
          do i = 1,nx
              u = wrk(i,j,k,1)
              v = wrk(i,j,k,2)
              w = wrk(i,j,k,3)
              write(filenumber, "(E16.10, A1, E16.10, A1, E16.10)") u, ",", v, ",", w
          end do
      end do
    end do

    flush(filenumber)

    meta_file = 620

    !write(sizefile, "('output/size_', i2.2, '.csv')") myid_world

    open(meta_file, file="output/size.csv")

    write(meta_file, "(A8)") "nx,ny,nz"
    write(meta_file, "(I5.5, A1, I5.5, A1, I5.5)") nx, ",", ny, ",", nz

end subroutine write_velocity_field

! create the file name and open the CSV files that we are using to write E(t) and h(t)
subroutine init_write_energy
    use m_parameters
    use m_data
    implicit none

    integer:: filenumber
    character(len=40) filename

    filenumber = 621
    call create_energy_filename(filename)

    open(filenumber, file=filename)
    write(filenumber, "('energy,helicity,advection_u,pressure_u,diffusion_u,divu,current_time,&
        solver_energy,solver_helicity,enstrophy,advection_omega,pressure_omega,diffusion_omega,alt_dh_dt')")

    filenumber = 622
    open(filenumber, file="output/velocity.csv")
    write(filenumber, "('u,v,w')")

    nz_global = numprocs * nz

    allocate(dUdX(nx,ny,nz_global));allocate(dUdY(nx,ny,nz_global));allocate(dUdZ(nx,ny,nz_global));
    allocate(dVdX(nx,ny,nz_global));allocate(dVdY(nx,ny,nz_global));allocate(dVdZ(nx,ny,nz_global));
    allocate(dWdX(nx,ny,nz_global));allocate(dWdY(nx,ny,nz_global));allocate(dWdZ(nx,ny,nz_global));
    allocate(OmgX(nx,ny,nz_global));allocate(OmgY(nx,ny,nz_global));allocate(OmgZ(nx,ny,nz_global));

    ! second order vars for Del^2(u)
    allocate(dUdX2(nx,ny,nz_global));allocate(dUdY2(nx,ny,nz_global));allocate(dUdZ2(nx,ny,nz_global));
    allocate(dVdX2(nx,ny,nz_global));allocate(dVdY2(nx,ny,nz_global));allocate(dVdZ2(nx,ny,nz_global));
    allocate(dWdX2(nx,ny,nz_global));allocate(dWdY2(nx,ny,nz_global));allocate(dWdZ2(nx,ny,nz_global));

    ! pressure terms
    allocate(dPxdX(nx,ny,nz_global));allocate(dPxdY(nx,ny,nz_global));allocate(dPxdZ(nx,ny,nz_global));
    allocate(dPydX(nx,ny,nz_global));allocate(dPydY(nx,ny,nz_global));allocate(dPydZ(nx,ny,nz_global));
    allocate(dPzdX(nx,ny,nz_global));allocate(dPzdY(nx,ny,nz_global));allocate(dPzdZ(nx,ny,nz_global));
    
    ! RHS curl terms
    allocate(dRHS1dX(nx,ny,nz_global));allocate(dRHS1dY(nx,ny,nz_global));allocate(dRHS1dZ(nx,ny,nz_global));
    allocate(dRHS2dX(nx,ny,nz_global));allocate(dRHS2dY(nx,ny,nz_global));allocate(dRHS2dZ(nx,ny,nz_global));
    allocate(dRHS3dX(nx,ny,nz_global));allocate(dRHS3dY(nx,ny,nz_global));allocate(dRHS3dZ(nx,ny,nz_global));

    ! global pressure / velocity arrays
    allocate(velocity_all(nx,ny,nz_global,3))
    allocate(pressure_all(nx,ny,nz_global,3))

end subroutine init_write_energy

! initalize the name of the file that we are writing to with E(t) and h(t) values
subroutine create_energy_filename(filename)
    use m_parameters
    character(len=40) :: filename
    write(filename, "('output/energy/energy_', i2.2, '_.csv')") myid_world
end subroutine create_energy_filename

! calculate and write the energy and helicity values for the current time step
! you must ensure that you call init_write_energy to open the correct files
! before calling this subroutine
subroutine write_energy(current_time)
    use m_work ! wrk
    use m_parameters
    implicit none

    integer::filenumber
    character(len=40) :: filename
    real*8 :: energy, helicity, advection_u, pressure_u, diffusion_u, divu, current_time, solver_energy, solver_helicity,enstrophy
    real*8 :: advection_omega, pressure_omega, diffusion_omega, alt_dh_dt

    call collect_mpi_arrays()

    if myid == 0 then
        ! initialize the name of the csv that this mpi process will write to
        call create_energy_filename(filename)

        ! ifft the RHS variables so we are in x-space
        call ifft_rhs

        filenumber = 621

        ! calculate the variables
        call calculate_energy(energy, solver_energy)
        call calculate_helicity(helicity,solver_helicity, enstrophy)
        ! all the terms in dE/dt
        call advection_dot(advection_u, advection_omega)
        call pressure_dot(pressure_u, pressure_omega)
        call diffusion_dot(diffusion_u, diffusion_omega)

        call calculate_divu(divu)

        call alternate_dh_dt_calculation(alt_dh_dt)

        ! The file will have already been opened by the init process
        ! now we write the calculated data to the file
        open(filenumber)
        write(filenumber, "(E16.10, ',', E16.10, ',', E16.10, ',', E16.10, &
            ',',E16.10, ',',E16.10, ',', E16.10, ',', E16.10, ',', E16.10, &
            ',',E16.10, ',',E16.10, ',', E16.10, ',', E16.10, ',', E16.10)") &
            energy, helicity, advection_u, pressure_u, diffusion_u, divu, current_time, solver_energy, solver_helicity, enstrophy,&
            advection_omega, pressure_omega, diffusion_omega,alt_dh_dt

        flush(filenumber)

        filenumber = 622
        open(filenumber)
        write(filenumber, "(E16.10, ',', E16.10, ',', E16.10)") wrk(32,32,32,1), wrk(32,32,32,2), wrk(32,32,32,3)
        flush(filenumber)
    end if
end subroutine

subroutine ifft_rhs
    use m_work !tmp_wrk + wrk + rhs_saved
    use x_fftw ! fft stuff
    use m_parameters ! kmax

    ! save the current wrk array
    tmp_wrk(:,:,:,1:3) = wrk(:,:,:,1:3)

    ! copy the RHS variables into wrk
    wrk(:,:,:,1:3) = rhs_saved(:,:,:,1:3)

    ! do the magic truncation stuff here
    rkmax2 = real(kmax,8)**2

    ! number of variables to write out
    ! =============================================================================
    ! brooks: used to be =3 but make =6 so that pressure could be modified as well
    ! =============================================================================
    ! putting all variables in wrk array
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+2
             ! magnitude of the wave number
             wmag2 = akx(i)**2 + aky(k)**2 + akz(j)**2

             if (wmag2 .gt. rkmax2) then
                wrk(i,j,k,1:3) = zip
             else
                wrk(i,j,k,1:3) = wrk(i,j,k,1:3)
             end if

          end do
       end do
    end do

    ! do the iffts using the wrk array
    call xFFT3d(-1,1)
    call xFFT3d(-1,2)
    call xFFT3d(-1,3)

    ! copy the data out of wrk and into the rhs
    rhs_saved(:,:,:,1:3) = wrk(:,:,:,1:3)

    ! copy the original data back into wrk
    wrk(:,:,:,1:3) = tmp_wrk(:,:,:,1:3)

end subroutine ifft_rhs

subroutine calculate_energy(energy, solver_energy)
    use m_work !wrk
    use m_parameters ! dx, dy, dz

    implicit none
    integer:: i, j, k
    real*8 :: energy, solver_energy, u, v, w

    energy = 0
    solver_energy = 0

    !write(*,*) "energy"

    do i =1,nx
        do j=1,ny
            do k=1,nz
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                energy = energy + u**2 + v**2 + w**2
                solver_energy = solver_energy + (rhs_saved(i,j,k,1) * u + rhs_saved(i,j,k,2) * v +rhs_saved(i,j,k,3) * w)
            end do
        end do
    end do

    energy = energy * dx * dy * dz * 0.5
    solver_energy = solver_energy * dx * dy * dz

end subroutine calculate_energy

subroutine calculate_helicity(helicity, solver_helicity, enstrophy)
    use m_work !wrk
    use m_data !allocatable arrays for vorticity calculation
    use m_parameters ! dx, dy, dz

    implicit none
    real*8 :: helicity, solver_helicity, tmp, enstrophy
    real*8 :: u,v,w
    real*8 :: omg_x, omg_y, omg_z
    integer :: i, j, k

    helicity = 0.
    solver_helicity = 0.
    enstrophy = 0.

    ! write(*,*) "helicity"

    ! start by calculating the vorticity,
    ! store the results in OmgX, OmgY, OmgZ
    call calculate_vorticity

    write(*,*) "got back from calculate vorticity, printing all values"
    do k =1,nz
        do j=1,ny
            do i=1,nx
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                ! if (u/=u) write(*,*) "u nan"
                ! if (v/=v) write(*,*) "v nan"
                ! if (w/=w) write(*,*) "w nan"

                omg_x = OmgX(i,j,k)
                omg_y = OmgY(i,j,k)
                omg_z = OmgZ(i,j,k)

                if (k ==  1 .and. j==1) then 
                    write(*,*) u,v,w,omg_x, omg_y, omg_z
                end if

                ! if (omg_x/=omg_x) write(*,*) "omgx nan"
                ! if (omg_y/=omg_y) write(*,*) "omgy nan"
                ! if (omg_z/=omg_z) write(*,*) "omgz nan"

                tmp = (u*omg_x) + (v * omg_y) + (w *omg_z)

                helicity = helicity + tmp
                solver_helicity = solver_helicity + &
                    (rhs_saved(i,j,k,1) * omg_x &
                   + rhs_saved(i,j,k,2) * omg_y &
                   + rhs_saved(i,j,k,3) * omg_z )

               enstrophy = enstrophy + (omg_x**2 + omg_y**2 + omg_z**2) * dx * dy  * dz
            end do
        end do
    end do

    if (helicity /= helicity) then
        write(*,*) "HELICITY IS NAN! This is a hard error - stopping"

        !write(*, *) "u", u
        !write(*, *) "v", v
        !write(*, *) "w", w

        !write(*, *) "omg x", omg_x
        !write(*, *) "omg y", omg_y
        !write(*, *) "omg z", omg_z

        !write(*,*) "u * omgx", u*omg_x
        !write(*,*) "v * omgy", v*omg_y
        !write(*,*) "w * omgz", w*omg_z

        !write(*,*) "tmp", tmp

        stop 1
    end if

    helicity = helicity * dx * dy * dz
    solver_helicity = solver_helicity * dx * dy * dz

end subroutine calculate_helicity

! \int_V  dot(u, div(u)) dV
subroutine advection_dot(udot, omegadot)
    use m_parameters
    use m_data ! d( )/d( ) values from vorticity
    use m_work ! wrk
    implicit none

    integer :: i, j, k
    real*8 :: a, b, c, divu, u, v,w, outer_dot_product
    real*8 :: omg_x, omg_y, omg_z

    real*8 ::udot ! argument, the advection term dotted with u
    real*8 ::omegadot ! argumnet, the advection term dotted with omega
    udot = 0.
    omegadot = 0.

    ! calculate_vortiticity has already been called so we know
    ! dU/dX,dU/dY,dU/dZ,
    ! dV/dX,dV/dY,dV/dZ,
    ! dW/dX,dW/dY,dW/dZ,
    ! have been calculated

    do i=1,nx
        do j=1,ny
            do k=1,nz
                ! \nabla \cdot u
                divu = dUdX(i,j,k) + dVdY(i,j,k) + dWdZ(i,j,k)
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                omg_x = OmgX(i,j,k)
                omg_y= OmgY(i,j,k)
                omg_z= OmgZ(i,j,k)

                ! calculating the whole term
                !outer_dot_product = -0.5 * divu* (u**2+ v**2 + w**2)
                !term = term + (outer_dot_product * dx * dy * dz)

                udot = udot + (-1*(u* dUdX(i,j,k) + v*dUdX(i,j,k) + w*dUdZ(i,j,k)) * dx * dy * dz * u)
                udot = udot + (-1*(u* dVdX(i,j,k) + v*dVdX(i,j,k) + w*dVdZ(i,j,k)) * dx * dy * dz * v)
                udot = udot + (-1*(u* dWdX(i,j,k) + v*dWdX(i,j,k) + w*dWdZ(i,j,k)) * dx * dy * dz * w)

                omegadot= omegadot + (-1*(u* dUdX(i,j,k) + v*dUdX(i,j,k) + w*dUdZ(i,j,k)) * dx * dy * dz * omg_x)
                omegadot= omegadot + (-1*(u* dVdX(i,j,k) + v*dVdX(i,j,k) + w*dVdZ(i,j,k)) * dx * dy * dz * omg_y)
                omegadot= omegadot + (-1*(u* dWdX(i,j,k) + v*dWdX(i,j,k) + w*dWdZ(i,j,k)) * dx * dy * dz * omg_z)
            end do
        end do
    end do
end subroutine advection_dot

! \int_V  dot(p,u) dV
subroutine pressure_dot(udot, omegadot)
    use m_parameters
    use m_work ! wrk
    use m_fields! pressure_field
    use m_data ! omega terms for dot product
    implicit none

    integer :: i, j, k
    real*8 a, b, c, u, v, w

    real*8 :: omg_x, omg_y, omg_z

    real*8 ::udot ! argument, the advection term dotted with u
    real*8 ::omegadot ! argumnet, the advection term dotted with omega
    udot = 0.
    omegadot = 0.

    do i=1,nx
        do j=1,ny
            do k=1,nz
                ! fetch variables from their terms
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                omg_x = OmgX(i,j,k)
                omg_y= OmgY(i,j,k)
                omg_z= OmgZ(i,j,k)

                ! calculate values

                a = u * pressure_field(i,j,k,1)
                b = v * pressure_field(i,j,k,2)
                c = w * pressure_field(i,j,k,3)

                udot = udot + (-1.* (a + b + c) * dx * dy * dz)

                a = omg_x * pressure_field(i,j,k,1)
                b = omg_y * pressure_field(i,j,k,2)
                c = omg_z * pressure_field(i,j,k,3)

                omegadot= omegadot+ (-1.* (a + b + c) * dx * dy * dz)
            end do
        end do
    end do
end subroutine pressure_dot

! nu * \int dot(del^2(u) , u) dV
subroutine diffusion_dot(udot, omegadot)
    use m_parameters
    use m_data ! d( )/d( ) values from vorticity
    use m_work ! wrk
    implicit none

    ! dUdX, dVdY, dWdZ have all ready been calculated in vortitcity

    integer :: i, j, k
    real*8 term, a, b, c, u,v,w
    real*8 :: omg_x, omg_y, omg_z

    real*8 ::udot ! argument, the advection term dotted with u
    real*8 ::omegadot ! argumnet, the advection term dotted with omega
    udot = 0.
    omegadot = 0.

    ! calculate second order gradients
    CALL gradient3D(nx,ny,nz,dUdX,dx,dy,dz,dUdX2,dUdY2,dUdZ2)
    CALL gradient3D(nx,ny,nz,dVdY,dx,dy,dz,dVdX2,dVdY2,dVdZ2)
    CALL gradient3D(nx,ny,nz,dWdZ,dx,dy,dz,dWdX2,dWdY2,dWdZ2)

    ! now we have the second order derivatives

    do i=1,nx
        do j=1,ny
            do k=1,nz
                !
                ! fetch the variables from their containers
                !
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                omg_x = OmgX(i,j,k)
                omg_y= OmgY(i,j,k)
                omg_z= OmgZ(i,j,k)

                !
                ! calculate DP for u
                !
                a = dUdX2(i,j,k) * u
                b = dVdY2(i,j,k) * v
                c = dWdZ2(i,j,k) * w

                udot = udot + (nu * (a + b + c) * dx * dy * dz)

                !
                ! calculate DP for omega
                !
                a = dUdX2(i,j,k) * omg_x
                b = dVdY2(i,j,k) * omg_y
                c = dWdZ2(i,j,k) * omg_z

                omegadot = omegadot + (nu * (a + b + c) * dx * dy * dz)
            end do
        end do
    end do

end subroutine diffusion_dot

subroutine calculate_divu(divu)
    use m_parameters
    use m_data
    use m_work
    implicit none

    real*8 :: divu, u, v, w
    integer:: i, j, k
    divu = 0

    do i = 1,nx
        do j = 1,ny
            do k = 1,nz
                divu = divu + (dUdX(i,j,k) + dVdY(i,j,k) + dWdZ(i,j,k))
            end do
        end do
    end do

    divu = divu * dx * dy * dz

end subroutine calculate_divu

subroutine calculate_vorticity
    use m_work
    use m_data
    implicit none

    real*8, dimension(nx, ny, nz) :: u,v,w
    integer::i

    !write(*,*) "voriticity"

    !u = wrk(:,:,:,1)
    !v = wrk(:,:,:,2)
    !w = wrk(:,:,:,3)

    ! calculate the gradients, store then in the d( )d( ) variables
    CALL gradient3D(nx,ny,nz,wrk(:,:,:,1),dx,dy,dz,dUdX,dUdY,dUdZ)
    CALL gradient3D(nx,ny,nz,wrk(:,:,:,2),dx,dy,dz,dVdX,dVdY,dVdZ)
    CALL gradient3D(nx,ny,nz,wrk(:,:,:,3),dx,dy,dz,dWdX,dWdY,dWdZ)

    ! use the derivatives to calculate the vorticity
    OmgX = dUdY-dVdX
    OmgY = dUdZ-dWdX
    OmgZ = dVdZ-dWdY

    !write(*,*) " from calculate vorticity, outputting the dUdX"
    !do i = 1, nx
        !write(*,*) dUdX(i,1,1), dUdX(i,1,1), dUdX(i,1,1)
    !end do

    !write(*,*) " from calculate vorticity, outputting the omegas"

    !do i = 1, nx
    !    write(*,*) OmgX(i,1,1), OmgY(i,1,1), OmgZ(i,1,1)
    !end do

end subroutine calculate_vorticity

! u . curl (rhs_saved) + \omega \dot rhs_saved
! which is an alternative way to calculate helicity derivative
subroutine alternate_dh_dt_calculation(dh_dt)

    use m_work
    use m_data
    implicit none
    
    integer :: i,j,k
    real*8 :: dh_dt, a, b, c, d,e,f
    real*8 :: omg_x, omg_y, omg_z, u, v,w

    dh_dt = 0.

    ! calculate curl(RHS
    CALL gradient3D(nx,ny,nz,rhs_saved(:,:,:,1),dx,dy,dz,dUdX,dUdY,dUdZ)
    CALL gradient3D(nx,ny,nz,rhs_saved(:,:,:,2),dx,dy,dz,dVdX,dVdY,dVdZ)
    CALL gradient3D(nx,ny,nz,rhs_saved(:,:,:,3),dx,dy,dz,dWdX,dWdY,dWdZ)

    do i =1,nx
        do j=1,ny
            do k=1,nz
                ! pull u/v/w/ omega terms out of the arrays
                u = wrk(i,j,k,1)
                v = wrk(i,j,k,2)
                w = wrk(i,j,k,3)

                omg_x = OmgX(i,j,k)
                omg_y= OmgY(i,j,k)
                omg_z= OmgZ(i,j,k)

                ! this is u. curl(RHS)
                a = dRHS1dX(i,j,k) * u 
                b = dRHS2dY(i,j,k) * v 
                c = dRHS3dZ(i,j,k) * w 

                ! this is \omega \dot rhs_saved
                d = omg_x * rhs_saved(i,j,k,1)
                e = omg_y * rhs_saved(i,j,k,2)
                f = omg_z * rhs_saved(i,j,k,3)

                ! add all the terms together
                dh_dt = dh_dt + a + b + c + d + e + f
                
            end do
        end do
    end do

    dh_dt = dh_dt * dx * dy * dz

end subroutine alternate_dh_dt_calculation

! Subroutine for gradient of a property in 3D domain - modified from MGM's postprocessessing work
subroutine gradient3D(m,n,o,f,hx,hy,hz,dfdx,dfdy,dfdz)
    INTEGER, INTENT(IN) :: m, n, o
    REAL*8, DIMENSION(m,n,o), INTENT(IN) :: f
    REAL*8, INTENT(IN) :: hx, hy, hz
    REAL*8, DIMENSION(m,n,o), INTENT(OUT) :: dfdx, dfdy, dfdz
    INTEGER :: i, j, k

    dfdx = 0.0; dfdy = 0.0; dfdz = 0.0;

    !dfdx
    do k = 1,o
        do j = 1,n
    !        forward difference at start point
            dfdx(1,j,k) = (f(2,j,k)-f(1,j,k)) / hx;
    !        central difference at middle region
            dfdx(2:m-1,j,k) = (f(3:m,j,k)-f(1:m-2,j,k)) / (2.*hx);
    !        backward difference at end point
            dfdx(m,j,k) = (f(m,j,k)-f(m-1,j,k)) / hx;
        end do
    end do

    !dfdy
    do k = 1,o
        do i = 1,m
    !        forward difference at start point
            dfdy(i,1,k) = (f(i,2,k)-f(i,1,k)) / hy;
    !        central difference at middle region
            dfdy(i,2:n-1,k) = (f(i,3:n,k)-f(i,1:n-2,k)) / (2.*hy);
    !        backward difference at end point
            dfdy(i,n,k) = (f(i,n,k)-f(i,n-1,k)) / hy;
        end do
    end do

    !dfdz
    do j = 1,n
        do i = 1,m
    !        forward difference at start point
            dfdz(i,j,1) = (f(i,j,2)-f(i,j,1)) / hz;
    !        central difference at middle region
            dfdz(i,j,2:o-1) = (f(i,j,3:o)-f(i,j,1:o-2)) / (2.*hz);
    !        backward difference at end point
            dfdz(i,j,o) = (f(i,j,o)-f(i,j,o-1)) / hz;
        end do
    end do

end subroutine gradient3D

subroutine collect_mpi_arrays
    use m_data
    use m_work
    implicit none

    integer :: point_count

    point_count = nx * ny * nz
    !
    ! start by sending / collecting wrk arrays to create a velocity array
    ! 


    ! if you are not the master process
    if (myid.ne.master) then
        id_to = master
        tag = myid
        call MPI_SEND(wrk(1:nx,1:ny,1:nz,1:3),point_count,MPI_REAL8,master,tag,MPI_COMM_TASK,mpi_err)

    ! if you ARE the master process then you collect everything together
    else
        ! we are master proc, so first store all of our information in the global array
        ! before its overwritten
        velocity_all(1:nx, 1:ny, 1:nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)

        do id_from=1,numprocs-1
            tag = id_from
            call MPI_RECV(wrk(1:nx,1:ny,1:nz,1:3),point_count,MPI_REAL8,id_from,tag,MPI_COMM_TASK,mpi_status,mpi_err)
            velocity_all(1:nx, 1:ny, (id_from*nz):(id_from+1)*nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)
        end do
    end if

    !
    ! do the same for the pressure arrays
    ! 

    if (myid.ne.master) then
        id_to = master
        tag = myid
        call MPI_SEND(wrk(1:nx,1:ny,1:nz,1:3),point_count,MPI_REAL8,master,tag,MPI_COMM_TASK,mpi_err)

    ! if you ARE the master process then you collect everything together
    else
        ! we are master proc, so first store all of our information in the global array
        ! before its overwritten
        velocity_all(1:nx, 1:ny, 1:nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)

        do id_from=1,numprocs-1
            tag = id_from
            call MPI_RECV(wrk(1:nx,1:ny,1:nz,1:3),point_count,MPI_REAL8,id_from,tag,MPI_COMM_TASK,mpi_status,mpi_err)
            velocity_all(1:nx, 1:ny, (id_from*nz):(id_from+1)*nz, 1:3) = wrk(1:nx, 1:ny, 1:nz, 1:3)
        end do
    end if

endsubroutine collect_mpi_arrays
