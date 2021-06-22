subroutine write_velocity_field(current_timestep)
    use m_work ! fields
    use m_parameters ! nx, ny, nz

    implicit none

    integer :: filenumber, meta_file, i, j, k, current_timestep
    !real*8, allocatable :: fields(:,:,:,:)
    real*8 :: u, v,w

    character(len=60) :: filename, sizefile

    write(*,*) "writing velocity field out"

    filenumber = 619

    write(filename, "('output/velocity_field/', i2.2, '_', i5.5, '.csv')") myid_world, current_timestep
    write(*,*) filename
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
