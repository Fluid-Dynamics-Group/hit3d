subroutine write_velocity_field
    use m_fields     ! fields
    use m_parameters ! nx, ny, nz
    implicit none

    integer :: filenumber, meta_file, i, j, k
    real*8 :: u, v,w
    character(len=30) :: filename, sizefile

    write(*,*) "writing velocity field out"

    filenumber = 619

    write(filename, "('output/velocity_field/', i2.2, '.csv')") myid_world
    write(*,*) filename
    open(filenumber, file=filename)

    write(filenumber, "(A5)") "u,v,w"

    do k = 1,nz
      do j = 1,ny
          do i = 1,nx
              u = fields(i,j,k,1)
              v = fields(i,j,k,2)
              w = fields(i,j,k,3)
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
