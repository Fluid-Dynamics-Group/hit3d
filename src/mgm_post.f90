!======================MGM's hit3d Postprocessing script===================
!Function:       Save u, v & w in my (organized .dat, ascii & tecplot) formats
!Description:    Save the resutls of hit3d (u, v & w) in proper organized manner of .dat, ascii & tecplot formats
!Details:        save data
!Date created:   07/06/2018
!Date edited:    07/09/2018; cleaned-up


program MgmTurbNetTools

implicit none

!========================
integer*4 :: nx, ny, nz, npoints, nscalars
integer*4 :: sizes(3)
real :: dx, dy, dz
real, parameter :: pi = 3.141592653589793
real, dimension(:), allocatable :: xx_v, yy_v, zz_v, XX, YY, ZZ
logical :: readforce
INTEGER :: tclock
REAL(4)    :: time, ttime
character(6) :: boxchar
character(2) :: boxchar1
character*80 :: fname
!real, dimension(:,:,:), allocatable :: u, v, w
real*4, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:), phi(:,:,:)
real, dimension(:,:,:), allocatable :: dUdX, dUdY, dUdZ, dVdX, dVdY, dVdZ, dWdX, dWdY, dWdZ, omgmag, strmag
real, dimension(:,:,:), allocatable :: OmgX, OmgY, OmgZ, omg1, omg2, omg3, gam1, gam2, gam3, sig1, sig2, sig3, qcrit
real, dimension(:), allocatable :: u_vec, v_vec, w_vec, OmgX_vec, OmgY_vec, OmgZ_vec, qcrit_vec, omgmag_vec, strmag_vec
real, dimension(:), allocatable :: phi_vec

real*4 :: a,b,c, energy

integer, parameter :: RP=4 ! Number of bytes for reals (single precision)

real(RP), dimension(:,:), allocatable :: points    ! input

! loop iterator
integer :: i, j , k, l, p, q, r, boxnum, itr

!============Load data====================================
!itr = 1000
write(*,*) 'Enter iteration number = '
read(*,*) itr
write(boxchar,"(I6.6)") itr

fname = "output/velocity/u."//boxchar
open(unit=100,file=fname,form="unformatted",access='stream')
read(100) sizes
nx = sizes(1); ny = sizes(2); nz = sizes(3);
npoints = nx*ny*nz;
print *, 'Sizes: ', nx, ny, nz
allocate(u(nx,ny,nz)); allocate(v(nx,ny,nz)); allocate(w(nx,ny,nz))
read(100) ( ( (u(i,j,k),i=1,nx) ,j=1,ny) ,k= 1,nz )   ! +(fact*(p-1))
close(100)

fname = "output/velocity/v."//boxchar
open(unit=100,file=fname,form="unformatted",access='stream')
read(100)
read(100) ( ( (v(i,j,k),i=1,nx) ,j=1,ny) ,k= 1,nz )
close(100)
print *, 'v(2,1,1) = ', v(2,1,1)

fname = "output/velocity/w."//boxchar
open(unit=100,file=fname,form="unformatted",access='stream')
read(100)
read(100) ( ( (w(i,j,k),i=1,nx) ,j=1,ny) ,k= 1,nz )
close(100)


write(*,*) 'Enter numer of passive scalars = '
read(*,*) nscalars
write(boxchar1,"(I2.2)") nscalars

fname = "output/sc"//boxchar1//"."//boxchar
allocate(phi(nx,ny,nz));
open(unit=100,file=fname,form="unformatted",access='stream')
read(100)
read(100) ( ( (phi(i,j,k),i=1,nx) ,j=1,ny) ,k= 1,nz )
close(100)

print *, 'Data loaded...'

!= ====


!=========Post-process data=================================
allocate(xx_v(nx)); allocate(yy_v(ny)); allocate(zz_v(nz));
allocate(XX(npoints)); allocate(YY(npoints)); allocate(ZZ(npoints));
allocate(points(3,npoints));
allocate(u_vec(npoints)); allocate(v_vec(npoints)); allocate(w_vec(npoints));
allocate(phi_vec(npoints));

dx = (2.*pi) / nx; dy = dx; dz = dx;

!Vorticity calculation
allocate(dUdX(nx,ny,nz));allocate(dUdY(nx,ny,nz));allocate(dUdZ(nx,ny,nz));
allocate(dVdX(nx,ny,nz));allocate(dVdY(nx,ny,nz));allocate(dVdZ(nx,ny,nz));
allocate(dWdX(nx,ny,nz));allocate(dWdY(nx,ny,nz));allocate(dWdZ(nx,ny,nz));
allocate(OmgX(nx,ny,nz));allocate(OmgY(nx,ny,nz));allocate(OmgZ(nx,ny,nz));
allocate(omg1(nx,ny,nz));allocate(omg2(nx,ny,nz));allocate(omg3(nx,ny,nz));
allocate(gam1(nx,ny,nz));allocate(gam2(nx,ny,nz));allocate(gam3(nx,ny,nz));
allocate(sig1(nx,ny,nz));allocate(sig2(nx,ny,nz));allocate(sig3(nx,ny,nz));
allocate(qcrit(nx,ny,nz)); allocate(omgmag(nx,ny,nz)); allocate(strmag(nx,ny,nz));
CALL gradient3D(nx,ny,nz,u,dx,dy,dz,dUdX,dUdY,dUdZ)
CALL gradient3D(nx,ny,nz,v,dx,dy,dz,dVdX,dVdY,dVdZ)
CALL gradient3D(nx,ny,nz,w,dx,dy,dz,dWdX,dWdY,dWdZ)

allocate(OmgX_vec(npoints));allocate(OmgY_vec(npoints));allocate(OmgZ_vec(npoints));
allocate(qcrit_vec(npoints)); allocate(omgmag_vec(npoints)); allocate(strmag_vec(npoints));

energy = 0.

do i =1,nx
    do j=1,ny
        do k=1,nz
            a = u(i,j,k)
            b = v(i,j,k)
            c = w(i,j,k)
            !write(*,*) a*a + b*b + c*c
            energy = energy + ((a**2 + b**2 + c**2) * dx * dy* dz)

            !write(*,*) energy

        end do
    end do
end do

write(*,*) "energy is ", energy

OmgX = dUdY-dVdX
OmgY = dUdZ-dWdX
OmgZ = dVdZ-dWdY

!Q-criterion
omg1 = -0.5 * (dUdY-dVdX)
omg2 = -0.5 * (dUdZ-dWdX)
omg3 = -0.5 * (dVdZ-dWdY)
gam1 = 0.5 * (dVdX+dUdY)
gam2 = 0.5 * (dWdX+dUdZ)
gam3 = 0.5 * (dWdY+dVdZ)
sig1 = dUdX
sig2 = dVdY
sig3 = dWdZ
qcrit = (omg1**2+omg2**2+omg3**2) - ((gam1**2+gam2**2+gam3**2) + 0.5*(sig1**2+sig2**2+sig3**2))
omgmag = omg1**2+omg2**2+omg3**2
strmag = (gam1**2+gam2**2+gam3**2) + 0.5*(sig1**2+sig2**2+sig3**2)
print *, 'Vorticity & Q-criterion calculated...'

!Mesh generation
xx_v(1) = 0.0; yy_v(1) = 0.0; zz_v(1) = 0.0;
do i = 2, nx
    xx_v(i) = xx_v(i-1) + dx
    yy_v(i) = yy_v(i-1) + dy
    zz_v(i) = zz_v(i-1) + dz
end do
l = 1
do k = 1, nz
    do j = 1, ny
        do i = 1, nx
            XX(l) = xx_v(i)
            YY(l) = yy_v(j)
            ZZ(l) = zz_v(k)
            u_vec(l) = u(i,j,k)
            v_vec(l) = v(i,j,k)
            w_vec(l) = w(i,j,k)
            phi_vec(l) = phi(i,j,k)
            OmgX_vec(l) = OmgX(i,j,k)
            OmgY_vec(l) = OmgY(i,j,k)
            OmgZ_vec(l) = OmgZ(i,j,k)
            qcrit_vec(l) = qcrit(i,j,k)
            omgmag_vec(l) = omgmag(i,j,k)
            strmag_vec(l) = strmag(i,j,k)
            l = l + 1
        end do
    end do
end do
points(1, :) = XX
points(2, :) = YY
points(3, :) = ZZ
print *, 'Cartesian mesh generated...'

!===============Write data (Binary)========================
inquire(file="post/veldata_"//boxchar//".dat",exist=readforce)
if (readforce) THEN
    open(unit=100,file="post/veldata_"//boxchar//".dat",&
        action="write",form="unformatted",status="replace")
else
    open(unit=100,file="post/veldata_"//boxchar//".dat",&
        action="write",form="unformatted",status="new")
end if
write(100) points
write(100) u_vec
write(100) v_vec
write(100) w_vec
write(100) OmgX_vec
write(100) OmgY_vec
write(100) OmgZ_vec
write(100) qcrit_vec
write(100) omgmag_vec
write(100) strmag_vec
write(100) phi_vec
close(100)
write(*,*) 'Data written (binary)'

!===============Write data (ASCII)========================
!inquire(file="post/veldata_"//boxchar//"ascii.dat",exist=readforce)
!if (readforce) THEN
!    open(unit=100,file="post/veldata_"//boxchar//"ascii.dat",&
!        action="write",form="formatted",status="replace")
!else
!    open(unit=100,file="post/veldata_"//boxchar//"ascii.dat",&
!        action="write",form="formatted",status="new")
!end if
!do i = 1,npoints
!    write(100,*) points(1,i), points(2,i), points(3,i),&
!    u_vec(i), v_vec(i), w_vec(i)
!end do
!close(100)
!write(*,*) 'Data written (ascii)'


!==================TECPLOT FILE==================
inquire(file="post/figs/veldata_"//boxchar//"plt.dat",exist=readforce)

if (readforce) THEN
    open(unit=100,file="post/figs/veldata_"//boxchar//"plt.dat",&
    action="write",form="formatted",status="replace")
else
    open(unit=100,file="post/figs/veldata_"//boxchar//"plt.dat",&
    action="write",form="formatted",status="new")
end if

WRITE(100,*) 'TITLE     = "3D Turbulent flow field (hit3d)"'
WRITE(100,*) 'VARIABLES = "X"'
WRITE(100,*) '"Y"'
WRITE(100,*) '"Z"'
WRITE(100,*) '"U"'
WRITE(100,*) '"V"'
WRITE(100,*) '"W"'
WRITE(100,*) '"VortX"'
WRITE(100,*) '"VortY"'
WRITE(100,*) '"VortZ"'
WRITE(100,*) '"Q"'
WRITE(100,*) '"Str"'
WRITE(100,*) '"Phi"'
WRITE(100,*) 'ZONE T="Rectangular zone"'
WRITE(100,*) 'I=',nx,', J=',ny,', K=',nz, ', ZONETYPE=Ordered'
WRITE(100,*) 'DATAPACKING=POINT'
WRITE(100,*) 'DT=(SINGLE SINGLE SINGLE)'

do i = 1,npoints
    write(100,*) points(1,i), points(2,i), points(3,i),&
    u_vec(i), v_vec(i), w_vec(i),&
    OmgX_vec(i), OmgY_vec(i), OmgZ_vec(i), qcrit_vec(i), strmag_vec(i),&
    phi_vec(i)
end do
close(100)
write(*,*) 'Figure plotted (TECPLOT)'


!============SUBROUTINES & FUNCTIONS=========================
CONTAINS
! Subroutine for clock start
SUBROUTINE tick(t)
INTEGER, INTENT(OUT) :: t

CALL system_clock(t)
END SUBROUTINE tick
! Function to return end time
REAL FUNCTION tock(t)
INTEGER, INTENT(IN) :: t
INTEGER :: now, clock_rate

call system_clock(now,clock_rate)
tock = real(now - t)/real(clock_rate)
END FUNCTION tock
! Function to perform cross product of 2 vectors
FUNCTION cross(a,b,n)
INTEGER, INTENT(IN) :: n
REAL, DIMENSION(n,3) :: cross
REAL, INTENT(IN), DIMENSION(n,3) :: a, b

cross(:,1) = a(:,2) * b(:,3) - a(:,3) * b(:,2)
cross(:,2) = a(:,3) * b(:,1) - a(:,1) * b(:,3)
cross(:,3) = a(:,1) * b(:,2) - a(:,2) * b(:,1)
END FUNCTION cross

!===========================
! Subroutine for gradient of a property in 3D domain
SUBROUTINE gradient3D(m,n,o,f,hx,hy,hz,dfdx,dfdy,dfdz)
INTEGER, INTENT(IN) :: m, n, o
REAL*4, DIMENSION(m,n,o), INTENT(IN) :: f
REAL, INTENT(IN) :: hx, hy, hz
REAL, DIMENSION(m,n,o), INTENT(OUT) :: dfdx, dfdy, dfdz
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

END SUBROUTINE gradient3D


end program MgmTurbNetTools

