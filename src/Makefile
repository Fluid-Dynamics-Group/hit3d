# Makefile with logic in it.  For different hosts runs different things

# Define the correct flags and compilers

MPIF90 = mpif90

# optimized flags - hide the issue
FCFLAGS = -O4 -c -fallow-argument-mismatch
LDFLAGS = -O4 -L/usr/local/lib -lfftw3

# unoptimized flags - issue is illustrated here
FCFLAGS = -O -c -fallow-argument-mismatch -ffree-line-length-none -ffixed-line-length-none -Wall -fcheck=all -g -fbacktrace 
LDFLAGS = -O -L/usr/local/lib -lfftw3

# Program name
PROG    = hit3d.x

# Modules
MODULES = m_openmpi.o\
	m_io.o\
	m_parameters.o\
	m_work.o\
	m_fields.o\
	m_timing.o\
	x_fftw.o\
	m_filter_xfftw.o\
	m_particles.o\
	m_stats.o\
	m_force.o\
	m_rand_knuth.o\
	m_les.o\
	RANDu.o

# Objects
OBJ     = main.o\
	begin_new.o\
	begin_restart.o\
	dealias_all.o\
	get_file_ext.o\
	init_velocity.o\
	init_scalars.o\
	io_write_4.o\
	my_dt.o\
	my_exit.o\
	pressure.o\
	restart_io.o\
	rhs_velocity.o\
	rhs_scalars.o\
	velocity_rescale.o\
	write_tmp4.o


# -------------------------------------------------------
# link

$(PROG):  $(MODULES) $(OBJ)
	$(MPIF90) $(MODULES) $(OBJ) -o $(PROG) $(LDFLAGS)
# -------------------------------------------------------
# compile
#$(MPIF90) $(FCFLAGS) $(FCFLAGS_F77)  $<

$(OBJ): $(MODULES) 

%.o: %.f
	$(MPIF90) $(FCFLAGS) -c $< -o $@

%.o: %.f90
	$(MPIF90) $(FCFLAGS) -c $< -o $@

clean:
	rm *.o *.mod $(PROG)

