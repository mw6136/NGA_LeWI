include ../../Makefile.in

F90FILES = volumeStats.f90 volumeStats_metric.f90 volumeStats_finitechem.f90 volumeStats_mpi.f90 volumeStats_average.f90 volumeStats_budgets.f90 volumeStats_velocity.f90 volumeStats_les.f90

BINFILE = volumeStats

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)
MECHFILE = HYDROGENF.o

.SUFFIXES: .f90 .o

$(BINFILE): $(OFILES) $(MECHFILE) $(LIBDIR)/liblibrary.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS) $(LAPACK_LIB) $(BLAS_LIB)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
