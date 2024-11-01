include ../../Makefile.in

F90FILES = planeStats.f90 planeStats_metric.f90 planeStats_mpi.f90 planeStats_average.f90 planeStats_budgets.f90 planeStats_budgets_velocity.f90 planeStats_budgets_scalar.f90 planeStats_scatter.f90 planeStats_cond.f90 planeStats_finitechem.f90

BINFILE = planeStats

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)
MECHFILE = HYDROGENF.o

.SUFFIXES: .f90 .o

$(BINFILE): $(OFILES) $(MECHFILE) $(LIBDIR)/liblibrary.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
