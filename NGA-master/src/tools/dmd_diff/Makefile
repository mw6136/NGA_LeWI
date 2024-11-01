include ../../Makefile.in

F90FILES = dmd_diff.f90

F77FILES = fdump.f i1mach.f j4save.f xercnt.f xerhlt.f xermsg.f xerprn.f xersve.f xgetua.f dsort.f

BINFILE = dmd_diff

OFILES   = $(F77FILES:.f=.o) $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .f .o

$(BINFILE): $(OFILES) $(LIBDIR)/liblibrary.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LAPACK_LIB) $(BLAS_LIB) $(FFTW_LIB) $(LDFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)

.f.o:
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f -o $(OBJDIR)/$*.o $(MODFLAGS)
