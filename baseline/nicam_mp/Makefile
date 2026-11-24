################################################################################
#
# Makefile grand central  (NICAM NSW6 cloud-microphysics mini-app)
#
# Same shape as nicam/NICAM/Makefile: pick the compiler set with NICAM_SYS,
# include Mkinclude, then include one 'dependency' file per source group which
# appends to VPATH / PRJS / LIBS.
#
#   make NICAM_SYS=Linux64-nvidia-ompi-time ENABLE_OPENACC=T
#
################################################################################

TOPDIR    = $(abspath .)
BUILDDIR  = ./workspace
SYSDEPDIR = ./sysdep

VPATH = $(BUILDDIR)

include $(SYSDEPDIR)/Makedef.$(NICAM_SYS)
include $(TOPDIR)/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makelib
	$(MAKE) makebin
	@echo; echo "##### [NICAM-MP] Build and install finished successfully #####"; echo

allclean: clean cleandir

# nhm libraries & programs
include ./nhm/dependency
VPATH +=:$(VPATH_NHM)
PRJS  += $(PRJ_NHM)
LIBS  += $(LIB_NHM)

# object-level dependencies (same format as nicam/NICAM/dependency)
include ./dependency

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)
	mkdir -p $(NICAM_LIBDIR)

makelib: makedir $(LIBS)

makebin: makedir $(PRJS)

nhm: makedir $(PRJ_NHM)

cleandir:
	-rm -f $(NICAM_LIBDIR)/*
	-rm -f $(NICAM_BINDIR)/*

clean:
	rm -f $(LIBS) $(PRJS)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILDDIR); rm -f $(LIBS) $(PRJS); rm -f *~ *.mod *.o *.lst *.L

check:
	$(MAKE) envlist
	$(MAKE) conflist

.SUFFIXES: .o .f .f90

.f.o:
	$(FC) $(FFLAGS) $(MODDIROPT) -o $(OBJDIR)/$@ -c $< $(F77FLAGS)
.f90.o:
	$(FC) $(FFLAGS) $(MODDIROPT) -o $(OBJDIR)/$@ -c $<

%.o: %.mod

$(BUILDDIR)/%.mod: %.f90
	$(make) $(patsubst %.f90,$(BUILDDIR)/%.o,$<)

.PHONY : all clean cleandir allclean nhm check
