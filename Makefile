PROCESSOR := $(shell uname -m)

#ifeq ($(PROCESSOR),ia64)
  F90=gfortran
  FFLAGS=-g -C -O0 -ffree-form -I/opt/local/include #-fcheck=all
  FFLAGS2=$(FFLAGS)
  LDFLAGS=-L/opt/local/lib -lnetcdf -lnetcdff
#else
#  ifeq ($(PROCESSOR),x86_64)
#    F90=/usr/local/pgi/linux86-64/5.2/bin/pgf90
#    FFLAGS=-Kieee -fastsse #-Kieee # 
#    LDFLAGS=-lnetcdf
#  else
#    F90=lf95
#    FFLAGS=-g --chk aesux --trap -I/usr/local/include #--o2  -I/usr/local/include #
#    LDFLAGS=-lnetcdf
#  endif
#endif

SOURCES= DGmod.f90 dgsweep.f90 ppmsweep_fct.f90 ppmsweep_pmod.f90 ppmwrap.f90 ppmsweep_pmod_clean_sl.f90 pcmsweep_semilagrangian.f90 tfcn.f90 
OBJECTS=$(SOURCES:.f90=.o)

all: $(SOURCES) test_skam_1d test_skam_2d

test_skam_1d: $(OBJECTS) test_advection_skam_1d.f90
	$(F90) $(FFLAGS) $(OBJECTS) test_advection_skam_1d.f90 -o $@ $(LDFLAGS) 

test_skam_2d: $(OBJECTS) test_advection_slskam_2d.f90
	$(F90) $(FFLAGS) $(OBJECTS) test_advection_slskam_2d.f90 \
	         -o $@ $(LDFLAGS) 

.PHONY:

clean:
	rm -f *.o test_skam_1d test_skam_2d

%.o : %.f90
	$(F90) -c $(FFLAGS) $<
