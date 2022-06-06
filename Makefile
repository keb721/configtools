#-*- mode: makefile; mode: font-lock; vc-back-end: RCS -*-
SHELL = /bin/sh

# Where you want the binary
prefix     = $(HOME)
bindir     = $(prefix)/bin

# Compiler and settings
F90       = gfortran
LD        = gfortran
FFLAGS    = -O3 -g
INCLUDE   = 

.PRECIOUS: %.o
.PHONY:  clean

%: %.o
%.o: %.f90
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.f95
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<
%: %.o
	$(F90) $(FFLAGS) $(INCLUDE)  -o $@ $^

all :  tip4p xmol2config

TIP4P_OBJECTS = tip4p.o

tip4p : $(TIP4P_OBJECTS)
	$(LD) -o $(bindir)/tip4p $(TIP4P_OBJECTS)


X2C_OBJECTS = xmol2config.o

xmol2config: $(X2C_OBJECTS)
	$(LD) -o $(bindir)/xmol2config $(X2C_OBJECTS)


clean : 

	rm -f *.mod *.d *.il *.o work.*
	rm -f $(bindir)/tip4p
	rm -f $(bindir)/xmol2config


