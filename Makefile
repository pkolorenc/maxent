#F90=gfortran
#OPTS=-march=native -O3 -Wopenmp $(DEFS)
#OPTS= -g -O0 -Wopenmp -fcheck=all -Wall -fbacktrace $(DEFS)
# 
F90=ifx
OPTS=-xHost -O3 -ipo -qopenmp -assume buffered_io $(DEFS)
#OPTS= -g -O0 -qopenmp -check all -warn all -traceback -debug full $(DEFS)
#
DEFS?=
TARG?=entropy.exe
#
#
#
OBJECTS=numerics_mod.o entropy_mod.o entropyMain.o
#
LIBS=
#
all: $(TARG) readEta.exe

$(TARG): $(OBJECTS) Makefile
	$(F90) $(OPTS) $(DEFS) $(LIBS) $(OBJECTS) -o $(TARG)

readEta.exe: readEtaLast.f90 Makefile
	$(F90) $(OPTS) $(DEFS) $(LIBS) readEtaLast.f90 -o $@

numerics_mod.o: numerics_mod.f90 Makefile
entropy_mod.o: entropy_mod.f90 Makefile
entropyMain.o: entropy_mod.o entropyMain.f90 Makefile

%.o: %.f90
	$(F90) $(OPTS) -c $<

clean:
	rm -f *.o *.mod *__genmod.f90

distclean:
	rm -f *.o *.mod *__genmod.f90 *.exe *~
