F90  = gfortran -static-libgfortran -O3
MPIF90 = ftn -DMPI -O3
OBJS_SQUARE  =  sstates-square.o lattice-square.o parameters.o constants.o utility.o

OBJS_HEX = sstates-hex.o lattice-hex.o parameters.o constants.o utility.o

OBJS_CUBIC = sstates-cubic.o lattice-cubic.o parameters.o constants.o utility.o

OBJS_TETRA = sstates-tetra.o lattice-tetra.o parameters.o constants.o utility.o

OBJS_J2_SQUARE = sstates-j2-square.o lattice-j2-square.o parameters.o constants.o utility.o

OBJS_DISC_CYLINDER = sstates-disc-cylinder.o lattice-disc-cylinder.o parameters-3D.o utility.o constants.o

OBJS_AC_CYLINDER = sstates-ac-cylinder.o lattice-ac-cylinder.o parameters-ac.o utility.o constants.o

default: square

square: main.f90 $(OBJS_SQUARE)
	$(F90) main.f90 $(OBJS_SQUARE) -o spinmc-square.x

hex: main.f90 $(OBJS_HEX)
	$(F90) main.f90 $(OBJS_HEX) -o spinmc-hex.x

cubic: main.f90 $(OBJS_CUBIC)
	$(F90) main.f90 $(OBJS_CUBIC) -o spinmc-cubic.x

tetra: main.f90 $(OBJS_TETRA)
	$(F90) main.f90 $(OBJS_TETRA) -o spinmc-tetra.x

j2-square: main.f90 $(OBJS_J2_SQUARE)
	$(F90) main.f90 $(OBJS_J2_SQUARE) -o spinmc-j2-square.x

disc-cylinder: main.f90 $(OBJS_DISC_CYLINDER)
	$(F90) main.f90 $(OBJS_DISC_CYLINDER) -o spinmc-disc-cylinder.x

ac-cylinder: main.f90 $(OBJS_AC_CYLINDER)
	$(F90) main.f90 $(OBJS_AC_CYLINDER) -o spinmc-ac-cylinder.x


######## tmp use
snapshot-q: snapshot-q.f90 $(OBJS_SQUARE)
	$(F90) snapshot-q.f90 $(OBJS_SQUARE) -o spinmc-snapshot-q.x
########

clean:
	rm -f *.o *.mod *.x

######## parellal shell
mpi: mpilink.f90
	$(MPIF90) mpilink.f90 -o spinmc_mpi.x
########

######## source

sstates-ac-cylinder.o: sstates-ac.f90 lattice-ac-cylinder.o parameters-ac.o utility.o constants.o
	$(F90) -c sstates-ac.f90 -o sstates-ac-cylinder.o

lattice-ac-cylinder.o: lattice-ac-cylinder.f90 parameters-ac.o utility.o constants.o
	$(F90) -c lattice-ac-cylinder.f90

parameters-ac.o: parameters-ac.f90 utility.o constants.o
	$(F90) -c parameters-ac.f90

sstates-disc-cylinder.o: sstates-disc.f90 lattice-disc-cylinder.o parameters-3D.o utility.o constants.o
	$(F90) -c sstates-disc.f90 -o sstates-disc-cylinder.o

lattice-disc-cylinder.o: lattice-disc-cylinder.f90 parameters-3D.o utility.o constants.o
	$(F90) -c lattice-disc-cylinder.f90

sstates-j2-square.o: sstates-j2.f90 lattice-j2-square.o parameters.o utility.o constants.o
	$(F90) -c sstates-j2.f90 -o sstates-j2-square.o

lattice-j2-square.o: lattice-j2-square.f90 parameters.o utility.o constants.o
	$(F90) -c lattice-j2-square.f90

sstates-cubic.o: sstates.f90 lattice-cubic.o parameters.o utility.o constants.o
	$(F90) -c sstates.f90 -o sstates-cubic.o

lattice-cubic.o: lattice-cubic.f90 parameters.o utility.o constants.o
	$(F90) -c lattice-cubic.f90

parameters-3D.o: parameters-3D.f90 constants.o
	$(F90) -c parameters-3D.f90

sstates-square.o: sstates.f90 lattice-square.o parameters.o utility.o constants.o
	$(F90) -c sstates.f90 -o sstates-square.o

lattice-square.o: lattice-square.f90 parameters.o utility.o constants.o
	$(F90) -c lattice-square.f90

sstates-tetra.o: sstates.f90 lattice-tetra.o parameters.o utility.o constants.o
	$(F90) -c sstates.f90 -o sstates-tetra.o

lattice-tetra.o: lattice-tetra.f90 parameters.o utility.o constants.o
	$(F90) -c lattice-tetra.f90

sstates-hex.o: sstates.f90 lattice-hex.o parameters.o utility.o constants.o
	$(F90) -c sstates.f90 -o sstates-hex.o

lattice-hex.o: lattice-hex.f90 parameters.o utility.o constants.o
	$(F90) -c lattice-hex.f90

parameters.o: parameters.f90 constants.o
	$(F90) -c parameters.f90

utility.o: utility.f90 constants.o
	$(F90) -c utility.f90

constants.o: constants.f90
	$(F90) -c constants.f90
########
