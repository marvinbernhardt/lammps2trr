PREFIX = $$HOME
OBJ = lammps2trr.c
PROGRAM = lammps2trr

CFLAGS = -Wall -O3 -std=gnu99
LDFLAGS = -lgromacs -llapacke -llapack -lcblas -lblas -lfftw3f -lm -lgfortran
INC = -I $${GROMACS_DIR}/include
LIB = -L $${GMXLDLIB}

$(PROGRAM): $(OBJ) verbPrintf.c
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJ) $(INC) $(LIB) $(LDFLAGS)

.PHONY: install
install: $(PROGRAM)
	mkdir -p $(PREFIX)/bin
	cp $(PROGRAM) $(PREFIX)/bin/$(PROGRAM)

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(PROGRAM)

.PHONY: clean
clean:
	rm -f $(PROGRAM)
