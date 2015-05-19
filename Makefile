IDIR =../include
CC=g++ -g
CFLAGS=-I$(IDIR) 
LFLAGS=-lm -DAXIAL

ODIR=obj
LDIR =../lib

_DEPS = md_sim.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =  md_sim_HS.o in_out.o ran2.o initial.o tools.o strucfac.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

hellomake: $(OBJ)
	$(CC) -o pp $^ $(CFLAGS) $(LFLAGS) 
	

	cp pp ../../../Testing/V_0.6.0/
	tar -zcf ../../code.tgz ../
	cp -f ../../code.tgz ../../../Testing/V_0.6.0/code.tgz

.PHONY: clean
	
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

