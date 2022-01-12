include SLmake.inc.intel

F90 = $(F77)

PMLIBS = $(BLACSLIB) $(BLASLIB)

#TESTDRV2 = OneTo2D_DoubleDataRedist5.o test_DataRedistDouble_timing.o 
#TESTDRV2 = test_dataredist.o 
#TESTDRV2 = TwoTo2D_DCplxRedistBLACS.o test_dataRedist2d.o
#TESTDRV2 = cardmpp.o OneTo1D_ComplexDataRedistBip.o test_DataRedist1DBip.o
TESTDRV2 = OneTo1D_ComplexDataRedist.o test_DataRedist1D.o

all: test2

test2: $(TESTDRV2) 
	$(F77LOADER) $(F77LOADFLAGS) -o $@ $(TESTDRV2) $(PMLIBS)
	mv $@ ../bin

%.o:%.f90
	$(F90) -c $(F77FLAGS) $<

%.o:%.f
	$(F90) -c $(F77FLAGS) $<

clean :
	rm -f *.o *.mod
	rm test1 test test2
