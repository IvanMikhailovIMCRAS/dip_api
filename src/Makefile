MF=      Makefile
FC=      gfortran
FFLAGS=  -O3 -fopenmp
LFLAGS=  $(FFLAGS)

# name of EXEcutable file:
EXE=     ../dip.exe

# list of SouRCe files:
SRC= \
       CommonParam.f90  \
       ErrorList.f90 \
       Lib.f90 \
       InputOutput.f90 \
       OnePoint.f90  \
       DIP.f90

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)
.f90.o:
	$(FC) $(FFLAGS) -c $<
all:	$(EXE)
$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)
$(OBJ):	$(MF)
tar:
	tar cvf $(EXE).tar $(MF) $(SRC)
clean:
	rm -f $(OBJ) $(EXE) *.mod core