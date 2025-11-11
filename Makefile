# Makefile for FULiou 90 simple version
FC = ifort 
FFLAGS = -cpp -w -O2 -auto -noalign -convert big_endian

FUCODE =  4S-FL.o Composit-FL.o gas-FL.o \
plank-others-FL.o rayley-FL.o thicks-FL.o rad-FL.o aero_cal-FL.o\
rain-FL.o ice-water-FL.o

all : FU-wang  FU-wang.o ${FUCODE}
#RadParams_flcode.o : RadParams_flcode.o
#	${FC} -c ${FFLAGS} RadParams_flcode.f90

4S-FL.o :  4S-FL.f90
	${FC} -c ${FFLAGS} 4S-FL.f90
Composit-FL.o :  Composit-FL.f90
	${FC} -c ${FFLAGS} Composit-FL.f90
gas-FL.o :  gas-FL.f90
	${FC} -c ${FFLAGS} gas-FL.f90
ice-water-FL.o:  ice-water-FL.f90
	${FC} -c ${FFLAGS} ice-water-FL.f90
plank-others-FL.o :  plank-others-FL.f90
	${FC} -c ${FFLAGS} plank-others-FL.f90
rain-FL.o :  rain-FL.f90
	${FC} -c ${FFLAGS} rain-FL.f90
rayley-FL.o :  rayley-FL.f90
	${FC} -c ${FFLAGS} rayley-FL.f90
thicks-FL.o:  thicks-FL.f90
	${FC} -c ${FFLAGS} thicks-FL.f90
aero_cal-FL.o:  aero_cal-FL.f90
	${FC} -c ${FFLAGS} aero_cal-FL.f90
rad-FL.o:  rad-FL.f90
	${FC} -c ${FFLAGS} rad-FL.f90
FU-wang.o : FU-wang.f90
	${FC} -c ${FFLAGS} FU-wang.f90
FU-wang : FU-wang.o ${FUCODE} 
	${FC} ${FFLAGS} FU-wang.f90 ${FUCODE} -o FU-wang
clean:
	rm  *.o 
