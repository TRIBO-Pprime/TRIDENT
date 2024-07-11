#!/bin/sh

cd ../TOOLIB
cd    UTILS ; make all
cd ../QSORT ; make all
cd ../CHOLE ; make all
cd ../ALGEN ; make all
cd ../LEAST ; make all
cd ../MINIM ; make all
cd ../SPLIN ; make all
cd ../TCHEV ; make all
cd ../DIGIS ; make all
cd ../FILES ; make all
cd ../GPLOT ; make all
cd ../FFTW3 ; make all
cd ../INTPL ; make all

cd ../TPGLIB
cd    STATS ; make all
cd ../ASFC2 ; make all
cd ../FILTR ; make all
cd ../ABBOT ; make all
cd ../DERIV ; make all
cd ../ANISO ; make all
cd ../MORPH ; make all

cd ../TRIDENT
make all

echo "OK"
read wait

