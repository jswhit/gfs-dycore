gfortran -O2 -fPIC -c ../src/shtns.f90
gfortran -O2 -fPIC -c ../src/kinds.f90
gfortran -O2 -fPIC -c ../src/sigio_module.f90
SHTNSDIR=/Users/${USER}
FFTWDIR=/opt/local
f2py -c read_sigma.f90 -m read_sigma --fcompiler=gnu95 shtns.o kinds.o sigio_module.o -L${SHTNSDIR}/lib -lshtns -L${FFTWDIR}/lib -lfftw3
