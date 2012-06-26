 setenv F77 gfortran
 setenv F90 gfortran
 setenv LIBDIR /lfs1/projects/globpsd/whitaker/EXP-hybens/nwprod/lib
 setenv W3LIB w3lib-2.0_d
 setenv CFLAGS LINUX
 setenv FINCM -I$LIBDIR/incmod/w31b-2.0_d
 setenv ARCHM 
 setenv PGSZM 
 setenv FRRM -ffree-form
 setenv FXXM 
 #setenv OPTSB "-O3 -convert big_endian -xHOST  -fp-model strict -heap-arrays"
 setenv OPTSB  "-O3 -march=native -frecursive -ffast-math -pipe -fomit-frame-pointer -fno-range-check"
 #setenv OPTSB "-O0 -pedantic -fimplicit-none -fbounds-check -fbacktrace -Wall -fcheck-array-temporaries -fno-range-check -g"
 #setenv OPTSB "-O0 -ffpe-trap=invalid,zero,overflow -fimplicit-none -fbounds-check -fbacktrace -Wall -fcheck-array-temporaries -fno-range-check -g"
 setenv OPTSBT "$OPTSB -fbacktrace"
 #setenv OPTSM "$OPTSBT -r8"
 #setenv OPTSIOM "$OPTSBT -r8 "
 #setenv OPTS_SERM "$OPTSBT -r8 $ARCHM"
 #setenv OPTS90M "$OPTSBT   -r8 $FRRM"
 #setenv OPTS90AM "$OPTSBT  -r8 $FRRM"
 setenv OPTSM "$OPTSBT"
 setenv OPTSIOM "$OPTSBT"
 setenv OPTS_SERM "$OPTSBT $ARCHM"
 setenv OPTS90M "$OPTSBT   $FRRM"
 setenv OPTS90AM "$OPTSBT  $FRRM"
 setenv LDFLAGSM ""
 setenv ESMFINC ""
# I had to add an explicit link to libmpichcxx on our system
# to avoid unresolved references to MPICH c++ symbols.
 setenv ESMFLIB ""
 setenv LDRM mpif90
 setenv LIBSM "-L${LIBDIR} -lbacio_d -lw3-2.0_d -lsp_d ${ESMFLIB}"
 #setenv LIBSM "-L${LIBDIR} -lbacio_d -lw3lib-1.9_d -lsp_d ${ESMFLIB}"
