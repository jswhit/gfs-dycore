 setenv F77 ifort
 setenv F90 ifort
 setenv W3LIB w3lib-2.0_d
 setenv CFLAGS LINUX
 setenv FINCM 
 setenv ARCHM 
 setenv PGSZM 
 setenv FRRM -FR
 setenv FXXM 
 #setenv OPTSB "-O3 -convert big_endian -xHOST  -fp-model strict -heap-arrays"
 setenv OPTSB "-O3 -convert big_endian -xHOST -recursive"
 setenv OPTSBT "$OPTSB -traceback"
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
 setenv LIBSM 
 #setenv LIBSM "-L${LIBDIR} -lbacio_d -lw3lib-1.9_d -lsp_d ${ESMFLIB}"
