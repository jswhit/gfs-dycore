for file in *.f; do
   file2=`basename $file`
   #/bin/cp -f /Volumes/Data1/emcsvn/gfs/global_fcst.fd/$file2 orig
   #diff -urN /Volumes/Data1/emcsvn/gfs/global_fcst.fd/$file2 $file
   diff -urN orig/$file $file
done
