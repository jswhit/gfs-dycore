#!/bin/tcsh 
##PBS -q batch
##PBS -l procs=12
##PBS -l walltime=8:00:00
##PBS -A fim
##PBS -N testgfs
##PBS -S /bin/tcsh
#$ -j y
#$ -cwd
#$ -l h_rt=06:00:00
#$ -A gfsenkf  
#$ -pe hfip 12  
#$ -N testgfs
#$ -o testgfs.out
#set date=2012060100
set homedir=$PWD
source analdate.csh
set datout=/pan2/projects/gfsenkf/whitaker/test/${analdate}
set incdate="/whome/whitaker/bin/incdate"
mkdir -p $datout
cd $datout
ln -fs /pan2/projects/gfsenkf/whitaker/test/init/gec00.${analdate}.sanl sanl.dat
ln -fs /pan2/projects/gfsenkf/whitaker/test/init/gec00.${analdate}.sfcanl sfcanl.dat
setenv HOMEGLOBAL /lfs1/projects/globpsd/whitaker/EXP-hybens
setenv EXECGLOBAL ${HOMEGLOBAL}/bin
setenv FIXGLOBAL ${HOMEGLOBAL}/gfs/fix_am
setenv FHMAX 120
setenv FHOUT 24
setenv JCAP 254
setenv LEVS 42
setenv LATB 384
setenv LONB 768
ln -fs ${FIXGLOBAL}/global_o3prdlos.f77 o3forcing.dat
ln -fs ${FIXGLOBAL}/global_mtnvar.t${JCAP}.${LONB}.${LATB}.f77 mtnvar.dat
ln -fs ${FIXGLOBAL}/global_orography.t${JCAP}.${LONB}.${LATB}.grb orography
ln -fs ${FIXGLOBAL}/global_orography_uf.t${JCAP}.${LONB}.${LATB}.grb orography_uf
ln -fs ${FIXGLOBAL}/global_climaeropac_global.txt     aerosol.dat
ln -fs ${FIXGLOBAL}/global_solarconstantdata.txt solarconstantdata.txt
set volcfiles=`ls $FIXGLOBAL | grep volcanic_aerosols`
foreach file ($volcfiles)
  set file2=`echo $file |sed -e "s/global_//g"`
  /bin/cp -f $FIXGLOBAL/$file $file2
end
set co2files=`ls $FIXGLOBAL | grep co2historicaldata`
foreach file ($co2files)
  set file2=`echo $file |sed -e "s/global_//g"`
  /bin/cp -f $FIXGLOBAL/$file $file2
end

cat  > gfs.nml <<EOF
 &nam_mrf
 initfile="sanl.dat",
 sfcinitfile="sfcanl.dat",
 fhmax=${FHMAX},
 fhout=${FHOUT},
 iaer=111,
 ico2=1,
 fhzer=${FHOUT},
 fhlwr=1,
 fhswr=1,
 timestepsperhr=6,
 explicit=.false.,
 postphys=.true.,
 ntrac=3,
 flgmin=0.220,
 ccnorm=.false.,
 ctei_rm=0.50,
 ICTM=1,
 /
 &soil_veg
  LPARAM = .FALSE.
 /
 &END
EOF

setenv OMP_NUM_THREADS 12
/bin/rm -f sig.f* sfc.f* flx.f* pgrb.f*
#export FCSTVARS="ras=.false.,nsout=0,lsm=1,tfiltc=0.85,liope=.true.,zhao_mic=.true.,ncw=50,150,crtrh=0.85,0.85,0.85,flgmin=0.220 ,IALB=$ialb,ccnorm=.false.,OUT_VIRTTEMP=.true.,LDIAG3D=.false.,mstrat=.false.,ctei_rm=0.50,MOM4ICE=.false.,cnvgwd=.false.,RUN_ENTHALPY=.false.,IOVR_SW=$IOVR_SW,zflxtvd=.true.,sashal=.true.,old_monin=.false.,newsas=.true.,ICTM=1,"
time /pan2/projects/gfsenkf/whitaker/gfs-dycore/bin/gfs < gfs.nml

setenv POSTGPSH /lfs1/projects/gfsenkf/hybda_realtime2012_scripts/global_postgpp.sh 
setenv POSTGPVARS "IDRT=0,IDRTC=4,IOC=$LONB,JOC=$LATB,MOO=255,MOOSLP=0"
setenv HOSTFILE $TMPDIR/machines
setenv OMP_NUM_THREADS 1
setenv nprocs $NSLOTS
setenv POSTGPLIST ${FIXGLOBAL}/global_kplist.1d.txt
setenv IO 360
setenv JO 181
set FH=0
while ($FH <= $FHMAX)
   set charfhr="f`printf %03i $FH`"
   setenv DATA postgp_fhr${charfhr}
   mkdir $DATA
   set PGBOUT=${datout}/pgrb.${charfhr}
   if ($FH == 0) then
     time sh $POSTGPSH ${datout}/sanl.dat ${datout}/flx.${charfhr} /dev/null $PGBOUT /dev/null $IO $JO 
   else
     time sh $POSTGPSH ${datout}/sig.${charfhr} ${datout}/flx.${charfhr} /dev/null $PGBOUT /dev/null $IO $JO 
   endif
   @ FH = $FH + $FHOUT
end

cd $homedir
setenv analdate `${incdate} $analdate 24`
echo "setenv analdate ${analdate}" >! analdate.csh
qsub testgfs.csh

exit
