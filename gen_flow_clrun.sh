#!/bin/bash

echo "Generation of clrun files for flow evolution"

Ns=32
Nt=14

tmax=1.0
eps=0.02

cores=4
mem=3.75
queue='mpi'

if [ $# -lt 4 ] 
then
  echo "./gen_flow_clrun.sh config start stop nrperrun"
  exit
fi

cnt=$2
runfilecnt=0
# for i in `seq $2 1 $3` ; do
while [ $cnt -le $3 ] ; do
  TIMERANDNR=${TIMESTAMP}_${RANDOM}
        TMPNAME=$TMPDIR/${TMPSTRING}_${TIMERANDNR}

  l=${1}${cnt}.bin
  runfile=clrun_mt_${l}

  echo "Generating $runfile ..."

  rm -f $runfile
  touch $runfile

  echo "#!/bin/bash" >> $runfile
  echo "#" >> $runfile
  echo "#$ -cwd" >> $runfile
  echo "#$ -S /bin/bash" >> $runfile
  echo "#$ -V" >> $runfile
  echo "#$ -q ${queue}" >> $runfile
  echo "#$ -pe threaded ${cores}" >> $runfile
  echo "#$ -l h_vmem=${mem}G" >> $runfile
  echo "#$ -N FLOW${l}" >> $runfile
  echo "export MV2_ENABLE_AFFINITY=0"  >> $runfile
  echo "export OMP_NUM_THREADS=\$NSLOTS"  >> $runfile
  echo "" >> $runfile
  echo "cd $PWD" >> $runfile
  for i in `seq 1 1 $4` ; do
    linside=${1}${cnt}.bin
    if [ ! -r $linside ]
    then
      echo "WARNING: $linside does not exist!"
    fi
  if [ $cnt -le $3 ] ; then
    echo "echo Calculation for ${linside} ..."  >> $runfile
    echo "/usr/bin/time ./flow -s $Ns -t $Nt -n 1 -c -w 5 -f $tmax -e $eps $linside -m &> ${linside}.status" >> $runfile
    echo "" >> $runfile
    echo "./genpollev_from_config.x $Ns $Nt ${linside} ${linside}.lpoll &>> ${linside}.status" >> $runfile
    echo "for ll in \`ls ${linside}_t*\` ; do" >> $runfile
    echo "./genpollev_from_config.x $Ns $Nt \$ll \${ll}.lpoll &>> ${linside}.status" >> $runfile
    echo "done" >> $runfile
  fi
  let cnt+=1
  done

  let runfilecnt+=1
done

num=$runfilecnt

echo "$num cluster runfiles generated!"
