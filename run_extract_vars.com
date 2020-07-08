#! /bin/sh
#reads *.a files from *.list file and extracts the data in the *.a files
#and appends it to a binary output file

#|: control operator to separate a sequence commands in a pipeline
#egrep: match pattern in file
#tail -7:  list last 7 lines
#wc -l: count lines
#>> extract_all.in: opens file extract_all.in and appends text
#> extract_all.in: opens file extract_all.in and (over)writes text
#awk prints text read in file test
#>& redirects both the stdout and the stderr of command to filename
#time: time a simple command or give resource usage

# run number
EXPT='21'
num='6'
RNNM=$EXPT$num
LISN=$RNNM'_archv.2019.lis'
echo $LISN

#range i:1-60, j:1-35
#jblks=1    # Surface stress/heat (i=1-60; j=1:35)
#jblke=35
#iblks=1
#iblke=60
jblks=15
jblke=19
iblks=23
iblke=25
stp=5

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
echo $PTHNM
echo $PTHHO

for j in $(seq $jblks $stp $jblke ); do
  echo $j

  i=$iblks
  je=$((j+stp-1))

  echo "along x: do ${i} to ${iblke}"
  echo "along y: do ${j} to ${je}"
  egrep "\.a" $LISN | tail -1000 | wc -l > test
  echo "$j    'jblks'  " > ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  echo "$je    'jblke'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  echo "$iblks    'iblks'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  echo "$iblke    'iblke'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  awk '{print $0"     \x27numfls\x27"}' test >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in
  egrep "\.a" $LISN | tail -1000 >> ${PTHNM}extract_${RNNM}_${j}_${i}a.in

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -N XTSj${j}_${i}"             >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -o XTSj${j}_${i}.out"         >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -e XTSj${j}_${i}.err"         >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -M miguel.solano@usm.edu" >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -m bea"                      >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -A ONRDC45592567"            >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -l walltime=03:00:00"        >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -q standard"                    >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "# Change to the specified directory" >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "aprun -b -n 1 ./extract_vars_TS_c004_v12.x  < ./extract_${RNNM}_${j}_${i}a.in > ./extract_${RNNM}_${j}_${i}a.out" >> ${PTHHO}pbsXTS_j${j}_${i}a

  qsub ${PTHHO}pbsXTS_j${j}_${i}a

done

