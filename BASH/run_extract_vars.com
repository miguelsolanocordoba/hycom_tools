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
EXPT='04'
num='3'
RNNM=$EXPT$num
LISN=$RNNM'_archs.0019.lis'
echo $LISN

#range i:1-52, j:1-38
iblks=1
iblke=52
jblks=1
jblke=38
stp=19

# path where the simulation is done
PTHNM='/p/work/tfohycom/hycom/ATLc0.02/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/ATLc0.02/expt_'$EXPT'.'$num'/'
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
  echo "#PBS -A ONRDC########"            >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -l walltime=05:00:00"        >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -q standard"                    >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "#PBS -l select=1:ncpus=44:mpiprocs=1"  >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "# Change to the specified directory" >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "cd /p/work/tfohycom/hycom/ATLc0.02/expt_${EXPT}.${num}" >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsXTS_j${j}_${i}a
  echo "aprun -b -n 1 ./extract_vars2D.x  < ./extract_${RNNM}_${j}_${i}a.in > ./extract_${RNNM}_${j}_${i}a.out" >> ${PTHHO}pbsXTS_j${j}_${i}a

  qsub ${PTHHO}pbsXTS_j${j}_${i}a

done

