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
EXPT='19'
num='0'
RNNM=$EXPT$num
maxobs=720
echo $RNNM

#range i:1-60, j:1-35
jblks=1
jblke=35
iblks=1
iblke=60
stp=5

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
echo $PTHNM
echo $PTHHO

#filter variables
ftnm='modes16'
#ftnm='surface'

for j in $(seq $jblks $stp $jblke ); do
  echo $j

  i=$iblks
  je=$((j+stp-1))

# energy balance
  echo "along x: do ${i} to ${iblke}"
  echo "along y: do ${j} to ${je}"
  echo "$j        'jblks'  "   >  ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "$je       'jblke'  "   >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  "   >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  "   >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "$maxobs   'maxobs'  "  >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in
  echo "$ftnm     'ftnm'  "    >> ${PTHNM}ebalI_${RNNM}_${j}_${i}.in

# make the pbs file EBAL
  echo "#! /bin/sh"                         > ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -N EBIj${j}_${i}"             >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -o EBIj${j}_${i}.out"         >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -e EBIj${j}_${i}.err"         >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -M miguel.solano@usm.edu"  >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -m bea"                       >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -A ONRDC45592567"             >> ${PTHHO}pbsEBI_j${j}_${i}
#  echo "#PBS -A NRLSS03755018"             >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -S /bin/bash"                 >> ${PTHHO}pbsEBI_j${j}_${i}
#  echo "#PBS -S /bin/csh"                 >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -l walltime=05:00:00"         >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -q standard"                  >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"          >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "# Change to the specified directory"            >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "# Execute the serial executable on 1 core"      >> ${PTHHO}pbsEBI_j${j}_${i}
  echo "./rotary_spectra_uvModes.x    < ./ebalI_${RNNM}_${j}_${i}.in   > ./ebalI_${RNNM}_${j}_${i}.out"   >> ${PTHHO}pbsEBI_j${j}_${i}

  qsub ${PTHHO}pbsEBI_j${j}_${i}
  
done

