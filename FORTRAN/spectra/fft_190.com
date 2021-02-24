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
  echo "$j        'jblks'  "   >  ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "$je       'jblke'  "   >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  "   >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  "   >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "$maxobs   'maxobs'  "  >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in
  echo "$ftnm     'ftnm'  "    >> ${PTHNM}fftuvp_${RNNM}_${j}_${i}.in

# make the pbs file EBIfft
  echo "#! /bin/sh"                         > ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -N EBIfftj${j}_${i}"             >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -o EBIfftj${j}_${i}.out"         >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -e EBIfftj${j}_${i}.err"         >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -M miguel.solano@usm.edu"  >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -m bea"                       >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -A ONRDC45592567"             >> ${PTHHO}pbsEBIfft_j${j}_${i}
#  echo "#PBS -A NRLSS03755018"             >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -S /bin/bash"                 >> ${PTHHO}pbsEBIfft_j${j}_${i}
#  echo "#PBS -S /bin/csh"                 >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -l walltime=05:00:00"         >> ${PTHHO}pbsEBIfft_j${j}_${i}  # HH:MM:SS
  echo "#PBS -q standard"                  >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"          >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "# Change to the specified directory"            >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "# Execute the serial executable on 1 core"      >> ${PTHHO}pbsEBIfft_j${j}_${i}
  echo "./fft_uvpModes.x    < ./fftuvp_${RNNM}_${j}_${i}.in   > ./fftuvp_${RNNM}_${j}_${i}.out"   >> ${PTHHO}pbsEBIfft_j${j}_${i}

  qsub ${PTHHO}pbsEBIfft_j${j}_${i}
  
done

