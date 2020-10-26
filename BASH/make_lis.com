#!/bin/sh
#make_lis.com, USM, MS, 2020/10/19 
#Makes a list of HYCOM output binaries (*.a)
#Modified from MCB 2016

#---------------------- INPUT ----------------------------#
FIN='archs.0019'   # file name root
RNM='043_'         # experiment number (e.g. 19.0 = 190)

is=183 # days start
ie=214 # days end
js=01  # hour start
je=23  # hour end

# HYCOM output directories
SD='/p/work/xbxu/hycom/ATLc0.02/expt_04.3/data'            # grid
SD2='/p/work/xbxu/hycom/ATLc0.02/expt_04.3/data/tars_019g' # variables

#--------------------------------------------------------#

#last day  minus 1
ee=$(( ie - 1 ))
#echo $ee

#first day  plus 1
gg=$(( is + 1 ))
#echo $gg

# difference in days
ff=$((ie-is))
#echo $ff

# first do less than one day ---------------------------
#echo $FIN.lis

# save regional.depth.a regional.grid.a
echo "$SD/regional.depth.a" > "$RNM$FIN.lis"
echo "$SD/regional.grid.a" >> "$RNM$FIN.lis"

# temp fix to test it
echo "$SD/tidal.rh.a" >> "$RNM$FIN.lis"
echo "$SD/cb.a" >> "$RNM$FIN.lis"

hh=0
if [ $ff -eq 0 ]; then

for i in $(seq $is $is); do
    for j in $(seq   $js $je); do
        hh=$((hh+1))
        echo "$SD2/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done
fi

# if more days -------------------------------------------

if [ $ff -ge 1 ]; then

# first day
for i in $(seq $is $is); do
    for j in $(seq   $js 23); do
        hh=$((hh+1))
        echo "$SD2/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done


# then do whole days first
# -w does the padding :-)

for i in $(seq $gg  $ee); do
    for j in $(seq   0 23); do
        hh=$((hh+1))
        echo "$SD2/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done

# last day
for i in $(seq $ie $ie); do
    for j in $(seq   0 $je); do
        hh=$((hh+1))
        echo "$SD2/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done
fi

echo "number of time steps: ${hh}"
