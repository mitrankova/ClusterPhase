#!/usr/bin/sh

# runs with a production file path, run number, and segment number you want to look at
# for example, sh make_rawhit_filelist.sh /sphenix/lustre01/sphnxpro/production/run2pp/physics/ana464_nocdbtag_v001/ 53877 0

productionpath=/sphenix/lustre01/sphnxpro/production/run3auau/physics/ana514_nocdbtag_v001
run=72592
#66522

runnumber=$(printf "%08d" $run)
lrun=$(( run / 100 * 100 ))
hrun=$(( run / 100 * 100 + 100 ))
lowrun=$(printf "%08d" $lrun)
hirun=$(printf "%08d" $hrun)

rundir="$productionpath/DST_STREAMING_EVENT_ebdc00_0/run_${lowrun}_${hirun}"

# Find all segments for this run
for filepath in $(ls $rundir/*-$runnumber-*.root ); do
  filename=$(basename $filepath)
  segment=00000 #$(echo $filename | sed -n "s/.*-${runnumber}-\([0-9]\+\)\.root/\1/p")

  # Call the original script with this run and segment
  ls $productionpath/DST_STREAMING_EVENT_*/run_$lowrun\_$hirun/*-$runnumber-$segment.root > filelists/filelist_${runnumber}_${segment}.list

  echo "Created filelist for run $run segment $segment"
done


