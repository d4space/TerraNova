#!/bin/bash -e

# Sets up and starts a run using the Condor system
# Usage: condor_run.sh w/z <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir> [ <which_sect> ]
#  Note: <which_sect> is optional and for advanced user only
#        <which_sect> specifies the one sector that user want to submit to condor system
#### Warning!  Single-sector Condor run broken and deprecated

condorusage(){
echo "Usage: `basename $0` w/z <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir>"
exit 1
}
#[ $# -lt 6 ] && condorusage \
#             || printf "Running Directory: $2\nInput Setting: $3\nHistogram File: $4\nOutput File: $5\nPDF Directory $6\n"
[ $# -lt 6 ] && condorusage

# check first argument, to make sure supported, and set executable
if [ $1 = "z" ] || [ $1 = "w" ]; then
   EXEC=fewz$1
else
   echo "Unrecognized argument; defaulting to neutral current."
   EXEC=fewzz
fi

RUNDIR=${2%/}
python scripts/create_parallel.py $1 $3 $RUNDIR
if ! [ -e $RUNDIR/$EXEC ] ; then
   cp $EXEC $RUNDIR
fi

### always copy input file, in case changed
cp $3 $RUNDIR
if ! [ -e $RUNDIR/$4 ] ; then
    python scripts/get_bin_files.py $4 $RUNDIR
fi
python scripts/create_condor_jobs.py $1 ${RUNDIR##*/} $3 ../$4 $5 ../../${6%/}
cd $2
for isec in `seq 1 $7`; do 
    echo ${PWD}
    cat job_desc | sed -e 's/..\/..\/\//\//g' -e 's/Y/'$(($isec-1))'/g' -e 's/X/'$2'/g' -e 's/sectm1/'$2$(($isec-1))'/g' > job_desc_$isec
    echo $isec
    cd $2$(($isec-1))
    mv ../job_desc_$isec ./
    chmod +x job_desc_$isec
    mv job_desc_$isec job_desc_$isec.sh
    qsub -q cms -l walltime=500:00:00,cput=500:00:00 job_desc_$isec.sh -N WpToEleNu_$isec
    cd ..
done;

echo "Run the following to post-process output files: finish.sh $2 <order>.$5"
