"""
Produce condor submission file "job_desc"
Run "condor_submit job_desc" to submit
Used by "condor_run.sh"
Usage: python $0 z/w <outputdir> <inputfile> <histofile> <outputfile_extension> <pdfdir>
Note: <outputdir> should be a directory name without a trailing '/' and a preceding path
      <pdfdir> should be a directory name or path without trailing '/'
"""

import sys
from defs import getsects

jobname = 'job_desc'
condor_out = 'condor_output.out'
condor_err = 'condor_error.err'
condor_log = 'condor_log.log'

try:
    if len(sys.argv) < 7:
        raise Exception

    else:
        boson = (sys.argv[1]).upper()
        sectors = getsects(sys.argv[3], boson)
        outputdir = sys.argv[2]
        inputfile = sys.argv[3]
        histofile = sys.argv[4]
        outputfile = sys.argv[5]
        pdfdir = sys.argv[6]

except Exception:
    print('Missing arguments.')
    raise

try:
    job_file = open(outputdir + '/' + jobname, 'w')
except IOError:
    print('Error creating job file.')
    raise

#Below, the single commented lines need to be changed into the PBS submission file
#
try:
    job_file.write('#!/bin/csh -f\n')
    job_file.write('#PBS -N fewzjob_Y\n')
    job_file.write('#PBS -o output.log\n')
    job_file.write('#PBS -e output.err\n')
    job_file.write('#PBS -m ae\n')
    job_file.write('#PBS -r n\n')
    job_file.write('#PBS -q cms\n')
    job_file.write('#PBS -W group_list=cms\n')
    #job_file.write('#PBS -l walltime=504:00:00\n')
    job_file.write('#PBS -l nodes=1:ppn=1\n')
    job_file.write('#PBS -l pmem=1000mb\n')
    job_file.write('\n')
    #job_file.write('set nonomatch\n')
    #job_file.write('echo ""\n')
    #job_file.write('echo "Job is running on \`uname -a\`"\n')
    #job_file.write('if ( \${OSTYPE} == "linux" ) then\n')
    #job_file.write('  set processor = \`sort /proc/cpuinfo | uniq | gawk -F: '(substr(\$1,1,10)=="model name"){print \$2}'\`\n')
    #job_file.write('  set rate = \`sort /proc/cpuinfo | uniq | gawk -F: '(substr(\$1,1,7)=="cpu MHz"){print substr(\$2,1,6)}'\`\n')
    #job_file.write('  echo "Processor info : " \$processor \$rate "MHz"\n')
    #job_file.write('endif\n')
    #job_file.write('set start = \`date\`\n')
    job_file.write('\n')
    #job_file.write('setenv WORKDIR `mktemp -d /tmp/asvyatko_XXXXXXXX`\n')
    job_file.write('setenv PDFDIR '+sys.argv[6]+'\n')
    #job_file.write('setenv FEWZDIR /scratch/lustreA/a/asvyatko/PBS_TESTS/FEWZ_BASE_2.1/bin\n')
    job_file.write('setenv WORKDIR ${PDFDIR}/bin/X/sectm1\n')
    job_file.write('setenv OUTDIR ${PDFDIR}/bin/X/sectm1\n')
    job_file.write('cd ${WORKDIR}\n')
    if (boson == 'W'):
        job_file.write('cp ../fewzw .\n')
        job_file.write('cp ../input_w.txt .\n')
        job_file.write('cp ../histograms.txt .\n')
        #job_file.write('./fewzz -i input_z.txt -h histograms.txt -o X.dat -p ${PDFDIR} -s Y\n')
        job_file.write('./fewzw -i input_w.txt -h histograms.txt -o X.dat -p ../../.. -s Y\n')
        #job_file.write('uptime\n')
        job_file.write('\n')
    else:
        if (boson != 'Z'):
            print('Warning: unrecognized parameter; defaulting to neutral current\n')
        job_file.write('cp ../fewzz .\n')
        job_file.write('cp ../input_z.txt .\n')
        job_file.write('cp ../histograms.txt .\n')
        job_file.write('\n')
        job_file.write('./fewzz -i input_z.txt -h histograms.txt -o X.dat -p ${PDFDIR} -s Y\n')
        #job_file.write('./fewzz -i input_z.txt -h histograms.txt -o X.dat -p ../../.. -s Y\n')
        #job_file.write('uptime\n')
        job_file.write('\n')
        #job_file.write('set rtime = \`tail -1 x | cut -f 1 -d " "\`\n')
        #job_file.write('set utime = \`tail -1 x | cut -f 2 -d " "\`\n')
        #job_file.write('set stime = \`tail -1 x | cut -f 3 -d " "\`\n')
        #job_file.write('set stat  = \`tail -1 x | cut -f 4 -d " "\`\n')
        #job_file.write('if ( ! -e  \${OUTDIR}/\${basename}.result ) cp outfile \${OUTDIR}/\${basename}.result\n')
        #job_file.write('if ( -e  \${OUTDIR}/\${basename}.result ) cat outfile >> \${OUTDIR}/\${basename}.result\n')
        job_file.write('echo "Current directory:"\n')
        job_file.write('pwd\n')
        job_file.write('ls -lrtAFh\n')
        #job_file.write('set end = \`date\`\n')
        #job_file.write('echo ""\n')
        #job_file.write('echo "Job end \`date\`"\n')
        job_file.write('echo ""\n')
        #job_file.write('echo \$PBS_JOBNAME \$stat \$start \$end \`uname -n | cut -f 1 -d .\` \$processor \$rate \$rtime \$utime \$stime  >> \${OUTDIR}/SUMMARY\n')
        #job_file.write('rm -r \${WORKDIR}\n')
        job_file.write('exit ${status}\n')
        job_file.close()

except IOError:
    print('Error writing to job file.')
    raise
