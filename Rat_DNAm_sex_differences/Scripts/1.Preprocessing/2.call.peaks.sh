#$ -S /bin/bash
#$ -cwd
#$ -N peaks
#$ -o peaks.out.$JOB_ID
#$ -j y
#$ -l h_rt=24:00:00
#$ -l h_vmem=12G

export PYTHONPATH=/programs/macs2-2.1.4/lib64/python2.7/site-packages
export PATH=/programs/macs2-2.1.4/bin:$PATH

cd location_of_data

for i in CCFR1ANXX_*;
do echo $i;
macs2 callpeak -t $i -n $i -f BAMPE -g 2.9e9 -q 0.05 -B --call-summits --verbose 1 --outdir location_of_data/dev2_seq/peaks/$i ;
done


