#$ -S /bin/bash
#$ -cwd
#$ -N filter
#$ -o filter.out.$JOB_ID
#$ -j y
#$ -l h_rt=24:00:00
#$ -l h_vmem=12G

cd location_of_data

for i in CCFR1ANXX*;
do echo $i; 
samtools view -b -q 10 -f 2 -F 1536 -o location_of_data/dev2_seq/filtered/$i $i;
done
