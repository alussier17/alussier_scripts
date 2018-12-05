#Courtesy of Alexandre Lussier, PhD
#September 22, 2017

import subprocess as sp
import os, sys

program = 'python'
code = '/home/alussier/Takao/BSseeker2/bs_seeker2-align.py'
import time
import glob

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))
# function to select only non-hidden files within a folder

forwarddir = '/home/alussier/Takao/fastq_files/forward/'
forwardfiles = listdir_nohidden(forwarddir)
print forwardfiles

revdir = '/home/alussier/Takao/fastq_files/reverse/'
revfiles = listdir_nohidden(revdir)
print revfiles


aligndir = '/home/alussier/Takao/aligned/'
#setting directory to save filtered bam files into
try:
    os.mkdir(aligndir)
#if it didn't exist before, create said directory. if it does exist, return error
except Exception as error:
    print error
    print "Warning: ALIGNDIR already exists!"

adaptdir = '/home/alussier/Takao/adapters/'
#setting directory to find adapter files - each file = 2 adapter sequences (paired-end)
adaptfiles = listdir_nohidden(adaptdir)
print adaptfiles

genome = '/home/alussier/Takao/mm10.fa'

## trying to reproduce this line of code -> python bs_seeker2-align.py -1 RRBS1.1.fastq -2 RRBS1.2.fastq --aligner=bowtie -o RRBS.bam -g genome.fa -a adapter.txt

##  python /home/alussier/Takao/BSseeker2/bs_seeker2-align.py --aligner=bowtie -g /home/alussier/Takao/mm10.fa -a /home/alussier/Takao/adapters/ACAGTG.txt -1 /home/alussier/Takao/fastq_files/forward/CBKY8ANXX_8_ACAGTG_1.fastq -2 /home/alussier/Takao/fastq_files/reverse/CBKY8ANXX_8_ACAGTG_2.fastq -o /home/alussier/Takao/aligned/ACAGTG.bam --low=50 --up=450 -r

n=1
count =1


for x in range(0,15):
    print n
    print 'Selecting files for alignment'
    adaptfile = sorted(adaptfiles)[x]
    print adaptfile
    forwardfile = sorted(forwardfiles)[x]
    print forwardfile
    revfile = sorted(revfiles)[x]
    print revfile
        #sanity check for files
    
    name = adaptfile.replace('/Users/Alex/Documents/test/adapters/','')
    name = name.replace('.adapter.txt','')
    print("Index ID == %s " % name)
        #shows which sample is being used
        
    outfile = adaptfile.replace('.txt','.bam')
    outfile = outfile.replace('/home/alussier/Takao/adapters/', aligndir)
    print("Output file == %s " % outfile)
        #shows where the new file is going
    print("This is sample # %s" %count)

    if n == 1:
        t= time.localtime()
        print "Aligning", name, "to genome. Start = %s" % time.asctime(t)
        print "python", code,"--aligner=bowtie", "-g", genome, "-a", adaptfile, "-1", forwardfile, "-2", revfile, "-o", outfile, "--low=50 --up=450 --bt-p 4 --temp_dir=/home/alussier/Takao/temp_files"
        
        #command = "python", code,"--aligner=bowtie", "-g", genome, "-a", adaptfile, "-1", forwardfile, "-2", revfile, "-o", outfile, "--low=50 --up=450", "-r"
        #proc1 = sp.Popen(command, stdout=sp.PIPE, shell=isinstance(command, str))
        
        proc1 = sp.Popen( ["python", code,"--aligner=bowtie", "-g", genome, "-a", adaptfile, "-1", forwardfile, "-2", revfile, "-o", outfile, "--low=50", "--up=450", "-r", "--bt-p 2", "--temp_dir=/home/alussier/Takao/temp_files"])

    if n == 2:
        t= time.localtime()
        print "Aligning", name, "to genome. Start = %s" % time.asctime(t)
        proc2 = sp.Popen( ["python", code,"--aligner=bowtie", "-g", genome, "-a", adaptfile, "-1", forwardfile, "-2", revfile, "-o", outfile, "--low=50", "--up=450", "-r", "--bt-p 2", "--temp_dir=/home/alussier/Takao/temp_files"    ])
        
        proc1.communicate() # now wait until the process is done before continuing.
        proc2.communicate()
        n=0
    
    n = n+1
    count = count+1


