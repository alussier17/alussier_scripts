#Courtesy of Alexandre Lussier, PhD
#September 22, 2017

import subprocess as sp
import os, sys

program = '/home/alussier/Takao/BSseeker2/bs_seeker2-call_methylation.py'
import time
import glob

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))
# function to select only non-hidden files within a folder

calldir = '/home/alussier/Takao/meth_call/'
#setting directory to save methylation call files into
try:
    os.mkdir(calldir)
#if it didn't exist before, create said directory. if it does exist, return error
except Exception as error:
    print error
    print "Warning:", calldir, "already exists!"

bamdir = '/home/alussier/Takao/aligned/*.bam'
#setting directory to find bam files
bamfiles = glob.glob(bamdir)
print bamfiles

genome = '--db=/home/alussier/Takao/BSseeker2/bs_utils/reference_genomes/mm10.fa_rrbs_50_450_bowtie/'

## trying to reproduce this line of code -> python bs_seeker2-call_methylation.py -x -i bamfile --db genome -o calldir

n=1

for bamfile in bamfiles:
    print bamfile
    #sanity check for files
    
    name = bamfile.replace('/home/Alex/Documents/test/aligned/','')
    name = name.replace('.bam','')
    print("Index ID == %s " % name)
    #shows which sample is being used
    
    
    if n == 1:
        t= time.localtime()
        print "Calling methylation scores for", name, "Start = %s" % time.asctime(t)
        proc = sp.Popen( ["python", program, "-i", bamfile,  "-o", calldir, "-x", genome])
    
    if n == 2:
        t= time.localtime()
        print "Calling methylation scores for", name, "Start = %s" % time.asctime(t)
        proc2 = sp.Popen( ["python", program, "-i", bamfile, "-o", calldir, "-x", genome] )
    
        proc.communicate()
        proc2.communicate()
        n=0
    
    n = n+1
