import subprocess as sp
import os, sys

program = 'macs2'
import time
import glob

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))
# function to select only non-hidden files within a folder

### PX0182 ###

bamdir = '/home/alussier/newdata/PX0182/merged/'
#setting directory to find unfiltered bam files
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0182/merged/')
#telling python that every non-hidden file here is one of the files to use - bamfile
print bamfiles


macsdir = '/home/alussier/newdata/PX0182/macs/'
#setting directory to save filtered bam files into
try:
    os.mkdir(macsdir)
#if it didn't exist before, create said directory. if it does exist, return error
except Exception as error:
    print error
    print "Warning: MACSDIR already exists!"

## trying to reproduce this line of code -> macs2 callpeak -t bamfile -f BAMPE -g 2.9e9 -q 0.05 -B --call-summits --outdir macsdir --verbose=1
n = 1
format = 'BAMPE'
for bamfile in bamfiles:
#for loop - will take files in bamfiles directory one at a time and run following code
    print bamfile
#shows which bamfile it is processing
    macsfile = bamfile.split('/')[-1].replace('PX0182_C5DN1ANXX_6_','')
    macsfile = macsfile.replace('.filtered.merged.bam','')
#creating macs file with path ~/bam/file.ID.bam -> ~/macs/ID
    print bamfile
    print macsfile
#shows which bamfile was used and which filtered file was created
    print n
    t= time.localtime()
    if n == 1:
        print "Start Time: %s " % time.asctime(t)
        proc1 = sp.Popen( [program, 'callpeak', '-t', bamfile, '-n', macsfile, '-f', format,'-g 2.9e9', '-q 0.05', '-B','--call-summits', '--outdir', macsdir, '--verbose=1'] )
    if n == 2:
        t= time.localtime()
        print "Start Time: %s " % time.asctime(t)
        proc2 = sp.Popen( [program, 'callpeak', '-t', bamfile, '-n', macsfile, '-f', format,'-g 2.9e9', '-q 0.05', '-B','--call-summits', '--outdir', macsdir, '--verbose=1'] )
        proc1.communicate() # now wait until the process is done before continuing.
        proc2.communicate()
        n=0
    n = n+1

# BROAD ANALYSIS OF PEAKS #
print "BEGINNING BROAD ANALYSIS"
macs_broaddir = '/home/alussier/newdata/PX0182/macs_broad/'
try:
    os.mkdir(macs_broaddir)
except Exception as error:
    print error
    print "Warning: MACS_BROADDIR already exists!"

n = 1
format = 'BAMPE'
for bamfile in bamfiles:
    #for loop - will take files in bamfiles directory one at a time and run following code
    print bamfile
    #shows which bamfile it is processing
    macsfile = bamfile.split('/')[-1].replace('PX0182_C5DN1ANXX_6_','')
    macsfile = macsfile.replace('.filtered.merged.bam','')
    #creating macs file with path ~/bam/file.ID.bam -> ~/macs/ID
    print bamfile
    print macsfile
    #shows which bamfile was used and which filtered file was created
    print n
    t= time.localtime()
    if n == 1:
        print "Start Time: %s " % time.asctime(t)
        proc1 = sp.Popen( [program, 'callpeak', '-t', bamfile, '-n', macsfile, '-f', format,'-g 2.9e9', '-q 0.05', '-B','--broad', '--outdir', macs_broaddir, '--verbose=1'] )
    if n == 2:
        t= time.localtime()
        print "Start Time: %s " % time.asctime(t)
        proc2 = sp.Popen( [program, 'callpeak', '-t', bamfile, '-n', macsfile, '-f', format,'-g 2.9e9', '-q 0.05', '-B','--broad', '--outdir', macs_broaddir, '--verbose=1'] )
        proc1.communicate() # now wait until the process is done before continuing.
        proc2.communicate()
        n=0
    n = n+1
