import subprocess as sp
import os, sys

program = 'samtools'
import time
import glob

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))
# function to select only non-hidden files within a folder


### PX0182 - C5DN1ANXX_6 ### - DONE -

bamdir = '/home/alussier/newdata/PX0182/C5DN1ANXX_6/125nt/Rnor_6.0/bwa-0.5.7/'
#setting directory to find unfiltered bam files
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0182/C5DN1ANXX_6/125nt/Rnor_6.0/bwa-0.5.7/')
#telling python that every non-hidden file here is one of the files to use - bamfile
print bamfiles


filtereddir = '/home/alussier/newdata/PX0182/C5DN1ANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
#setting directory to save filtered bam files into
try:
    os.mkdir(filtereddir)
#if it didn't exist before, create said directory. if it does exist, return error
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam

for bamfile in bamfiles:
#for loop - will take files in bamfiles directory one at a time and run following code
    print bamfile
#shows which bamfile it is processing
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
#creating filtered file with path ~/bam/file.bam -> ~/filtered/file.filtered.bam
    print bamfile
    print filteredfile
#shows which bamfile was used and which filtered file was created

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
#reproducing this line of code in terminal -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam
# -b = bam file, -q 10 = quality score >10, -f 2 = keep only matching pair (paired-end reads), -F 1536 = remove duplicates (1024) + remove QC failures (512)
    proc.communicate() # now wait until the process is done before continuing.
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)


## NEXT ONES ARE THE EXACT SAME, BUT WITH DIFFERENT INPUT FOLDERS. 

### PX0182 - C5DWPANXX_5 ###

bamdir = '/home/alussier/newdata/PX0182/C5DWPANXX_5/125nt/Rnor_6.0/bwa-0.5.7/'
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0182/C5DWPANXX_5/125nt/Rnor_6.0/bwa-0.5.7/')
print bamfiles

filtereddir = '/home/alussier/newdata/PX0182/C5DWPANXX_5/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
try:
    os.mkdir(filtereddir)
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam
for bamfile in bamfiles:
    print bamfile
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
    print bamfile
    print filteredfile

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
    proc.communicate()
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)



### PX0183 - C5DN1ANXX_7 ### - DONE -

bamdir = '/home/alussier/newdata/PX0183/C5DN1ANXX_7/125nt/Rnor_6.0/bwa-0.5.7/'
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0183/C5DN1ANXX_7/125nt/Rnor_6.0/bwa-0.5.7/')
print bamfiles


filtereddir = '/home/alussier/newdata/PX0183/C5DN1ANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
try:
    os.mkdir(filtereddir)
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam

for bamfile in bamfiles:
    print bamfile
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
    print bamfile
    print filteredfile

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
    proc.communicate()
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)


### PX0183 - C5DWPANXX_6 ###

bamdir = '/home/alussier/newdata/PX0183/C5DWPANXX_6/125nt/Rnor_6.0/bwa-0.5.7/'
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0183/C5DWPANXX_6/125nt/Rnor_6.0/bwa-0.5.7/')
print bamfiles


filtereddir = '/home/alussier/newdata/PX0183/C5DWPANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
try:
    os.mkdir(filtereddir)
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam

for bamfile in bamfiles:
    print bamfile
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
    print bamfile
    print filteredfile

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
    proc.communicate()
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)


### PX0184 - C5DN1ANXX_8 ###

bamdir = '/home/alussier/newdata/PX0184/C5DN1ANXX_8/125nt/Rnor_6.0/bwa-0.5.7/'
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0184/C5DN1ANXX_8/125nt/Rnor_6.0/bwa-0.5.7/')
print bamfiles


filtereddir = '/home/alussier/newdata/PX0184/C5DN1ANXX_8/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
try:
    os.mkdir(filtereddir)
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam

for bamfile in bamfiles:
    print bamfile
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
    print bamfile
    print filteredfile

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
    proc.communicate()
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)



### PX0184 - C5DWPANXX_7 ###

bamdir = '/home/alussier/newdata/PX0184/C5DWPANXX_7/125nt/Rnor_6.0/bwa-0.5.7/'
bamfiles = listdir_nohidden('/home/alussier/newdata/PX0184/C5DWPANXX_7/125nt/Rnor_6.0/bwa-0.5.7/')
print bamfiles

filtereddir = '/home/alussier/newdata/PX0184/C5DWPANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
try:
    os.mkdir(filtereddir)
except Exception as error:
    print error
    print "Warning: FILTEREDDIR already exists!"

## trying to reproduce this line of code -> samtools view -b -q 10 -f 2 -F 1536 -o out.bam in.bam

for bamfile in bamfiles:
    print bamfile
    filteredfile = filtereddir + bamfile.split('/')[-1].replace('.bam','.filtered.bam')
    print bamfile
    print filteredfile

    t= time.localtime()
    print "Start Time: %s " % time.asctime(t)
    proc = sp.Popen( [program, 'view', '-b', '-q 10', '-f 2', '-F 1536', '-o', filteredfile,bamfile] )
    proc.communicate()
    t= time.localtime()
    print "End Time: %s " % time.asctime(t)

