import subprocess as sp
import os, sys
import glob
import time

program = 'samtools'

def listdir_nohidden(path):
    return glob.glob(os.path.join(path,'*'))


### PX0182 ###

filtereddir1 = '/home/alussier/newdata/PX0182/C5DN1ANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles1 = listdir_nohidden('/home/alussier/newdata/PX0182/C5DN1ANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles1

filtereddir2 = '/home/alussier/newdata/PX0182/C5DWPANXX_5/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles2 = listdir_nohidden('/home/alussier/newdata/PX0182/C5DWPANXX_5/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles2

mergeddir = '/home/alussier/newdata/PX0182/merged/'
try:
    os.mkdir(mergeddir)
except Exception as error:
    print error
    print "Warning: MERGEDDIR already exists!"

## Trying to reproduce the following samtools code:
#   samtools merge -r out.bam in1.bam in2.bam

for filteredfile1 in filteredfiles1:
    print filteredfile1
    mergedfile = mergeddir + filteredfile1.split('/')[-1].replace('.bam','.merged.bam')
    print filteredfile1
    print mergedfile
    filt1name = filteredfile1.split('/')[-1].replace('PX0182_C5DN1ANXX_6_','')
    print filt1name
    for filteredfile2 in filteredfiles2:
        if filt1name in filteredfile2:
            print filteredfile2
            t= time.localtime()
            print "Start Time: %s " % time.asctime(t)
            print (program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2)
            proc=sp.Popen ( [program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2 ])
            proc.communicate()
    t= time.localtime()
    print "End time: %s " % time.asctime(t)

print "HERE MARKS THE END OF PX0182"

### PX0183 ###

filtereddir1 = '/home/alussier/newdata/PX0183/C5DN1ANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles1 = listdir_nohidden('/home/alussier/newdata/PX0183/C5DN1ANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles1

filtereddir2 = '/home/alussier/newdata/PX0183/C5DWPANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles2 = listdir_nohidden('/home/alussier/newdata/PX0183/C5DWPANXX_6/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles2


mergeddir = '/home/alussier/newdata/PX0183/merged/'
try:
    os.mkdir(mergeddir)
except Exception as error:
    print error
    print "Warning: MERGEDDIR already exists!"

## Trying to reproduce the following samtools code:
#   samtools merge -r out.bam in1.bam in2.bam

for filteredfile1 in filteredfiles1:
    print filteredfile1
    mergedfile = mergeddir + filteredfile1.split('/')[-1].replace('.bam','.merged.bam')
    print filteredfile1
    print mergedfile
    filt1name = filteredfile1.split('/')[-1].replace('PX0183_C5DN1ANXX_7_','')
    print filt1name
    for filteredfile2 in filteredfiles2:
        if filt1name in filteredfile2:
            print filteredfile2
            t= time.localtime()
            print "Start Time: %s " % time.asctime(t)
            print (program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2)
            proc=sp.Popen ( [program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2 ])
            proc.communicate()
    t= time.localtime()
    print "End time: %s " % time.asctime(t)

print "HERE MARKS THE END OF PX0183"

### PX0184 ###

filtereddir1 = '/home/alussier/newdata/PX0184/C5DN1ANXX_8/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles1 = listdir_nohidden('/home/alussier/newdata/PX0184/C5DN1ANXX_8/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles1

filtereddir2 = '/home/alussier/newdata/PX0184/C5DWPANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/'
filteredfiles2 = listdir_nohidden('/home/alussier/newdata/PX0184/C5DWPANXX_7/125nt/Rnor_6.0/bwa-0.5.7/filtered/')
print filteredfiles2

mergeddir = '/home/alussier/newdata/PX0184/merged/'
try:
    os.mkdir(mergeddir)
except Exception as error:
    print error
    print "Warning: MERGEDDIR already exists!"

## Trying to reproduce the following samtools code:
#   samtools merge -r out.bam in1.bam in2.bam

for filteredfile1 in filteredfiles1:
    print filteredfile1
    mergedfile = mergeddir + filteredfile1.split('/')[-1].replace('.bam','.merged.bam')
    print filteredfile1
    print mergedfile
    filt1name = filteredfile1.split('/')[-1].replace('PX0184_C5DN1ANXX_8_','')
    print filt1name
    for filteredfile2 in filteredfiles2:
        if filt1name in filteredfile2:
            print filteredfile2
            t= time.localtime()
            print "Start Time: %s " % time.asctime(t)
            print (program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2)
            proc=sp.Popen ( [program, 'merge', '-r', mergedfile, filteredfile1, filteredfile2 ])
            proc.communicate()
    t= time.localtime()
    print "End time: %s " % time.asctime(t)

print "HERE MARKS THE END OF PX0184"
