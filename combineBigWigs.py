#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
import argparse

parser = argparse.ArgumentParser(description='This program takes a BW file as input and applies one or more provided mathematical operations to the data, in the order provided.')
parser.add_argument('-i1',dest='inBW1',	metavar='<inBigWig1>',help='Input bigwig file1 to operate on', required=True);
parser.add_argument('-i2',dest='inBW2',	metavar='<inBigWig2>',help='Input bigwig file2 to operate on', required=True);
parser.add_argument('-f',dest='function',	metavar='<function>',help='Function with which to combine bigwigs', required=True);
parser.add_argument('-c',dest='chrsFile',	metavar='<chromSizesFile>',help='Input file containing the chromsome sizes: chr\\tsize\\n', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePrefix>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output messages', required=False);
parser.add_argument('-k',dest='chunks', metavar='<chunkSize>',help='The size of the chunks (in bp) to read/write when outputting [default=500000]', default=500000, required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

#import itertools
import os
import errno
import time
import warnings
import subprocess
import MYUTILS
import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter;
#from scipy import linalg
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from sklearn import mixture
import re
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
import threading

# this is for dubugging in interactive mode - comment out normally
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

args.chunks = int(args.chunks)

def applyFunction(data1,data2,f):
	if f=="^":
		return data1**data2;
	elif f=="*":
		return data1*data2;
	elif f=="+":
		return data1+data2;
	elif f=="/":
		return data1/data2;
	elif f=="-":
		return data1-data2;
	raise Exception("Unknown function provided: %s"%f)


chromSizesFile = MYUTILS.smartGZOpen(args.chrsFile,'r');
chromSizes = {};
for line in chromSizesFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	chromSizes[data[0]]=int(data[1]);

curBW1 = BigWigFile(open(args.inBW1))
curBW2 = BigWigFile(open(args.inBW2))


outStream = MYUTILS.smartGZOpen("%s.wig.gz"%(args.outFPre),"w");
outStream.write("track type=wiggle_0\n")

for chr in chromSizes.keys():
	last = 0;
	final = chromSizes[chr];
	sys.stderr.write("Outputting data for %s:\n"%(chr));
	while last!=final: # this breaks it up into chunks so that I'm not piping entire (human) chromosomes at once
		if args.verbose>0: sys.stderr.write("  Section %i - %i:\n"%(last,curLast));
		curLast = np.min([last+args.chunks,final]);
		values1 = curBW1.get_as_array( chr, last, curLast )
		values2 = curBW2.get_as_array( chr, last, curLast )
	#print(chr);
		if values1 is not None and values2 is not None: # what if only a chunk of a chromosome is missing? then I will get errors
			values = applyFunction(values1,values2,args.function);
			if last==0:
				outStream.write("fixedStep chrom=%s start=1 step=1\n"%(chr))
			outStream.write("\n".join(map(str,values)));
			outStream.write("\n");
			outStream.flush();
		last=curLast;


outStream.close();

toBW = subprocess.Popen(["wigToBigWig","%s.wig.gz"%(args.outFPre),args.chrsFile,"%s.bw"%(args.outFPre)])
temp = toBW.communicate()
if temp[0] is not None:
	sys.stderr.write("wigToBigWig: %s"%(temp[0]));
if temp[1] is not None:
	sys.stderr.write("wigToBigWig: %s"%(temp[1]));
if temp[0] is None and temp[1] is None and os.path.isfile("%s.bw"%(args.outFPre)): # if no errors, delete the original
	os.remove("%s.wig.gz"%(args.outFPre))
else:
	sys.stderr.write("Left wig.gz because of error.");
	
sys.stderr.write("Output all data.\n");

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
if (args.logFP is not None):
	logFile.close();


