#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
import argparse

parser = argparse.ArgumentParser(description='This program takes a BW file as input and applies one or more provided mathematical operations to the data, in the order provided.')
parser.add_argument('-i',dest='inBW',	metavar='<inBigWig>',help='Input bigwig file to translate', required=True);
parser.add_argument('-gt',dest='greaterThan',	metavar='<greaterThan>',help='A float that when data points are above, BED entries will be made.', required=False);
parser.add_argument('-lt',dest='lessThan',	metavar='<lessThan>',help='A float that when data points are below, BED entries will be made.', required=False);
parser.add_argument('-c',dest='chrsFile',	metavar='<chromSizesFile>',help='Input file containing the chromsome sizes: chr\\tsize\\n', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output messages', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
parser.add_argument('-k',dest='chunks', metavar='<chunkSize>',help='The size of the chunks (in bp) to read/write when outputting [default=500000]', default=500000, required=False);
args = parser.parse_args();

#import itertools
import os
import warnings
import subprocess
import MYUTILS
import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter;
from scipy.ndimage.filters import uniform_filter1d;
#from scipy import linalg
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from sklearn import mixture
import re
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

# this is for dubugging in interactive mode - comment out normally
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

args.chunks = int(args.chunks);
if args.lessThan is not None:
	args.lessThan = float(args.lessThan);
if args.greaterThan is not None:
	args.greaterThan = float(args.greaterThan);
if args.lessThan is None and args.greaterThan is None:
	raise Exception("must specify one of -gt or -lt!");





chromSizesFile = MYUTILS.smartGZOpen(args.chrsFile,'r');
chromSizes = {};
for line in chromSizesFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	chromSizes[data[0]]=int(data[1]);

curBW = BigWigFile(open(args.inBW))

outFile = MYUTILS.smartGZOpen(args.outFP,"w");

gtI=1;
ltI=1;

for chr in chromSizes.keys():
	last = 0;
	final = chromSizes[chr];
	sys.stderr.write("Outputting data for %s:\n"%(chr));
	gtStart = None;
	ltStart = None;
	while last!=final: # this breaks it up into chunks so that I'm not piping entire (human) chromosomes at once
		curLast = np.min([last+args.chunks,final]);
		if args.verbose>0: sys.stderr.write("  Section %i - %i:\n"%(last,curLast));
		curEnd = np.min([curLast, final]);
		curSt = np.max([last,0]);
		values = curBW.get_as_array( chr, curSt, curEnd )
		if values is not None:
			if args.lessThan is not None:
				i=curSt;
				for v in values:
					if v>=args.lessThan:
						if ltStart is not None:
							if args.verbose>1: sys.stderr.write("  	Found LT end %i - %i:\n"%(ltStart,i-1));
							outFile.write("%s\t%i\t%i\tLT%i\n"%(chr, ltStart, i-1, ltI))
							ltStart=None;
							ltI+=1;
					elif ltStart is None: # v is less than
						if args.verbose>1: sys.stderr.write("  	Found LT start %i:\n"%(i));
						ltStart = i;
					i+=1;
			if args.greaterThan is not None:
				i=curSt;
				for v in values:
					if v<=args.greaterThan:
						if gtStart is not None:
							if args.verbose>1: sys.stderr.write("  	Found GT end %i - %i:\n"%(gtStart,i-1));
							outFile.write("%s\t%i\t%i\tGT%i\n"%(chr, gtStart, i-1, gtI))
							gtStart=None;
							gtI+=1;
					elif gtStart is None: # v is greater than
						if args.verbose>1: sys.stderr.write("  	Found GT start %i:\n"%(i));
						gtStart = i;
					i+=1;
		last = curLast;
	if ltStart is not None:
		if args.verbose>1: sys.stderr.write("  	Found final LT end %i - %i:\n"%(ltStart,i-1));
		outFile.write("%s\t%i\t%i\tLT%i\n"%(chr, ltStart, i-1, ltI))
		ltStart = None;
		ltI+=1;
	if gtStart is not None:
		if args.verbose>1: sys.stderr.write("  	Found final GT end %i - %i:\n"%(gtStart,i-1));
		outFile.write("%s\t%i\t%i\tGT%i\n"%(chr, gtStart, i-1, gtI))
		gtStart = None;
		gtI+=1;
outFile.close();

if (args.logFP is not None):
	logFile.close();


