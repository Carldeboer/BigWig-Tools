#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
import argparse
parser = argparse.ArgumentParser(description='This program takes a BW file as input and a file containing a mapping of chromsome names and outputs a new BW file with the new names.')
parser.add_argument('-i',dest='inBW',	metavar='<inBigWig>',help='Input bigwig file to translate', required=True);
parser.add_argument('-m',dest='inMap',	metavar='<inMap>',help='Input file containing a mapping: chrNew\\tchrOld[\\tchrOld...]\\n', required=True);
parser.add_argument('-c',dest='chrsFile',	metavar='<chromSizesFile>',help='Input file containing the chromsome sizes: chrNew\\tsize\\n', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFilePath>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output messages', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

#import itertools
import warnings
import subprocess
import MYUTILS
import numpy as np
import scipy as sp
#from scipy.ndimage.filters import gaussian_filter;
#from scipy import linalg
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from sklearn import mixture
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile


# this is for dubugging in interactive mode - comment out normally
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

oldToNew = {};
transChrs = [];
inFile=MYUTILS.smartGZOpen(args.inMap,'r');
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	for i in range(0,len(data)):
		oldToNew[data[i]] = data[0];
		transChrs.append(data[i]);

inFile.close();

chromSizesFile = MYUTILS.smartGZOpen(args.chrsFile,'r');
chromSizes = {};
for line in chromSizesFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	chromSizes[data[0]]=int(data[1]);

curBW = BigWigFile(open(args.inBW))

toBW = subprocess.Popen(["wigToBigWig","stdin",args.chrsFile,args.outFP], stdin=subprocess.PIPE)
toBW.stdin.write("track type=wiggle_0\n")

for chr in transChrs:
	values = curBW.get_as_array( chr, 0, chromSizes[oldToNew[chr]] )
	#print(chr);
	if values is not None:
		sys.stderr.write("Adding %s -> %s\n"%(chr, oldToNew[chr]));
		toBW.stdin.write("fixedStep chrom=%s start=1 step=1\n"%(oldToNew[chr]))
		toBW.stdin.write("\n".join(map(str,values)));
		toBW.stdin.write("\n");

temp = toBW.communicate();

if temp[0] is not None:
	sys.stderr.write(temp[0]);

if temp[1] is not None:
	sys.stderr.write(temp[1]);

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
if (args.logFP is not None):
	logFile.close();


