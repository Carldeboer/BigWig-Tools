#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
#import itertools
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
import argparse
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
parser = argparse.ArgumentParser(description='This program takes a BW file as input and a file containing a mapping of chromsome names and outputs a new BW file with the new names.')
parser.add_argument('-i',dest='inBW',	metavar='<inBigWig>',help='Input bigwig file to translate', required=True);
parser.add_argument('-f',dest='functions',	metavar='<functions>',help='Functions to apply to data\n  One or more of {log(2|10)|[+-^*/]<num>|smooth<SD>|exp|pow<num>|default<num>}; separated by commas.', required=True);
parser.add_argument('-c',dest='chrsFile',	metavar='<chromSizesFile>',help='Input file containing the chromsome sizes: chr\\tsize\\n', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFilePath>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output messages', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

# this is for dubugging in interactive mode - comment out normally
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

allFunctions = args.functions.split(",");

def applyFunction(data,f):
	if re.match("^log$",f):
		return np.log(data);
	elif re.match("^log10$",f):
		return np.log10(data);
	elif re.match("^log2$",f):
		return np.log2(data);
	elif re.match("^exp$",f):
		return np.exp(data);
	else: 
		m = re.match("^([\^\*\+\-\/])([.-eE0-9]+)$",f);
		if m:
			if m.group(1)=="^":
				return data**float(m.group(2));
			elif m.group(1)=="*":
				return data*float(m.group(2));
			elif m.group(1)=="+":
				return data+float(m.group(2));
			elif m.group(1)=="/":
				return data/float(m.group(2));
			elif m.group(1)=="-":
				return data-float(m.group(2));
		m = re.match("^default([.-eE0-9]*)$",f);
		if m:
			invalid = np.isnan( data )
			data2 = data;
			data2[ invalid ] = float(m.group(1));
			return data2;
		m = re.match("^smooth([.-eE0-9]*)$",f);
		if m:
			return gaussian_filter(data, float(m.group(1)), mode='reflect', truncate=4.0)
		m = re.match("^pow([.-eE0-9]*)$",f);
		if m:
			return np.power(float(m.group(1)),data)
		raise Exception("Unknown function provided: %s"%f)
data = np.array(range(1,10));
for f in allFunctions:
	data = applyFunction(data,f);
	print(data);

print(allFunctions);
exit;

chromSizesFile = MYUTILS.smartGZOpen(args.chrsFile,'r');
chromSizes = {};
for line in chromSizesFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	chromSizes[data[0]]=int(data[1]);

curBW = BigWigFile(open(args.inBW))

toBW = subprocess.Popen(["wigToBigWig","stdin",args.chrsFile,args.outFP], stdin=subprocess.PIPE)
toBW.stdin.write("track type=wiggle_0\n")

for chr in chromSizes.keys():
	values = curBW.get_as_array( chr, 0, chromSizes[chr] )
	#print(chr);
	if values is not None:
		for f in allFunctions:
			values = applyFunction(values,f);
		toBW.stdin.write("fixedStep chrom=%s start=1 step=1\n"%(chr))
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


