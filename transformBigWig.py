#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
import argparse

parser = argparse.ArgumentParser(description='This program takes a BW file as input and applies one or more provided mathematical operations to the data, in the order provided.')
parser.add_argument('-i',dest='inBW',	metavar='<inBigWig>',help='Input bigwig file to translate', required=True);
parser.add_argument('-f',dest='functions',	metavar='<functions>',help='Functions to apply to data\n  One or more of {log(2|10)|[+-^*/]<num>|(>=|<=|<|>|==|!=)<num>|smooth<SD>|exp|pow<num>|default<num>}; separated by commas. \n  log(2|10|) = apply log2, log10, or natural log;\n  [+-*/^]<num> = <data> <operation> <num>, where ^ is the exponent operator;\n  (>=|<=|<|>|==|!=)<num> = inequalities, thereby turning the data into 1s and 0s;\n  smoothG<num> = gaussian smooth, SD = <num>;\n  exp = e^data; \n  smoothU<num> = uniform (sliding window) smoothing with window = <num>;\n  smoothZ<num> = transform to Z scores, with SD and mean defined over window <num>/2 to either side of datum;\n  pow<num> = <num>^data;\n  default<num> = replace missing values with <num>;', required=True);
parser.add_argument('-c',dest='chrsFile',	metavar='<chromSizesFile>',help='Input file containing the chromsome sizes: chr\\tsize\\n', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFilePath>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output messages', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
parser.add_argument('-k',dest='chunks', metavar='<chunkSize>',help='The size of the chunks (in bp) to read/write when outputting [default=500000]', default=500000, required=False);
args = parser.parse_args();

#import itertools
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
allFunctions = args.functions.split(",");

additionalFlankSize = 0;
for f in allFunctions:
		m = re.match("^smoothG([.-eE0-9]*)$",f);
		if m:
			additionalFlankSize = np.max([additionalFlankSize,float(m.group(1)*4)])
		m = re.match("^smoothU([.-eE0-9]*)$",f);
		if m:
			additionalFlankSize = np.max([additionalFlankSize,int(m.group(1))])
		m = re.match("^smoothZ([.-eE0-9]*)$",f);
		if m:
			additionalFlankSize = np.max([additionalFlankSize,int(m.group(1))])
additionalFlankSize = int(additionalFlankSize);


def z_scoreize_slow(data,window):
	window = int(window/2)# window in either direction
	newData = np.zeros(data.shape);
	for i in range(0,np.size(data)):
		newData[i] = (data[i]-np.mean(data[np.max([0,i-window]):np.min([np.size(data),i+window])]))/np.std(data[np.max([0,i-window]):np.min([np.size(data),i+window])]);
	return newData;

def z_scoreize_online(data,window):
	window = int(window/2)# window in either direction
	mean=0;
	var=0;
	newData = np.zeros(data.shape);
	mean = np.mean(data[0:window]);
	var = np.var(data[0:window]);
	lastMean=mean;
	n=window;
	for i in range(0,window+1): #add things until i==window
		n=n+1
		mean = mean + (data[i+window]-mean)/(n);
		var = ((n-1)*var + (data[i+window]-lastMean)*(data[i+window]-mean))/(n);
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
	if n!=window*2+1: raise Exception("n not equal to 2*window+1!!: %i\n"%(n));
	for i in range(window+1,np.size(data)-window): #in the middle, add and subtract until i==len-window-1
		mean = mean + (data[i+window]-data[i-window-1])/(n);
		var = var + ((data[i+window]**2 - data[i-window-1]**2 + (lastMean + mean)*(data[i-window-1] - data[i+window]))/(n)) 
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
	for i in range(np.size(data)-window,np.size(data)): # in the end, subtract until i= len-1
		n=n-1;
		mean = mean - (data[i-window-1]-mean)/(n);
		var = ((n+1)*var - (data[i-window-1]-lastMean)*(data[i-window-1]-mean))/(n);
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
		
	if n!=window+1: raise Exception("n not equal to window+1!!: %i\n"%(n));
	return newData;

def z_scoreize_online2(data,window):
	window = int(window/2)# window in either direction
	mean=0;
	var=0;
	newData = np.zeros(data.shape);
	mean = np.mean(data[0:window]);
	var = np.var(data[0:window]);
	lastMean=mean;
	for i in range(0,window+1): #add things until i==window
		mean = mean + (data[i+window]-mean)/(i+window+1);
		var = ((i+window-1)*var + (data[i+window]-lastMean)*(data[i+window]-mean))/(i+window);
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
	for i in range(window+1,np.size(data)-window): #in the middle, add and subtract until i==len-window-1
		mean = mean + (data[i+window]-data[i-window-1])/(2*window+1);
		var = var + ((data[i+window]**2 - data[i-window-1]**2 + (lastMean + mean)*(data[i-window-1] - data[i+window]))/(2*window+1)) 
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
	for i in range(np.size(data)-window,np.size(data)): # in the end, subtract until i= len-1
		mean = mean - (data[i-window-1]-mean)/(window+(np.size(data)-i));
		var = ((np.size(data)-i+window+1)*var - (data[i-window-1]-lastMean)*(data[i-window-1]-mean))/(np.size(data)-i+window);
		lastMean=mean;
		newData[i] = (data[i]-mean)/np.sqrt(var);
		
	return newData;

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
		m = re.match("^([<>=!]+)([-.eE0-9]+)$",f);
		if m:
			if m.group(1)=="<=":
				return (data<=float(m.group(2))).astype(int);
			elif m.group(1)==">=":
				return (data>=float(m.group(2))).astype(int);
			elif m.group(1)=="<":
				return (data<float(m.group(2))).astype(int);
			elif m.group(1)==">":
				return (data>float(m.group(2))).astype(int);
			elif m.group(1)=="!=":
				return (data!=float(m.group(2))).astype(int);
			elif m.group(1)=="==":
				return (data==float(m.group(2))).astype(int);
		m = re.match("^([\^\*\+\-\/])([-.eE0-9]+)$",f);
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
		m = re.match("^default([-.eE0-9]*)$",f);
		if m:
			invalid = np.isnan( data )
			data2 = data;
			data2[ invalid ] = float(m.group(1));
			return data2;
		m = re.match("^smoothG([-.eE0-9]*)$",f);
		if m:
			return gaussian_filter(data, float(m.group(1)), mode='reflect', truncate=4.0)
		m = re.match("^smoothU([-.eE0-9]*)$",f);
		if m:
			return uniform_filter1d(data, int(m.group(1)), mode='reflect')
		m = re.match("^smoothZ([-.eE0-9]*)$",f);
		if m:
			return z_scoreize_online(data, int(m.group(1)))
		m = re.match("^pow([-.eE0-9]*)$",f);
		if m:
			return np.power(float(m.group(1)),data)
		raise Exception("Unknown function provided: %s"%f)
#data = np.array(range(1,10));
#for f in allFunctions:
#	data = applyFunction(data,f);
#	print(data);

#print(allFunctions);
#exit;

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
	last = 0;
	final = chromSizes[chr];
	sys.stderr.write("Outputting data for %s:\n"%(chr));
	while last!=final: # this breaks it up into chunks so that I'm not piping entire (human) chromosomes at once
		if args.verbose>0: sys.stderr.write("  Section %i - %i:\n"%(last,curLast));
		curLast = np.min([last+args.chunks,final]);
		curEnd = np.min([curLast+additionalFlankSize, final]);
		curSt = np.max([last-additionalFlankSize,0]);
		values = curBW.get_as_array( chr, curSt, curEnd )
	#print(chr);
		if values is not None:
			for f in allFunctions:
				values = applyFunction(values,f);
			values = values[(last - curSt):(curLast-last + (last-curSt))];# set them only to the middle part of this data so that the additionalFlankSize regions are not output.
			#print(values.shape);
			try:
				if last==0:
					toBW.stdin.write("fixedStep chrom=%s start=1 step=1\n"%(chr))
				toBW.stdin.write("\n".join(map(str,values)));
				toBW.stdin.write("\n");
				toBW.stdin.flush();
			except IOError as e:
				sys.stderr.write("IOError on %s, component %i\n"%(chr,i))
				sys.stderr.write("If this is a broken pipe, try reducing -k <chunks> from %i\n"%(args.chunks))
				raise(e);
		last=curLast;

toBW.stdin.close();
temp = toBW.communicate();

if temp[0] is not None:
	sys.stderr.write(temp[0]);

if temp[1] is not None:
	sys.stderr.write(temp[1]);

sys.stderr.write("Output all data!\n");
#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
if (args.logFP is not None):
	logFile.close();


