#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3


#import sys; sys.argv= "findCopyNumber.py  -i example_input_CopyNumber.txt -o test_CNV -c /home/unix/cgdeboer/genomes/sc/20110203_R64/chrom.sizes -s 1 -v -v -v -w".split();

import argparse
parser = argparse.ArgumentParser(description='This program takes GB tracks as input and decomposes them using either PCA or NMF, outputting the components.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file containing a list of files, one per line with the following format, tab separated: \n<ID>\t<filePath>\t<gaussianSmoothingSD>\t<isOpen>\t<defaultValue>\t<doLog?>\n  where <isOpen> is 1 for openness tracks (e.g. FAIRE, DHS, ATAC-seq) and -1 for occupancy (e.g. nucleosomes)\n  and <doLog> is the number you want added to the data before log transform, or a negative number  otherwise', required=True);
parser.add_argument('-c',dest='chrsFile', metavar='<chrom.sizes>',help='A chrom.sizes file contiaining the sizes for all the chromosomes you want output', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePre>',help='Where to output results, prefix [default=stdout]\n  Multiple files will be created, including a BW track, and diagnostics showing the clustering of data', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-s',dest='sample', metavar='<sampleSpacing>',help='The spacing of data points to use (e.g. take one every X bp) [default=100]', required=False, default = 100);
parser.add_argument('-t',dest='iterations', metavar='<numIterations>',help='The number of iterations of ML-reestimation [default=10]', required=False, default = 10);
parser.add_argument('-p',dest='ploidy', metavar='<defaultPloidy>',help='The default ploidy (usually 1 for haploid, 2 for diploid) [default=1]', required=False,default=1);
parser.add_argument('-x',dest='transition', metavar='<log10TransitionP>',help='The log10(P(transition)) between states [default=-10]', required=False,default=-10);
parser.add_argument('-d',dest='dim', metavar='<graphDim>',help='Inches per graph [default=3]', required=False, default = 3);
parser.add_argument('-e',dest='eliminateMissing', action='count',help='eliminate base pairs with missing values; often, these represent 0s, so this is not always appropriate', required=False, default=0);
parser.add_argument('-w',dest='wigs', action='count',help='Output to wigs first? Useful if RAM is limiting.  Sometimes wigToBigWig gets killed for it and its less RAM intensive to convert everything at the end', required=False, default=0);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

import os
import itertools
import warnings
import subprocess
import MYUTILS
from  hmmlearn.hmm import GaussianHMM;
import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter;
#from scipy.stats import multivariate_normal;
#from scipy.stats import boxcox;
#from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import decomposition
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

args.sample = int(args.sample);
args.dim = float(args.dim);
args.iterations = int(args.iterations);
args.transition = int(args.transition);


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;


IDs =[];
files = [];
smoothings = [];
defaultVal = [];
doLog = [];
inFile=MYUTILS.smartGZOpen(args.inFP,'r');
for line in inFile:
	if line is None or line == "" or line.rstrip()=="" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	if len(data)!=5: raise Exception("Incorrect number of fields in input: %s\n"%(line));
	IDs.append(data[0]);
	files.append(data[1]);
	smoothings.append(float(data[2]));
	defaultVal.append(float(data[3]));
	doLog.append(float(data[4]));

inFile.close();

#get loci of interest and their sizes
chromSizesFile = MYUTILS.smartGZOpen(args.chrsFile,'r');
chromSizes = {};
chrOrder = [];
for line in chromSizesFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	chromSizes[data[0]]=int(data[1]);
	chrOrder.append(data[0]);
	

#determine what positions will be used in the data
useThese = {};
totalLength = 0
allData = {};
if args.eliminateMissing>0:
	keepThese = {};

for chr in chrOrder:
	#sample positions
	useThese[chr] = np.arange(0,chromSizes[chr]-1,args.sample,int)
	totalLength = totalLength + useThese[chr].shape[0];
	allData[chr] = np.empty([useThese[chr].shape[0],len(IDs)]);
	if args.eliminateMissing>0:
		keepThese[chr] = np.ones([useThese[chr].shape[0]]).astype(bool); # onlt those for which data was observed in all tracks


#make a matrix of the data

for i in range(0,len(IDs)):
	#input GB tracks
	curBW = BigWigFile(open(files[i]))
	curTot = 0;
	if args.verbose>1: sys.stderr.write("Inputting data for %s.\n"%(IDs[i]));
	for chr in chrOrder:
		if args.verbose>1: sys.stderr.write("  Inputting data for %s.\n"%(chr));
		if args.verbose>2: sys.stderr.write("    Getting data from BW.\n");
		values = curBW.get_as_array( chr, 0, chromSizes[chr] )
		if values is None:
			sys.stderr.write("%s is missing %s... skipping it for all\n"%(IDs[i],chr));
			chrOrder.remove(chr)
			totalLength = totalLength -  np.sum(useThese[chr]);
			del allData[chr]
			if args.eliminateMissing>0:
				del keepThese[chr];
			del useThese[chr]
			del chromSizes[chr]
			continue
		if args.verbose>2: sys.stderr.write("    Checking for missing data.\n");
		if args.eliminateMissing>0:
			#keepThese[curTot:(curTot+sum(useThese[chr]))] = np.logical_and(keepThese[curTot:(curTot+sum(useThese[chr]))], np.logical_not(np.isnan( values ))[useThese[chr]]);
			keepThese[chr][np.nonzero(np.isnan( values[useThese[chr]]))] = False;
		#print(np.add(curTot,np.nonzero(np.isnan( values[useThese[chr]]))))
		#fill in blanks; these are only used for smoothing anyways as keepThese will filter them out.
		if args.verbose>2: sys.stderr.write("    Filling in missing data.\n");
		invalid = np.isnan( values )
		values[ invalid ] = defaultVal[i];
		#log transform
		if doLog[i]>=0:
			if args.verbose>2: sys.stderr.write("    Log transforming data.\n");
			values =np.log10(values + doLog[i])
		#smooth data
		if smoothings[i]!=0:
			#reflect = at edge of array, the data will be mirrored
			#truncate = don't incorporate data more than X away - probably a great speed increase since the distrib goes out to infinity
			if args.verbose>2: sys.stderr.write("    Smoothing data.\n");
			values = gaussian_filter(values, smoothings[i], mode='reflect', truncate=4.0)
		#sample data
		if args.verbose>2: sys.stderr.write("    Sampling data.\n");
		#TODO: instead, calculate the mean of each sample range and take that
		values = values [useThese[chr]]
		#append data;
		if args.verbose>2: sys.stderr.write("    Appending data.\n");
		allData[chr][:,i] = values;
		curTot = curTot + len(values);

allDataCat = np.empty([curTot,len(IDs)]);
curTot=0;
for chr in chrOrder:
	if args.eliminateMissing>0:
		curLen = np.sum(keepThese[chr]);
		allDataCat[curTot:(curTot+curLen),:] = allData[chr][keepThese[chr],:];
		useThese[chr] = useThese[chr][keepThese[chr]];
	else:
		curLen = allData[chr].shape[0];
		allDataCat[curTot:(curTot+curLen),:] = allData[chr];
	curTot=curTot+curLen;

if args.eliminateMissing>0:
	allDataCat = allDataCat[0:curTot,:];
	

if 0: #don't do this because it is non-linear and I can't assume that CNVs will scale linearly within the transformation
	## perform box-cox transformation of data
	if args.verbose>0: sys.stderr.write("Performing BoxCox transformations.\n");
	allDataCat = allDataCat+1;
	bcLambdas = np.zeros(len(IDs));
	for i in range(0, len(IDs)):
		allDataCat[:,i], bcLambdas[i]= boxcox(allDataCat[:,i]);
		#TODO: print a histogram of the transformed data

if args.verbose>0: sys.stderr.write("Building HMM.\n");

#1. Define PDF of signal with genome-wide data
covAll = np.cov(allDataCat,rowvar=0);
meanAll = np.mean(allDataCat,axis=0);

#2. Define three (or two for haploid) state HMM with means equal to mean, mean/2, and mean*1.5 (or just mean, mean*2 for haploid), user-defined transition probability
numStates = args.ploidy+2
stateMeans = [0]*numStates;
stateCovs = [0]*numStates;
stateIsToCNVs = [];
cnvsToStateIs = {};
for i in range(0,numStates):
	stateMeans[i] = float(i)/args.ploidy *meanAll;
	stateCovs[i] = covAll;
	stateIsToCNVs.append(i);
	cnvsToStateIs[i]=i


if len(IDs)==1:
	stateCovs = np.expand_dims(stateCovs,1)

### make transmat
transitionMatrix = np.add(np.eye(numStates)*(1-(numStates-1)*10**args.transition),(1-np.eye(numStates))*10**args.transition);

#model = GaussianHMM(len(states),covariance_type="full",n_iter=1);
model = GaussianHMM(numStates,covariance_type="full", n_iter=1);
###insert my own params
model.means_ = stateMeans;
model.covars_ = stateCovs;
model.transmat=transitionMatrix;


meanNormal = meanAll;
normalState = cnvsToStateIs[args.ploidy];
lastClass = {};
for chr in chrOrder:
	lastClass[chr] = np.tile(args.ploidy,allData[chr].shape[0]);


for i in range(0,args.iterations):
	warned=False;
	if args.verbose>1: sys.stderr.write("  Iteration %i.\n"%(i));
	viterbi = {}
	noChange=False;
	viterbiCat =np.zeros(allDataCat.shape[0]);
	curTot=0;
	curNumCNVs
	for chr in chrOrder:
		#3. Calculate Viterbi path given data
		if args.verbose>2: sys.stderr.write("    i=%i; Calculating Viterbi path for %s.\n"%(i,chr));
		viterbi[chr] = model.predict(allData[chr]);
		curLen = len(viterbi[chr]);
		#4. For each non-standard state, calculate the mean in that state and add a state with a mean representing that ploidy
		changeStart=-1
		viterbi[chr] = np.insert(viterbi[chr],[0,curLen],[normalState,normalState]); # add initial and terminal normalStates so that telomeres in CNV will be detected.
		for j in range(1,len(viterbi[chr])):
			if viterbi[chr][j]!=viterbi[chr][j-1]:#there was a change
				if changeStart==-1:
					changeStart=j;
				else: #from changeStart to j-1
					curNumCNVs+=1;
					#calculate the means of this region
					localMean = np.mean(allData[chr][changeStart:j,:],axis=0);
					#figure out the local CN as the local means divided by the global means, rounded to the nearest logical ploidy
					meanRatio = np.divide(localMean,meanNormal);
					localPloidy = np.mean(meanRatio*args.ploidy);
					if args.verbose>2: sys.stderr.write("    Local ploidy for %s (%i:%i)*%i equal to %f (rounded to %i).\n"%(chr,changeStart,j-1,args.sample,localPloidy,int(np.round(localPloidy))));
					#re-assign viterbi to have be this (potentially new) state
					localPloidy = int(np.round(localPloidy));
					if localPloidy==args.ploidy and not warned:
						warned=True;
						sys.stderr.write("WARNING: Local ploidy rounds to global ploidy; transition log-probability (%f) is probably too high\n"%(args.transition));
					if localPloidy not in cnvsToStateIs:
						cnvsToStateIs[localPloidy]=numStates;
						stateIsToCNVs.append(localPloidy);
						numStates+=1;
					viterbi[chr][changeStart:j]=cnvsToStateIs[localPloidy];
					changeStart=-1;	
		viterbi[chr] = viterbi[chr][1:len(viterbi[chr])-1] #remove the initial and terminal normalState
		viterbiCat[curTot:(curTot+curLen)] = viterbi[chr];
		curTot=curTot+curLen;
		noChange = noChange and np.all(viterbi[chr]==lastClass[chr]);
		lastClass[chr]=viterbi[chr]#update
	if noChange: #no chr was updated
		break;
	#5. Re-calculate the means and SDs of each state given the viterbi path
	meanNormal = np.mean(allDataCat[viterbiCat==cnvsToStateIs[args.ploidy],:],axis=0)
	covNormal = np.cov(allDataCat[viterbiCat==cnvsToStateIs[args.ploidy],:],rowvar=0)
	stateMeans = [0]*numStates;
	stateCovs = [0]*numStates;
	for s in range(0,numStates):
		stateMeans[s] = meanNormal*(float(stateIsToCNVs[s])/args.ploidy);
		if np.sum(viterbiCat==s)<3:
			stateCovs[s] = covNormal;
		else:
			stateCovs[s] = np.cov(allDataCat[viterbiCat==s,:],rowvar=0)
	model = GaussianHMM(numStates,covariance_type="full", n_iter=1);
	transitionMatrix = np.add(np.eye(numStates)*(1-(numStates-1)*10**args.transition),(1-np.eye(numStates))*10**args.transition);
	model.transmat=transitionMatrix;
	model.means_ = stateMeans;
	model.covars_ = stateCovs;
	#6. Repeat steps 3-6 until convergence (E-M)




sys.stderr.write("Finished in %i/%i iterations, yielding %i CNVs\n",%(i+1, args.iterations, curNumCNVs));

sys.stderr.write("Making output streams\n");

outBW = []
outStream = [];
if args.wigs>0:
	if args.verbose>1: sys.stderr.write("Making wig output file:  %s.wig.gz.\n"%(args.outFPre));
	outStream = MYUTILS.smartGZOpen("%s.wig.gz"%(args.outFPre),"w");
else:
	if args.verbose>1: sys.stderr.write("Making wig->bigwig pipe for output file:  %s.bw.\n"%(args.outFPre));
	outBW = subprocess.Popen(["wigToBigWig","stdin",args.chrsFile,"%s.bw"%(args.outFPre)], stdin=subprocess.PIPE)
	outStream = toBW.stdin;

nbytes = outStream.write("track type=wiggle_0\n");

for chr in chrOrder:
	sys.stderr.write("Outputting inferred ploidy for %s:\n"%(chr));
	#get data for this chromosome
	try:		
		nBytes = outStream.write("variableStep chrom=%s span=%i\n"%(chr,args.sample))
		for i in range(0,len(lastClass[chr])):
			nBytes = outStream.write("%i\t%i\n"%(useThese[chr][i]+1,stateIsToCNVs[lastClass[chr][i]]));
		outStream.flush();
	except IOError as e:
		sys.stderr.write("IOError on %s, component %i\n"%(chr,i))
		sys.stderr.write("If this is a broken pipe, try using the -w option\n")
		raise(e);


sys.stderr.write("Output all data.\n");

if args.wigs>0:
	outStream.close()
	toBW = subprocess.Popen(["wigToBigWig","%s.wig.gz"%(args.outFPre),args.chrsFile,"%s.bw"%(args.outFPre)])
	temp = toBW.communicate()
	if temp[0] is not None:
		sys.stderr.write("wigToBigWig: %s"%(temp[0]));
	if temp[1] is not None:
		sys.stderr.write("wigToBigWig: %s"%(temp[1]));
	if temp[0] is None and temp[1] is None: # if no errors, delete the original
		os.remove("%s.wig.gz"%(args.outFPre))
else: #streams are pipes to wigToBigWig
	outBW.stdin.close();
	temp = outBW.communicate()
	if temp[0] is not None:
		sys.stderr.write("wigToBigWig: %s"%(temp[0]));
	if temp[1] is not None:
		sys.stderr.write("wigToBigWig: %s"%(temp[1]));
	if args.verbose>0: sys.stderr.write("Successfully closed wigToBigWig processes.\n");

if (args.logFP is not None):
	logFile.close();
