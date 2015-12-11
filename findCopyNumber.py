#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3


#import sys; sys.argv= "findCopyNumber.py  -i example_input_CopyNumber3.txt -o test_MNase3-5_rawPDF_x300 -c /home/unix/cgdeboer/genomes/sc/20110203_R64/chrom.sizes -x -300 -b 0.001 -s 1 -v -v -v -w".split();

import argparse
parser = argparse.ArgumentParser(description='This program takes GB tracks as input and decomposes them using either PCA or NMF, outputting the components.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file containing a list of files, one per line with the following format, tab separated: \n<ID>\t<filePath>\t<gaussianSmoothingSD>\t<defaultValue>\t<doLog?>\n where <doLog> is the number you want added to the data before log transform, or a negative number  otherwise', required=True);
parser.add_argument('-c',dest='chrsFile', metavar='<chrom.sizes>',help='A chrom.sizes file contiaining the sizes for all the chromosomes you want output', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePre>',help='Where to output results, prefix [default=stdout]\n  Multiple files will be created, including a BW track, and diagnostics showing the clustering of data', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-s',dest='sample', metavar='<sampleSpacing>',help='The spacing of data points to use (e.g. take one every X bp) [default=100]', required=False, default = 100);
parser.add_argument('-t',dest='iterations', metavar='<numIterations>',help='The number of iterations of ML-reestimation [default=10]', required=False, default = 10);
parser.add_argument('-p',dest='ploidy', metavar='<defaultPloidy>',help='The default ploidy (usually 1 for haploid, 2 for diploid) [default=1]', required=False,default=1);
parser.add_argument('-x',dest='transition', metavar='<log10TransitionP>',help='The log10(P(transition)) between states [default=-100]; more negative is more conservative (calls fewer CNVs) - more resolution (lower -s <sample>) requires a more negative transition', required=False,default=-100);
parser.add_argument('-b',dest='fractionBG', metavar='<backgroundAlignmentPct>',help='the fraction of the covariance of the standard ploidy state to use for the null state (no copies) - shouldn\'t be 0 since reads may align by chance [default=0.001]', required=False,default=0.001);
parser.add_argument('-r',dest='standardPrior', metavar='<standardPloidyPrior>',help='The prior (in log probability) preference for the default (as specified by -p) state [default=-1]', required=False,default=-0.1);
#parser.add_argument('-d',dest='dim', metavar='<graphDim>',help='Inches per graph [default=3]', required=False, default = 3);
parser.add_argument('-z',dest='scalePDF', action='count',help='Scale the state PDFs to have a max of 1; this has the effect of making variance less important than means and often results in a larger number of more specific states.', required=False, default=0);
parser.add_argument('-e',dest='eliminateMissing', action='count',help='eliminate base pairs with missing values; often, these represent 0s, so this is not always appropriate', required=False, default=0);
parser.add_argument('-ns',dest='notSumSamples', action='count',help='Do not sum adjacent data points when sampling data [e.g. if sample is 10, the first base will sum 1:10]', required=False, default=0);
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
from scipy.stats import multivariate_normal;
#from scipy.stats import boxcox;
#from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from sklearn import decomposition
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

args.sample = int(args.sample);
#args.dim = float(args.dim);
args.iterations = int(args.iterations);
args.fractionBG = float(args.fractionBG);
args.transition = float(args.transition);
args.standardPrior = float(args.standardPrior);
args.ploidy = int(args.ploidy);


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
	useThese[chr] = np.arange(0,np.floor(chromSizes[chr]/args.sample)*args.sample,args.sample,int) #skip the last position since it will have incomplete data
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
			sys.stderr.write("%s is missing %s... skipping it for all and removing %i from total length (now %i)\n"%(IDs[i],chr, len(useThese[chr]), totalLength -  useThese[chr].shape[0]));
			chrOrder.remove(chr)
			#totalLength = totalLength -  np.sum(useThese[chr]);
			totalLength = totalLength -  useThese[chr].shape[0];
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
		if np.max(useThese[chr])>=len(values):
			raise Exception("useThese[chr] contains indeces greater than len(values)! (%i vs %i)"%(np.max(useThese[chr]), len(values)));
		#TODO: instead, calculate the mean of each sample range and take that
		mergedValues = values [useThese[chr]]
		if args.notSumSamples==0:
			for si in range(1,args.sample):
				mergedSamples = values[useThese[chr]+si];
		#append data;
		if args.verbose>2: sys.stderr.write("    Appending data.\n");
		if allData[chr].shape[0]!=len(mergedValues):
			raise Exception("mergedValues and allData[chr] are not of equal extent! (%i vs %i)"%(len(mergedValues), allData[chr].shape[0]));
		allData[chr][:,i] = mergedValues;
		curTot = curTot + len(mergedValues);

if curTot!=totalLength:
	sys.stderr.write("curTot and totalLength are not of equal extent! (%i vs %i)\n"%(curTot,totalLength));
allDataCat = np.empty([totalLength,len(IDs)]);
curTot=0;
for chr in chrOrder:
	if args.eliminateMissing>0:
		curLen = np.sum(keepThese[chr]);
		allDataCat[curTot:(curTot+curLen),:] = allData[chr][keepThese[chr],:];
		useThese[chr] = useThese[chr][keepThese[chr]];
	else:
		curLen = allData[chr].shape[0];
		sys.stderr.write("Cating allData: curLen=%i, curTot=%i, allData[chr].shape=%s, allDataCat.shape = %s\n"%(curLen, curTot, str(allData[chr].shape), str(allDataCat.shape)));
		allDataCat[curTot:(curTot+curLen),:] = allData[chr];
	curTot=curTot+curLen;

sys.stderr.write("curTot=%i, totalLength=%i; resizing to curTot\n"%(curTot,totalLength));

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
stateIsToCNVs = range(0,(args.ploidy+2))
numStates = len(stateIsToCNVs)
stateMeans = [0]*numStates;
stateCovs = [0]*numStates;
statePDFMaxima = [0]*numStates;
cnvsToStateIs = {};
for i in range(0,numStates): #assume poisson, where mean is lambda and variance is also lambda
	stateMeans[i] = float(i)/args.ploidy *meanAll;
	if i==0:
		stateCovs[i] = covAll * args.fractionBG/args.ploidy; # since if the variance is 0, the probability of observing anything but the mean (0) is 0
	else:
		stateCovs[i] = covAll * float(i)/args.ploidy;
	cnvsToStateIs[i]=i
	statePDFMaxima[i]=np.log(multivariate_normal.pdf(x=stateMeans[i],mean=stateMeans[i],cov=stateCovs[i]))

cnvsToStateIs0=cnvsToStateIs;
stateIsToCNVs0 = stateIsToCNVs;
if len(IDs)==1:
	stateCovs = np.expand_dims(stateCovs,1)


#model = GaussianHMM(len(states),covariance_type="full",n_iter=1);
model = GaussianHMM(numStates,covariance_type="full", n_iter=1);
###insert my own params
model.means_ = stateMeans;
model.covars_ = stateCovs;

### make transmat
if args.transition <= -100:
	transitionMatrix = (1-np.eye(numStates))*args.transition*np.log(10);
	model._log_transmat =transitionMatrix;
else:
	transitionMatrix = np.add(np.eye(numStates)*(1-(numStates-1)*10**args.transition),(1-np.eye(numStates))*10**args.transition);
	model._set_transmat(transitionMatrix);

if args.verbose>0: sys.stderr.write(np.array_str(model._log_transmat)+"\n");

#exit;
meanNormal = meanAll;
normalState = cnvsToStateIs[args.ploidy];
lastClass = {};
for chr in chrOrder:
	lastClass[chr] = np.tile(args.ploidy,allData[chr].shape[0]);


for i in range(0,args.iterations):
	warned=False;
	if args.verbose>1: sys.stderr.write("  Iteration %i.\n"%(i));
	viterbi = {}
	noChange=0;
	viterbiCat =np.zeros(allDataCat.shape[0]);
	curTot=0;
	curNumCNVs=0;
	for chr in chrOrder:
		#3. Calculate Viterbi path given data
		if args.verbose>2: sys.stderr.write("    i=%i; Calculating Viterbi path for %s.\n"%(i,chr));
		framelogprob = model._compute_log_likelihood(allData[chr])
		#sys.stderr.write("framelogprob dim: "+str(framelogprob.shape)+"\n");
		framelogprob[:,cnvsToStateIs[args.ploidy]] = np.subtract(framelogprob[:,cnvsToStateIs[args.ploidy]], args.standardPrior); #add log(prior)
		if args.scalePDF>0:
			framelogprob = np.subtract(framelogprob,statePDFMaxima) #### This requires some explanation.  See Note 1 below. 
		logprob, viterbi[chr] = model._do_viterbi_pass(framelogprob);
		curLen = len(viterbi[chr]);
		#4. For each non-standard state, calculate the mean in that state and add a state with a mean representing that ploidy
		changeStart=-1
		viterbi[chr] = np.insert(viterbi[chr],[0,curLen],[normalState,normalState]); # add initial and terminal normalStates so that telomeres in CNV will be detected.
		for j in range(1,len(viterbi[chr])):
			if viterbi[chr][j]!=viterbi[chr][j-1]:#there was a change
				if changeStart==-1:
					if viterbi[chr][j]==normalState:
						raise Exception("new state is normal ploidy state");
					changeStart=j;
				else: #from changeStart to j-1
					#calculate the means of this region
					localMean = np.mean(allData[chr][changeStart:j,:],axis=0);
					#figure out the local CN as the local means divided by the global means, rounded to the nearest logical ploidy
					meanRatio = np.divide(localMean,meanNormal);
					localPloidy = np.mean(meanRatio*args.ploidy);
					if args.verbose>2: sys.stderr.write("      Local ploidy for %s:%i-%i i=(%i-%i)*%i (len=%ibp) better fit by %i, empirically equal to %f (rounded to %i).\n"%(chr,changeStart*args.sample,(j-1)*args.sample,changeStart,j-1,args.sample,(j-changeStart)*args.sample,stateIsToCNVs[viterbi[chr][j-1]],localPloidy,int(np.round(localPloidy))));
					#re-assign viterbi to have be this (potentially new) state
					localPloidy = int(np.round(localPloidy));
					if localPloidy==args.ploidy and not warned:
						warned=True;
					#	sys.stderr.write("WARNING: Local ploidy rounds to global ploidy; transition log-probability (%f) is probably too high\n"%(args.transition));
					else: 
						curNumCNVs+=1;
					if localPloidy not in cnvsToStateIs:
						cnvsToStateIs[localPloidy]=numStates;
						stateIsToCNVs.append(localPloidy);
						numStates+=1;
					viterbi[chr][changeStart:j]=cnvsToStateIs[localPloidy];
					changeStart=-1;	
					if viterbi[chr][j]!=normalState:
						changeStart=j;
		viterbi[chr] = viterbi[chr][1:len(viterbi[chr])-1] #remove the initial and terminal normalState
		#if args.verbose>1: sys.stderr.write("    merging revised viterbi with the rest: curTot=%i, curLen=%i, viterbi[chr].shape=%i, replacing %i:%i/%i.\n"%(curTot, curLen, viterbi[chr].shape[0], curTot, (curTot+curLen),allDataCat.shape[0]));
		viterbiCat[curTot:(curTot+curLen)] = viterbi[chr];
		curTot=curTot+curLen;
		noChange = noChange + np.sum(viterbi[chr]!=lastClass[chr]);
		lastClass[chr]=viterbi[chr]#update
	#5. Re-calculate the means and SDs of each state given the revised viterbi path
	#print(viterbiCat.dtype);
	viterbitCat = viterbiCat.astype(int);
	usedStates = np.unique(np.append(np.array([0,1]),np.unique(viterbiCat))).astype(int);
	numStates = len(usedStates);
	sys.stderr.write("Last iteration had %i CNVs and changed at %i positions; now have %i CNV states\n"%(curNumCNVs, noChange, numStates));
	if noChange==0: #no chr was updated
		break;
	#redefine standard ploidy
	meanNormal = np.mean(allDataCat[viterbiCat==cnvsToStateIs[args.ploidy],:],axis=0)
	covNormal = np.cov(allDataCat[viterbiCat==cnvsToStateIs[args.ploidy],:],rowvar=0)
	#reset states and remove unused states
	cnvsToStateIsNew = {};
	stateIsToCNVsNew =[0]*numStates;
	for s in range(0,numStates):
		cnvsToStateIsNew[stateIsToCNVs[usedStates[s]]]=s;
		stateIsToCNVsNew[s]=stateIsToCNVs[usedStates[s]];
	cnvsToStateIs=cnvsToStateIsNew;
	stateIsToCNVs =stateIsToCNVsNew;
	stateMeans = [0]*numStates;
	stateCovs = [0]*numStates;
	statePDFMaxima = [0]*numStates;
	for s in range(0,numStates):
		if args.verbose>1: sys.stderr.write("  Calculating mean and covariance for state i=%i: CN-%i; have %i examples\n"%(s,stateIsToCNVs[s], np.sum(viterbiCat==usedStates[s])));
		stateMeans[s] = meanNormal*(float(stateIsToCNVs[s])/args.ploidy);
		if np.sum(viterbiCat==usedStates[s])>=3: # too few examples to estimate covar ## sometimes estimating the covariance emperically leads to some weird states defined more by covariance than anything else.
			stateCovs[s] = np.cov(allDataCat[viterbiCat==usedStates[s],:],rowvar=0)
		redoCV = False;
		try: #test if positive definite
			temp = np.linalg.cholesky(stateCovs[s]);
		except np.linalg.LinAlgError as e:
			redoCV=True;
		if np.sum(viterbiCat==usedStates[s])<3 or redoCV: # too few examples to estimate covar or variance was 0 for at least one ### sometimes estimating the covariance emperically leads to some weird states defined more by covariance than anything else.

			# it appears to be better to estimate this from the empirical covariance (else, below) rather than the CNV* standard covariance as here
			if float(stateIsToCNVs[s])==0:
				stateCovs[s] = covNormal * args.fractionBG/args.ploidy; # since if the variance is 0, the probability of observing anything but the mean (0) is 0
			else:
				stateCovs[s] = covNormal * float(stateIsToCNVs[s])/args.ploidy;
			stateCovs[s] = covNormal;
		statePDFMaxima[s]=np.log(multivariate_normal.pdf(x=stateMeans[s],mean=stateMeans[s],cov=stateCovs[s]))
	model = GaussianHMM(numStates,covariance_type="full", n_iter=1);
	if args.transition <= -100:
		transitionMatrix = (1-np.eye(numStates))*args.transition*np.log(10);
		model._log_transmat =transitionMatrix;
	else:
		transitionMatrix = np.add(np.eye(numStates)*(1-(numStates-1)*10**args.transition),(1-np.eye(numStates))*10**args.transition);
		model._set_transmat(transitionMatrix);
	model.means_ = stateMeans;
	try:
		model.covars_ = stateCovs;
	except ValueError as e:
		print(stateCovs);
		raise(e);


statesOldToNew = {}
for s in range(0,numStates):
	statesOldToNew[usedStates[s]]=s;


sys.stderr.write("Finished in %i/%i iterations, yielding %i CNVs\n"%(i+1, args.iterations, curNumCNVs));

sys.stderr.write("Making output streams\n");

outStream = [];
if args.verbose>1: sys.stderr.write("Making wig output file:  %s.wig.gz.\n"%(args.outFPre));
outStream = MYUTILS.smartGZOpen("%s.wig.gz"%(args.outFPre),"w");

nbytes = outStream.write("track type=wiggle_0\n");

for chr in chrOrder:
	sys.stderr.write("Outputting inferred ploidy for %s:\n"%(chr));
	#get data for this chromosome
	try:		
		nBytes = outStream.write("variableStep chrom=%s span=%i\n"%(chr,args.sample))
		for i in range(0,len(lastClass[chr])):
			nBytes = outStream.write("%i\t%i\n"%(useThese[chr][i]+1,stateIsToCNVs[statesOldToNew[lastClass[chr][i]]])); #lastClass -> new stateIs -> CNVs
		outStream.flush();
	except IOError as e:
		sys.stderr.write("IOError on %s, component %i\n"%(chr,i))
		raise(e);



outStream.close()
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

if (args.logFP is not None):
	logFile.close();


### Note 1
# A PDF retruns the probability of observing the data for a given distribution.  Indeed, this is what is returned to framelogprob.
# This means that the maximum value of the function depends on the variance.  For example, if two states had the same mean, but different variance, 
#  what state is a value equal to the mean more likely to come from? Because of the increased variance of the one state, it is less likely to emit 
#  values close to its mean.  This kind of makes sense in the case where the variance of a state is informative.  However, take a second example:
#  state A has mean 0, SD=1, and state B has mean 1, SD=20.  Here, P(x=1) is much greater for state A than state B because of the vast difference in SD.
# Since the values of adjacent bases are not independent of one another, a whole string of 1s will greatly favour state A. Since adjacent bases are not
#  independent of one another (because reads span more than one base and few, if any, methods are truly single BP resolution, this is a likely scenario.
# Thus, dividing by the max of the PDF (subtracting the log(max(PDF(x)))) is equivalent to dividing the PDF by 1/(sigma*sqroot(2)*pi), making the PDF
#  exp(-(x-mu)^2/(2*sigma^2). Note that this transformation means that the AUPDF now increases with increasing variance.
# However, this also means that states with extreme variance can dominate predictions because there is little penalty from being far from the mean.
# I tried both ways and using the standard PDF seemed to work better. Scaling results in a large number of states being created.
