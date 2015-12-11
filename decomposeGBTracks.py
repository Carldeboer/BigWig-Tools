#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3


#import sys; sys.argv= "decomposeGBTracks.py -a ICA  -i test_samples.txt -o test_ICA -c /home/unix/cgdeboer/genomes/sc/20110203_R64/chrom.sizes -x 2 -v -v -v".split();

import argparse
parser = argparse.ArgumentParser(description='This program takes GB tracks as input and decomposes them using either PCA, ICA, or NMF, outputting the components.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file containing a list of files, one per line with the following format, tab separated: \n<ID>\t<filePath>\t<gaussianSmoothingSD>\t<defaultValue>\t<doLog?>\n where <doLog> is the number you want added to the data before log transform, or a negative number  otherwise', required=True);
parser.add_argument('-c',dest='chrsFile', metavar='<chrom.sizes>',help='A chrom.sizes file contiaining the sizes for all the chromosomes you want output', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePre>',help='Where to output results, prefix [default=stdout]\n  Multiple files will be created, including a BW track, and diagnostics showing the clustering of data', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-s',dest='sample', metavar='<sampleFrac>',help='The fraction of data points to use [default=0.01]', required=False, default = 0.01);
parser.add_argument('-a',dest='approach', metavar='<approach>',help='Approach to use; on eof {PCA,NMF,ICA} [default=PCA]', required=False, default = 'PCA');
parser.add_argument('-t',dest='iterations', metavar='<numIterations>',help='The number of iterations of NMF [default=1000]', required=False, default = 1000);
parser.add_argument('-r',dest='scaleTo', metavar='<scaleRangeTo>',help='How to scale the range of each track (none|SD|mean|max|median) [default=SD]', required=False, default="SD");
parser.add_argument('-x',dest='components', metavar='<numComponents>',help='The number of components to output [for PCA, defaults to all significant; required for NMF and ICA]', required=False,default=-1);
parser.add_argument('-k',dest='chunks', metavar='<chunkSize>',help='The size of the chunks (in bp) to read/write when outputting [default=500000]', default=500000, required=False);
parser.add_argument('-d',dest='dim', metavar='<graphDim>',help='Inches per graph [default=3]', required=False, default = 3);
parser.add_argument('-e',dest='eliminateMissing', action='count',help='eliminate base pairs with missing values; often, these represent 0s, so this is not always appropriate', required=False, default=0);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

import os
import itertools
import warnings
import subprocess
import MYUTILS
import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter;
#from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import decomposition
import sys
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile

args.sample = float(args.sample);
args.dim = float(args.dim);
args.iterations = int(args.iterations);
args.components = int(args.components);
args.chunks = int(args.chunks);
args.scaleTo = args.scaleTo.upper();

if args.approach!="PCA" and args.components==-1:
	raise Exception("components has no default for NMF/ICA - must specify with -x");

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
for chr in chrOrder:
	#sample positions
	useThese[chr] = np.random.random_sample((chromSizes[chr]))<args.sample;
	totalLength = totalLength + np.sum(useThese[chr]);


#make a matrix of the data
allData = np.empty([totalLength,len(IDs)]);
if args.eliminateMissing>0:
	keepThese = np.ones([totalLength]).astype(bool); # onlt those for which data was observed in all tracks

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
			allData = np.delete(allData, [range(curTot, (curTot+np.sum(useThese[chr])))],0);
			if args.eliminateMissing>0:
				keepThese = np.delete(keepThese, [range(curTot, (curTot+np.sum(useThese[chr])))],0);
			totalLength = totalLength -  np.sum(useThese[chr]);
			del useThese[chr]
			del chromSizes[chr]
			continue
		if args.verbose>2: sys.stderr.write("    Checking for missing data.\n");
		if args.eliminateMissing>0:
			#keepThese[curTot:(curTot+sum(useThese[chr]))] = np.logical_and(keepThese[curTot:(curTot+sum(useThese[chr]))], np.logical_not(np.isnan( values ))[useThese[chr]]);
			keepThese[np.add(curTot,np.nonzero(np.isnan( values[useThese[chr]])))] = False;
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
		values = values [useThese[chr]]
		#append data;
		if args.verbose>2: sys.stderr.write("    Appending data.\n");
		allData[curTot:(curTot+len(values)),i] = values;
		curTot = curTot + len(values);


if args.eliminateMissing>0:
	totalLength=np.sum(keepThese);
	allData = allData[keepThese,:];
	if args.verbose>0: sys.stderr.write("Selected %i, ended up with %i, will only keep %i for missing data.\n"%(totalLength,allData.shape[0], np.sum(keepThese)));

dataMedians = np.median(allData,axis=0);
sys.stderr.write("Decomposing data based on %i sampled values\n"%(allData.shape[0]))

#scale allData such that the tracks are on a similar scale
if args.scaleTo=="MEAN":
	scaleTo = np.mean(allData,axis=0);
elif args.scaleTo=="MAX":
	scaleTo = np.max(allData,axis=0);
elif args.scaleTo=="MEDIAN":
	scaleTo = np.median(allData,axis=0);
elif args.scaleTo=="NONE":
	scaleTo = [1]*len(IDs);
elif args.scaleTo=="SD":
	scaleTo = np.std(allData,axis=0);
else:
	raise Exception("Unrecognized scale: %s\n"%(args.scaleTo));


if args.approach=="ICA":
	np.divide(allData, np.array([scaleTo,]*allData.shape[0]), allData);# -> allData
	myDecomp = decomposition.FastICA(n_components=args.components, whiten=True, max_iter=args.iterations);
	myDecomp.fit(allData);
	W = myDecomp.mixing_;
	W = np.transpose(W);
	compFile = open("%s_ica_mix.txt"%(args.outFPre),"w");
	compFile.write("track\tcomponent_"+"\tcomponent_".join(map(str, range(0,args.components)))+"\n");
	for i in range(0,len(IDs)):
		compFile.write(IDs[i]+"\t"+"\t".join(map(str,W[:,i]))+"\n");
	
	compFile.close()

elif args.approach=="NMF":
	np.divide(allData, np.array([scaleTo,]*allData.shape[0]), allData);# -> allData
	#TODO: check for negative values and die if they are present, telling user which data have negative values
	#n_components = number of components to keep
	#init= how to initialize
	#sparseness {data|components|None]: where to enforce sparsity
	#beta: degree of sparseness (larger is more sparse, default=1
	#eta: degree of correctness to maintain; smaller=more error; default=0.1
	#tol: tolerance in stopping
	#max_iter: number of iterations
	#random_state: seed
	myDecomp = decomposition.NMF(sparseness='components', n_components=args.components, max_iter=args.iterations);
	myDecomp.fit(allData);
	
	#output stats to do with the fit
	sys.stderr.write("Done NMF: reconstruction error = %g\n"%(myDecomp.reconstruction_err_)); #Frobenius norm of the matrix differences between the training data and the reconstructed data; something like a sqrt(sum of squared error)
	# not really sure what the reconstruction error means in real terms since it depends on the amplitudes of the tracks
	W = myDecomp.components_ #non-negative components of the data; these correspond to the weights of each of the components 
	# the components is a (components x tracks) matrix, corresponding to W in W*H =~ V, where V is the input GB tracks
	#Therefore, the entries of W indicate how much each component contributes to each track.

	compFile = open("%s_component_weights.txt"%(args.outFPre),"w");
	compFile.write("track\tcomponent_"+"\tcomponent_".join(map(str, range(0,args.components)))+"\n");
	for i in range(0,len(IDs)):
		compFile.write(IDs[i]+"\t"+"\t".join(map(str,W[:,i]))+"\n");
	
	compFile.close()


else: #PCA
	#perform PCA to find axis of maximum variation
	#center and scale data
	dataMeans = np.mean(allData,axis=0);
	allData = np.divide(np.subtract(allData,[dataMeans,]*totalLength), [scaleTo,]*totalLength);
	
	myDecomp = decomposition.PCA()
	myDecomp.fit(allData);
	#myDecomp.components_
	if args.components<=0:
		args.components = myDecomp.n_components_;
	
	W = myDecomp.components_[0:args.components,:];
	sys.stderr.write("Done PCA: num components = %i. first PC explains %f%% of the variance \n"%(myDecomp.n_components_, myDecomp.explained_variance_ratio_[0]*100))
	
	#PVE plot
	plt.figure(figsize=(8,8),dpi=300);
	plt.bar(range(1,len(myDecomp.explained_variance_ratio_)+1),myDecomp.explained_variance_ratio_, 0.5, color='b',label=None);
	plt.xlabel('Component');
	plt.ylabel('PVE');
	plt.savefig('%s_PCA-PVE.png'%(args.outFPre))
	plt.savefig('%s_PCA-PVE.pdf'%(args.outFPre))
	
transformedData = myDecomp.transform(allData);
componentNames = np.char.mod('comp_%d',range(1,transformedData.shape[1]+1));

colour = "black";
gridSizeC = [np.min((32767,args.components*args.dim)), np.min((32767,args.components*args.dim))];

plt.figure(figsize=gridSizeC,dpi=300);
for x in range(0,args.components):
	for y in range(0,args.components): # each pair of inputs
		splot = plt.subplot(args.components,args.components,x*args.components+y+1)
		xUp = args.components-x-1;
		plt.scatter(W[y,:], W[xUp,:],color = colour, alpha=1)
		for i, lab in enumerate(IDs):
			plt.gca().annotate(lab,(W[y,i],W[xUp,i]));
		if y==0:
			plt.ylabel(componentNames[xUp]);
		if x==(args.components-1):
			plt.xlabel(componentNames[y]);

plt.savefig('%s_track_weight_scatters.png'%(args.outFPre))
plt.savefig('%s_track_weight_scatters.pdf'%(args.outFPre))

colour = "black";
plt.figure(figsize=gridSizeC,dpi=300);
for x in range(0,args.components):
	for y in range(0,args.components): # each pair of inputs
		xUp = args.components-x-1;
		splot = plt.subplot(args.components,args.components,x*args.components+y+1)
		if xUp==y:
			plt.hist(transformedData[:,y],100,facecolor = colour);
		else:
			plt.scatter(transformedData[:,y], transformedData[:,xUp],color = colour, alpha=0.1);
		if y==0:
			plt.ylabel(componentNames[xUp]);
		if x==(args.components-1):
			plt.xlabel(componentNames[y]);

plt.savefig('%s_all_comp_scatter.png'%(args.outFPre))

#plot first two components
plt.figure(figsize=(8,8),dpi=300);
plt.scatter(transformedData[:,0],transformedData[:,1],color = "b", alpha=0.1)
plt.xlabel(componentNames[0]);
plt.ylabel(componentNames[1]);
plt.savefig('%s_top2Components.png'%(args.outFPre))


corMat = np.empty([allData.shape[1],transformedData.shape[1]]);
for i in range(0,allData.shape[1]):
	for j in range(0,transformedData.shape[1]):
		corMat[i,j] = np.corrcoef(allData[:,i],transformedData[:,j])[1,0];

plt.figure(figsize=(10,8),dpi=300);
plt.pcolor(corMat, cmap=plt.cm.RdBu, vmin=-np.max(np.abs(corMat)), vmax=np.max(np.abs(corMat)));
plt.axis([0, transformedData.shape[1], 0, allData.shape[1]]) #x x, y y
plt.colorbar()
ax = plt.gca()
ax.set_xticks(np.arange(corMat.shape[1])+0.5, minor=False)
ax.set_yticks(np.arange(corMat.shape[0])+0.5, minor=False)
#ax.invert_yaxis() # because labels 
#ax.xaxis.tick_top()
ax.set_yticklabels(IDs, minor=False)
ax.set_xticklabels(componentNames, minor=False)
plt.xticks(rotation=90)
#plt.gca().tight_layout()
plt.gcf().subplots_adjust(left=0.3,bottom=0.15)
plt.savefig('%s_track_corrs.png'%(args.outFPre))
plt.savefig('%s_track_corrs.pdf'%(args.outFPre))

allData=None;
transformedData=None;

allInBWs = [];
for i in range(0,len(IDs)):
	#input GB tracks
	curBW = BigWigFile(open(files[i]))
	allInBWs.append(curBW);

sys.stderr.write("Making output streams\n");

outStreams = [];
for i in range(0, args.components):
	if args.verbose>1: sys.stderr.write("Making wig output file:  %s_comp%i.wig.gz.\n"%(args.outFPre,i+1));
	outStreams.append( MYUTILS.smartGZOpen("%s_comp%i.wig.gz"%(args.outFPre,i+1),"w") );

for i in range(0, args.components):
	outStreams[i].write("track type=wiggle_0\n")

for chr in chrOrder:
	sys.stderr.write("Outputting components for %s:\n"%(chr));
	last = 0;
	final = chromSizes[chr];
	while last!=final: # this breaks it up into chunks so that I'm not piping entire (human) chromosomes at once
		curLast = np.min([last+args.chunks,final]);
		chrData  = np.empty([curLast-last,len(IDs)]);
		if args.verbose>0: sys.stderr.write("  %s: Section %i - %i:\n"%(chr, last,curLast));
		for i in range(0,len(IDs)):
			if args.verbose>1: sys.stderr.write("    Input %i: %s.\n"%(i,IDs[i]));
			values = allInBWs[i].get_as_array( chr, last, curLast )
			#fill in blanks
			invalid = np.isnan( values )
			values[ invalid ] = dataMedians[i]; # at this stage, replace missing values with median
			#log transform
			if doLog[i]>=0:
				values =np.log10(values + doLog[i])
			#smooth data
			if smoothings[i]!=0:
				#reflect = at edge of array, the data will be mirrored
				#truncate = don't incorporate data more than X away - probably a great speed increase since the distrib goes out to infinity
				values = gaussian_filter(values, smoothings[i], mode='reflect', truncate=4.0)
			#add data;
			chrData[:,i] = np.divide(values,scaleTo[i]); # put on same scale as before
		transed = myDecomp.transform(chrData);
		for i in range(0,args.components):
			if args.verbose>1: sys.stderr.write("    Printing component %i.\n"%(i+1));
			try:
				if last==0:
					outStreams[i].write("fixedStep chrom=%s start=1 step=1\n"%(chr))
				outStreams[i].write("\n".join(map(str,transed[:,i])));
				outStreams[i].write("\n");
				outStreams[i].flush();
			except IOError as e:
				sys.stderr.write("IOError on %s, component %i\n"%(chr,i))
				sys.stderr.write("If this is a broken pipe, try reducing -k <chunks> from %i or increasing the amount of available RAM\n"%(args.chunks))
				raise(e);
		last=curLast;


sys.stderr.write("Output all data.\n");

for i in range(0,args.components):
	outStreams[i].close()
	toBW = subprocess.Popen(["wigToBigWig","%s_comp%i.wig.gz"%(args.outFPre,i+1),args.chrsFile,"%s_comp%i.bw"%(args.outFPre,i+1)])
	temp = toBW.communicate()
	if temp[0] is not None:
		sys.stderr.write("component %i : %s"%(i+1,temp[0]));
	if temp[1] is not None:
		sys.stderr.write("component %i : %s"%(i+1,temp[1]));
	if temp[0] is None and temp[1] is None and os.path.isfile("%s_comp%i.bw"%(args.outFPre,i+1)): # if no errors, delete the original
		os.remove("%s_comp%i.wig.gz"%(args.outFPre,i+1))
	else:
		sys.stderr.write("Left %s_comp%i.bw because of error "%(args.outFPre,i+1));

if args.verbose>0: sys.stderr.write("Successfully closed wigToBigWig processes.\n");

if (args.logFP is not None):
	logFile.close();
