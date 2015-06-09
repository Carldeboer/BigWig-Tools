#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
#doesn't work with:
#!/home/unix/cgdeboer/bin/python3
import itertools
import warnings
import subprocess
import MYUTILS
import numpy as np
import scipy as sp
from scipy.ndimage.filters import gaussian_filter;
from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import decomposition
import sys
import argparse
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
parser = argparse.ArgumentParser(description='This program takes GB tracks as input and performs a NMF on it, outputting the components.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file containing a list of files, one per line with the following format, tab separated: \n<ID>\t<filePath>\t<gaussianSmoothingSD>\t<isOpen>\t<defaultValue>\t<doLog?>\n  where <isOpen> is 1 for openness tracks (e.g. FAIRE, DHS, ATAC-seq) and -1 for occupancy (e.g. nucleosomes)\n  and <doLog> is the number you want added to the data before log transform, or a negative number  otherwise', required=True);
parser.add_argument('-c',dest='chrsFile', metavar='<chrom.sizes>',help='A chrom.sizes file contiaining the sizes for all the chromosomes you want output', required=True);
parser.add_argument('-o',dest='outFPre', metavar='<outFilePre>',help='Where to output results, prefix [default=stdout]\n  Multiple files will be created, including a BW track, and diagnostics showing the clustering of data', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-s',dest='sample', metavar='<sampleFrac>',help='The fraction of data points to use [default=0.01]', required=False, default = 0.01);
parser.add_argument('-x',dest='components', metavar='<numComponents>',help='The number of components to output [default=0.01]', required=True);
parser.add_argument('-d',dest='dim', metavar='<graphDim>',help='Inches per graph [default=3]', required=False, default = 3);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);
args = parser.parse_args();

# this is for dubugging in interactive mode - comment out normally
args = lambda: None
setattr(args,"sample","0.01")
setattr(args,"dim","3")
setattr(args,"inFP","calcNMFTestData.txt")
setattr(args,"outFPre","calcNMFTestData_out")
setattr(args,"chrsFile","/home/unix/cgdeboer/genomes/sc/20110203_R64/chrom.sizes");
setattr(args,"components","3");
setattr(args,"logFP",None);
	

args.sample = float(args.sample);
args.dim = float(args.dim);
args.components = int(args.components);

if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;


IDs =[];
files = [];
smoothings = [];
isOpenness = [];
defaultVal = [];
doLog = [];
inFile=MYUTILS.smartGZOpen(args.inFP,'r');
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
	data=line.rstrip().split("\t");
	IDs.append(data[0]);
	files.append(data[1]);
	smoothings.append(float(data[2]));
	isOpenness.append(int(data[3]));
	defaultVal.append(float(data[4]));
	doLog.append(float(data[5]));

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
for i in range(0,len(IDs)):
	#input GB tracks
	curBW = BigWigFile(open(files[i]))
	curTot = 0;
	for chr in chrOrder:
		values = curBW.get_as_array( chr, 0, chromSizes[chr] )
		#fill in blanks
		invalid = np.isnan( values )
		values[ invalid ] = defaultVal[i];
		#log transform
		if doLog[i]>=0:
			values =np.log10(values + doLog[i])
		#smooth data
		if smoothings[i]!=0:
			#reflect = at edge of array, the data will be mirrored
			#truncate = don't incorporate data more than X away - probably a great speed increase since the distrib goes out to infinity
			values = gaussian_filter(values, smoothings[i], mode='reflect', truncate=4.0)
		#sample data
		values = values [useThese[chr]]
		#append data;
		allData[curTot:(curTot+len(values)),i] = values;
		curTot = curTot + len(values);


dataMedians = np.median(allData,axis=0);

#n_components = number of components to keep
#init= how to initialize
#sparseness {data|components|None]: where to enforce sparsity
#beta: degree of sparseness (larger is more sparse, default=1
#eta: degree of correctness to maintain; smaller=more error; default=0.1
#tol: tolerance in stopping
#max_iter: number of iterations
#random_state: seed
myDecomp = decomposition.NMF(sparseness='components', n_components=args.components);
myDecomp.fit(allData);

allInBWs = [];
for i in range(0,len(IDs)):
	#input GB tracks
	curBW = BigWigFile(open(files[i]))
	allInBWs.append(curBW);


allOutBWs = []
for i in range(0, len(myDecomp.components_)):
	toBW = subprocess.Popen(["wigToBigWig","stdin",args.chrsFile,"%s_comp%i_pOpen.bw"%(args.outFPre,i)], stdin=subprocess.PIPE)
	toBW.stdin.write("track type=wiggle_0\n")
	allOutBWs.append(toBW);


for chr in chrOrder:
	chrData  = np.empty([chromSizes[chr],len(IDs)]);
	for i in range(0,len(IDs)):
		values = allInBWs[i].get_as_array( chr, 0, chromSizes[chr] )
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
		chrData[:,i] = values;
	transed = myDecomp.transform(chrData);
	for i in range(0,len(myDecomp.components_)):
		allOutBWs[i].stdin.write("fixedStep chrom=%s start=1 step=1\n"%(chr))
		allOutBWs[i].stdin.write("\n".join(map(str,transed[:,i])));
		allOutBWs[i].stdin.write("\n");


for i in range(0,len(myDecomp.components_)):
	temp = allOutBWs[i].communicate()
	if temp[0] is not None:
		sys.stderr.write("component %i : %s"%(i,temp[0]));
	if temp[1] is not None:
		sys.stderr.write("component %i : %s"%(i,temp[1]));

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
if (args.logFP is not None):
	logFile.close();
