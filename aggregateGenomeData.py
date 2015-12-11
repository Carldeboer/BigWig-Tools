#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python2.7
import argparse


parser = argparse.ArgumentParser(description='This program retrieves the values of genome tracks that correspond to specific loci and outputs a value that is an aggregation function of the scores in that region.')
parser.add_argument('-c',dest='lociFile',	metavar='<inGFF>',help='Input GFF file of loci to scan', required=True);
parser.add_argument('-b',dest='inBED',	action='count', help='Input loci file is in BED format', required=False,default=0);
parser.add_argument('-i',dest='inFile',	metavar='<inFile1>',help='Genomic track to collect data from', required=True);
parser.add_argument('-i2',dest='inFile2',	metavar='<inFile2>',help='Genomic track to collect data from for reverse strand', required=False);
parser.add_argument('-t',dest='format',	metavar='<format>',help='inFile format: one of {BEDGR,WIG,BIGWIG,BIGBED}', required=False,default="BW");
parser.add_argument('-f',dest='flank',	metavar='<flank>',help='Number of flanking bases', required=False,default=0);
parser.add_argument('-r',dest='rev',	action='count', help='Scan reverse complements of loci (still output in + orientation)', required=False,default=0);
parser.add_argument('-e',dest='exclusive', action='count', help='Coordinates of loci are exclusive of ends', required=False, default=0);
parser.add_argument('-k',dest='keep', action='count', help='Keep loci on absent chromosomes', required=False, default=0);
parser.add_argument('-d',dest='default', metavar='<default>', help='Use this as the default value for missing data', required=False, default="NaN");
parser.add_argument('-a',dest='agFtn', metavar='<aggregateFtn>', help='Use this as the aggregation function. One of {min,max,mean,sum}', required=False, default="mean");
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results', required=True);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-cc',dest='correctChr', action='count',help='add/remove "chr" to/from chromosome if chr absent?', required=False, default=0);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();

from numpy import *
 
from bx.intervals.io import GenomicIntervalReader
import warnings
import MYMATH
import MYUTILS
import GENOMEDATA
import sys

if args.inBED>0:
	scanThese = GENOMEDATA.BED(args.lociFile);
else: #GFF
	scanThese = GENOMEDATA.GFF(args.lociFile);

inclusive=1;
if args.exclusive>0:
	inclusive=0;


if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

outFile = MYUTILS.smartGZOpen(args.outFP,'w');

#reverse loci if req
if args.rev>0:
	scanThese.flipStrands_();

if scanThese.length()==0:
	outFile.write("");
	outFile.close();
	quit();
#calculate average and max length of GFF entries
scanThese.coord_to_i_();
lengthMax = 0;
lengthSum = 0;
lengthCount = 0;
numStrandless=0;
for i in range(0, scanThese.length()):
	if scanThese[i][GENOMEDATA.EN]<scanThese[i][GENOMEDATA.ST]:
		sys.stderr.write("Locus %s has end less than start - swapping\n"%(scanThese[i][GENOMEDATA.NAME]));
		tmp = scanThese[i][GENOMEDATA.EN];
		scanThese[i][GENOMEDATA.EN]= scanThese[i][GENOMEDATA.ST];
		scanThese[i][GENOMEDATA.ST] = tmp;
	if scanThese[i][GENOMEDATA.ST]<0:
		raise Exception("Locus %s start before 0"%(scanThese[i][GENOMEDATA.NAME]));
	curLen = scanThese[i][GENOMEDATA.EN]-scanThese[i][GENOMEDATA.ST]+inclusive;
	lengthMax = max([lengthMax,curLen]);
	lengthCount+=1;
	lengthSum+=curLen;
	if scanThese[i][GENOMEDATA.STR]!="+" and scanThese[i][GENOMEDATA.STR]!="-":
		numStrandless+=1
avgLength = int(round(lengthSum/lengthCount))

if numStrandless>0 and args.inFile2 is not None:
	raise Exception("Error: loci contain strandless entries, but genome tracks provided for each strand!!");
if numStrandless>0:
	sys.stderr.write("Warning: Strandless loci detected; assuming forward orientation.\n");
padding = int(args.flank);

#read in track file(s)
if args.format=="BIGWIG" or args.format=="BW":
	from bx.bbi.bigwig_file import BigWigFile
	inFile1 = BigWigFile(open(args.inFile))
	if args.inFile2 is not None:
		inFile2 = BigWigFile(open(args.inFile2))
elif args.format=="BIGBED" or args.format=="BB":
	from bx.bbi.bigwig_file import BigBedFile
	inFile1 = BigBedFile(open(args.inFile))
	if args.inFile2 is not None:
		inFile2 = BigBedFile(open(args.inFile2))
elif args.format=="WIG" or args.format=="W":
	from bx.arrays.wiggle import WiggleReader
	inFile1 = WiggleReader(open(args.inFile))
	if args.inFile2 is not None:
		inFile2 = WiggleReader(open(args.inFile2))
elif args.format=="BEDGR" or args.format=="BG":
	from bx.arrays.bed import BedReader
	inFile1 = BedReader(open(args.inFile))
	if args.inFile2 is not None:
		inFile2 = BedReader(open(args.inFile2))
else:	
	raise Exception("Unrecognized format!");



for locus in scanThese:
	if args.verbose>0:
		print("Scanning %s"%(locus[GENOMEDATA.NAME]))
	stF = max(locus[GENOMEDATA.ST] - padding,0);
	enF = locus[GENOMEDATA.EN] + padding + inclusive;
	try:
		if (locus[GENOMEDATA.STR]=="-" and args.inFile2 is not None):
			if args.format=="BIGWIG" or args.format=="BW" or args.format=="BIGBED" or args.format=="BB":
				values = inFile2.get_as_array( locus[GENOMEDATA.CHR], stF, enF )
		else:
			if args.format=="BIGWIG" or args.format=="BW" or args.format=="BIGBED" or args.format=="BB":
				values = inFile1.get_as_array( locus[GENOMEDATA.CHR], stF, enF )
	except OverflowError as e:
		sys.stderr.write("OverflowError at '%s'; st=%d, en=%d\n"%(locus[GENOMEDATA.NAME],locus[GENOMEDATA.ST],locus[GENOMEDATA.EN]));
		raise(e);
	if values is None and args.correctChr>0:
		#try again adding chr or taking it away
		if locus[GENOMEDATA.CHR][:3]=="chr":
			 locus[GENOMEDATA.CHR] =  locus[GENOMEDATA.CHR][3:]
		else:
			 locus[GENOMEDATA.CHR]="chr"+ locus[GENOMEDATA.CHR];
		try:
			if (locus[GENOMEDATA.STR]=="-" and args.inFile2 is not None):
				if args.format=="BIGWIG" or args.format=="BW" or args.format=="BIGBED" or args.format=="BB":
					values = inFile2.get_as_array( locus[GENOMEDATA.CHR], stF, enF )
			else:
				if args.format=="BIGWIG" or args.format=="BW" or args.format=="BIGBED" or args.format=="BB":
					values = inFile1.get_as_array( locus[GENOMEDATA.CHR], stF, enF )
		except OverflowError as e:
			sys.stderr.write("OverflowError at '%s'; st=%d, en=%d\n"%(locus[GENOMEDATA.NAME],locus[GENOMEDATA.ST],locus[GENOMEDATA.EN]));
			raise(e);
	if values is None and args.keep==0:
		print("Skipping %s because there is no data for it"%(locus[GENOMEDATA.NAME]))
		continue;
	if values is None:
		curAgg="NaN"
	else:
		if stF>locus[GENOMEDATA.ST] - padding: #add nans if off chr start
			values = concatenate([ [ NAN ] * ( stF-(locus[GENOMEDATA.ST] - padding) ), values ]);
		#remove nans
		invalid = isnan( values )
		if args.default=="NaN":
			del values[ invalid ];#remove them
		else:
			values[ invalid ] = float(args.default);
		if len(values)==0:
			curAgg="NaN"
		else:
			if args.agFtn=="mean" or args.agFtn=="sum":
				curAgg = sum(values);
				if args.agFtn=="mean":
					curAgg = curAgg/len(values);
			elif args.agFtn=="max":
				curAgg = max(values)
			elif args.agFtn=="min":
				curAgg = min(values)
			else: 
				raise Exception("Unrecognized aggrgate function!");
				
	outFile.write("%s\t%i\t%i\t%s\t%g\n"%(locus[GENOMEDATA.CHR], locus[GENOMEDATA.ST], locus[GENOMEDATA.EN], locus[GENOMEDATA.NAME], curAgg));
 
#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));


outFile.close();
if (args.logFP is not None):
	logFile.close();
