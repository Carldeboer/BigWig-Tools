import MYUTILS

CHR=0
ST = 1
EN = 2
NAME=3
STR=4
SCORE=5
RF = 6
SOURCE=7
TYPE=8


GFF_CHR=0;
GFF_SOURCE=1;
GFF_TYPE=2;
GFF_ST=3;
GFF_EN=4;
GFF_SCORE=5;
GFF_STR=6;
GFF_RF=7;
GFF_NAME=8;

BED_CHR=0;
BED_ST=1;
BED_EN=2;
BED_NAME=3;
BED_SCORE=4;
BED_STR=5;

def leg(a,b):
	if a>b:
		return 1;
	if b<a:
		return -1;
	return 0;

def locusGTLT_noStrand(a,b): #-1= a<b, 1= a>b, 0= a=b
	if a[CHR]==b[CHR]:
		return leg(a[ST],b[ST]);
	else:
		return leg(a[CHR],b[CHR]);
	

def locusGTLT(a,b): #-1= a<b, 1= a>b, 0= a=b
	if a[CHR]==b[CHR]:
		if a[STR]==b[STR]:
			return leg(a[ST],b[ST]);
		else:
			return leg(a[STR],b[STR]);
	return leg(a[CHR],b[CHR]);
	

def locusOverlap(a,b): #-1= a<b, 1= a>b, 0= a overlaps b
	if a[CHR]==b[CHR]:
		if a[STR]==b[STR]:
			return leg(a[ST],b[ST])
			if (a[ST]>b[EN]):
				return 1;
			elif (a[EN]<b[ST]):
				return -1;
			else:
				return 0;
		else:
			return leg(a[STR],b[STR]);
	return leg(a[CHR],b[CHR])

class LOCI:
	
	def readGFF(self,gffFP): #readGFF
		self.allData = [];
		for line in MYUTILS.smartGZForeach(gffFP):
			line = line.rstrip();
			if line is None  or line=="" or line[0]=="#" or line[0]==">" or line[0:11]=="track name=":
			  continue;
			data = line.split("\t");
			data2 = [data[GFF_CHR], data[GFF_ST], data[GFF_EN], data[GFF_NAME], data[GFF_STR], data[GFF_SCORE], data[GFF_RF], data[GFF_SOURCE], data[GFF_TYPE]];
			self.allData.append(data2);
	
	def readBED(self,bedFP): #readBED
		self.allData = [];
		for line in MYUTILS.smartGZForeach(bedFP):
			line = line.rstrip();
			if line is None  or line=="" or line[0]=="#" or line[0]==">" or line[0:11]=="track name=":
			  continue;
			data = line.split("\t");
			data2 = [data[BED_CHR], data[BED_ST], data[BED_EN]];
			if len(data)>=(BED_NAME+1):
				data2.append(data[BED_NAME]);
			else:
				data2.append("");
			if len(data)>=(BED_STR+1):
				data2.append(data[BED_STR]);
			else:
				data2.append("");
			if len(data)>=(BED_SCORE+1):
				data2.append(data[BED_SCORE]);
			else:
				data2.append("");
			if len(data)>=(BED_SCORE+2):
				data2.append(data[6:]); # the rest
			self.allData.append(data2);
	
	def writeBED(self,bedFP):
		outFile = smartGZOpen(bedFP, "w");
		for e in self.allData:
			e2=[e[CHR], e[ST], e[EN], e[NAME], e[SCORE], e[STR]];
			if length(e)>6:
				e2 = e2 + e[6:];
			outFile.write("\t".join(e2)+"\n");
		outFile.close();

	def writeGFF(self,gffFP,source="SOURCE",type="TYPE",rf="0"):
		outFile = smartGZOpen(gffFP, "w");
		for e in self.allData:
			e2=[e[CHR], source, type, e[ST], e[EN], e[SCORE], e[STR], rf, e[NAME]];
			if TYPE in e:
				e2[GFF_TYPE] = e[TYPE]
				e2[GFF_SOURCE] = e[SOURCE]
				e2[GFF_RF] = e[RF]
			elif SOURCE in e:
				e2[GFF_SOURCE] = e[SOURCE]
				e2[GFF_RF] = e[RF]
			elif RF in e:
				e2[GFF_RF] = e[RF]
			outFile.write("\t".join(e2)+"\n");
		outFile.close();
		
	def __getitem__(self,i):
		return self.allData[i];
	
	def __setitem__(self,i, vals):
		self.allData[i]=vals;
	
	def __add__(self,other):
		return GFF.new(self.allData.dup()+other.allData);
	
	def length(self):
		return len(self.allData);
	
	def append(self,newOne):
		self.allData.append(newOne);
	
	def coord_to_i_(self):
		for i in range(0,self.length()):
			self.allData[i][ST]=int(self.allData[i][ST]);
			self.allData[i][EN]=int(self.allData[i][EN]);
	
	def dup(self):
		theDup = GFF.new(None);
		theDup.allData=copy.deepcopy(self.allData);
		return theDup;
	
#	def assemblePieces():
#		uniqueHash = Hash.new();
#		for entry in self.allData:
#			if uniqueHash[entry[NAME]] is None:
#				uniqueHash[entry[NAME]]=[];
#			
#			uniqueHash[entry[NAME]].append(entry);
#		newData=[];
#		for name in uniqueHash:
#			newEntry = uniqueHash[name][0];
#			newEntry[ST]=newEntry[ST].to_i();
#			newEntry[EN]=newEntry[EN].to_i();
#	
#			1.upto(uniqueHash[name].length()-1){|i|
#				newEntry[ST] = [newEntry[ST],uniqueHash[name][i][ST].to_i()].min();
#				newEntry[EN] = [newEntry[EN],uniqueHash[name][i][EN].to_i()].max();
#			}
#			newData.append(newEntry);
#		self.allData=newData;
	
	
	def flipStrands_(self): #MirrorStrands
		for i in range(0,self.length()):
			if self.allData[i][STR]=="-":
				self.allData[i][STR]="+";
			elif self.allData[i][STR]=="+": 
				self.allData[i][STR]="-";
	
	def mirrorStrands_(self): #MirrorStrands
		newData=[];
		for i in range(0,self.length()):
			newEntry=self.allData[i].dup;
			newEntry[STR]="+";
			newData.append(newEntry);
			newEntry=self.allData[i].dup;
			newEntry[STR]="-";
			newData.append(newEntry);
		self.allData=newData;
	
	
	def splitStrandless(self): #SplitStrandless
		temp = self.dup();
		temp.splitStrandless_();
		return temp;
	
	def splitStrandless_(self): #SplitStrandless
		newData=[];
		for i in range(0,self.length()):
			if self.allData[i][STR]==".":
				newEntry=self.allData[i].dup;
				newEntry[STR]="+";
				newData.append(newEntry);
				newEntry=self.allData[i].dup;
				newEntry[STR]="-";
				newData.append(newEntry);
			else:
				newData.append(self.allData[i].dup);
		self.allData=newData;
	
	
	def sort_noStrand_(self):
		self.coord_to_i_();
		sort(self.allData,cmp=locusGTLT_noStrand)
	
	def sort_(self):
		self.coord_to_i_();
		sort(self.allData,cmp=locusGTLT)
	
	def removeRedundant(self):
		temp = self.dup();
		temp.removeRedundant_();
		return temp;
	
	def removeRedundant_(self):
		self.sort_();
		for i in range(self.length()-2,-1,-1): 
			#check if the one before it is the same
			if locusOverlap(self.allData[i],self.allData[i+1])==0:
				self.allData.delete_at(i);
			
	
	def splitBiStranded_(self):
		newAllData=[];
		for entry in self.allData:
			if entry[STR]==".":
				tempN = entry;
				tempP = entry.dup();
				tempN[STR]="-";
				tempP[STR]="+";
				newAllData.append(tempP);
				newAllData.append(tempN);
			else:
				newAllData.append(entry);
			
		self.allData=newAllData;
	
	def dropOverlappingSelf(self):
		temp = self.dup();
		temp.dropOverlappingSelf_();
		return temp;
	
	def dropOverlappingSelf_(self):
		self.sort_();
		dropThese=[false]*self.length();
		i=1;
		while (i<self.length()):
			#p([i,u]);
			overlap=locusOverlap(self.allData[i],self.allData[i-1]);
			if overlap==0:
				dropThese[i]=true;
				dropThese[i-1]=true;
			
			i+=1;
		
		for i in range(self.length()-1,-1,-1):
			if dropThese[i]:
				del self.allData[i] 
	
	def dropOverlapping(self,other):
		temp = self.dup();
		temp.dropOverlapping_(other);
		return temp;
	
	def dropOverlapping_(self,otherGFF):
		otherGFF = otherGFF.dup();
		otherGFF.splitBiStranded_();
		otherGFF.sort_();
		self.splitBiStranded_();
		self.sort_();
		i=0;
		u=0;
		newData=[];
		while (i<self.length()  and  u<otherGFF.length()):
			#p([i,u]);
			overlap=locusOverlap(self.allData[i],otherGFF.allData[u]);
			if overlap>0:
				u+=1;
			elif overlap<0:
				newData.append(self.allData[i].dup());
				i+=1;
			else:
				#the entries overlap
				i+=1;
			
		
		while (i<self.length()):
			newData.append(self.allData[i].dup());
			i+=1;
		
		self.allData=newData;
	
	def notIn(self,other):
		temp = self.dup();
		temp.notIn_(other);
		return temp;
	
	def notIn_(self,otherGFF):
		otherGFF = otherGFF.dup();
		otherGFF.splitBiStranded_();
		otherGFF.sort_();
		self.splitBiStranded_();
		self.sort_();
		i=0;
		u=0;
		newData=[];
		while (i<self.length()  and  u<otherGFF.length()):
			#p([i,u]);
			overlap=locusOverlap(self.allData[i],otherGFF.allData[u]);
			if overlap>0: #other happens before mine
				u+=1;
			elif overlap<0: #mine happens before other; keep it
				newEntry = self.allData[i].dup();
				newData.append(newEntry);
				i+=1;
			else: #they overlap
				i+=1;
			
		
		while (i<self.length()):
				newEntry = self.allData[i].dup();
				newData.append(newEntry);
				i+=1;
		
		self.allData=newData;
	
	def ands(self,other):
		temp = self.dup();
		temp.and_(other);
		return temp;
	
	def and_(self,otherGFF):
		otherGFF = otherGFF.dup();
		otherGFF.splitBiStranded_();
		otherGFF.sort_();
		self.splitBiStranded_();
		self.sort_();
		i=0;
		u=0;
		newData=[];
		while (i<self.length()  and  u<otherGFF.length()):
			#p([i,u]);
			overlap=locusOverlap(self.allData[i],otherGFF.allData[u]);
			if overlap>0:
				u+=1;
			elif overlap<0:
				i+=1;
			else:
				#the entries overlap
				#find overlap
				newEntry = self.allData[i].dup();
				newEntry[ST] = [self.allData[i][ST],otherGFF.allData[u][ST]].max();
				newEntry[EN] = [self.allData[i][EN],otherGFF.allData[u][EN]].min();
				newEntry[NAME]+=";"+otherGFF.allData[u][NAME];
				#p(newEntry);
				newData.append(newEntry);
				if self.allData[i][EN]>otherGFF.allData[u][EN]:
					u+=1;
				else:
					i+=1;
		self.allData=newData;
	
	def ors(self):
		temp = self.dup();
		temp.or_();
		return temp;
	
	def or_(self): #Or
		self.sort_();
		newEntries = [];
		newEntries.append(self.allData[0]);
		for i in range(1,self.length()):
			last=newEntries.length()-1;
			if self.allData[i][CHR]!=newEntries[last][CHR]  or  self.allData[i][STR]!=newEntries[last][STR]: #new Chr or strand
				newEntries.append(self.allData[i]);
			elif self.allData[i][ST]<=newEntries[last][EN]: #next one starts before end of last
				newEntries[last][EN] = [newEntries[last][EN],self.allData[i][EN]].max();
				newEntries[last][NAME]=newEntries[last][NAME]+";"+self.allData[i][NAME];
			else: # non overlapping
				newEntries.append(self.allData[i]);
		self.allData=newEntries;
	
	
	def nots(self,chrSizes):
		temp = self.dup();
		temp.not_(chrSizes);
		return temp;
	
	def not_(self,chrSizes): #Not
		self.splitStrandless_();
		self.or_();
		self.sort_();

		strands =['-','+'];
		allDataHash=Hash.new();
		for s in strands:
			allDataHash[s]=Hash.new();
			for k in chrSizes:
				allDataHash[s][k]=[];
		for i in range(0,allData.length()):
			if chrSizes[self.allData[i][CHR]] is None:
				print("Skipping #%i because %s not in chrom sizes...\n"%(i,self.allData[i][CHR]));
			elif self.allData[i] is None:
				raise Exception("Why is #%i None?"%(i));
			else:
				allDataHash[self.allData[i][STR]][self.allData[i][CHR]].append(self.allData[i]);
		newEntries = [];
	
		for chr in allDataHash["+"]:
			for strand in allDataHash:
				last=1;
				for i in range(0, allDataHash[strand][chr].length()):
					newEntry = allDataHash[strand][chr][i].dup;
					newEntry[ST]=last;
					newEntry[EN]=allDataHash[strand][chr][i][ST]-1;
					if (newEntry[ST]<=newEntry[EN]):
						newEntries.append(newEntry);
					last=allDataHash[strand][chr][i][EN]+1;
				newEntry = [];
				newEntry[CHR]=chr;
				newEntry[STR]=strand;
				newEntry[ST]=last;
				newEntry[EN]=chrSizes[chr];
				newEntry[NAME]='full';
				if (newEntry[ST]<=newEntry[EN]):
					newEntries.append(newEntry);
				
					
		self.allData=newEntries;
	
	def mergeNonOverlapping(self,other):
		allData1=self.allData;
		allData2=other.allData;
		self.allData=[];
		
		sort(allData1,cmp=locusGTLT_noStrand)
		sort(allData2,cmp=locusGTLT_noStrand)
		
		for e in allData1:
			self.allData.append(e);
		i=0;
		j=0;
		while (i<allData1.length()  and  j<allData2.length()):
			if allData1[i][CHR]>allData2[j][CHR]:
				self.allData.append(allData2[j]);
				j+=1;
			elif allData1[i][CHR]<allData2[j][CHR]:
				i+=1;
			elif allData1[i][ST]>allData2[j][EN]: # 2 ends before 1 starts
				self.allData.append(allData2[j]);
				j+=1
			elif allData1[i][EN]<allData2[j][ST]: # 1 ends before 2 starts
				i+=1;
			else: #they overlap, so skip 2
				j+=1;
			
		
		# if j runs out first, all j are added.
		# if i runs out first
		

class BED (LOCI):
	def __init__(self,bedFile=None):
		if bedFile is None:
			self.allData=[];
		elif isinstance(bedFile,str):
			self.readBED(bedFile);
		else:
			self.allData=bedFile;

class GFF (LOCI):
	def __init__(self,gffFile=None):
		if gffFile is None:
			self.allData=[];
		elif isinstance(gffFile,str):
			self.readGFF(gffFile);
		else:
			self.allData=gffFile;
	
	

class ChrSizes:
	def __init__(self,file):
		self.chrs={};
		inFile = MYUTILS.smartGZOpen(file,"r")
		for l in inFile:
			l = l.rstrip();
			c,len = l.split("\t");
			self.chrs[c]=int(len);
	def __iter__(self):
		return iter(sorted(self.chrs.keys()));
	def __getitem__(self,c):
		return self.chrs[c];



#
#def operate(*data):
#	raise Exception("No data provided!") if data.length()==0;
#	dataNew= makeNewAllData(data[0],nil);
#	data[0].keys.each(){|key|
#		dataNew[key]=[nil]*data[0][key].length();
#		0.upto(data[0][key].length()-1){|i|
#			if data.length==1:
#				dataNew[key][i]=yield(data[0][key][i]);
#			else:
#				curData=[];
#				0.upto(data.length()-1){|j|
#					curData.append(data[j][key][i]);
#				}
#				dataNew[key][i]=yield(curData);
#			
#		}
#	}
#	return dataNew;
#
#
#def operate_(*data):
#	raise Exception("No data provided!") if data.length()==0;
#	data[0].keys.each(){|key|
#		data[0][key]=[nil]*data[0][key].length();
#		0.upto(data[0][key].length()-1){|i|
#			if data.length==1:
#				data[0][key][i]=yield(data[0][key][i]);
#			else:
#				curData=[];
#				0.upto(data.length()-1){|j|
#					curData.append(data[j][key][i]);
#				}
#				data[0][key][i]=yield(curData);
#			
#		}
#	}
#
#
#
#def readGPRFile(allData,gprFP,parseString, parseFeats, idPos, subFromCoords):
#	#allData = Hash.new();
#	allCounts = Hash.new();#to divide by at the end for overlapping probes
#	allData.keys.each(){|chr|
#		allCounts[chr]=[0]*allData[chr].length();
#	}
#	skip=2;
#	theRegex=Regexp.new(/#{parseString}/);
#	chrPos=nil;
#	stPos=nil;
#	enPos=nil;
#	0.upto(parseFeats.length()-1){|i|
#		if parseFeats[i]=="C"#:chr
#			chrPos=i+1;#+1 because the first match of a regex is the whole string
#		elif parseFeats[i]=="S":#start
#			stPos=i+1;
#		elif parseFeats[i]=="E":
#			enPos=i+1;
#		
#	}
#	if enPos is None:
#		enPos=stPos;
#	
#	#counter=0;
#	for line in MYUTILS.smartGZForeach(gprFP):
#		#print(line);
#		line = e.rstrip();
#		if skip>0:
#			skip-=1;
#			continue
#		
#		continue if line is None  or  line=="" or line[0]=="\"";
#		stuff = line.split("\t");
#		curData = stuff[26].to_f();#26=ratio of medians, 27=ratio of means, 28=median of ratios, 29=mean of ratios
#		rawName = stuff[idPos];
#		#p(rawName)
#		#theRegex.match(rawName){|match| #executes only if match
#		match=/#{parseString}/.match(rawName);
#		if match:
#			#print("string matched ");
#			st=match[stPos].to_i()-subFromCoords;
#			en=match[enPos].to_i()-subFromCoords;
#			chr=match[chrPos];
#			continue if allData[chr] is None;
#			#allData[chr]=[] if allData[chr] is None;
#			#allCounts[chr]=[] if allCounts[chr] is None;
#
#			st.upto(en){|i|
#				if allData[chr][i] is None:
#					allData[chr][i]=curData;
#					allCounts[chr][i]=1;
#				else:
#					allData[chr][i]+=curData;
#					allCounts[chr][i]+=1;
#				
#			}
#		
#		#counter+=1;
#		#return if counter>100;
#	}
#	allData=operate(allData,allCounts){|sum,count| 
#		if (sum and count) is None:
#			nil;
#		else:
#			sum/count;
#		end};
#	return allData;
#
#
#def readBedGraphFile(allData,bedgrFP):
#	skipChr = nil;
#	MYUTILS.smartGZForeach(bedgrFP){|line|
#		line = e.rstrip();
#		if line is None  or  line==""  or  line[0:19]=="track type="  or  line[0]=="#:"
#			continue
#		
#		chr,st,en,val = line.split("\t");
#		continue if skipChr==chr;
#		if allData[chr] is None:
#			print("skipping #{chr}\n");
#			skipChr=chr;
#			continue
#		
#		st.to_i().upto(en.to_i()){|i|
#			allData[chr][i]=val.to_f();
#		}
#	}
#	return allData;
#
#		
#def readMapFile(allData,inFP):
#	curData=ChrVals.new(inFP);
#	if !allData is None:
#		allData.keys.each(){|k|
#			if !curData[k] is None:
#				allData[k]=curData[k];
#			
#		}
#	else:
#		allData=curData;
#	
#	return allData;
#
#
#def shiftDataBy(allData,by, defVal):
#	allData2=Hash.new()
#	allData.keys.each{|chr|
#		allData2[chr] = [defVal]*allData[chr].length();
#		if by<0:
#			allData2[chr][0,allData[chr].length()+by] = allData[chr][-by,allData[chr].length()+by];
#		else:
#			allData2[chr][by,allData[chr].length()-by] = allData[chr][0,allData[chr].length()-by];
#		
#	}
#	return allData2;
#
#
#class BigArray:
#	def initialize(size):
#		self.mySize=0;
#		self.full=size;
#		self.array = Array.new(size);
#	
#	def trim():
#		self.array=self.array[0,self.mySize];
#	
#	def [](ind,len):
#		#return nil if ind+len>=self.mySize;
#		return self.array[ind,len];
#	
#	def [](ind):
#		return nil if ind>=self.mySize;
#		return self.array[ind];
#	
#	def []=(i,len,vals):
#		raise Exception("Assignment array not same length #{vals.length().to_s()} as assignment length #{len.to_s}") if vals.length!=len;
#		0.upto(len-1){|j|
#			self[i+j]=vals[j];
#		}
#	
#	def []=(i,val):
#		self.mySize=[i+1,self.mySize].max();
#		if i>self.full:
#			self.full=self.full*2;
#			self.array[self.full]=nil;
#		
#		self.array[i]=val;
#	
#	def to_array():
#		return self.array[0,self.mySize];
#	
#	def from_array(a):
#		self.array=a;
#		self.mySize=a.length();
#		self.full=a.length();
#	
#
#
#class GenomeTrack:
#	def initialize():
#		self.data = Hash.new();
#	
#	def []=(i,v):
#		self.data[i]=BigArray.new(1000);
#		if !v.empty?:
#			self.data[i].from_array(v);
#		
#	
#	def [](i):
#		self.data[i];
#	
#	def keys():
#		return self.data.keys;
#	
#	def has?(chr):
#		return !self.data[chr] is None;
#	
#	def add(chr):
#		self.data[chr]=BigArray.new(100000);
#	
#	def to_hash():
#		self.data.keys.each{|k|
#			self.data[k]=self.data[k].to_array();
#		}
#		return self.data;
#	
#	def from_hash(h):
#		h.keys.each{|k|
#			self.data[k] = BigArray.new(h[k].length());
#			self.data[k][0,h[k].length()]=h[k];
#		}
#	
#
#
#def readWigFile(wigFP):
#	return readWigFile(nil,wigFP);
#
#def readWigFile(allData,wigFP):
#	readAll=false;
#	if allData is None :
#		readAll=true;
#		allData = GenomeTrack.new();
#	
#	curChr=nil;
#	curI=nil;
#	vs=nil;
#	step=1;
#	span=1;
#	skip=false;
#	isBG = false;
#	MYUTILS.smartGZForeach(wigFP){|line|
#		#print(line);
#		continue if line is None  or  line=="";
#		if line[0]>=?A:
#			if line[0:19]=="track type=wiggle_0":
#				#print("Converting #{line}\n");
#				continue
#			elif line[0:19]=="track type=bedGraph":
#				isBG=true;
#				break;
#			elif line=~/^track type=.*wiggle_0/:
#				continue
#			
#			if line=~/^fixedStep.*\s+chrom=([._a-zA-Z0-9]+)/:
#				curChr = $1; 
#				if allData[curChr] is None:
#					if readAll:
#						allData[curChr]=[];
#						skip=false;
#					else:
#						#print("skipping #{curChr}\n");
#						skip=true;
#						continue
#					
#				else:
#					skip=false;
#				
#				if line=~/start=([0-9]+)/:
#					curI=$1.to_i()-1;
#				else:
#					curI=0;
#				
#				if line=~/step=([0-9]+)/:
#					step=$1.to_i();
#				else:
#					step=1;
#				
#				if line=~/span=([0-9]+)/:
#					span=$1.to_i();
#				else:
#					span=1;
#				
#	
#				#print("Doing FS #{curChr}...\n");
#				vs=false;
#				continue
#			
#			if line=~/^variableStep.*\s+chrom=([a-zA-Z0-9]+)/:
#				curChr = $1; 
#				if line=~/span=([0-9]+)/:
#					span=$1.to_i();
#				else:
#					span=1;
#				
#				if allData[curChr] is None:
#					if readAll:
#						allData[curChr]=[];
#						skip=false;
#					else:
#						#print("skipping #{curChr}\n");
#						skip=true;
#						continue
#					
#				else:
#					skip=false;
#				
#				#print("Doing VS #{curChr}...\n");
#				vs=true;
#				continue
#			
#		
#		continue if skip;
#		if vs:
#			stuff = line.split("\t");
#			curI=stuff[0].to_i()-1;
#			if stuff[1]=="NaN\n" or stuff[1]=="nil\n"  or  stuff[1]=="NA\n":
#				val=nil;
#			else:
#				val=stuff[1].to_f();
#			
#			
#		else:
#			if line=="NaN\n" or line=="nil\n"  or  line=="NA\n":
#				val=nil;
#			else:
#				val=line.to_f();
#			
#		
#		#print(curChr+"\n") if allData[curChr] is None;
#		
#		0.upto(span-1){|i|
#			allData[curChr][curI+i]=val;
#		}
#		curI+=step;
#	}
#	if isBG:
#		return readBedGraphFile(allData,wigFP);
#	
#	if allData.kind_of?(GenomeTrack):
#		return allData.to_hash();
#	
#	return allData;
#
#
#class ChrSepWig(ChrSeq):
#  def loadKey(key):
#    temp = readWigFile(nil,self.fileHash[key]);
#    self.seqHash[key] = temp[temp.keys()[0]];
#
#
#
#def readBedFile(allData, bedFP):
#	#allData = Hash.new();
#	MYUTILS.smartGZForeach(bedFP){|line|
#		line = e.rstrip();
#		continue if line is None  or  line=="";
#		chr, st, en, val = line.split("\t");
#		continue if allData[chr] is None;
#		#if allData[chr] is None
#		#	allData[chr]=[];
#		#
#		en=en.to_i()-1;
#		st=st.to_i()-1;#-1 for index from 0 change
#		val = val.to_f();
#		st.upto(en){|i|
#			allData[chr][i]=val;
#		}
#	}
#	return allData;
#
#
#def readNDBedFile(allData,bedFP):#goddamn morons... Stands for "no data bed" because they don't have data in it, the reads are mapped, but not added up
#	#allData = Hash.new();
#	MYUTILS.smartGZForeach(bedFP){|line|
#		line = e.rstrip();
#		continue if line is None  or  line=="";
#		chr, st, en = line.split("\t");
#		continue if allData[chr] is None;
#		#if allData[chr] is None
#		#	allData[chr]=[];
#		#
#		en=en.to_i()-1;
#		st=st.to_i()-1;#-1 for index from 0 change
#		st.upto(en){|i|
#			allData[chr][i]=0.0 if allData[chr][i] is None;
#			allData[chr][i]+=1;
#		}
#	}
#	return allData;
#
#
#def readGFFFile(allData, inFP):
#	return readCustomFile(allData, inFP, 3, 4, 5, 0, 1, 1, 0, "#", "\t");
#
#
#def makeNewAllData(chrMap, defaultValue):
#	allData=Hash.new();
#	chrMap.keys.each(){|chr|
#		allData[chr]=[defaultValue]*chrMap[chr].length();
#	}
#	return allData;
#
#
#def readCustomFile(allData, inFP, startPos, endPos, valPos, chrPos, subS, subE, skipFirst, skipPre, delim):
#	#allData = Hash.new();
#	MYUTILS.smartGZForeach(inFP){|line|
#		continue if skipPre!=nil  and  line[0:skipPre.length()]==skipPre;
#		if skipFirst>0:
#			skipFirst-=1;
#			continue
#		
#		line = e.rstrip();
#		continue if line is None  or  line=="";
#		data = line.split(delim);
#		continue if allData[data[chrPos]] is None;
#		st = data[startPos].to_i()-subS;
#		en= data[endPos].to_i()-subE;
#		data[valPos]=data[valPos].to_f();
#		#if allData[data[chrPos]] is None
#		#	allData[data[chrPos]]=[];
#		#
#		st.upto(en){|i|
#			allData[data[chrPos]][i]=data[valPos];
#		}
#	}
#	return allData;
#
#
#def smoothAll(allData, chrMap, defaultValue, numDP):
#	allDataNew = Hash.new();
#	chrMap.keys.each{|chr|
#		allDataNew[chr] = [defaultValue]*chrMap[chr].length();
#		runningSum = 0.0;
#		numEntries = 0;
#		posQ=[];
#		0.upto(numDP-1){|i| #don't subtract anything, since we're building it up.
#			if allData[chr][i] is None:
#				runningSum=runningSum+defaultValue;
#				posQ.append(defaultValue);
#			else:
#				runningSum=runningSum+allData[chr][i];
#				posQ.append(allData[chr][i]);
#			
#			numEntries+=1;
#			if i-((numDP-1)/2)>=0:
#				allDataNew[chr][i-((numDP-1)/2)]=runningSum/numEntries;
#			
#				
#		}
#		(numDP).upto(chrMap[chr].length()-1){|i|
#			if allData[chr][i] is None:
#				runningSum=runningSum+defaultValue;
#				posQ.append(defaultValue);
#			else:
#				runningSum=runningSum+allData[chr][i];
#				posQ.append(allData[chr][i]);
#			
#			runningSum=runningSum-posQ.shift();
#			allDataNew[chr][i-((numDP-1)/2)]=runningSum/numEntries;
#		}
#		#don't add anything since we're running off the
#		(chrMap[chr].length()-((numDP-1)/2)).upto(chrMap[chr].length()-1){|i|
#			runningSum=runningSum-posQ.shift();
#			numEntries-=1;
#			allDataNew[chr][i]=runningSum/numEntries;
#		}
#	}
#	return allDataNew;
#
#
#def fillBlankData(allData, chrMap, defaultValue):
#	allDataNew = Hash.new();
#	chrMap.keys.each{|chr|
#		allDataNew[chr] = [defaultValue]*chrMap[chr].length();
#		lastVal=defaultValue;
#		lastPos=0;
#		if allData[chr] is None:
#			continue
#		
#		(0).upto(chrMap[chr].length()-1){|i|
#			if !allData[chr][i] is None:
#				allDataNew[chr][i]=allData[chr][i];
#				#fill blanks till here
#				step = (allData[chr][i]-lastVal) / (0.0+i-lastPos);
#				total=lastVal;
#				(lastPos+1).upto(i-1){|j|
#					total+=step;
#					allDataNew[chr][j]=total;
#				}
#				#update vars
#				lastVal=allData[chr][i];
#				lastPos=i;
#			
#		}
#		step = (defaultValue-lastVal) / (0.0+chrMap[chr].length()-1-lastPos);
#		total=lastVal;
#		(lastPos+1).upto(chrMap[chr].length()-1){|j|
#			total+=step;
#			allDataNew[chr][j]=total;
#		}
#		
#	}
#	return allDataNew;
#
#
#def readElandFile(allData, inFP):
#	#allData=Hash.new();
#	chrPos=6;
#	startPos=7;
#	orientationPos=8;
#	oligoPos=1;
#	File.foreach(inFP){|line|
#		line = e.rstrip();
#		continue if line is None   or  line=="";
#		stuff=line.split("\t");
#		continue if stuff[chrPos]==""  or  allData[stuff[chrPos]] is None;
#		#allData[stuff[chrPos]]=[] if allData[stuff[chrPos]] is None;
#		thePos = stuff[startPos].to_i()+(stuff[oligoPos].length()/2)-1;
#		if allData[stuff[chrPos]][thePos] is None:
#			allData[stuff[chrPos]][thePos]=0;
#		
#		allData[stuff[chrPos]][thePos]+=1;
#	}
#	return allData;
#
#
#def readElandFileSS(allDataF,allDataR,inFP):
#	#allDataF=Hash.new();
#	#allDataR=Hash.new();
#	chrPos=6;
#	startPos=7;
#	orientationPos=8;
#	oligoPos=1;
#	File.foreach(inFP){|line|
#		line = e.rstrip();
#		continue if line is None   or  line=="";
#		stuff=line.split("\t");
#		continue if stuff[chrPos]==""  or  allData[stuff[chrPos]] is None;
#		#allDataF[stuff[chrPos]]=[] if allDataF[stuff[chrPos]] is None;
#		#allDataR[stuff[chrPos]]=[] if allDataR[stuff[chrPos]] is None;
#		thePos = stuff[startPos].to_i()+(stuff[oligoPos].length()/2)-1;
#		if stuff[orientationPos]=="F":
#			if allDataF[stuff[chrPos]][thePos] is None:
#				allDataF[stuff[chrPos]][thePos]=0;
#			
#			allDataF[stuff[chrPos]][thePos]+=1;
#		else:
#			if allDataR[stuff[chrPos]][thePos] is None:
#				allDataR[stuff[chrPos]][thePos]=0;
#			
#			allDataR[stuff[chrPos]][thePos]+=1;
#		
#	}
#	return [allDataF, allDataR];
#
#
#
#def GDtoS(allData):
#	theStr="";
#	allData.keys.each(){|chr|
#		0.upto(allData[chr].length()-1){|i|
#			theStr+=chr+"\t"+(i+1).to_s()+"\t"+allData[chr][i].to_s()+"\n";
#		}
#	}
#	return theStr;
#
#def toSGR(allData,outFP):
#	outFile = File.open(outFP, "w");
#	allData.keys.each(){|chr|
#		0.upto(allData[chr].length()-1){|i|
#			outFile.write([chr,(i+1).to_s(),allData[chr][i]].join("\t")+"\n");
#		}
#	}
#
#def toGZWig(allData,outFP):
#	outFile = outGZ(outFP,"w");
#	toWigStream(allData,outFile);
#	outFile.close();
#
#
#def toWig(allData,outFP):
#	outFile = smartGZOut(outFP, "w");
#	toWigStream(allData,outFile);
#	outFile.close();
#
#
#def toWigStream(allData,outFile):
#	outFile.write("track type=wiggle_0\n");
#	allData.keys.each(){|chr|
#		outFile.write("fixedStep chrom=#{chr} start=1 step=1\n");
#		outFile.write(allData[chr].join("\n")+"\n");
#	}
#
#
