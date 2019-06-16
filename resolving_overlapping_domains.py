#************  THE PROGRAM TO RESOLVE OVERLAPPING BETWEEN DOMAINS IN PROTEINS  *************
#           Author: Elzbieta Janowska-Sedja, 16/06/2019
#		        As input file hmmer file with domain signatures is used. 
#		        The format of the file is described below    
#*******************************************************************************************   
# The general rules for solving overlapping problem were adopted from previous study (Seidl et al., 2011)
#*******************************************************************************************

#!/usr/bin/python
import sys

#*********WORKING DIRECTORY*******************************************************************
workDir1 = "your_working_directory"

#************** INPUT FILE *******************************************************************
file1 = "%s/your_input_file_with_hmmer_output.txt" % (workDir1)
#
#  hmmer output is a tab-delimited file wih a following format including 9 columns:
#  protein_id	protein_description	pfam_domain_id	pfam_domain_description	start_domain_coordinates	end_domain_coordinates	score	e-value	-log10(e-value)
#**********************************************************************************************
#************** OUTPUT FILES ******************************************************************
file2 = "%s/output_file_with_overlapping_domains_to_be_resolved_manualy_based_on_the_sore_in_matrix_included.txt" % (workDir1)
#file3 includes all proteins even the one where manual solving of overlapping is required. Thus, once resolve the manual overlapping, the domains for protein not resolved automaticaly (proteins from output file2) need to be updated in file3
file3 = "%s/output_for_all_domains_per_protein.txt" % (workDir1)
#**********************************************************************************************
repDom = dict()
fh1 = open (file1, "r")
fh2 = open (file2, "a")
fh3 = open (file3, "w")
t = tuple()
set1=set()
set2= set()
set3 = set()
#********************************************************************************************
class domainInstance:
	def __init__(self, did, descr, start, end, evalue):
		self.domainId = did
		self.domainDescr = descr
		self.start = float(start)
		self.end = float(end)
		self.length = float(end) - float(start)
		self.evalue = float(evalue)
		self.score = float(score)
	def __str__(self):
		return self.domainId
#********************************************************************************************	
def findOverlappingGroups(group, value):
	ovearlappingGroups = list()
	for dom, rules in group.items():
		print dom.domainId, rules
		fh2.write("\n%s\t%s\n" %(dom.domainId, rules))
		if 0 in rules:
			processed = False
			for ovearlappingGroup in ovearlappingGroups:
				if ovearlappingGroup.__contains__(dom):
					processed = True
					break
			if processed == True:
				continue
			ovearlappingGroup = depthFirstSeardch(dom, group, value, list())
			ovearlappingGroups.append(ovearlappingGroup)
	return ovearlappingGroups
#*********************************************************************************************			
def depthFirstSeardch(dom, group, value, result):
	if result.__contains__(dom):
		return result
	result.append(dom)
	rules = group[dom]
	for i in range (0, len(rules)):
		if rules[i] != -1:
			depthFirstSeardch(value[i], group, value, result)
	return result
#*********************************************************************************************	
def start_compare(a, b):
	if a.start == b.start :
		return 0
	elif a.start < b.start:
		return -1
	return 1

#**********************************************************************************************
def resolveOverlap(group, value):

	remove = []
	removeZ = []
	
	for dom, rules in group.items():
		if 0 not in rules:
			for z in range(0, len(rules)):
				if rules[z] == 1:
					remove.append(value[z])
					removeZ.append(z)
					

			
	removeZ = list(set(removeZ))			
	removeZ.sort()		
	removeZ.reverse()
	for r in remove:
		try:
			del group[r]
		except KeyError:
			pass
	nextPass = False
	for z in removeZ:
		value.pop(z)
	for dom, rules in group.items():
		for z in removeZ:
			rules.pop(z)
		if 0 in rules:
			nextPass = True
	retain = [x for x in value if x not in remove]
	
	if nextPass == True:
		if len(remove) == 0:
			print " found a loop",key
			fh2.write("\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
			fh2.write("\nfound a loop %s\n" %key)
			fh2.write("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n")
			for d in retain:
				print d.domainId,
				fh2.write("%s\t" %d.domainId)
			fh2.write("              \n")
			print " "
			ovearlappingGroups = findOverlappingGroups(group, value)
			for ovearlappingGroup in ovearlappingGroups:
				print key, " ",
				#fh2.write("\n*****\n")
				#fh2.write("%s\t\n" %key)
				#for d in ovearlappingGroup:
					#print d.domainId,
					#fh2.write("\t%s" %d.domainId) 
				print ""
				for i in range (1, len(ovearlappingGroup)):
					retain.remove(ovearlappingGroup[i])
 
			return retain
			
		else:
			retain = resolveOverlap(group, retain)
			#fh3.write("\n+++++%s\n" %group)
	return retain
	
		         
def applyRules(l,k):
	if (l.end <= k.end):
		if l.end < k.start:
			return -1
		else:
			overlap = l.end - k.start
                        if overlap/l.length < 0.1 and overlap/k.length < 0.1:
				return -1
	if (k.end <= l.end):
		if k.end < l.start:

			return -1
		else:
			overlap = k.end - l.start
			#print overlap, overlap/l.length, overlap/k.length
                        if overlap/l.length < 0.1 and overlap/k.length < 0.1:
				return -1
	
#evalue rule
	if abs(l.evalue - k.evalue) > 5: 
		if l.evalue > k.evalue:
			return 1
		else:
			return 0
#lenght rule				 
        if l.length != k.length:
		if l.length > k.length:
			return 1
		else:
			return 0
#score rule
	
	if l.length == k.length:
		if l.score > k.score:
			return 1
		else:
			return 0
	return 1

#************************************************************************************************	
#******************************MAIN PROGRAM******************************************************
#************************************************************************************************
for line in fh1.readlines():
    line= line.split('\t')
    fgId = line[0]#first column in the input file is FG id
    pfamAcc = line[2]#third column in the input file is pfam accession number
    pfamDescr = line[3]
    start = int(line[4])#fith column in the input file is a start coordinate of the domain
    end = int(line[5])#sixth column in the input file is a end coordinate of the domain
    if start > end:
	temp = start
	start = end
	end = temp
    eValue= float(line[8].replace('\n',''))#column 9 in the input file is -log10(e-value)
    score = float(line[6])#column 7 in the input file is a score
    t = domainInstance(pfamAcc,pfamDescr, start, end, eValue) # t is an instance of domainInstance class
    repDom.setdefault(fgId, list())#assigning key and value to dictionary repDom where value is a list
    repDom[fgId].append(t)#filling the value(list)with list of domains
    
c = 0
for key, value in repDom.items():
   
    if len(value) >= 1:
	group = dict()
	for i in range (0, len(value)):
		rules = list()
		group.setdefault(value[i], rules)
		for j in range (0, len(value)):
			if value[i] == value [j]:
				rules.append(-1)
			else:
				rules.append(applyRules(value[i], value[j]))
		
	retain = resolveOverlap(group, value)	
	result = list()
	retain = sorted(retain, cmp=start_compare)
	for d in retain:
		result.append(d.domainId)
		#c = c+1
		#print c, key, result
	
	#fh3.write("%s\t%s\t%s\n" %(c, key, result))
		fh3.write("%s\t%s\t%s\t%s\t%s\n" %(key, d.domainId, d.domainDescr, d.start, d.end))
fh1.close()  
fh2.close()
fh3.close() 
#*****************************************************************************************
	

