#!/usr/bin/env python

'''Accepts fasta files containing amino acid sequence, reading them in as
amino acid sequence objects.  Reads influenza clade defintions (i.e. amino 
acids at certain positions) from .csv file into dictionary structure. Searches
each of the amino acid sequence objects for a list of matching clades, assigns
the most 'evolved' (i.e. child as opposed to parent clade) to the sequence. Appends
"_cladename" to the Sequence name and generates a fasta file of original sequences with 
modified names.'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, Oct 2017'''

import sys,string,os, time, Bio
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
inFileHandle1 = sys.argv[1] #batch fasta file with sequences to be parsed
inFileHandle2 = sys.argv[2] # .csv file containing clade definitions and "depth"
outFileHandle = sys.argv[3] #user-specified name for output file of aa seq's with clade suffixes
outFile= open(outFileHandle,'w') #open a writable, appendable output file
seqList = [] #list of aa sequence objects to parse for clade definitions
cladeList = [] #empty list to hold clade tuples i.e. ("3C.3a", 1 ,{"3":"I", "9":"V"..})

'''Searches record for required amino acids at defined positions. If found, assigns
clade name to sequence name by appending underscore and clade name to record id.'''
def call_clade(record):
    print "---------------------------------------------------------------------"
    print "Parsing %s for matching flu clade definitions..." % (record.id)
    matchList = [] #empty list to hold clades that match 100%
    #iterate over each tuple in the clade list
    for clade in cladeList:
        cladeName = clade[0] #temp variable for name
        depth = clade[1] #temp variable for depth
        sites = clade[2] #temp variable for aa def dictionary
        shouldFind = len(sites) #number of sites that should match
        found = 0 #a counter to hold matches to antigenic sites
        #iterate over each position in sites dictionary
        for pos, aa in sites.iteritems():
            #translate pos to corresponding index in target sequence
            index = int(pos) - 1
            #if record at index has same amino acid as 'aa', increment 'found'
            if record[index] == aa:
                found += 1
        if (found == shouldFind):
            #add the matching clade tuple to the list of matches
            matchList.append(clade)
    return matchList

'''Compares depth level of clades in a list and returns the most granular one'''
def decide_clade_by_depth(matchList):
    #empty variable for maximum depth encountered
    max_depth = 0
    best_match_name = '' #variable to hold most granular clade
    #for each matching clade, check depth of the corresponding tuple
    for clade in matchList:
        #if the current clade is 'deeper' than the one before it
        if clade[1] > max_depth:
            #store this depth
            max_depth = clade[1]
            #store name of the clade
            best_match_name = clade[0]
    return best_match_name

'''opens the .csv file of clade definitions and clade "depth" '''
with open (inFileHandle2, 'r') as clade_file:
    #remove whitespace from the end of each line and split elements at commas
    for line in clade_file:
        #print "Current Line in File:" + line
        sites={} #initialize a dictionary for clade
        elementList = line.rstrip().split(',')
        new_list = [] #start a new list to put non-empty strings into
        #remove empty stings in list
        for item in elementList:
            if item != '':
                new_list.append(item)
        name = new_list.pop(0) #move 1st element to name field
        depth = int(new_list.pop(0)) #move 2nd element to depth field
        #read remaining pairs of non-null elements into clade def dictionary
        for i in range(0, len(new_list), 2):
            #move next 2 items from the list into the dict
            pos = new_list[i]
            aa = new_list[i + 1]
            sites[pos] = aa
        #add the clade info as a tuple to the cladeList[]
        oneClade =(name, depth, sites)
        cladeList.append(oneClade)
    print "The List of Clades:"
    for clade in cladeList:
        print "Clade Name: %s Depth: %i Antigenic Sites: %i" % (clade[0], clade[1], len(clade[2]))
        for pos, aa in clade[2].iteritems():
            print "Pos: %s\tAA: %s" % (pos,aa)

'''opens readable input file of sequences to parse using filename from cmd line,
    instantiates as AA Sequence objects, with ppercase sequences'''
with open(inFileHandle1,'r') as inFile:
    #read in Sequences from fasta file, uppercase and add to seqList
    for record in SeqIO.parse(inFile, "fasta", alphabet=IUPAC.protein):
        record = record.upper()
        seqList.append(record) #add Seq to list of Sequences
    print "\n%i flu HA sequences will be compared to current clade definitions..." % len(seqList)
    #parse each target sequence object
    for record in seqList:
        clade_call = '' #empty variale for final clade call on sequence
        matchingCladeList = call_clade(record) #holds matching clade tuples
        #if there is more than one clade match
        if len(matchingCladeList) > 1:
            #choose the most granular clade based on depth
            clade_call = decide_clade_by_depth(matchingCladeList)
        #if there is only one clade call
        elif len(matchingCladeList) > 0:
            clade = matchingCladeList[0] #take the first tuple in the list
            clade_call = clade[0] #clade name is the first item in the tuple
        #empty list return, no matches
        else:
            clade_call = "No_Match"
        print clade_call
        seq_name = record.id
        mod_name = seq_name + "_" + clade_call
        print "New Sequence Name: " + mod_name
        record.id = mod_name


#output fasta file with clade calls appended to sequence names
SeqIO.write(seqList,outFile,"fasta")

#print("\n%i Sequences Extracted to Output file: %s"  % ((len(extractedSeqList),outFileHandle)))
inFile.close()
clade_file.close()
outFile.close()

