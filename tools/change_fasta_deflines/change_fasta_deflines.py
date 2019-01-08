#!/usr/bin/env python
import sys, argparse
'''Accepts either csv (default) or tab-delimited files with old/new sequence names, creating a dictionary of 
respective key:value pairs. Parses an input fasta file for 'old' names, replacing them with 'new' names, writing
renamed sequences to a fasta file. NOTE: use of tab-delim text file for renaming requires '-t' on cmd line.'''
#USAGE EXAMPLE 1: python change_fasta_def_lines.py csv_rename_file.csv fasta_2_rename.fasta renamedSequences.fasta
#USAGE EXAMPLE 2: python change_fasta_def_lines.py tab_delim_rename_file.txt -t fasta_2_rename.fasta renamedSequences.fasta

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Sept 2017'''

#parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument ("-t", "--tab_delim", help = "name fasta definition lines from tab-delim file", action = "store_true")
parser.add_argument("inFileHandle") #csv file with current fasta file names in column 1 and desired names in col 2
parser.add_argument("inFileHandle2") #fasta file containing sequences requiring name replacement
parser.add_argument("outFileHandle") #user-specified output filename
args = parser.parse_args()

#open a writable output file that will be over-written if it already exists
outfile= open(args.outFileHandle,'w')
dict = {} #dictionary to hold old_name:new_name key:value pairs
splitter = ',' #default char to split lines at
#determine if input naming file is csv (default) or tab delim text
if args.tab_delim:
    splitter = '\t' #change splitter to tab if comd line args contain '-t'

#create dictionary using key/value pairs from csv file of old/new names
with open(args.inFileHandle,'r') as inputFile:
#read in each line and split at comma into key:value pairs
    for line in inputFile:
        #remove whitespace from end of lines, split at comma, assigning to key:value pairs
        line2 = line.rstrip()
        splitLine = line2.split(splitter)
        old_name = splitLine[0]
        new_name = splitLine[1]
        dict[old_name] = new_name

#parse fasta deflines for 'old' names and, if found, replace with new names
with open(args.inFileHandle2,'r') as inputFile2:
    for line in inputFile2:
        #find the definition lines, remove trailing whitespace & '>'
        if ">" in line:
            originalDefline = line.rstrip().replace(">","",1)
            #check for a match to any of the dict key
            if dict.has_key(originalDefline):
                #find the index of that item in the list
                newDefline= dict[originalDefline]
                #print("the new name"), newDefline
                # print each item to make sure the right name is being entered
                outfile.write(">" + newDefline + "\n")
            else:
                #write out the original defline sequence name
                print ("Defline not in dictionary: "), originalDefline
                outfile.write(">" + originalDefline + "\n")
        else:
        #in lines without ">", write out sequence as it was
            seq = line.rstrip()
            outfile.write(seq+"\n")

inputFile.close()
inputFile2.close()
outfile.close()
