#!/usr/bin/env python
import sys, re

'''Accepts a sequence-collapsed fasta output from USEARCH (drive5) software and reformats the fasta 
definition lines by replacing occurences of ';size=N;' with '_xN' and writing output to fasta 
(N = number of identical sequences represented by the collapsed sequence). If N is not greater
than 1 (i.e. only 1 sample with that sequence), replaces ';size=N;' with ''. For example, 
'>sequence_A;size=2;' is replaced with '>sequence_A_x2', whereas '>sequence_B;size=1;' is
replaced with 'sequence_B'.
#USAGE EXAMPLE: python reformat_usearch_collapsed_fasta.py usearch_collapsed_sequences.fasta output.fasta

Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Feb 2018'''

inFileHandle = sys.argv[1] #input fasta filename
outFileHandle = sys.argv[2] #output fasta filename
outFile = open(outFileHandle,'w') #open a writable output file

separator = "_x" #the string separating sequence name from number of sequences, N
regex = re.compile(";size=[0-9]{0,};") #regex snippet from debuggex

#parse fasta definition lines for pattern matching regex
with open(inFileHandle,'r') as inFile:
    for line in inFile:
        if ">" in line:
            #look for regex pattern in fasta definition line
            matchArray = regex.findall(line)
            if len(matchArray) > 0: #replace the matching substring
                substringToReplace = matchArray[0]
                endIndex = len(substringToReplace) 
                digits = substringToReplace[6:endIndex -1] #digits between ';size=' and ';'
                if int(digits) > 1: #show number of sequences if greater than 1
                    replacementString = separator + digits
                else:
                    replacementString = "" #otherwise, just display sequence name
                newDefline = line.rstrip().replace(substringToReplace, replacementString)
                outFile.write(newDefline + "\n")
        else: #in lines without ">", write out sequence unmodified
            seq = line.rstrip()
            outFile.write(seq+"\n")

inFile.close()
outFile.close()
