#!/usr/bin/env python
'''Reads in a fasta file of extracted antigenic sites and one containing a 
reference flu antigenic map, reading them in as protein SeqRecords. Compares each amino
acid of each sample antigenic map to corresponding sites in the reference and replaces
identical amino acids with dots. Writes headers (including amino acid position numbers
read in from the respective index array), the reference amino acid sequence and column
headings required for non-aggregated line lists. Outputs headers and modified (i.e. dotted)
sequences to a csv file.'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, Nov 2017'''

import sys,string,os, time, Bio, re, argparse
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

inputAntigenicMaps = sys.argv[1] #batch fasta file with antigenic map sequences
refAntigenicMap = sys.argv[2] #fasta file of reference antigenic map sequence
antigenicSiteIndexArray = sys.argv[3] #antigenic site index array csv file
cladeDefinitionFile = sys.argv[4] #clade definition csv file
outFileHandle = sys.argv[5] #user-specifed output filename

lineListFile = open(outFileHandle,'w') #open a writable output file

indicesLine = "" #comma-separated antigenic site positions
cladeList = [] #list of clade names read from clade definition file
ref_seq = "" #reference antigenic map (protein sequence)
seqList = [] #list of aa sequences to compare to reference

BC_list = [] #empty list for BC samples
AB_list = [] #empty list for AB samples
ON_list = [] #empty list for ON samples
QC_list = [] #empty list for QC samples
nonprov_list = [] #empty list for samples not in above 4 provinces
#dictionary for location-separated sequence lists
segregated_lists = {'1_BC':BC_list,'2_AB':AB_list,'3_ON':ON_list,'4_QC': QC_list, '5_nonprov': nonprov_list}
uniqueSeqs = {} #empty dict with unique seqs as keys and lists of SeqRecords as values

def replace_matching_aa_with_dot(record):
    """Compare amino acids in record to reference sequence, replace matching symbols
    with dots, and return record with modified amino acid sequence."""
    orig_seq = str(record.seq) #get sequence string from SeqRecord
    mod_seq = ""
    #replace only those aa's matching the reference with dots
    for i in range(0, len(orig_seq)):
        if (orig_seq[i] == ref_seq[i]):
            mod_seq = mod_seq  + '.'
        else:
            mod_seq  = mod_seq  + orig_seq[i]
    #assign modified sequence to new SeqRecord and return it
    rec = SeqRecord(Seq(mod_seq,IUPAC.protein), id = record.id, name = "", description = "")
    return rec

def extract_clade(record):
    """Extract clade name (or 'No_Match') from sequence name and return as clade name. """
    if record.id.endswith('No_Match'):
        clade_name = 'No_Match'
        end_index = record.id.index(clade_name)
        record.id = record.id[:end_index -1]
        return clade_name
    else: #
        for clade in cladeList:
            if record.id.endswith(clade):
                clade_name = clade
                end_index = record.id.index(clade)
                record.id = record.id[:end_index -1]
                return clade_name
    
def sort_by_location(record):
    """Search sequence name for province name or 2 letter province code and add SeqRecord to
    province-specific dictionary."""
    seq_name = record.id
    if ('-BC-' in seq_name) or ('/British_Columbia/' in seq_name):
        BC_list.append(record) #add Sequence record to BC_list
    elif ('-AB-' in seq_name) or ('/Alberta/' in seq_name):
        AB_list.append(record) #add Sequence record to AB_list
    elif ('-ON-' in seq_name) or ('/Ontario/' in seq_name):
        ON_list.append(record) #add Sequence record to ON_list
    elif ('-QC-' in seq_name) or ('/Quebec/' in seq_name):
        QC_list.append(record) #add Sequence record to QC_list
    else:
        nonprov_list.append(record) #add Sequence record to nonprov_list
    return

def extract_province(record):
    """Search sequence name for province name or 2 letter province code and return province."""
    seq_name = record.id
    if ('-BC-' in seq_name) or ('/British_Columbia/' in seq_name):
        province = 'British Columbia'
    elif ('-AB-' in seq_name) or ('Alberta' in seq_name):
        province = '/Alberta/'
    elif ('-ON-' in seq_name) or ('/Ontario/' in seq_name):
        province = 'Ontario'
    elif ('-QC-' in seq_name) or ('/Quebec/' in seq_name):
        province = 'Quebec'
    else:
        province = "other"
    return province

def get_sequence_length(record):
    """Return the length of a sequence in a Sequence record."""
    sequenceLength = len(str((record.seq)))
    return sequenceLength

def get_antigenic_site_substitutions(record):
    """Count number of non-dotted amino acids in SeqRecord sequence and return as substitutions."""
    sequenceLength = get_sequence_length(record)
    seqString = str(record.seq)
    matches = seqString.count('.')
    substitutions = sequenceLength - matches
    return substitutions

def calculate_percent_id(record, substitutions):
    """Calculate sequence identity to a reference (based on substitutions and sequence length) and return percent id."""
    sequenceLength = get_sequence_length(record)
    percentID = (1.00 - (float(substitutions)/float(sequenceLength)))
    return percentID

def output_linelist(sequenceList):
    """Output a list of SeqRecords to a non-aggregated line list in csv format."""
    for record in sequenceList:
        #get province, clade from sequence record
        province = extract_province(record)
        clade = extract_clade(record)
        #calculate number of substitutions and % id to reference
        substitutions = get_antigenic_site_substitutions(record)
        percentID = calculate_percent_id(record,substitutions)
        name_part = (record.id).rstrip() + ','
        clade_part = clade + ','
        substitutions_part = str(substitutions) + ','
        percID_part = str(percentID) + ','
        col = " ," #empty column
        sequence = str(record.seq).strip()
        csv_seq = ",".join(sequence) +","
        #write linelisted antigenic maps to csv file
        comma_sep_line = name_part + col + clade_part + col + csv_seq + substitutions_part + percID_part + "\n"
        lineListFile.write(comma_sep_line)
    return
	
with open (antigenicSiteIndexArray,'r') as siteIndices:
    """Read amino acid positions from antigenic site index array and print as header after one empty row."""
    col = "," #empty column
    #read items separated by comma's to position list
    for line in siteIndices:
        #remove whitespace from the end of each line
        indicesLine = line.rstrip()
    row1 = "\n"
    #add comma-separated AA positions to header line
    row2 = col + col + col + col + indicesLine + "\n"
    #write first (empty) and 2nd (amino acid position) lines to linelist output file
    lineListFile.write(row1)
    lineListFile.write(row2)

with open (refAntigenicMap,'r') as refMapFile:
    """Read reference antigenic map from fasta and output amino acids, followed by column headers."""
    #read sequences from fasta to SeqRecord, uppercase, and store sequence string to ref_seq
    record = SeqIO.read(refMapFile,"fasta",alphabet=IUPAC.protein)
    record = record.upper()
    ref_seq = str(record.seq).strip() #store sequence in variable for comparison to sample seqs
    col = "," #empty column
    name_part = (record.id).rstrip() + ','
    sequence = str(record.seq).strip()
    csv_seq = ",".join(sequence)
    #output row with reference sequence displayed above sample sequences
    row3 = name_part + col + col + col + csv_seq + "\n"
    lineListFile.write(row3)
    #replaces digits in the indicesLine with empty strings
    positions = indicesLine.split(',')
    numPos = len(positions)
    empty_indicesLine = ',' * numPos
    #print column headers for sample sequences
    row4 = "Sequence Name,N,Clade,Extra Substitutions," + empty_indicesLine + "Number of Amino Acid Substitutions in Antigenic Sites,% Identity of Antigenic Site Residues\n"
    lineListFile.write(row4)
    print ("\nREFERENCE ANTIGENIC MAP: '%s' (%i amino acids)" % (record.id, len(record)))

with open(cladeDefinitionFile,'r') as cladeFile:
    """Read clade definition file and store clade names in a list."""
    #remove whitespace from the end of each line and split elements at commas
    for line in cladeFile:
        elementList = line.rstrip().split(',')
        name = elementList[0] #move 1st element to name field
        cladeList.append(name)

with open(inputAntigenicMaps,'r') as extrAntigMapFile:
    """Read antigenic maps as protein SeqRecords and add to list."""
    #read Sequences from fasta file, uppercase and add to seqList
    for record in SeqIO.parse(extrAntigMapFile, "fasta", alphabet=IUPAC.protein):
        record = record.upper()
        seqList.append(record) #add Seq to list of Sequences

#print number of sequences to be process as user check
print "\nCOMPARING %i flu antigenic map sequences to the reference..." % len(seqList)
#parse each antigenic map sequence object
for record in seqList:
    #assign Sequence to dictionaries according to location in name
    sort_by_location(record)
#sort dictionary keys that access province-segregated lists
sorted_segregated_list_keys = sorted(segregated_lists.keys())
print "\nSequence Lists Sorted by Province: "
#process each province-segregated SeqRecord list
for listname in sorted_segregated_list_keys:
    #acesss list of sequences by the listname key
    a_list = segregated_lists[listname]
    # sort original SeqRecords by record id (i.e. name)
    a_list = [f for f in sorted(a_list, key = lambda x : x.id)]
    mod_list = [] # empty temporary list
    for record in a_list:
        #replace matching amino acid symbols with dots
        rec = replace_matching_aa_with_dot(record)
        mod_list.append(rec) #populate a list of modified records
    segregated_lists[listname] =  mod_list
    print "\n'%s' List (Amino Acids identical to Reference Masked): " % (listname)
    #output the list to csv as non-aggregated linelist
    output_linelist(segregated_lists[listname])

extrAntigMapFile.close()
refMapFile.close()
lineListFile.close()