#!/usr/bin/env python
'''Reads in a fasta file of antigenic maps and one with the reference antigenic map as 
protein SeqRecords. Compares amino acids of sample antigenic maps to corresponding sites
in the reference and masks identical amino acids with dots. Writes headers (including
amino acid position numbers read from the respective index array), the reference amino
acid sequence and column headings required for both non-aggregated and aggregated line lists. 
Outputs all headers and modified (i.e. dotted) sample sequences to a csv file.'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, Jan 2018'''

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

agg_lineListFile = open(outFileHandle,'w') #open a writable output file

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
prov_lists = {'1_BC':BC_list,'2_AB':AB_list,'3_ON':ON_list,'4_QC': QC_list, '5_nonprov': nonprov_list}

def replace_matching_aa_with_dot(record):
    """Compare amino acids in record to reference, mask identical symbols with dots, and return modified record."""
    orig_seq = str(record.seq) #sequence string from SeqRecord
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
    else: #
        for clade in cladeList:
            if record.id.endswith(clade):
                clade_name = clade
    return clade_name
    
def extract_sample_name(record, clade):
    """Extract sample name from sequence name and return sample name. """
    end_index = record.id.index(clade)
    sample_name = record.id[:end_index -1]
    #return sample name as sequence name minus underscore and clade name
    return sample_name

def sort_by_location(record):
    """Search sequence name for province name or 2-letter province code and add SeqRecord to
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
    """Search sequence name for province name or 2-letter province code and return province."""
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
    """Return length of sequence in a SeqRecord."""
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
    """Calculate percent sequence identity to reference sequence, based on substitutions
and sequence length and return percent id as a ratio (i.e. 0.90 no 90%)."""
    sequenceLength = get_sequence_length(record)
    percentID = (1.00 - (float(substitutions)/float(sequenceLength)))
    return percentID

def output_aggregated_linelist(a_list):
    """Output aggregated line list of SeqRecords in csv format."""
    sequevars = {} #dict of sequevar: SeqRecord list
    firstRecordID = None
    #examine dotted/masked sequences in list and assign unique ones as dict keys
    for rec in a_list:
        rec = replace_matching_aa_with_dot(rec)
        sequence =str(rec.seq)
        #if the sequence is a key in the dict, add SeqRecord to list
        if sequence in sequevars:
            #if sequence already in dict as a key, increment the value
            sequevars[sequence].append(rec)
        else:
            #if sequence not in dict, add is as new key with list of 1 SeqRecord
            sequevars[sequence] = [rec]
    #get list of sorted unique sequence keys
    sorted_unique_seq_keys = sorted(sequevars.keys())
    #process each list of SeqRecords sharing a unique sequence
    for u in sorted_unique_seq_keys:
        #access list of sequences by unique sequence
        listOfSeqs = sequevars[u]
        #sort this list of SeqRecords by record.id (i.e. name)
        listOfSeqs = [f for f in sorted(listOfSeqs, key = lambda x : x.id)]
        N = len(listOfSeqs)
        #output details of first SeqRecord to csv
        firstRecord = listOfSeqs[0]
        province = extract_province(firstRecord)
        clade = extract_clade(firstRecord)
        substitutions = get_antigenic_site_substitutions(firstRecord)
        percentID = calculate_percent_id(firstRecord,substitutions)
        name = extract_sample_name(firstRecord, clade)
        name_part = name.rstrip() + ','
        N_part = str(N) + ','
        clade_part = clade + ','
        substitutions_part = str(substitutions) + ','
        percID_part = str(percentID) + ','
        col = " ," #empty column
        sequence = str(firstRecord.seq).strip()
        csv_seq = ",".join(sequence) +","
        comma_sep_output = name_part + N_part + clade_part + col + csv_seq + substitutions_part + percID_part + "\n"
        #write first member of unique sequence list to csv
        agg_lineListFile.write(comma_sep_output)
        #print sequence records in sequevar to console
        print "\n\t\t%i SeqRecords matching Sequevar: %s" % (len(listOfSeqs), u)

        #to uncollapse sequevar group, print each member of the sequevar list to csv output
        '''for i in range(1,len(listOfSeqs)):
            currentRec = listOfSeqs[i]
            province = extract_province(currentRec)
            clade = extract_clade(currentRec)
            substitutions = get_antigenic_site_substitutions(currentRec)
            percentID = calculate_percent_id(currentRec,substitutions)
            name_part = (currentRec.id).rstrip() + ','
            N_part = "n/a" + ','
            clade_part = clade + ','
            substitutions_part = str(substitutions) + ','
            percID_part = str(percentID) + ','
            col = " ," #empty column
            sequence = str(currentRec.seq).strip()
            csv_seq = ",".join(sequence) +","
            comma_sep_output = name_part + N_part + clade_part + col + csv_seq + substitutions_part + percID_part + "\n"
            agg_lineListFile.write(comma_sep_output)   '''
    return
	
with open (antigenicSiteIndexArray,'r') as siteIndices:
    """Read amino acid positions from antigenic site index array and print as header after one empty row."""
    col = "," #empty column
    #read amino acid positions and remove trailing whitespace
    for line in siteIndices:
        #remove whitespace from the end of each line
        indicesLine = line.rstrip()
    row1 = "\n"
    #add comma-separated AA positions to header line
    row2 = col + col + col + col + indicesLine + "\n"
    #write first (empty) and 2nd (amino acid position) lines to output file
    agg_lineListFile.write(row1)
    agg_lineListFile.write(row2)

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
    agg_lineListFile.write(row3)
    positions = indicesLine.split(',')
    numPos = len(positions)
    empty_indicesLine = ',' * numPos
    #print column headers for sample sequences
    row4 = "Sequence Name,N,Clade,Extra Substitutions," + empty_indicesLine + "Number of Amino Acid Substitutions in Antigenic Sites,% Identity of Antigenic Site Residues\n"
    agg_lineListFile.write(row4)
    print ("\nREFERENCE ANTIGENIC MAP: '%s' (%i amino acids)" % (record.id, len(record)))

with open(cladeDefinitionFile,'r') as cladeFile:
    """Read clade definition file and store clade names in list."""
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

#print number of sequences to be processed as user check
print "\nCOMPARING %i flu antigenic map sequences to the reference..." % len(seqList)
for record in seqList:
    #assign SeqRecords to province-specific dictionaries
    sort_by_location(record)

#access prov segregated lists in order
sorted_prov_keys = sorted(prov_lists.keys())
print "\nSequence Lists Sorted by Province: "
for prov in sorted_prov_keys:
    current_list = prov_lists[prov]
    #mask AA's identical to reference sequence with dot
    masked_list = [] # empty temporary list to park masked sequences
    for record in current_list:
        masked_rec = replace_matching_aa_with_dot(record)
        masked_list.append(masked_rec)
    prov_lists[prov] =  masked_list #replace original SeqRecord list with masked list

#group sequences in province-sorted list into clades
for prov in sorted_prov_keys:
    prov_list = prov_lists[prov]
    by_clades_dict = {} #empty dict for clade:seqRecord list groups
    print "\n'%s' List (Amino Acids identical to Reference are Masked): " % (prov)
    for rec in prov_list:
        clade = extract_clade(rec)
        if clade in by_clades_dict:
            #if clade already in dict as key, append record to list (value)
            by_clades_dict[clade].append(rec)
        else: #add clade as key to dict, value is list of 1 SeqRecord
            by_clades_dict[clade] = [rec]
    #get list of alphabetically sorted clade keys
    sorted_clade_keys = sorted(by_clades_dict.keys())
    print "\tNumber of clades: ", len(by_clades_dict)
    #group each list of sequences in clade by sequevars
    for key in sorted_clade_keys:
        print "\n\tCLADE: %s Number of Members: %i" % (key, len(by_clades_dict[key]))
        a_list = by_clades_dict[key]
        for seqrec in a_list:
            print "\t %s: %s" %(seqrec.id,str(seqrec.seq))
        #output the list to csv as aggregated linelist
        output_aggregated_linelist(a_list)
    
print("Aggregated Linelist written to file: '%s\n'"  % (outFileHandle))
extrAntigMapFile.close()
refMapFile.close()
agg_lineListFile.close()
