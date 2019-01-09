# Flu Classification Suite
Influenza viruses continually evolve to evade population immunity. We have developed a publicly-available Galaxy workflow Flu Analysis Suite, for rapid clade-mapping of sequenced influenza viruses. This suite provides rapid, high-resolution understanding of circulating influenza strain evolution to inform influenza vaccine effectiveness and the need for potential vaccine reformulation. 

# Installation

# Tools

## Aggregate Line List
Transforms fasta files of flu antigenic site amino acids into aggregated line lists, comparing antigenic maps to that of a reference sequence and collapsing and enumerating identical sequences.

## Antigenic Site Extraction
Extracts antigenic amino acids from flu sequence, using a specific index array (i.e. for H3, H1 etc.).

## Assign Clades
Assign clade designations to influenza amino acid fasta files.

## Change Fasta Deflines
Renames definition lines in fasta files. Requires a fasta file requiring sequence name changes and a 2-column renaming file (either tab-delimited text or csv). Searches for fasta definition lines matching column 1 and, if found, replaces fasta definition line with string specified in column 2 of the renaming file.

## Line List
Transforms fasta files of flu antigenic site amino acids into line lists, comparing antigenic maps to that of a reference sequence.

## Reformat USearch-Collapsed Fasta
Parses format of USearch-collapsed fasta output files and outputs fasta with customized definition line formatting.
