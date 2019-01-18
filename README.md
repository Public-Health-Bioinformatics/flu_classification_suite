# Flu Classification Suite
Influenza viruses continually evolve to evade population immunity. We have developed a publicly-available Galaxy workflow Flu Classification Suite, for rapid clade-mapping of sequenced influenza viruses. This suite provides rapid, high-resolution understanding of circulating influenza strain evolution to inform influenza vaccine effectiveness and the need for potential vaccine reformulation. 

# Installation
 1. In the Galaxy Admin panel, select 'Install New Tools'
 2. Select the 'Galaxy Test Toolshed'
 3. Search for `flu_classification_suite`
 4. Click the button labeled `flu_classification_suite` and select 'Preview and install'
 5. Click 'Install to Galaxy'
 6. Select a tool panel section to install the tools under, or create a new section. We recommend creating a section called 'Flu Classification Suite'
 7. Click 'Install' 

# Tools 
Template files and sample input and output files can be found in the 'test-data' folder for each respective tool.

Each tool can be selected from the “Flu Classification Suite” menu and used individually or chained together to create a workflow. 

## Aggregate Line List
Transforms fasta files of flu antigenic site amino acids into aggregated line lists, comparing antigenic maps to that of a reference sequence and collapsing and enumerating identical sequences.

Input - Antigenic Site Extraction output (fasta), amino acid index array (csv), clade definition file (csv)  
Output - csv

## Antigenic Site Extraction
Extracts antigenic amino acids from influenza hemagglutinin (HA) sequences, using a flu type-specific array of amino acid positions to be extracted (i.e. for H3, H1 etc.), and outputs as a fasta file.

Input - Assign Clades output (fasta), amino acid index array (csv)  
Output - fasta, csv

## Assign Clades
Assigns clade designations to influenza HA amino acid fasta files.

Input - Sequence files (fasta), clade definition file (csv)  
Output - fasta

## Change Fasta Deflines
Renames definition lines in fasta files. Requires a fasta file requiring sequence name changes and a 2-column renaming file (either tab-delimited text or csv). Searches for fasta definition lines matching column 1 and, if found, replaces fasta definition line with string specified in column 2 of the renaming file.

Input - Sequence file to be renamed (fasta), 2-column renaming file (txt or csv)  
Output - fasta

## Line List
Transforms fasta files of flu antigenic site amino acids into line lists, comparing antigenic maps to that of a reference sequence.

Input - Antigenic Site Extraction output (fasta), amino acid index array (csv), clade definition file (csv)  
Output - csv

## Reformat USearch-Collapsed Fasta
Parses format of USearch-collapsed fasta output files and outputs fasta with customized definition line formatting.

# Workflows
While each tool could be selected from the “Flu Classification Suite” menu and used individually, a workflow was created by chaining tools in a pipeline to automate a series of tasks in a standardized, user-friendly manner.

## Assign clades and extract antigenic maps to csv
Input - Sequence files (fasta), clade definition file (csv), amino acid index array (csv)  
Output - csv
