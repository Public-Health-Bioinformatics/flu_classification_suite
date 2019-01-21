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
 
# General How-To

Task | Action
------------ | -------------
Upload a fasta file | <ul><li> Select **Get Data** from the **Tools** menu </li><li> Select **Upload File from your computer**</li><li> Drag the fasta file(s) into the window </li><li> Select **Start** to upload the file(s) </li><li> Select **Close** once each file's green progress bar reads 100% </li><li> Collapse the **Get Data** menu by selecting it </li></ul>	
Upload a comma-separated value (csv) file | <ul><li> Select **Get Data** from the **Tools** menu </li><li> Select **Upload File from your computer** </li><li> Drag the csv file(s) into the window </li><li> Select "csv" under the **Type** column </li><li> Select **Start** to upload the file(s) </li><li> Select **Close** once each file's green progress bar reads 100% </li><li> Collapse the **Get Data** menu by selecting it </li></ul>
View the contents of a file | <ul><li> Select the eye icon by the filename </li></ul>
Remove a file | <ul><li> Select the pencil icon by the filename </li></ul>
Download a file | <ul><li> Select the computer disk icon by the filename </li></ul>
View the metadata of a file | <ul><li> Select the information icon by the filename </li></ul>
Use an individual tool | <ul><li> Select the tool under **Flu Classification Suite** in the **Tools** pane </li><li> Input the required files </li><li>	Execute the operation </li><li>	Download the output file </li></ul>	
Use a workflow | <ul><li> Select the workflow under **Flu Classification Suite** in the **Tools** pane </li><li> Input the required files </li><li> Select **Run workflow** to execute the operations </li><li> Wait for each workflow step to complete (highlighted in green) </li><li> Download or delete output files from each step as desired</li></ul>
Determine if a tool is running | The operation will display as highlighted in grey while it is waiting to start on the server, as yellow during execution, and as green when complete

# Tools 
Template files and sample input and output files can be found in the 'test-data' folder for each respective tool.

Each tool can be selected from the “Flu Classification Suite” menu and used individually or chained together to create a workflow. 

## Aggregate Line List
Transforms fasta files of flu antigenic site amino acids into aggregated line lists, comparing antigenic maps to that of a reference sequence and collapsing and enumerating identical sequences.

Input - Antigenic Site Extraction output (fasta), amino acid index array (csv), clade definition file (csv)  
Output - csv

1.	Select the **Aggregate Line List** tool
2.	Input the fasta file with extracted antigenic maps and clade calls under **Sample Sequences fasta**
3.	Input the fasta file with the reference sequence antigenic map under **Reference Sequence fasta**
4.	Select the index array under **Antigenic Site Index Array File** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the clade definition file under **Clade Definition File** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
6.	Execute the operation
7.	Download the output file

## Antigenic Site Extraction
Extracts antigenic amino acids from influenza hemagglutinin (HA) sequences, using a flu type-specific array of amino acid positions to be extracted (i.e. for H3, H1 etc.), and outputs as a fasta file.

Input - Assign Clades output (fasta), amino acid index array (csv)  
Output - fasta, csv

1.	Select the **Antigenic Site Extraction** tool
2.	Select the fasta file with protein sequences to extract under **input_fasta**
3.	Select the antigenic site index array file under **index_array** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Choose **Yes** to output results in csv format. (*Note: The default **No** selection outputs results in fasta*)
5.	Execute the operation
6.	Download the output file

## Assign Clades
Assigns clade designations to influenza HA amino acid fasta files.

Input - Sequence files (fasta), clade definition file (csv)  
Output - fasta

1. Select the **Assign Clades** tool
2.	Select the fasta file containing the amino acid sequences under **input_fasta**
3.	Select the clade definition file under **clade_definitions** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Execute the operation
5.	Download the fasta file with the clade-called sequences

## Change Fasta Deflines
Renames definition lines in fasta files. Requires a fasta file requiring sequence name changes and a 2-column renaming file (either tab-delimited text or csv). Searches for fasta definition lines matching column 1 and, if found, replaces fasta definition line with string specified in column 2 of the renaming file.

Input - Sequence file to be renamed (fasta), 2-column renaming file (txt or csv)  
Output - fasta

1. Create a renaming file in Excel with current sequence names in column 1 and desired names in column 2 and export in csv or tab-delimited text format
2. Upload this file into Galaxy
3. *If the renaming file is tab-delimited text, select the pencil icon beside the file name, select datatypes, and ensure the datatype displayed is "csv"*
4. Select the **Change Fasta Deflines** tool
5.	Choose a fasta file under the **“input_fasta”** parameter
6.	Choose a renaming file under the **“key_value_pairs”** parameter
7.	Select **Yes/No** under the **“Names file is tab-delimited”** (*Note: The default renaming file format is csv and the default selection is set to “No”*)
8.	Press **Execute** to start the operation

## Line List
Transforms fasta files of flu antigenic site amino acids into line lists, comparing antigenic maps to that of a reference sequence.

Input - Antigenic Site Extraction output (fasta), amino acid index array (csv), clade definition file (csv)  
Output - csv

1.	Select the **Line List** tool
2.	Input the fasta file with extracted antigenic maps and clade calls under **Sample Sequences fasta**
3.	Input the fasta file with the reference sequence antigenic map under **Reference Sequence fasta**
4.	Select the index array under **Antigenic Site Index Array File** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the clade definition file under **Clade Definition File** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
6.	Execute the operation
7.	Download the output file


## Reformat USearch-Collapsed Fasta
Parses format of USearch-collapsed fasta output files and outputs fasta with customized definition line formatting.

# Workflows
While each tool could be selected from the “Flu Classification Suite” menu and used individually, a workflow was created by chaining tools in a pipeline to automate a series of tasks in a standardized, user-friendly manner.

## Assign clades and extract antigenic maps to csv
Input - Sequence files (fasta), clade definition file (csv), amino acid index array (csv)  
Output - csv
