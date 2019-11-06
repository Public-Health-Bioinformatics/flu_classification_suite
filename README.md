# Influenza Classification Suite
Influenza viruses continually evolve to evade population immunity. We have developed a publicly-available Galaxy workflow Influenza Classification Suite, for rapid clade-mapping of sequenced influenza viruses. This suite provides rapid, high-resolution understanding of circulating influenza strain evolution to inform influenza vaccine effectiveness and the need for potential vaccine reformulation. 

# Table of Contents
 * [Installation](#installation)
 * [General How-To](#general-how-to)
 * [Tools](#tools)
    * [Change Fasta Deflines](#change-fasta-deflines)
    * [Assign Clades](#assign-clades)
    * [Antigenic Site Extraction](#antigenic-site-extraction)
    * [Line List](#line-list)
    * [Aggregate Line List](#aggregate-line-list)
    * [Reformat USearch-Collapsed Fasta](#reformat-usearch-collapsed-fasta)
 * [Workflows](#workflows)
    * [Assign clades and extract antigenic maps](#assign-clades-and-extract-antigenic-maps)
    * [Assign clades, extract antigenic maps and output to line list](#assign-clades-extract-antigenic-maps-and-output-to-line-list)
    * [Assign clades, extract antigenic maps and output to aggregated line list](#assign-clades-extract-antigenic-maps-and-output-to-aggregated-line-list)
 
# Installation
 1. In the Galaxy Admin panel, select 'Install New Tools'
 2. Select the 'Galaxy Test Toolshed'
 3. Search for `flu_classification_suite`
 4. Click the button labeled `flu_classification_suite` and select 'Preview and install'
 5. Click 'Install to Galaxy'
 6. Select a tool panel section to install the tools under, or create a new section. We recommend creating a section called 'Influenza Classification Suite'
 7. Click 'Install' 
 
# General How-To

Task | Action
------------ | -------------
[Upload a fasta file](https://github.com/Public-Health-Bioinformatics/flu_classification_suite/blob/master/doc/images/uploading_data_to_galaxy.png) | <ul><li> Select **Get Data** from the **Tools** menu </li><li> Select **Upload File from your computer**</li><li> Drag the fasta file(s) into the window </li><li> Select **Start** to upload the file(s) </li><li> Select **Close** once each file's green progress bar reads 100% </li><li> Collapse the **Get Data** menu by selecting it </li></ul>	
Upload a comma-separated value (csv) file | <ul><li> Select **Get Data** from the **Tools** menu </li><li> Select **Upload File from your computer** </li><li> Drag the csv file(s) into the window </li><li> Select "csv" under the **Type** column </li><li> Select **Start** to upload the file(s) </li><li> Select **Close** once each file's green progress bar reads 100% </li><li> Collapse the **Get Data** menu by selecting it </li></ul>
[View the contents of a file](https://github.com/Public-Health-Bioinformatics/flu_classification_suite/blob/master/doc/images/viewing_data_in_galaxy.png) | <ul><li> Select the eye icon by the filename </li></ul>
Edit a file | <ul><li> Select the pencil icon by the filename </li></ul>
Remove a file | <ul><li> Select the `X` icon by the filename </li></ul>
Download a file | <ul><li> Select the computer disk icon by the filename </li></ul>
View the metadata of a file | <ul><li> Select the information icon by the filename </li></ul>
Use an individual tool | <ul><li> Select the tool under **Flu Classification Suite** in the **Tools** pane </li><li> Input the required files </li><li>	Execute the operation </li></ul>	
Use a workflow | <ul><li> Select the workflow under **Flu Classification Suite** in the **Tools** pane </li><li> Input the required files </li><li> Select **Run workflow** to execute the operations </li></ul>
Determine if a tool is running | The operation will display as highlighted in grey while it is waiting to start on the server, as yellow during execution, and as green when complete

New to Galaxy? Try [Galaxy 101 tutorial](https://galaxyproject.org/tutorials/g101/)

# Tools 
Template files and sample input and output files can be found in the 'test-data' folder for each respective tool.

Each tool can be selected from the “Influenza Classification Suite” menu and used individually or chained together to create a workflow. 

## Change Fasta Deflines
Renames definition lines in fasta files. Requires a fasta file requiring sequence name changes and a 2-column renaming file (either [tab-delimited text](doc/images/change_fasta_deflines-tsv_renaming-file.tsv) or [csv](doc/images/change_fasta_deflines-csv_renaming-file.csv)). Searches for fasta definition lines matching column 1 and, if found, replaces fasta definition line with string specified in column 2 of the renaming file.

![](doc/images/change-fasta-deflines-tool-usage-with-csv.png)

**Input** - Sequence file to be renamed (fasta), 2-column renaming file ([csv](tools/change_fasta_deflines/test-data/csv_rename_file.csv) or [txt](tools/change_fasta_deflines/test-data/tab_delim_rename_file.txt))  
**Output** - fasta

**Command line usage**

   Using csv renaming file:
   
      python change_fasta_def_lines.py csv_rename_file.csv fasta_2_rename.fasta renamedSequences.fasta

   Using tsv renaming file:
   
      python change_fasta_def_lines.py tab_delim_rename_file.txt -t fasta_2_rename.fasta renamedSequences.fasta

**Galaxy tool usage**
1. Create a renaming file in Excel with current sequence names in column 1 and desired names in column 2 and export in [csv](doc/images/change-fasta-deflines-csv-seq-rename-file.png) or [tsv](doc/images/change-fasta-deflines-tsv-seq-rename-file.png) format
2. Upload this file into Galaxy
3. *If the renaming file is tab-delimited text, select the pencil icon beside the file name, select datatypes, and ensure the datatype displayed is "csv"*
4. Select the **Change Fasta Deflines** tool
5.	Choose a fasta file under the **“input_fasta”** parameter
6.	Choose a renaming file under the **“key_value_pairs”** parameter
7.	Select **Yes/No** under the **“Names file is tab-delimited”** (*Note: The default renaming file format is csv and the default selection is set to “No”*)
8.	Press **Execute** to start the operation

![](doc/images/change_fasta_deflines-fasta_file_after_changing_definition_lines.png)

## Assign Clades
Assigns clade designations to influenza HA amino acid fasta files.

![](doc/images/assign_clades-tool_image.png)

**Input** - Sequence files (fasta), [clade definition file (csv)](tools/assign_clades/test-data/clades.csv)  
**Output** - fasta

**Command line usage**

    python assign_clades.py input_sequences.fasta clade_definitions.csv clade-assigned-output-sequences.fasta

**Galaxy tool usage**
1. Select the **Assign Clades** tool
2.	Select the [fasta file (without clades assigned)](doc/images/assign_clades-fasta_file_before_assigning_clades.png) containing the amino acid sequences under **input_fasta**
3.	Select the clade definition file under **clade_definitions** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Execute the operation
5.	Download the output file

![](doc/images/assign_clades-fasta_file_after_assigning_clades.png)

## Antigenic Site Extraction
Extracts antigenic amino acids from influenza hemagglutinin (HA) sequences, using a flu type-specific array of amino acid positions to be extracted (i.e. for H3, H1 etc.), and outputs as a fasta file.

![](doc/images/antigenic-site-extraction_tool_usage.png)

**Input** - Assign Clades output (fasta), [amino acid index array (csv)](tools/antigenic_site_extraction/test-data/)  
**Output** - fasta, csv

**Command line usage**

   Output extracted antigenic amino acids to fasta:
   
      python antigenic_site_extraction.py input.fasta index-array.csv extracted-antigenic-sites.fasta
   
   Output extracted antigenic amino acids to csv:
   
      python antigenic_site_extraction.py input.fasta index-array.csv -c extracted-antigenic-sites.csv

**Galaxy tool usage**
1.	Select the **Antigenic Site Extraction** tool
2.	Select the fasta file with protein sequences to extract under **input_fasta**
3.	Select the antigenic site index array file under **index_array** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Choose **Yes** to output results in csv format (*Note: The default **No** selection outputs results in fasta*)
5.	Execute the operation
6.	Download the output file

![](doc/images/antigenic-site-extracted-fasta.png)

## Line List
Transforms fasta files of flu antigenic site amino acids into line lists, comparing antigenic maps to that of a reference sequence.

![](doc/images/line-list-tool-usage.png)

**Input** - Antigenic Site Extraction output (fasta), [amino acid index array (csv)](tools/linelisting/test-data/FluA_H3_antigenic_aa_indices.csv), [clade definition file (csv)](tools/linelisting/test-data/Flu_Clade_Definitions_H3_20171121.csv)  
**Output** - csv

**Command line usage**

  `python linelisting.py input.fasta reference-antigenic-sites.fasta index-array.csv clade-definitions.csv output.csv`

**Galaxy tool usage**
1.	Select the **Line List** tool
2.	Input the fasta file with extracted antigenic maps and clade calls under **Sample Sequences fasta**
3.	Input the fasta file with the reference sequence antigenic map under **Reference Sequence fasta**
4.	Select the index array under **Antigenic Site Index Array File** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the clade definition file under **Clade Definition File** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
6.	Execute the operation
7.	Download the output file
8. View the line list as a [spreadsheet](doc/images/spreadsheet_view_linelist_and_agg_linelist.pdf) to compare antigenic amino acid sequences of samples to the reference.

![](doc/images/lline_list-output.png)

## Aggregate Line List
Transforms fasta files of flu antigenic site amino acids into aggregated line lists, comparing antigenic maps to that of a reference sequence and collapsing and enumerating identical sequences.

![](doc/images/aggLL-align_sample_antigenic_sites_to_ref.png)

**Input** - Antigenic Site Extraction output (fasta), [amino acid index array (csv)](tools/aggregate_linelisting/test-data/FluA_H3_antigenic_aa_indices.csv), [clade definition file (csv)](tools/aggregate_linelisting/test-data/Flu_Clade_Definitions_H3_20171121.csv)   
**Output** - csv

**Command line usage**

  `python aggregate_linelisting.py input.fasta reference-antigenic-sites.fasta index-array.csv clade-definitions.csv output.csv`


**Galaxy tool usage**
1.	Select the **Aggregate Line List** tool
2.	Input the fasta file with extracted antigenic maps and clade calls under **Sample Sequences fasta**
3.	Input the fasta file with the reference sequence antigenic map under **Reference Sequence fasta**
4.	Select the index array under **Antigenic Site Index Array File** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the clade definition file under **Clade Definition File** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
6.	Execute the operation
7.	Download the output file
8. View the line list as a [spreadsheet](doc/images/spreadsheet_view_linelist_and_agg_linelist.pdf) to compare antigenic amino acid sequences of samples to the reference.

![](doc/images/aggregate_line_list-output.png)

## Reformat USearch-Collapsed Fasta
Parses format of USearch-collapsed fasta output files and applies a custom format to the fasta definition lines.

![](doc/images/reformat-usearch-collapsed-fasta-usage.png)

**Input** - [USearch-outputted sequence files (fasta)](tools/reformat_usearch_collapsed_fasta/test-data/10_usearch_collapsed_sequences.fasta)

**Output** - [sequence files with custom-formatted definition lines (fasta)](tools/reformat_usearch_collapsed_fasta/test-data/output.fasta)

**Command line usage**

    python reformat_usearch_collapsed_fasta.py usearch_collapsed_sequences.fasta output.fasta

**Galaxy tool usage**
1.	Select the **Reformat USearch-Collapsed Fasta** tool
2.	Select the [USearch-collapsed fasta file to reformat](doc/images/usearch-collapsed-fasta-before-reformat.png) under **input_fasta**
3.	Execute the operation
4.	Download the output file

![](doc/images/reformat-usearch-collapsed-fasta-output.png)

# Workflows
While each tool could be selected from the “Influenza Classification Suite” menu and used individually, a workflow was created by chaining tools in a pipeline to automate a series of tasks in a standardized, user-friendly manner.

## Assign clades and extract antigenic maps

This workflow assigns and appends clade names to fasta definition lines of flu HA amino acid sequences. It then extracts the antigenic amino acids and outputs the resulting antigenic maps in fasta format.

![](doc/images/assign-clades-and-extract-antigenic-maps-workflow.png)

**Input** - Sequence files (fasta), clade definition file (csv), amino acid index array (csv) (*Note: Use the provided clade definition and amino acid index array files or provide your own respective versions of these files*)

**Output** - csv

**Galaxy tool usage**
1.	Select the workflow **Assign clades and extract antigenic maps**
2.	Select whether to send the results to a new history (*Note: This is not required but facilitates convenient tracking and deletion of files within an analysis run*)
3.	Select the clade definition file under **Clade Definitions** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Select the index array under **Antigenic Amino Acid Index Array** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the fasta file to perform all operations as the **input_fasta** under **Assign Clades**
6.	Select **Run workflow** to execute the operations
7.	Wait for each workflow step to complete (*highlighted in green*)
8. Download or delete output files from each step as desired

![](doc/images/Assign-clades-and-extract-antigenic-maps-workflow-usage.png)

## Assign clades, extract antigenic maps and output to line list

This workflow assigns and appends clade names to fasta definition lines of flu HA amino acid sequences, extracts the antigenic sites and displays the resulting sequences in relation to a reference (e.g. vaccine strain) sequence in a csv file.

![](doc/images/assign_clades_extract_antigenic_sites_output_to_line_list-workflow.png)

**Input** - Sequence files (fasta), reference antigenic map (e.g. vaccine influenza strain) (fasta), clade definition file (csv), amino acid index array (csv) (*Note: Use the provided reference antigenic map, clade definition and amino acid index array files or provide your own respective versions of these files*)

**Output** - csv

**Galaxy tool usage**
1.	Select the workflow **Assign clades, extract antigenic maps and output to line list**
2.	Select whether to send the results to a new history (*Note: This is not required but facilitates convenient tracking and deletion of files within an analysis run*)
3.	Select the clade definition file under **Clade Definitions** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Select the index array under **Antigenic Amino Acid Index Array** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the reference antigenic map under **Reference Antigenic Map** (*Note: The reference antigenic map file must be the extracted antigenic site amino acids in fasta format and of the correct version for the respective flu type to obtain accurate results*)
6. Select the fasta file to perform all operations as the **input_fasta** under **Assign Clades**
7.	Select **Run workflow** to execute the operations
8.	Wait for each workflow step to complete (*highlighted in green*)
9. Download or delete output files from each step as desired

![](doc/images/assign_clades_extract_antigenic_sites_output_to_line_list-usage.png)

## Assign clades, extract antigenic maps and output to aggregated line list

This workflow assigns and appends clade names to fasta definition lines of flu HA amino acid sequences, extracts antigenic sites and displays the resulting sequences in relation to a reference (e.g. vaccine strain) sequence in a csv file. In addition, the aggregated view collapses and enumerates identical antigenic map sequences among samples.

![](doc/images/assign_clades_extract_antigenic_sites_output_to_aggregated_line_list-workflow.png)

**Input** - Sequence files (fasta), reference antigenic map (e.g. vaccine influenza strain) (fasta), clade definition file (csv), amino acid index array (csv) (*Note: Use the provided reference antigenic map, clade definition and amino acid index array files or provide your own respective versions of these files*)

**Output** - csv

**Galaxy tool usage**
1.	Select the workflow **Assign clades, extract antigenic maps and output to aggregated line list**
2.	Select whether to send the results to a new history (*Note: This is not required but facilitates convenient tracking and deletion of files within an analysis run*)
3.	Select the clade definition file under **Clade Definitions** (*Note: The clade definition file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
4.	Select the index array under **Antigenic Amino Acid Index Array** (*Note: The index array file must be in csv format and of the correct version for the respective flu type to obtain accurate results*)
5.	Select the reference antigenic map under **Reference Antigenic Map** (*Note: The reference antigenic map file must be the extracted antigenic site amino acids in fasta format and of the correct version for the respective flu type to obtain accurate results*)
6. Select the fasta file to perform all operations as the **input_fasta** under **Assign Clades**
7.	Select **Run workflow** to execute the operations
8.	Wait for each workflow step to complete (*highlighted in green*)
9. Download or delete output files from each step as desired

![](doc/images/assign_clades_extract_antigenic_sites_output_to_aggregated_line_list-workflow-usage.png)
