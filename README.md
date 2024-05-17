# BMI8950Deployment
## Install Instructions (if applicable)
### Requirements
* Python

The following accession numbers were used in this project:

|Variant |Subvariant ID |NCBI Virus Accession Number |Protein Data Bank Accession Number |
|--------|-----------|-----------------------------|----------------------------------- |
|Alpha |AY30 |WLY63243 | |
|Delta |B.1.617.2 |WWZ28257 |7SBK |
|Delta |B.1.1.7 |WWQ76894 |7SBK |
|Omicron |BA.1 |WWZ28966 |7TNW |
|Omicron |BA.2 |WXB41833 |7XIW |

### Installation Instructions

Download protein sequence files from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:2697049. For each accession number, refine results by accession, enter the accession number into the text box, and hit submit. Select the link to the accession number to open the GenBank page for the submission. Select the "Send to" dropdown in the corner of the page, select 'File', change the format to FASTA, and select 'Create File' to download the sequence. 

Download Spike protein-ACE2 receptor binding complex reference structure from https://www.rcsb.org/structure/6LZG. Select the 'Download files' dropdown menu and select 'PDB format'. To access the SARS-CoV-2 protein data bank files, search for the listed Protein Data Bank Accession Number.

### Getting started
Prior to running the python programs, blosum, Python math, and Python Bio must be installed on the console through the pip command.
