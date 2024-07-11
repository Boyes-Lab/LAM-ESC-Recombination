LAM Recombination analysis script

README

Time required to download script and dependencies: Approximately 5 minutes
Time required to run on a typical sample: <30 minutes

Required: python3 with multiprocessing module: install using 'python3 -m pip install multiprocessing' if required

Required dependencies:

FLASH version 1.2.11:

https://ccb.jhu.edu/software/FLASH/

Also available to download from: https://sourceforge.net/projects/flashpage/files/
Place FLASH binary into same folder as this script, named 'flash' (location can be changed by modifying line 23 of the script)

blastn:

Available to download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/
Tested with version 2.12.0
Place blastn binary into same folder as this script, named 'blastn' (location can be changed by modifying line 21 of the script)


Custom BLAST database:
The custom BLAST database folder (blast_db) should be located in the same folder as the script. This is a BLASTN database containing all combinations of V-J recombination events that can be generated at the immunoglobulin human light chain loci. This database is available in this Github repository along with the Python script.


System requirements:
This script has been run and tested on an Ubuntu 22.04 machine with an Intel i7-6700K processor and 16 GB RAM, using versions of FLASH and blastn compiled for linux (Intel).


Workflow:
1. R1 and R2 reads from Fastq files are joined using FLASH to generate a single overlapping read (if possible; if not, the R1 read is analysed)
2. Reads are demultiplexed based on each J gene segment (or KDE) primer based on a simple sequence match
3. Reads from unrecombined J gene segments are filtered out via BLAST match to custom databases containing J gene segments
4. Remaining reads are BLASTed against custom BLAST databases of V gene segment RSSs
5. Matches are then BLASTed against custom BLAST databases of V-J SJs

*****
Usage:
python3 lam_rec_script.py SAMPLE_NAME ANALYSIS_FOLDER_PATH

Argument 1 = SAMPLE_NAME
This will be prefixed to all output files

Argument 2 = ANALYSIS_FOLDER_PATH
Provide the full path to the analysis folder containing the two paired end .fastq files (or fastq.gz)
This can be found by changing to the analysis folder and running the 'pwd' command 
*****

Output files:
Individual tsv (tab separated values) files are produced for each J (or KDE) RSS:

*_blast_result_top_two.tsv
These contain the top two hits for each read (can be useful in the case of equal matches to inversional
and deletional V gene segments at the kappa locus)

*_blast_result_top_hit.tsv
This contains the top recombination blast hit for each read (based on e value score) *careful that some V kappa gene segments can match identically to deletional and inversional segments and therefore have the same e value

*_blast_result_top_hit.tsv
This contains the top recombination blast hit for each read (based on e value score) *careful that some V kappa gene segments can match identially to deletional and inversional segments and therefore have the same e value

*_Rec_summary.txt
A summary file with counts of each recombination event detected, generated from the *_top_hit.tsv files

*_blast_result_top_hit_e_val.tsv
This contains the top SJ blast hit for each read (based on e value score); reads where the e value is identical for the top two reads are removed

*_Rec_summary_unique_e_val.txt
A summary file with counts of each recombination event detected, generated from the *_blast_result_top_hit_e_val.tsv files; useful to filter out ambiguous deletional/inversional SJs



