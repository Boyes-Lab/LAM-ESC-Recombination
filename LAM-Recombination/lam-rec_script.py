#!/usr/bin/env python3
# coding: utf-8

import os
import re
import subprocess
import time
import glob
import shlex
from itertools import zip_longest
from multiprocessing import Pool
from collections import Counter
import sys

system_platform = sys.platform

current_dir = os.getcwd()


#path to blastn binary:
blastn = f'{current_dir}/blastn'
#path to FLASH binary:
flash_path = f'{current_dir}/flash'

if sys.platform == 'win32':
    print('Please run this script on linux')
    exit()


cpu_threads = 8

notebook_path = os.getcwd()




#path to blast db
#put the blast_db folder in the same dir that you run the script from
blast_db_path = current_dir


# In[2]:


def run_shell_command(command):
    print(f'Running command: {command} \n')
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output = process.communicate()
    print(output[1].decode("utf-8"))
    print('\n')


# In[39]:


print('Welcome to the LAM-Recombination analysis script')
print('')
print('Requirements: flash https://ccb.jhu.edu/software/FLASH/ ')
print('Required: blastn https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/')
print('Place the above binaries in the same folder as the python script')
print('')
print('Please make a new folder for each sample to be analysed, and run the script from that folder')
print('')
print('Make sure the folder only contains the forward and reverse FASTQ files')
print('Please give the replicate name as the first script argument; otherwise this will default to "sample"')
print('')
print('Please give the full folder path as the second script argument; otherwise this will run from the current directory')
print('')
print('Example: python3 lam-rec.py sample_name folder_path')

sample_name = ''
folder_path = ''

try:
	sample_name = sys.argv[1]
	folder_path = sys.argv[2]
except:
	pass

if not sample_name:
    sample_name = 'sample'

if folder_path:
    os.chdir(folder_path)



print('Running LAM-Recombination analysis! This might take a while....')
print('')


if sys.argv[1] == '--help':
    print('Help shown above; script exiting')
    exit()


# In[13]:




barcodes = {}

barcodes['IGKJ1'] = 'ACTGAGGAAGCAAAG'
barcodes['IGKJ2'] = 'ATCTCCAGCTTGGTC'
barcodes['IGKJ3'] = 'TTACGTTTGATATCC'
barcodes['IGKJ4'] = 'ACGTTTGATCTCCAC'
barcodes['IGKJ5'] = 'CTCCAGTCGTGTCCC'
barcodes['KDE'] = 'AGCTGCAGACTCATG'
barcodes['IGLJ1'] = 'AGAGAGAGAAAACAG'
barcodes['IGLJ23'] = 'ACTCACCTAGGACGG'

J_segments = ['IGKJ1','IGKJ2','IGKJ3','IGKJ4','IGKJ5','KDE','IGLJ1','IGLJ23']


#Demultiplex fastq based on Nested/Red primers for each J RSS


for folder , sub_folders , files in os.walk(os.getcwd()):
    
    for gzipped_fastq in files:
        match = re.search(r'.+fastq.gz$',gzipped_fastq)
        if match:
            command = f'gunzip {gzipped_fastq}'
            run_shell_command(command)

            

for folder , sub_folders , files in os.walk(os.getcwd()):
    
    for fastq_file in files:
        match = re.search(r'.+R1_001.fastq$',fastq_file)
        if match:
            print(fastq_file)
            fastq_r1 = folder+'/'+fastq_file.split('R1')[0]+'R1_001.fastq'
            fastq_r2 = folder+'/'+fastq_file.split('R1')[0]+'R2_001.fastq'
            
            command = f'{flash_path} {fastq_r1} {fastq_r2} -M 250'
            run_shell_command(command)
            
            command = f'mv out.extendedFrags.fastq {sample_name}.extendedFrags.fastq'
            run_shell_command(command)
            
            for j in J_segments:
                fastq_file = f'{sample_name}.extendedFrags.fastq'
                barcode = barcodes[j]

                output = []

                with open(fastq_file,'r') as f:
                    for line in f:
                        fastq_1 = line
                        fastq_2 = f.readline()
                        fastq_3 = f.readline()
                        fastq_4 = f.readline()
                        if barcode in fastq_2:
                            output.append(fastq_1)
                            output.append(fastq_2)
                            output.append(fastq_3)
                            output.append(fastq_4)
                            
                #also check fastq reads that didn't combine (R1 notCombined file):
                with open('out.notCombined_1.fastq','r') as f:
                    for line in f:
                        fastq_1 = line
                        fastq_2 = f.readline()
                        fastq_3 = f.readline()
                        fastq_4 = f.readline()
                        if barcode in fastq_2:
                            output.append(fastq_1)
                            output.append(fastq_2)
                            output.append(fastq_3)
                            output.append(fastq_4)

                with open(f'{j}.fastq','w') as f:
                    for item in output:
                        f.write(item)


# In[33]:


def run_analysis(gene):
    print(f'running for {gene}')
    Rec_blast_result = []
    unmatched = []
    sequences = []
    V_blast_result = []
    with open(f'{gene}.fastq','r') as f:
        for line in f:
            seq_id = '>' + line.split(' ')[0]
            seq = f.readline().strip()
            if len(seq) > 100:
                sequences.append(seq_id)
                sequences.append(seq)
            f.readline()
            f.readline()

    check_unmatched_output = []
    with open(f'{gene}_tmp.fasta','w') as f:
        for item in sequences:
            f.write(f'{item}\n')
            
    evalue = 1e-5
    max_target_seqs = 5
    outfmt = "'6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'"
    command = f'{blastn} -query {gene}_tmp.fasta -db {blast_db_path}/blast_db/{gene}_RSS.fasta -outfmt {outfmt} -evalue {evalue} -max_target_seqs {max_target_seqs} -task blastn -num_threads {cpu_threads} > {gene}_tmp_blast_output.txt'
    os.system(command)
    
    sequence_id_match = []
    with open(f'{gene}_tmp_blast_output.txt','r') as f:
        for line in f:
            seq_id = line.split('\t')[0]
            sequence_id_match.append(seq_id)
    
    sequence_id_match = set(sequence_id_match)
    with open(f'{gene}_tmp.fasta','r') as f:
        for line in f:
            seq_id = line
            seq = f.readline()
            if seq_id.split('>')[1] not in sequence_id_match:
                check_unmatched_output.append(f'{seq_id} {seq}')
    
    os.remove(f'{gene}_tmp.fasta')
    os.remove(f'{gene}_tmp_blast_output.txt')
  

    V_blast_result = []
    fasta = []
    fasta_output = []
    for item in check_unmatched_output:
        seq_id = item.split(' ')[0]
        seq = item.split(' ')[1]
        fasta.append(f'{seq_id}{seq}')
    
    with open(f'{gene}_tmp2.fasta','w') as f:
        for item in fasta:
            f.write(item)
    
    if 'IGK' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGKV_gene_segments.fasta'
    if 'KDE' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGKV_gene_segments.fasta'
    if 'IGL' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGLV_gene_segments.fasta'
    
    
    evalue = 1e-5
    max_target_seqs = 5
    outfmt = "'6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'"
    command = f'{blastn} -query {gene}_tmp2.fasta -db {blast_db} -outfmt {outfmt} -evalue {evalue} -max_target_seqs {max_target_seqs} -task blastn -num_threads {cpu_threads} -mt_mode 1 > {gene}_tmp_blast_output2.txt'
    os.system(command)
    
    sequence_id_match = []
    with open(f'{gene}_tmp_blast_output2.txt','r') as f:
        for line in f:
            seq_id = line.split('\t')[0]
            sequence_id_match.append(seq_id)
            
    with open(f'{gene}_tmp2.fasta','r') as f:
        for line in f:
            fasta_id = line.strip()
            fasta_seq = f.readline().strip()
            if fasta_id.split('>')[1] in sequence_id_match:
                fasta_output.append(f'{fasta_id}\t{fasta_seq}')
                
    os.remove(f'{gene}_tmp2.fasta')
    os.remove(f'{gene}_tmp_blast_output2.txt')
    
            
    Rec_blast_result = []
    fasta = []
    for item in fasta_output:
        seq_id = item.split('\t')[0]
        seq = item.split('\t')[1]
        fasta.append(seq_id)
        fasta.append(seq)
        
    with open(f'{gene}_tmp3.fasta','w') as f:
        for item in fasta:
            f.write(f'{item}\n')
    
    
    evalue = 1e-21
    max_target_seqs = 10
    outfmt = "'6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'"
    
    
    if 'IGK' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGKV_{gene}.fasta'
    if 'KDE' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGKV_{gene}.fasta'
    if 'IGLJ1' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGLV_{gene}.fasta'
    if 'IGLJ23' in gene:
        blast_db = f'{blast_db_path}/blast_db/IGLV_{gene}.fasta'
    
    command = f'{blastn} -query {gene}_tmp3.fasta -db {blast_db} -outfmt {outfmt} -evalue {evalue} -max_target_seqs {max_target_seqs} -task blastn -num_threads {cpu_threads} > {gene}_tmp_blast_output3.txt'
    os.system(command)
    
    command = "sort -t '\t' -k11,11 -g " + f'{gene}_tmp_blast_output3.txt > {gene}_tmp_blast_output_eval_sorted.txt'
    os.system(command)
    
    processed_ids_r1 = []
    processed_ids_r2 = []
    output = []
    with open(f'{gene}_tmp_blast_output_eval_sorted.txt','r') as f:
        for line in f:
            seq_id = line.split('\t')[0]
            if seq_id not in set(processed_ids_r1):
                output.append(line)
                processed_ids_r1.append(seq_id)
            elif seq_id in set(processed_ids_r1):
                if seq_id not in set(processed_ids_r2):
                    output.append(line)
                    processed_ids_r2.append(seq_id)        

                                
    os.remove(f'{gene}_tmp3.fasta')
    os.remove(f'{gene}_tmp_blast_output3.txt')
    
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_two_tmp.tsv','w') as f:
        for item in output:
            f.write(f'{item}')
            
    command = "sort -t '\t' -k1,1 -k11g,11 " + f'{sample_name}_{gene}_Rec_blast_result_top_two_tmp.tsv > {sample_name}_{gene}_Rec_blast_result_top_two.tsv'
    os.system(command)
    os.remove(f'{sample_name}_{gene}_Rec_blast_result_top_two_tmp.tsv')


# In[34]:


from time import time
t0 = time()

genes = ['IGLJ1','IGLJ23','IGKJ1','IGKJ2','IGKJ3','IGKJ4','IGKJ5','KDE']
#genes = ['KDE']

pool = Pool(processes=cpu_threads)
pool.map(run_analysis, genes)

t1 = time()
print (f'Run time: threads={cpu_threads} time: {t1-t0} seconds')


#Generate output files with just the top hit for each sequence


for gene in genes:
    #need to control for only one blast hit in top_two.tsv file; therefore make a list of all IDS and use count
    #to exclude where != 2
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_two.tsv','r') as f:
        ids = [line.split('\t')[0] for line in f]
    output = []
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_two.tsv','r') as a:
        for line in a:
            ident = line.split('\t')[0]
            if ids.count(ident) !=2:
                output.append(line)
            else:
                output.append(line)
                a.readline()
        with open(f'{sample_name}_{gene}_Rec_blast_result_top_hit.tsv','w') as b:
            for item in output:
                b.write(item)

#Count and sort recombination events and print output:

recombination_events = []
for gene in genes:
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_hit.tsv','r') as f:
        for line in f:
            try:
                recombination_event = line.split('\t')[1]
                # E value filter 1e-21
                if float(line.split('\t')[10]) < 1e-21:
                    recombination_events.append(recombination_event)
            except:
                pass
            
x = Counter(recombination_events)
print(x.most_common())

with open(f'{sample_name}_Rec_summary.txt','w') as f:
    for item in x.most_common():
        f.write(f'{item}\n')


#Generate output files with just the top hit for each sequence; 
#removing BOTH entries if the E value is the same = ambiguous

for gene in genes:
    
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_two.tsv','r') as f:
        ids = [line.split('\t')[0] for line in f]
    
    output = []
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_two.tsv','r') as a:
        for line in a:
            ident = line.split('\t')[0]
            if ids.count(ident) !=2:
                output.append(line)
            else:
                match_1 = line
                match_2 = a.readline()
                e_value_1 = match_1.split('\t')[10]
                try:
                    e_value_2 = match_2.split('\t')[10]
                except:
                    pass
                if e_value_1 != e_value_2:
                    output.append(match_1)
        with open(f'{sample_name}_{gene}_Rec_blast_result_top_hit_e_val.tsv','w') as b:
            for item in output:
                b.write(item)

#Count and sort recombination events and print output (same E values/ambiguous matches removed)

recombination_events = []
for gene in genes:
    with open(f'{sample_name}_{gene}_Rec_blast_result_top_hit_e_val.tsv','r') as f:
        for line in f:
            try:
                recombination_event = line.split('\t')[1]
                # E value filter 1e-21
                if float(line.split('\t')[10]) < 1e-21:
                    recombination_events.append(recombination_event)
            except:
                pass
            
x = Counter(recombination_events)
print(x.most_common())

with open(f'{sample_name}_Rec_summary_unique_e_val.txt','w') as f:
    for item in x.most_common():
        f.write(f'{item}\n')
