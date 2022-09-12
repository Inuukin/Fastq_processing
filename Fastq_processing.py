#!/usr/bin/env python

### Fastq processing

from Bio import SeqIO ### I am using biopython here
import statistics as stat
import pandas as pd
import matplotlib as mappet ## :D
import matplotlib.pyplot as plt
import numpy as np


def Fastq_proc():
    
    ## create dataframes
    
    DF = pd.DataFrame()
    DF_stat = pd.DataFrame()
    count = 1
    
    ## read the fasta file in
    
    for record in SeqIO.parse("test.fastq", "fastq"):
            x = record.letter_annotations["phred_quality"]
            df = pd.DataFrame (x, columns = [count])
            DF[count] = list(x)
            count += 1
            
            ## create a dataframe as required
    
    DF_stat['read_position'] = np.arange(1, len(DF) + 1) 
    DF_stat['mean_Phred_qual'] = DF.mean(axis=1) 
    DF_stat['standard_deviation_Phred_qual'] = DF.std(axis =1)
   
    
    ## Graph

    Graph = DF_stat.plot(kind = "barh", x = 'read_position', y = "mean_Phred_qual", legend = False, title = "Average read quality",
        xerr = "standard_deviation_Phred_qual", figsize=(10,20))

    Graph.set(ylabel="Base position", xlabel="Phred")
    plt.savefig('Barplot_read_quality.pdf', dpi=300)  

    ## save to tsv file

    DF_stat.to_csv('output.tsv', sep="\t", index = False)
    
if __name__ == "__main__":
    Fastq_proc()
