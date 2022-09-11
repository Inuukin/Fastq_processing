#!/usr/bin/env python

### Fastq processing

from Bio import SeqIO ### I am using biopython here
import statistics as stat
import pandas as pd
import matplotlib as mappet ## :D
import matplotlib.pyplot as plt
import numpy as np
import unittest

def Fastq_proc():
    DF = pd.DataFrame()
    DF_stat = pd.DataFrame()
    count = 1

    for record in SeqIO.parse("test.fastq", "fastq"):
            x = record.letter_annotations["phred_quality"]
            df = pd.DataFrame (x, columns = [count])
            DF[count] = list(x)
            count += 1

    ## create a dataframe as required
    DF_stat['read_position'] = np.arange(len(DF)) 
    DF_stat['mean_Phred_qual'] = DF.mean(axis=1) 
    DF_stat['standard_deviation_Phred_qual'] = DF.std(axis =1)


    Graph = DF_stat.plot(kind = "barh", y = "mean_Phred_qual", legend = False, title = "Average read quality",
        xerr = "standard_deviation_Phred_qual")

    Graph.set(ylabel="Base position", xlabel="Phred")
    plt.savefig('Barplot_read_quality.pdf')  

    ## save to tsv file

    DF_stat.to_csv('output.tsv', sep="\t", index = False)
    
if __name__ == "__main__":
    Fastq_proc()
