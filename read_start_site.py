# libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from Bio.Seq import Seq

def TIR_extraction(gff_genes, gene_lst, fasta_str, output, basename):
    file = open(gene_lst, 'r')
    genes = file.read()
    file.close()
    genes = genes.split(', ')
    genes[-1] = genes[-1].replace('\n', '')
    f = open(output + '/' + basename + '_windows.fasta', 'w')
    
    for x in genes:
        for y in gff_genes:
            
            if (y['strand'] == '+') and (y['gene'] == x):
                start = y['start'] - 31
                stop = y['start'] + 89
                for i in range(start, stop, 1):
                    end = i + 30
                    region_1 = i - y['start'] +1
                    region_2 = end - y['start'] +1
                    if region_2 > 90:
                        break
                    f.write('>'+y['gene']+'_'+str(region_1)+'_'+str(region_2) +'\n')
                    f.write(fasta_str[i:end] +'\n')
                for i in range(start, stop, 1):
                    end = i + 40
                    region_1 = i - y['start'] +1
                    region_2 = end - y['start'] +1
                    if region_2 > 90:
                        break
                    f.write('>'+y['gene']+'_'+str(region_1)+'_'+str(region_2) +'\n')
                    f.write(fasta_str[i:end] +'\n')
                for i in range(start, stop, 1):
                    end = i + 50
                    region_1 = i - y['start'] +1
                    region_2 = end - y['start'] +1
                    if region_2 > 90:
                        break
                    f.write('>'+y['gene']+'_'+str(region_1)+'_'+str(region_2) +'\n')
                    f.write(fasta_str[i:end] +'\n')
                for i in range(start, stop, 1):
                    end = i + 60
                    region_1 = i - y['start'] +1
                    region_2 = end - y['start'] +1
                    if region_2 > 90:
                        break
                    f.write('>'+y['gene']+'_'+str(region_1)+'_'+str(region_2) +'\n')
                    f.write(fasta_str[i:end] +'\n')

            if (y['strand'] == '-') and (y['gene'] == x):
                    
                start = y['stop'] + 30
                stop = y['stop'] - 90
                for i in range(stop, start, 1):
                    end = i + 30
                    fasta = Seq(fasta_str[i:end])
                    fasta = fasta.reverse_complement()
                    region_1 = -1*(i - y['stop'])
                    region_2 = -1*(end - y['stop'])
                    if (region_1 > 90) or (region_2 < -30):
                        break
                    f.write('>'+y['gene']+'_'+str(region_2)+'_'+str(region_1) +'\n')
                    f.write(str(fasta) +'\n')
                for i in range(stop, start, 1):
                    end = i + 40
                    region_1 = -1*(i - y['stop'])
                    region_2 = -1*(end - y['stop'])
                    if (region_1 > 90) or (region_2 < -30):
                        break
                    fasta = Seq(fasta_str[i:end])
                    fasta = fasta.reverse_complement()
                    f.write('>'+y['gene']+'_'+str(region_2)+'_'+str(region_1) +'\n')
                    f.write(str(fasta) +'\n')
                    #print('>'+y['gene']+'_'+str(region_2)+'_'+str(region_1) +'\n')
                    #print(str(fasta) +'\n')
                for i in range(stop, start, 1):
                    end = i + 50
                    region_1 = -1*(i - y['stop']) 
                    region_2 = -1*(end - y['stop']) 
                    if (region_1 > 90) or (region_2 < -30):
                        break
                    fasta = Seq(fasta_str[i:end])
                    fasta = fasta.reverse_complement()
                    f.write('>'+y['gene']+'_'+str(region_2)+'_'+str(region_1) +'\n')
                    f.write(str(fasta) +'\n')
                for i in range(stop, start, 1):
                    end = i + 60
                    region_1 = -1*(i - y['stop'])
                    region_2 = -1*(end - y['stop'])
                    if (region_1 > 90) or (region_2 < -30):
                        break
                    fasta = Seq(fasta_str[i:end])
                    fasta = fasta.reverse_complement()
                    f.write('>'+y['gene']+'_'+str(region_2)+'_'+str(region_1) +'\n')
                    f.write(str(fasta) +'\n')
    f.close()