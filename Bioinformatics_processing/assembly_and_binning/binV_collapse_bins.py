#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd
import os

###function for concatenating fasta files of species
def concat_species_bins(sample_name,species_name,genome_length):

    fasta_out = open('%s_%s_len%s_concat.fa' %(sample_name,species_name,genome_length), "w")

    for filename in file_list:

        if filename.endswith('.fa'):

            file_split = filename.split("_")
            if file_split[4] == 'BVAB1':
                scientific_name = file_split[4]
            elif file_split[4] == 'g':
                scientific_name = file_split[5]
            else:
                scientific_name = file_split[4] + "_" + file_split[5]
        
            if scientific_name==species_name:
                
                for rec in SeqIO.parse(open(filename, "r"),"fasta"):
                    fasta_out.write(">" + str(rec.id) + '\n'+str(rec.seq) + '\n')

    fasta_out.close()

###reading in the genome bin file
genome_size = pd.read_csv("../virgo_genome_sizes.csv",sep=",")

###
file_list = []
for file in os.listdir('.'):
    #check if the file is an excel file
    if file.endswith('.fa'):
        file_list.append(file)

bins_by_species = {}

for filename in file_list:
    file_split = filename.split("_")
    if file_split[4] == 'BVAB1':
        scientific_name = file_split[4]
        ength = int(file_split[6].split('n')[1])
    elif file_split[4] == 'g':
        scientific_name = file_split[5]
        ength = int(file_split[7].split('n')[1])
    else:
        scientific_name = file_split[4] + "_" + file_split[5]
        length = int(file_split[7].split('n')[1])
    bins_by_species.setdefault(scientific_name, []).append(length)

species_bins = {}
sample_name = file_split[0] + '_' + file_split[1] + '_' + file_split[2] + '_' + file_split[3]
#summing bin lengths of species
for k, v in bins_by_species.items():

    if len(v) > 1:
        
        species_bins[k] = sum(v)
        print(k)

        species_genome_size = float(genome_size.loc[genome_size['Species'] == k, 'GenomeSize'])

        if species_genome_size*1000000*1.2 > species_bins[k] > species_genome_size*1000000*.8:

            concat_species_bins(sample_name,k,species_bins[k])

            