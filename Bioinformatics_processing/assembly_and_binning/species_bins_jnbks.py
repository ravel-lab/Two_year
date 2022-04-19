#!/usr/bin/python3

#loading in required python modules
from Bio import SeqIO
import sys
import matplotlib.pyplot as pyplot
from statistics import median
import pandas as pd
import os
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans
import numpy as np
from scipy import stats
import jenkspy

np.seterr(divide='ignore', invalid='ignore')

#defining a function that pulls out all of the entries belonging to a taxa
def output_bin(taxa_dataframe,taxa,sample_name,coverage_bin):

    #creating a list of needed contigs
    needed_contigs = list(taxa_dataframe['Contig'].values)

    contig_name_key = {}
    #cleaning up the contig names and generating variables to use in bin name
    total_length = 0
    total_cov = 0
    contig_count = 0

    for contig_name in needed_contigs:

        components = contig_name.split("_")
        total_length += int(components[3])
        total_cov += int(taxa_dataframe.loc[taxa_dataframe['Contig'] == contig_name,'Coverage'].values[0])
        contig_count+=1

        new_name = components[0] + "_" + components[1] + "_Len" + components[3] + "_Cov" + str('{0:.2f}'.format(taxa_dataframe.loc[taxa_dataframe['Contig'] == contig_name,'Coverage'].values[0])) + "_Per" + str('{0:.2f}'.format(taxa_dataframe.loc[taxa_dataframe['Contig'] == contig_name,'Taxa_perc_map'].values[0]))
        contig_name_key[contig_name] = new_name

    ave_cov = int(float(total_cov)/contig_count)

    species_bin_fa = open('%s_%s_%s_len%s_cov%s.fa' %(sample_name,taxa,coverage_bin,total_length,ave_cov), "w")

    for rec in SeqIO.parse(open(sys.argv[2], "r"),"fasta"):

        if str(rec.id) in needed_contigs:

            species_bin_fa.write(">" + str(contig_name_key[rec.id]) + '\n'+str(rec.seq) + '\n')

    species_bin_fa.close()

def coverage_clustering(taxa_data):

    #test cluster numbers 2 to 10
    cluster_range = [2,3,4,5,6,7,8,9,10]
    taxa_data = taxa_data.sort_values('Coverage',axis=0,ascending=False)
    #use only contigs at least 10kb in the clustering
    taxa_data_trim = taxa_data[taxa_data['Contig_len']>=10000]

    cov = np.array(taxa_data_trim['Coverage'])
    cov_r = cov.reshape(-1,1)    
    
    silhouette_scores = []
    
    #identify breaks in coverage for each cluster number using jenks natural breaks, calculate silhouette score and output in array
    for cluster in cluster_range:

        breaks = jenkspy.jenks_breaks(taxa_data_trim['Coverage'], nb_class=cluster)
        labels  = pd.cut(taxa_data_trim['Coverage'],bins=breaks,labels=list(range(cluster)),include_lowest=True).to_numpy()

        silhouette_scores.append(metrics.silhouette_score(cov_r, labels,metric='braycurtis'))

    #take the max silhouette score
    max_score = max(silhouette_scores)

    equitable_scores = []

    for clust_score in silhouette_scores:

        pos = silhouette_scores.index(clust_score)

        if np.abs((max_score-clust_score)/max_score) <= 0.01:

           equitable_scores.append(pos)

    max_score_n = min(equitable_scores) + 2

    breaks = jenkspy.jenks_breaks(taxa_data_trim['Coverage'], nb_class=max_score_n)

    breaks[0] = taxa_data['Coverage'].min()
    breaks[-1] = taxa_data['Coverage'].max()

    taxa_data['cluster'] = pd.cut(taxa_data['Coverage'],labels=list(range(max_score_n)),bins=breaks,include_lowest=True)

    return taxa_data


contig_data = pd.read_csv(sys.argv[1],sep="\t")
sample_name = sys.argv[3]


os.makedirs(sample_name,exist_ok=True)
os.chdir(sample_name)

for taxa in sorted(set(contig_data['Taxa'].values)):

    if taxa == 'unknown':

        continue

    #filtering out contigs that don't completely map to one species in VIRGO or have low coverage
    taxa_data = contig_data[contig_data['Taxa']==taxa]
    #contig must have at least 90% of taxa id'd genes going to the same taxa
    taxa_data = taxa_data[taxa_data['Taxa_perc_map'] >= 0.90]
    #contig must have at least 5x coverage
    taxa_data = taxa_data[taxa_data['Coverage']>= 5]

    #if there are more than 15 contigs check for outliers and save in a phage category
    if len(taxa_data.index) > 15:
        #filting out outlier contigs with extremely high coverage, labeling them as maybe phage
        
        pot_phage_data = taxa_data[(np.abs(stats.zscore(taxa_data['Coverage']))>6)]
        taxa_data = taxa_data[(np.abs(stats.zscore(taxa_data['Coverage']))<6)]

        if len(pot_phage_data.index) > 0:

            output_bin(pot_phage_data,taxa,sample_name,'phage')



    #calculate max difference in coverage
    cov_dif = taxa_data['Coverage'].max() - taxa_data['Coverage'].min()
    cov_difper = (taxa_data['Coverage'].max() - taxa_data['Coverage'].min())/taxa_data['Coverage'].max()

    #if the length of the contig is more than 1Mb and there are at least 15 10kb contigs and the percent difference in coverage is at least 50% and the total difference in coverage is at least 20
    if (taxa_data['Contig_len'].sum() >= 1000000) & (len(taxa_data[taxa_data['Contig_len']>=10000].index)>15) & (cov_difper >.5) & (cov_dif>20):     

        taxa_data = coverage_clustering(taxa_data)

        for coverage_bin in sorted(set(taxa_data['cluster'].values)):

            coverage_bin_data = taxa_data[taxa_data['cluster']==coverage_bin]
        
            output_bin(coverage_bin_data,taxa,sample_name,coverage_bin)

    #if the combined length is at least 500kb
    elif taxa_data['Contig_len'].sum() >= 500000:

        output_bin(taxa_data,taxa,sample_name,0)
    
