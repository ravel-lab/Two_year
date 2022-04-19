#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys
import matplotlib
import seaborn as sns
import random
import matplotlib.patches as mpatches


#taxa color dictionary
taxa_color_scheme = {'Lactobacillus_crispatus':'#ff0000','Gardnerella_vaginalis':'#20b2aa','g_Lactobacillus':'#eef06c',
                     'Lactobacillus_iners':'#ff8c00','Lactobacillus_gasseri':'#7fff00','Lactobacillus_jensenii':'#333333',
                     'Enterococcus_faecalis':'#bbffff','Raoultella_planticola':'#ffc2ff','g_Peptoniphilus':'#CCCC00',
                     'Sneathia_sanguinegens':'#c7aa8f','Atopobium_vaginae':'#0000cd','g_Atopobium':'#0000cd','Prevotella_sp.':'#70005e',
                     'Lacotbacillus_helveticus':'#00ccff','Mageeibacillus_indolicus':'#3cb371','g_Anaerococcus':'#87cefa','g_Gardnerella':'#20b2aa',
                     'Megasphaera_genomosp.':'#38b9ff','Streptococcus_agalactiae':'#ff738a','g_Megasphaera':'#008B45',
                     'Streptococcus_oralis':'#ff66ff','Prevotella_bivia':'#bfbfbf','Aerococcus_christensenii':'#bebebe',
                     'Anaerococcus_tetradius':'#87cefa','g_Gemella':'#daa520','Prevotella_genogroup_1':'#b0b0b0','Prevotella_disiens':'#9bbade',
                     'Lactobacillus_vaginalis':'#ffffff','other':'#808080','Bifidobacterium_longum':'#c1ffc1','Prevotella_buccalis':'#debdff',
                     'Bifidobacterium_breve':'#329ba8','Eggerthella':'#DE7710','Mycoplasma_hominis':'#10DE4E',
                     'Porphyromonas_uenonis':'#DE4310','Eubacterium_saphenum':'#8F10DE','Fusobacterium_nucleatum':'#CD853F',
                     'Fusobacterium_gonidiaformans':'#CD853F','Streptococcus_anginosus':'#ffc0cb','Peptostreptococcus_anaerobius':'#DEDB10',
                     'Arcanobacterium_phocae':'#8c10de','Bacteroides_uniformis':'#de1058','Ureaplasma_parvum':'#9999ff','Prevotella_corporis':'#ad8bab',
                     'Peptoniphilus_harei':'#CCCC00','Mobiluncus_mulieris':'#f08080','Megasphaera_sp._type_2':'#008B45',
                     'BVAB1':'#b31900','BVAB2':'#3bb16f','g_Escherichia.Shigella':'#12456b','Peptoniphilus_lacrimalis':'#803849',
                     'Veillonella_montpellierensis':'#ff8c69','Prevotella_genogroup_3':'#b36200','Parvimonas_micra':'#cdcd00',
                     'Corynebacterium_accolens':'#ffff00','Finegoldia_magna':'#800080','Prevotella_genogroup_2':'#b0b0b0','Propionibacterium':'#5f7362',
                     'Staphylococcus_epidermidis':'#ffffff','Prevotella_timonensis':'#f8fca9','g_Streptococcus':'#ffc0cb','g_Bifidobacterium':'#c1ffc1','g_Enterococcus':'#bbffff'
                     ,'g_Staphylococcus':'#800080','g_Finegoldia':'#800080','g_Prevotella':'#ffd9b3','g_Sneathia':'#d1e8eb','g_Aerococcus':'#e8e1ba'
                     ,'g_Leptotrichia':'#7ca386','g_Veillonella':'#499996','g_Dialister':'#997dbd','g_Corynebacterium_1':'#ffff00'
                     ,'g_Varibaculum':'#094717','g_Delftia':'#090c47','g_Corynebacterium_1':'#ffff00','Prevotella_amnii':'#ac8fc7','Sneathia_amnii':'#6963ff'}

def get_color(taxa):

    chars = '0123456789ABCDEF'

    if taxa in taxa_color_scheme:
        
        taxa_color = taxa_color_scheme[taxa]

    else:
        taxa_color_scheme[taxa] = '#'+''.join(random.sample(chars,6))
        taxa_color = taxa_color_scheme[taxa]

    return taxa_color


matplotlib.rc('font', serif='Helvetica')
matplotlib.rc('text', usetex=True)

data_all = pd.read_csv("MU01B_comp_glcorr_081521.csv")

data_all = data_all[data_all.columns[0:19]]

#subjects = list(set(data_all['SID']))
subjects = [1,3,5,7,9,2,4,6,8,10]

#setting up the plot to be a 6x8 matrix 
stacked_fig, stacked_axs = plt.subplots(2,5, figsize=(22,9), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.05,right=0.87,bottom=0.05,top=0.95,hspace=0.3, wspace=0.3)

stacked_axs = stacked_axs.ravel()


bar_width = 0.3
plot_count=0
legend_entries = {}

for subject in subjects:

    print(subject)
        
    #subsetting the data to just have this subject
    data_sub_rel = data_all.loc[data_all['Woman'] == subject]

    data_sub_micro = data_sub_rel[data_sub_rel.columns[4:]]

    #subsetting the dataframe to only include the columns with the top 15 summed relative abundance for the subject's data

    data_sub_micro['other'] = data_sub_micro.apply(lambda row: 1.0 - row.sum(), axis=1)

    data_sub_rel = pd.concat([data_sub_rel[data_sub_rel.columns[0:4]],data_sub_micro],axis=1,sort=False)

    data_sub_rel['bottom_count'] = pd.Series([0.0 for x in range(len(data_sub_rel.index))], index=data_sub_rel.index)

    for taxa in range(4,20):

        taxa = data_sub_rel.columns[taxa]

        taxa_color = get_color(taxa)

        if taxa not in legend_entries:

            legend_entries[taxa] = taxa_color_scheme[taxa]
        
        stacked_axs[plot_count].bar(data_sub_rel['Timepoint'],data_sub_rel[taxa],width=bar_width,bottom=data_sub_rel['bottom_count'],color=taxa_color,clip_on=False)
        
        data_sub_rel['bottom_count'] = data_sub_rel['bottom_count'] + data_sub_rel[taxa]

    stacked_axs[plot_count].set_title("W%s" %(subject),fontsize=24)
    stacked_axs[plot_count].set_yticks([0.0,0.2,0.4,0.6,0.8,1.0], minor=False)
    stacked_axs[plot_count].set_ylim(0,1)
    stacked_axs[plot_count].set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'],minor=False,fontsize=20)

    #if plot_count in [0,1,2,3,4]:
    #    stacked_axs[plot_count].set_xlabel('Relative abundance',fontsize=16)

    stacked_axs[plot_count].set_xticks(data_sub_rel.Timepoint.tolist())
    xlabs = ['T1','T2']
    #xlabs = [w.replace('.0','') for w in xlabs]
    stacked_axs[plot_count].set_xticklabels(xlabs,fontsize=24)

    #stacked_axs[plot_count].yaxis.grid(color='gray')
    stacked_axs[plot_count].set_axisbelow(True)
    plot_count += 1 

patch_list = []
for taxa in legend_entries:
    data_key = mpatches.Patch(color=legend_entries[taxa],label=taxa)
    patch_list.append(data_key)

stacked_fig.text(0.01,0.5,s='Relative abundance',fontsize=26,va='center',rotation='vertical',fontname='serif')
taxa_labels = [r'\textit{L. crispatus}',r'\textit{L. iners}',r'\textit{L. jensenii}',r'\textit{Gardnerella}',r'\textit{L. gasseri}',r'\textit{B. longum}',r'\textit{Ca.} L. vaginae',r'\textit{A. vaginae}',r'\textit{S. agalactiae}',r'\textit{P. buccalis}',r'\textit{Megasphaera}',r'\textit{Prevotella}',r'\textit{P. amnii}',r'\textit{S. sanguinegens}',r'\textit{P. harei}',r'Other']
stacked_fig.legend(handles=patch_list,labels=taxa_labels,loc=(0.875,0.2),ncol=1,fontsize=16,columnspacing=0.5)
stacked_fig.savefig('Figure1_MU01B_barcharts_all3.pdf')

