#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import matplotlib.patches as mpatches

data = pd.read_csv("MU01B_popANI.csv",sep=",")

stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(4.5,6), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.14,right=0.95,bottom=0.1,top=0.9,hspace=0.3, wspace=0.3)

taxa_palette = {'Lcrispatus':'#ff0000','Gardnerella':'#20b2aa','Liners':'#ff8c00',
				'Lgasseri':'#7fff00','Ljensenii':'#333333','Blongum':'#c1ffc1',
				'BVAB1':'#b31900','Prevotella':'#70005e','Avaginae':'#0000cd','Pbuccalis':'#debdff',
				'Megasphaera':'#38b9ff'}
sns.histplot(data=data, x="popANI_plotting",bins=[0.9978,0.998,0.9982,0.9984,0.9986,0.9988,0.9990,0.9992,0.9994,0.9996,0.9998,1.0],
			 hue='Species',palette=taxa_palette,axes=stacked_axs,legend=False,multiple="stack")

stacked_axs.set_ylim(0,17)
stacked_axs.set_xlabel('Population ANI',fontsize=16)
stacked_axs.set_ylabel('Count',fontsize=16)
stacked_axs.set_yticks([0,4,8,12,16])
stacked_axs.set_yticklabels(['0','4','8','12','16'],fontsize=14)
stacked_axs.set_xticks([0.9979,0.9985,0.9990,0.9995,1.000])
stacked_axs.set_xticklabels(['<99.7%','99.85%','99.9%','99.95%','100%'],fontsize=14)
stacked_axs.axvline(x=0.9990,linewidth=4,c='#000000',linestyle='dotted')

stacked_fig.savefig('MU01B_popANI.pdf')
