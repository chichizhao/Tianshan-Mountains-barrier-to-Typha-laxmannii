#!/bin/bash
# function: plot the frequency of each type of typha_lax

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import log
import numpy as np

# usage example: python plot_freq.py -1 north.txt -2 south.txt -o output.png
# or useage example: python plot_freq.py -i north_south.txt -o output.png

chr_len ={"chr1": 24094684, "chr2": 19094038, "chr3": 18910458, "chr4": 18888309, "chr5": 17861205,
          "chr6": 15065052, "chr7": 14536590, "chr8": 13397664, "chr9": 13040483,
          "chr10": 12597730, "chr11": 11722036, "chr12": 9996585, "chr13": 9906942,
          "chr14": 9656787, "chr15": 9218693}
chr_start_pos = {"chr1": 0, "chr2": 24094684, "chr3": 43088722, "chr4": 61999180, "chr5": 80887589,
                "chr6": 98748794, "chr7": 113013846, "chr8": 127550436, "chr9": 141948100,
                "chr10": 154988583, "chr11": 167586313, "chr12": 179308349, "chr13": 189304934,
                "chr14": 199211876, "chr15": 208868563}

"""
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot frequency of Typha Lax types by region.')
    parser.add_argument('-1', '--north', required=True, help='File containing North region data')
    parser.add_argument('-2', '--south', required=True, help='File containing South region data')
    parser.add_argument('-o', '--output', required=True, help='Output file for the plot')
    return parser.parse_args()
"""
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Plot frequency of Typha Lax types by region.')
    parser.add_argument('-i', '--north_south', required=True, help='File containing North and South region data')
    parser.add_argument('-o', '--output', required=True, help='Output file for the plot')
    return parser.parse_args()

def plot_freq1(north_file, south_file, output_file):
    # Load data
    north_data = pd.read_csv(north_file, sep='\t', header=None, names=['type', 'freq'])
    south_data = pd.read_csv(south_file, sep='\t', header=None, names=['type', 'freq'])

def plot_freq2():
    # Load data
    #data = pd.read_csv(north_south_file, sep='\t')
    #print(data.head())
    #data = data.head(1000)
    
    #  here we need to plot the frequency of 
    """
    the data should be like this:
        chromosome  position  N_heterozygous  S_heterozygous
    0       chr1     28122            0.02        0.076923
    1       chr1     28140            0.01        0.000000
    
    """
    # we should calculate the new position based on the chromosome and chr_start_pos position
    #data['new_position'] = data.apply(lambda row: row['position'] + chr_start_pos[row['chromosome']], axis=1)

    #print(data.head())
    # read the fst from the file
    fst_data = pd.read_csv('/home/chichi/data2/chichi/typha_lax/rasid/fst_10kb.windowed.weir.fst', sep='\t')
    """
    the data 
    CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
    chr1    1       500000  5894    0.0947266       0.0639279
    chr1    500001  1000000 7029    0.101329        0.0685913
    chr1    1000001 1500000 8035    0.114601        0.0707064
    """
    # data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    fst_data['new_position'] = fst_data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    # read the tajimaD from the file
    south_tajimaD = pd.read_csv('/home/chichi/data2/chichi/typha_lax/rasid/S_TajimaD_50kb.Tajima.D', sep='\t')
    north_tajimaD = pd.read_csv('/home/chichi/data2/chichi/typha_lax/rasid/N_tajima.Tajima.D', sep='\t')
    """
    CHROM	BIN_START	N_SNPS	TajimaD
    chr1	0	167	0.412421
    chr1	50000	438	1.20234
    chr1	100000	351	0.0222918
    """
    #data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    south_tajimaD['new_position'] = south_tajimaD.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    north_tajimaD['new_position'] = north_tajimaD.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    # read the pi from the file
    south_pi = pd.read_csv('/home/chichi/data2/chichi/typha_lax/rasid/S_pi.windowed.pi', sep='\t')
    north_pi = pd.read_csv('/home/chichi/data2/chichi/typha_lax/rasid/N_pi.windowed.pi', sep='\t')
    """
    CHROM	BIN_START	BIN_END	N_VARIANTS	PI
    chr1	1	50000	192	0.000826651
    chr1	50001	100000	543	0.0027542
    chr1	100001	150000	415	0.00159258
    chr1	150001	200000	434	0.0019616
    chr1	200001	250000	639	0.00334921
    """
    #data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    # here we need facilitate the plots 
    # we select the first 1000 row for adjust the plot
    # data = data.head(1000)
    south_pi['new_position'] = south_pi.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    north_pi['new_position'] = north_pi.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)

    fig, ax = plt.subplots(figsize=(8, 5))
    #sns.scatterplot(data=data, x='new_position', y='N_heterozygous', label='North', facecolors='blue', edgecolor='none', alpha=0.5, ax=ax)
    # modified the S_heterozygous to be negitive for better visualization
    #sns.scatterplot(data=data, x='new_position', y=-data['S_heterozygous'], label='South', facecolors='red', edgecolor='none', alpha=0.5, ax=ax)
    # may we can use the lineplot to show the frequency
    #lineplot(data=data, x='new_position', y='N_heterozygous', label='North', color='blue', ax=ax, alpha=0.2)
    #sns.lineplot(data=data, x='new_position', y=data['S_heterozygous'], label='South', color='red', ax=ax, alpha=0.2)
    #sns.lineplot(data=data, x='new_position', y='MEAN_FST', label='Mean FST', color='blue', ax=ax, alpha=0.8)
    # plot with the chromesomes

    # plot the manhattan plot for var of each chr, Fst
    ax.set_xlim(-3500000, 208868563+9218693)
    ax.set_ylim(-0.2, 7.8)
    box = plt.Rectangle((0, -8), 208868563+9218693, 2, fc='#34587F', ec='#34587F', lw=1)
    ax.add_patch(box)
    num = 0
    ax.axis('off')
    chr_number = 0
    data = fst_data
    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        x = []
        y = []
        x_significant = []
        y_significant = []
        if chr_number % 2 == 0:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['MEAN_FST'].iloc[i]/0.6*1.6)
                if chr_data['MEAN_FST'].iloc[i] > 0.2:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['MEAN_FST'].iloc[i]/0.6*1.6)
            ax.scatter(x, y, s=0.5, label=f'{chr}', color='green', alpha=0.2)
            ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.5)

        else:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['MEAN_FST'].iloc[i]/0.6*1.6)
                if chr_data['MEAN_FST'].iloc[i] > 0.2:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['MEAN_FST'].iloc[i]/0.6*1.6)

            ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.5)
        
        chr_number += 1
    

    for i in range(1, 16):
        if i < 10:
            chr = 'chr' + str(i)
        else:
            chr = 'chr' + str(i)
        ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -0.06/0.6*2, str(i), ha='center', va='center', fontsize=8, color='black')
        # add the chromosome region line 
        #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
        if i % 2 == 1:
            # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
            box = plt.Rectangle((chr_start_pos[chr], -0.06), chr_len[chr], 8, fc='grey', ec='none', lw=0.5, alpha=0.2)
            ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
    #ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
    #ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')
    box = plt.Rectangle((-3500000, 0), 1, 120, fc='black', ec='black', lw=1)
    ax.add_patch(box)

    box = plt.Rectangle((-3500000, 0), 1, 120, fc='black', ec='black', lw=1)
    ax.add_patch(box)
    #ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)
    ax.axhline(y=0.2/0.6*1.8, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    ax.text(2000000, 0.2/0.6*1.6, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-12500000, 0.3/0.6*1.6, " Fst \n (North vs South)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    ax.text(-1000000,  -0.06/0.6*2, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 0, "0 ", ha='right', va='center', fontsize=7, color='black')
    box = plt.Rectangle((-3500000, 0), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 0.2/0.6*1.6, "0.2", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 0.2/0.6*1.6), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 0.6/0.6*1.6, "0.6", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 0.6/0.6*1.6), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 0.4/0.6*1.6, "0.4", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 0.4/0.6*1.6), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    #ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)
    box = plt.Rectangle((0, -0.02/0.6*1.6), 208868563+9218693, 0.02, fc='grey', ec='none', lw=0.2, alpha=0.7)
    ax.add_patch(box)

    ax.set_xlabel('Position on Chromosome')
    ax.set_ylabel('Frequency')
    #ax.set_title('Frequency of Typha Lax Types by Region')
    #ax.legend()
    #plt.xticks(rotation=45)
    plt.tight_layout()
    # set the x-axis limits based on chromosome lengths
    #ax.set_xlim(0, sum(chr_len.values()))
    # set the y-axis limits
    #ax.set_ylim(-1, 1)
    # add the margin ot the Fst plot
    box = plt.Rectangle((-3500000, 0.7/0.6*1.6), 208868563+9218693+4000000, 0.1, fc='white', ec='none', lw=0.2)
    ax.add_patch(box)

    
    data = south_tajimaD
    # we need adjust the TajimaD value to locate the plots of TajimaD

    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        x = []
        y = []
        x_significant = []
        y_significant = []
        if chr_number % 2 == 0:
            for i in range(len(chr_data)):
                # here we use the bin to plot TajimaD
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 2.5), 50000, chr_data['TajimaD'].iloc[i]/4, fc='green', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)

        else:
            for i in range(len(chr_data)):
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 2.5), 50000, chr_data['TajimaD'].iloc[i]/4, fc='blue', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)


        chr_number += 1
    

    for i in range(1, 16):
        if i < 10:
            chr = 'chr' + str(i)
        else:
            chr = 'chr' + str(i)
        #ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -0.06, str(i), ha='center', va='center', fontsize=8, color='black')
        # add the chromosome region line 
        #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
        #if i % 2 == 1:
            # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
            #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 1, fc='grey', ec='none', lw=0.5, alpha=0.2)
            #ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
    #ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
    #ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')

    ax.axhline(y=2.5, color='black', alpha=1, linewidth=0.4)
    #ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    #ax.text(2000000, 0.2, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-9500000, 2.5, "Sorth", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    #ax.text(-1000000,  -0.06, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 2.5, "0", ha='right', va='center', fontsize=7, color='black')
    box = plt.Rectangle((-3500000, 0/3+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 1/4+2.5, "1", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 1/4+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 2/4+2.5, "2", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 2/4+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 3/4+2.5, "3", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 3/4+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, -2/4+2.5, "-2", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, -2/4+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, -1/4+2.5, "-1", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, -1/4+2.5), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    # add the margin to the TajimaD plot
    box = plt.Rectangle((-3500000,3.38), 208868563+9218693+4000000, 0.1, fc='white', ec='none', lw=0.2)
    ax.add_patch(box)
    box = plt.Rectangle((-3500000,4.88), 208868563+9218693+4000000, 0.62, fc='white', ec='none', lw=0.2)
    ax.add_patch(box)
    data = north_tajimaD

    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        x = []
        y = []
        x_significant = []
        y_significant = []
        if chr_number % 2 == 0:
            for i in range(len(chr_data)):
                # here we use the bin to plot TajimaD
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 4), 50000, chr_data['TajimaD'].iloc[i]/4, fc='green', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)

        else:
            for i in range(len(chr_data)):
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 4), 50000, chr_data['TajimaD'].iloc[i]/4, fc='blue', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)


        chr_number += 1
    

    for i in range(1, 16):
        if i < 10:
            chr = 'chr' + str(i)
        else:
            chr = 'chr' + str(i)
        #ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -0.06, str(i), ha='center', va='center', fontsize=8, color='black')
        # add the chromosome region line 
        #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
        #if i % 2 == 1:
            # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
            #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 1, fc='grey', ec='none', lw=0.5, alpha=0.2)
            #ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
    #ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
    #ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')

    ax.axhline(y=2.5, color='black', alpha=1, linewidth=0.4)
    #ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    #ax.text(2000000, 0.2, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-9500000, 4, "North", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    #ax.text(-1000000,  -0.06, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 4, "0", ha='right', va='center', fontsize=7, color='black')
    box = plt.Rectangle((-3500000, 0/3+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 1/4+4, "1", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 1/4+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 2/4+4, "2", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 2/4+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, 3/4+4, "3", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 3/4+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, -2/4+4, "-2", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, -2/4+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)
    ax.text(-3700000, -1/4+4, "-1", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, -1/4+4), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    ax.add_patch(box)

    # add the TajimaD label to the margin plot
    ax.text(-13500000, 3.4, " Tajima's D", ha='center', va='center', fontsize=7, color='black', rotation= 90)

    # plot the pi for each chr
    data = south_pi
    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        x = []
        y = []
        x_significant = []
        y_significant = []
        if chr_number % 2 == 0:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['PI'].iloc[i])
                # here we use the bin to plot PI
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 5.5), 50000, chr_data['PI'].iloc[i]/0.01, fc='green', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)
                if chr_data['PI'].iloc[i] > 0:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['PI'].iloc[i])
            #ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            #ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.2)
            # here we plot the line with green color
            #x.plot(x, y, label=f'{chr}', color='green', alpha=0.6, linewidth=0.4)
            

        else:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['PI'].iloc[i])
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 5.5), 50000, chr_data['PI'].iloc[i]/0.01, fc='blue', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)
                if chr_data['PI'].iloc[i] > 0:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['PI'].iloc[i])

            #ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            #ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.2)
            # here we plot the line with red color
            #ax.plot(x, y, label=f'{chr}', color='blue', alpha=0.6, linewidth=0.4)


        chr_number += 1
    


    #box = plt.Rectangle((-3500000, -5), 1, 120, fc='black', ec='black', lw=1)
    #ax.add_patch(box)
    ax.axhline(y=5.5, color='black', alpha=1, linewidth=0.4)
    #ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    #ax.text(2000000, 0.2, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-7500000, 0.005/0.01*1+5.5, "South", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    #ax.text(-1000000,  -0.06, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 0/0.01*1+5.5, "0", ha='right', va='center', fontsize=7, color='black')
    #box = plt.Rectangle((-3500000, 0), 10000000, 0.01, fc='black', ec='none', lw=0.4)
    #ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 0.01/0.01*1+5.5, "0.01", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 1+5.5), 1000000, 0.01, fc='black', ec='none', lw=1)
    ax.add_patch(box)

    # add the margin to the pi plot
    box = plt.Rectangle((-3500000, 6.60), 208868563+9218693+4000000, 0.15, fc='white', ec='none', lw=0.2)
    ax.add_patch(box)

    data = north_pi
    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        x = []
        y = []
        x_significant = []
        y_significant = []
        if chr_number % 2 == 0:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['PI'].iloc[i])
                # here we use the bin to plot PI
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 6.75), 50000, chr_data['PI'].iloc[i]/0.01, fc='green', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)
                if chr_data['PI'].iloc[i] > 0:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['PI'].iloc[i])
            #ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            #ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.2)
            # here we plot the line with green color
            #x.plot(x, y, label=f'{chr}', color='green', alpha=0.6, linewidth=0.4)
            

        else:
            for i in range(len(chr_data)):
                x.append(chr_data['new_position'].iloc[i])
                y.append(chr_data['PI'].iloc[i])
                box = plt.Rectangle((chr_data['new_position'].iloc[i], 6.75), 50000, chr_data['PI'].iloc[i]/0.01, fc='blue', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)
                if chr_data['PI'].iloc[i] > 0:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['PI'].iloc[i])

            #ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            #ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.2)
            # here we plot the line with red color
            #ax.plot(x, y, label=f'{chr}', color='blue', alpha=0.6, linewidth=0.4)


        chr_number += 1
    


    #box = plt.Rectangle((-3500000, -5), 1, 120, fc='black', ec='black', lw=1)
    #ax.add_patch(box)
    ax.axhline(y=6.75, color='black', alpha=1, linewidth=0.4)
    #ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    #ax.text(2000000, 0.2, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-7500000, 0.005/0.01*1+6.75, "North", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    #ax.text(-1000000,  -0.06, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 0/0.01*1+6.75, "0", ha='right', va='center', fontsize=7, color='black')
    #box = plt.Rectangle((-3500000, 0), 10000000, 0.01, fc='black', ec='none', lw=0.4)
    #ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 0.01/0.01*1+6.75, "0.01", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 1+6.75), 1000000, 0.01, fc='black', ec='none', lw=1)
    ax.add_patch(box)

    ax.text(-13500000, 6.8, "π", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    


    plt.savefig("vcf_result.png", dpi=300, bbox_inches='tight')


    plt.close()

if __name__ == "__main__":
    #args = get_args()
    #plot_freq2(args.north_south, args.output)
    plot_freq2()
