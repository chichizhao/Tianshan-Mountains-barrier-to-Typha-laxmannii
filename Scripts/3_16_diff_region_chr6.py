
#!/bin/bash
# function: plot the frequency of each type of typha_lax

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

def plot_freq2(north_south_file, output_file):
    # Load data
    data = pd.read_csv(north_south_file, sep='\t')
    print(data.head())
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

    """
    the data 
    CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
    chr1    1       500000  5894    0.0947266       0.0639279
    chr1    500001  1000000 7029    0.101329        0.0685913
    chr1    1000001 1500000 8035    0.114601        0.0707064
    """
    # data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)

    """
    CHROM	BIN_START	N_SNPS	TajimaD
    chr1	0	167	0.412421
    chr1	50000	438	1.20234
    chr1	100000	351	0.0222918
    """
    #data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    """
    CHROM	BIN_START	BIN_END	N_VARIANTS	PI
    chr1	1	50000	192	0.000826651
    chr1	50001	100000	543	0.0027542
    chr1	100001	150000	415	0.00159258
    chr1	150001	200000	434	0.0019616
    chr1	200001	250000	639	0.00334921
    """
    data['new_position'] = data.apply(lambda row: row['BIN_START'] + chr_start_pos[row['CHROM']], axis=1)
    # here we need facilitate the plots 
    # we select the first 1000 row for adjust the plot
    #data = data.head(1000)

    fig, ax = plt.subplots(figsize=(5, 1))
    #sns.scatterplot(data=data, x='new_position', y='N_heterozygous', label='North', facecolors='blue', edgecolor='none', alpha=0.5, ax=ax)
    # modified the S_heterozygous to be negitive for better visualization
    #sns.scatterplot(data=data, x='new_position', y=-data['S_heterozygous'], label='South', facecolors='red', edgecolor='none', alpha=0.5, ax=ax)
    # may we can use the lineplot to show the frequency
    #lineplot(data=data, x='new_position', y='N_heterozygous', label='North', color='blue', ax=ax, alpha=0.2)
    #sns.lineplot(data=data, x='new_position', y=data['S_heterozygous'], label='South', color='red', ax=ax, alpha=0.2)
    #sns.lineplot(data=data, x='new_position', y='TajimaD', label='Mean FST', color='blue', ax=ax, alpha=0.8)
    # plot with the chromesomes

    # plot the manhattan plot for var of each chr
    ax.set_xlim(-3500000, 16065052)
    ax.set_ylim(-0.0001,0.01)
    #box = plt.Rectangle((0, -8), 16065052, 2, fc='#34587F', ec='#34587F', lw=1)
    #ax.add_patch(box)
    num = 0
    ax.axis('off')
    chr_number = 0
    for chr, length in chr_len.items():
        chr = 'chr6'
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
                box = plt.Rectangle((chr_data['new_position'].iloc[i]-chr_start_pos[chr], 0), 50000, chr_data['PI'].iloc[i], fc='green', ec='none', lw=0.5, alpha=0.5)
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
                box = plt.Rectangle((chr_data['new_position'].iloc[i]-chr_start_pos[chr], 0), 50000, chr_data['PI'].iloc[i], fc='blue', ec='none', lw=0.5, alpha=0.5)
                ax.add_patch(box)
                if chr_data['PI'].iloc[i] > 0:
                    x_significant.append(chr_data['new_position'].iloc[i])
                    y_significant.append(chr_data['PI'].iloc[i])

            #ax.scatter(x, y, s=0.5, label=f'{chr}', color='blue', alpha=0.2)
            #ax.scatter(x_significant, y_significant, s=1, label=f'{chr} significant', color='red', alpha=0.2)
            # here we plot the line with red color
            #ax.plot(x, y, label=f'{chr}', color='blue', alpha=0.6, linewidth=0.4)


        chr_number += 1
    

    #for i in range(1, 16):
    #    if i < 10:
    #        chr = 'chr' + str(i)
    #    else:
    #        chr = 'chr' + str(i)
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
    box = plt.Rectangle((-3500000, 0), 1, 120, fc='black', ec='black', lw=0.4)
    ax.add_patch(box)

    #box = plt.Rectangle((-3500000, -5), 1, 120, fc='black', ec='black', lw=1)
    #ax.add_patch(box)
    ax.axhline(y=0, color='black', alpha=1, linewidth=0.4)
    #ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=0.4)
    # add the y axis label
    #ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
    #ax.text(2000000, 0.2, "0.2", ha='right', va='bottom', fontsize=8, color='red')
    ax.text(-5000000, 0.005, "π", ha='center', va='center', fontsize=7, color='black', rotation= 90)
    #ax.text(-1000000,  -0.06, "Chr", ha='center', va='center', fontsize=7, color='black')
    ax.text(-3700000, 0, "0", ha='right', va='center', fontsize=7, color='black')
    #box = plt.Rectangle((-3500000, 0), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    #ax.add_patch(box)
    #ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
    ax.text(-3700000, 0.01, "0.01", ha='right', va='center', fontsize=7)
    box = plt.Rectangle((-3500000, 1), 100000, 0.01, fc='black', ec='none', lw=1)
    ax.add_patch(box)
    #ax.text(-3700000, 0.005, "π", ha='right', va='center', fontsize=7)
    #box = plt.Rectangle((-3500000, 2), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    #ax.add_patch(box)
    #ax.text(-3700000, 0.03, "0.03", ha='right', va='center', fontsize=7)
    #box = plt.Rectangle((-3500000, 0.03), 1000000, 0.01, fc='black', ec='none', lw=0.4)
    #ax.add_patch(box)

    #ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)
    #box = plt.Rectangle((0, -0.02), 208868563+9218693, 0.02, fc='grey', ec='none', lw=0.2, alpha=0.7)
    #ax.add_patch(box)


    """
    
    chr_number = 0
    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        if chr_number % 2 == 0:
            sns.lineplot(data=chr_data, x='new_position', y='PI', label=f'{chr}', color='blue', ax=ax, alpha=0.5)
        else:
            sns.lineplot(data=chr_data, x='new_position', y='Tajima', label=f'{chr}', color='red', ax=ax, alpha=0.5)
    """
    """
    chr_number = 0
    for chr, length in chr_len.items():
        chr_data = data[data['CHROM'] == chr]
        if chr_number % 2 == 0:
            sns.scatterplot(data=chr_data, x='new_position', y='PI', label=f'{chr}', facecolors='blue', edgecolor='none', alpha=0.5, ax=ax)
        else:
            sns.scatterplot(data=chr_data, x='new_position', y='PI', label=f'{chr}', facecolors='red', edgecolor='none', alpha=0.5, ax=ax)
        chr_number += 1
    """
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

    plt.savefig(output_file, dpi=300, bbox_inches='tight')


    plt.close()

if __name__ == "__main__":
    args = get_args()
    plot_freq2(args.north_south, args.output)

# the core division region is 5550000-8550000 

# we first read the postive selectio signal on chr5
# gff file
gff="/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.gene.gff"
# annotation file
anno="/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"

# read the gff file to get the gene on chr5
import pandas as pd
gff_file = pd.read_csv(gff, sep="\t", header=None, comment="#")
# filter the gff file to get the gene on chr5
chr6_genes = gff_file[(gff_file[0] == "chr6") & (gff_file[2] == "mRNA")]
# get the gene names
chr6_gene_names = chr6_genes[8].str.split("ID=").str[1].str.split(";").str[0]
# get the gene start and end positions
chr6_gene_starts = chr6_genes[3]
chr6_gene_ends = chr6_genes[4]
# create a dataframe with the gene names, starts and ends
chr6_genes_df = pd.DataFrame({
    "gene_name": chr6_gene_names,
    "start": chr6_gene_starts,
    "end": chr6_gene_ends
})
# filter the gene in the sweep region
sweep_start = 5550000
sweep_end = 8550000
sweep_genes = chr6_genes_df[(chr6_genes_df["start"] >= sweep_start) & (chr6_genes_df["end"] <= sweep_end)]

print("Significant genes on chr6 in the sweep region:")
print(sweep_genes)

# read the annotation file to get the gene functions
anno_file = pd.read_csv(anno, sep="\t", header=0)
# filter the annotation file to get the gene functions for the significant genes
# write the annotation into a file
with open("chr6_significant_genes_annotation.txt", "w") as f:
    for gene in sweep_genes["gene_name"]:
        gene_info = anno_file[anno_file["SeqID"] == gene]
        if not gene_info.empty:
            gene_info.to_csv(f, sep="\t", header=False, index=False)
        else:
            f.write(f"{gene}\tNo annotation found\n")
