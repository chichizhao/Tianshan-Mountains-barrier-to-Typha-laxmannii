#!/bin/python
# function: plot the manhanttan plot of the gwas data and annotate the genes
# usage: python plot_sig_and_annotegene.py
#

#!/bin/python
# function: plot the manhanttan plot of the gwas data
# usage: python plot_snp.py

# step 1: read the gwas data with p value
# file: sig_SNPs_bio_5_mod_sub_env_data4.part5_merge3_filter_variants.csv
# read the data with pandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

chr_len ={"chr1": 24094684, "chr2": 19094038, "chr3": 18910458, "chr4": 18888309, "chr5": 17861205,
          "chr6": 15065052, "chr7": 14536590, "chr8": 13397664, "chr9": 13040483,
          "chr10": 12597730, "chr11": 11722036, "chr12": 9996585, "chr13": 9906942,
          "chr14": 9656787, "chr15": 9218693}
chr_start_pos = {"chr1": 0, "chr2": 24094684, "chr3": 43088722, "chr4": 61999180, "chr5": 80887589,
                "chr6": 98748794, "chr7": 113013846, "chr8": 127550436, "chr9": 141948100,
                "chr10": 154988583, "chr11": 167586313, "chr12": 179308349, "chr13": 189304934,
                "chr14": 199211876, "chr15": 208868563}


file = '/home/chichi/data2/chichi/typha_lax/gwas/pos_mod_sub_N_S_126_raw_filter_variants_maf0.01_miss0.9.assoc.txt'
data = pd.read_csv(file, sep='\t',index_col= False)
#print(data.p_value)
print(data.head())
# remove the line for the chr key has the Contig1 and Contig2
# Remove the rows where 'chr' contains 'Contig' or is NaN
#data = data[data['chr'].notna() & ~data['chr'].str.contains('Contig')]

# step 2: plot the manhattan plot
fig = plt.figure(figsize=(16, 2.5))
ax = fig.add_subplot(111)
# off the axis
ax.axis('off')
ax.set_xlim(0, 208868563+9218693+5000000)
ax.set_ylim(0, 55)
num = 0
#print(len(data))
#print(len(data))
# remove the data with p_wald < 1e-4
# Create a new column 'chr_str' with formatted chromosome numbers
data['chr_str'] = 'chr' + data['chr'].apply(lambda x: f'{x:01}')
print('new column chr_str created')
# Create a new column 'color' based on the 'p_wald' values
data['color'] = np.where(data['p_wald'] >= 7.46e-08, 'grey', 'blue')

# Create a new column 'x' for the x-coordinates of the scatter plot
data['x'] = data.apply(lambda row: chr_start_pos[row['chr_str']] + row['ps'], axis=1)

print('new column x created')
# Create a new column 'y' for the y-coordinates of the scatter plot
data['y'] = -np.log10(data['p_wald'])
print('new column y created')

# Plot all rows
#for _, row in data.iterrows():
#    ax.scatter(row['x'], row['y'], color=row['color'], s=2,alpha=0.5,rasterized=True)

# plot the data with with the chromosomes 
for i in range(1, 16):
    if i < 10:
        chr = 'chr' + str(i)
    else:
        chr = 'chr' + str(i)
    chr_data_x = data[data['chr_str'] == chr]['x']
    chr_data_y = data[data['chr_str'] == chr]['y']
    chr_data_color = data[data['chr_str'] == chr]['color']
    if i%2 == 0:
        # replace the grey color with black
        chr_data_color = chr_data_color.replace('grey', 'black')
    ax.scatter(chr_data_x, chr_data_y, color=chr_data_color, s=2,alpha=0.5,rasterized=True)
#for _, row in data.iloc[:5000].iterrows():
#    ax.scatter(row['x'], row['y'], color=row['color'], s=1)
print('plot the manhattan plot')

# annotate the genes 
## read the gene annotation file
gff ="/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.gene.gff.cds"
gff_data = pd.read_csv(gff, sep='\t', header=None, comment='#')
gff_data.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'info']
#print(gff_data.head())

## read the gene annotation file 
annotation = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"
annotation_data = pd.read_csv(annotation, sep='\t',index_col= False)
annotation_dict = {}
for i in range(len(annotation_data)):
    if str(annotation_data.loc[i, 'Pfam']).startswith('PF'):
        # SeqID	KEGG	Pathway	Nr	Uniprot	GO	KOG	Pfam
        annotation_dict[annotation_data.loc[i, 'SeqID']] = [annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'Pathway'], annotation_data.loc[i, 'Nr'], annotation_data.loc[i, 'Uniprot'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KOG'], annotation_data.loc[i, 'Pfam'], annotation_data.loc[i, 'Interpro']]
        #annotation_dict[annotation_data.loc[i, 'SeqID']] = [annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'Pathway'], annotation_data.loc[i, 'Nr'], annotation_data.loc[i, 'Uniprot'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KOG'], annotation_data.loc[i, 'Pfam'], annotation_data.loc[i, 'Interpro']]
    annotation_dict[annotation_data.loc[i, 'SeqID']] = [annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'Pathway'], annotation_data.loc[i, 'Nr'], annotation_data.loc[i, 'Uniprot'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KOG'], annotation_data.loc[i, 'Pfam'], annotation_data.loc[i, 'Interpro']]
    if  annotation_data.loc[i, 'SeqID'] == 'TlaxChr11G00166340.1':
        print(annotation_data.loc[i, 'SeqID'], annotation_data.loc[i, 'KEGG'], annotation_data.loc[i, 'Pathway'], annotation_data.loc[i, 'Nr'], annotation_data.loc[i, 'Uniprot'], annotation_data.loc[i, 'GO'], annotation_data.loc[i, 'KOG'], annotation_data.loc[i, 'Pfam'], annotation_data.loc[i, 'Interpro'])
# print(annotation_dict.keys())

file = '/home/chichi/data2/chichi/typha_lax/gwas/p_wald_pos_mod_sub_N_S_126_raw_filter_variants_maf0.01_miss0.9_top0.01.csv'
data = pd.read_csv(file, sep=',',index_col= False)
# keep the data with the first 100 rows
# data = data.iloc[:100]

print('gene annotation file read')
gene_list_info = []


# here we need check the snp in CDS region or not
for i in range(len(data)):
    # print the processing infomation
    if i % 100 == 0:
        print(f'Processing {i}th SNP')
    chr = data.loc[i, 'chr']
    pos = data.loc[i, 'ps']
    if chr < 10:
        chr = 'chr' + str(chr)
    else:
        chr = 'chr' + str(chr)


    #print(chr)
    # check the snp in the CDS region or not
    print(gff_data.columns)
    gff_chr = gff_data[gff_data['chr'] == chr]
    #print(gff_chr)
    if -np.log10(data.loc[i, 'p_wald']) >= 7.46:
        for j in range(len(gff_chr)):
            if pos >= gff_chr.iloc[j, 3] and pos <= gff_chr.iloc[j, 4]:
                #ax.text(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald']), gff_chr.iloc[j, 8], ha='center', va='center', fontsize=8, color='black')
                id = gff_chr.iloc[j, 8].split('Parent=')[1]
                id = str(id)
                if id in annotation_dict:
                    # ax.text(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald'])+1, annotation_dict[id], ha='center', va='center', fontsize=12, color='black')
                    # add the gene information to the gene_list_info
                    gene_list_info.append([chr, pos, str(id), -np.log10(data.loc[i, 'p_wald']), annotation_dict[id][0], annotation_dict[id][1], annotation_dict[id][2], annotation_dict[id][3], annotation_dict[id][4], annotation_dict[id][5], annotation_dict[id][6]])
                    # highlight the gene point
                    ax.scatter(chr_start_pos[chr] + pos, -np.log10(data.loc[i, 'p_wald']), color='red', s=8,rasterized=True)
                break

#sort the gene_list_info based on the chromosome and position
#gene_list_info = sorted(gene_list_info, key=lambda x: (x[0], x[1]))
# sort the gene_list_info based on the p_wald value
gene_list_info = sorted(gene_list_info, key=lambda x: (x[3], x[0], x[1]), reverse=True)


with open('bio16_top0.01_annotated_genes.xls', 'w') as f:  
    #f.write("chr\tpos\tgene_id\tPFAM_family\tPFAM\tGO\tKEGG\tNR\tSwissprot\n")
    #f.write("chr,pos,gene_id,PFAM_family,PFAM,GO,KEGG,NR,Swissprot\n")
    f.write("chr\tpos\tgene_id\tp_wald\tKEGG\tPathway\tNr\tUniprot\tGO\tKOG\tPfam\tInterpro\n")
   
    for item in gene_list_info:
        item = item[:-1]
        #f.write("%s\n" % str(item).replace('[','').replace(']','').replace('\'',''))
        #f.write("%s\n" % '\t'.join(map(str, item)))
        f.write("%s\n" % '\t'.join(map(str, item)))

#print( gene_list_info)


# remove the snp on the same gene
new_gene_list_info = []
for i in range(len(gene_list_info)):
    if i == 0:
        new_gene_list_info.append(gene_list_info[i])
    else:
        if gene_list_info[i][2] != gene_list_info[i-1][2]:
            new_gene_list_info.append(gene_list_info[i])


# add the line at the significant threshold of 7.46e-08
#ax.axhline(y=-np.log10(7.46e-08), color='b', linestyle='--', alpha=0.5)
plt.plot([0, 216599865], [-np.log10(7.46e-08), -np.log10(7.46e-08)], color='blue', linestyle='--', alpha=0.5)
# add the mark of the 0
#ax.axhline(y=0, color='black', linestyle='-')
plt.plot([-4000000, 223599865], [-0.1, -0.1], color='black', linestyle='-')
box = plt.Rectangle((-4000000, 0), 100000, 20, fc='none', ec='black', lw=1)
ax.add_patch(box)


# add the mark of the chromosome
for i in range(1, 16):
    if i < 10:
        chr = 'chr' + str(i)
    else:
        chr = 'chr' + str(i)
    ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -3.4, str(i), ha='center', va='center', fontsize=14, color='black')
    # add the chromosome region line 
    #if i%2 == 1:
    #    box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 20, fc='grey', ec='none', lw=1, alpha=0.5)
    #    ax.add_patch(box)
#ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
ax.axvline(x=chr_start_pos["chr1"], color='black', linestyle='-')
# add the y axis label
#ax.text(-0.1, 7.46, '-log10(7.46e-08)', ha='center', va='center', fontsize=12, color='black')
ax.text(-13000000, 25, '-log10(p_wald)', ha='center', va='center', fontsize=14, color='black', rotation=90) 
ax.text(-3700000, 0, '0', ha='right', va='center', fontsize=14, color='black')
ax.text(-3700000, 10, '10', ha='right', va='center', fontsize=14, color='black')
#ax.text(-3000000, 7.46, '7.46', ha='center', va='center', fontsize=12, color='blue')
ax.text(-3700000, 20, '20', ha='right', va='center', fontsize=14, color='black')
ax.text(-3700000, 30, '30', ha='right', va='center', fontsize=14, color='black')
ax.text(-3700000, 40, '40', ha='right', va='center', fontsize=14, color='black')
ax.text(-3700000, 50, '50', ha='right', va='center', fontsize=14, color='black')
ax.text(0, -3.4, 'chr', ha='center', va='center', fontsize=14, color='black')
#box = plt.Rectangle((0, 0), 216599865, 20, fc='none', ec='black', lw=3)
#ax.add_patch(box)
box = plt.Rectangle((208868563+9218693, 0), 5000000, 55, fc='grey', ec='black', lw=1, alpha=0.1)
ax.add_patch(box)
# add the bio2 label
ax.text(216599865+4500000, 25, 'North vs South', ha='center', va='center', fontsize=14, color='black',rotation=90)
# add the x-axis and y-axis label
# ax.set_xlabel('Position')
# ax.set_ylabel('-log10(p-value)')
# ax.set_title('Manhattan plot of GWAS data')

# here we want add the top 10 genes to the manhattan plot
for i in range(len(new_gene_list_info)):
    # add the gene id to the manhattan plot
    if i <5:
        annotation_pos = new_gene_list_info[i][1]
        annotation_chr = new_gene_list_info[i][0]
        annotation_id = new_gene_list_info[i][2].replace('Tlax','')
        annotation_p_wald = new_gene_list_info[i][3]
        ax.text(chr_start_pos[annotation_chr] + annotation_pos, annotation_p_wald+5, annotation_id, ha='center', va='center', fontsize=12, color='black')
        # add the connection line
        plt.plot([chr_start_pos[annotation_chr] + annotation_pos, chr_start_pos[annotation_chr] + annotation_pos], [annotation_p_wald+5, annotation_p_wald], color='black', linewidth=0.5, alpha=0.5)    
    elif i == 5:
        annotation_pos = new_gene_list_info[i][1]
        annotation_chr = new_gene_list_info[i][0]
        annotation_id = new_gene_list_info[i][2].replace('Tlax','')
        annotation_p_wald = new_gene_list_info[i][3]
        ax.text(chr_start_pos[annotation_chr] + annotation_pos, annotation_p_wald+20, annotation_id, ha='center', va='center', fontsize=12, color='black')
        # add the connection line
        plt.plot([chr_start_pos[annotation_chr] + annotation_pos, chr_start_pos[annotation_chr] + annotation_pos], [annotation_p_wald+20, annotation_p_wald], color='black', linewidth=0.5, alpha=0.5)
    elif i == 7:
        annotation_pos = new_gene_list_info[i][1]
        annotation_chr = new_gene_list_info[i][0]
        annotation_id = new_gene_list_info[i][2].replace('Tlax','')
        annotation_p_wald = new_gene_list_info[i][3]
        ax.text(chr_start_pos[annotation_chr] + annotation_pos, annotation_p_wald+12, annotation_id, ha='center', va='center', fontsize=12, color='black')
        # add the connection line
        plt.plot([chr_start_pos[annotation_chr] + annotation_pos, chr_start_pos[annotation_chr] + annotation_pos], [annotation_p_wald+12, annotation_p_wald], color='black', linewidth=0.5, alpha=0.5)
    elif i == 6:
        annotation_pos = new_gene_list_info[i][1]
        annotation_chr = new_gene_list_info[i][0]
        annotation_id = new_gene_list_info[i][2].replace('Tlax','')
        annotation_p_wald = new_gene_list_info[i][3]
        ax.text(chr_start_pos[annotation_chr] + annotation_pos, annotation_p_wald+12, annotation_id, ha='center', va='center', fontsize=12, color='black')
        # add the connection line
        plt.plot([chr_start_pos[annotation_chr] + annotation_pos, chr_start_pos[annotation_chr] + annotation_pos], [annotation_p_wald+12, annotation_p_wald], color='black', linewidth=0.5, alpha=0.5)
plt.savefig('manhattan_gwas.png', dpi=300, bbox_inches='tight')
