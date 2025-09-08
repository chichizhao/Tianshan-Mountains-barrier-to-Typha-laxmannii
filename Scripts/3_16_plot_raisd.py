
#!/bin/py
# -*- coding: utf-8 -*-
# function: plot the RAISD results with manhattan plot
# usage: python plot_raisd.py -i raisd_result -o output_prefix
# author: chichi

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


#chr_len ={'Chr01': 24370473,'Chr02': 18857030,'Chr03': 18714887,'Chr04': 18467694,'Chr05': 17512124,'Chr06': 14553481,'Chr07': 14503471,'Chr08': 13751786,'Chr09': 11865621,'Chr10': 11839088,'Chr11': 11614203,'Chr12': 10714389,'Chr13': 10242974,'Chr14': 9831573,'Chr15': 9461071,'Contig1': 452046,'Contig2': 203240}
#chr_start_pos= {'Chr01': 0,'Chr02': 24370473,'Chr03': 43227503,'Chr04': 61942390,'Chr05': 80410084,'Chr06': 97922208,'Chr07': 112775689,'Chr08': 127279160,'Chr09': 141030946,'Chr10': 152896567,'Chr11': 164735655,'Chr12': 176349858,'Chr13': 187064247,'Chr14': 197307221,'Chr15': 207138794,'Contig1': 216599865,'Contig2': 261846911}
chr_len ={"chr1": 24094684, "chr2": 19094038, "chr3": 18910458, "chr4": 18888309, "chr5": 17861205,
          "chr6": 15065052, "chr7": 14536590, "chr8": 13397664, "chr9": 13040483,
          "chr10": 12597730, "chr11": 11722036, "chr12": 9996585, "chr13": 9906942,
          "chr14": 9656787, "chr15": 9218693}
chr_start_pos = {"chr1": 0, "chr2": 24094684, "chr3": 43088722, "chr4": 61999180, "chr5": 80887589,
                "chr6": 98748794, "chr7": 113013846, "chr8": 127550436, "chr9": 141948100,
                "chr10": 154988583, "chr11": 167586313, "chr12": 179308349, "chr13": 189304934,
                "chr14": 199211876, "chr15": 208868563}

#def plot_raisd(input_file, output_prefix):  
data ={}
with open("/home/chichi/data2/chichi/typha_lax/rasid/RAiSD_Report.N_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())

# read the gene annotation file for the genome
file = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"
df2 = pd.read_csv(file, sep="\t")




# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 5.635
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 2.345
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 1.611
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 0.5776
print(snp_values[int(len(snp_values)*0.9999)])
# 99.99% snp cutoff line 18.15
print(snp_values[int(len(snp_values)*0.9995)])




# plot the manhattan plot for var of each chr
fig = plt.figure(figsize=(8,2))
ax = fig.add_subplot(111)
ax.set_xlim(-3500000, 208868563+9218693+4000000)
ax.set_ylim(-500, 500)
box = plt.Rectangle((0, 490), 208868563+9218693, 10, fc='#34587F', ec='#34587F', lw=1)
ax.add_patch(box)
num = 0
ax.axis('off')

# here we need mark the significant snp with their name 
#for chr in data:
gff_file = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.gene.gff.mrna"
gff_data = pd.read_csv(gff_file, sep="\t", header=None)
print(gff_data.head())

# read the gene annotation file for the genome
annotation = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"
annotation_data = pd.read_csv(annotation, sep='\t',index_col= False)
annotation_dict = {}
for i in range(len(annotation_data)):
    if str(annotation_data.loc[i, 'Pfam']).startswith('PF'):
        annotation_dict[annotation_data.loc[i, 'SeqID']] = annotation_data.loc[i, 'Pfam'].split('(')[0]
# print(annotation_dict.keys())

with open("S_sig_snp.txt", "w") as f1:
    positive_select_pos = []
    x = []
    y = []
    x_sig = []
    y_sig = []
    x_sig_marker = []
    y_sig_marker = []
    sig_gene = []
    for chr in data:

        gff_data_chr = gff_data[gff_data[0] == chr]
        print(gff_data_chr.head())
        for i in range(len(data[chr])):
            x.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
            if float(data[chr][i][6]) > 10:
                x_sig.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                positive_select_pos.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                y_sig.append(float(data[chr][i][6]))
                # to determine the significant snp in which gene or not
                if float(data[chr][i][6]) > 10:
                    for index, row in gff_data_chr.iterrows():
                        if int(row[3]) <= int(data[chr][i][1]) <= int(row[4]) or int(row[3]) >= int(data[chr][i][1]) + 1 >= int(row[4]):
                            #print(row[3], row[4], data[chr][i][1])
                            x_sig_marker.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                            y_sig_marker.append(float(data[chr][i][6]))
                            id = row[8].split(";Parent=")[0].replace("ID=", "")
                            id = str(id)
                            print(id)
                            print(data[chr][i][1], id)
                            #if id in annotation_dict.keys():
                            sig_gene.append([int(data[chr][i][1])+int(chr_start_pos[chr]), float(data[chr][i][6]), id])
                            #    sig_gene.append(annotation_dict[id])
                            #print(sig_gene)
                            break
            y.append(float(data[chr][i][6]))
        #print(chr, positive_select_pos)
        f1.write(chr+"\n")
        for i in positive_select_pos:
            f1.write(str(i)+",")
        f1.write("\n")

        num += 1# add the chromosome length scale
        #if chr == "Chr02":
        #    break
    # add the significant line 10.04
    # give the y value multiply 4 
    y = [i*4 for i in y]
    y_sig = [i*4 for i in y_sig]
    y_sig_marker = [i*4 for i in y_sig_marker]
    ax.scatter(x, y, s=1, label=chr, color = "grey")
    ax.scatter(x_sig, y_sig, s=1, color = "blue")
    ax.scatter(x_sig_marker, y_sig_marker, s=3, color = "red")        
    #ax.axhline(y=10.04, color='b', linestyle='--', linewidth=0.5)
    #ax.axhline(y=30.32, color='r', linestyle='--', linewidth=0.5)
# print the number of significant snp
# print(len(positive_select_pos))


print(sig_gene)
# add the gene annotation
# write the sig_gene to the file
# unique the gene within the same gene
if len(sig_gene) != 0:
    unique_sig_gene = [sig_gene[0]]
    for i in range(1, len(sig_gene)):
        print(sig_gene[i])
        if sig_gene[i][2] != sig_gene[i-1][2]:
            unique_sig_gene.append(sig_gene[i])
            #ax.text(sig_gene[i][0], sig_gene[i][1]+5, sig_gene[i][2], fontsize=5)
        else:
            continue

    print("unique_sig_gene", unique_sig_gene)


    gene_name_pos = unique_sig_gene
    print(gene_name_pos)
    print(gene_name_pos[0])
    print(gene_name_pos[0][0])
    gene_count = 0
    for new_gene in gene_name_pos:
        #new_num = new_num + 1
        # replace the Tlax in the new_gene[0] with the blank space
        new_gene[2] = new_gene[2].replace("Tlax", "")
        gene_count += 1 
        if gene_count == 2:
            arrow = plt.Arrow(new_gene[0], (new_gene[1]+ 5)*4, 0, 20*4, width=1, color='black',alpha=0.3, lw=0.3)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0] , (new_gene[1]+30)*4, new_gene[2], fontsize=7, va='center', ha='center')


        elif gene_count == 3:
            arrow = plt.Arrow(new_gene[0], (new_gene[1]+5)*4, 0, 40*4, width=7, color='black',alpha=0.3, lw=0.3)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0]+ 2000000, (new_gene[1]+50)*4, new_gene[2], fontsize=7, va='center', ha='center')

        elif gene_count == 4:
            arrow = plt.Arrow(new_gene[0], (new_gene[1]+40)*4, 0, 140*4, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0] -4000000, (new_gene[1]+200)*4, new_gene[2], fontsize=7, va='center', ha='left')
        else:
            arrow = plt.Arrow(new_gene[0], (new_gene[1]+5)*4, 0, 60*4, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], (new_gene[1]+ 70)*4, new_gene[2], fontsize=7, va='center', ha='center')
    # add the mark of the chromosome

df2_columns = df2.columns.tolist()
with open("N_sig_gene.txt", "w") as f2:
    f2.write("GeneID\tchr\tpos\tμ\t"+"\t".join(df2_columns[1:])+"\n")
    for gene in sig_gene:
        gene_id = gene[2]
        chr = gene_id.split("G")[0].replace("Chr", "chr").replace("Tlax", "")
        pos = gene[0] - chr_start_pos[chr]
        mu = gene[1]
        if gene_id in df2['SeqID'].values:
            f2.write(gene_id + "\t" + chr + "\t" + str(pos) + "\t" + str(mu) +"\t"+ "\t".join([str(df2.loc[df2['SeqID'] == gene_id, col].values[0]) for col in df2_columns[1:]]) + "\n")
        else:
            f2.write(gene_id + "\t" + chr + "\t" + str(pos) + "\t" + str(mu) + "\t" + "No function annotated" + "\n")
    
        

for i in range(1, 16):
    if i < 10:
        chr = 'chr' + str(i)
    else:
        chr = 'chr' + str(i)
    #ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -20, str(i), ha='center', va='center', fontsize=8, color='black')
    # add the chromosome region line 
    #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
    if i % 2 == 1:
        # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
        box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 870, fc='grey', ec='none', lw=0.5, alpha=0.2)
        ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
#ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
#ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')
box = plt.Rectangle((-3500000, 0), 1, 500, fc='black', ec='black', lw=1)
ax.add_patch(box)
#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)
ax.axhline(y=10*4, color='blue', linestyle='--', alpha=0.5, linewidth=0.4)
# add the y axis label
#ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
ax.text(1800000, 10*4, "10", ha='right', va='bottom', fontsize=7, color='blue')
#ax.text(-12500000, 60*4, "-log10(µ)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
#ax.text(-1000000, -20, "Chr", ha='center', va='center', fontsize=7, color='black')
ax.text(-3700000, 0, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
ax.text(-3700000, 50*4, "50", ha='right', va='center', fontsize=7)
#ax.text(-200000, 60, "60", ha='right', va='center', fontsize=10)
ax.text(-3700000, 100*4, "100", ha='right', va='center', fontsize=7)
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)

box = plt.Rectangle((208868563+9218693, 0), 4000000, 500, fc='#34587F', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(208868563+9218693 + 2000000, 60*4, "North", ha='center', va='center', fontsize=7, color='black', rotation= 90)
ax.axhline(y=0, color='black', linewidth=0.6)

data ={}
with open("/home/chichi/data2/chichi/typha_lax/rasid/RAiSD_Report.S_clean") as f:
    for line in f:
        line = line.strip()

        if line.startswith("// ch"):
            chr = line.split( )[1]
            print(chr)
            data[line.split( )[1]] = []
            print(chr)
        elif line.startswith("// Contig"):
            break
        else:
            data[chr].append(line.split( ))
print(data.keys())

# here need the 99% snp cutoff line
snp_values = []
for chr in data:
    for i in range(len(data[chr])):
        snp_values.append(float(data[chr][i][6]))
snp_values.sort()
print(snp_values[int(len(snp_values)*0.999)])
# 99.9% snp cutoff line 11.32
print(snp_values[int(len(snp_values)*0.995)])
# 99.5% snp cutoff line 5.605
print(snp_values[int(len(snp_values)*0.99)])
# 99% snp cutoff line 4.333
print(snp_values[int(len(snp_values)*0.95)])
# 95% snp cutoff line 2.32
print(snp_values[int(len(snp_values)*0.9999)])
# 99.99% snp cutoff line 37.11
print(snp_values[int(len(snp_values)*0.9995)])





#ax.set_xlim(-3500000, 208868563+9218693+4000000)
#ax.set_ylim(-25, 500)
box = plt.Rectangle((0, -500), 208868563+9218693, 10, fc='#68AC56', ec='#68AC56', lw=1)
ax.add_patch(box)
num = 0
ax.axis('off')

# here we need mark the significant snp with their name 
#for chr in data:
gff_file = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.gene.gff.mrna"
gff_data = pd.read_csv(gff_file, sep="\t", header=None)
print(gff_data.head())

# read the gene annotation file for the genome
annotation = "/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"
annotation_data = pd.read_csv(annotation, sep='\t',index_col= False)
annotation_dict = {}
for i in range(len(annotation_data)):
    if str(annotation_data.loc[i, 'Pfam']).startswith('PF'):
        annotation_dict[annotation_data.loc[i, 'SeqID']] = annotation_data.loc[i, 'Pfam'].split('(')[0]
# print(annotation_dict.keys())

with open("S_sig_snp.txt", "w") as f1:
    positive_select_pos = []
    x = []
    y = []
    x_sig = []
    y_sig = []
    x_sig_marker = []
    y_sig_marker = []
    sig_gene = []
    for chr in data:

        gff_data_chr = gff_data[gff_data[0] == chr]
        print(gff_data_chr.head())
        for i in range(len(data[chr])):
            x.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
            if float(data[chr][i][6]) > 30:
                x_sig.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                positive_select_pos.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                y_sig.append(float(data[chr][i][6]))
                # to determine the significant snp in which gene or not
                if float(data[chr][i][6]) > 30:
                    for index, row in gff_data_chr.iterrows():
                        if int(row[3]) <= int(data[chr][i][1]) <= int(row[4]) or int(row[3]) >= int(data[chr][i][1]) + 1 >= int(row[4]):
                            #print(row[3], row[4], data[chr][i][1])
                            x_sig_marker.append(int(data[chr][i][1])+int(chr_start_pos[chr]))
                            y_sig_marker.append(float(data[chr][i][6]))
                            id = row[8].split(";Parent=")[0].replace("ID=", "")
                            id = str(id)
                            print(data[chr][i][1], id)
                            #if id in annotation_dict.keys():
                            sig_gene.append([int(data[chr][i][1])+int(chr_start_pos[chr]), float(data[chr][i][6]), id])
                            #    sig_gene.append(annotation_dict[id])
                            #print(sig_gene)
                            break
            y.append(float(data[chr][i][6]))
        #print(chr, positive_select_pos)
        f1.write(chr+"\n")
        for i in positive_select_pos:
            f1.write(str(i)+",")
        f1.write("\n")

        num += 1# add the chromosome length scale
        #if chr == "Chr02":
        #    break
    # add the significant line 10.04
    # give the y value the negative value
    y = [-i for i in y]
    y_sig = [-i for i in y_sig]
    y_sig_marker = [-i for i in y_sig_marker]

    ax.scatter(x, y, s=1, label=chr, color = "grey")
    ax.scatter(x_sig, y_sig, s=1, color = "blue")
    ax.scatter(x_sig_marker, y_sig_marker, s=3, color = "red")        
    #ax.axhline(y=10.04, color='b', linestyle='--', linewidth=0.5)
    #ax.axhline(y=30.32, color='r', linestyle='--', linewidth=0.5)
# print the number of significant snp
# print(len(positive_select_pos))


print(sig_gene)
# add the gene annotation
# write the sig_gene to the file 
# with open("S_sig_gene.txt", "w") as f2:
#     for gene in sig_gene:
#f2.write("\t".join([str(i) for i in gene]) + "\n")
# unique the gene within the same gene
if len(sig_gene) != 0:
    unique_sig_gene = [sig_gene[0]]
    for i in range(1, len(sig_gene)):
        print(sig_gene[i])
        if sig_gene[i][2] != sig_gene[i-1][2]:
            unique_sig_gene.append(sig_gene[i])
            #ax.text(sig_gene[i][0], sig_gene[i][1]+5, sig_gene[i][2], fontsize=5)
        else:
            continue

    print("unique_sig_gene", unique_sig_gene)


    gene_name_pos = unique_sig_gene
    print(gene_name_pos)
    print(gene_name_pos[0])
    print(gene_name_pos[0][0])
    gene_count = 0
    for new_gene in gene_name_pos:
        #new_num = new_num + 1
        # replace the Tlax in the new_gene[0] with the blank space
        new_gene[2] = new_gene[2].replace("Tlax", "")
        gene_count += 1 
        if gene_count == 2 or gene_count == 5  :
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -100, width=1, color='black',alpha=0.3, lw=0.3)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+260), new_gene[2], fontsize=7, va='center', ha='center')


        elif gene_count == 3:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -240, width=7, color='black',alpha=0.3, lw=0.3)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+350), new_gene[2], fontsize=7, va='center', ha='center')

        elif gene_count == 4:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -140, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+200), new_gene[2], fontsize=7, va='center', ha='center')

        elif gene_count == 7:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -140, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+200), new_gene[2], fontsize=7, va='center', ha='center')
        elif gene_count == 8:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -140, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+230), new_gene[2], fontsize=7, va='center', ha='center')

        elif gene_count == 9:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -100, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0]-4000000, -(new_gene[1]+150), new_gene[2], fontsize=7, va='center', ha='left')
        elif gene_count == 6:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -140, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+100), new_gene[2], fontsize=7, va='center', ha='right')
                               
        else:
            arrow = plt.Arrow(new_gene[0], -(new_gene[1]+40), 0, -140, width=1, color='black',alpha=0.3, lw=0.5)
            ax.add_patch(arrow)
            print(new_gene)               
            ax.text(new_gene[0], -(new_gene[1]+200), new_gene[2], fontsize=7, va='center', ha='center')
    # add the mark of the chromosome
with open("S_sig_gene.txt", "w") as f2:
    f2.write("GeneID\tchr\tpos\tμ\t"+"\t".join(df2_columns[1:])+"\n")
    for gene in sig_gene:
        gene_id = gene[2]
        chr = gene_id.split("G")[0].replace("Chr", "chr").replace("Tlax", "")
        pos = gene[0] - chr_start_pos[chr]
        mu = gene[1]
        if gene_id in df2['SeqID'].values:
            f2.write(gene_id + "\t" + chr + "\t" + str(pos) + "\t" + str(mu) +"\t" + "\t".join([str(df2.loc[df2['SeqID'] == gene_id, col].values[0]) for col in df2_columns[1:]]) + "\n")
        else:
            f2.write(gene_id + "\t" + chr + "\t" + str(pos) + "\t" + str(mu) + "\t" + "No function annotated" + "\n")


for i in range(1, 16):
    if i < 10:
        chr = 'chr' + str(i)
    else:
        chr = 'chr' + str(i)
    ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -550, str(i), ha='center', va='center', fontsize=8, color='black')
    # add the chromosome region line 
    #ax.axvline(x=chr_start_pos[chr], color='black', linestyle='-', alpha=0.2)
    if i % 2 == 1:
        # ax.text(chr_start_pos[chr] + chr_len[chr] / 2, -10, str(i), ha='center', va='center', fontsize=7, color='black')
        box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], -870, fc='grey', ec='none', lw=0.5, alpha=0.2)
        ax.add_patch(box)
    #box = plt.Rectangle((chr_start_pos[chr], 0), chr_len[chr], 100, fc='grey', ec='none', lw=0.5, alpha=0.2)
#ax.axvline(x=chr_start_pos['Contig1'], color='black', linestyle='-')
#ax.axvline(x=chr_start_pos["Chr01"], color='black', linestyle='-')
box = plt.Rectangle((-3500000, 0), 1, -500, fc='black', ec='black', lw=1)
ax.add_patch(box)
#ax.axhline(y=14, color='black', linestyle='--', alpha=0.5)
ax.axhline(y=-30, color='blue', linestyle='--', alpha=0.5, linewidth=0.4)
# add the y axis label
#ax.text(-1800000, 10.5, "99%", ha='right', va='center', fontsize=10, color='blue')
ax.text(1800000, -30, "30", ha='right', va='top', fontsize=7, color='blue')
ax.text(-15500000, 0, "-log10(µ)", ha='center', va='center', fontsize=7, color='black', rotation= 90)
ax.text(-1000000, -550, "Chr", ha='center', va='center', fontsize=7, color='black')
ax.text(-3700000, 0, "0", ha='right', va='center', fontsize=7, color='black')
#ax.text(-200000, 20, "20", ha='right', va='center', fontsize=7)
ax.text(-3700000, -200, "200", ha='right', va='center', fontsize=7)
#ax.text(-200000, 60, "60", ha='right', va='center', fontsize=10)
ax.text(-3700000, -400, "400", ha='right', va='center', fontsize=7)
#ax.text(-200000, 100, "100", ha='right', va='center', fontsize=10)

box = plt.Rectangle((208868563+9218693, 0), 4000000, -500, fc='#68AC56', ec='none', lw=0.5, alpha=0.7)
ax.add_patch(box)
ax.text(208868563+9218693 + 2000000, -200, "South", ha='center', va='center', fontsize=7, color='black', rotation= 90)





plt.savefig("N_S.png", dpi=1200, bbox_inches='tight')
