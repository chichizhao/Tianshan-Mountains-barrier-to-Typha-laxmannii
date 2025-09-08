#!/bin/bash
# function: plot the chromosome 5 in typha_lax
# here we plot the pi value, Tajima's D value, Fst and Rasid value in one plot

import matplotlib.pyplot as plt
import pandas as pd

# the chr5 length: "chr5": 17861205,

# read the Pi in the north and south

# north pi
north_pi = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/N_pi.windowed.pi", sep="\t", header=None)

# print the head 
# print(north_pi.head())
"""
       0          1        2           3            4
0  CHROM  BIN_START  BIN_END  N_VARIANTS           PI
1   chr1          1    50000         192  0.000826651
2   chr1      50001   100000         543    0.0027542
3   chr1     100001   150000         415   0.00159258
"""

# we should plot the chr5 only
length = 17861205
north_pi_chr5 = north_pi[north_pi[0] == "chr5"]

south_pi = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/S_pi.windowed.pi", sep="\t", header=None)
# south pi
south_pi_chr5 = south_pi[south_pi[0] == "chr5"]

# read the Tajima's D value in the north and south

# north Tajima's D
north_tajima = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/N_tajima.Tajima.D", sep="\t", header=None)
# print the head
# print(north_tajima.head())
"""
       0          1       2          3
0  CHROM  BIN_START  N_SNPS    TajimaD
1   chr1          0     167   0.412421
2   chr1      50000     438    1.20234
3   chr1     100000     351  0.0222918
4   chr1     150000     345   0.651549
"""

# the Tajima's D value is in the second column
north_tajima_chr5 = north_tajima[north_tajima[0] == "chr5"]

# south Tajima's D
south_tajima = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/S_tajima.Tajima.D", sep="\t", header=None)
south_tajima_chr5 = south_tajima[south_tajima[0] == "chr5"]

# read the Fst value in the north and south
# Fst
Fst = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/fst_50kb.windowed.weir.fst", sep="\t", header=None)

# print the head
# print(Fst.head())

"""
       0          1        2           3             4          5
0  CHROM  BIN_START  BIN_END  N_VARIANTS  WEIGHTED_FST   MEAN_FST
1   chr1      20001    30000          24     0.0163011  0.0135761
2   chr1      30001    40000          83     0.0403032  0.0300312
3   chr1      40001    50000          87     0.0319616   0.020519
4   chr1      50001    60000          61     0.0461447  0.0340955
"""

# the Fst value is in the second column
Fst_chr5 = Fst[Fst[0] == "chr5"]


# get the Rasid value for the north and south
# north Rasid
north_rasid = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/RAiSD_Report.N_clean_modified", sep="\t", header=0)
# give the head of the first column

# print the head
#print(north_rasid.head())

# get the Rasid value for chr5
north_rasid_chr5 = north_rasid[north_rasid["chr"] == 5]

#print(north_rasid_chr5.head())

# south Rasid
south_rasid = pd.read_csv("/home/chichi/data2/chichi/typha_lax/rasid/RAiSD_Report.S_clean_modified", sep="\t", header=0)
# get the head of the first column

# print the head
# print(south_rasid.head())

# get the Rasid value for chr5
south_rasid_chr5 = south_rasid[south_rasid["chr"] == 5]

# first we need to plot the pi value, Tajima's D value, Fst and Rasid value in one plot

fig, ax1 = plt.subplots(figsize=(10, 10))
# set the xlimits
ax1.set_xlim(-1500000, 18500000)

# mark the core sweep region 
box = plt.Rectangle((6871824-100000, -0.01), 2469404+200000, 0.99, facecolor='white', edgecolor='blue', alpha=0.5, linewidth=2, label='Core Sweep Region')
ax1.add_patch(box)
box = plt.Rectangle((6871824-100000, -0.01), 2469404+200000, 0.99, facecolor='blue', edgecolor='blue', alpha=0.1, linewidth=2, label='Core Sweep Region')
ax1.add_patch(box)

# mark the extended sweep region
box = plt.Rectangle((5650000-100000, -0.02), 9750000-5650000+200000, 1.02, facecolor='yellow', edgecolor='red', alpha=0.15, linewidth=2, label='Extended Sweep Region')
ax1.add_patch(box)
box = plt.Rectangle((5650000-100000, -0.02), 9750000-5650000+200000, 1.02, facecolor='white', edgecolor='red', alpha=0.5, linewidth=2, label='Extended Sweep Region')
ax1.add_patch(box)
# plot the north pi value
# set the start point
box = plt.Rectangle((5850000-100000, 0.895), 9600000-5850000+200000, 0.09, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box) 
# add the y label
box = plt.Rectangle((-90000,0.9),20000,0.09, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.9), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.9, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.94), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.94, "4" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.98), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.98, "8" , fontsize=12, ha='right', va='center')
ax1.text(-1600000, 0.9, "π x 0.001" , fontsize=12, ha='center', va='center', rotation=90)
box = plt.Rectangle((-1350000,0.81), 10000, 0.18, color='black')
ax1.add_patch(box)


start_point_y1 = 0.9
# plot the north pi value
for start, end, pi in zip(north_pi_chr5[1], north_pi_chr5[2], north_pi_chr5[4]):
    box = plt.Rectangle((float(start), float(start_point_y1)), float(end) - float(start), float(pi)*10, facecolor='#34587F', edgecolor='#34587F', linewidth=0.5)
    ax1.add_patch(box)
north_pi_sweep_region = []
# here, we need get the core sweep region from the pi value
# for start, end, pi in zip(north_pi_chr5[1], north_pi_chr5[2], north_pi_chr5[4]):
#    if float(pi) < 0.002:
#        north_pi_sweep_region.append(float(start))
# print("Pi core sweep region start positions:", north_pi_sweep_region)


# plot the south pi value
box = plt.Rectangle((5850000-100000, 0.795), 9600000-5850000+200000, 0.09, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0.8),20000,0.09, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.8), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.8, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.84), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.84, "4" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.88), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.88, "8" , fontsize=12, ha='right', va='center')
#ax1.text(-1700000, 0.85, "π x0.001" , fontsize=12, ha='center', va='center', rotation=90)


start_point_y2 = 0.8
for start, end, pi in zip(south_pi_chr5[1], south_pi_chr5[2], south_pi_chr5[4]):
    box = plt.Rectangle((float(start), float(start_point_y2)), float(end) - float(start), float(pi)*10, facecolor='#68AC56', edgecolor='#68AC56', linewidth=0.5)
    ax1.add_patch(box)
# plot the Tajima's D value
#pi_sweep_region = []
# set the bundaries for the pi value
#for start, end, pi in zip(north_pi_chr5[1], north_pi_chr5[2], north_pi_chr5[4]):
#    if float(pi) < 0.002:
#        pi_sweep_region.append(float(start))
#print("Pi core sweep region start positions:", pi_sweep_region)

# plot the north Tajima's D value
box = plt.Rectangle((5700000-50000-100000, 0.68), 9550000-5700000-50000+200000, 0.08, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0.66),20000,0.09, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.7), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.7, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.73), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.73, "+ " , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.67), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.67, "- " , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-1350000,0.56), 10000, 0.18, color='black')
ax1.add_patch(box)
ax1.text(-1600000, 0.65, "Tajima D" , fontsize=12, ha='center', va='center', rotation=90)


start_point_y3 = 0.7
for start, tajima in zip(north_tajima_chr5[1], north_tajima_chr5[3]):
    bin= 50000
    box = plt.Rectangle((float(start), float(start_point_y3)), bin, float(tajima)*0.03, facecolor='#34587F', edgecolor='#34587F', linewidth=0.5)
    ax1.add_patch(box)
south_tajima_sweep_region = []
# here, we need get the core sweep region from the Tajima's D value
#k = 0
#for start, tajima in zip(south_tajima_chr5[1], south_tajima_chr5[3]):
#    if float(tajima) > 0:
#        k =1
#    if float(tajima) < 0:
#        if k == 1:
#            k = 0
#        else:
#            south_tajima_sweep_region.append(float(start))
#print("Tajima's D core sweep region start positions:", south_tajima_sweep_region)


# plot the south Tajima's D value
box = plt.Rectangle((5700000-50000-100000, 0.56), 9550000-5700000+200000, 0.1, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0.56),20000,0.09, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.6), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.6, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.63), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.63, "+ " , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.57), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.57, "- " , fontsize=12, ha='right', va='center')
#ax1.text(-1700000, 0.6, "Tajima D" , fontsize=12, ha='center', va='center', rotation=90)


start_point_y4 = 0.6
for start, tajima in zip(south_tajima_chr5[1], south_tajima_chr5[3]):
    bin = 50000
    box = plt.Rectangle((float(start), float(start_point_y4)), bin, float(tajima)*0.03, facecolor='#68AC56', edgecolor='#68AC56', linewidth=0.5)
    ax1.add_patch(box)
tajima_sweep_region = []
# here, we need get the core sweep region from the Tajima's D value
#k = 0
#for start, tajima in zip(north_tajima_chr5[1], north_tajima_chr5[3]):
#    if float(tajima) > 0:
#        k =1
#    if float(tajima) < 0:
#        if k == 1:
#            k = 0
#        else:
#            tajima_sweep_region.append(float(start))
#print("Tajima's D core sweep region start positions:", tajima_sweep_region)

    



# plot the Fst value
box = plt.Rectangle((6050000-100000, 0.39), 9750000-6050000+200000, 0.15, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0.4),20000,0.15, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.4), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.4, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.47), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.47, "0.07" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.54), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.54, "0.14" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-1350000,0.4), 10000, 0.15, color='black')
ax1.add_patch(box)
ax1.text(-1600000, 0.5, "Fst" , fontsize=12, ha='center', va='center', rotation=90)

start_point_y5 = 0.4
for start, end, fst in zip(Fst_chr5[1], Fst_chr5[2], Fst_chr5[5]):
    box = plt.Rectangle((float(start), float(start_point_y5)), float(end) - float(start), float(fst), color='green', alpha=0.5)
    ax1.add_patch(box)
Fst_sweep_region = []
# here, we need get the get the core sweep region from the Fst
#for start, end, fst in zip(Fst_chr5[1], Fst_chr5[2], Fst_chr5[5]):
#    if float(fst) <0.03:
#        Fst_sweep_region.append(float(start))
#print("Fst core sweep region start positions:", Fst_sweep_region)

# plot the north Rasid value
box = plt.Rectangle((6871824-100000, 0.19), 2469404+200000, 0.18, facecolor='white', edgecolor='red', alpha=0.5, linewidth=1, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0.2),20000,0.19, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0.2), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.2, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.29), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.29, "30" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.38), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.38, "60" , fontsize=12, ha='right', va='center')
#ax1.text(-1300000, 0.3, "-log(u)" , fontsize=12, ha='center', va='center', rotation=90)

start_point_y6 = 0.2
x=[]
y=[]
x_significant = []
y_significant = []
for pos, rasid in zip(north_rasid_chr5["pos"], north_rasid_chr5["RAiSD"]):
    x.append(float(pos))
    y.append(float(rasid)*0.003+start_point_y6)
    if float(rasid) > 5.635:
        x_significant.append(float(pos))
        y_significant.append(float(rasid)*0.003+start_point_y6)

ax1.scatter(x, y, facecolor='#34587F', edgecolor='#34587F', label='North RAiSD', alpha=0.5, s=1)
ax1.scatter(x_significant, y_significant, facecolor='red', edgecolor='red', label='North RAiSD Significant', alpha=0.5, s=1)
# 99.9% snp cutoff line 5.635
# 99.5% snp cutoff line 2.345
# 99% snp cutoff line 1.611
# 95% snp cutoff line 0.5776
# 99.99% snp cutoff line 18.15
north_core_sweep_start = []

# we need get the core sweep region, which is the u value over 5.635 (99.9% snp cutoff line)
#for pos, rasid in zip(north_rasid_chr5["pos"], north_rasid_chr5["RAiSD"]):
#    if float(rasid) > 5.635:
#        north_core_sweep_start.append(float(pos))
#print("North core sweep region start positions:", north_core_sweep_start)
# the start region is 6871824
# the end region is 9341228



# plot the south Rasid value
# add the core sweep region mark
box = plt.Rectangle((6871984-100000, -0.01), 2460754+200000, 0.19, facecolor='white', edgecolor='red', alpha=0.5, label='Core Sweep Region')
ax1.add_patch(box)
# add the y label
box = plt.Rectangle((-90000,0),20000,0.19, color='black')
ax1.add_patch(box)
box = plt.Rectangle((-250000,0), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0, "0" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.09), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.09, "60" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-250000,0.18), 160000, 0.001, color='black')
ax1.add_patch(box)
ax1.text(-300000, 0.18, "120" , fontsize=12, ha='right', va='center')
box = plt.Rectangle((-1350000,0.02), 10000, 0.36, color='black')
ax1.add_patch(box)
ax1.text(-1600000, 0.2, "RAiSD   -log(u)" , fontsize=12, ha='center', va='center', rotation=90)


start_point_y7 = 0
x_south = []
y_south = []
x_south_significant = []
y_south_significant = []
for pos, rasid in zip(south_rasid_chr5["pos"], south_rasid_chr5["RAiSD"]):
    x_south.append(float(pos))
    y_south.append(float(rasid)*0.0015+start_point_y7)
    if float(rasid) > 11.32:
        x_south_significant.append(float(pos))
        y_south_significant.append(float(rasid)*0.0015+start_point_y7)

ax1.scatter(x_south, y_south, facecolor='#68AC56', edgecolor='#68AC56', label='South RAiSD', alpha=0.5, s=1)
ax1.scatter(x_south_significant, y_south_significant, facecolor='red', edgecolor='red', label='South RAiSD Significant', alpha=0.5, s=1)
# set the ylimits

# here, we need get the core sweep regin of the south
core_sweep_start = []
# 99.9% snp cutoff line 11.32
# 99.5% snp cutoff line 5.605
# 99% snp cutoff line 4.333
# 95% snp cutoff line 2.32
# 99.99% snp cutoff line 37.11
# here we need record the core sweep region, whch is the u value overe  11.32 (99.9% snp cutoff line)
#for pos, rasid in zip(south_rasid_chr5["pos"], south_rasid_chr5["RAiSD"]):
#    if float(rasid) > 11.32:
#        core_sweep_start.append(float(pos))
#print("Core sweep region start positions:", core_sweep_start)
# the start region is 6871984
# the end region is 9332738


ax1.set_ylim(-0.05, 1)

# off the axis
ax1.axis('off')
# save the north pi value in the figure
# add the x label
box = plt.Rectangle((0, -0.025), 18000000, 0.001, color='black')
ax1.add_patch(box)

for i in range(0, 18500000, 2000000):
    box = plt.Rectangle((i, -0.032), 10000, 0.007, color='black')
    ax1.add_patch(box)
    ax1.text(i, -0.06, str(i/1000000), fontsize=12, ha='center')
# add the chr5 label
ax1.text(9000000, -0.08, "Chromosome 5 (mb)", fontsize=12, ha='center', va='top')
# here we print the  sweep region 
print("pi north core sweep region start positions:", 5850000, 9600000)
print("pi south core sweep region start positions:", 5850000, 9600000)
print("Tajima's D north core sweep region start positions:", 5650000, 9750000)
print("Tajima's D south core sweep region start positions:", 5650000, 9750000)
print("Fst core sweep region start positions:", 6050000, 9750000)
print("North Rasid core sweep region start positions:", 6871824, 9341228)
print("South Rasid core sweep region start positions:", 6871984, 9332738)

# the core sweep region is 6,871,984 – 9,332,738
# and the extended sweep region is 	5,650,000 – 9,750,000
plt.savefig("/home/chichi/data2/chichi/typha_lax/rasid/chr5.png", dpi=600, bbox_inches='tight')

# the sweep region
# the core sweep region is 6,871,984 – 9,332,738

# we first read the postive selectio signal on chr5
# gff file
gff="/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.gene.gff"
# annotation file
anno="/home/chichi/data2/chichi/typha_lax/ref/Typha_laxmannii.all_annotation.xls"

# read the gff file to get the gene on chr5
import pandas as pd
gff_file = pd.read_csv(gff, sep="\t", header=None, comment="#")
# filter the gff file to get the gene on chr5
chr5_gff = gff_file[gff_file[0] == "chr5"]
# filter the gene type
chr5_gene = chr5_gff[chr5_gff[2] == "mRNA"]
# get the gene id
chr5_gene_ids = chr5_gene[8].str.split("ID=").str[1].str.split(";").str[0]
# get the gene start and end position
chr5_gene_starts = chr5_gene[3]
chr5_gene_ends = chr5_gene[4]
# create a dataframe
chr5_gene_df = pd.DataFrame({"gene_id": chr5_gene_ids, "start": chr5_gene_starts, "end": chr5_gene_ends})
# filter the gene in the sweep region
sweep_start = 6871984
sweep_end = 9332738
sweep_genes = chr5_gene_df[(chr5_gene_df["start"] >= sweep_start) & (chr5_gene_df["end"] <= sweep_end)]

print("Genes in the sweep region on chr5:")
print(sweep_genes)

# read the annotation file to get the gene function
anno_file = pd.read_csv(anno, sep="\t", header=0)
# filter the annotation file to get the genes in the sweep region
# write the annotation into a file
with open("chr5_sweep_genes_annotation.txt", "w") as f:
    header = anno_file.columns.tolist()
    f.write("chr"\t +"start" + "\t" + "end" + "\t" + "\t".join(header) + "\n")
    for gene_id in sweep_genes["gene_id"]:
        chr = "chr5"
        start = sweep_genes[sweep_genes["gene_id"] == gene_id]["start"].values[0]
        end = sweep_genes[sweep_genes["gene_id"] == gene_id]["end"].values[0]
        gene_anno = anno_file[anno_file["SeqID"] == gene_id]
        if not gene_anno.empty:
            gene_anno.to_csv(f, sep="\t", header=False, index=False)
        else:
            f.write(f"No annotation found for gene {gene_id}\n")
