#!/bin/python
# -*- coding: utf-8 -*-
# function: plot the population bio data from world climate data

import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt


# DATA
# data set 
# N populations
Xu0856 = [1.9875, 13.1250, 24.3507, 1518.2690, 26.9000, -27.0000, 53.9000, 19.1000, -14.6500, 19.1000, -18.2667, 160.0000, 27.0000, 6.0000, 41.2571, 61.0000, 19.0000, 61.0000, 22.0000]
Xu2534 = [3.2417, 12.0167, 26.1801, 1273.8130, 24.4000, -21.5000, 45.9000, 17.7833, -13.6000, 17.7833, -13.6000, 170.0000, 34.0000, 2.0000, 73.7344, 91.0000, 8.0000, 91.0000, 8.0000]
Xu2542 = [5.5750, 11.0667, 23.0076, 1398.8610, 29.6000, -18.5000, 48.1000, 21.7500, -10.1167, 21.7500, -12.4167, 192.0000, 31.0000, 7.0000, 38.6953, 72.0000, 24.0000, 72.0000, 31.0000]
Xu2565 = [11.1917, 15.5000, 32.3591, 1199.5220, 33.0000, -14.9000, 47.9000, 24.2167, -3.7167, 24.2167, -4.8167, 56.0000, 14.0000, 1.0000, 77.9576, 33.0000, 3.0000, 33.0000, 3.0000]
Xu4091 = [10.7750, 15.4333, 30.1432, 1328.2600, 35.2000, -16.0000, 51.2000, 26.0000, -2.4500, 26.0000, -6.6667, 38.0000, 9.0000, 1.0000, 62.1084, 21.0000, 3.0000, 21.0000, 3.0000]
Xu4097 = [8.5333, 11.2667, 21.7503, 1513.2510, 33.0000, -18.8000, 51.8000, 24.0500, -7.9000, 25.7167, -11.4667, 151.0000, 24.0000, 5.0000, 42.5188, 61.0000, 18.0000, 58.0000, 20.0000]
Xu4124 = [5.6458, 12.2917, 23.7291, 1471.8590, 30.3000, -21.5000, 51.8000, 20.6333, -11.0000, 22.3500, -13.6833, 169.0000, 23.0000, 5.0000, 37.8707, 57.0000, 19.0000, 56.0000, 21.0000]
Xu4132 = [4.8667, 13.2000, 24.7191, 1499.6580, 30.3000, -23.1000, 53.4000, 20.1167, -12.0000, 21.9167, -14.8167, 162.0000, 22.0000, 6.0000, 31.7405, 53.0000, 21.0000, 52.0000, 25.0000]
Xu4141 = [1.3292, 12.9250, 24.1589, 1519.0910, 25.5000, -28.0000, 53.5000, 18.1667, -15.3000, 18.1667, -19.1667, 154.0000, 27.0000, 6.0000, 44.6508, 61.0000, 18.0000, 61.0000, 21.0000]
Xu4158 = [3.2417, 13.9833, 25.7520, 1518.6710, 28.8000, -25.5000, 54.3000, 18.5500, -13.7667, 20.4500, -16.7500, 149.0000, 22.0000, 5.0000, 36.3003, 50.0000, 18.0000, 50.0000, 24.0000]
Xu4166 = [5.3333, 12.6167, 24.3096, 1465.3150, 30.1000, -21.8000, 51.9000, 20.3167, -11.2833, 21.9833, -13.8167, 154.0000, 22.0000, 6.0000, 32.3042, 51.0000, 21.0000, 51.0000, 28.0000]
Xu4216 = [4.9083, 10.9500, 23.7527, 1354.4270, 28.3000, -17.8000, 46.1000, 18.9000, -10.4333, 20.5833, -12.3833, 250.0000, 39.0000, 9.0000, 38.3120, 92.0000, 31.0000, 91.0000, 39.0000]
Xu4224 = [3.5708, 10.8417, 24.5843, 1276.4480, 26.2000, -17.9000, 44.1000, 16.8333, -10.6167, 18.7667, -12.2833, 280.0000, 44.0000, 9.0000, 42.5761, 109.0000, 32.0000, 103.0000, 35.0000]
Xu4232 = [3.3625, 11.4083, 25.4650, 1269.2970, 25.8000, -19.0000, 44.8000, 16.5333, -10.7333, 18.3333, -12.6500, 253.0000, 41.0000, 7.0000, 48.8270, 104.0000, 25.0000, 100.0000, 26.0000]
Xu4245 = [7.8333, 13.1667, 26.9809, 1315.4840, 32.4000, -16.4000, 48.8000, 21.2667, -6.8167, 23.2333, -9.0333, 215.0000, 29.0000, 9.0000, 36.6566, 85.0000, 31.0000, 73.0000, 35.0000]
Xu4248 = [6.5875, 11.8083, 24.0496, 1406.8790, 30.4000, -18.7000, 49.1000, 20.8000, -9.2667, 22.6667, -11.7667, 203.0000, 31.0000, 7.0000, 46.1150, 87.0000, 23.0000, 76.0000, 25.0000]
Xu4260 = [5.4792, 11.0417, 23.6438, 1350.1460, 29.2000, -17.5000, 46.7000, 19.2333, -11.5333, 21.3167, -11.5333, 231.0000, 38.0000, 6.0000, 52.1822, 102.0000, 22.0000, 93.0000, 22.0000]
Xu4272 = [4.0167, 13.1667, 29.0015, 1201.7910, 26.2000, -19.2000, 45.4000, 16.0833, -9.8000, 17.9000, -11.5667, 314.0000, 43.0000, 14.0000, 33.4372, 116.0000, 50.0000, 103.0000, 50.0000]
Xu4290 = [9.0625, 13.1583, 27.6436, 1238.0650, 31.0000, -16.6000, 47.6000, 17.8167, -4.8500, 22.7667, -7.7333, 225.0000, 27.0000, 12.0000, 25.3739, 74.0000, 43.0000, 61.0000, 47.0000]
Xu4305 = [8.6375, 13.3917, 27.5549, 1267.2680, 30.8000, -17.8000, 48.6000, 17.6500, -5.6667, 22.5500, -8.6500, 238.0000, 27.0000, 13.0000, 23.8541, 77.0000, 46.0000, 62.0000, 51.0000]
Xu4310 = [5.7708, 13.1250, 28.0449, 1221.7810, 27.5000, -19.3000, 46.8000, 17.6333, -10.7500, 19.1500, -10.7500, 258.0000, 39.0000, 9.0000, 44.5050, 107.0000, 31.0000, 101.0000, 31.0000]
Xu4322 = [6.5792, 12.7917, 27.9294, 1202.9420, 27.5000, -18.3000, 45.8000, 18.3833, -9.7333, 19.7167, -9.7333, 232.0000, 36.0000, 9.0000, 45.3339, 97.0000, 29.0000, 92.0000, 29.0000]
Xu4328 = [1.2583, 11.3167, 26.2568, 1186.7560, 21.5000, -21.6000, 43.1000, 14.3500, -14.5667, 14.3500, -14.5667, 223.0000, 42.0000, 3.0000, 79.8076, 118.0000, 9.0000, 118.0000, 9.0000]
Xu4336 = [10.0042, 13.1250, 29.1020, 1192.7500, 30.5000, -14.6000, 45.1000, 23.2000, -4.4333, 23.2000, -5.8333, 90.0000, 21.0000, 1.0000, 88.9631, 54.0000, 3.0000, 54.0000, 3.0000]
Xu4343 = [9.7917, 11.2667, 25.0370, 1219.8580, 30.1000, -14.9000, 45.0000, 23.3500, 2.5000, 23.3500, -6.5500, 76.0000, 18.0000, 1.0000, 82.5040, 45.0000, 4.0000, 45.0000, 4.0000]
Xu4397 = [11.3208, 15.5750, 32.3805, 1212.4400, 33.0000, -15.1000, 48.1000, 24.5833, -3.6500, 24.5833, -4.8333, 41.0000, 11.0000, 0.0000, 84.4183, 26.0000, 0.0000, 26.0000, 1.0000]
Xu4420 = [11.4292, 13.5083, 29.7540, 1191.5560, 32.1000, -13.3000, 45.4000, 24.4667, -4.5833, 24.4667, -4.5833, 133.0000, 30.0000, 2.0000, 88.6054, 81.0000, 6.0000, 81.0000, 6.0000]
Xu4425 = [7.6625, 10.9083, 21.2638, 1504.3440, 31.4000, -19.9000, 51.3000, 22.6333, -9.0667, 24.3167, -12.5833, 169.0000, 23.0000, 7.0000, 36.1431, 63.0000, 23.0000, 54.0000, 24.0000]
Xu4438 = [8.0792, 11.2750, 21.8085, 1520.9960, 31.8000, -19.9000, 51.7000, 23.2167, -8.7500, 24.7167, -12.5500, 159.0000, 22.0000, 6.0000, 41.7985, 63.0000, 20.0000, 56.0000, 20.0000]
Xu4440 = [8.3250, 10.8833, 20.9295, 1530.3390, 32.3000, -19.7000, 52.0000, 23.5000, -8.5833, 25.1167, -12.4667, 149.0000, 20.0000, 5.0000, 43.9755, 59.0000, 17.0000, 54.0000, 17.0000]
Xu4453 = [8.4625, 11.1083, 21.7384, 1511.4280, 31.6000, -19.5000, 51.1000, 23.4667, -8.1333, 24.8000, -12.1667, 139.0000, 20.0000, 5.0000, 47.4961, 59.0000, 16.0000, 54.0000, 16.0000]
# List of populations



# Step 2: Combine all the lists into a single data frame
data = np.array([Xu0856, Xu2534, Xu2542, Xu2565, Xu4091, Xu4097, Xu4124, Xu4132, Xu4141, Xu4158, Xu4166, Xu4216, Xu4224, Xu4232, Xu4245, Xu4248, Xu4260, Xu4272, Xu4290, Xu4305, Xu4310, Xu4322, Xu4328, Xu4336, Xu4343, Xu4397, Xu4420, Xu4425, Xu4438, Xu4440, Xu4453])   
df = pd.DataFrame(data)
# add the group name to the data frame
#df['group'] = ['N4', 'S', 'N4', 'S', 'S', 'N3', 'N4', 'N4', 'N4', 'N4', 'N4', 'N4', 'N4', 'N3', 'N3', 'N3', 'N3', 'N2', 'N2', 'N2', 'N2', 'N2', 'N2', 'S', 'S', 'S', 'S', 'N1', 'N1', 'N1', 'N1']
df['group'] = ['N', 'S', 'N', 'S', 'S', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'S', 'S', 'S', 'S', 'N', 'N', 'N', 'N']



#print(df)          
# add the column name to the data frame
df.columns = ['bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19','group']

print(df)

# unique group names
groups = df['group'].unique()
print(groups)
#groups = ['N1', 'N2', 'N3', 'N4', 'S']
#colors = {  'N1':'#2F5597', 'N2':'#7E8BDC','N3':'#00B0F0','N4':'#08C2BE','S':'#68AC56'}
groups = ['N', 'S']
colors = {  'N':'#00B0F0','S':'#68AC56'}
#  we need first to identify the similar bio data 
# here we use the R code to identify the similar bio data
# here is the result
unique_group_bios =[["bio1","bio6","bio9","bio11"],
                    ["bio12","bio13","bio16","bio18"],
                    ["bio8","bio5","bio10"],
                    ["bio14","bio17","bio19"],
                    ["bio4","bio7"],    
                    ["bio2"],
                    ["bio3"],
                    ["bio15"]]

# set the font into arail
plt.rcParams['font.sans-serif'] = ['Arial']

fig, ax = plt.subplots()
fig.set_size_inches(14,8)
ax.set_xlim(-0.5,15.5)
ax.set_ylim(2.2,12.5)
t =0
for i in range(8):
    for t in range(len(unique_group_bios[i])):
        if (i%2 == 0 and t %2 == 0) or (i%2 == 1 and t %2 == 1):
            print(i)
            box = plt.Rectangle((-0.5+i*2,11-t*2.5-1.1),2,2.5,fc='grey',ec='k',lw=0.2,alpha=0.1)
            ax.add_artist(box)
            box = plt.Rectangle((-0.5+i*2,11-t*2.5-1.1),2,2.2,fc='none',ec='k',lw=0.5)
            ax.add_artist(box)
            box = plt.Rectangle((-0.5+i*2,11-t*2.5-1.1),2,2.5,fc='none',ec='k',lw=0.5)
            ax.add_artist(box)
            ax.text(-0.5+i*2+0.96 ,11-t*2.5-1+2.15,str(unique_group_bios[i][t]),fontsize=14,horizontalalignment='center')
        else:
            box = plt.Rectangle((-0.5+i*2,11-t*2.5-1.1),2,2.2,fc='none',ec='k',lw=0.5)
            ax.add_artist(box)
            box = plt.Rectangle((-0.5+i*2,11-t*2.5-1.1),2,2.5,fc='none',ec='k',lw=0.5)
            ax.add_artist(box)
            ax.text(-0.5+i*2+0.96 ,11-t*2.5-1+2.15,str(unique_group_bios[i][t]),fontsize=14,horizontalalignment='center')


for i in range(19):
    print(i)
    bio_name = 'bio'+str(i+1)
    print(bio_name)

    max_value = np.max(df[bio_name])
    min_value = np.min(df[bio_name])
    t = 0
    # match the position of the bio data
    for k1 in range(len(unique_group_bios)):
        for k2 in range(len(unique_group_bios[k1])):
            if bio_name == unique_group_bios[k1][k2]:
                i = k1
                t = k2
                print(i,t)
                for j in range(len(groups)):
                    bio_data = df[df['group']==groups[j]]
                    bio_data = bio_data[bio_name]
                    bio_data = (bio_data-min_value)/(max_value-min_value)*2*0.8-t*2.5+10
                    print(bio_data)
                    print("this is the bio data for "+bio_name+'_'+groups[j])
                    if j == 0:
                        bios_data_1 = bio_data
                    elif j == 1:
                        bios_data_2 = bio_data


                    pos = [2*(i+0.6)-0.6*j-0.6]*len(bio_data)
                    plt.plot(pos,bio_data,'o',color=colors[groups[j]],markersize=3.5,alpha=0.9)

                    
                    
                    plt.boxplot(bio_data, positions=[2*(i+0.6)-0.6*j-0.6], showfliers=False, widths=0.2)
                    """
                    parts1 =plt.violinplot(bio_data, showmedians=True,positions=[2*(i+0.6)-0.6*j-0.6], showextrema=False, widths=0.3)
                    for pc in parts1['bodies']:
                        pc.set_facecolor(colors[groups[j]])
                        pc.set_edgecolor('black')
                        pc.set_alpha(0.5)
                    """
                step = (max_value-min_value)/5
                for j in range(5):
                    box = plt.Rectangle((-0.5+i*2+1.95,11-t*2.5-1.1+j*0.5*0.8+0.1),0.05,0.01,fc='k',ec='k',lw=0.2,alpha=0.5)
                    ax.add_artist(box)
                    ax.text(-0.65+i*2+2.1 ,11-t*2.5-1+j*0.5*0.8,str(round(min_value+j*step,1)),fontsize=8,horizontalalignment='right',verticalalignment='center')
                # the len(groups) is 2, so we can calculate the significance level between the two groups
    
                from scipy.stats import ttest_ind
                t_stat, p_value = ttest_ind(bios_data_1, bios_data_2, equal_var=False)
                print(f' {bio_name} for {groups[0]} and {groups[1]}: t-statistic = {t_stat}, p-value = {p_value}')
                if p_value < 0.001:
                    ax.text(2*(i+0.6)-0.6*0-0.9,12-t*2.5-0.15,'***',fontsize=18,horizontalalignment='center',verticalalignment='center')
                elif p_value < 0.01:
                    ax.text(2*(i+0.6)-0.6*0-0.9,12-t*2.5-0.15,'**',fontsize=18,horizontalalignment='center',verticalalignment='center')
                elif p_value < 0.05:
                    ax.text(2*(i+0.6)-0.6*0-0.9,12-t*2.5-0.15,'*',fontsize=18,horizontalalignment='center',verticalalignment='center')
                else:
                    ax.text(2*(i+0.6)-0.6*0-0.9,12-t*2.5-0.1,'ns',fontsize=12,horizontalalignment='center',verticalalignment='center')
                
                # add the connecting line
                ax.plot([2*(i+0.6)-0.6*0-0.6,2*(i+0.6)-0.6*1-0.6],[12-t*2.5-0.2,12-t*2.5-0.2],color='k',lw=1,alpha=0.5)
                ax.plot([2*(i+0.6)-0.6*0-0.6,2*(i+0.6)-0.6*0-0.6],[12-t*2.5-0.3,12-t*2.5-0.2],color='k',lw=1,alpha=0.5)
                ax.plot([2*(i+0.6)-0.6*1-0.6,2*(i+0.6)-0.6*1-0.6],[12-t*2.5-0.3,12-t*2.5-0.2],color='k',lw=1,alpha=0.5)
                #ax.text(2*(i+0.6)-0.6*0-0.6,12-t*2.5-0.1,bio_name,fontsize=14,horizontalalignment='center',verticalalignment='center')
                #ax.text(2*(i+0.6)-0.6*0-0.6,9.5,bio_name,fontsize=14,horizontalalignment='center',verticalalignment='center')

ax.axis('off')

"""
BIO1 = Annual Mean Temperature
BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
BIO3 = Isothermality (BIO2/BIO7) (×100)
BIO4 = Temperature Seasonality (standard deviation ×100)
BIO5 = Max Temperature of Warmest Month
BIO6 = Min Temperature of Coldest Month
BIO7 = Temperature Annual Range (BIO5-BIO6)
BIO8 = Mean Temperature of Wettest Quarter
BIO9 = Mean Temperature of Driest Quarter
BIO10 = Mean Temperature of Warmest Quarter
BIO11 = Mean Temperature of Coldest Quarter
BIO12 = Annual Precipitation
BIO13 = Precipitation of Wettest Month
BIO14 = Precipitation of Driest Month
BIO15 = Precipitation Seasonality (Coefficient of Variation)
BIO16 = Precipitation of Wettest Quarter
BIO17 = Precipitation of Driest Quarter
BIO18 = Precipitation of Warmest Quarter
BIO19 = Precipitation of Coldest Quarter         
"""

# add the annotation for the bios
plt.text(3.8, 4.4, 'bio1: Annual Mean Temperature', fontsize=10)
plt.text(3.8, 4.0, 'bio2: Mean Diurnal Range', fontsize=10)
plt.text(3.8, 3.6, 'bio3: Isothermality', fontsize=10)
plt.text(3.8, 3.2, 'bio4: Temperature Seasonality', fontsize=10)
plt.text(3.8, 2.8, 'bio5: Max Temperature of Warmest Month', fontsize=10)
plt.text(3.8, 2.4, 'bio6: Min Temperature of Coldest Month', fontsize=10)
plt.text(7.8, 5.6, 'bio7: Temperature Annual Range', fontsize=10)
plt.text(7.8, 5.2, 'bio8: Mean Temperature of Wettest Quarter', fontsize=10)
plt.text(7.8, 4.8, 'bio9: Mean Temperature of Driest Quarter', fontsize=10)
plt.text(7.8, 4.4, 'bio10: Mean Temperature of Warmest Quarter', fontsize=10)
plt.text(7.8, 4.0, 'bio11: Mean Temperature of Coldest Quarter', fontsize=10)
plt.text(7.8, 3.6, 'bio12: Annual Precipitation', fontsize=10)
plt.text(7.8, 3.2, 'bio13: Precipitation of Wettest Month', fontsize=10)
plt.text(7.8, 2.8, 'bio14: Precipitation of Driest Month', fontsize=10)
plt.text(7.8, 2.4, 'bio15: Precipitation Seasonality', fontsize=10)
plt.text(12, 3.6, 'bio16: Precipitation of Wettest Quarter', fontsize=10)
plt.text(12, 3.2, 'bio17: Precipitation of Driest Quarter', fontsize=10)
plt.text(12, 2.8, 'bio18: Precipitation of Warmest Quarter', fontsize=10)
plt.text(12, 2.4, 'bio19: Precipitation of Coldest Quarter', fontsize=10)

# add the significance level
plt.text(7.8, 7.0, 'Significance Level', fontsize=12, fontweight='bold')
plt.text(7.8, 6.6, '*** p < 0.001', fontsize=12)
plt.text(9.4, 6.6, '** p < 0.01', fontsize=12)
plt.text(7.8, 6.2, '* p < 0.05', fontsize=12)
plt.text(9.4, 6.2, 'ns p > 0.05', fontsize=12)
plt.savefig('8pop_bio2_2.png',dpi=300, bbox_inches='tight')

features = ['bio2', 'bio3', 'bio4', 'bio15', 'bio1', 'bio14', 'bio12','bio8']

# write the data to a csv file
df.to_csv('PCA_East_Asia_climate_TOP_3_5pops.csv', index=False)
# standardize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df[features])
#print(df_scaled)
# apply PCA
pca = PCA()
# fit PCA with the data
components = pca.fit_transform(df[features])
#print(components)
# Get the feature importance
importance = np.abs(pca.components_[0])


# Sort the feature labels based on the importance
feature_names = np.array(features)
sorted_idx = np.argsort(importance)
sorted_feature_names = feature_names[sorted_idx]

# Create labels with PC number, variance and original variable name
labels = {
    str(i): f"PC {i+1} ({var:.2f}%) - {sorted_feature_names[i]}"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
}
print(labels)
#print(components)
# here we want plot the pca result with matplotlib
import matplotlib.pyplot as plt
# plot the pca result with matplotlib

colors = {'N':'#00B0F0','S':'#68AC56'}

# create the lagend plot for the pca result
fig, ax = plt.subplots()
fig.set_size_inches(2.5,2.5)
ax.set_xlim(0.2,7.61)
ax.set_ylim(8.2,15.61)
# set the  font in arail bold
plt.rcParams['font.sans-serif'] = ['Arial']
# plt.rcParams['font.weight'] = 'bold'
# aff the ticks of the subplots
ax.set_xticks([])
ax.set_yticks([])
# off the margin of the subplots
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
# add the pca result
for i in range(3):
    for j in range(3):
        if i < j: 
        # in the origial plot, we have 19*19 subplots, so we need to change the position of the subplots
            standard_i = max(components.__abs__()[:,i])
            #print(standard_i)
            standard_j = max(components.__abs__()[:,j])
            print(standard_j)
            # according to the standard_i and standard_j, we can change the position of the subplots into this plot
            new_i = (14-4*i)+ 1.5*components[:,i]/standard_i
            new_j = (2+4*(j-i-1))+ 1.5*components[:,j]/standard_j
            print(new_i)
            print(new_j)
            #print(df['group'])
            ax.scatter(new_j, new_i, c=df['group'].apply(lambda x: colors[x]), s=5)
            # add box
            box = plt.Rectangle((0.4+4*(j-i-1),12.4-4*i),3.2,3.2,fc='none',ec='k',lw=1)
            ax.add_artist(box)
            # add the 0 line
            line = plt.Line2D([0.4+4*(j-i-1),3.6+4*(j-i-1)],[14-4*i,14-4*i],color='k',linestyle='--',linewidth=0.5,alpha=0.4)
            ax.add_artist(line)
            line = plt.Line2D([2+4*(j-i-1),2+4*(j-i-1)],[12.4-4*i,15.6-4*i],color='k',linestyle='--',linewidth=0.5,alpha=0.4)
            ax.add_artist(line)
            # add the label of the pca
            ax.text(2+4*(j-i-1),12.2-4*i-0.1,labels[str(j)],fontsize=6,horizontalalignment='center',verticalalignment='center')
            ax.text(0.2+4*(j-i-1)-0.1,14-4*i,labels[str(i)],fontsize=6,horizontalalignment='center',verticalalignment='center',rotation=90)

#box = plt.Rectangle((4.4,8.4),3.2,3.2,fc='none',ec='k',lw=1)
#ax.add_artist(box)
for i in range(2):
    if i%2 == 0:
        circle = plt.Circle((4.6,10.3+(i//2)*0.8), 0.1, color=list(colors.values())[i])
        ax.add_artist(circle)
        ax.text(4.45,11+(i//2)*0.8,"North",fontsize=6,verticalalignment='center',rotation=90)
    else:
        circle = plt.Circle((4.6,9.0+(i//2)*0.8), 0.1, color=list(colors.values())[i])
        ax.add_artist(circle)
        ax.text(4.45,9.68+(i//2)*0.8,"South",fontsize=6,verticalalignment='center',rotation=90)
 
# save the figure
plt.savefig('PCA_East_Asia_climate_TOP_3_5_2pops.png', dpi=600, bbox_inches='tight')
