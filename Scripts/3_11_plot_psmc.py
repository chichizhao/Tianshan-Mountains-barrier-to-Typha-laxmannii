
#!/bin/bash
# function: plot the psmc result and niche in one figure

# Path: data/plot_psmc/script/plot_psmc.py
# function: plot psmc result with chichi's understanding
# reference: easy_psmc_plot.py
# usage python plot_psmc.py -namelist popmap.txt -psmcdir psmc_result -out psmc_plot

import sys
import os
import argparse
import click
import matplotlib.pyplot as plt
import numpy as np

# we need set some parameters 
mutation_rate = 5.13e-9 # we estimate from the whole genome data
generation_time = 1 # we set it as 1 year for one generation
size = 100 # we set the bin size as 100 in default in psmc analysis
def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-namelist', '--namelist', type=str, help='the namelist of the samples')
    parser.add_argument('-psmcdir', '--psmcdir', type=str, help='the directory of the psmc result')
    parser.add_argument('-out', '--out', type=str, help='the output directory of the psmc plot')
    return parser



# here is the idea of the psmc plot
# we want to plot three regions of the psmc result in one figure 
# and mark the three regions with three different colors 
# and for the samples, from the same site, we may want to get the average plot of this site
# while the above may not that true, it is just a initial idea now.
def read_name_and_group(filename):
    name_group_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            name = line.split('\t')[0]
            region = line.split('\t')[1]
            print(name, region)
            name_group_dict[name] = [region]
    return name_group_dict



# the following function is get the psmc result from the psmcfa file
def psmc_fun(filename, size, mutation, generation_time):
    # Read the raw file
    psmc_rlt = open(filename, "r")
    result = psmc_rlt.read()
    psmc_rlt.close()
    # Getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    #print(last_block)
    # this is the last block of psmc result. for example, the last block of psmc result is like this: the rpund of 50 
    last_block = last_block.split('\n')
    #print(last_block)
    # convert the last block to a list
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2] =="RS":
            #print(line)
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))
    #print(time_windows)
    #print(estimated_lambdas)        
    # based on the result of the time_windows and estimated_lambdas, get the information of the last round of psmc
    '''
    rs period   time    estimated_lambda  
    RS	0	0.000000	1407962.503178	0.000456	0.000000	0.000000
    '''
    # Getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the estimated parameters
    #print(result)
    result = result[-1].split('\n')[0]
    #print(result)
    # comments by chichi 20230-09-14, this result use the last round of psmc
    # it is the last pa line, so it is the last round of psmc
    result = result.split(' ')
    #print(result)
    theta = float(result[1])
    #print(theta)
    N0 = theta/(4*mutation)/size
    # here the size is the bin size, here it is 100 in default
    # mutation is the mutation rate we difined
    # theta is pa line in the last block of psmc result, and for the NO, the theta is the first element of the pa line 
    # Scalling times and sizes
    # the method one we use to estimate the time and size, we know the mutation rate.
    times = [generation_time * 2 * N0 * i for i in time_windows]
    # set the scale of the last round of psmc, here the generation time is 1.
    sizes = [N0 * i for i in estimated_lambdas]
    # we get the real population size that we estimate.
    #print(times)
    #print(sizes)
    # Remove the false positive result
    raw_dict = {}
    false_result = sizes[-1]
    #print(false_result)
    for i in range(len(sizes)):
        raw_dict[times[i]] = sizes[i]
    times = []
    sizes = []
    for k, v in raw_dict.items():
        if str(v) != str(false_result):
            times.append(k)
            sizes.append(v)
        else:
            break
    #print(times)
    #print(sizes)
    return(times, sizes)

def plot_result(result_dict, out):
    fig, ax = plt.subplots(figsize=(18, 8), dpi=600)
    # set font style in arail
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.xscale('log')
    plt.xlabel('Years before present', fontsize=24)
    plt.ylabel('Effective population size (x 100000)', fontsize=24)
    plt.xlim(500, 4e6)
    plt.ylim(0, 5.5e5)

    
    # add the climate date of last  800 kya
    import pandas as pd
    import numpy as np

    data = pd.read_csv('/home/chichi/data/china/china3/momi/EDC_dD_temp_estim.tab', sep='\t', header=0)
    print(data.head())
    # remove the rows with missing values
    data = data.dropna()
    print(data.head())
    age = data['Age model [ka]']*1000
    T = data['delta T [°C]'] 
    # turn in to exponential scale in 10
    T = T*7000 + 420000
    print(T)
    line = plt.plot(age, T, color='black', linewidth=1,alpha=0.5)
    ax = plt.gca().add_line(line[0])
    

    box = plt.Rectangle((800000, 340000), 10000, 120000, facecolor='black', edgecolor='black', linewidth=2)
    plt.gca().add_patch(box)
    for i in range(1,5):
        martT =  (i*5-15)*7000+ 420000
        Rectangle = plt.Rectangle((800000, martT), 50000, 1e2, facecolor='none', edgecolor='black', linewidth=1)
        plt.gca().add_patch(Rectangle)
        plt.text(1000000, martT, str(i*5-15), fontsize=14, va='center', ha='center')
    
    plt.text(1.4e6, 420000, 'Antarctic Temperature (°C) \n relative to present', fontsize=14, va='center', ha='center', rotation=90)
    # add the line of 0 degree
    line = plt.plot([500,8e5], [420000, 420000], color='black', linewidth=1,alpha=0.5, linestyle='--')
    plt.gca().add_line(line[0])

    # here we need set the color of the regions
    colors = { 'N':'#34587F', 'S':'#68AC56'}
    N = {}
    S = {}
    N_times = []
    S_times = []

    for name in result_dict.keys():
        region = result_dict[name][0]
        times = result_dict[name][1]
        sizes = result_dict[name][2]
        # here we plot the psmc result with the times and sizes, and it is steps plot
        #print(times)
        #print(sizes)
        #print(region, len(times), len(sizes))
        sizes = [i/2.5 for i in sizes if i > 250000]+ [i for i in sizes if i <= 250000]
        #plt.step(times, sizes, color=colors[region], alpha=0.5, label=region,linewidth=0.2)
        
        if region == 'N':
            N[name] = [times, sizes]
            N_times = N_times + times
        elif region == 'S':
            S[name] = [times, sizes]
            S_times = S_times + times

        # here we need modify the every element in sizes, if the sizes over 200000, we divide it in to 2
        #print(region, len(times), len(sizes))
        #sizes = [i/2 for i in sizes if i > 200000]+ [i for i in sizes if i <= 200000]
        plt.step(times, sizes, color=colors[region], alpha=0.4, label=region,linewidth=0.2)
        
        # here we need to add the region of the samples
        #times_over_10000 = [t for t in times if t > 10000]
        #sizes_over_10000 = [sizes[i] for i in range(len(times)) if times[i] > 10000]
        #sizes_over_10000 = [i *4 + 300000 for i in sizes_over_10000]
        #plt.step(times_over_10000, sizes_over_10000, color=colors[region], alpha=0.5, linewidth=0.2)

    # sort the N times
    N_times = sorted(N_times)
    # so we have the time points of the N region
    # so we can calculate the each time region of the average size
    N_average_size = calculate_average_zise(N, N_times)
    #print(average_size)
    N_times, N_average_size = smooth_carve(N_times, N_average_size)
    # add the average steps plot of the N region
    plt.step(N_times, N_average_size, color=colors['N'], label='N', linewidth=5, alpha=0.7)    

    # sort the S times
    S_times = sorted(S_times)
    S_average_size = calculate_average_zise(S, S_times)
    S_times, S_average_size = smooth_carve(S_times, S_average_size)
    S_times= S_times[:-5]
    S_average_size = S_average_size[:-5]

    plt.step(S_times, S_average_size, color=colors['S'], label='S', linewidth=5, alpha=0.7)

  

    # set x axis and y axis font size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # change the axis mark inside the plot
    plt.tick_params(axis='x', which='major', direction='in', length=8, width=1.5)
    plt.tick_params(axis='y', which='major', direction='in', length=8, width=1.5)
    # off the minor tick
    plt.tick_params(axis='x', which='minor', bottom=False, top=False)
    plt.tick_params(axis='y', which='minor', left=False, right=False)
    # set axis line width
    plt.rcParams['axes.linewidth'] = 1.5
    
    # add the lengend
    box = plt.Rectangle((6e5, 5e4), 4.5e5, 1.5e4, facecolor=colors['N'])
    plt.text(11e5, 5.2e4, 'North', fontsize=14)
    plt.gca().add_patch(box)
    box = plt.Rectangle((6e5, 2e4), 4.5e5, 1.5e4, facecolor=colors['S'])
    plt.text(11e5, 2.2e4, 'South', fontsize=14)
    plt.gca().add_patch(box)

    
######################################################
    # add three important time region# 
    # 6000 years ago : Mid-Holocene
    # add the region of the mid-holocene from 5000 to 7000 years ago
    #box = plt.Rectangle((5e3, 1e2), 2e3, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2) 
    #plt.gca().add_patch(box)
    arrow = plt.arrow(6e3, 1e2, 0, 410000-2e4, head_width=1e3, head_length=2e4, fc='grey', ec='grey',width=500,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(6e3, 3.4e5, '~6 Kya', fontsize=14, ha='right')
    plt.text(6e3, 3.2e5, 'Mid-Holocene', fontsize=14, ha='right')

    # here we need add the time line from 12 kya to 14 kya
    # 13000 years ago : Younger Dryas (YD; ~12,900 - 11,700 years BP)
    # add the region of the younger dryas from 12000 to 14000 years ago
    #box = plt.Rectangle((1.2e4, 1e2), 2e3, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
    
    # Last glacial maximum (LGM; ~21,000 years BP)
    # add the region of the last glacial maximum from 20000 to 22000 years ago
    #box = plt.Rectangle((1.8e4, 1e2), 4e3, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
    
    arrow = plt.arrow(2.0e4,0, 0, 340000-2e4, head_width=4e3, head_length=2e4, fc='grey', ec='grey',width=2000,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(2.1e4, 3.0e5, '~21 Kya', fontsize=14, ha='left')
    plt.text(2.1e4, 2.8e5, 'LGM', fontsize=14, ha='left')
    
    # 130000 years ago : Last inter-glacial (LIG; ~120,000 - 140,000 years BP)
    # add the region of the last inter-glacial from 120000 to 140000 years ago
    arrow = plt.arrow(1.3e5,1000, 0, 340000-2e4, head_width=2e4, head_length=1e4, fc='grey', ec='grey',width=10000,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(1.4e5, 3.0e5, '~130 Kya', fontsize=14, ha='left')
    plt.text(1.4e5, 2.8e5, 'LIG', fontsize=14, ha='left')
    #box = plt.Rectangle((1.2e5, 1e2), 2e4, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
        # add the cited date
    arrow = plt.arrow(2.78e5, 1000, 0, 330000-1e2, head_width=4e4, head_length=2e4, fc='grey', ec='grey',width=20000,alpha=0.1)
    plt.gca().add_patch(arrow)
    plt.text(2.8e5, 2.8e5, '~280 Kya', fontsize=14, ha='left')
    plt.text(2.8e5, 2.6e5, 'MIS 8', fontsize=14, ha='left')
    

    plt.text(4e4, 440000, 'Jouzel et al., 2007', fontsize=14, ha='center')
    


    # 400 000 years ago : mark region
    #box = plt.Rectangle((4e5, 1e2), 1e5, 3e6-1e2, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.2)
    #plt.gca().add_patch(box)
######################################################
    # set the y mark
    plt.yticks([50000, 100000, 150000,200000,250000,330000,420000,500000,580000], ['0.5', '1.0', '1.5','2.0','2.5','4.5','6.5','8.5','10.5'], fontsize=14)

##########################################################
    plt.savefig("psmc_8pop7.png", dpi=600, bbox_inches='tight')
    plt.close()

    ## here we need plot this figure in normal scale
    fig, ax = plt.subplots(figsize=(13, 11.3), dpi=600)
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.xlabel('Years before present', fontsize=16)
    plt.ylabel('Effective population size (x 100000)', fontsize=16)
    plt.xlim(10, 140)
    plt.ylim(22, 135)
    box = plt.Rectangle((13, 51.2), 38, 5, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.3)
    plt.gca().add_patch(box)

    trapezoid = plt.Polygon([[13, 56.2], [51, 56.2], [84, 62], [64, 62]], closed=True, facecolor='grey', edgecolor='none', alpha=0.3)
    plt.gca().add_patch(trapezoid)

    box = plt.Rectangle((55, 51.2), 38, 5, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.3)
    plt.gca().add_patch(box)

    trapezoid = plt.Polygon([[55, 56.2], [93, 56.2], [93.7, 62], [85, 62]], closed=True, facecolor='grey', edgecolor='none', alpha=0.3)
    plt.gca().add_patch(trapezoid)

    box = plt.Rectangle((97, 51.2), 38, 5, facecolor='grey', edgecolor='none', linewidth=2, alpha=0.3)
    plt.gca().add_patch(box)

    trapezoid = plt.Polygon([[97, 56.2], [135, 56.2], [109, 62], [95, 62]], closed=True, facecolor='grey', edgecolor='none', alpha=0.3)
    plt.gca().add_patch(trapezoid)
    # add the photo of the 130kya
    photo_130kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/presence_130kya.png')
    # and we need cut the photo to the special region
    photo_130kya = photo_130kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_130kya, extent=(95, 135, 104, 134), aspect='auto')

    # add the photo of the 21kya
    photo_21kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/presence_21kya.png')
    # and we need cut the photo to the special region
    photo_21kya = photo_21kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_21kya, extent=(53, 93, 104, 134), aspect='auto')

    # add the photo of the 6kya
    photo_6kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/presence_6kya.png')
    # and we need cut the photo to the special region
    photo_6kya = photo_6kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_6kya, extent=(11, 51, 104, 134), aspect='auto')

    # add the photo of the present
    photo_present = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/presence_present.png')
    # and we need cut the photo to the special region
    photo_present = photo_present[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_present, extent=(10.5, 58.5, 61, 97), aspect='auto')

    # add the photo of the present_minus_6kya 
    photo_present_minus_6kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/presence_minus_6kya.png')
    # and we need cut the photo to the special region
    photo_present_minus_6kya = photo_present_minus_6kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_present_minus_6kya, extent=(11, 51, 23, 53), aspect='auto')

    # add the photo of the 6kya_minus_21kya
    photo_6kya_minus_21kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/6kya_minus_21kya.png')
    # and we need cut the photo to the special region
    photo_6kya_minus_21kya = photo_6kya_minus_21kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_6kya_minus_21kya, extent=(53, 93, 23, 53), aspect='auto')

    # add the photo of the 21kya_minus_130kya
    photo_21kya_minus_130kya = plt.imread('/home/chichi/data2/chichi/typha_lax/psmc/21kya_minus_130kya.png')
    # and we need cut the photo to the special region
    photo_21kya_minus_130kya = photo_21kya_minus_130kya[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_21kya_minus_130kya, extent=(95, 135, 23, 53), aspect='auto')

    # add the photo of the psmc_8pop7.png
    photo_psmc_8pop7 = plt.imread('psmc_8pop7.png')
    # and we need cut the photo to the special region
    # photo_psmc_8pop7 = photo_psmc_8pop7[0:3162, 400:4562]
    # add the photo to the box
    ax.imshow(photo_psmc_8pop7, extent=(60, 138, 58, 96), aspect='auto')
    # off the x and y axis
    plt.axis('off')

    # here we need add the time line to the figure

    box = plt.Rectangle((22, 99), 83,1.5, facecolor='blue', edgecolor='none', linewidth=2, alpha=0.5)
    plt.gca().add_patch(box)

    box = plt.Rectangle((107, 99), 5, 1.5, facecolor='blue', edgecolor='none', linewidth=2, alpha=0.3)
    plt.gca().add_patch(box)

    box = plt.Rectangle((114, 99), 6, 1.5, facecolor='blue', edgecolor='none', linewidth=2, alpha=0.3)
    plt.gca().add_patch(box)
    box = plt.Rectangle((122, 99), 7, 1.5, facecolor='blue', edgecolor='none', linewidth=2, alpha=0.1)
    plt.gca().add_patch(box)

    # add the end arrow
    arrow = plt.arrow(21, 100.5, 0, -2.5, width=2, head_width=5, head_length=2, fc='blue', ec='none',alpha=0.3)
    plt.gca().add_patch(arrow)

    # add the time line text marks
    plt.text(25,59, 'Present', fontsize=14, va='center', ha='center')
    plt.text(30, 102,'Mid-Holocene ~ 6 Kya', fontsize=14, va='center', ha='center')
    plt.text(73, 102, 'Last Glacial Maximum ~ 21 Kya', fontsize=14, va='center', ha='center')
    plt.text(115, 102, 'Last Inter-glacial ~ 130 Kya', fontsize=14, va='center', ha='center')

    plt.text(30,54.5, 'Present <= Mid-Holocene', fontsize=14, va='center', ha='center')
    plt.text(73, 54.5, 'Mid-Holocene <= LGM', fontsize=14, va='center', ha='center')
    plt.text(115, 54.5, 'LGM <= LIG', fontsize=14, va='center', ha='center')
    
    # add the legend of Precence probability

    color_set= ['#4a9c80','#88d190','#d8f3c1','#e6f6df']
    box = plt.Rectangle((136, 105), 1, 7, fc="#4a9c80",ec='black',linewidth=0.2)
    ax.add_patch(box)
    box = plt.Rectangle((136, 112), 1, 7, fc="#88d190",ec='black',linewidth=0.2)
    ax.add_patch(box)
    box = plt.Rectangle((136, 119), 1, 7, fc="#d8f3c1",ec='black',linewidth=0.2)
    ax.add_patch(box)
    box = plt.Rectangle((136, 126), 1, 7, fc="white",ec='black',linewidth=0.2)
    ax.add_patch(box)

    plt.text(145, 118, 'Presence probability', fontsize=14, va='center', ha='center',rotation=90)
    plt.text(137.5, 105, '1.00', fontsize=14, va='center', ha='left')
    plt.text(137.5, 112, '0.75', fontsize=14, va='center', ha='left')
    plt.text(137.5, 119, '0.50', fontsize=14, va='center', ha='left')
    plt.text(137.5, 126, '0.25', fontsize=14, va='center', ha='left')
    plt.text(137.5, 133, '0.00', fontsize=14, va='center', ha='left')
    # here we need add the color set of the psmc result, it is the size of the population
    
    # and we need add the lable to show the in crease  and the decrease
    color_set=['#ff0303','#ff7a46','#ffcc8d','#ffcfcf','#ffffd4','#d4ffd4','#f9fad5','#d2ed8e','#8fd147','#3aa903']
    box = plt.Rectangle((136, 23), 1, 8, fc="#ff0303",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31), 1, 0.2/0.5*8, fc="#ff7a46",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8), 1, 0.15/0.5*8, fc="#ffcc8d",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8), 1, 0.10/0.5*8, fc="#ffcfcf",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8+0.10/0.5*8), 1, 0.10/0.5*8, fc="#ffffd4",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8+0.20/0.5*8), 1, 0.10/0.5*8, fc="#d2ed8e",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8), 1, 0.15/0.5*8, fc="#d2ed8e",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8), 1, 0.2/0.5*8, fc="#8fd147",ec='black',linewidth=0.5)
    ax.add_patch(box)
    box = plt.Rectangle((136, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.2/0.5*8), 1, 0.5/0.5*8, fc="#3aa903",ec='black',linewidth=0.5)
    ax.add_patch(box)

    plt.text(136, 21.5, 'Decrease', fontsize=14, va='center', ha='left')
    plt.text(137.5, 23, '1.00', fontsize=14, va='bottom', ha='left')
    plt.text(137.5, 31, '0.50', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8, '0.30', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8, '0.15', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8, '0', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8, '0.15', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8, '0.30', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.2/0.5*8, '0.50', fontsize=14, va='center', ha='left')
    plt.text(137.5, 31+0.2/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.15/0.5*8+0.2/0.5*8+0.5/0.5*8, '1.00', fontsize=14, va='top', ha='left')
    plt.text(136, 57.5, 'Increase', fontsize=14, va='center', ha='left')

    plt.text(145,40, 'Presence probability change', fontsize=14, va='center', ha='center', rotation=90)

    plt.savefig("psmc_8pop8.png", dpi=150, bbox_inches='tight')
    plt.close()

def calculate_average_zise(dict, times):
    # base one each some samples in the dict has it period of the psmc result.
    # the data has such character:
    # 1. it goes with step plot with corresponding times and sizes
    # 2. each samples's time line points are not the same
    # the solutions is gather all the time points of the samples, sort them.
    # from the very first time point, we calculate the average size of the this first time period.
    # and keep the thoughts going on. for the end of this study line, it may not that all the samples reach the end of the time line.
    # so we depends on how many sample go there and calculate the average size of this period to represent the average size of this period, also this population.
    # we get the sorted time points line: times
    # we get the dict of the samples: dict
    summary_times = times
    average_size = [0] * len(summary_times)
    average_sample_number = [0]* len(summary_times)
    #print(len(summary_times))
    #print(len(average_size))
    for sample in dict.keys():
        sample_times = dict[sample][0]
        sample_sizes = dict[sample][1]
        j = 0
        k = 0
        for i in range(len(summary_times)):
            #print(i)
            #print(j)
            if summary_times[i] == sample_times[j]:
                for t in range(k, i+1):
                    average_size[t] = average_size[t] + sample_sizes[j]
                    average_sample_number[t] = average_sample_number[t] + 1
                k = i
                j = j + 1
                if j == len(sample_times):
                    break
            elif summary_times[i] > sample_times[j]:
                break
    #print(average_size)
    #print(average_sample_number)
    average_size = [average_size[i]/average_sample_number[i] for i in range(len(average_size))]
    return average_size
def smooth_carve(time, average_zise):
    # here we want to smooth the carve of the psmc result
    # here is my solution: gather 4 time points,in to one and give the new time point the average size of the 4 time points
    step = 5
    smooth_time = []
    smooth_size = []
    for i in range(0, len(time), step):
        smooth_time.append(time[i])
        smooth_size.append(sum(average_zise[i:i+step])/step)
    return smooth_time, smooth_size

if __name__ == '__main__':
    args = get_argparser().parse_args()
    namelist = args.namelist
    psmcdir = args.psmcdir
    out = args.out
    name_group_dict = read_name_and_group(namelist)
    result_dict = {}
    for name in name_group_dict.keys():
        region = name_group_dict[name][0]
        psmcfile = os.path.join(psmcdir, name + '.psmc')
        print(psmcfile)
        times, sizes = psmc_fun(psmcfile, size, mutation_rate, generation_time)
        result_dict[name] = [region, times, sizes]
    #print(result_dict)
    plot_result(result_dict, out)
