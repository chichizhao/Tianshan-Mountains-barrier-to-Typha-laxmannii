#!/bin/python
# function: plot momi results
# author: chichi
# function, plot the momi results and psmc results

#!/bin/python
# -*- coding: utf-8 -*-
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
    fig = plt.figure(figsize=(12, 10))
    # set font style in arail
    plt.rcParams['font.sans-serif'] = ['Arial']
    
    # add the climate date of last  800 kya
    import pandas as pd
    import numpy as np


    # here we need set the color of the regions
    colors = { 'N':'#34587F', 'S':'#68AC56', 'N1':'#9BB2E2','N2':'#BACEC9','N3':'#9EAFC5', 'N4':'#34587F','EA':'#4A7DB4','MA':'#68AC56'}

    N = {}
    S = {}
    N1 = {}
    N2 = {}
    N3 = {}
    N4 = {}

    N_times = []
    S_times = []
    N1_times = []
    N2_times = []
    N3_times = []
    N4_times = []



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
        elif region == 'N1':
            N1[name] = [times, sizes]
            N1_times = N1_times + times
        elif region == 'N2':
            N2[name] = [times, sizes]
            N2_times = N2_times + times
        elif region == 'N3':
            N3[name] = [times, sizes]
            N3_times = N3_times + times
        elif region == 'N4':
            N4[name] = [times, sizes]
            N4_times = N4_times + times

        # here we need modify the every element in sizes, if the sizes over 200000, we divide it in to 2
        #print(region, len(times), len(sizes))
        #sizes = [i/2 for i in sizes if i > 200000]+ [i for i in sizes if i <= 200000]
        #plt.step(times, sizes, color=colors[region], alpha=0.4, label=region,linewidth=0.2)
        
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
    #plt.step(N_times, N_average_size, color=colors['N'], label='N', linewidth=5, alpha=0.7)    
    # here we need plot the some area of the N region which is the N_TIME OVER 10000
    #N_times_over_10000 = [t for t in N_times if t > 8200]
    #N_average_size_over_10000 = [N_average_size[i] for i in range(len(N_times)) if N_times[i] >8200]
    #print(len(N_times_over_10000), len(N_average_size_over_10000))
    # ZOOM N_average_size_over_10000 * 10 + 1000
    #N_average_size_over_10000 = [i *4 + 300000 for i in N_average_size_over_10000]

    # add the area of the N region which is the N_TIME OVER 10000
    #plt.step(N_times_over_10000, N_average_size_over_10000, color=colors['N'], linewidth=5, alpha=0.7)
    # plot the N region with the area of the N region which is the N_TIME below 8300
    #N_times_below_8300 = [t for t in N_times if t <= 8200]
    #N_average_size_below_8300 = [N_average_size[i] for i in range(len(N_times)) if N_times[i] <= 8200]
    #N_average_size_below_8300 = [i /2.5 for i in N_average_size_below_8300]
    # add the area of the N region which is the N_TIME below 8300
    #plt.step(N_times_below_8300, N_average_size_below_8300, color=colors['N'], linewidth=5, alpha=0.7)

    # sort the S times
    S_times = sorted(S_times)
    S_average_size = calculate_average_zise(S, S_times)
    S_times, S_average_size = smooth_carve(S_times, S_average_size)
    S_times= S_times[:-5]
    S_average_size = S_average_size[:-5]


    """
    for i in range(1,5):
        martT =  (i*5-15)*20000+ 1500000
        Rectangle = plt.Rectangle((800000, martT), 50000, 1e2, facecolor='none', edgecolor='black', linewidth=1)
        plt.gca().add_patch(Rectangle)
        plt.text(1000000, martT, str(i*5-15), fontsize=14, va='center', ha='center')
    
    plt.text(1.4e6, 1500000, 'Antarctic Temperature (°C) \n relative to present', fontsize=14, va='center', ha='center', rotation=90)
    # add the line of 0 degree
    line = plt.plot([500,8e5], [1500000, 1500000], color='black', linewidth=1,alpha=0.5, linestyle='--')
    plt.gca().add_line(line[0])
    """




    
    # sort the N1 times
    N1_times = sorted(N1_times)
    N1_average_size = calculate_average_zise(N1, N1_times)
    N1_times, N1_average_size = smooth_carve(N1_times, N1_average_size)
    #plt.step(N1_times, N1_average_size, color=colors['N1'], label='N1', linewidth=5, alpha=0.7)
    
    # sort the N2 times
    N2_times = sorted(N2_times)
    N2_average_size = calculate_average_zise(N2, N2_times)
    N2_times, N2_average_size = smooth_carve(N2_times, N2_average_size)
    #plt.step(N2_times, N2_average_size, color=colors['N2'], label='N2', linewidth=5, alpha=0.7)
    
    # sort the N3 times
    N3_times = sorted(N3_times)
    N3_average_size = calculate_average_zise(N3, N3_times)
    N3_times, N3_average_size = smooth_carve(N3_times, N3_average_size)
    #plt.step(N3_times, N3_average_size, color=colors['N3'], label='N3', linewidth=5, alpha=0.7)
    
    # sort the N4 times
    N4_times = sorted(N4_times)
    N4_average_size = calculate_average_zise(N4, N4_times)
    N4_times, N4_average_size = smooth_carve(N4_times, N4_average_size)
    #plt.step(N4_times, N4_average_size, color=colors['N4'], label='N4', linewidth=5, alpha=0.7)
    from math import log
    # plot the S
    # the S and N2 split time is 64000
    for i in range(len(S_times)-1):
        if S_times[i+1] < 64000:   
            #box = plt.Rectangle((100000- S_average_size[i]/2, S_times[i]),S_average_size[i], S_times[i+1]-S_times[i], facecolor=colors['S'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((30- S_average_size[i]/20000/2, S_times[i]),S_average_size[i]/20000, S_times[i+1]-S_times[i], facecolor=colors['S'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)
    # the N2 and N1 split time is 73500
    for i in range(len(N2_times)-1):
        if N2_times[i+1] < 64000:   
            #box = plt.Rectangle((400000- N2_average_size[i]/2, N2_times[i]),N2_average_size[i], N2_times[i+1]-N2_times[i], facecolor=colors['N2'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((50- N2_average_size[i]/20000/2, N2_times[i]),N2_average_size[i]/20000, N2_times[i+1]-N2_times[i], facecolor=colors['N2'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)
    
    # plot the 64000 to the 73500 time period
    for i in range(len(N2_times)-1):
        if N2_times[i+1] >= 64000 and N2_times[i] < 72000:   
            #box = plt.Rectangle((400000- N2_average_size[i]/2, N2_times[i]),N2_average_size[i], N2_times[i+1]-N2_times[i], facecolor=colors['N2'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((40- N2_average_size[i]/20000/2, N2_times[i]),N2_average_size[i]/20000, N2_times[i+1]-N2_times[i], facecolor=colors['N2'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)
    # add the connection line between S and N2
    box = plt.Rectangle((28, 64000),26, 100, facecolor='none', edgecolor='black', linewidth=1, alpha=0.7)
    plt.gca().add_patch(box)
    plt.text(40, 56000, "64 kya", fontsize=16, va='center', ha='center', color='black')

    # the N1 and N3 split time is 37600

    for i in range(len(N3_times)-1):
        if N3_times[i+1] < 37600:   
            #box = plt.Rectangle((600000- N3_average_size[i]/2, N3_times[i]),N3_average_size[i], N3_times[i+1]-N3_times[i], facecolor=colors['N3'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((90- N3_average_size[i]/20000/2, N3_times[i]),N3_average_size[i]/20000, N3_times[i+1]-N3_times[i], facecolor=colors['N3'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)

    # plot the N3 time period 37600
    for i in range(len(N1_times)-1):
        if N1_times[i+1] < 37600: 
            #box = plt.Rectangle((800000- N1_average_size[i]/2, N1_times[i]),N1_average_size[i], N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((130- N1_average_size[i]/20000/2, N1_times[i]),N1_average_size[i]/20000, N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)
    
    # add the connection line between N1 and N3
    box = plt.Rectangle((90, 37600),40, 100, facecolor='none', edgecolor='black', linewidth=1, alpha=0.7)
    plt.gca().add_patch(box)
    plt.text(110, 41600, "37.6 kya", fontsize=16, va='center', ha='center', color='black')


    # plot the time period  from 37600 to 73500
    for i in range(len(N1_times)-1):
        if N1_times[i+1] >= 37600 and N1_times[i] < 73500:   
            #box = plt.Rectangle((800000- N1_average_size[i]/2, N1_times[i]),N1_average_size[i], N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((115- N1_average_size[i]/20000/2, N1_times[i]),N1_average_size[i]/20000, N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)

    # plot the time period  from 73500 to 85500
    for i in range(len(N1_times)-1):
        if N1_times[i+1] > 78000 and N1_times[i] < 85500:   
            #box = plt.Rectangle((800000- N1_average_size[i]/2, N1_times[i]),N1_average_size[i], N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((78- N1_average_size[i]/20000/2, N1_times[i]),N1_average_size[i]/20000, N1_times[i+1]-N1_times[i], facecolor=colors['N1'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)

    # add the connection line between N1 and N2
    box = plt.Rectangle((35, 77000),80,100, facecolor='black', edgecolor='black', linewidth=1, alpha=0.7)
    plt.gca().add_patch(box)
    plt.text(75, 68000, "77 kya", fontsize=16, va='center', ha='center', color='black')
    # the N1 and N4 split time is 85500
    for i in range(len(N4_times)-1):
        if N4_times[i+1] < 85500:   
            #box = plt.Rectangle((1100000- N4_average_size[i]/2, N4_times[i]),N4_average_size[i], N4_times[i+1]-N4_times[i], facecolor=colors['N4'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((160- N4_average_size[i]/20000/2, N4_times[i]),N4_average_size[i]/20000, N4_times[i+1]-N4_times[i], facecolor=colors['N4'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)

    # plot the N4 time period 85500
    for i in range(len(N4_times)-1):
        if N4_times[i+1] >= 85500 and N4_times[i] < 400000:
            #box = plt.Rectangle((1100000- N4_average_size[i]/2, N4_times[i]),N4_average_size[i], N4_times[i+1]-N4_times[i], facecolor=colors['N4'], edgecolor='none', alpha=0.7)
            #plt.gca().add_patch(box)
            box = plt.Rectangle((117- N4_average_size[i]/20000/2, N4_times[i]),N4_average_size[i]/20000, N4_times[i+1]-N4_times[i], facecolor=colors['N4'], edgecolor='none', alpha=0.7)
            plt.gca().add_patch(box)
    # add the connection line between N1 and N4
    box = plt.Rectangle((80, 85500),80, 100, facecolor='none', edgecolor='black', linewidth=1, alpha=0.7)
    plt.gca().add_patch(box)
    plt.text(120,91500, "85.5 kya", fontsize=16, va='center', ha='center', color='black')
    # add the arrow to mark the migration 
    # model.add_pulse_param("p_mig_N2_N1", lower=0.72, upper=0.735)
    # model.add_time_param("t_mig_N2_N1", lower=19000, upper=20000)
    arrow = plt.arrow(54.5, 19000, 67, 0, head_width=2000, head_length=3, fc='grey', ec='grey',width=1000,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(82,21000, "0.27", fontsize=16, va='center', ha='center', color='black')
    plt.text(82, 16000, "19 kya", fontsize=16, va='center', ha='center', color='black')
    #model.add_pulse_param("p_mig_N3_N1", lower=0.28, upper=0.30)
    #model.add_time_param("t_mig_N3_N1", lower=29000, upper=31000)

    arrow = plt.arrow(124, 30000, -24.5, 0, head_width=2000, head_length=3, fc='grey', ec='grey',width=1000,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(112, 32000, "0.28", fontsize=16, va='center', ha='center', color='black')
    plt.text(112, 26000, "30 kya", fontsize=16, va='center', ha='center', color='black')

    # model.add_pulse_param("p_mig_N4_N1", lower=0.78, upper=0.79)
    # model.add_time_param("t_mig_N4_N1", lower=13000, upper=14000)

    arrow = plt.arrow(156.5, 13000, -18, 0, head_width=2000, head_length=3, fc='grey', ec='grey',width=600,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(147, 14500, "0.22", fontsize=16, va='center', ha='center', color='black')
    plt.text(147, 11000, "13 kya", fontsize=16, va='center', ha='center', color='black')

    # model.add_pulse_param("p_mig_N4_N2", lower=0.71, upper=0.73)
    # model.add_time_param("t_mig_N4_N2", lower=6000, upper=7000)

    arrow = plt.arrow(120.5, 6000, -62, 0, head_width=1000, head_length=3, fc='grey', ec='grey',width=300,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(67, 7000, "0.28", fontsize=16, va='center', ha='center', color='black')
    plt.text(67, 5000, "6 kya", fontsize=16, va='center', ha='center', color='black')

    # model.add_pulse_param("p_mig_N4_N3", lower=0.78, upper=0.79)
    # model.add_time_param("t_mig_N4_N3", lower=14000, upper=15000)

    arrow = plt.arrow(156.5, 14000, -58, 0, head_width=2000, head_length=3, fc='grey', ec='grey',width=600,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(111, 16000, "0.22", fontsize=16, va='center', ha='center', color='black')
    plt.text(111, 12000, "14 kya", fontsize=16, va='center', ha='center', color='black')

    # model.add_pulse_param("p_mig_N3_N2", lower=0.53, upper=0.54)
    # model.add_time_param("t_mig_N3_N2", lower=25000, upper=26000)
    
    arrow = plt.arrow(84.5, 25000, -26, 0, head_width=2000, head_length=3, fc='grey', ec='grey',width=1000,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(68, 28000, "0.47", fontsize=16, va='center', ha='center', color='black')
    plt.text(68, 22000, "25 kya", fontsize=16, va='center', ha='center', color='black')


    # set the x and y limit
    plt.xlim(18, 180)
    plt.ylim(700, 650000)
    # set the y axis to log scale
    plt.yscale('log') 

    # add the white block
    box = plt.Rectangle((20, 100), 160, 900, facecolor='white', edgecolor='none')
    plt.gca().add_patch(box)
    # add the group name to the legend
    plt.text(30, 850,"South", fontsize=14, va='center', ha='center')
    plt.text(50, 850,"North_1", fontsize=14, va='center', ha='center')
    plt.text(90, 850,"North_3", fontsize=14, va='center', ha='center')
    plt.text(130, 850,"North_2", fontsize=14, va='center', ha='center')
    plt.text(160, 850,"North_4", fontsize=14, va='center', ha='center')

    # add the time mark lengend
    box = plt.Rectangle((20, 1000), 0.1,600000, facecolor='black', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    
    # add the lengend of the effective population size
    box = plt.Rectangle((50, 420000), -10000/20000, 10, facecolor='black', edgecolor='grey', linewidth=15)
    plt.gca().add_patch(box)
    plt.text(50, 510000, "Effective Population Size", fontsize=14, va='center', ha='center', color='black')
    plt.text(52, 410000, "10,000", fontsize=14, va='center', ha='left', color='black')
    box = plt.Rectangle((50, 320000), -100000/20000, 10, facecolor='blue', edgecolor='grey', linewidth=15)
    plt.gca().add_patch(box)
    plt.text(52, 310000, "100,000", fontsize=14, va='center', ha='left', color='black')

    box = plt.Rectangle((50, 220000), -500000/20000, 10, facecolor='green', edgecolor='grey', linewidth=15)
    plt.gca().add_patch(box)
    plt.text(52, 210000, "500,000", fontsize=14, va='center', ha='left', color='black')

    arrow = plt.arrow(50, 150000, -20, 0, head_width=20000, head_length=3, fc='grey', ec='grey',width=10000,alpha=1)
    plt.gca().add_patch(arrow)
    plt.text(40, 170000, "0.47", fontsize=16, va='center', ha='center', color='black')
    plt.text(40, 130000, "25 kya", fontsize=16, va='center', ha='center', color='black')

    plt.text( 51, 150000, ": Migration", fontsize=14, va='center', ha='left', color='black')

    # add the legend box 
    box = plt.Rectangle((24, 110000), 60, 480000, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)

    box = plt.Rectangle((74, 110000), 10, 480000, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)

    plt.text(79, 240000, "Legend", fontsize=18, va='center', ha='center', color='black',rotation=90)

    # add the box of time mark
    box = plt.Rectangle((20, 1000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 1000, "1000", fontsize=14, va='center', ha='right', color='black')

    box = plt.Rectangle((20, 5000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 5000, "5000", fontsize=14, va='center', ha='right', color='black')

    box = plt.Rectangle((20, 10000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 10000, "10000", fontsize=14, va='center', ha='right', color='black')

    box = plt.Rectangle((20, 50000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 50000, "50000", fontsize=14, va='center', ha='right', color='black')

    box = plt.Rectangle((20, 100000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 100000, "100000", fontsize=14, va='center', ha='right', color='black')

    box = plt.Rectangle((20, 500000), 2, 1, facecolor='none', edgecolor='black', linewidth=1)
    plt.gca().add_patch(box)
    plt.text(19, 500000, "500000", fontsize=14, va='center', ha='right', color='black')

    plt.text( 20, 680000, "Time (years)", fontsize=14, va='center', ha='center', color='black')


    plt.axis('off')   
    plt.savefig("psmc_8pop9.png", dpi=600, bbox_inches='tight')
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


