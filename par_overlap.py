# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np 
import pandas as pd 
import re
import math
import matplotlib.pyplot as plt
import random
from numpy import array
import itertools
cn_summary = []
event_list = []
for chr in sorted(cn_summary):
    cn_summary_ = cn_summary[chr]
    for each_summary in sorted(cn_summary_):
        start, end = each_summary.strip('.')
        event_list.append([chr, int(start), int(end)])

def check_overlap(e1, e2):
    start = 1
    end = 2
    flag = False
    if e1[0] == e2[0]:
        if e1[start] > e2[start] and e1[start] < e2[end] :
            flag = True
        elif e1[end] > e2[start] and e1[end] < e2[end]:
            flag = True
        elif e1[start]> e2[start] and e1[end] < e2[end]:
            flag = True
        elif e2[start]> e1[start] and e2[end] < e1[end]:
            flag = True
    return flag
            
    

def compute_score(events):
    count = 0
    n = len(events)
    for i in range(0, len(events)):
        for j in range(0, len(events)): 
            if i == j : 
                pass 
            else: 
                if check_overlap(events[i], events[j]):
                    count = count + 1
    return count/pow(n, 2)


def plot_distr(length, freq, amp):
    x = np.linspace(start = 0, stop = 1, num = length)
    y = np.zeros(length)
    for i in range(0, 30):
        temp = random.uniform(-math.pi, math.pi)
        temp2 = random.uniform(0, amp)
        temp3 = random.uniform(1,freq)
        for j in range(0, length):
            y[j] = y[j] + math.sin(x[j]*temp3 + temp)*temp2
    y = np.exp(y)
    plt.plot(x,y)
    return x,y


#event.append()
def check_overlap_test(e1, e2):
    start = 0
    end = 1
    flag = False
    if e1[start] > e2[start] and e1[start] < e2[end] :
        flag = True
    elif e1[end] > e2[start] and e1[end] < e2[end]:
        flag = True
    elif e1[start]> e2[start] and e1[end] < e2[end]:
        flag = True
    elif e2[start]> e1[start] and e2[end] < e1[end]:
        flag = True
    return flag

def compute_score(events):
    count = 0
    n = len(events)
    for i in range(0, len(events)):
        for j in range(0, len(events)):
            if i == j :
                pass
            else:
                if check_overlap_test(events[i], events[j]) is True:
                    count = count + 1
                    break
                    #print("here")

    print("overlap count is: ", count)
    print("number of events: ", n)
    score = count*1.0/n*1.0
    print("percet of events overlap:", score)
    #score = count*1.0/pow(n-1,2)
    return score

def overlap_cover_chr(events):
    ind_list = []
    for event in events:
        ind_list.append(event[0])
        ind_list.append(event[1])
    ind_list = sorted(list(set(ind_list)))
    #print(ind_list)
    #ind_list = list(set(ind_list))
    #print(ind_list)
    seg_count = [0] * (len(ind_list)-1)
    for event in events:
        s_index = ind_list.index(event[0])
        e_index = ind_list.index(event[1])
        for i in range(s_index, e_index):
            seg_count[i] = seg_count[i] + 1
    seg = array(seg_count)
    #find where event count is larger than 1
    seg = np.where(seg>1)
    #compute the overlap length 
    ol= 0
    print(seg)
    print(type(seg))
    seg = list(seg[0])
    print(type(seg))
    print(seg)
    for j in seg:
       ol = ol + ind_list[j + 1] - ind_list[j]
    #print(ol)
    return ol

def overlap_cover(events):
    sorted_event = [list(item[1]) for item in itertools.groupby(sorted(events), key=lambda x: x[0])]
    events = []
    for chr_ in sorted_event:
        temp_list = []
        for event in chr_:
            temp_list.append(event[1:])
        events.append(temp_list)
    sum_ = 0
    for chr_ in events:
        sum_ = sum_ + overlap_cover_chr(chr_)    
    return sum_
        
scaling = 100000
x, y = plot_distr(scaling, 1000, 4)
# event = []
# start = np.random.choice(np.arange(scaling), size = 10000, p = y/sum(y))
# ave_length = 0.014
# for n in range(0, 30):
#     end = start[n] + np.random.choice([1,-1])* ave_length * scaling
#     event.append([start[n], end])
# compute_score(event)
# print(overlap_cover_chr(event)/scaling)


# L = [[0,1,2], [1,2,3], [0,2,3]]
# L = [list(item[1]) for item in itertools.groupby(sorted(L), key=lambda x: x[0])]
# L_ =[]
# for items in L:
#     items_ = []
#     for item in items:
#         items_.append(item[1:])
#     L_.append(items_)
# print(L_)
    
from dendropy.simulate import treesim

def BD_tree(birth_rate, dir, ntips):
    #global birth_rate
    print("birth rate is:", birth_rate)
    try:
        tree = treesim.birth_death_tree(birth_rate, death_rate=0.75 * birth_rate, gsa_ntax=ntips*2, num_extant_tips=ntips)
        print(tree)
        t = treesim.birth_death_tree(birth_rate, birth_rate*0.75, gsa_ntax=2*ntips, num_extant_tips=ntips)
    except TypeError:
        print("dendropy tree simulation failed, try one more time")
        return BD_tree(birth_rate, dir, ntips)


#tree = treesim.birth_death_tree(1, 0.75, gsa_ntax=100, num_extant_tips=200)
