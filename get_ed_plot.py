# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:49:01 2019

@author: liuqingshan
"""
from itertools import groupby, cycle
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pandas import DataFrame
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm 
from argparse import ArgumentParser
lowess=sm.nonparametric.lowess
def get_args():
    parser=ArgumentParser(description='get ed data')
    parser.add_argument('--i',help='input ed value file')
    parser.add_argument('--o',help='output_file')
    parser.add_argument('--ed',help='ed value')
    args=parser.parse_args()
    return args
def plot_ed(infile,ED,outfile):
    data=pd.read_csv(infile,sep='\t')
    data['ED2']=data['ED']**2
    data['ED4']=data['ED']**4
    data['ED6']=data['ED']**6
    data['position']=range(len(data))
    length=len(set(data['chr']))
    if length >20:
        chrlist=set(data['chr'])[0:20]
        data=data[data['chr'].isin(chrlist)]
    else:
        data=data
    plt.figure(figsize=(10,4))
    plt.style.use('ggplot')
    xs_by_id = {} # use for collecting chromosome's position on x-axis
    name=['blue','red','green','yellowgreen', 'orangered', 'deepskyblue', 'deeppink', 'seagreen','maroon']*10
    colors=iter(name)
    ed_value=[]
    pos=[]
    for seqid, rlist in data.groupby('chr', sort=True):  
        pos.extend(rlist['position'])
        ed_value.extend(rlist[ED])
        x=rlist['position']
        y=rlist[ED]
        c=next(colors)
        xs_by_id[seqid] =(min(rlist['position']) +max(rlist['position']))/ 2#中间位置
        plt.scatter(x=x, y=y,alpha=0.5,c=c,s=5,label=seqid)
    xticks=list(xs_by_id.keys())
    xpos=list(xs_by_id.values())
    plt.xticks(xpos,xticks)
    plt.xlabel('position')
    plt.ylabel(ED)
    #plt.legend()
    z=lowess(ed_value,pos,frac=0.01)#lowess函数拟合
    plt.plot(z[:,0],z[:,1],color='black',lw=1)
    plt.title('lowess')
    plt.savefig(outfile,format="pdf",dpi=800)
def main():
    args=get_args()
    infile=args.i
    outfile=args.o
    ED=args.ed
    plot_ed(infile,ED,outfile)   
if __name__=="__main__":
    main()      
