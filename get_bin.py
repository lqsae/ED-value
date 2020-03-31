#usr/bin/python3
#author liuqingshan
'''caclulate the mean ed value in windown '''
import numpy as np
import pandas as pd
from argparse import ArgumentParser
def get_args():
    parser=ArgumentParser(description='get ed data')
    parser.add_argument('--i',help='input ed value file')
    parser.add_argument('--step',type=int,help='length of step')
    parser.add_argument('--o',help='output_file')
    parser.add_argument('--ed',type=int,help='ed value')
    args=parser.parse_args()
    return args
def get_bin_ed(infile,step,ed,outfile):
    df=pd.read_csv(infile,sep='\t')
    f1=open(outfile,'w+')
    dict_window1={}
    for chr,data in df.groupby('chr'):
        pos_max=data['pos'].max()
        all_qujian=range(0,pos_max,step)
        s={}
        for m in all_qujian :
            s[m]=0
        dict_window1[chr]=s
        for i in range(len(data)):
            qujian=int(int(data.iloc[i,1])/step)*step
            dict_window1[chr][qujian]+=1
    dict_window2={}
    for chr,data in df.groupby('chr'):
        pos_max=data['pos'].max()
        all_qujian=range(0,pos_max,step)
        s={}
        for m in all_qujian :
            s[m]=0
        dict_window2[chr]=s
        for i in range(len(data)):
            qujian=int(int(data.iloc[i,1])/step)*step
            value=data.iloc[i,8]**ed
            dict_window2[chr][qujian]+=value
    head="chr"+'\t'+"start"+"\t"+"end"+"\t"+"ED"+"\n"
    f1.write(head)
    for m in dict_window1.keys():
        for j in dict_window1[m]:
            if dict_window1[m][j]==0:
                line=str(m)+'\t'+str(j)+'\t'+str(j+step)+'\t'+str(0)+'\n'
                f1.write(line)
            else:
                line=str(m)+'\t'+str(j)+'\t'+str(j+step)+'\t'+str(dict_window2[m][j]/dict_window1[m][j])+'\n'
                f1.write(line)
def main():
    args=get_args()
    infile=args.i
    outfile=args.o
    ed=args.ed
    step=args.step
    get_bin_ed(infile,step,ed,outfile)
if __name__=='__main__':
    main()  
