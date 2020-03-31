import numpy as np
from argparse import ArgumentParser

def parse_namelist():#解析命令函数参数namelist
    args=get_args()
    name_list=args.name_list
    name_list=name_list.split(',')
    return name_list


def parse_vcf(infile):
    with open(infile) as f:
        for line in f:
            if  not line.strip().startswith('##'):
                yield line

def print_head(type):#打印表头
    name_list=parse_namelist()
    if type=='yes':
        s1,s2=name_list[2:]
    elif type=='no':
        s1,s2=name_list
    head='\t'.join(['chr','pos','ref','alt','depth_ref_'+s1,'depth_ref_'+s2,'depth_alt_'+s1,'depth_alt_'+s2,'ED'])
    print(head)



def get_dict_index(infile):#获取亲本和子代的位置索引
    dict_index={}
    name_list=parse_namelist()
    f=parse_vcf(infile)
    for line in f:
        if line.strip().startswith('#'):
            head=head=line.strip().split('\t')
            break
    for name in name_list:
        dict_index[name]=head.index(name)
    return dict_index


def get_tag(line):##过滤双等位基因
    length=len(line.strip().split('\t')[4])
    if length ==1:
        tag=line.strip().split('\t')
        return tag



def get_GeneType_GeneDepth(dict_index,name_list,tag):#获取不同个体基因型和深度信息
    dict_depth={}
    for i,j in dict_index.items():
        try:
            dict_depth['genetype'+i]=tag[j].split(':')[0]
            if dict_depth['genetype'+i]=='./.':
                dict_depth['depth_ref'+i]=0
                dict_depth['depth_alt'+i]=0
            elif dict_depth['genetype'+i]=='0/0':
                dict_depth['depth_ref'+i]=int(tag[j].split(':')[1].split(',')[0])
                dict_depth['depth_alt'+i]=0
            elif dict_depth['genetype'+i]=='1/1':
                dict_depth['depth_ref'+i]=0
                dict_depth['depth_alt'+i]=int(tag[j].split(':')[1].split(',')[1])
            elif dict_depth['genetype'+i]=='0/1':
                dict_depth['depth_ref'+i]=int(tag[j].split(':')[1].split(',')[0])
                dict_depth['depth_alt'+i]=int(tag[j].split(':')[1].split(',')[1])
        except:
            pass

    return dict_depth

def filter_depth(name_list,dict_depth):#过滤深度
    ref_all_depth=[dict_depth['depth_ref'+i] for i in name_list]
    alt_all_depth=[dict_depth['depth_alt'+i] for i in name_list]
    all_depth=np.array(ref_all_depth)+np.array(alt_all_depth)
    min_depth=min(all_depth)
    flag=min_depth>=7
    return flag


def filter_genetype(type,dict_depth,name_list):#有亲本就过纯和差异位点，无亲本就不过滤
    p1=name_list[0]
    p2=name_list[1]
    alle_p1 = dict_depth.get('genetype'+p1)
    alle_p2 = dict_depth.get('genetype'+p2)
    if type=='yes':
        m1=(alle_p1=="0/0")&(alle_p2=="1/1")
        m2=(alle_p1=="1/1")&(alle_p2=="0/0")
        flag=any([m1,m2])
    elif type=='no':
        flag=True
    return flag

def get_ed(name_list,type,dict_depth):#计算ed值
    if type=='yes':
        s1,s2=name_list[2:]
    elif type=='no':
        s1,s2=name_list
    s1_ref_freq=dict_depth['depth_ref'+s1]/(dict_depth['depth_ref'+s1]+dict_depth['depth_alt'+s1])
    s1_alt_freq=dict_depth['depth_alt'+s1]/(dict_depth['depth_ref'+s1]+dict_depth['depth_alt'+s1])
    s2_ref_freq=dict_depth['depth_ref'+s2]/(dict_depth['depth_ref'+s2]+dict_depth['depth_alt'+s2])
    s2_alt_freq=dict_depth['depth_alt'+s2]/(dict_depth['depth_ref'+s2]+dict_depth['depth_alt'+s2])
    ED=np.sqrt(np.square(s1_ref_freq-s2_ref_freq)+np.square(s1_alt_freq-s2_alt_freq))
    return ED


def get_depth_S(name_list,type,dict_depth):#获取子代的深度信息最后打印出来
    if type=='yes':
        s1=name_list[2]
        s2=name_list[3]
    elif type=='no':
        s1=name_list[0]
        s2=name_list[1]
    m=[ dict_depth.get('depth_ref'+s1),dict_depth.get('depth_ref'+s2),dict_depth.get('depth_alt'+s1),dict_depth.get('depth_alt'+s2)]
    m=[str(i) for i in m]
    return m

def get_line(tag,m,ed):#格式化输出
    tag_list=[str(tag[i]) for i in [0,1,3,4]]+m
    line='\t'.join(tag_list)+'\t'+str(ed)
    return line


def format(infile,type):#主函数
    dict_index=get_dict_index(infile)
    name_list=parse_namelist()
    f=parse_vcf(infile)
    for line in f:
        if not line.strip().startswith('#'):
            tag=get_tag(line)#过滤到双等位基因
            dict_depth=get_GeneType_GeneDepth(dict_index,name_list,tag)
            flag1=filter_genetype(type,dict_depth,name_list)
            if flag1:#过滤纯合差异位点
                flag2=filter_genetype(type,dict_depth,name_list)
                if flag2:#过滤深度
                    ed=get_ed(name_list,type,dict_depth)
                    m=get_depth_S(name_list,type,dict_depth)
                    format_data=get_line(tag,m,ed)
                    print(format_data) 


def get_args():
    parser=ArgumentParser(description='get ed data')
    parser.add_argument('--i',help='input vcf file')
    parser.add_argument('--type',help='type of parent[yes or no]')
    parser.add_argument('--name_list',type=str,help='the name of parent and generation[p1,p2,s1,s2]')
    args=parser.parse_args()
    return args


def main():
    args=get_args()
    type=args.type
    infile=args.i
    print_head(type)
    format(infile,type)

if __name__=="__main__":
    main()
