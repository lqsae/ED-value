import parsevcf
import collections
import math
import os


import pandas as pd
from scipy.stats import chi2
from scipy.special import comb
import argparse

all_info = collections.namedtuple('info', 'chr pos ref alt s1_ref_ad s1_ad_alt s2_ref_ad s2_ad_alt ed lod')


class ED(parsevcf.VCF):
    alle_frq = collections.namedtuple('alle_freq', 'ref alt')
    def __init__(self, vcf_file, type, p1=None, p2=None, s1=None, s2=None,):
        # 继承vcf方法
        super(ED, self).__init__(vcf_file)
        # 执行VCF初始化
        parsevcf.VCF.__init__(self, vcf_file)
        self.type = type
        self.p1 = p1
        self.p2 = p2
        self.s1 = s1
        self.s2 = s2

    def __iter__(self):
        for line in parsevcf.VCF.get_line(self.vcf):
            if not line.startswith('#'):
                tags = line.split()
                alt = tags[4]
                if len(alt.split(',')) == 1:
                    record = parsevcf.Record(tags, self.map_index, self.samples)
                    miss = record.miss()
                    if miss == 0:
                        all_ad_s1 = record.get_sample_ad(self.s1)
                        all_ad_s2 = record.get_sample_ad(self.s2)
                         #两个亲本选择和纯和差异的位点
                        if self.is_dobule_p():
                            p1_alle = record.get_sample_gt(self.p1)
                            p2_alle = record.get_sample_gt(self.p2)
                            if self.is_hom_diff(p1_alle, p2_alle):
                                ed = self.get_ed(all_ad_s1, all_ad_s2)
                                lod = self.get_lod_value(all_ad_s1, all_ad_s1)
                                info = [record.chrom(), record.pos(), record.ref(), record.alt(), all_ad_s1[0],
                                        all_ad_s1[1], all_ad_s2[0], all_ad_s2[1], ed, lod]
                                yield all_info._make(info)
                        #一个亲本选择纯和位点
                        elif self.is_only_one_p():
                            p1_alle = record.get_sample_gt(self.p1)
                            if self.is_hom(p1_alle):
                                ed = self.get_ed(all_ad_s1, all_ad_s2)
                                lod = self.get_lod_value(all_ad_s1, all_ad_s1)
                                info = [record.chrom(), record.pos(), record.ref(), record.alt(), all_ad_s1[0],
                                        all_ad_s1[1], all_ad_s2[0], all_ad_s2[1], ed, lod]
                                yield all_info._make(info)
                        #无亲本不过滤
                        elif self.is_only_s():
                            ed = self.get_ed(all_ad_s1, all_ad_s2)
                            lod = self.get_lod_value(all_ad_s1, all_ad_s1)
                            info = [record.chrom(), record.pos(), record.ref(), record.alt(),
                                    all_ad_s1[0], all_ad_s1[1], all_ad_s2[0], all_ad_s2[1], ed, lod]
                            yield all_info._make(info)

    def head(self):
        """
        返回表头
        :return: str
        """
        head = 'chr\tpos\tref\talt\tdepth_ref_{0}\tdepth_ref_{1}\tdepth_alt_{0}' \
               '\tdepth_alt_{1}\tED\tLOD\n'.format(self.s1, self.s2)
        return head

    def is_dobule_p(self):
        """
        判断是否为双亲本
        :return: bool
        """
        return self.type == 'twop'

    def is_only_one_p(self):
        """
        判断是否为单亲本
        :return: bool
        """
        return self.type == 'onep'

    def is_only_s(self):
        """
        判断是否仅为子代
        :return: bool
        """
        return self.type == 'onlys'

    def is_hom(self, p_alle):
        """
        判断是否为纯和基因型
        :param p_alle: str 0/1, 1/1, 0/0
        :return: bool
        """
        p = set(p_alle.split('/'))
        return len(p) == 1

    def is_hom_diff(self, p1_alle, p2_alle):
        """
        判断两个亲本是否纯和差异
        :param p1_alle: str 0/1, 1/1, 0/0
        :param p2_alle: str 0/1, 1/1, 0/0
        :return: bool
        """
        p1 = set(p1_alle.split('/'))
        p2 = set(p2_alle.split('/'))
        return (len(p1) == 1) and (len(p2) == 1) and (p1_alle != p2_alle)

    def get_freq(self, ad_s1):
        ref_dp, alt_dp = ad_s1
        # alt_dp = ad_s1.ref
        # alt_dp = ad_s1.alt
        all_dp = ref_dp + alt_dp
        return ED.alle_frq._make([ref_dp / all_dp, alt_dp / all_dp])

    def get_ed(self, ad_s1, ad_s2):
        freq_s1 = self.get_freq(ad_s1)
        freq_s2 = self.get_freq(ad_s2)
        sum_pow = math.pow((freq_s1.ref - freq_s2.ref), 2) + math.pow((freq_s1.alt - freq_s2.alt), 2)
        return math.sqrt(sum_pow)

    def get_prob(self, nL, nAL, naL, pL):
        return comb(nL, nAL)*math.pow(pL, nAL) * math.pow(1-pL, naL)

    def get_lod_value(self, L, H):
        #for low
        nAL, naL = L
        # naL = L.alt
        nL = nAL + naL
        pL = nAL/nL
        #for high
        nAH, naH = H
        nH = nAH + naH
        pH = nAH/nH
        unexist_QTN = comb(nL, nAL)*comb(nH, nAH)*math.pow(0.5, (nL+nH))
        exist_QTN = self.get_prob(nL, nAL, naL, pL) * self.get_prob(nH, nAH, naH, pH)
        #exist_QTN = comb(nL, nAL)*math.pow(pL, nAL) * math.pow(1-pL, naL) *comb(nH, nAH)*math.pow(pH, nAH) * math.pow(1-pH,naH)
        lod_value = exist_QTN/unexist_QTN
        return math.log(lod_value, 10)


def smooth(infile, w, oufile, pos='Pos', value='LOD'):
    df = pd.read_csv(infile, sep='\t')
    x = df[pos]
    y = df[value]
    # xfirst = x[0]
    olen = len(x)
    # beginning intervals
    xfirst = x[0]  # start site
    startband = (x <= xfirst + w)
    xstart = xfirst - (x[startband] - xfirst)
    ystart = y[startband]

    # end of intervals
    xlast = x.iloc[-1]  # last site
    endband = x >= xlast - w
    xend = xlast + (xlast - x[endband])
    yend = y[endband]

    a = list(xstart[::-1])  # reversal
    a.pop(-1)  # remove the last site of a

    b = list(xend[::-1])
    b.pop(0)  # remove the first site of b

    c = list(ystart[::-1])
    c.pop(-1)  # remove the last site of c

    d = list(yend[::-1])
    d.pop(0)  # remove the first site of d

    X = []
    X.extend(a)
    X.extend(x)
    X.extend(b)

    Y = []
    Y.extend(c)
    Y.extend(y)
    Y.extend(d)

    df_X = pd.Series(X)
    df_Y = pd.Series(Y)
    ywins = []
    for i in range(len(df_X)):
        c = df_X[i]
        inband = ((c - w) <= df_X) & (df_X >= (c - w))

        xfrac = abs((df_X[inband]) - c) / w

        xwt = pow((1.0 - pow(abs(xfrac), 3)), 3)  # Weights
        xwt[abs(xfrac) >= 1] = 0
        ywin = sum(df_Y[inband] * xwt) / sum(xwt)
        ywins.append(ywin)
    h1 = len(a) + 1
    h2 = len(a) + olen
    df['smooth{}'.format(value)] = ywins[h1:h2+1]
    df['p'] = df['smooth{}'.format(value)].apply(lambda x: chi2.sf(x*math.log(10)*2, 1))
    df.to_csv(oufile, sep='\t', index=None)


def slide(step, windsize, length):
    wind = collections.namedtuple('wid', 'start end')
    start = 0
    i = 0
    end = 0
    while end <= length:
        start = step*i + 1
        end = start + windsize
        i += 1
        yield wind._make([start, end])
    if end > length:
        end = length
        yield wind._make([start, end])


def win_stat(info_data, number):
    info = collections.namedtuple('win', 'chr start end pos ed')
    flag = 0
    sum = 0
    yield 'chrt\tstart\tend\tPos\ted'.split('\t')
    for i in info_data:
        if flag == 0:
            start = i.pos
            flag += 1
            sum += i.ed
        elif flag < number:
            flag += 1
            sum += i.ed
        elif flag == number:
            end = i.pos
            pos = int((start + end)/2)
            chr_name = i.chr
            sum += i.ed
            mean = sum /flag
            yield info._make([chr_name, start, end, pos, mean])
            sum = 0
            flag = 0


def winsize(step, winsie, info):
    dict_info = collections.defaultdict(list)
    yield 'chrt\tstart\tend\tPos\ted'.split('\t')
    key_info = collections.namedtuple('win', 'chr start end pos ed')
    for i in info:
        chr_name = i.chr
        pos = int(i.pos)
        ed = i.ed
        win = pos//step
        win_start = win*step + 1
        win_end = win_start + winsie
        key = (chr_name, win_start, win_end)
        dict_info[key].append(ed)
    for key, value in dict_info.items():
        mean = sum(value)/len(value)
        chr_name, win_start, win_end = key
        pos = int((win_start + win_end)/2)
        yield key_info._make([chr_name, win_start, win_end, pos, mean])


def get_args():
    parser = argparse.ArgumentParser(description="calculate the ed value")
    parser.add_argument("--v", "--vcf", help="vcf file")
    parser.add_argument("--p1", help="the name of p1")
    parser.add_argument("--p2", type=str, help="the name of p2")
    parser.add_argument("--s1", type=str, help="the name of s1")
    parser.add_argument("--s2", type=str, help="the name of s2")
    parser.add_argument("--type", type=str, help="the type of the data [twop ,onep, onlys]")
    parser.add_argument("--w", type=int, help="g static of the window", default=100000)
    parser.add_argument("--o", type=str, help="out file")
    parser.add_argument("--ws", type=int, help="the length of the ed window")
    parser.add_argument("--step", type=int, help="the length of the ed step ")
    return parser.parse_args()


def main():
    args = get_args()
    vcf_file = args.v
    p1 = args.p1
    p2 = args.p2
    s1 = args.s1
    s2 = args.s2
    ty = args.type
    win_size_ed = args.ws
    step = args.step
    out = '{}_{}.ed.value.data'.format(s1, s2)

    if not os.path.exists(out):
        out_file = open(out, 'w')
        eds = ED(vcf_file, ty, p1, p2, s1, s2)
        #gethead
        out_file.write(eds.head())
        for ed in eds:
            line = '\t'.join(map(str, ed)) + '\n'
            out_file.write(line)
        #write ed value
        out_file.close()

    #get win ed
    win_ed_name = '{0}{1}.win{2}.ed.data'.format(s1, s2, win_size_ed)
    if not os.path.exists(win_ed_name):
        eds = ED(vcf_file, ty, p1, p2, s1, s2)
        win_ed_file = open(win_ed_name, 'w')
        for info in winsize(step, win_size_ed, eds):
            line = '\t'.join(map(str, info)) + '\n'
            win_ed_file.write(line)
    #write win ed value
        win_ed_file.close()

    #get smooth ed ad write
    if os.path.exists(out):
        smooth(infile=out, w=args.w, oufile=args.o, pos='pos', value='LOD')


main()
