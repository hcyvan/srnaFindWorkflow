#!/usr/bin/python3
import sys
import os
import argparse

import numpy as np


def get_candidates(bedgraph, cutoff=10, p=False):
    result=[]
    with open(bedgraph) as f:
        stack=[]
        for l in f.readlines():
            pos = l.strip().split('\t')
            if float(pos[2]) >= cutoff:
                stack.append(pos)
            else:
                if len(stack)>0:
                    if sum([float(x[2]) for x in stack])/len(stack) >= cutoff*2 and len(stack)>=40 and len(stack)<=500:
                        r=[stack[0][0], str(int(stack[0][1])+1), str(int(stack[-1][1])+1)]
                        if p:
                            print('\t'.join(r))
                        else:
                            result.append(r)
                stack=[]
    if not p:
        return result

def get_candidates_from_file(candidate_file):
    with open(candidate_file) as f:
        candidates = [x.strip().split('\t') for x in f.readlines()]
    return [[x[0],int(x[1]),int(x[2])] for x in candidates]

def find_srna(input1,input2,gff,cutoff=10):
    candidates1=get_candidates(input1, cutoff=cutoff)
    print('Finding {} candidates in + strand...'.format(len(candidates1)))
    candidates2=get_candidates(input2, cutoff=cutoff)
    print('Finding {} candidates in - strand...'.format(len(candidates2)))

    g1=[]
    g2=[]
    with open(gff) as f:
        for l in f.readlines():
            if l.startswith('#'):
                continue
            items=l.strip().split('\t')
            if items[6]=='+':
                g1.append(items)
            elif items[6]=='-':
                g2.append(items)

    def _intersect(g, candidates):
        result=[]
        for c in candidates:
            for i in range(len(g)):
                if i==1:
                    continue
                left=int(g[i-1][4])+100 if g[i-1][2] in('gene','CDS') else int(g[i-1][4])
                right=int(g[i][3])-60 if g[i-1][2] in('gene','CDS') else int(g[i][3])

                if int(c[1]) >= left and int(c[2]) <= right:
                    result.append(c)
        return result

    r1 = _intersect(g1, candidates1)
    print('Get {} sRNA in + strand...'.format(len(r1)))
    r2 = _intersect(g2, candidates2)
    print('Get {} sRNA in - strand...'.format(len(r2)))

    return r1, r2

def find(args):
    def _save(r, output_file):
        rna = []
        for i in r[0]:
            i.append('+')
            rna.append(i)
        for i in r[1]:
            i.append('-')
            rna.append(i)
        rna = sorted(rna,key=lambda x:x[1])
        with open(output_file,'w') as f:
            f.writelines(['\t'.join(x)+'\n' for x in rna])

    if args.input_and_output_files is None:
        r=find_srna(args.input1, args.input2, args.gff, int(args.cutoff))
        _save(r, args.output)
    else:
        with open(args.input_and_output_files) as f:
            for l in f.readlines():
                if l.startswith('#') or l.strip()=='':
                    continue
                k = l.strip()
                input_1 = '{}.1.scale.bedgraph'.format(k)
                input_2 = '{}.2.scale.bedgraph'.format(k)
                output = '{}.o.scale.bedgraph'.format(k)

                print('parsing input1: {}, input2: {}'.format(input_1, input_2))
                r=find_srna(os.path.join(args.input_file_dir, input_1),
                            os.path.join(args.input_file_dir, input_2),
                            args.gff,
                            int(args.cutoff))

                _save(r, os.path.join(args.output_file_dir, output))

def scale_and_merge(files_and_score, input_dir):
    m=None
    for fs in files_and_score:
        a_path=os.path.join(input_dir, fs[0])
        print(' reading {}'.format(a_path))
        a=np.genfromtxt(a_path)[:,2]*1000000/int(fs[1])
        if m is not None:
            m=np.vstack((m,a))
        else:
            m=a
    if m.ndim==1:
        return m.tolist()
    else:
        return np.mean(m,axis=0).tolist()

def parse_input_file(file_path):
    size_map=dict()
    with open(file_path) as f:
        for l in f.readlines():
            if l.startswith('#'):
                continue
            item = l.strip().split('\t')
            if size_map.get(item[2]):
                size_map[item[2]].append(item[0:2])
            else:
                size_map[item[2]]=[item[0:2]]
    return size_map


def scale(args):
    size_map=parse_input_file(args.bedgraph_files)
    for k,v in size_map.items():
        print('scale group: {}'.format(k))
        scores=scale_and_merge(v,args.input_file_dir)
        print('printing {}.bedgraph'.format(k))
        with open(os.path.join(args.input_file_dir, v[0][0])) as fi:
                with open(os.path.join(args.output_file_dir, '{}.bedgraph'.format(k)), 'a') as fo:
                    i=0
                    for li in fi.readlines():
                        lis=li.strip().split('\t')
                        lis[2]=str(scores[i])
                        fo.write('\t'.join(lis)+'\n')
                        i=i+1

def merge_two(s1, s2):
    def split_s(s):
        return [x for x in s if x[3]=='+'],[x for x in s if x[3]=='-']  
    def merge_s(s1_x, s2_x):
        merge=[]
        merge.extend(s1_x)
        merge.extend(s2_x)
        merge=sorted(merge,key=lambda x:x[1])
        ret=[]
        tmp=None
        while len(merge)!=0:
            if tmp is None:
                tmp=merge.pop(0)
            else:
                if tmp[2] <= merge[0][1]:
                    ret.append(tmp)
                    tmp=merge.pop(0)
                elif tmp[2] >= merge[0][2]:
                    tmp[1]=merge[0][1]
                    tmp[2]=merge[0][2]
                    ret.append(tmp)
                    merge.pop(0)
                    if len(merge) != 0:
                        tmp=merge.pop(0)
                    else:
                        tmp=None
                else:
                    tmp[1]=merge[0][1]
                    ret.append(tmp)
                    merge.pop(0)
                    if len(merge) != 0:
                        tmp=merge.pop(0)
                    else:
                        tmp=None

        if tmp is not None:
            ret.append(tmp)
        return ret

    s1_1, s1_2=split_s(s1)
    s2_1, s2_2=split_s(s2)

    ret_1 = merge_s(s1_1, s2_1)
    ret_2 = merge_s(s1_2, s2_2)
    ret=[]
    ret.extend(ret_1)
    ret.extend(ret_2)

    return sorted(ret, key=lambda x: x[1])


def merge(args):
    """Merge bedgraph files
    """
    def read_file(file_name,label=''):
        ret=[]
        with open(file_name) as f:
            ret=[x.strip().split() for x in list(f.readlines())]
        ret=[[x[0],int(x[1]),int(x[2]),x[3]] for x in ret]
        ret=sorted(ret, key=lambda x:x[1])
        final=[]
        i=1
        for r in ret:
            r.append('{}.{}'.format(label,i))
            final.append(r)
            i=i+1
        return final
    merge=[]
    i=0
    with open(args.config_file) as f:
        files=['{}.o.scale.bedgraph'.format(x.strip()) for x in f]
    for fn in files:
        if args.label:
            label=args.label.split(',')[i]
        else:
            label=i
        merge=merge_two(merge, read_file(os.path.join(args.input_file_dir, fn),label=i))
        i=i+1
    merge=sorted(merge,key=lambda x:x[1])
    for m in merge:
        if args.format=='gff':
            print('{}\tSRNA\tncRNA\t{}\t{}\t0\t{}\t0\tID=srna.{};Name=srna.{}'.format(m[0],m[1],m[2],m[3],m[4],m[4]))
        else:
            print('{}\t{}\t{}\t{}'.format(m[0],m[1],m[2],m[3]))

def merge2anno(args):
    def read_gff(gff):
        ret=[]
        with open(gff) as f:
            for line in f.readlines():
                line=line.strip()
                if not line.startswith('#'):
                    i=line.split('\t')
                    ret.append([i[0],i[1],i[2],int(i[3]),int(i[4]),i[5],i[6],i[7],i[8]])
        return ret

    def parse(item):
        attrs=item[8].split(';')
        ret=dict()
        for attr in attrs:
            kv=attr.split('=')
            ret[kv[0]]=kv[1]
        return ret
    def get_name_or_gene(item):
        attrs = parse(item)
        if attrs.get('Name'):
            return attrs.get('Name')
        return attrs['gene']

    def get_flag(rna, f_1, f_2):
        """
        f_1 strand and rna strand is the same
        """
        flag=None
        for i2 in f_2:
            if rna[4] > i2[3] and rna[3] < i2[4]:
                flag='anti.{}'.format(get_name_or_gene(i2))
                break
        if flag is None:
            if f_1[0][3] > rna[4]:
                flag='inter.{}.{}'.format('start',parse(f_1[1])['Name'])
            elif f_1[-1][4] < rna[3]: 
                flag='inter.{}.{}'.format(parse(f_1[1])['Name'],'end')
            else:
                for k in range(len(f_1)-1):
                    if f_1[k][4]<rna[3] and f_1[k+1][3]>rna[4]:
                        flag='inter.{}.{}'.format(parse(f_1[k])['Name'],parse(f_1[k+1])['Name'])
                        break

        return flag

    f1=read_gff(args.file)
    f2=read_gff(args.anno)
    f2_1=[]
    f2_2=[]
    for i in f2:
        if i[2]=='gene':
            if i[6]=='+':
                f2_1.append(i)
            else:
                f2_2.append(i)
    for i,rna in enumerate(f1):
        if rna[6]=='+':
            flag=get_flag(rna,f2_1,f2_2)
        else:
            flag=get_flag(rna,f2_2,f2_1)
        f1[i][8]=rna[8]+';ncType='+flag
    f2.extend(f1)
    f2=sorted(f2,key=lambda x: x[3])
    for i in f2:
        print('\t'.join([str(x) for x in i]))

if __name__=='__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-i','--input', metavar='input file dir',default='./',
                        dest='input_file_dir', action='store', help='input file dir')
    parent_parser.add_argument('-o','--output', metavar='output file dir',default='./',
                        dest='output_file_dir', action='store', help='output file dir')
    parent_parser.add_argument('-co','--config', metavar='configure file',default='./',
                        dest='config_file', action='store', help='configure file')

    parser = argparse.ArgumentParser(description='Find Small non-coding RNAs')
    subparsers = parser.add_subparsers(help='Find sRNAs')

    # find
    parser_0 = subparsers.add_parser('find', help='a help', parents=[parent_parser])
    parser_0.add_argument('-1', metavar='input file 1',
                        dest='input1', action='store', help='The + strand bedgraph file')
    parser_0.add_argument('-2', metavar='input file 2',
                        dest='input2', action='store', help='The - strand bedgraph file')
    parser_0.add_argument('-s', metavar='sRNA bed file',
                        dest='output', action='store', help='output file')
    parser_0.add_argument('-g', '--genome', metavar='genome gff file', required=True,
                        dest='gff', action='store', help='The genome GFF file')
    #parser_0.add_argument('-c', '--candidate', metavar='bed file',
    #                   dest='candidate', action='store', help='The candidate bed file')
    parser_0.add_argument('-c', '--cutoff', metavar='cutoff', default=20,
                       dest='cutoff', action='store', help='The cutoff of bedgraph curve for sRNA')
    parser_0.add_argument('-b', metavar='input and output batch files',
                        dest='input_and_output_files', action='store', help='<input1>\t<input2>\t<output>')
    parser_0.set_defaults(func=find)

    # merge
    parser_merge = subparsers.add_parser('merge', help='merge bedgraph files', parents=[parent_parser])
    parser_merge.add_argument('-f', metavar='input files',
                        dest='files', action='append', help='Bedgraph files to be merged')
    parser_merge.add_argument('-r', '--format', metavar='format output', default='gff',
                        dest='format', action='store', help='format output')
    parser_merge.add_argument('-l', '--label', metavar='output labels',
                        dest='label', action='store', help='output labels')

    parser_merge.set_defaults(func=merge)

    #merge2anno
    parser_merge2anno = subparsers.add_parser('merge2anno', help='merge gff to annotation file', parents=[parent_parser])
    parser_merge2anno.add_argument('-f', metavar='input file',required=True,
                        dest='file', action='store', help='input_file')
    parser_merge2anno.add_argument('-a', metavar='annotation file',required=True,
                        dest='anno', action='store', help='annotation file')
    parser_merge2anno.set_defaults(func=merge2anno)

    # scale
    parser_a = subparsers.add_parser('scale', help='Scale bedgraph files', parents=[parent_parser])
    parser_a.add_argument('-f', metavar='bedgraph file list', required=True,
                        dest='bedgraph_files', action='store', help='<file>\t<library size>')
    parser_a.set_defaults(func=scale)
    if len(sys.argv)==1:
        args_list=['-h']
    else:
        args_list=sys.argv[1:]

    args = parser.parse_args(args_list)
    args.func(args)

