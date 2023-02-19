#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 18:40:28 2023

@author: haoran
"""

'''This is another detection method I came up with: first transcribe the DNA base sequence into RNA base sequence, 
then translate it into a protein peptide chain and compare the amino acid sequence of the peptide chain.'''
#这是个我想出来的另一种检测方法：将DNA碱基序列先转录成RNA碱基序列，再翻译成蛋白质的肽链后比对肽链的氨基酸顺序。
'''Now, it is still in a state of not being able to fully run normally. I have tested several of the core functions. 
With a small amount of DNA sequence samples, 
the comparison can be basically completed and the correct results can be output. 
However, after reading the 'fasta' file in this program, it can't run normally.'''
#但是目前仍处于无法完全正常运行的状态，我测试过其中的几个核心函数，用少量的DNA序列样本可以基本完成比对并输出正确结果，但是在这个程序里面读取'fasta'文件之后就无法正常运行
'If someone visits this program, please feel free to enlighten me!'
#若有人莅临查看，望不吝赐教！

import Bio
from Bio import SeqIO
import re
import os


#读取基因报告文件函数
def get_report_data(file_name):
    return Bio.SeqIO.read('./gene_report/'+str(file_name),'fasta')


#读取文件夹下全部fasta文件并以碱基序列的方式输出的函数
def read_all_fasta(folder_path):
    fasta_list = []
    seq_list=[]
    #循环指定文件夹
    for file in os.listdir(folder_path):
        #fasta结尾的话
        if file.endswith('.fasta'):
            #将文件名整合在一起
            fasta_file = os.path.join(folder_path, file)
            #读取文件并把每个基因分别放入list中
            for gene in SeqIO.parse(fasta_file, 'fasta'):
                fasta_list.append(gene)
            #把碱基序列提出来
            for geneseq in fasta_list:
                seq_list.append(geneseq)
    return seq_list

disease=read_all_fasta('./disease')
for i in range(len(disease[:])):
    disease[i]=str(disease[:][i].seq)


#翻译函数
def rna_trans_proteip(seq):
    codonTable = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'',
        'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W',
    }
    proteinSeq=''
    for codonStart in range(0,len(seq),3):
        codon=seq[codonStart:codonStart+3]
        if codon in codonTable:
            proteinSeq+=codonTable[codon]
    return proteinSeq

def dna_trans_rna(seq):
    rnaSeq=re.sub('T','U',seq)
    return rnaSeq

#将DNA中的密码子提取并转化成蛋白质的函数
def dna_translate_protein(dna):
    #处理DNA部分
    dna_seqs=[]
    for i in range(len(dna)):
        dna=str(dna)
        #找到初始密码子，把初始密码子前面掐掉
        if dna[i:i+3]=='ATG':
            #切出掐头部分用于循环
            headless=dna[i:]
            for j in range(len(headless)):
                #找到中止密码子，把终止密码子后面去掉
                if headless[j:j+3] in ['TGA','TAA','TAG']:
                    product=headless[:j+3]
                    #成品需要可以被3整除，不然会编码错误
                    if len(product)%3==0:
                        #满足条件后收入list
                        dna_seqs.append(product)
    #处理RNA部分，将DNA的list跑一遍转录的for循环，生成RNA的list
    rna_seqs=[]
    for k in range(len(dna_seqs)):
        rna=re.sub('T','U',dna_seqs[k])
        rna_seqs.append(rna)
    #处理蛋白质部分，将RNA的list跑一遍翻译的for循环，生成含有肽链的list
    pro_seqs=[]
    for l in range(len(rna_seqs)):
        pro=rna_trans_proteip(rna_seqs[l])
        pro_seqs.append(pro)
    result=[]
    #去除残链部分，若一个蛋白质中包含多个甲硫氨酸(ATG->AUG->M)可能会包含多段不完整的肽链，用下面的代码消除
    #先循环基因中的每段肽链
    for m in range(len(pro_seqs)):
        flag=True
        #循环其中一个之后再次循环进行依次比对
        for n in range(len(pro_seqs)):
            #如果检测出来残肽链（一个肽链存在于另一个肽链之中）就不添加进结果
            if m!=n and pro_seqs[m] in pro_seqs[n]:
                flag=False
        if flag:
            result.append(pro_seqs[m]) 
    print(result)
    return result


#比较函数
def compare_dna_by_protein(data,sample):
    pro_data=dna_translate_protein(data)
    pro_sample=dna_translate_protein(sample)
    same_list=[]
    for i in range(len(pro_data)):
        same=0
        if len(pro_data[i])>len(pro_sample[i]):
            for j in range(len(pro_sample[i])):
                if pro_data[i][j]==pro_sample[i][j]:
                    same+=1
            similarity=same/len(pro_sample[i])
        else:
            for j in range(len(pro_data[i])):
                if pro_data[i][j]==pro_sample[i][j]:
                    same+=1
            similarity=same/len(pro_data[i])
        print(similarity)
        same_list.append(similarity)
        
    #print(same_list)
    if same_list==[]:
        return 0
    else:
        return sum(same_list)/len(same_list)


#输出结果函数
def output_result(data,threshold=0.95):
    result=[]
    for seq in disease:
        similarity=compare_dna_by_protein(data, seq)
        result.append(similarity)
    risk=0
    for i in result:
        if i<=threshold:
            risk+=1
    if risk<len(disease):
        return '依此阀值计算，受验个体大概率不携带相关遗传病基因'
    else:
        return '依此阀值计算，受验个体很可能携带相关遗传病基因'
        


    
#print(compare_dna_by_protein('ATGACCACGTAG', 'ATGACGACCTGA'))




def main():
    report_data=get_report_data('gene_seq.fasta').seq
    pro_report_data=dna_translate_protein(report_data)
    print(output_result(pro_report_data,0.95))#<=set the threshold here
    



if __name__=='__main__':
    main()
