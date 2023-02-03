#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 03:42:46 2023

@author: haoran
"""

import Bio
from Bio import SeqIO
import os




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
#疾病数据目录
disease=read_all_fasta('./disease')


#读取基因报告文件函数
def get_report_data(file_name):
    return Bio.SeqIO.read('./gene_report/'+str(file_name),'fasta')


#碱基序列对比函数     
def compare_gene_seq(reference_seq, test_seq):
    #利用一个数据统计相同的数量
    same = 0
    #计算相同的数量
    if len(reference_seq)>len(test_seq):
        for i in range(len(test_seq)):
            if reference_seq[i]==test_seq[i]:
                same+=1
    else:
        for i in range(len(reference_seq)):
            if reference_seq[i]==test_seq[i]:
                same+=1
    #计算相似率
    similarity_rate = same / len(reference_seq)
    return similarity_rate


#处理报告数据
def process_report_data(file_data,threshold=0.8):
    #可能进行多次对比所以先用一个list收
    result=[]
    for gene_seq in disease:
        #调用碱基序列对比函数，计算出每个的相似值
        similarity=compare_gene_seq(file_data, gene_seq)
        result.append(similarity)
    risk=0
    #风险是用来统计结果的，处理的数据报告若是与样本相似率中的一条都达不到阀值，就会报危险
    for i in result:
        if i<=threshold:
            risk+=1
        else:
            pass
    if risk<len(disease):
        return '依此阀值计算，受验个体大概率不携带相关遗传病基因'
    else:
        return '依此阀值计算，受验个体很可能携带相关遗传病基因'
    
    
    





#主函数
def main():
    #获取基因报告数据
    report_data = get_report_data('gene_seq.fasta')
    
    #分析基因信息，确定是否存在遗传病风险
    disease_risk = process_report_data(report_data,0.7)#<=set the threshold here
    print(disease_risk)



if __name__ == '__main__':
    main()
    
    
    
    
    
    
    
    
    
    
