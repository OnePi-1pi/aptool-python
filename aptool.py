# -*- encoding: utf-8 -*-
import re
import sys
import os
import fileinput
import copy
import analysis
from time import time
################################
def main ():
    fqpaths=[]
    #fqpaths用于存储fastq文件路径
    primerfile=""
    #primerfile用于存储primer文件的路径
    pathin=input("输入路径")
    #输入路径程序开始
    start=time()
    #开始计时
    pathin=pathin.replace("\\","/")
    pathin=os.path.abspath(pathin)+"\\"
    #用于windows系统替换路径
    files = os.listdir(pathin)
    for fr in files:
        #将文件名和后缀名赋在两个值中
        name = os.path.splitext(fr)[0]
        suff = os.path.splitext(fr)[-1]
        if suff ==".fq"or suff ==".fastq":
        #如果找到符合的文件。
            fqpaths.append(pathin+fr)
            filename = name
        if name == "primer":
            primerfile=pathin+fr
    if fqpaths == []:
        print("没有找到合适的fastq文件,请检查文件路径是否正确")
        os.system("pause")
        sys.exit()
    #没有fastq文件则抛出错误。
    if primerfile == "":
        print("没有找到合适的primer.csv文件,请检查文件路径是否正确")
        os.system("pause")
        sys.exit()
    #没有primer文件则抛出错误。
    #print(fqpaths)
    #print(primerfile)
    pnum=0
    primernames=[]
    primers=[]
    oppoprimers=[]
    array0pnum=[0]
    global results
    results={}
    global resultd
    resultd={}
    stop1=time()
    print("第一阶段用时%s秒"%str(stop1-start))
    with open(primerfile, 'r', encoding='UTF-8') as primerf:
        for pr in primerf:
                    
                    if pnum==0:
                        numn=int(pr.split(',')[0])
                        rightprimer=((pr.split(',')[1]).replace("\r","").replace("\n","").replace(" ",""))
                        rightprimer=rightprimer.upper()
                        dicl1="序列/随机区长度%d" %numn
                        results[dicl1]=["总数"]
                    else:
                        pname=pr.split(',')[0]
                        primernames.append(pname)
                        strpri=(pr.split(',')[1]).replace("\r","").replace("\n","").replace(" ","").upper()
                        primers.append(strpri)
                        strpri=strpri[::-1]
                        strpri=strpri.replace("A","t").replace("C","g").replace("G","c").replace("T","a")
                        strpri=strpri.upper()
                        results[dicl1].append(pname)
                        array0pnum.append(0)
                    pnum+=1
    resultd=(copy.deepcopy(results))
    #print("primers=%s"%primers)
    numsing=numn+len(rightprimer)
    print ("随机区+末端引物个数为%s" %numsing)
    print ("末端序列为%s" %rightprimer)
    seqnums=1
    seqnumd=1
    stop2=time()
    print("第二阶段用时%s秒"%str(stop2-stop1))
    for fqpath in fqpaths:
        s=1
        with open(fqpath,'r',-1)as file_all:
            for line in file_all.readlines():
                s+=1
                #if s % 1000000 == 3:
                        #print("行数%d" %(s-3))
                if s%4==3:
                    row=0
                    for primer in primers:
                        ma=re.search(primer+"(.*)",line,re.I)
                        if ma:
                            if len(ma.group(1))>=numsing:
                                    #print("第%d个序列" %seqnums)
                                    temps=(re.match('.{%d}' %numsing,ma.group(1),re.I).group(0))
                                    #temps:挑选后的符合的序列。      
                                    #seqnums是已经有的不重复单端序列
                                    #查找序列是否已经存在
                                    if temps in results:
                                            results[temps][0]+=1
                                            results[temps][row+1]+=1

                                    ####序列不存在，新建一行，并把对应地方加1
                                    else:
                                        #sls.add(temps)
                                        seqnums+=1
                                        results[temps]=copy.deepcopy(array0pnum)
                                        results[temps][0]=1
                                        results[temps][row+1]=1

                            if re.match('(.*)%s' %rightprimer,ma.group(1),re.I):
                                    tempd=(re.match('(.*)%s' %rightprimer,ma.group(1),re.I).group(1))
                                    #tempd:挑选后的符合的双端序列。      
                                    #seqnumd是已经有的不重复双端序列
                                    #查找序列是否已经存在
                                    if tempd in resultd:
                                            resultd[tempd][0]+=1
                                            resultd[tempd][row+1]+=1

                                    ####序列不存在，新建一行，并把对应地方加1
                                    else:
                                        #sld.add(tempd)
                                        seqnumd+=1
                                        resultd[tempd]=copy.deepcopy(array0pnum)
                                        resultd[tempd][0]=1
                                        resultd[tempd][row+1]=1

                        row+=1
    stop3=time()
    print("第三阶段用时%s秒"%str(stop3-stop2))
    results[dicl1][0]="每条总数/全部%d条符合序列" %seqnums
    resultd[dicl1][0]="每条总数/全部%d条符合序列" %seqnumd
    results=analysis.analysis(dicl1,results)
    resultd=analysis.analysis(dicl1,resultd)
    ####输出序列总数
    filesingle=open(pathin+filename+"单端.csv","w+",newline='\n')
    for i in results:
        if i!=dicl1 and results[i][0] >= 0:
            #单端筛选器的个数限制在此处更改
            filesingle.write(str(i)+",")
            for j in results[i]:
                    filesingle.write(str(j)+",")
                    
            filesingle.write('\n')
        else:
            filesingle.write(str(i)+",")
            for j in results[i]:
                    filesingle.write(str(j)+",")
                    
            filesingle.write('\n')
    filesingle.close()
    #单端数据带筛选器
    filedouble=open(pathin+filename+"双端.csv","w+",newline='\n')
    for i in resultd:

        if  i!=dicl1 and resultd[i][0] >=0:
            #双端筛选器的个数限制在此处更改
            filedouble.write(str(i)+",")
            for j in resultd[i]:
                    filedouble.write(str(j)+",") 
            filedouble.write('\n')
        else:
            filedouble.write(str(i)+",")
            for j in results[i]:
                    filedouble.write(str(j)+",")
                    
            filedouble.write('\n')
    filedouble.close()
    stop4=time()
    print("总计时"+str(stop4-start) + "秒")
    #双端数据带筛选器

main()
print("程序已运行完毕")
os.system("pause")
