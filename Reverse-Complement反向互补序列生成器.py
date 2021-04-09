# -*- encoding: utf-8 -*-
import re
import os
import fileinput
import copy
from time import time
################################

def main ():
    fqpaths=[]
    pathin=input("输入路径")
    #start
    start=time()
    #
    pathin=pathin.replace("\\","/")
    pathin=os.path.abspath(pathin)+"\\"
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

    pnum=0
    primernames=[]
    primers=[]
    oppoprimers=[]
    array0pnum=[0]
    results={}
    resultd={}
    stop1=time()
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
                        oppoprimers.append(strpri)
                        results[dicl1].append(pname)
                        array0pnum.append(0)
                    pnum+=1
    resultd=(copy.deepcopy(results))
    print("primerf=%s"%primerf)
    print("primers=%s"%primers)
    print("oppoprimers=%s"%oppoprimers)
    with open(primerfile, 'rb') as pf:
        pf.seek(-1,2)
        qwer=pf.read(1)
        print("qwer=%s"%qwer)
    files=open(pathin+"primer.csv","a+",newline='\n')
    j=0
    if qwer != b"\n":
        files.write('\n')
    for i in primernames:
        files.write(str(i)+"_opposite,"+oppoprimers[j])     
        files.write('\n')
        j=j+1
    files.close()
main()
