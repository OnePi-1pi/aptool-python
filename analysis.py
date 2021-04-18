# -*- encoding: utf-8 -*-
def analysis(dicl1,result):
    score={}
    sumrow=[]
    result[dicl1].append("序列得分")
    for i in range(len(result[dicl1])-1):
        sumrow.append(0)
    for i in result:
        if i!=dicl1 :
            #i不在第一行则
            numrow=0
            for j in result[i]:
                sumrow[numrow]+=j
                #print("sumrow[%d]=%s"%(numrow,sumrow[numrow]))
                numrow+=1
    #上面得到了sumrow
    if numrow > 2:
        for i in result:
            if i!=dicl1 :
                score[i]=0
                j=0
                for n in result[i]:
                    #print("result[i][j-1]=%s|result[i][j]=%s|n=%d|i=%s|j=%s"%(result[i][j-1],result[i][j],n,i,j))
                    if n != 0 and j > 1 and sumrow[j]!= 0 and sumrow[j-1]!=0:
                        diff=(result[i][j]-result[i][j-1])
                        if diff > 0:
                            booln=1
                        else:
                            booln=(result[i][j]/sumrow[j])-(result[i][j-1]/sumrow[j-1])
                        if booln > 0:
                            score[i]+=int(100/(numrow-2)+diff/sumrow[0]*10000*numrow)
                    else:
                        score[i]+=0
                    j+=1
                result[i].append(score[i])
    #print(sumrow)
    return(result)


































