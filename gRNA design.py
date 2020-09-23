import csv
import random
#This is the function that sorts a given list of crRNAs based on different criteria described in the paper
#These criteria include number of purines next to PAM, G/C content, position relative to ATG
def sorting(crRNAlist):
    score=[]
    answer=[]
    z=1
    for i in range(0,len(crRNAlist)-1):
        z=int(crRNAlist[i][1])       
        crRNA=crRNAlist[i][0]        
        length=len(crRNAlist)
        purine=[]
        if crRNA[21:23]=='GG':
            if crRNA[0:2]=='CC':
                isitplus=ifplus(crRNA)
            else:
                isitplus=True
        else:
            isitplus=False      
        if isitplus:
            c=crRNA[16:20]
            GC=gccounter(crRNA[0:20])
            AGS=0.25*c.count('G')+0.2*c.count('A')+0.15*c.count('C')
        else:
            c=crRNA[3:7]
            GC=gccounter(crRNA[3:23])
            AGS=0.25*c.count('C')+0.2*c.count('T')+0.15*c.count('G')
        if i>=length*0.8:
            closetobegining=0
        elif i>= length*0.6:
            closetobegining=0.85
        else:
            closetobegining=1
        if GC>15:
            GCS=0
        elif GC>6:
            GCS=1
        else:
            GCS=0
        
        s=[100*(0.33333*GCS+0.33333*AGS+0.33333*closetobegining)/(z**2), crRNA]        
        score.append(s)
        
    score.sort(key=lambda tup: tup[0])
    for i in score:
        answer.append(i)
    answer.reverse()
    return answer

#Calculates the minus strand when given the plus strand
def minusstrand(seq):
    answer=''
    for i in seq:
        if i =='A':
            answer=answer+'T'
        elif i =='T':
            answer=answer+'A'
        elif i =='C':
            answer=answer+'G'
        elif i =='G':
            answer=answer+'C'
    return answer
#GC Counter

def gccounter(seq):
    GC=0
    for i in range(0,len(seq)-1):
        if seq[i]=='C':
            GC+=1            
        if seq[i]=='G':
            GC+=1            
    return GC
#GC Content Calculator
def gccontent(seq):
    GC=0
    for i in range(0,len(seq)-1):
        if seq[i]=='C':
            GC+=1            
        if seq[i]=='G':
            GC+=1            
    return GC*100/len(seq)
#Reading the genome of SC
file=open('Refer.txt','r')
gene= file.read()
file.close
chromosome=[]
loc=gene.find('chromosome')

for i in range(1,17):    
    loc2=gene.find('chromosome',loc+100)
    for i in range(0,20):
        if gene[loc+i]==']':
            chrom=gene[loc+i+2:loc2]
            chromosome.append(''.join(chrom.split('\n')))           
            
    loc=loc2
#Reading the non-essential genes
with open('All.txt','r') as nsg:
    reader = csv.reader(nsg, delimiter="\t")
    non = list(reader)

def ifplus(seq):
    for i in allorfs:
        if seq==i[3]:
            if i[2]=='+':
                answer=True
            if i[2]=='-':
                answer=False
    return answer
    
    
#Reverse complement function
def revcom(seq):
    answer=''
    for i in range(0,len(seq)):       
        if seq[i]== 'A':
            answer ='T'+answer
        elif seq[i]=='C':
            answer ='G'+answer
        elif seq[i]=='G':
           answer ='C'+answer
        elif seq[i]=='T':
            answer ='A'+answer
    return answer
#gets the sequence of a crRNA and the name of the gene corresponding to it and
#returns 100 bp upstream and downstream of the crRNA in the genome but in forward direction
#Returns empty string if an essential gene
def HA(seq,seqname,plus):
    a=0
    answer=''
    b=0
    t=0
    chromenum=-1
    for i in non:        
        if seqname=='#'+i[0].split(' ')[0]:
            chromelet=i[1]
            if chromelet=='I':
                chromenum=0
            elif chromelet=='II':
                chromenum=1
            elif chromelet=='III':
                chromenum=2
            elif chromelet=='IV':
                chromenum=3
            elif chromelet=='V':
                chromenum=4
            elif chromelet=='VI':
                chromenum=5
            elif chromelet=='VII':
                chromenum=6
            elif chromelet=='VIII':
                chromenum=7
            elif chromelet=='IX':
                chromenum=8
            elif chromelet=='X':
                chromenum=9
            elif chromelet=='XI':
                chromenum=10
            elif chromelet=='XII':
                chromenum=11
            elif chromelet=='XIII':
                chromenum=12
            elif chromelet=='XIV':
                chromenum=13
            elif chromelet=='XV':
                chromenum=14
            elif chromelet=='XVI':
                chromenum=15               
            if plus:                
                a=chromosome[chromenum].find(i[4])
            else:
                a=chromosome[chromenum].find(revcom(i[4]))  
            break
    if a<0:
        for j in range(0,15):            
            if plus:                
                a=chromosome[j].find(item2[1])
                if a>-1:
                    chromenum=j                
            else:
                a=chromosome[j].find(item2[1])
                if a>-1:
                    chromenum=j        
    if a>-1:
        if plus:
            b=chromosome[chromenum][a-2000:a+2000].find(seq)
            answer=chromosome[chromenum][a-2000+b-100:a-2000+b+100]
        else:            
            b=chromosome[chromenum][a-2000:a+2000].find(revcom(seq))
            answer =revcom(chromosome[chromenum][a-2000+b-100:a-2000+b+100])
    
        
        
    else:
        answer=''
    
    
    return answer

def nonrep(length, number):
    b=[]
    k=0
    while len(b)<number:
        seq=''
        a=0
        k+=1
        ifpresent=False
        for i in range(1,length+1):
            seq=seq+random.choice('ACGT')
        
        if 'TTTTT' in seq or 'GGTCTC' in seq or 'GAGACC' in seq or 'GGGGGG' in seq or 'AAAAAA' in seq or 'CC' in seq or 'GG' in seq:
            ifpresent=True                
        else:
            if gccontent(seq)>20:
                if gccontent(seq)<80:
                   b.append(seq)
    return b
initialseq='TATCTACACG'+'GGTCTC'+'ACCAAAAC'
finalseq='GTTTTAGAGA'+'GAGACC'+'AGCGTAACTC'                       

def controlgen(number):
    answer=[]
    nexttopam=[]
    arms=nonrep(50,10*number)
    crRNA=nonrep(16,15*number)
    con=0
    for j in range(0,2*number):
        bp4=''
        for i in range(1,5):
            bp4=bp4+random.choice('ACG')
        nexttopam.append(bp4)
    while len(answer)<number:
        con=con+1
        const=arms[con*2]+arms[2*con+1]+crRNA[con]+nexttopam[con]
        if 'TTTTT' in const or 'GGTCTC' in const or 'GAGACC' in const or 'GGGGGG' in const or 'AAAAAA' in const:
            length=2
        else:            
            answer.append(initialseq+const+finalseq)
    return answer
###parsing the file into a list with each gene name like #PGK followed by its corresponding crRNAs
d=[]

with open('orf4 copy.txt','r') as f:
    reader = csv.reader(f, delimiter="\t")
    temp = list(reader)
f.close
d=d+temp

f.close
allorfs=d

counter=0
summation=0
falsepositive=[]
b=[]

for i in range (0, len(d)-1):
    first_string= d[i][0]
    if first_string=='# sequence_name:':
        b.append([])
print len(b)


j=-1
for i in range(0,len(d)-1):
    first_string= d[i][0]
    if first_string=='# sequence_name:':               
        j=j+1
        n=[]
        b[j].append('#'+ d[i][1])
    if first_string[0]!="#" and int(d[i][7])>0:
        b[j].append([d[i][3],d[i][8]])               

            
#The list with no BsaI and TTTTT >0 20mer hit and <2 12mer hit is made at this point
b2empty=0
bempty=0
for i in range(0,len(b)-1):    
    if len(b[i])<2:
        bempty+=1
    for j in range(len(b[i])-1,0-1,-1):    
                    
        if 'TTTTT' in b[i][j][0]:
            b[i].pop(j)                
        elif 'GGTCTC' in b[i][j][0]:
            b[i].pop(j)
        elif 'GAGACC' in b[i][j][0]:
            b[i].pop(j)
        elif 'GGGGGG' in b[i][j][0]:
            b[i].pop(j)                
        elif 'AAAAAA' in b[i][j][0]:
            b[i].pop(j)
    if len(b[i])<2:
        b2empty+=1

###Sorting the crRNAs based on the scoring system
sortedscore=[]
for i in range (0,len(b)-1):
    k=[]
    score=[]
    score.append(b[i][0])
    k.append(b[i][0])
    a=sorting(b[i][1:])
    for j in a:
        k.append(j[1])
        score.append(j[0])
    b[i]=k
    sortedscore.append(score)
        
sequence=''        
empty=0        
#remove replicates from the list of crRNAs
for i in b:
    for j in range(len(i)-1,0-1,-1):
        if i.count(i[j])>1:
            i.pop(j)


##### Now the combining part. The crRNAs are all in a list called b.
###print 'len bc:  ' +str(len(bc))
with open("boutput.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(b)
f.close 
with open('scgenes.txt','r') as f:
    reader = csv.reader(f, delimiter="\t")
    d = list(reader)
p=[]
a=0
f.close

for i in range(0,len(d)-1):
    if d[i][0][0]=='>':
        p.append([])
j=-1
Firstone=True 
for i in range(0,len(d)-1):
    first_string= d[i][0][0]
    if first_string=='>':
        j=j+1
        p[j].append('#'+ d[i][0][1:13])
        Firstone=True
    elif Firstone:
        p[j].append(d[i][0])
        Firstone=False
    else:
        p[j][1]=p[j][1]+d[i][0]
counter=0
e=0
q=-1
k=-1
revcr=''
construct=[]
nf=0
testing=[]
for kos in b:
    construct.append([])
for item in range(0,len(b)-1):
    q=q+1
    construct[q].append(b[item][0])
    for item2 in p:       
                    
        if b[item][0][0:11]==item2[0][0:11]:
          
            for i in range(1,len(b[item])-1):
                if len(construct[q])<35:
                    if b[item][i] in item2[1]:
                        a=item2[1].find(b[item][i])
                        if b[item][i][21:23]=='GG':
                            if b[item][i][0:2]=='CC':
                                isitplus=ifplus(b[item][i])
                            else:
                                isitplus=True
                        else:
                            isitplus=False
                        
                            
                        if isitplus:
                            if a>32 and len(item2[1])-a>75:
                                k+=1 
                                construct[q].append([b[item][0],b[item][i] ,sortedscore[item][i],initialseq+item2[1][a-33:a+17]+item2[1][a+25:a+75]+b[item][i][0:20]+finalseq])
                                

                            elif len(construct[q])<10:
                                k+=1
                                raw_HA=HA(b[item][i][0:20],b[item][0].split(' ')[0],True)
                                if raw_HA!='':
                                    loc_gen=raw_HA.find(b[item][i][0:20])
                                    LHA=raw_HA[loc_gen-33:loc_gen+17]
                                    RHA=raw_HA[loc_gen+25:loc_gen+75]
                                    construct[q].append([b[item][0],b[item][i],sortedscore[item][i] ,initialseq+LHA+RHA+b[item][i][0:20]+finalseq])
                                    
                               
                        else:                             
                            if a>51 and len(item2[1])-a>56:
                                k+=1
                                revcr=revcom(b[item][i])
                                construct[q].append([b[item][0],b[item][i],sortedscore[item][i] ,initialseq+minusstrand(item2[1][a-52:a-2])+minusstrand(item2[1][a+6:a+56])+revcr[0:20]+finalseq])
                                
                            elif len(construct[q])<10:
                                k+=1
                                revcr=revcom(b[item][i])
                                raw_HA=HA(b[item][i][3:23],b[item][0].split(' ')[0],False)
                                if raw_HA!='':
                
                                    loc=raw_HA.find(b[item][i][0:20])
                                    LHA=raw_HA[loc-52:loc-2]
                                    RHA=raw_HA[loc+6:loc+56]
                                    construct[q].append([b[item][0],b[item][i],sortedscore[item][i] ,initialseq+minusstrand(LHA)+minusstrand(RHA)+revcr[0:20]+finalseq])              
                        
##bar=k                            
##
for i in range(0,len(construct)-1):    
    for j in range(len(construct[i])-1,-1,-1):    
                   
        if 'TTTTT' in construct[i][j][3][23:150]:
            construct[i].pop(j)
               
        elif 'GGTCTC' in construct[i][j][3][23:150]:
            construct[i].pop(j)
                
        elif 'GAGACC' in construct[i][j][3][23:150]:
            construct[i].pop(j)
               
        elif 'GGGGGG' in construct[i][j][3][23:150]:
            construct[i].pop(j)
               
        elif 'AAAAAA' in construct[i][j][3][23:150]:
            construct[i].pop(j)

                
                 



for k in range(0,len(construct)-1):
    if len(construct[k])<2:
        fails.append(construct[k][0])
               
summation=0            
for k in range (0,len(construct)-1):
    summation =summation+len(construct[k])
print (str(len(fails))+'empty')
k=b[2]

test=[['Gene Name','crRNA','score', 'Construct']]
###Testing for quality control

control=controlgen(100)
testcontrol=[]
for i in control:
    testcontrol.append(['Control',i[124:144]+'GGG', 'NA' , i])
counter=0
q=0
for i in range(0,len(construct)-1):
    for j in range(1,5):
        if len(construct[i])>j:
            test.append(construct[i][j])

k=1
for j in range(0,len(test)):
    if test[j][0]==test[j-k][0]:
        test[j][0]=' '
        k+=1
    else:
        k=1
test=test+testcontrol
with open("output.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(test)
f.close 




