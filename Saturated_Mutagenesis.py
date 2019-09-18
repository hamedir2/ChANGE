import csv
import random
import math
from difflib import SequenceMatcher
#Reverse complement function
initialseq='TATCTACACG'+'GGTCTC'+'ACCAAAAC'
finalseq='GTTTTAGAGT'+'GAGACC'+'AGCGTAACTC'   
all_codons=['Ala',
'Arg',
'Asn',
'Asp',
'Cys',
'Gln',
'End',
'Glu',
'Gly',
'His',
'Ile',
'Leu',
'Lys',
'Met',
'Phe',
'Pro',
'Ser',
'Thr',
'Trp',
'Tyr',
'Val'
]
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

#Gets the DNA sequence of the codon and retunrs its amino acid 
def find_codon(seq):
	codon_dict= {"GCT" : "Ala",
				"GCA" : "Ala",
				"GCC" : "Ala",
				"GCG" : "Ala",
				"AGA" : "Arg",
				"AGG" : "Arg",
				"CGT" : "Arg",
				"CGA" : "Arg",
				"CGC" : "Arg",
				"CGG" : "Arg",
				"AAT" : "Asn",
				"AAC" : "Asn",
				"GAT" : "Asp",
				"GAC" : "Asp",
				"TGT" : "Cys",
				"TGC" : "Cys",
				"TAA" : "End",
				"TGA" : "End",
				"TAG" : "End",
				"CAA" : "Gln",
				"CAG" : "Gln",
				"GAA" : "Glu",
				"GAG" : "Glu",
				"GGT" : "Gly",
				"GGA" : "Gly",
				"GGC" : "Gly",
				"GGG" : "Gly",
				"CAT" : "His",
				"CAC" : "His",
				"ATT" : "Ile",
				"ATA" : "Ile",
				"ATC" : "Ile",
				"TTG" : "Leu",
				"TTA" : "Leu",
				"CTA" : "Leu",
				"CTT" : "Leu",
				"CTG" : "Leu",
				"CTC" : "Leu",
				"AAA" : "Lys",
				"AAG" : "Lys",
				"ATG" : "Met",
				"TTT" : "Phe",
				"TTC" : "Phe",
				"CCA" : "Pro",
				"CCT" : "Pro",
				"CCC" : "Pro",
				"CCG" : "Pro",
				"TCT" : "Ser",
				"TCA" : "Ser",
				"AGT" : "Ser",
				"TCC" : "Ser",
				"AGC" : "Ser",
				"TCG" : "Ser",
				"ACT" : "Thr",
				"ACA" : "Thr",
				"ACC" : "Thr",
				"ACG" : "Thr",
				"TGG" : "Trp",
				"TAT" : "Tyr",
				"TAC" : "Tyr",
				"GTT" : "Val",
				"GTA" : "Val",
				"GTG" : "Val",
				"GTC" : "Val",}
	return codon_dict[seq]
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
#gets init codon (in AA) and returns the mostly codon with  nucleotides
def change_codon(init,highest):
	highest_usage={"Ala" : "GCT",
					"Arg" : "AGA",
					"Asn" : "AAT",
					"Asp" : "GAT",
					"Cys" : "TGT",
					"End" : "TAA",
					"Gln" : "CAA",
					"Glu" : "GAA",
					"Gly" : "GGT",
					"His" : "CAT",
					"Ile" : "ATT",
					"Leu" : "TTG",
					"Lys" : "AAA",
					"Met" : "ATG",
					"Phe" : "TTT",
					"Pro" : "CCA",
					"Ser" : "TCT",
					"Thr" : "ACT",
					"Trp" : "TGG",
					"Tyr" : "TAT",
					"Val" : "GTT",}
	second_highest={"Ala":"GCA",
					"Arg":"AGG",
					"Asn":"AAC",
					"Asp":"GAC",
					"Cys":"TGC",
					"End":"TGA",
					"Gln":"CAG",
					"Glu":"GAG",
					"Gly":"GGA",
					"His":"CAC",
					"Ile":"ATA",
					"Leu":"CTA",
					"Lys":"AAG",
					"Met":"ATG",
					"Phe":"TTC",
					"Pro":"CCT",
					"Ser":"TCA",
					"Thr":"ACA",
					"Trp":"TGG",
					"Tyr":"TAC",
					"Val":"GTA"}
	third_highest={"Asn": 'AAC',
					"Lys": 'AAG',
					"Thr": 'ACC',
					"Ser": 'AGT',
					"Ile": 'ATC',
					"Met": 'ATG',
					"His": 'CAC',
					"Gln": 'CAG',
					"Pro": 'CCC',
					"Arg": 'CGT',
					"Leu": 'CTA',
					"Asp": 'GAC',
					"Glu": 'GAG',
					"Ala": 'GCC',
					"Gly": 'GGC',
					"Val": 'GTG',
					"Tyr": 'TAC',
					"End": 'TAG',
					"Cys": 'TGC',
					"Trp": 'TGG',
					"Phe": 'TTC',}
	if highest==1:
		answer=highest_usage[init]
	if highest==2:
		answer=second_highest[init]
	if highest==3:
		answer=third_highest[init]
	return answer
def mutate(init):
	with open('Mutate.csv','r') as f:
	    reader = csv.reader(f, delimiter=",")
	    temp = list(reader)
	
	mutation_list=temp
	Amino=find_codon(init)
	for i in mutation_list:
		min_similarity=1
		if i[0]==Amino:
			
			for j in range(1,len(i)):
				if i[j]!='':					
					similarity=similar(init,i[j])
					if similarity<=min_similarity:
						min_similarity=similarity
						mutated=i[j]
	return mutated
print mutate('ATG')
#this function gets a gene and removes all BsaI sites from it
def BsaI_removal(seq):
	still_BsaI=True
	while (still_BsaI):
		BsaI_loc=seq.find('GGTCTC')
		# print BsaI_loc
		if BsaI_loc<0:
			BsaI_loc=seq.find(revcom('GGTCTC'))
		if BsaI_loc>0:
			mutation_loc=int(BsaI_loc/3+1)*3
			mutated_seq=mutate(seq[mutation_loc:mutation_loc+3])
			seq=seq[:mutation_loc]+mutated_seq+seq[mutation_loc+3:]
		else:
			still_BsaI=False
	return seq


def poly_removal(seq):
	is_poly=True
	fixed_location=0
	while is_poly:
		poly_loc=seq.find('TTTTT',fixed_location)
		if poly_loc>0:
			if seq[poly_loc-90:poly_loc+90].find('AAAAA')>0:
				randomizer=random.choice([0,1])
				seq=seq[:int(poly_loc/3+randomizer)*3]+mutate(seq[int(poly_loc/3+randomizer)*3:int(poly_loc/3+1+randomizer)*3])+seq[int(poly_loc/3+1+randomizer)*3:]
			else:
				fixed_location=poly_loc+5
		else:
			is_poly=False
	return seq

#this function finds the PAM sequences that are closer than 17 bp to the codon, if none is found, it looks for closer than 30 and if none is found, it'll look for closer than 40
def find_all_PAMS(location,gRNAlist):
	eligible_pam_list=[]
	for i in range(1,len(gRNAlist)-1):
		if gRNAlist[i][2]=='+':
			PAM_loc=int(gRNAlist[i][1])-3			
		else:
			PAM_loc=int(gRNAlist[i][0])-1
		if abs(PAM_loc-location)<17:
			eligible_pam_list.append(gRNAlist[i])
	if len(eligible_pam_list)==0:
		for i in range(1,len(gRNAlist)-1):
			if gRNAlist[i][2]=='+':
				PAM_loc=int(gRNAlist[i][1])-3			
			else:
				PAM_loc=int(gRNAlist[i][0])-1
			if abs(PAM_loc-location)<30:
				eligible_pam_list.append(gRNAlist[i])
		
	if len(eligible_pam_list)==0:
		for i in range(1,len(gRNAlist)-1):
			if gRNAlist[i][2]=='+':
				PAM_loc=int(gRNAlist[i][1])-3			
			else:
				PAM_loc=int(gRNAlist[i][0])-1
			if abs(PAM_loc-location)<40:
				eligible_pam_list.append(gRNAlist[i])	
	
	return eligible_pam_list

def sorting(gRNAlist,location):   
    global cRNA
    unsortedlist=[]
    for i in range(0,len(gRNAlist)):
    	z=int(gRNAlist[i][9])
    	GC_content=float(gRNAlist[i][4])/100
    	#calculate purine score
    	if gRNAlist[i][2]=='+':
    		next_to_PAM=gRNAlist[i][3][16:20]
    		AGS=0.25*next_to_PAM.count('G')+0.2*next_to_PAM.count('A')+0.15*next_to_PAM.count('C')
    		PAM_loc=int(gRNAlist[i][1])-1
    		dist_score=(17-abs(PAM_loc-location))/17
    		PAM_seq=gRNAlist[i][3][20:]
    		cRNA=gRNAlist[i][3][0:20]
    	if gRNAlist[i][2]=='-':
    		next_to_PAM=gRNAlist[i][3][3:7]
    		AGS=0.25*next_to_PAM.count('C')+0.2*next_to_PAM.count('T')+0.15*next_to_PAM.count('G')
    		PAM_loc=int(gRNAlist[i][0])+1
    		dist_score=(17-abs(PAM_loc-location))/17
    		PAM_seq=gRNAlist[i][3][0:3]
    		cRNA=gRNAlist[i][3][3:]
    	score=100*(0.33333*GC_content+0.33333*AGS+0.33333*dist_score)/(z**2)
    	unsorted_gRNA=gRNAlist[i]
    	unsorted_gRNA.append(score)
    	unsorted_gRNA.append(PAM_loc)
    	unsorted_gRNA.append(PAM_seq)
    	unsorted_gRNA.append(cRNA)
    	unsortedlist.append(unsorted_gRNA)
    unsortedlist.sort(key=lambda tup: tup[11])
    answer=[]
    for i in unsortedlist:
        answer.append(i)
    answer.reverse()
    return answer
#Main code

f=open("flankin_UBC4.txt","r") 
gene_seq=f.read()
cRNA=''
f=open("UBC4_short.txt","r") 
gene_seq_short=f.read()

with open('UBC4_gRNA.csv','r') as f:
    reader = csv.reader(f, delimiter=",")
    temp = list(reader)
gene_seq= poly_removal(gene_seq)
No_PAM=0
gene_seq=BsaI_removal(gene_seq)
finallist=[['cRNA', 'donor','codon location','PAM location','full gRNA', 'old codon', 'new codon']]

for codon in range (20,int(len(gene_seq_short)/3)+20):
	for AA in all_codons:
		count=0
		not_confirmed=True
		current_location=codon*3
		if_no_PAM=False
		
		while not_confirmed:			
			
			all_PAMS=find_all_PAMS(current_location+151,temp)
			if len(all_PAMS)>0:
				sorted_PAMS=sorting(all_PAMS,current_location+151)
				cRNA=sorted_PAMS[0][14]
			else:
				print 'No PAMS found'
				No_PAM+=1
				if_no_PAM=True
				break
			current_codon=gene_seq[current_location:current_location+3]
			
			if sorted_PAMS[0][2]=='+':
				PAM_location=int(sorted_PAMS[0][1])-3+60-210
				print gene_seq[PAM_location:PAM_location+3]
			else:
				PAM_location=int(sorted_PAMS[0][0])-1+60-210
			distance=PAM_location-current_location
			if distance>0:
				if sorted_PAMS[0][2]=='+':				
					donor=gene_seq[current_location-40:current_location]+change_codon(AA,1)+gene_seq[current_location+3:PAM_location+43]
					off_set= PAM_location % 3
					# PAM= donor[distance+40:distance+43]
					if off_set==0:
						donor=donor[:40+int((distance)/3)*3]+mutate(donor[40+int((distance)/3)*3:40+int((distance)/3+1)*3])+donor[40+int((distance)/3+1)*3:]					
						# print 40+int((distance)/3)*3
						if donor[distance+41:distance+43]=='GG' or donor[distance+41:distance+43]=='AG':
							donor=donor[:40+int((distance)/3-1)*3]+mutate(donor[40+int((distance)/3-1)*3:40+int((distance)/3)*3])+donor[40+int((distance)/3)*3:]						
					else:
						donor=donor[:40+int((distance)/3)*3]+mutate(donor[40+int((distance)/3)*3:40+int((distance)/3+1)*3])+mutate(donor[40+int((distance)/3+1)*3:40+int((distance)/3+2)*3])+donor[40+int((distance)/3+2)*3:]
						if donor[distance+41:distance+43]=='GG' or donor[distance+41:distance+43]=='AG':
							donor=donor[:40+int((distance)/3-1)*3]+mutate(donor[40+int((distance)/3-1)*3:40+int((distance)/3)*3])+donor[40+int((distance)/3)*3:]		
									
				else:			
					donor=gene_seq[current_location-40:current_location]+change_codon(AA,1)+gene_seq[current_location+3:PAM_location+43]
					off_set= PAM_location % 3
					if off_set==0:
						donor=donor[:40+int((distance)/3)*3]+mutate(donor[40+int((distance)/3)*3:40+int((distance)/3+1)*3])+donor[40+int((distance)/3+1)*3:]
						if donor[distance+40:distance+42]=='CC' or donor[distance+40:distance+42]=='CT':
							donor=donor[:40+int((distance)/3+1)*3]+mutate(donor[40+int((distance)/3+1)*3:40+int((distance)/3+2)*3])+donor[40+int((distance)/3+2)*3:]	
					else:
						donor=donor[:40+int((distance)/3)*3]+mutate(donor[40+int((distance)/3)*3:40+int((distance)/3+1)*3])+mutate(donor[40+int((distance)/3+1)*3:40+int((distance)/3+2)*3])+donor[40+int((distance)/3+2)*3:]
						if donor[distance+40:distance+42]=='CC' or donor[distance+40:distance+42]=='CT':
							donor=donor[:40+int((distance)/3+2)*3]+mutate(donor[40+int((distance)/3+2)*3:40+int((distance)/3+3)*3])+donor[40+int((distance)/3+3)*3:]		
			#Finding the donor if the codon is downstream of the PAM
			else:
				if sorted_PAMS[0][2]=='+':
					donor=gene_seq[PAM_location-40:current_location]+change_codon(AA,1)+gene_seq[current_location+3:current_location+43]
					off_set= PAM_location % 3
					PAM_codon=len(donor)+int(distance/3)*3-40
					if off_set==0:					
						donor=donor[:40]+mutate(donor[40:43])+donor[43:]
						if donor[41:43]=='GG' or donor[41:43]=='AG':
							donor=donor[:37]+mutate(donor[37:40])+donor[40:]							
					else:
						donor=donor[:PAM_codon-3]+mutate(donor[PAM_codon-3:PAM_codon])+mutate(donor[PAM_codon:PAM_codon+3])+donor[PAM_codon+3:]					
						if donor[41:43]=='GG' or donor[41:43]=='AG':						
							donor=donor[:PAM_codon-6]+mutate(donor[PAM_codon-6:PAM_codon-3])+donor[PAM_codon-3:]							
				else:
					donor=gene_seq[PAM_location-40:current_location]+change_codon(AA,1)+gene_seq[current_location+3:current_location+43]
					off_set= PAM_location % 3 
					PAM_codon=len(donor)+int(distance/3)*3-40
					if off_set==0:
						donor=donor[:40]+mutate(donor[40:43])+donor[43:]
						if donor[40:42]=='CC' or donor[40:42]=='CT':									
							donor=donor[:43]+mutate(donor[43:46])+donor[46:]					
					else:
						donor=donor[:PAM_codon-3]+mutate(donor[PAM_codon-3:PAM_codon])+mutate(donor[PAM_codon:PAM_codon+3])+donor[PAM_codon+3:]	
						if donor[40:42]=='CC' or donor[40:42]=='CT':						
							donor=donor[:PAM_codon+3]+mutate(donor[PAM_codon+3:PAM_codon+6])+donor[PAM_codon+6:]	
					if len(donor)>100:
						trim_length=int(math.ceil((len(donor)-100)*1.0/2))
						donor=donor[trim_length:len(donor)-trim_length]
			#revcom if negative strand
			if sorted_PAMS[0][2]=='-':
				donor=revcom(donor)
				cRNA=revcom(cRNA)
			if 	'TTTTT' in donor:
				donor=revcom(donor)
			donor=BsaI_removal(donor)
			if 'GGTCTC' in donor or revcom('GGTCTC') in donor:
				not_confirmed=True
				sorted_PAMS.pop(0)
			else:
				not_confirmed=False
			count+=1			
			if count>5:
				not_confirmed=False
				print 'oh noooo'
		#making grna with variable length
		
		partial_gRNA=initialseq+donor+cRNA+finalseq
		seq=''
		#adding random nucleotides to the 3' side to make the length homogenous
		for i in range(1,170-len(partial_gRNA)+1):
			seq=seq+random.choice('ACGT')
		gRNA=partial_gRNA+seq
		
		if if_no_PAM==False:
			finallist.append([cRNA,donor,current_location-60+211, PAM_location-60+211,gRNA,find_codon(current_codon), AA])
		else:
			print 'poogh'
with open("output.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(finallist)