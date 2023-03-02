#!/usr/bin/python3.6.4
# -*- coding: utf-8 -*-

#############################
##          Import         ##
#############################
import pprint
import os
import sys
from optparse import OptionParser
import argparse
import subprocess

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
usage = "usage: python3 eve_mapping.py -f 1000 -g /mnt/65To/Carole/1-data/genomes_transcriptomes/Apis_GT.fna -pb $Apis.bed -sB /mnt/65To/Carole/2-mapping/smalls/s_i_BAM/Apis_FG.bam #-lB /mnt/65To/Carole/2-mapping/longs/s_i_BAM/Apis_FG.bam -o Apis_FG.tab  --h\' : usage parametre informations"
if ((len(sys.argv) == 2) and (sys.argv[1] == '-h')) :
    print(usage)
    sys.exit(0) 

parser.add_argument('-f', '--flancSize', type=int, help='size of EVE side  \n', default=1000)
parser.add_argument('-g', '--genome', type=str, help='genome (.fasta file) \n')
parser.add_argument('-pb', '--ori_bed', type=str, help='Eve annotation file (.bed file) \n')
#parser.add_argument('-lB', '--lB', type=str, help='mapping of long RNA (.bam file) \n')
parser.add_argument('-sB', '--sB', type=str, help='mapping of small RNA (.bam file) \n' )
parser.add_argument('-o', '--output', type=str, help=' output (.tab file) \n', default="output.table")

args = parser.parse_args()   

args.lenC=args.genome+'.len'

# ###########################
# ##       FUNCTIONS       ##
# ###########################   
# gestion erreur
def errorArgsBed(bed_init,bed_fin):
    b_bedCreate=False # flag contigLen
    try: ## Bed annotation EVEs
        with open(bed_init): pass
        print( 'Fichier '+bed_init+' existe OK!\n creation fichier'+bed_fin+' possible... ' )
        b_bedCreate=True ## flag le fichier bed ( avec flanq n'existe pas ) mais le fichier ori_bed ( sans flanq ) existe 
    except IOError :
        print( 'Fichier '+bed_init+'n existe pas; erreur fatale !!!!!!\n======> INDIQUER LE CHEMIN DU FICHIER BED INDIQUAND LES REGIONS GENOMIQUES DES EVEs\n' )
        sys.exit(0)
    return b_bedCreate

# contig size
def f_LenContig (lenC,genome,ori_bed) :
    LenContig={}
    try: ## recup info taille des contigs
        with open(lenC): pass
        print( 'Fichier '+lenC+' existe déjà OK!****' )
    except IOError :
        print ('**Le fichier :'+lenC+' existe pas et va etre... \n')
        try: ## Creation fichier taille contig à partir du fichier genome taille des contigs
            with open(genome): pass # test si fichier genome existe
            print('**Fichier '+genome+' existe OK!\n creation fichier'+lenC+' possible... ')
            os.system('awk \'/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}\' '+genome+' |  tr \"\n$\" \"\t\" | tr -s \'>\' \'\n>\' |sed "1d" > '+args.lenC)
            #os.system("sed -i '1d' "+lenC )
        except IOError :
            print( 'Fichier '+genome+'n existe pas FATAL ERROR !!!!!!\n======> INDIQUER LE CHEMIN DU FICHIER FASTA CONTENANT LES SEQUENCES GENOMIQUES\n' )
            sys.exit(0)
          
    #-- stock contig length in dico
    os.system('cut -f 1 ' +ori_bed+' > '+ori_bed+'.cut ; grep -w -f '+ori_bed+'.cut '+lenC+' >'+ lenC+'.cut' )        
    LenC = open(lenC+'.cut', "r")    
    for ligne in LenC :
        ligne=ligne.rstrip() # suppr. saut de ligne
        l=ligne.split("\t")
        LenContig[l[0]]=int(l[1])
    LenC.close()
    os.remove(ori_bed+'.cut')
    os.remove(lenC+'.cut')
    return LenContig;

# create mpileup file with bam file
def f_mpileup (bam,bed,genome) :
	try:## Bed annotation EVEs + Flanq
		with open(bam): pass
		os.system( '/opt/samtools-1.9/bin/samtools mpileup -f '+genome+' -l '+bed+' '+bam+' > '+bam+'.mpileup')
		os.system("cut -f 1,2,5 " +bam+'.mpileup'+' > '+bam+'.mpileup.cut' ) ##
		print('small mpileup ok! ')
	except IOError :
		print ('Le fichier :'+bam+'.n existe pas FATAL ERROR !!!!!!\n======> INDIQUER LE CHEMIN DU FICHIER BAM, alignement de  ARN \n')
		sys.exit(0)


 # Recup info  mpileup
def f_mpileupInfo(Dico,bam,typ) :
    Mpileup = open(bam+'.mpileup.cut' , "r")     
    for ligne in Mpileup :
       # print(ligne)
        ligne=ligne.rstrip() # suppr. saut de ligne
        l=ligne.split("\t")
        ## stockage des infos  
        if (len(l) > 2 ) :
            if (l[0], int(l[1]), typ+'_forw') not in Dico: 
                Dico[l[0], int(l[1]), typ+'_forw'] = str(l[2].count('.') + l[2].count('>'))
            if(l[0], int(l[1]), typ+'_rev') not in Dico :
                Dico[l[0], int(l[1]), typ+'_rev']  = str(l[2].count(',') + l[2].count('<'))
        else :
            if (l[0], int(l[1]), typ+'_forw') not in Dico: 
                Dico[l[0], int(l[1]), typ+'_forw'] = '0'
            if (l[0], int(l[1]), typ+'_rev') not in Dico :  
                Dico[l[0], int(l[1]), typ+'_rev']  = '0'
      #  print(str(Dico[l[0], int(l[1]), typ+'_forw']) + ' '+str(Dico[l[0], int(l[1]), typ+'_rev']))
    Mpileup.close()
    return (Dico)

# small and long mpileup file info in dico 
def f_table(bed,output,Dico):
    test={}
    Bed = open(bed, "r")  
    OutFile = open(output, "w")
    OutFile.write('contig \t position \t fill \t count \n')        
    info={}
    for ligne in Bed : ## parcours le fichier bed
        ligne=ligne.rstrip() # suppr. saut de ligne
        l=ligne.split("\t")
        contig=l[0]
        r= range(int(l[1]), int(l[2])) # liste des position de l eve
        for posi in r :
            
            ##long
            #if (l[0],posi) not in test :
             #   test[l[0],posi]=''
              #  if (l[0],posi,'long_forw')  in Dico.keys() :
               #     info['l_forw'] = Dico[l[0],posi,'long_forw']
             #
              #  else :
               #     info['l_forw'] = '0'

                #if (l[0],posi,'long_rev')  in Dico.keys() :
                 #   info['l_rev'] = Dico[l[0],posi,'long_rev']
                #else :
                 #   info['l_rev'] = '0'
                ##small
                if (l[0],posi,'small_forw')  in Dico.keys() :
                    info['s_forw'] = Dico[l[0],posi,'small_forw']
                else :
                    info['s_forw'] = '0'
                
                if (l[0],posi,'small_rev')  in Dico.keys() :
                    info['s_rev'] = Dico[l[0],posi,'small_rev']
                else :
                    info['s_rev'] = '0'
               # print('contig :'+contig+' posi '+str(posi)+' '+str(l[1])+' '+str(info['s_rev']) + ' '+str(info['s_forw']))

               # OutFile.write(l[0] +'\t'+ str(posi) +'\t'+ info['l_forw'] +'\t'+ info['l_rev']+'\t'+ info['s_forw'] +'\t'+ info['s_rev'] +'\n')
                info['contig'] = l[0]
                info['position'] = str(posi)
                OutFile.write(contig +'\t'+ str(posi)+'\t'+ 's_forw'+'\t'+ info['s_forw']+'\n')
                OutFile.write(contig +'\t'+ str(posi)+'\t'+ 's_rev'+'\t-'+ info['s_rev']+'\n')
                #OutFile.write(contig +'\t'+ str(posi)+'\t'+ 'l_forw'+'\t'+  info['l_forw']+'\n')
                #OutFile.write(contig +'\t'+ str(posi)+'\t'+ 'l_rev'+'\t-'+ info['l_rev'] +'\n')
    print('==> mpileup.table OK !!! \n outputFile :: '+ args.output)
    Bed.close()
    OutFile.close()
    #return(Dico)

    
# Main
def f_bamToTable(sB,lB,bed,genome,output) :
    f_mpileup (sB,bed,genome) ## test input mpileup
    #f_mpileup (lB,bed,genome)
    Mpileup_dicoo ={}
    Mpileup_dicoo=f_mpileupInfo (Mpileup_dicoo,sB,'small')
    
#    Mpileup_dicoo=f_mpileupInfo (Mpileup_dicoo,lB,'long') # mpileup in dico
    f_table(bed,output,Mpileup_dicoo)
    

##################################################################################################################################################################
#########################                                          ***  TABLES    ***                                                  ###########################
##################################################################################################################################################################

###########################
##       EACH  EVE       ##
###########################
b_bedCreate=errorArgsBed(args.ori_bed,args.ori_bed+'.flanc')
if b_bedCreate == True : 
	LenContig = f_LenContig(args.lenC,args.genome,args.ori_bed)

# creation un fichier  bed avec flanq    PreBed = open(ori_bed, "r") 
PreBed = open(args.ori_bed, "r")
Bed = open(args.ori_bed+'.flanc', "w")
for ligne in PreBed :
	ligne=ligne.rstrip() # suppr. saut de ligne
	l=ligne.split("\t")
	if int(l[1]) > int(l[2]) :
		debF=int(l[2])-args.flancSize
		finF=int(l[1])+args.flancSize
		l[3]='-'
	else :
		debF=int(l[1])-args.flancSize
		finF=int(l[2])+args.flancSize
	if(debF < 0) : debF = 0
	if (finF > LenContig[l[0]]) and (debF < LenContig[l[0]]): finF = LenContig[l[0]]
	Bed.write(l[0]+'\t'+str(debF)+'\t'+str(finF)+'\t'+l[3]+'\t'+l[4]+'\t'+l[5]+'\t'+l[6]+'\n')
PreBed.close()
Bed.close()
args.bed=args.ori_bed+'.flanc'
f_bamToTable(args.sB,args.lB,args.bed,args.genome,args.output)


