import pprint   # debugg -- affichage Dump des structures de données
import os       # manip systeme
import sys      # gestion arguments et opt
import pandas as pd 

df=pd.read_csv('Classified-goodSpeciesNames.txt', sep = '\t')
df = df.rename(columns={'2':'species','3':'accViral','16':'ScaffoldName','12':'starts','13':'stops', '24':'virus'})

print(df.shape)


## -- virus family
print('start')
distance = 50
output=df[df.species ==''].copy() # copy empty tab

for Scaffold in df.ScaffoldName.unique(): # pour chaque scaffold
    EVEscontig=df.loc[df.ScaffoldName == Scaffold].reset_index(drop=True)
    EVEscontig=EVEscontig.sort_values(by=['starts'])
    
    if len(EVEscontig.ScaffoldName) == 1 : # une seule occurence eve id dans le tableau d'origine 
        output=output.append(EVEscontig)#.reset_index(drop=True) # copy line dans tab de sortie
        next
        
    else: # si plusieurs eve avec meme id
        #print('else1')
        
        line1=EVEscontig.iloc[0,:]
        output=output.append(line1).reset_index(drop=True) # ajoute première ligne dans la table 
        EVEscontig=EVEscontig.drop(0,axis=0) # remove first line 
        
        for i in EVEscontig.index : # EVEs fichier In 

            lineI = EVEscontig.loc[i,:]
            startI = lineI['starts']
            stopI = lineI['stops']
            virAccI = lineI['accViral']
            virusI = lineI['virus']
            outEvesI = output.loc[output['ScaffoldName'] == Scaffold]
            statI = 'not'
            
            
            for i_outEve in outEvesI.index : # EVEs fichier out
                if statI == 'yes' : break

                line_output = outEvesI.loc[i_outEve,:]
                
                startO = line_output['starts']
                stopO = line_output['stops']
                virAccO = line_output['accViral']
                virusO = line_output['virus']
                
                if virusI == virusO : # meme virus identitque
                    if ((startI   > startO) and  ((startI ) <= stopO  + distance)) :#((stopO   > startL) and 

                        indd = output[(output['ScaffoldName'] == Scaffold) & (output['starts'] == startO) ]
                        ind=indd.index
                        output.loc[ind,['stops']] = stopI
                        statI = 'yes'
                        
                      #  print(Scaffold,'****',startI,stopI)
                       # print(Scaffold,'****',startO,stopO)
                        #print('_________________________________')
                        
                        
                        break
                        #print(indd)
            if statI == 'not' : output=output.append(lineI).reset_index(drop=True)      
             
print('end')
#output.to_csv('fusionEve_50pb_virFamily.csv', sep ='\t', index=None)
#output.to_csv('fusionEve_10pb_virFamily.csv', sep ='\t', index=None)


 ## -- Virus accession number
  print('start')
distance = 50
output=df[df.species ==''].copy() # copy empty tab

for Scaffold in df.ScaffoldName.unique(): # pour chaque scaffold
    EVEscontig=df.loc[df.ScaffoldName == Scaffold].reset_index(drop=True)
    EVEscontig=EVEscontig.sort_values(by=['starts'])
    
    if len(EVEscontig.ScaffoldName) == 1 : # une seule occurence eve id dans le tableau d'origine 
        output=output.append(EVEscontig)#.reset_index(drop=True) # copy line dans tab de sortie
        next
        
    else: # si plusieurs eve avec meme id
        #print('else1')
        
        line1=EVEscontig.iloc[0,:]
        output=output.append(line1).reset_index(drop=True) # ajoute première ligne dans la table 
        EVEscontig=EVEscontig.drop(0,axis=0) # remove first line 
        
        for i in EVEscontig.index : # EVEs fichier In 

            lineI = EVEscontig.loc[i,:]
            startI = lineI['starts']
            stopI = lineI['stops']
            virAccI = lineI['accViral']
            virusI = lineI['virus']
            outEvesI = output.loc[output['ScaffoldName'] == Scaffold]
            statI = 'not'
            
            
            for i_outEve in outEvesI.index : # EVEs fichier out
                if statI == 'yes' : break

                line_output = outEvesI.loc[i_outEve,:]
                
                startO = line_output['starts']
                stopO = line_output['stops']
                virAccO = line_output['accViral']
                virusO = line_output['virus']
                
                if virAccI == virAccO : # meme virus identitque
                    if ((startI   > startO) and  ((startI ) <= stopO  + distance)) :#((stopO   > startL) and 

                        indd = output[(output['ScaffoldName'] == Scaffold) & (output['starts'] == startO) ]
                        ind=indd.index
                        output.loc[ind,['stops']] = stopI
                        statI = 'yes'
                        
                      #  print(Scaffold,'****',startI,stopI)
                       # print(Scaffold,'****',startO,stopO)
                        #print('_________________________________')
                        
                        
                        break
                        #print(indd)
            if statI == 'not' : output=output.append(lineI).reset_index(drop=True)      
             
print('end')

output.to_csv('fusionEve_50pb_virAcc.csv', sep ='\t', index=None)
#output.to_csv('fusionEve_10pb_virAcc.csv', sep ='\t', index=None)



print(df.shape, output.shape)
