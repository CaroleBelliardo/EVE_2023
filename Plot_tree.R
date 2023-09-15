#install.packages("ape")
#install.packages("phangorn")
#install.packages("phytools")
#install.packages("geiger")
#install.packages("devtools")
#install_github("liamrevell/phytools")
install.packages('dplyr')
library(ggplot2)
library(dplyr)
library(devtools)
library(ape)
library(RColorBrewer)

## read tree from string
text.string<-
  "((Limulus_polyphemus,(Centruroides_exilicauda,Parasteatoda_tepidariorum)),((Strigamia_maritima),(((Armadillidium_vulgare),((Locusta_migratoria),(((Oncopeltus_fasciatus), Acyrthosiphon_pisum),(((Bombus_terrestris), Apis_mellifera),((Nicrophorus_vespilloides, (Tribolium_castaneum, Diabrotica_virgifera)),
  ((Plutella_xylostella,(Heliconius_melpomene)),
  (Aedes_aegypti,(Musca_domestica,((Drosophila_virilis), Drosophila_melanogaster))))))))))));"

vert.tree<-read.tree(text=text.string)
# 
c=brewer.pal(n = 4, name = "PRGn")

## ---
vec=vert.tree$tip.label

d1=read.table('Bureau/clement_Avril2021/Len_eve_onlyG/sumEVE_exp-Noexp.txt2', sep='\t', stringsAsFactors=FALSE) #row.names = 1
names(d1)<-c( 'species', 'No expressed', 'Expressed' )
lab=c( 'No expressed', 'Expressed' )

data_new1 <- d1[match(vec, d1$species), ]

#barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,       legend =lab, beside=FALSE)



## plot type
d2_fg=read.table('Bureau/eve_mai2021/FG_sipi_formated.txt', sep='\t',
              header = TRUE, stringsAsFactors=FALSE) #row.names = 1
d2_fs=read.table('Bureau/eve_mai2021/FS_sipi_formated.txt', sep='\t',
              header = TRUE, stringsAsFactors=FALSE) #row.names = 1


fg <- d2_fg[match(vec, d2_fg$species), ]
fs <- d2_fs[match(vec, d2_fs$species), ]

#barplot(rbind(fg$piRNA, fg$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,       legend = colnames(fg[,4:5]), beside=FALSE,main = 'FG')

#barplot(rbind(fs$piRNA, fs$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
 #       legend = colnames(fs[,4:5]), beside=FALSE,main = 'FS')

##--- DRAWING

#--FS
par(mfrow=c(1,4))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)
barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE,  space=rep(0.8,19))


barplot(rbind(fs$piRNA, fs$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
       legend = colnames(fs[,4:5]), beside=FALSE,main = 'FS')

barplot(t(m), main='Viral families',horiz=TRUE, xpd=TRUE,las=2,cex.names=0.2,  ylab = '', xlab = '',axisnames=F,
        col = sc, beside=FALSE, space=rep(0.8,19))
legend("bottomright", colnames(m),pch=15,
       col=sc, cex=1.09)

#--FG
par(mfrow=c(1,4))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)
barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE,  space=rep(0.8,19))

barplot(rbind(fg$piRNA, fg$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(fg[,4:5]), beside=FALSE,main = 'FG')

barplot(t(m), main='Viral families',horiz=TRUE, xpd=TRUE,las=2,cex.names=0.2,  ylab = '', xlab = '',axisnames=F,
        col = sc, beside=FALSE, space=rep(0.8,19))
legend("bottomright", colnames(m),pch=15,
       col=sc, cex=1.09)

# -- FG FS
par(mfrow=c(1,4))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)
barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE,  space=rep(0.8,19))

barplot(rbind(fs$piRNA, fs$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(fs[,4:5]), beside=FALSE,main = 'FS')

barplot(rbind(fg$piRNA, fg$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(fg[,4:5]), beside=FALSE,main = 'FG')




#-- viral families 
#long version for ggplot # d1=read.table('Bureau/eve_mai2021/familles/FG_merged_longV.txt', sep='\t', header=TRUE, stringsAsFactors=FALSE)[2:4,] #row.names = 1

d3=read.table('Bureau/eve_mai2021/familles/FG_merged_sum_manual_completName.txt', sep='\t', header=TRUE) #row.names = 1
head(d3)
dim(d3)
data_new3 <- d3[match(vec, d3$species), ]

m=data.matrix(data_new3, rownames.force = NA)[,3:21]

# DRAW species
library(RColorBrewer)
n <- 21
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sc=sample(col_vector, n)

par(mfrow=c(1,3))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)

barplot(rbind(data_new1$`No expressed`, data_new1$Expressed),
        col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        beside=FALSE,  space=rep(0.8,19))
legend("bottomright", lab, pch=15,
       col=c( "grey", c[4]), cex=1.09)

barplot(t(m), horiz=TRUE, xpd=TRUE,las=2,cex.names=0.2,  ylab = '', xlab = '',axisnames=F,
         col = sc, beside=FALSE, space=rep(0.8,19))
legend("bottomright", colnames(m),pch=15,
       col=sc, cex=1.09)



#all
par(mfrow=c(1,5))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)

barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE,  space=rep(0.8,19))

barplot(rbind(fg$piRNA, fg$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(fg[,4:5]), beside=FALSE,main = 'FG')

barplot(rbind(fs$piRNA, fs$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(fs[,4:5]), beside=FALSE,main = 'FS')

barplot(t(m), main='Viral families',horiz=TRUE, xpd=TRUE,las=2,cex.names=0.2,  ylab = '', xlab = '',axisnames=F,
        col = sc, beside=FALSE, space=rep(0.8,19))
legend("bottomright", colnames(m),pch=15,
       col=sc, cex=1.09)



