#!/usr/bin/Rscript
install.packages("reshape");
install.packages("extrafont")
install.packages("wesanderson");
install.packages("ggplot2")
install.packages("gridExtra")
vignette("xkcd-intro")
install.packages('jpeg')
#setwd('Bureau/clement_Avril2021/')

require(jpeg)

df <- data.frame(stringsAsFactors=FALSE)
dl <- data.frame(stringsAsFactors=FALSE)

library(reshape); library(extrafont); library(wesanderson); library("ggplot2");library("gridExtra");library(xkcd)


# load("temp.RData") 

sevEves = c()
flanq = 1000
args <- commandArgs(TRUE)

#args=c('4_small_long_expTAB/FS/Avul_FS.tab',  'bed_GT/strandedBed_new_familyAdded/Avul_stranded.txt','bed_GT/strandedBed_new_familyAdded/Avul_stranded.txt', '4_genomes_transcriptomes_len/Avul_GT.fna.len', 'Plots_2021_v4/test.pdf', '4_ET_out_filtre/Avul.out', 'Avul', 'FS', 'bed_exp_eve/small/Avul.txt' ,'bed_exp_eve/long/Av')

if (length(args) < 5) {
print('Rscript --vanilla Figs.R *deep_table *bed_combin *bed *genomeLen *outFile *teFile ') }
if (length(args) >= 5) {
#imports
	table_FG <- read.table(args[1], header = TRUE, sep = "\t", dec = ".", strip.white = TRUE)
	EVE_combin <- read.table(args[2],header = TRUE, sep = "\t", strip.white = TRUE)
	 names(EVE_combin) <- c("Contig", "start", "stop","stranded","id%",'evalue', "lineage","family")
	EVE <- read.table(args[3],header = TRUE, sep = "\t", strip.white = TRUE)
	 names(EVE) <- c("Contig", "start", "stop","stranded","id%",'evalue', "lineage","family")
	ContigLen <- read.table(args[4],header = FALSE, sep = "\t", strip.white = TRUE)[,1:2]
	 names(ContigLen) <- c("Contig", "len")
	outFile <- args[5]
	te = F
	if ((length(args) >= 6) && (file.exists(args[6]))) {
	  te = T
	  TE <- read.table(args[6], header = F, sep = "\t")
	  names(TE) <- c("Contig", "start", "stop", "sens", "type")
	  te = T
	  arrowCol = c( 'red3', 'cyan2', 'skyblue1',  'blue', 'orange', 'yellow')
	  names(arrowCol) <- c("ClasseII","LTR" ,"nonLTR", "Retro-other", "Satellite","Unknown") #"Low_complexity","Simple_repeat"
	}
	
	speciesAbbrev = args[7]
	condition = args[8]
	
	legende <- 0
	pdf(outFile)
	# par(family='CM Roman')
	par(mfrow = c(2,1));

	mf = 2
	par( xpd = T) # 2 images
	par(mar = c(4, 4, 4, 4))
	gb2 <- c(alpha("blueviolet", 0.6 ),alpha("blueviolet", 0.6),alpha("green2", 0.45),alpha("green2", 0.45))
	names(gb2) <- c('s_forw','s_rev','l_forw','l_rev')

#---
	EVE_combin = EVE_combin[order(EVE_combin$family, EVE_combin$Contig,  EVE_combin$start, decreasing = F),]
	for (i in 1:length(EVE_combin$Contig)) { ## parcours tous les Contigs 
		if (legende == 1) {par(mfrow = c(2,1)); mf = 2}
	  legende = legende + 1
		if (legende == 1 ) {  
		plot.new()
		legend(0.225,1.2,legend = c('small','long'),
		       fill = c(gb2[1], gb2[3]), ncol = 1,horiz = FALSE,xpd = TRUE, bty = "n",
		       cex = 1)
		if (te == T) {
		  legend(0.25,0.7,legend = names(arrowCol),lty = 1,col = arrowCol,
		         ncol = 1,horiz = FALSE,xpd = TRUE, bty = "n",
		         cex = 1)}
		}
		
		contig = as.character(EVE_combin[i,1])
		Start = EVE_combin[i,2] 
		Stop = EVE_combin[i,3] 
		startR = EVE_combin[i,2] - 1000 
		stopR = EVE_combin[i,3] + 1000
		vir = EVE_combin$lineage[i] ## pb
		
		
		
		## test
		contig_depth = table_FG[which(as.character(table_FG$contig) == contig
				& table_FG$position > startR 
				& table_FG$position < stopR),] # pro/locus de l'EVE
		
		eve_depth = table_FG[which(as.character(table_FG$contig) == contig
		                             & table_FG$position > Start 
		                             & table_FG$position < Stop),]
		
		eve_small = eve_depth[which(as.character(eve_depth$contig) == contig  
		                             & (as.character(eve_depth$fill) == 's_forw' 
		                             | as.character(eve_depth$fill) == 's_rev')),] # pro/locus de l'EVE
		
	
		eve_long = eve_depth[which(as.character(eve_depth$contig) == contig  
		                            & (as.character(eve_depth$fill) == 'l_forw' 
		                               | as.character(eve_depth$fill) == 'l_rev')),] # pro/locus de l'EVE
		
		contig_depth = contig_depth[order(contig_depth$position, decreasing = F),]
	  	## test
	  
		### DISPLAY SELON NB EVE DANS CONTIG
			if (is.data.frame(contig_depth) && nrow(eve_small) < 100) {next} #
		  	if (any(abs(eve_small$count) > 10) == TRUE) {
        
		  	  df = rbind(df,EVE_combin[i,])	  	
		  	  
		  	  if (any(abs(eve_long$count) > 3) == TRUE) {
		  	    dl = rbind(dl,EVE_combin[i,])
		  	  }
		  	  
				starts = EVE$start[which(as.character(EVE$Contig) == contig 
		      & EVE$start >= Start & EVE$stop <= Stop) ];starts = sort(starts, decreasing = FALSE) # recupere ts les starts des EVE du contig
				stops = EVE$stop[which(as.character(EVE$Contig) == contig
		      & EVE$start >= Start & EVE$stop <= Stop)]; stops = sort(stops, decreasing = FALSE)    # recupere ts les stops des EVE du contig
	legende = legende + 1
	vec = c('s_forw','s_rev','l_forw','l_rev') #,
	xl = c(startR ,stopR)
    
	if (length(unique(starts)) > 1) {
	    sevEves = c(sevEves,contig)
    	next}
	
	
	cexLab = 1
	cexAxis = 1
	
    if (xl[1] < 0 ) {xl[1] = 0}
	if (xl[2] < 0 ) {xl[2] = 0} # gestion erreur plot size vs genome size
    
    m <- (min(contig_depth$count) + 1)*110/100
    if (min(contig_depth$count) > -(max(contig_depth$count)/3)) { m <- -(max(contig_depth$count)/2) }
	yl = c(min(contig_depth$count) + m,max(contig_depth$count) + (max(contig_depth$count)*5/100))
  	y = 0; x = 0
## PLOTS
  	yy = max(abs(contig_depth$count))
 		plot( y ~ x,
			las = 1, type = 'n', xpd = FALSE,
			frame = F, xlab = "",ylab = "", ylim = c(-(yy + (yy*0.4)), yy + (yy*0.4)),
			cex.lab = cexLab,cex.axis = cexAxis,
			xlim = xl #,ylim = c(yl[1],yl[2] + (yl[2]*0.4))
			)
 		

 		for (v in 1:length(vec)) {	# pour chaque type arn plots
 		  x = contig_depth$position[which(contig_depth$fill == vec[v]  )]
 		  y = contig_depth$count[which(contig_depth$fill == vec[v] )]
 		  y = c(0,y,0) 
 		  x = c(x[1] - 1,x,x[length(x)] + 1)
 		  polygon( y ~ x,
 		           las = 1, col = gb2[v], border = NA ,
 		           cex.lab = 2,cex.axis = cexAxis,lwd = 0.7, main = paste(contig,'\n'),
 		           xlim = xl, ylim = yl)
 		}

 		if (te == T) {
 		  contig_TE = TE[which(as.character(TE$Contig) == contig
 		                      & ((TE$start > (startR)  & (TE$start < (stopR) ) # &  TE$stop > (startR)
 		                            | ((TE$stop > (startR) &  TE$stop < (stopR)  ) )
 		                      ))),] # pro/locus de l'EVE
 		  
 		  if (te == T) {
 		    if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) { 
 		      contig_TE$start[contig_TE$start < startR & contig_TE$stop > startR  ] <- startR #+1000
 		      contig_TE$stop[contig_TE$stop > stopR & contig_TE$start < stopR] <- stopR
 		    }
 		  }
 		}
 		
 		
		for (j in 1:length(starts)) { ## pour chaque EVEs du contig	
		
			st =	starts[j]		
			sp = stops[j]
			nameEVE = EVE$EVE[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
				& EVE$start == starts[j]
				& EVE$stop == stops[j] )]
			lineageEVE = EVE$family[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
			                      & EVE$start == starts[j]
			                      & EVE$stop == stops[j] )]
			
			eSens = EVE$stranded[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
				& EVE$start == st
				& EVE$stop == sp )] 
			  tmp = sp
			
			  
      			  
		  if (te == T) {
		    if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) { 
		      te_todel = contig_TE[contig_TE$start >= st  & contig_TE$stop <= sp,  ]
		      contig_TE = contig_TE[!(contig_TE$start %in% te_todel$start & contig_TE$stop %in% te_todel$stop ),]
		      
		      te_todel = contig_TE[contig_TE$start <= st  & contig_TE$stop >= sp,  ]
		      contig_TE = contig_TE[!(contig_TE$start %in% te_todel$start & contig_TE$stop %in% te_todel$stop ),]
		      
		      contig_TE$start[contig_TE$start > st  & contig_TE$start < sp] <- sp + 2
		      contig_TE$stop[contig_TE$stop > st & contig_TE$stop < sp] <- st - 3
		    }
		  }
			  
			#if (as.character(eSens) == '-'){ xl=rev(xl);contig_depth$count<- - contig_depth$count }
			if (as.character(eSens) == '-') {
			  tmp = sp
			  sp = st
			  st = tmp
			}
			  
			  
			 ## distrib taille *****************************************************************************************
			 
			   distrib = paste('/home/cbelliardo/Bureau/clement_Avril2021/distrib_jpg/',speciesAbbrev,'/',contig,'__',starts[j],'__',stops[j],'_',condition,'_run1_5prime.jpg', sep = '')
			  
			 if (file.exists(distrib)){
			  img <- readJPEG(distrib)
			  gauche = stops[j] 
			  bas = max(contig_depth$count)/10
			  droite = stopR
			  haut = yy*2
			  rasterImage(img,gauche,bas,droite,haut)
			 }
			  #addImg(img, x =stopR, y = yy , width = 1400)
			  

			  ## distrib taille end  *****************************************************************************************
			
# 			arrows(st, yl[1] - yl[1]*10/100, x1 = sp, y1 = yl[1] - yl[1]*10/100, length = 0.045, angle = 40,
#        			code = 2, col = par("fg"), lty = par("lty"),
#        			lwd = par("lwd"))
			  
			  arrows(st, -yy - yy*10/100, x1 = sp, y1 = -yy - yy*10/100, length = 0.045, angle = 40,
			         code = 2, col = par("fg"), lty = par("lty"),
			         lwd = par("lwd"))
			
			segments(st,(max(contig_depth$count)*10/100), # vertical start
		             x1 = st,
		             y1 = min(contig_depth$count)*10/100,
		             col = 'black', lty = 1, lwd = 1.5)
		    segments(sp,(max(contig_depth$count)*10/100), # vertical stop
		             x1 = sp,
		             y1 = min(contig_depth$count)*10/100,
		             col = 'black', lty = 1, lwd = 1.5)
		    
		    AT_lineage = st + ((sp - st)/2)      
		    if (AT_lineage < 300) {AT_lineage = 500}
		    #mtext(nameEVE,cex=0.45,side = 1, at=st+((sp-st)/2), line=-1.8, col = 'grey')
		    mtext(lineageEVE, cex = 0.50 ,side = 1, at = AT_lineage, line = -1.5, col = 'gray48')
		    
		    
		  	#}
		}
 		
    		lines(x = xl, y = c(0,0), col = "black")
        #mtext(vir,cex=0.65,line=0, col = 'grey')
        title(contig)

        if ((length(args) >= 6) && (file.exists(args[6]))) {
          if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) {
          
          # segments(contig_TE$start,yl[1] - yl[1]*10/100, # trait horiz
          #          x1 = contig_TE$stop,lwd = 2,col = arrowCol[as.character(contig_TE$type)],
          #          y1 = yl[1] - yl[1]*10/100, lty = 1 )
            
            segments(contig_TE$start,-yy - yy*10/100, # trait horiz
                     x1 = contig_TE$stop,lwd = 2,col = arrowCol[as.character(contig_TE$type)],
                     y1 = -yy - yy*10/100, lty = 1 )
            }
        }
		}
	}

	###https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r
	
	
	print (nrow(EVE))
	print (paste('nbExp:',nrow(df), args[1]))
	
	######################################"###################"###################"
	###################"###################"###################"###################"###################"
	###################"###################"###################"###################"###################"###################"
	###################"###################"###################"###################"###################"###################"###################"
	EVE_combin = EVE_combin[EVE_combin$Contig %in% sevEves,]
    
    par(mfrow = c(3,1))

    for (i in 1:length(EVE_combin$Contig)) { ## parcours tous les Contigs
      contig = as.character(EVE_combin[i,1])
      Start = EVE_combin[i,2]
      Stop = EVE_combin[i,3]
      startR = EVE_combin[i,2] - 1000
      stopR = EVE_combin[i,3] + 1000
      vir = EVE_combin$lineage[i] ## pb

      contig_depth = table_FG[which(as.character(table_FG$contig) == contig
                                   & table_FG$position > startR
                                   & table_FG$position < stopR),] # pro/locus de l'EVE

      eve_depth = table_FG[which(as.character(table_FG$contig) == contig
                                & table_FG$position > Start
                                & table_FG$position < Stop),]

      eve_small = eve_depth[which(as.character(eve_depth$contig) == contig
                                 & (as.character(eve_depth$fill) == 's_forw'
                                    | as.character(eve_depth$fill) == 's_rev')),] # pro/locus de l'EVE



      contig_depth = contig_depth[order(contig_depth$position, decreasing = F),]
      ## test

      ### DISPLAY SELON NB EVE DANS CONTIG
      if (is.data.frame(contig_depth) && nrow(eve_small) < 100 ) {next} #
      if (any(eve_small$count > 10) == TRUE) {
  
        starts = EVE$start[which(as.character(EVE$Contig) == contig
        & EVE$start >= Start & EVE$stop <= Stop) ];starts = sort(starts, decreasing = FALSE) # recupere ts les starts des EVE du contig
        stops = EVE$stop[which(as.character(EVE$Contig) ==  contig
                              & EVE$start >= Start & EVE$stop <= Stop)]; stops = sort(stops, decreasing = FALSE)    # recupere ts les stops des EVE du contig
        legende = legende + 1
        vec = c('s_forw','s_rev','l_forw','l_rev') #,
        xl = c(startR ,stopR)
  
  
        if (xl[1] < 0 ) {xl[1] = 0}; if (xl[2] < 0 ) {xl[2] = 0} # gestion erreur plot size vs genome size
  
        m <- (min(contig_depth$count) + 1)*110/100
        if (min(contig_depth$count) > -(max(contig_depth$count)/3)) {m <- -(max(contig_depth$count)/2) }
        yl = c(min(contig_depth$count) + m,max(contig_depth$count) + (max(contig_depth$count)*5/100))
        y = 0; x = 0
        ## PLOTS
        yy = max(abs(contig_depth$count))
        plot( y ~ x,
              las = 1, type = 'n', xpd = FALSE,
              frame = F, xlab = "",ylab = "", ylim = c(-(yy + (yy*0.4)), yy + (yy*0.4)),
              cex.lab = cexLab,cex.axis = cexAxis,
              xlim = xl #,ylim = c(yl[1],yl[2] + (yl[2]*0.4))
        )
        for (v in 1:length(vec)) {	# pour chaque type arn plots
          x = contig_depth$position[which(contig_depth$fill == vec[v]  )]
          y = contig_depth$count[which(contig_depth$fill == vec[v] )]
          y = c(0,y,0) 
          x = c(x[1] - 1,x,x[length(x)] + 1)
          polygon( y ~ x,
                    las = 1, col = gb2[v], border = NA ,
                    cex.lab = 2,cex.axis = cexAxis,lwd = 0.7, main = paste(contig,'\n'),
                    xlim = xl, ylim = yl)
        }
  
        if (te == T) {
          contig_TE = TE[which(as.character(TE$Contig) == contig
                              & ((TE$start > (startR)  & (TE$start < (stopR)) # &  TE$stop > (startR)
                                    | ((TE$stop > (startR) &  TE$stop < (stopR)))
                              ))),] # pro/locus de l'EVE
  
          if (te == T) {
            if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) {
              contig_TE$start[contig_TE$start < startR & contig_TE$stop > startR  ] <- startR #+1000
              contig_TE$stop[contig_TE$stop > stopR & contig_TE$start < stopR] <- stopR
            }
          }
        }
  
  
        for (j in 1:length(starts)) { ## pour chaque EVEs du contig
          st =	starts[j]
          sp = stops[j]
          nameEVE = EVE$EVE[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
                                & EVE$start == starts[j]
                                & EVE$stop == stops[j] )]
          lineageEVE = EVE$lineage[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
                                       & EVE$start == starts[j]
                                       & EVE$stop == stops[j] )]
  
          eSens = EVE$stranded[which(as.character(EVE$Contig) == contig   # SOUSTABLE EVE_flanq
                                   & EVE$start == st
                                   & EVE$stop == sp )]
          
          if (te == T) {
            if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) { 
                te_todel = contig_TE[contig_TE$start >= st  & contig_TE$stop <= sp,  ]
                contig_TE = contig_TE[!(contig_TE$start %in% te_todel$start & contig_TE$stop %in% te_todel$stop ),]
                
                te_todel = contig_TE[contig_TE$start <= st  & contig_TE$stop >= sp,  ]
                contig_TE = contig_TE[!(contig_TE$start %in% te_todel$start & contig_TE$stop %in% te_todel$stop ),]
                
                contig_TE$start[contig_TE$start > st  & contig_TE$start < sp] <- sp + 3
                contig_TE$stop[contig_TE$stop > st & contig_TE$stop < sp] <- st - 3
              }
            }
          tmp = sp
          if (as.character(eSens) == '-') {
            tmp = sp
            sp = st
            st = tmp
          }


          
          # arrows(st, yl[1] - yl[1]*10/100, x1 = sp, y1 = yl[1] - yl[1]*10/100, length = 0.045, angle = 40,
          #        code = 2, col = par("fg"), lty = par("lty"),
          #        lwd = par("lwd"))
          # 
          arrows(st, yy - yy*10/100, x1 = sp, y1 = yy - yy*10/100, length = 0.045, angle = 40,
                 code = 2, col = par("fg"), lty = par("lty"),
                 lwd = par("lwd"))
          
          segments(st,(max(contig_depth$count)*10/100), # vertical start
                   x1 = st,
                   y1 = min(contig_depth$count)*10/100,
                   col = 'black', lty = 1, lwd = 1.5)
          segments(sp,(max(contig_depth$count)*10/100), # vertical stop
                   x1 = sp,
                   y1 = min(contig_depth$count)*10/100,
                   col = 'black', lty = 1, lwd = 1.5)

          
          
          mtext(lineageEVE,cex = 0.45,side = 1, at = st + ((sp - st)/2), line = -1.8, col = 'grey')
          #}
        }
        lines(x = xl, y = c(0,0), col = "black")
        title(contig)
  
  
        if ((length(args) >= 6) && (file.exists(args[6]))) {
          if (is.data.frame(contig_TE) && nrow(contig_TE) > 0) {
  
            # segments(contig_TE$start,yl[1] - yl[1]*10/100, # trait horiz
            #          x1 = contig_TE$stop,lwd = 2,col = arrowCol[as.character(contig_TE$type)],
            #          y1 = yl[1] - yl[1]*10/100, lty = 1 )
            
            segments(contig_TE$start,yy - yy*10/100, # trait horiz
                     x1 = contig_TE$stop,lwd = 2,col = arrowCol[as.character(contig_TE$type)],
                     y1 = yy - yy*10/100, lty = 1 )
          }
        }
      }
    }
    write.table(df, file = args[9], col.names = FALSE, row.names = FALSE,sep = "\t", quote = FALSE)
    write.table(df, file = args[10], col.names = FALSE, row.names = FALSE,sep = "\t", quote = FALSE)
    ###################
	print(paste('plot OK, path : ',outFile))
	dev.off()
# save.image('temp.RData')
}


for e in $(cat l); do for i in {FG,FS,MG,MS}; do Rscript Figs2020.R files/tabs/${i}/${e}_${i}.tab files/bed-2020/${e}.txt files/bed-2020/${e}.txt files/len/${e}_GT.fna.len figures2020_v1/${e}_${i} ;done ;done
#Rscript --vanilla eve_mapping_plots.R 4_small_long_expTAB/FG/Aaeg_FG.tab strandedBed_new/Aaeg_stranded.txt strandedBed_new/Aaeg_stranded.txt 4_genomes_transcriptomes_len/Aaeg_GT.fna.len Aaeg_FG.pdf 4_ET_out/Aaeg.out 
