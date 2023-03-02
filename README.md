# EVEs_2021
  - Plot_tree.R à éxecuter dans R studio pour generer l'arbre des espèces et distrib EVE

  - remove_closeEve.py : pour générer le fichier le fichier bed_combin, possibilité de parametrer le nombre de pb séparant les EVE mergés avec la variable "distance". Par default 50 (input = bed file)
  - eve_mapping.py : parsing bam files for extract EVE mapping info (input : bam file -> output : deep_table file)
  - eve_mapping_plots_family.R : drawing plots (input : deep_table, output : plots)
  - reads_len.py : compute si/pi RNA %

 
## absolut location start and stop
for i in *; do Rscript --vanilla ./absolutPosi.R $i '.' ;done # retrouver start et stop des EVE sur les scafffolds

####################
##     Mapping    ##
####################
#-- indexes fastq***
for i in *.fna; do bowtie2-build  --large-index -f${genome} ${genome}.bwt; done #smalls
for i in *.fna; do hisat2-build --large-index -p 25  ${genome} ${genome}.hst; done  #longs

#-- Align
for cond1 in F_G F_S M_G M_S; do 
cond2=$(echo $cond1| sed 's/_//')
bowtie2 --threads 25 -x ${GT_path}/${genome}_GT.fna.bwt -U ${fastq_path}/${genome}_${cond1}_sRNA_f_t.fq -S ${sam_path}/${genome}_${cond2}_small.SAM --no-unal 
done 

for cond1 in F_G F_S M_G M_S; do 
cond2=$(echo $cond1| sed 's/_//')
/mnt/65To/bin/hisat2-2.1.0/hisat2 -p 25 -x ${GT_path}/${genome}_GT.fna.hst -q -1 ${fastq_path}/${genome}_${cond1}_RNAseq_f_t.1.fq -2 ${fastq_path}/${genome}_${cond1}_RNAseq_f_t.2.fq --rna-strandness FR -S ${sam_path}/${genome}_${cond2}_long.SAM --no-unal 
done 


#-- SAM to BAM
for l in small long; do  
  for cond1 in F_G F_S M_G M_S; do 
    /mnt/65To/bin/samtools/bin/samtools view -h -bS ${sam_path}/${genome}_${cond2}_${l}.SAM > ${bam_path}/${genome}_${cond2}_${l}.bam
    /mnt/65To/bin/samtools/bin/samtools sort ${bam_path}/${genome}_${cond2}_${l}.bam -o ${sorted_bam_path}/${genome}_${cond2}_${l}.bam
    /mnt/65To/bin/samtools/bin/samtools index ${sorted_bam_path}/${genome}_${cond2}_${l}.bam 
  done
done

####################
##     Plots      ##
####################

#-- One plot

Rscript --vanilla eve_mapping_plots.R Apis_FG.tab /mnt/65To/Carole/1-data/genomes_transcriptomes/Apis_GT.fna.len Apis_FG.pdf

#-- All plots
for cond1 in FG FS MG MS; do 
python3 eve_mapping.py -f 80000 -g  -pb ${bed_path}/${genome}.bed -sB ${sorted_bam_path}/${genome}_${cond2}_small.bam -lB ${sorted_bam_path}/${genome}_${cond2}_long.bam -o ${tab_path}/${genomes}_${cond1}.tab
done

for cond1 in FG FS MG MS; do 
for i in $(cat list_18_species.txt); do Rscript --vanilla eve_mapping_plots.R ${cond1}_${genome}.tab /mnt/65To/Carole/1-data/genomes_transcriptomes/${cond1}.fna.len ${cond1}_${genome}.pdf ;done ;done' &

## nohup sh -c 'for ...'&