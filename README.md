# EVEs_expression
### Purpose 
This repository contains scripts for the analysis of EVEs (Endogenous Viral Elements) and their mapping in various genomes. The provided scripts facilitate data processing, alignment, visualization, and analysis of EVEs. Please follow the outlined steps below for efficient execution of the analysis pipeline


### Table of Contents
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Directory contents](#Directorycontents)
- [Contributing](#contributing)
- [License](#license)

## Installation
Clone this repository to your local machine and provide necessary permissions to the scripts using the following commands:

```bash
git clone https://github.com/carolebelliardo/EVE_expression.git
chmod +x EVE_expression
```

### Dependencies 
In order to successfully utilize the scripts within this repository, certain external tools are required. Please ensure you have the following tools installed and accessible in your environment:

- [Bowtie2](https://github.com/BenLangmead/bowtie2): A tool for aligning sequencing reads to a reference genome. It provides fast and accurate alignment, suitable for both small and large reads[^1].
- [Hisat2](https://github.com/DaehwanKimLab/hisat2): A fast and sensitive tool for aligning RNA-seq reads to genomes. It's optimized for mapping RNA-seq data and supports both short and long reads[^2].
- [Samtools](https://github.com/samtools/): A powerful tool for working with SAM (Sequence Alignment/Map) files. It provides functions for manipulating, processing, and querying SAM and BAM files[^3].
- [R](https://www.r-project.org/): A programming language for statistical computing and graphics[^4.]
- [Python 3](https://www.python.org/download/releases/3.0/): A programming language required for executing the provided Python scripts[^5].



The following R librairies are required:
- *dplyr*
- *ggplot2*
- *dplyr*
- *devtools*
- *ape*
- *RColorBrewer*
- *reshape*
- *extrafont*
- *wesanderson*
- *ggplot2*
- *gridExtra*


The following python librairies are required:
```
* *pprint*
* *os*
* *sys*
* *optparse*
* *argparse*
* *subprocess*
```



--------
Repository Contents

*Plot_tree.R: Run this script in R Studio to generate a species distribution tree depicting EVEs.
*remove_closeEve.py: This script generates the combined BED file ("bed_combin") by merging EVEs with the ability to set a threshold distance parameter ("distance"). Default value is 50. Input required is a BED file.
*eve_mapping.py: This script parses BAM files to extract EVE mapping information. Input <BAM file> / Output :deep_table <tab file>.
*eve_mapping_plots_family.R: This script generate visualisation related to EVE mapping profils. Input: deep_table <tab file> / Output <PDF file>.
*reads_len.py: Compute the percentage of siRNA/piRNA


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
bowtie2 --threads 25 -x ${GT_path}/${genome}_GT.fna.bwt -U ${fastq_path}/${genome}_${cond1}_sRNA_f_t.fq -S ${sam_path}/${genome}_${cond2}_small.SAM --no-20232023
Private
Privateunal 
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



[^1]: R Core Team (2021). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. [URL: https://www.R-project.org/](https://www.R-project.org/)
[^2]: Python Software Foundation. (2021). Python Language Reference, version 3.9.6. [URL: https://www.python.org/](https://www.python.org/)

