# Pipeline to process 10x atac-seq
## process_10xatac.snakemake
Input files: 
1. bam files outputed by 10x in bams folder
2. list.txt with names of bam files without .bam extension
3. config_file: "atac_config.yaml"
4. Scripts needed
 - Rscripts: getinsertsize.R, findknee_default_fitmidRank.R, filtercells_2wide.R, plotPreseq.R
 - BASH: mapping_stats.sh

Software Requirements:
picard, samtools, BEDtools, deeptools, MACS2, featurCounts,umi-tools, R, python3
