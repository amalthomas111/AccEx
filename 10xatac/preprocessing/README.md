# Pipeline to process 10x atac-seq
## process_10xatac.snakemake
Input files: 
1. bam files outputed by 10x in bams folder
2. list.txt with names of bam files without .bam extension
3. config_file: "atac_config.yaml"
4. Rscripts: getinsertsize.R, findknee_default_fitmidRank.R, filtercells_2wide.R

Software Requirements:
picard, samtools, BEDtools, deeptools, MACS2, featurCounts,umi-tools, R, python3
