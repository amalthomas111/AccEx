#Author: AT
#Script to filter reads based on whitelist and write to wide format
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
      stop("Usage: Rscript script.R <whitelist_file> <counts-long-file>", call.=FALSE)
}

suppressPackageStartupMessages({
      library(reshape2)
        library(Matrix)
})
whitelist_file = args[1]
if(!file.exists(whitelist_file)) stop("\nWhitelist file not found!",call.=FALSE)
#"SINAA2_default_mod_selected.cells.tsv"
#barcode    reads
#CACAGATGTTCTTAGG-1  120521
#TGATTTCGTACGGAGT-1  109707
#AGGCGAACATCACAGT-1  106708
#ACCCAAACAAATAGTG-1  105128

df_whitelist = read.table(whitelist_file,sep="\t",header=T)
head(df_whitelist)
cat("# unique whitelist:",nrow(df_whitelist),"\n")

read_file = args[2]
if(!file.exists(read_file)) stop("\nReads file not found!",call.=FALSE)
#"SINAA2.default.long.txt"
#gene   cell    count
#r10000_chr10.116949677.116950785    AAAGGATAGTCTCGAT-1  1
#r10000_chr10.116949677.116950785    AAAGGATCAATGAAAC-1  2
#r10000_chr10.116949677.116950785    AACGGGACAGAACGAC-1  1

df_read = read.table(read_file,sep="\t",header=T)
head(df_read)
dim(df_read)

outputname = gsub(".long.txt","",basename(read_file))
outputname = gsub(".txt","",outputname)
cat("\nOutputname:",outputname,"\n")

cat("\nReading long counts file\n")
df_read = df_read[df_read$cell %in% df_whitelist$barcode,]
cat("\nAfter filtering using whitelist\n")
dim(df_read)

if(nrow(df_whitelist) != length(unique(df_read$cell))){
    stop("\nLength of unique Whitelist barcode not same! Exiting\n")
}
cat("\nconverting to long")
df = dcast(df_read,gene~cell,value.var = "count")
row.names(df) = df$gene
df$gene = NULL
df[is.na(df)] = 0
#write.table(df,file = paste0("countstable/",outputname,".tsv"),sep ="\t",quote = F,col.names = NA)

cat("\nConverting to matrix")
df = as.matrix(df)
dim(df)

mainDir=getwd()
subDir="countstable"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
cat("\nconverting to sparse")
df.sparse = Matrix(df,sparse = T)

#write sparse matrix files

writeMM(obj = df.sparse,file = paste0("countstable/",outputname,".mtx"))
write(x=row.names(df),file = paste0("countstable/",outputname,"_regions.tsv"))
write(x=colnames(df),file = paste0("countstable/",outputname,"_barcodes.tsv"))
