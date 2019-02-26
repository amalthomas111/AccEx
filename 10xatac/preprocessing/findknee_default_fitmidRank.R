#Author: AT
#Script to find the inflection point
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1)
{
  stop("Usage: Rscript script.R <counts>", call.=FALSE)
}
suppressPackageStartupMessages({
 library(Matrix)
library(stats)
})
file=args[1]
outputname= basename(file)
outputname = gsub("_barcode.counts.sorted.txt","",outputname)
outputname= gsub(".barcode.counts.sorted.txt","",outputname)
outputname = gsub("_barcode.counts.txt","",outputname)
outputname= gsub(".barcode.counts.txt","",outputname)
cat("outputname:",outputname,"\n")
reads_per_barcode = read.table(file,col.names = c("barcode","reads"))
o = order(reads_per_barcode[,2],decreasing=TRUE)
ordered.freq = reads_per_barcode[,2][o]
freq.sum = rle(ordered.freq)
mid.rank = cumsum(freq.sum$lengths) - (freq.sum$lengths - 1)/2
barcode.values = freq.sum$values

#log.barcode.values = log10(freq.sum$values)
#log.barcodes.rank = log10(mid.rank)
#lower.counts = 10
#keep = barcode.values > lower.counts
y = log10(barcode.values)
x = log10(mid.rank)

#fit smooth cubic spline
fit1 <- smooth.spline(x=x, y=y, cv = TRUE)

x.p =predict(fit1)$x
y.p = predict(fit1)$y
d1 = predict(fit1,deriv = 1)$y
inflection.g = which.min(d1)
inflection.g
x.p[inflection.g]


threshold = 10^y.p[inflection.g]
selected = reads_per_barcode[reads_per_barcode[,2] > threshold,]
no.of.cells = nrow(selected)
median.reads = median(selected$reads)
mean.reads = round(mean(selected$reads),2)
lowest = tail(selected,n=1)[,2]

mainDir=getwd()
subDir="midrank_inflection"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="selected"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

png(filename = paste0("midrank_inflection/",outputname,".png") )
plot(x=x.p,y=y.p,xlab = "log10(barcode Rank)",ylab= "log10(# barcodes)")
abline(v=x.p[inflection.g], col="magenta", lwd=2,lty=2)
dev.off()

out1 = "fitmidrank_stats.tsv"
line1 =  paste(outputname, no.of.cells, median.reads,mean.reads,lowest,sep="\t")
if(!file.exists(out1)){
write(paste("name","no.of.cells","median_reads","mean_reads","lowest",sep="\t"),file=out1)
}
write(line1,file = out1,append=T)

write.table(selected,file=paste0("selected/",outputname,"_selected.cells.tsv"),
	sep="\t",row.names=F,quote=F)
