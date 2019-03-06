args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1)
{
  stop("Usage: Rscript script <bam>", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input file not found!\
Usage: Rscript peakannot.R <bamfile>", call.=FALSE)
}
attach(mtcars)
mainDir=getwd()
subDir="preseq"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="preseq/plots"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="preseq/plots/subset"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

name = basename(args[1])
#name=gsub(".sorted.bam","",name)
#name=gsub(".potential.bam","",name)
name=gsub(".bam","",name)
cat("\noutputname:",name)

command1=paste0("preseq lc_extrap -o preseq/",name,".extrap.txt -B ",args[1])
command1
command2=paste0("preseq c_curve -o preseq/",name,".txt -B ",args[1])
command2

print("preseq extrap started")
if (!file.exists(paste0("preseq/",name,".extrap.txt"))){
system(command1)
}else{
cat(paste0("preseq/",name,".extrap.txt")," exists")
}

print("preseq extrap done! Curve started")
if (!file.exists(paste0("preseq/",name,".txt"))){
	system(command2)
}else{
cat(paste0("preseq/",name,".txt")," exists")
}
print("curve done")

plotcomplexity = function(name,x,y,x1,y1,present_size){

png(filename = paste0("preseq/plots/",name,"_combined.png"))#,width=400,height=350,res=300)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plot(x,y,xlab="Total read in M",ylab="Unique reads in M",main="Present")
plot(x1,y1,xlim=c(0,1000),ylim=c(0,400),xlab="Total read in M",ylab="Expected Unique reads in M",main="Expected")
title(paste0("Library Complexity ",name),outer = T)
abline(v=present_size,col = "blue")
abline(h=30,col="red")
dev.off()

pdf(file = paste0("preseq/plots/",name,"_combined.pdf"))#,width=400,height=350,res=300)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plot(x,y,xlab="Total read in M",ylab="Unique reads in M",main="Present")
plot(x1,y1,xlim=c(0,1000),ylim = c(0,400),xlab="Total read in M",ylab="Expected Unique reads in M",main="Expected")
title(paste0("Library Complexity ",name),outer = T)
abline(v=present_size,col = "blue")
abline(h=30,col="red")
dev.off()

pdf(file = paste0("preseq/plots/",name,"_present.pdf"))#,width=400,height=350,res=300)
plot(x,y,xlab="Total read in M",ylab="Unique reads in M",main="Present")
dev.off()

pdf(file = paste0("preseq/plots/",name,"_expected.pdf"))#,width=400,height=350,res=300)
plot(x1,y1,xlim=c(0,1000),ylim = c(0,400),xlab="Total read in M",ylab="Expected Unique reads in M",main="Expected")
title(paste0("Library Complexity ",name),outer = T)
abline(v=present_size,col = "blue",lty=2)
#abline(h=30,col="red")
dev.off()


}
oneM=1000000
complexity= read.table(paste0("preseq/",name,".txt"),header = T)
x=complexity$total_reads/oneM
y=complexity$distinct_reads/oneM
present_size=system(paste0("samtools view -F 0x904 -c ",args[1]),intern = T)
present_size=as.numeric(present_size)/oneM
extrap=read.table(paste0("preseq/",name,".extrap.txt"),header = T)#,nrows = 1000)

x1=extrap$TOTAL_READS/oneM
y1=extrap$EXPECTED_DISTINCT/oneM
plotcomplexity(paste0(name,"_complexity"),x,y,x1,y1,present_size)

extrap1=read.table(paste0("preseq/",name,".extrap.txt"),header = T,nrows = 1000)
x1=extrap1$TOTAL_READS/oneM
y1=extrap1$EXPECTED_DISTINCT/oneM
plotcomplexity(paste0("subset/",name,"_complexity_1Krows"),x,y,x1,y1,present_size)

extrap1=read.table(paste0("preseq/",name,".extrap.txt"),header = T,nrows = 300)
x1=extrap1$TOTAL_READS/oneM
y1=extrap1$EXPECTED_DISTINCT/oneM
plotcomplexity(paste0("subset/",name,"_complexity_300rows"),x,y,x1,y1,present_size)
