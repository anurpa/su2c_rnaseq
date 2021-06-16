#Pavana Anur
#R Script : Plots for Baylor & Geogetown data 

#Format files to include UCSC gene symbols 
ucsc<- read.table("UCSC_Names.txt",sep="\t",header=T)

#BT474AZ_PvsELR
fpkm<-read.table("BT474AZ_PvsELR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_PvsELR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_PvsELR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values from Tumor and Normal and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_P FPKM (log2)", ylab="BT474AZ_ELR FPKM (log2)", main="BT474AZ_P vs BT4
74AZ_ELR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N  and display them on the plot
topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2
#Order the output by fold-change or p-value
sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]
write.csv(output, file="BT474AZ_PvsELR_sigDe.csv", row.names=FALSE)

##################################################################
#BT474AZ_PvsLLR
fpkm<-read.table("BT474AZ_PvsLLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_PvsLLR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_PvsLLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values from Tumor and Normal and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_P FPKM (log2)", ylab="BT474AZ_LLR FPKM (log2)", main="BT474AZ_P vs BT4
74AZ_LLR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]
write.csv(output, file="BT474AZ_PvsLLR_sigDe.csv", row.names=FALSE)

#############################################################################
#BT474AZ_PvsLTR

fpkm<-read.table("BT474AZ_PvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_PvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_PvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_P FPKM (log2)", ylab="BT474AZ_LTR FPKM (log2)", main="BT474AZ_P vs BT4
74AZ_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)


sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]
write.csv(output, file="BT474AZ_PvsLTR_sigDe.csv", row.names=FALSE)

#####################################################################
#BT474AZ_PvsTR

fpkm<-read.table("BT474AZ_PvsTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_PvsTR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_PvsTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_P FPKM (log2)", ylab="BT474AZ_TR FPKM (log2)", main="BT474AZ_P vs BT47
4AZ_TR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)


sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="BT474AZ_PvsLTR_sigDe.csv", row.names=FALSE)
#####################################################################
#BT474AZ_TRvsLLR

fpkm<-read.table("BT474AZ_TRvsLLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_TRvsLLR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_TRvsLLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_TR FPKM (log2)", ylab="BT474AZ_LLR FPKM (log2)", main="BT474AZ_TR vs B
T474AZ_LLR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)



sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="BT474AZ_TRvsLLR_sigDe.csv", row.names=FALSE)
#####################################################################
#BT474AZ_ELRvsLLR

fpkm<-read.table("BT474AZ_ELRvsLLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_ELRvsLLR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_ELRvsLLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_ELR FPKM (log2)", ylab="BT474AZ_LLR FPKM (log2)", main="BT474AZ_ELR vs
 BT474AZ_LLR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:10]
topn = order(de[,"q_value"])[1:10]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="BT474AZ_ELRvsLLR_sigDe.csv", row.names=FALSE)
#####################################################################
#BT474AZ_TRvsLTR

fpkm<-read.table("BT474AZ_TRvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"BT474AZ_TRvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("BT474AZ_TRvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="BT474AZ_TR FPKM (log2)", ylab="BT474AZ_LTR FPKM (log2)", main="BT474AZ_TR vs B
T474AZ_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a  table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="BT474AZ_TRvsLTR_sigDe.csv", row.names=FALSE)

#####################################################################
#####################################################################

#SKBR3_PvsLTR

fpkm<-read.table("SKBR3_PvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"SKBR3_PvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("SKBR3_PvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="SKBR3_P FPKM (log2)", ylab="SKBR3_LTR FPKM (log2)", main="SKBR3_P vs SKBR3_LTR
 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)


sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="SKBR3_PvsLTR_sigDe.csv", row.names=FALSE)

#####################################################################
#SKBR3_PvsTR

fpkm<-read.table("SKBR3_PvsTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"SKBR3_PvsTR_de.txt",row.names=F,sep='\t')

de<-read.table("SKBR3_PvsTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="SKBR3_P FPKM (log2)", ylab="SKBR3_TR FPKM (log2)", main="SKBR3_P vs SKBR3_TR F
PKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a  table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="SKBR3_PvsTR_sigDe.csv", row.names=FALSE)

#####################################################################

#SKBR3_TRvsLR

fpkm<-read.table("SKBR3_TRvsLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"SKBR3_TRvsLR_de.txt",row.names=F,sep='\t')

de<-read.table("SKBR3_TRvsLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="SKBR3_TR FPKM (log2)", ylab="SKBR3_LR FPKM (log2)", main="SKBR3_TR vs SKBR3_LR
 FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="SKBR3_TRvsLR_sigDe.csv", row.names=FALSE) 

#####################################################################

#SKBR3_TRvsLTR

fpkm<-read.table("SKBR3_TRvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"SKBR3_TRvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("SKBR3_TRvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="SKBR3_TR FPKM (log2)", ylab="SKBR3_LTR FPKM (log2)", main="SKBR3_TR vs SKBR3_L
TR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="SKBR3_TRvsLTR_sigDe.csv", row.names=FALSE) 
#####################################################################################
#####################################################################################

#UACC812_PvsLTR

fpkm<-read.table("UACC812_PvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"UACC812_PvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("UACC812_PvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="UACC812_P FPKM (log2)", ylab="UACC812_LTR FPKM (log2)", main="UACC812_P vs UAC
C812_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="UACC812_PvsLTR_sigDe.csv", row.names=FALSE)

#####################################################################

#UACC812_PvsTR

fpkm<-read.table("UACC812_PvsTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"UACC812_PvsTR_de.txt",row.names=F,sep='\t')

de<-read.table("UACC812_PvsTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="UACC812_P FPKM (log2)", ylab="UACC812_TR FPKM (log2)", main="UACC812_P vs UACC
812_TR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:10]
topn = order(de[,"q_value"])[1:10]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="UACC812_PvsTR_sigDe.csv", row.names=FALSE)
#####################################################################

#UACC812_PvsLR

fpkm<-read.table("UACC812_PvsLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"UACC812_PvsLR_de.txt",row.names=F,sep='\t')

de<-read.table("UACC812_PvsLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="UACC812_P FPKM (log2)", ylab="UACC812_LR FPKM (log2)", main="UACC812_P vs UACC
812_LR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:10]
topn = order(de[,"q_value"])[1:10]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="UACC812_PvsLR_sigDe.csv", row.names=FALSE)


#####################################################################

#UACC812_TRvsLR

fpkm<-read.table("UACC812_TRvsLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"UACC812_TRvsLR_de.txt",row.names=F,sep='\t')

de<-read.table("UACC812_TRvsLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="UACC812_TR FPKM (log2)", ylab="UACC812_LR FPKM (log2)", main="UACC812_TR vs UA
CC812_LR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="UACC812_TRvsLR_sigDe.csv", row.names=FALSE) 

#####################################################################

#UACC812_TRvsLTR

fpkm<-read.table("UACC812_TRvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"UACC812_TRvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("UACC812_TRvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="UACC812_TR FPKM (log2)", ylab="UACC812_LTR FPKM (log2)", main="UACC812_TR vs U
ACC812_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)



sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="UACC812_TRvsLTR_sigDe.csv", row.names=FALSE) 

#####################################################################
#####################################################################

#HCC1954_PvsLTR

fpkm<-read.table("HCC1954_PvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"HCC1954_PvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("HCC1954_PvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="HCC1954_P FPKM (log2)", ylab="HCC1954_LTR FPKM (log2)", main="HCC1954_P vs HCC
1954_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="HCC1954_PvsLTR_sigDe.csv", row.names=FALSE)

#####################################################################

#HCC1954_PvsTR

fpkm<-read.table("HCC1954_PvsTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"HCC1954_PvsTR_de.txt",row.names=F,sep='\t')

de<-read.table("HCC1954_PvsTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="HCC1954_P FPKM (log2)", ylab="HCC1954_TR FPKM (log2)", main="HCC1954_P vs HCC1
954_TR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:10]
topn = order(de[,"q_value"])[1:10]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a  table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="HCC1954_PvsTR_sigDe.csv", row.names=FALSE)
#####################################################################

#HCC1954_PvsLR

fpkm<-read.table("HCC1954_PvsLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"HCC1954_PvsLR_de.txt",row.names=F,sep='\t')

de<-read.table("HCC1954_PvsLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="HCC1954_P FPKM (log2)", ylab="HCC1954_LR FPKM (log2)", main="HCC1954_P vs HCC1
954_LR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:10]
topn = order(de[,"q_value"])[1:10]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a  table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="HCC1954_PvsLR_sigDe.csv", row.names=FALSE)


#####################################################################

#HCC1954_TRvsLR

fpkm<-read.table("HCC1954_TRvsLR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"HCC1954_TRvsLR_de.txt",row.names=F,sep='\t')

de<-read.table("HCC1954_TRvsLR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="HCC1954_TR FPKM (log2)", ylab="HCC1954_LR FPKM (log2)", main="HCC1954_TR vs HC
C1954_LR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)

#Write a table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2

sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="HCC1954_TRvsLR_sigDe.csv", row.names=FALSE) 

#####################################################################

#HCC1954_TRvsLTR

fpkm<-read.table("HCC1954_TRvsLTR.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"HCC1954_TRvsLTR_de.txt",row.names=F,sep='\t')

de<-read.table("HCC1954_TRvsLTR_de.txt",header=T,sep="\t")
sig = which(de[,"significant"]=="yes")
#Plot the grand expression values  and mark those that are differentially expressed
x=log2(de[,"value_1"])
y=log2(de[,"value_2"])
plot(x=x, y=y, pch=16, cex=0.25, xlab="HCC1954_TR FPKM (log2)", ylab="HCC1954_LTR FPKM (log2)", main="HCC1954_TR vs H
CC1954_LTR FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

topn = order(abs(de[,"log2.fold_change."]), decreasing=TRUE)[1:15]
topn = order(de[,"q_value"])[1:15]
text(x[topn], y[topn], de[topn,"hg19.kgXref.geneSymbol"], col="black", cex=0.75, srt=45)



sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_c
hange.","q_value")]

write.csv(output, file="HCC1954_TRvsLTR_sigDe.csv", row.names=FALSE) 

