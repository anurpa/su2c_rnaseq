#Author: Pavana Anur
#R code for formatting cuffdiff output files for SU2C 

#Read gene names from ucsc.txt file
ucsc<- read.table("ucsc.txt",sep="\t",header=T)
names(ucsc)<-c("hg19.knownGene.chrom","gene_id","hg19.kgXref.geneSymbol")
str(ucsc)

#Add gene names to fpkm files
fpkm<-read.table("outputFile.fpkm",header=T,sep='\t')
mergeFile<-merge(fpkm,ucsc)
write.table(mergeFile,"outputFile_de.txt",row.names=F,sep='\t')

#Write out de table
de<-read.table("outputFile_de.txt",header=T,sep="\t")

de$fold_change_computed<-((log2(de$value_2+1))-(log2(de$value_1+1)))

write.table(de,file="outputFile_de_v2.txt",sep="\t",row.names=FALSE)

#Identify top significant differentially expressed genes and write to a file. p <0.05 , log2.fold_change>2
sig = which(de[,"significant"]=="yes" & abs(de[,"log2.fold_change."]) >= 2)
sig_de = de[sig,]
o = order(abs(sig_de[,"q_value"]), decreasing=FALSE)
output = sig_de[o,c("gene_id","hg19.kgXref.geneSymbol","locus","sample_1","sample_2","value_1","value_2","log2.fold_change.","fold_change_
computed","q_value")]
str(output)
write.table(output, file="outputFile_sigDe_v2.txt", row.names=FALSE,sep='\t')
