#read in count table and biotypes 
args = commandArgs(trailingOnly =T)
type = args[1]
value = args[2]
base_dir = args[3]
temp_dir = args[4]
anno_dir = args[5]


countsDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="longRNA")))
countsDirExists <- length(countsDir) == 1
write(paste("countsDirExists: ", countsDirExists), stderr())

rpkmDir <- unique(dirname(list.files(base_dir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE, pattern="gene_counts_rpkmtpm.txt")))
rpkmDirExists <- length(rpkmDir) == 1
write(paste("rpkmDirExists: ", rpkmDirExists), stderr())

biotypes <- read.table(paste(anno_dir,"gene_biotypes.tsv", sep="/"), sep="\t", header=T)
countLong <- read.table(paste(countsDir, "gene_counts_longRNA", sep="/"), sep="\t", header=T, col.names=c("Gene", "Chr", "Start", "End", "Strand", "Length", "Count"))
#countLong <- left_join(countLong, biotypes, by = c("Gene" = "gene_id"))
#counts with biotypes as one table 
colnames(biotypes)[1] = "Gene"
countsWbiotype = merge(countLong, biotypes[biotypes[,"Gene"] %in% countLong[,"Gene"],], by = c("Gene"), all=T)
normal <- read.table(paste(rpkmDir,"gene_counts_rpkmtpm.txt", sep="/"), sep="\t", header=T) 
colnames(normal)[1] = "Gene"
together = merge(countsWbiotype, normal, by="Gene")

write.table(together, "Full_count_table.csv", sep=",", row.names=F)

if(type =="reads"){
	together = together[together$Count >value,]
	write.table(together, "Filter_count_table.csv", sep=",", row.names=F)
}else if( type == "RPKM"){
	together = together[together$RPKM >value,]
	write.table(together, "Filter_count_table.csv", sep=",", row.names=F)
}else if( type == "TPM"){
	together = together[together$TPM >value,]
	write.table(together, "Filter_count_table.csv", sep=",",  row.names=F)
}else{
	write.table(together, "Filter_count_table.csv", sep=",",  row.names=F)
}

