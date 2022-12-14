argv <- commandArgs(TRUE)
start <- as.numeric(argv[1])
####haplotype
# hap.cs = read.table('case.haplotype_output_file',fill = T,col.names = c(1:12))
# hap.pa = read.table('case_pa.haplotype_output_file',fill = T,col.names = c(1:12))
# hap.ma = read.table('case_ma.haplotype_output_file',fill = T,col.names = c(1:12))
# hap.cs = read.table('/home/ccwu/gpu/2021_downsy/haplotype/case/case.haplotype_output_file',fill = T,col.names = c(1:12))
# hap.pa = read.table('/home/ccwu/gpu/2021_downsy/haplotype/case_paternal_2/case_pa.haplotype_output_file',fill = T,col.names = c(1:12))
# hap.ma = read.table('/home/ccwu/gpu/2021_downsy/haplotype/case_maternal_2/case_ma.haplotype_output_file',fill = T,col.names = c(1:12))
##############################################################################################################################################################################
# ch = c(1:22,'X','Y')
# i = start
# ci = ch[i]
# chr = paste0('chr',ch[i])
# temp.cs = hap.cs[which(hap.cs[,4] == chr),]
# temp.pa = hap.pa[which(hap.pa[,4] == chr),]
# temp.ma = hap.ma[which(hap.ma[,4] == chr),]
# temp.pos = unique(c(temp.cs[,5],temp.pa[,5],temp.ma[,5]))
# bp = read.table(paste0('/home/ccwu/haplotype/haplotype/whole_genome/chr/chr',ci,'.uniq.bedpe'))
###################
# n1 = apply(bp,1,function(x) length(intersect(temp.pos,c((as.numeric(x[2]) + 1):(as.numeric(x[3]))))))
# n2 = apply(bp,1,function(x) length(intersect(temp.pos,c((as.numeric(x[5]) + 1):(as.numeric(x[6]))))))
#########
# n1 = as.numeric(n1)
# n2 = as.numeric(n2)
# save(n1,n2,file = paste0('pair_hit_chr',ci,'.Rdata'))
########
####pair
# pa <- file("/home/ccwu/gpu/2021_downsy/fetal/dlohic/case_alter/02-bedpe/1.uniq.bedpe", "r")
# ep = strsplit(readLines(pa,1),'\t')
###############################################################################################################################################################################
ch = c(1:22,'X','Y')
# for (i in 1:24){
i = start
ci = ch[i]
# chr = paste0('chr',ch[i])
# temp.cs = hap.cs[which(hap.cs[,4] == chr),]
# temp.pa = hap.pa[which(hap.pa[,4] == chr),]
# temp.ma = hap.ma[which(hap.ma[,4] == chr),]
# temp.pos = unique(c(temp.cs[,5],temp.pa[,5],temp.ma[,5]))
bp = read.table(paste0('/home/ccwu/haplotype/haplotype/whole_genome/chr/chr',ci,'.uniq.bedpe'))
##################################
# load(paste0('pair_hit_chr',ci,'.Rdata'))
####n1 n2
# num1 = which(n1 > 0)
# num2 = which(n2 > 0)
# num3 = union(num1,num2)
bt = bp
btseq = as.data.frame(array(,dim = c(nrow(bt),2)))
######################################################################################################
#/home/ccwu/haplotype/haplotype/whole_genome##########################################################
name = paste0('/home/ccwu/haplotype/haplotype/whole_genome/pet/chr',ci,'.pet1.uniq.sam')
sys = paste('wc','-l',name)
long = system(sys,intern = T)
hg = as.numeric(substr(long,1,(nchar(long) - nchar(name) - 1)))
#
pet1 <- file(name,"r")
for (j in 1:hg){
	ep=strsplit(readLines(pet1,1),'\t')
	ll = ep[[1]][1]
	num = which(bt[,7] == ll)
	if (length(num) > 0){
		# btseq[num,1] = ep[[1]][10]
		fin = paste(j,num,ep[[1]][10])
		write.table(fin,file = paste0('chr',ci,'_pet1.txt'),append = T,quote =F,row.names =F,col.names =F)
	}
}
close(pet1)
##################################################
name = paste0('/home/ccwu/haplotype/haplotype/whole_genome/pet/chr',ci,'.pet2.uniq.sam')
sys = paste('wc','-l',name)
long = system(sys,intern = T)
hg = as.numeric(substr(long,1,(nchar(long) - nchar(name) - 1)))
#
pet2 <- file(name,"r")
for (j in 1:hg){
	ep=strsplit(readLines(pet2,1),'\t')
	ll = ep[[1]][1]
	num = which(bt[,7] == ll)
	if (length(num) > 0){
		# btseq[num,2] = ep[[1]][10]
	}
}
close(pet2)
#
save(bt,btseq,file = paste0('hap_pair_chr',ci,'.Rdata'))
# }

#############################################################
# chr = read.table('/home/ccwu/haplotype/haplotype/whole_genome/homo_sapiens/chr6.fa')
# chr6 = chr[-1,]
#source('pickup.sequence.R')
























